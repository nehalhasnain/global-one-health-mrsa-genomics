#!/usr/bin/env python3
"""
04_source_stats.py — Phase 3: inferential comparison of AMR / virulence / MRSA
markers across One Health source categories.

Reads the LOCAL processed feasibility dataset (no network) built by
01_feasibility_ncbi_s_aureus.py:
    data/processed/s_aureus_primary_collection_year_2015_2025_feasibility_rows.csv.gz

For a locked panel of AMR and virulence genes it computes, across the five
informative source categories (human, livestock, companion animal, food,
environment; 'unknown' reported but excluded from tests):

  * prevalence + Wilson 95% CI per gene per source              -> table 26
  * omnibus gene x source test (Pearson chi-square) + BH-FDR q  -> table 27
  * pairwise odds ratio vs human reference (Fisher, Wald 95% CI)-> table 28
  * AMR / virulence gene-burden distribution + Kruskal-Wallis   -> table 29
  * narrative report                                            -> reports/04_source_stats.md

Design choices (defensible for peer review):
  - 'unknown' source excluded from hypothesis tests (sensitivity: reported separately).
  - Multiple-testing correction (Benjamini-Hochberg) applied across the omnibus
    p-values (one per gene) -> q-values.
  - mecA / mecC use the pre-computed boolean columns; all other genes use
    normalized token membership (suffixes like =PARTIAL / =MISTRANSLATION /
    =POINT stripped, matching 01's parse_gene_set).
  - Zero-cell odds ratios use a 0.5 Haldane-Anscombe correction for the OR/CI
    point estimate; the p-value is the (uncorrected) Fisher exact p.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.proportion import proportion_confint

ROOT = Path(__file__).resolve().parent.parent
PROCESSED = ROOT / "data" / "processed" / "s_aureus_primary_collection_year_2015_2025_feasibility_rows.csv.gz"
TABLES = ROOT / "tables"
REPORTS = ROOT / "reports"

# Source categories used in tests, in display order. 'unknown' is excluded.
TEST_SOURCES = ["human_or_clinical", "livestock", "companion_animal", "food", "environment"]
ALL_SOURCES = TEST_SOURCES + ["unknown"]
REF = "human_or_clinical"

# ---------------------------------------------------------------------------
# Gene panel. Each entry: (label, field, matcher)
#   field: 'amr' (AMR_genotypes + AMR_genotypes_core), 'vir' (virulence_genotypes),
#          or 'bool' (a pre-computed boolean column named by `match`).
#   matcher: for 'amr'/'vir' a callable(token_set) -> bool; for 'bool' the column name.
# ---------------------------------------------------------------------------

def _exact(*tokens):
    toks = {t.lower() for t in tokens}
    return lambda s: bool(s & toks)

def _prefix(*prefixes):
    pre = tuple(p.lower() for p in prefixes)
    return lambda s: any(t.startswith(pre) for t in s)

AMR_PANEL = [
    # (label, matcher) — matched against normalized AMR token set
    ("mecA (methicillin)", "bool:mecA_positive"),
    ("mecC (methicillin)", "bool:mecC_positive"),
    ("blaZ (penicillinase)", _exact("blaZ")),
    ("tet(K)", _exact("tet(K)")),
    ("tet(M)", _exact("tet(M)")),
    ("tet(L)", _exact("tet(L)")),
    # NB: tet(38) (99.6%) is intrinsic/chromosomal and deliberately excluded;
    # "acquired tetracycline" = tet(K)/tet(L)/tet(M)/tet(O) only.
    ("any acquired tetracycline [tet(K/L/M/O)]", _exact("tet(K)", "tet(L)", "tet(M)", "tet(O)")),
    ("erm(A)", _exact("erm(A)")),
    ("erm(B)", _exact("erm(B)")),
    ("erm(C)", _exact("erm(C)")),
    ("any MLSb [erm(*)]", _prefix("erm(")),
    ("msr(A)", _exact("msr(A)")),
    ("mph(C)", _exact("mph(C)")),
    ("aac(6')-Ie/aph(2'')-Ia (gentamicin)", _exact("aac(6')-ie/aph(2'')-ia")),
    ("any aminoglycoside [aac/ant/aph]", _prefix("aac(", "ant(", "aph(")),
    ("dfrG (trimethoprim)", _exact("dfrG")),
    ("dfrS1 (trimethoprim)", _exact("dfrS1")),
    ("any trimethoprim [dfr*]", _prefix("dfr")),
    ("fusB (fusidic acid)", _exact("fusB")),
    ("fusC (fusidic acid)", _exact("fusC")),
    ("gyrA/parC QRDR mutation (fluoroquinolone)", _prefix("gyra_", "parc_", "grla_")),
]

VIR_PANEL = [
    ("PVL [lukS-PV/lukF-PV]", _exact("lukS-PV", "lukF-PV")),
    ("scn (IEC / staph. complement inhibitor)", _exact("scn")),
    ("sak (IEC / staphylokinase)", _exact("sak")),
    # chp (IEC chemotaxis inhibitor) is not emitted by AMRFinderPlus in this
    # snapshot (0 tokens) and is therefore not measurable here; excluded.
    ("tst (TSST-1)", _exact("tst")),
    ("sea (enterotoxin A)", _exact("sea")),
    ("seb (enterotoxin B)", _exact("seb")),
    ("sec (enterotoxin C)", _exact("sec", "sec1", "sec2", "sec3", "sec4")),
    ("sed (enterotoxin D)", _exact("sed")),
    ("see (enterotoxin E)", _exact("see")),
    ("eta (exfoliative toxin A)", _exact("eta")),
    ("etb (exfoliative toxin B)", _exact("etb")),
    ("cna (collagen adhesin)", _exact("cna")),
]


def parse_tokens(s) -> set[str]:
    """Normalized lowercase token set; strips =PARTIAL/=MISTRANSLATION/=POINT suffixes."""
    if pd.isna(s):
        return set()
    out = set()
    for part in str(s).replace('"', "").split(","):
        g = part.strip()
        if not g:
            continue
        out.add(g.split("=", 1)[0].strip().lower())
    return out


def build_presence(df: pd.DataFrame) -> pd.DataFrame:
    """Return a DataFrame of boolean gene-presence columns aligned to df.index."""
    core = df["AMR_genotypes_core"] if "AMR_genotypes_core" in df.columns else pd.Series([""] * len(df), index=df.index)
    amr_tokens = pd.Series(
        [parse_tokens(a) | parse_tokens(b) for a, b in zip(df["AMR_genotypes"], core)],
        index=df.index)
    vir_tokens = pd.Series([parse_tokens(v) for v in df["virulence_genotypes"]], index=df.index)

    cols = {}
    for label, matcher in AMR_PANEL:
        if isinstance(matcher, str) and matcher.startswith("bool:"):
            col = matcher.split(":", 1)[1]
            cols[label] = df[col].astype(str).str.lower().isin({"true", "1", "1.0", "yes"})
        else:
            cols[label] = amr_tokens.map(matcher)
    for label, matcher in VIR_PANEL:
        cols[label] = vir_tokens.map(matcher)
    return pd.DataFrame(cols, index=df.index)


def wilson(count: int, n: int):
    if n == 0:
        return (np.nan, np.nan)
    lo, hi = proportion_confint(count, n, alpha=0.05, method="wilson")
    return (round(100 * lo, 2), round(100 * hi, 2))


def odds_ratio(a, b, c, d):
    """OR for [[a,b],[c,d]] with Haldane 0.5 correction if any zero cell. Returns (or, lo, hi)."""
    if min(a, b, c, d) == 0:
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    orr = (a * d) / (b * c)
    se = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    lo = np.exp(np.log(orr) - 1.96 * se)
    hi = np.exp(np.log(orr) + 1.96 * se)
    return round(orr, 3), round(lo, 3), round(hi, 3)


def main():
    if not PROCESSED.exists():
        sys.exit(f"ERROR: processed data not found: {PROCESSED}\nRun 01_feasibility_ncbi_s_aureus.py first.")

    df = pd.read_csv(PROCESSED, low_memory=False)
    df = df[df["source_category"].isin(ALL_SOURCES)].copy()
    presence = build_presence(df)
    df = pd.concat([df, presence], axis=1)
    gene_cols = list(presence.columns)

    n_by_source = df.groupby("source_category").size().reindex(ALL_SOURCES).fillna(0).astype(int)
    print("Source Ns:", n_by_source.to_dict())

    # ---- Table 26: prevalence + Wilson CI per gene per source -------------
    prev_rows = []
    for gene in gene_cols:
        row = {"gene": gene}
        for src in ALL_SOURCES:
            sub = df[df["source_category"] == src]
            k = int(sub[gene].sum())
            n = len(sub)
            lo, hi = wilson(k, n)
            row[f"{src}_n"] = n
            row[f"{src}_pos"] = k
            row[f"{src}_pct"] = round(100 * k / n, 2) if n else np.nan
            row[f"{src}_ci"] = f"{lo}-{hi}" if n else ""
        prev_rows.append(row)
    prev = pd.DataFrame(prev_rows)
    prev.to_csv(TABLES / "26_gene_prevalence_by_source.csv", index=False)

    # ---- Table 27: omnibus gene x source test (5 informative sources) + BH-FDR
    test_df = df[df["source_category"].isin(TEST_SOURCES)]
    omni_rows = []
    for gene in gene_cols:
        ct = pd.crosstab(test_df["source_category"], test_df[gene])
        ct = ct.reindex(index=TEST_SOURCES, columns=[False, True], fill_value=0).astype(int)
        table = ct.values
        min_exp = np.nan
        try:
            chi2, p, dof, expected = stats.chi2_contingency(table)
            min_exp = float(expected.min())
        except ValueError:
            chi2, p, dof = np.nan, np.nan, np.nan
        omni_rows.append({
            "gene": gene, "chi2": round(chi2, 3) if pd.notna(chi2) else np.nan,
            "dof": dof, "p_value": p, "min_expected_count": round(min_exp, 2) if pd.notna(min_exp) else np.nan,
        })
    omni = pd.DataFrame(omni_rows)
    valid = omni["p_value"].notna()
    omni["q_value_bh"] = np.nan
    if valid.any():
        omni.loc[valid, "q_value_bh"] = multipletests(omni.loc[valid, "p_value"], method="fdr_bh")[1]
    omni["significant_q<0.05"] = omni["q_value_bh"] < 0.05
    omni = omni.sort_values("q_value_bh", na_position="last")
    omni.to_csv(TABLES / "27_gene_omnibus_tests.csv", index=False)

    # ---- Table 28: pairwise OR vs human reference -------------------------
    or_rows = []
    ref_df = df[df["source_category"] == REF]
    for gene in gene_cols:
        c = int(ref_df[gene].sum()); d = len(ref_df) - c   # human pos/neg
        for src in [s for s in TEST_SOURCES if s != REF]:
            sub = df[df["source_category"] == src]
            a = int(sub[gene].sum()); b = len(sub) - a       # source pos/neg
            orr, lo, hi = odds_ratio(a, b, c, d)
            try:
                _, p = stats.fisher_exact([[a, b], [c, d]])
            except ValueError:
                p = np.nan
            or_rows.append({
                "gene": gene, "source": src, "source_pos": a, "source_n": a + b,
                "human_pos": c, "human_n": c + d, "odds_ratio_vs_human": orr,
                "ci95_low": lo, "ci95_high": hi, "fisher_p": p,
            })
    ordf = pd.DataFrame(or_rows)
    ordf.to_csv(TABLES / "28_gene_or_vs_human.csv", index=False)

    # ---- Table 29: gene-burden distribution + Kruskal-Wallis --------------
    burden_rows = []
    for metric in ["number_amr_genes", "number_virulence_genes"]:
        vals = pd.to_numeric(df[metric], errors="coerce")
        groups = [vals[df["source_category"] == s].dropna() for s in TEST_SOURCES]
        try:
            H, p = stats.kruskal(*[g for g in groups if len(g)])
        except ValueError:
            H, p = np.nan, np.nan
        for src, g in zip(TEST_SOURCES, groups):
            burden_rows.append({
                "metric": metric, "source": src, "n": int(len(g)),
                "median": float(g.median()) if len(g) else np.nan,
                "q1": float(g.quantile(0.25)) if len(g) else np.nan,
                "q3": float(g.quantile(0.75)) if len(g) else np.nan,
                "mean": round(float(g.mean()), 2) if len(g) else np.nan,
                "kruskal_H": round(H, 2) if pd.notna(H) else np.nan,
                "kruskal_p": p,
            })
    burden = pd.DataFrame(burden_rows)
    burden.to_csv(TABLES / "29_gene_burden_by_source.csv", index=False)

    # ---- Narrative report -------------------------------------------------
    n_sig = int(omni["significant_q<0.05"].sum())
    lines = []
    lines.append("# Phase 3 — AMR / virulence / MRSA comparison across One Health source categories\n")
    lines.append(f"Snapshot input: `{PROCESSED.name}` (collection-year 2015-2025, n={len(df):,}).\n")
    lines.append("Statistics: per-gene omnibus Pearson chi-square across the five informative "
                 "source categories (human, livestock, companion animal, food, environment), "
                 "Benjamini-Hochberg FDR across genes; pairwise odds ratios vs human reference "
                 "(Fisher exact p, Wald 95% CI, Haldane 0.5 correction for zero cells). "
                 "'unknown' source excluded from tests.\n")
    lines.append("## Source denominators\n")
    lines.append("| source | n |\n|---|---|")
    for s in ALL_SOURCES:
        lines.append(f"| {s} | {n_by_source[s]:,} |")
    lines.append(f"\n**{n_sig}/{len(gene_cols)} genes differ significantly by source (q<0.05).**\n")

    lines.append("## Prevalence by source (%), selected markers\n")
    lines.append("| gene | " + " | ".join(TEST_SOURCES) + " |")
    lines.append("|" + "---|" * (len(TEST_SOURCES) + 1))
    focus = ["mecA (methicillin)", "mecC (methicillin)", "PVL [lukS-PV/lukF-PV]",
             "scn (IEC / staph. complement inhibitor)", "sak (IEC / staphylokinase)",
             "tst (TSST-1)", "blaZ (penicillinase)",
             "any acquired tetracycline [tet(K/L/M/O)]", "any MLSb [erm(*)]",
             "aac(6')-Ie/aph(2'')-Ia (gentamicin)"]
    for gene in focus:
        r = prev[prev["gene"] == gene].iloc[0]
        cells = " | ".join(f"{r[f'{s}_pct']}" for s in TEST_SOURCES)
        lines.append(f"| {gene} | {cells} |")

    lines.append("\n## Strongest source differences (top 15 by q-value)\n")
    lines.append("| gene | chi2 | p | q(BH) | min exp |\n|---|---|---|---|---|")
    for _, r in omni.head(15).iterrows():
        lines.append(f"| {r['gene']} | {r['chi2']} | {r['p_value']:.2e} | "
                     f"{r['q_value_bh']:.2e} | {r['min_expected_count']} |")

    lines.append("\n## Gene burden (median [IQR]) by source\n")
    lines.append("| metric | " + " | ".join(TEST_SOURCES) + " | Kruskal p |")
    lines.append("|" + "---|" * (len(TEST_SOURCES) + 2))
    for metric in ["number_amr_genes", "number_virulence_genes"]:
        b = burden[burden["metric"] == metric]
        cells = " | ".join(
            f"{b[b['source']==s]['median'].iloc[0]:.0f} [{b[b['source']==s]['q1'].iloc[0]:.0f}-{b[b['source']==s]['q3'].iloc[0]:.0f}]"
            for s in TEST_SOURCES)
        kp = b["kruskal_p"].iloc[0]
        lines.append(f"| {metric} | {cells} | {kp:.2e} |")

    lines.append("\n## Guardrails\n")
    lines.append("- Prevalence reflects **deposited genomes**, not population prevalence; "
                 "convenience/deposition bias applies, especially to small non-human arms.\n"
                 "- Odds ratios contrast source categories within the public repository; "
                 "they are associations in deposited data, not transmission or causal effects.\n"
                 "- Small denominators (environment n<500, companion animal n<1000) widen CIs; "
                 "interpret those ORs cautiously.\n")

    (REPORTS / "04_source_stats.md").write_text("\n".join(lines))
    print("Wrote tables 26-29 and reports/04_source_stats.md")
    print(f"{n_sig}/{len(gene_cols)} genes significant at q<0.05")


if __name__ == "__main__":
    main()
