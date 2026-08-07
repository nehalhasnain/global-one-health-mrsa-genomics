#!/usr/bin/env python3
"""
05_inequity.py — Phase 4: global inequity in public MRSA/S. aureus genomic
sequence VISIBILITY (the paper's headline), plus AMR-gene profile by income.

Reads the LOCAL processed dataset built by 01, and pulls a small set of World
Bank country covariates live (SSL via certifi). Produces:

  * table 30  visibility by income group (records, records/million, non-human/million)
  * table 31  metadata completeness by income group
  * table 32  submission lag (deposition - collection year) by income group + Kruskal-Wallis
  * table 33  country-level dataset (deliverable): counts, population, covariates, rates
  * table 34  negative-binomial regression: country genome count ~ income (+covariates), offset log(pop)
  * table 35  AMR/virulence gene prevalence by income group (human isolates only) + Wilson CI
  * table 36  country COVERAGE by income group (share of countries with >=1 genome)
  * table 37  contemporary sensitivity (2020-2025): headline visibility metrics by income
  * reports/05_inequity.md

Guardrails baked into the report: visibility != disease burden/prevalence;
gene-by-income reflects deposited-genome composition (lineage + assembly-quality
confounded), not AMR burden by healthcare strength.
"""
from __future__ import annotations

import json
import ssl
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats
from statsmodels.stats.proportion import proportion_confint

try:
    import certifi
    _SSL = ssl.create_default_context(cafile=certifi.where())
except Exception:
    _SSL = ssl.create_default_context()

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed" / "s_aureus_primary_collection_year_2015_2025_feasibility_rows.csv.gz"
TABLES = ROOT / "tables"
REPORTS = ROOT / "reports"

INCOME_ORDER = ["High income", "Upper middle income", "Lower middle income", "Low income"]

WB_INDICATORS = {
    "gdp_pc_usd": "NY.GDP.PCAP.CD",           # GDP per capita (current US$)
    "health_exp_pc_usd": "SH.XPD.CHEX.PC.CD",  # current health expenditure per capita
    "physicians_per_1k": "SH.MED.PHYS.ZS",     # physicians per 1,000
    "hosp_beds_per_1k": "SH.MED.BEDS.ZS",      # hospital beds per 1,000
}


def wb_fetch(url: str, timeout: int = 60):
    with urllib.request.urlopen(url, timeout=timeout, context=_SSL) as r:
        return json.loads(r.read())


def wb_countries() -> pd.DataFrame:
    # Offline-reproducible: prefer the frozen cache (00_freeze_raw_sources.py).
    cache = RAW / "worldbank_countries.csv"
    if cache.exists():
        c = pd.read_csv(cache, dtype=str)
        c = c.rename(columns={"region": "region_wb", "income_group": "income_group_wb"})
        return c[["iso3", "income_group_wb", "region_wb"]]
    data = wb_fetch("https://api.worldbank.org/v2/country/all?format=json&per_page=400")
    rows = []
    for r in data[1]:
        if r.get("region", {}).get("id") == "NA":   # aggregates, not countries
            continue
        rows.append({"iso3": r.get("id"), "income_group_wb": r.get("incomeLevel", {}).get("value"),
                     "region_wb": r.get("region", {}).get("value")})
    return pd.DataFrame(rows)


def wb_indicator(indicator: str, value_name: str, start=2015, end=2023) -> pd.DataFrame:
    # Offline-reproducible: prefer the frozen cache (00_freeze_raw_sources.py).
    cache = RAW / "worldbank_indicators.csv"
    if cache.exists():
        c = pd.read_csv(cache)
        c = c[(c["indicator"] == indicator) & (c["year"] >= start) & (c["year"] <= end)]
        if c.empty:
            return pd.DataFrame(columns=["iso3", value_name])
        c = c.sort_values(["iso3", "year"]).groupby("iso3", as_index=False).tail(1)
        out = c.rename(columns={"value": value_name})[["iso3", value_name]]
        if value_name.startswith("population"):
            out[value_name] = out[value_name].astype("int64")
        return out
    url = (f"https://api.worldbank.org/v2/country/all/indicator/{indicator}"
           f"?format=json&per_page=20000&date={start}:{end}")
    try:
        data = wb_fetch(url, timeout=90)
        rows = [{"iso3": r["countryiso3code"], "year": int(r["date"]), value_name: r["value"]}
                for r in data[1] if r.get("value") is not None]
        df = pd.DataFrame(rows)
        if df.empty:
            return pd.DataFrame(columns=["iso3", value_name])
        df = df.sort_values(["iso3", "year"]).groupby("iso3", as_index=False).tail(1)
        return df[["iso3", value_name]]
    except Exception as e:
        print(f"WARNING: WB indicator {indicator} failed: {e}", file=sys.stderr)
        return pd.DataFrame(columns=["iso3", value_name])


def parse_tokens(s) -> set[str]:
    if pd.isna(s):
        return set()
    return {p.strip().split("=", 1)[0].strip().lower()
            for p in str(s).replace('"', "").split(",") if p.strip()}


def wilson(k, n):
    if not n:
        return (np.nan, np.nan)
    lo, hi = proportion_confint(k, n, alpha=0.05, method="wilson")
    return round(100 * lo, 1), round(100 * hi, 1)


def add_gene_flags(df):
    A = df["AMR_genotypes"].map(parse_tokens)
    V = df["virulence_genotypes"].map(parse_tokens)
    df = df.copy()
    df["mecA"] = df["mecA_positive"].astype(str).str.lower().isin({"true", "1", "1.0", "yes"})
    df["acq_tet"] = A.map(lambda s: bool(s & {"tet(k)", "tet(l)", "tet(m)", "tet(o)"}))
    df["any_erm"] = A.map(lambda s: any(t.startswith("erm(") for t in s))
    df["msrA"] = A.map(lambda s: "msr(a)" in s)
    df["fq_qrdr"] = A.map(lambda s: any(t.startswith(("gyra_", "parc_", "grla_")) for t in s))
    df["PVL"] = V.map(lambda s: bool(s & {"luks-pv", "lukf-pv"}))
    df["scn"] = V.map(lambda s: "scn" in s)
    return df


def income_visibility(df, pop_by_income):
    """Records + records/million + non-human/million by income group."""
    d = df[df["income_group"].isin(INCOME_ORDER)]
    rows = []
    for grp in INCOME_ORDER:
        sub = d[d["income_group"] == grp]
        nonhuman = sub[sub["source_category"].isin(["livestock", "companion_animal", "food", "environment"])]
        pop = pop_by_income.get(grp, np.nan)
        rows.append({
            "income_group": grp, "records": len(sub),
            "mrsa_marker_positive": int(sub["mrsa_marker_positive"].astype(str).str.lower().isin({"true", "1", "1.0", "yes"}).sum()),
            "non_human_records": len(nonhuman),
            "unique_countries": sub["iso3"].nunique(),
            "pop_million": round(pop / 1e6, 1) if pd.notna(pop) else np.nan,
            "records_per_million": round(len(sub) / pop * 1e6, 3) if pd.notna(pop) and pop else np.nan,
            "non_human_per_million": round(len(nonhuman) / pop * 1e6, 4) if pd.notna(pop) and pop else np.nan,
        })
    return pd.DataFrame(rows)


def main():
    if not PROCESSED.exists():
        sys.exit(f"ERROR: {PROCESSED} not found; run 01 first.")
    df = pd.read_csv(PROCESSED, low_memory=False)
    df["collection_year"] = pd.to_numeric(df["collection_year"], errors="coerce")
    df["target_creation_year"] = pd.to_numeric(df["target_creation_year"], errors="coerce")

    print("Fetching World Bank covariates (certifi SSL)...")
    wbc = wb_countries()
    pop = wb_indicator("SP.POP.TOTL", "population")
    cov = pop.copy()
    for name, ind in WB_INDICATORS.items():
        cov = cov.merge(wb_indicator(ind, name), on="iso3", how="left")
    cov = cov.merge(wbc, on="iso3", how="left")
    print(f"  covariates for {len(cov)} countries; pop non-null {cov['population'].notna().mean():.2f}")

    # population per income group (WB universe, not repository)
    pop_by_income = cov.groupby("income_group_wb")["population"].sum().to_dict()

    # ---- country-level aggregation (primary window 2015-2025) --------------
    def country_table(dsub):
        agg = dsub.groupby("iso3").agg(
            records=("iso3", "size"),
            income_group=("income_group", "first"),
            country=("country_parsed", "first"),
        ).reset_index()
        nh = (dsub[dsub["source_category"].isin(["livestock", "companion_animal", "food", "environment"])]
              .groupby("iso3").size().rename("non_human_records"))
        agg = agg.merge(nh, on="iso3", how="left")
        agg["non_human_records"] = agg["non_human_records"].fillna(0).astype(int)
        agg = agg.merge(cov[["iso3", "population"] + list(WB_INDICATORS)], on="iso3", how="left")
        agg["records_per_million"] = (agg["records"] / agg["population"] * 1e6).round(3)
        agg["log_pop"] = np.log(agg["population"])
        return agg

    prim = df  # 2015-2025 primary
    ctab = country_table(prim)
    ctab.sort_values("records", ascending=False).to_csv(TABLES / "33_country_level_dataset.csv", index=False)

    # ---- Table 30: visibility by income ------------------------------------
    vis = income_visibility(prim, pop_by_income)
    vis.to_csv(TABLES / "30_visibility_by_income.csv", index=False)

    # ---- Table 36: coverage (share of countries with >=1 genome) -----------
    present = set(prim["iso3"].dropna())
    cov_rows = []
    for grp in INCOME_ORDER:
        universe = cov[cov["income_group_wb"] == grp]["iso3"].nunique()
        have = cov[(cov["income_group_wb"] == grp) & (cov["iso3"].isin(present))]["iso3"].nunique()
        cov_rows.append({"income_group": grp, "countries_total": universe,
                         "countries_with_genome": have,
                         "coverage_pct": round(100 * have / universe, 1) if universe else np.nan})
    coverage = pd.DataFrame(cov_rows)
    coverage.to_csv(TABLES / "36_country_coverage_by_income.csv", index=False)

    # ---- Table 31: metadata completeness by income -------------------------
    comp_fields = {"host": "host", "source": "isolation_source", "collection_date": "collection_date"}
    comp_rows = []
    d = prim[prim["income_group"].isin(INCOME_ORDER)]
    for grp in INCOME_ORDER:
        sub = d[d["income_group"] == grp]
        row = {"income_group": grp, "n": len(sub)}
        for lab, col in comp_fields.items():
            row[f"{lab}_pct"] = round(100 * sub[col].notna().mean(), 1) if len(sub) else np.nan
        comp_rows.append(row)
    completeness = pd.DataFrame(comp_rows)
    completeness.to_csv(TABLES / "31_metadata_completeness_by_income.csv", index=False)

    # ---- Table 32: submission lag by income --------------------------------
    prim = prim.copy()
    prim["lag_years"] = prim["target_creation_year"] - prim["collection_year"]
    lag_rows = []
    lag_groups = []
    for grp in INCOME_ORDER:
        g = prim.loc[prim["income_group"] == grp, "lag_years"].dropna()
        g = g[g >= 0]
        lag_groups.append(g)
        lag_rows.append({"income_group": grp, "n": int(len(g)),
                         "median_lag_yr": float(g.median()) if len(g) else np.nan,
                         "q1": float(g.quantile(.25)) if len(g) else np.nan,
                         "q3": float(g.quantile(.75)) if len(g) else np.nan,
                         "mean_lag_yr": round(float(g.mean()), 2) if len(g) else np.nan})
    try:
        Hk, pk = stats.kruskal(*[g for g in lag_groups if len(g)])
    except ValueError:
        Hk, pk = np.nan, np.nan
    lag = pd.DataFrame(lag_rows)
    lag["kruskal_H"] = round(Hk, 2) if pd.notna(Hk) else np.nan
    lag["kruskal_p"] = pk
    lag.to_csv(TABLES / "32_submission_lag_by_income.csv", index=False)

    # ---- Table 34: negative binomial regression ----------------------------
    import statsmodels.formula.api as smf
    reg = ctab.dropna(subset=["population", "records"]).copy()
    reg = reg[reg["income_group"].isin(INCOME_ORDER)]
    reg["income_group"] = pd.Categorical(reg["income_group"], categories=INCOME_ORDER, ordered=False)
    formula = "records ~ C(income_group, Treatment(reference='High income'))"
    nb_out = []
    nb_alpha = np.nan
    try:
        # Step 1: Poisson GLM with log(pop) offset.
        pois = smf.glm(formula, data=reg, family=sm.families.Poisson(),
                       offset=reg["log_pop"]).fit()
        mu = pois.mu
        # Step 2: Cameron-Trivedi dispersion alpha (aux OLS through origin).
        aux_y = ((reg["records"] - mu) ** 2 - reg["records"]) / mu
        nb_alpha = float(sm.OLS(aux_y, mu).fit().params.iloc[0])
        nb_alpha = max(nb_alpha, 1e-6)
        # Step 3: NB GLM with fixed alpha -> proper Wald CIs.
        m = smf.glm(formula, data=reg, family=sm.families.NegativeBinomial(alpha=nb_alpha),
                    offset=reg["log_pop"]).fit()
        ci = m.conf_int()
        for name in m.params.index:
            if name == "Intercept":
                continue
            nb_out.append({"term": name, "IRR": round(np.exp(m.params[name]), 3),
                           "ci_low": round(np.exp(ci.loc[name, 0]), 3),
                           "ci_high": round(np.exp(ci.loc[name, 1]), 3),
                           "p_value": m.pvalues[name]})
    except Exception as e:
        print(f"WARNING: NB regression failed: {e}", file=sys.stderr)
    nb = pd.DataFrame(nb_out)
    nb.to_csv(TABLES / "34_negbin_visibility_by_income.csv", index=False)

    # ---- Table 35: AMR/virulence by income (human isolates only) -----------
    h = add_gene_flags(prim[prim["source_category"] == "human_or_clinical"])
    genes = ["mecA", "acq_tet", "any_erm", "msrA", "fq_qrdr", "PVL", "scn"]
    g_rows = []
    for grp in INCOME_ORDER:
        sub = h[h["income_group"] == grp]
        row = {"income_group": grp, "n": len(sub)}
        for gene in genes:
            k = int(sub[gene].sum()); n = len(sub)
            lo, hi = wilson(k, n)
            row[f"{gene}_pct"] = round(100 * k / n, 1) if n else np.nan
            row[f"{gene}_ci"] = f"{lo}-{hi}" if n else ""
        g_rows.append(row)
    gene_income = pd.DataFrame(g_rows)
    gene_income.to_csv(TABLES / "35_amr_by_income_human.csv", index=False)

    # ---- Table 37: contemporary sensitivity (2020-2025) --------------------
    sub2020 = df[df["collection_year"] >= 2020]
    vis2020 = income_visibility(sub2020, pop_by_income)
    vis2020.to_csv(TABLES / "37_visibility_by_income_2020_2025.csv", index=False)

    # ---- Narrative report --------------------------------------------------
    hi = vis[vis.income_group == "High income"].iloc[0]
    low = vis[vis.income_group == "Low income"].iloc[0]
    ratio = round(hi["records_per_million"] / low["records_per_million"], 1) if low["records_per_million"] else np.nan
    L = []
    L.append("# Phase 4 — Global inequity in public MRSA/S. aureus genomic sequence visibility\n")
    L.append(f"Primary window: collection-year 2015-2025 (n={len(df):,}). "
             f"World Bank covariates fetched live (population, GDP pc, health exp pc, physicians/1k, beds/1k).\n")
    L.append("## Visibility by income group (per-million normalised)\n")
    L.append("| income group | records | records/million | non-human/million | unique countries |")
    L.append("|---|---|---|---|---|")
    for _, r in vis.iterrows():
        L.append(f"| {r['income_group']} | {r['records']:,} | {r['records_per_million']} | "
                 f"{r['non_human_per_million']} | {r['unique_countries']} |")
    L.append(f"\n**High-income vs low-income records/million ratio ≈ {ratio}:1.**\n")

    L.append("## Country coverage (share of countries with ≥1 deposited genome)\n")
    L.append("| income group | with genome | total countries | coverage % |\n|---|---|---|---|")
    for _, r in coverage.iterrows():
        L.append(f"| {r['income_group']} | {r['countries_with_genome']} | {r['countries_total']} | {r['coverage_pct']} |")

    L.append("\n## Metadata completeness by income group (%)\n")
    L.append("| income group | n | host | source | collection date |\n|---|---|---|---|---|")
    for _, r in completeness.iterrows():
        L.append(f"| {r['income_group']} | {r['n']:,} | {r['host_pct']} | {r['source_pct']} | {r['collection_date_pct']} |")

    L.append("\n## Submission lag (deposition − collection year) by income\n")
    L.append(f"Kruskal-Wallis p = {pk:.2e}.\n")
    L.append("| income group | n | median lag (yr) | IQR |\n|---|---|---|---|")
    for _, r in lag.iterrows():
        L.append(f"| {r['income_group']} | {r['n']:,} | {r['median_lag_yr']} | {r['q1']}-{r['q3']} |")

    if len(nb):
        L.append("\n## Negative-binomial regression — genome count ~ income, offset log(population)\n")
        L.append("Reference = High income. IRR<1 means fewer deposited genomes per capita.\n")
        L.append("| term | IRR | 95% CI | p |\n|---|---|---|---|")
        for _, r in nb.iterrows():
            term = r["term"].replace("C(income_group, Treatment(reference='High income'))", "").strip("[T.]")
            L.append(f"| {term} | {r['IRR']} | {r['ci_low']}-{r['ci_high']} | {r['p_value']:.2e} |")

    L.append("\n## AMR/virulence gene profile by income (human isolates only)\n")
    L.append("| income group | n | mecA | acq.tet | erm | msr(A) | FQ-QRDR | PVL | scn |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for _, r in gene_income.iterrows():
        L.append(f"| {r['income_group']} | {r['n']:,} | {r['mecA_pct']} | {r['acq_tet_pct']} | "
                 f"{r['any_erm_pct']} | {r['msrA_pct']} | {r['fq_qrdr_pct']} | {r['PVL_pct']} | {r['scn_pct']} |")

    L.append("\n## Contemporary sensitivity (collection-year 2020-2025)\n")
    L.append("| income group | records | records/million | non-human/million |\n|---|---|---|---|")
    for _, r in vis2020.iterrows():
        L.append(f"| {r['income_group']} | {r['records']:,} | {r['records_per_million']} | {r['non_human_per_million']} |")

    L.append("\n## Guardrails\n")
    L.append("- **Visibility ≠ burden.** Per-million rates measure public *sequencing/deposition* intensity, "
             "not MRSA incidence or prevalence.\n"
             "- The gene-by-income profile reflects **which genomes were deposited** (lineage mix + assembly "
             "quality, both income-correlated), not AMR burden by healthcare strength; interpret descriptively.\n"
             "- The regression is over countries **with ≥1 genome**; the ~half of low-income countries with "
             "zero genomes (coverage table) are the strongest inequity signal and sit outside the model.\n")

    (REPORTS / "05_inequity.md").write_text("\n".join(L))
    print("Wrote tables 30-37 and reports/05_inequity.md")
    print(f"High:low records/million ratio ≈ {ratio}:1")


if __name__ == "__main__":
    main()
