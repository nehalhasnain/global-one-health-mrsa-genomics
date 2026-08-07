#!/usr/bin/env python3
"""
07_appendix.py — Technical Appendix: supporting analyses, tables and figures
that document the methods and substantiate the limitations of the main study.

Reads the LOCAL processed dataset built by 01 and the derived tables from 03/04/05,
and pulls the same small World Bank covariate set live (SSL via certifi) used by 05
so that regional per-capita rates use the full regional population universe.

Produces (tables/):
  * table 38  source-classification QA: assignment basis + exemplar raw values
  * table 39  geographic concentration: Gini + cumulative country shares +
              per-One-Health-arm single-country dominance
  * table 40  regional visibility (World Bank region): records, per-million,
              non-human, coverage
  * table 41  window & date-definition sensitivity (source composition stability)
  * table 42  master gene table: prevalence + Wilson CI by source, omnibus chi2/q

Produces (figures/), Okabe-Ito, 300 dpi, .png + .pdf, matching 06_figures.py:
  * figS1_geographic_concentration   Lorenz curve + per-arm single-country dominance
  * figS2_temporal_dynamics          collection- vs deposition-year deposition + non-human share
  * figS3_regional_visibility        genomes per million by World Bank region

Guardrails identical to the main analysis: visibility != burden; concentration and
dominance are properties of the public repository (a convenience sample), not of MRSA
biology; gene prevalence is over deposited genomes only.
"""
from __future__ import annotations

import json
import ssl
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import certifi
    _SSL = ssl.create_default_context(cafile=certifi.where())
except Exception:
    _SSL = ssl.create_default_context()

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROCESSED = ROOT / "data" / "processed" / "s_aureus_primary_collection_year_2015_2025_feasibility_rows.csv.gz"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
REPORTS = ROOT / "reports"
FIGS.mkdir(exist_ok=True)

INCOME_ORDER = ["High income", "Upper middle income", "Lower middle income", "Low income"]
NONHUMAN = ["livestock", "companion_animal", "food", "environment"]
ARMS = ["livestock", "companion_animal", "food", "environment"]

# ---- house style (identical to 06_figures.py) -----------------------------
plt.rcParams.update({
    "font.family": "Arial", "font.size": 9, "axes.titlesize": 10,
    "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False,
    "axes.linewidth": 0.8, "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
})
OI = {"blue": "#0072B2", "sky": "#56B4E9", "green": "#009E73", "orange": "#E69F00",
      "vermillion": "#D55E00", "purple": "#CC79A7", "yellow": "#F0E442", "grey": "#999999"}
ARM_COL = {"livestock": OI["vermillion"], "companion_animal": OI["green"],
           "food": OI["orange"], "environment": OI["sky"]}
ARM_LABEL = {"livestock": "Livestock", "companion_animal": "Companion", "food": "Food",
             "environment": "Environment"}


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(FIGS / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote figures/{name}.png / .pdf")


# ---- World Bank (same pattern as 05_inequity.py) --------------------------
def wb_fetch(url, timeout=90):
    with urllib.request.urlopen(url, timeout=timeout, context=_SSL) as r:
        return json.loads(r.read())


def wb_countries():
    # Offline-reproducible: prefer the frozen cache (00_freeze_raw_sources.py).
    cache = RAW / "worldbank_countries.csv"
    if cache.exists():
        c = pd.read_csv(cache, dtype=str)
        c = c.rename(columns={"region": "region_wb", "income_group": "income_group_wb"})
        return c[["iso3", "region_wb", "income_group_wb"]]
    data = wb_fetch("https://api.worldbank.org/v2/country/all?format=json&per_page=400")
    rows = []
    for r in data[1]:
        if r.get("region", {}).get("id") == "NA":
            continue
        rows.append({"iso3": r.get("id"),
                     "region_wb": r.get("region", {}).get("value"),
                     "income_group_wb": r.get("incomeLevel", {}).get("value")})
    return pd.DataFrame(rows)


def wb_population(start=2015, end=2023):
    # Offline-reproducible: prefer the frozen cache (00_freeze_raw_sources.py).
    cache = RAW / "worldbank_indicators.csv"
    if cache.exists():
        c = pd.read_csv(cache)
        c = c[(c["indicator"] == "SP.POP.TOTL") & (c["year"] >= start) & (c["year"] <= end)]
        c = c.sort_values(["iso3", "year"]).groupby("iso3", as_index=False).tail(1)
        out = c.rename(columns={"value": "population"})[["iso3", "population"]]
        out["population"] = out["population"].astype("int64")
        return out
    url = ("https://api.worldbank.org/v2/country/all/indicator/SP.POP.TOTL"
           f"?format=json&per_page=20000&date={start}:{end}")
    data = wb_fetch(url)
    rows = [{"iso3": r["countryiso3code"], "year": int(r["date"]), "population": r["value"]}
            for r in data[1] if r.get("value") is not None]
    df = pd.DataFrame(rows).sort_values(["iso3", "year"]).groupby("iso3", as_index=False).tail(1)
    return df[["iso3", "population"]]


def gini(x):
    """Gini coefficient of a non-negative array (0 = equal, 1 = maximally concentrated)."""
    x = np.sort(np.asarray(x, dtype=float))
    n = x.size
    if n == 0 or x.sum() == 0:
        return np.nan
    return float((2.0 * np.sum(np.arange(1, n + 1) * x) / (n * x.sum())) - (n + 1) / n)


def lorenz_points(counts):
    """Return (cum_frac_units, cum_frac_mass) including the origin, units sorted ascending."""
    x = np.sort(np.asarray(counts, dtype=float))
    cum = np.cumsum(x)
    cum_mass = np.concatenate([[0.0], cum / cum[-1]])
    cum_units = np.concatenate([[0.0], np.arange(1, x.size + 1) / x.size])
    return cum_units, cum_mass


def parse_tokens(s):
    if pd.isna(s):
        return set()
    return {p.strip().split("=", 1)[0].strip().lower()
            for p in str(s).replace('"', "").split(",") if p.strip()}


# ===========================================================================
def main():
    if not PROCESSED.exists():
        sys.exit(f"ERROR: {PROCESSED} not found; run 01 first.")
    df = pd.read_csv(PROCESSED, low_memory=False)
    df["collection_year"] = pd.to_numeric(df["collection_year"], errors="coerce")
    df["target_creation_year"] = pd.to_numeric(df["target_creation_year"], errors="coerce")
    for m in ["number_amr_genes", "number_virulence_genes"]:
        df[m] = pd.to_numeric(df[m], errors="coerce")
    n_total = len(df)

    print("Fetching World Bank region + population (certifi SSL)...")
    wbc = wb_countries()
    pop = wb_population()
    wb = wbc.merge(pop, on="iso3", how="left")
    wb["region_wb"] = wb["region_wb"].str.strip()  # WB emits trailing spaces on some region names
    region_by_iso = wb.set_index("iso3")["region_wb"].to_dict()
    pop_by_region = wb.groupby("region_wb")["population"].sum().to_dict()
    countries_by_region = wb.groupby("region_wb")["iso3"].nunique().to_dict()

    # ---- Table 38: source-classification QA --------------------------------
    qa_rows = []
    for src in ["human_or_clinical", "livestock", "companion_animal", "food", "environment", "unknown"]:
        sub = df[df["source_category"] == src]
        n = len(sub)
        has_host = sub["host"].notna()
        has_src = sub["isolation_source"].notna()
        has_type = sub["source_type"].notna()
        has_any = has_host | has_src | has_type
        qa_rows.append({
            "source_category": src, "n": n,
            "pct_of_total": round(100 * n / n_total, 2),
            "pct_with_host": round(100 * has_host.mean(), 1) if n else np.nan,
            "pct_with_isolation_source": round(100 * has_src.mean(), 1) if n else np.nan,
            "pct_with_source_type": round(100 * has_type.mean(), 1) if n else np.nan,
            "pct_with_any_informative_field": round(100 * has_any.mean(), 1) if n else np.nan,
            "n_with_no_informative_field": int((~has_any).sum()),
        })
    qa = pd.DataFrame(qa_rows)
    qa.to_csv(TABLES / "38_source_classification_qa.csv", index=False)

    # exemplar raw values per category (validation that classification is faithful)
    ex_rows = []
    for src in ["human_or_clinical", "livestock", "companion_animal", "food", "environment"]:
        sub = df[df["source_category"] == src]
        for field in ["host", "isolation_source"]:
            vc = sub[field].dropna().astype(str).str.strip().value_counts().head(6)
            for val, cnt in vc.items():
                ex_rows.append({"source_category": src, "field": field,
                                "raw_value": val, "records": int(cnt)})
    pd.DataFrame(ex_rows).to_csv(TABLES / "38b_source_classification_exemplars.csv", index=False)

    # ---- Table 39: geographic concentration + per-arm dominance -------------
    cc = df.groupby("country_parsed").size().sort_values(ascending=False)
    total = int(cc.sum())
    conc_rows = [{
        "stratum": "all_sources", "n_genomes": total, "n_countries": int(cc.size),
        "gini_across_countries": round(gini(cc.values), 3),
        "top1_country": cc.index[0], "top1_share_pct": round(100 * cc.iloc[0] / total, 1),
        "top3_share_pct": round(100 * cc.head(3).sum() / total, 1),
        "top5_share_pct": round(100 * cc.head(5).sum() / total, 1),
        "top10_share_pct": round(100 * cc.head(10).sum() / total, 1),
    }]
    # MRSA-marker-positive subset
    mr = df[df["mrsa_marker_positive"].astype(str).str.lower().isin({"true", "1", "1.0", "yes"})]
    ccm = mr.groupby("country_parsed").size().sort_values(ascending=False)
    tm = int(ccm.sum())
    conc_rows.append({
        "stratum": "mrsa_marker_positive", "n_genomes": tm, "n_countries": int(ccm.size),
        "gini_across_countries": round(gini(ccm.values), 3),
        "top1_country": ccm.index[0], "top1_share_pct": round(100 * ccm.iloc[0] / tm, 1),
        "top3_share_pct": round(100 * ccm.head(3).sum() / tm, 1),
        "top5_share_pct": round(100 * ccm.head(5).sum() / tm, 1),
        "top10_share_pct": round(100 * ccm.head(10).sum() / tm, 1),
    })
    # per-One-Health-arm single-country dominance
    for arm in ARMS:
        sub = df[df["source_category"] == arm]
        vc = sub["country_parsed"].value_counts()
        n = len(sub)
        conc_rows.append({
            "stratum": arm, "n_genomes": n, "n_countries": int(vc.size),
            "gini_across_countries": round(gini(vc.values), 3),
            "top1_country": vc.index[0], "top1_share_pct": round(100 * vc.iloc[0] / n, 1),
            "top3_share_pct": round(100 * vc.head(3).sum() / n, 1),
            "top5_share_pct": round(100 * vc.head(5).sum() / n, 1),
            "top10_share_pct": round(100 * vc.head(10).sum() / n, 1),
        })
    conc = pd.DataFrame(conc_rows)
    conc.to_csv(TABLES / "39_geographic_concentration.csv", index=False)

    # ---- Table 40: regional visibility -------------------------------------
    df["region_wb"] = df["iso3"].map(region_by_iso)
    reg_rows = []
    present_iso = set(df["iso3"].dropna())
    for region in sorted(pop_by_region, key=lambda r: -df[df["region_wb"] == r].shape[0]):
        sub = df[df["region_wb"] == region]
        n = len(sub)
        nh = int(sub["source_category"].isin(NONHUMAN).sum())
        popr = pop_by_region.get(region, np.nan)
        have = wb[(wb["region_wb"] == region) & (wb["iso3"].isin(present_iso))]["iso3"].nunique()
        tot_c = countries_by_region.get(region, 0)
        reg_rows.append({
            "region": region, "records": n,
            "pct_of_total": round(100 * n / n_total, 2),
            "population_million": round(popr / 1e6, 1) if pd.notna(popr) else np.nan,
            "records_per_million": round(n / popr * 1e6, 3) if pd.notna(popr) and popr else np.nan,
            "non_human_records": nh,
            "non_human_per_million": round(nh / popr * 1e6, 4) if pd.notna(popr) and popr else np.nan,
            "countries_with_genome": int(have), "countries_total": int(tot_c),
            "coverage_pct": round(100 * have / tot_c, 1) if tot_c else np.nan,
        })
    reg = pd.DataFrame(reg_rows)
    reg.to_csv(TABLES / "40_regional_visibility.csv", index=False)

    # ---- Table 41: window & date-definition sensitivity --------------------
    w = pd.read_csv(TABLES / "18_window_comparison_2020_2025_vs_2015_2025.csv")
    keep = ["window", "year_definition", "start", "end", "records", "mrsa_pct",
            "nonhuman", "nonhuman_pct", "livestock", "companion_animal", "food",
            "environment", "unique_countries"]
    w[keep].to_csv(TABLES / "41_window_sensitivity.csv", index=False)

    # ---- Table 42: master gene table (prevalence+CI by source + omnibus) ----
    prev = pd.read_csv(TABLES / "26_gene_prevalence_by_source.csv")
    omni = pd.read_csv(TABLES / "27_gene_omnibus_tests.csv")[
        ["gene", "chi2", "dof", "p_value", "q_value_bh", "min_expected_count", "significant_q<0.05"]]
    srcs = ["human_or_clinical", "livestock", "companion_animal", "food", "environment"]
    master = prev[["gene"]].copy()
    for s in srcs:
        master[f"{s}"] = (prev[f"{s}_pct"].map(lambda v: f"{v:.1f}")
                          + " (" + prev[f"{s}_ci"].astype(str) + ")")
        master[f"{s}_pos_n"] = prev[f"{s}_pos"].astype(str) + "/" + prev[f"{s}_n"].astype(str)
    master = master.merge(omni, on="gene", how="left").sort_values("q_value_bh", na_position="last")
    master.to_csv(TABLES / "42_gene_master_by_source.csv", index=False)

    # =======================================================================
    # FIGURE S1 — geographic concentration
    # =======================================================================
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 3.4))
    u, m = lorenz_points(cc.values)
    a.plot([0, 1], [0, 1], color=OI["grey"], ls="--", lw=1.0)
    a.plot(u, m, color=OI["blue"], lw=2.0)
    a.fill_between(u, m, u, color=OI["blue"], alpha=0.12)
    a.set_xlabel("Cumulative share of countries")
    a.set_ylabel("Cumulative share of deposited genomes")
    a.set_title("A  Concentration across countries")
    a.set_xlim(0, 1); a.set_ylim(0, 1)
    a.text(0.05, 0.86, f"Gini = {gini(cc.values):.3f}\nTop 10 countries = "
           f"{100*cc.head(10).sum()/total:.0f}%\n{cc.size} countries with ≥1 genome",
           transform=a.transAxes, fontsize=8, va="top")

    dom = conc[conc["stratum"].isin(ARMS)].set_index("stratum").reindex(ARMS)
    y = np.arange(len(ARMS))[::-1]
    b.barh(y, dom["top1_share_pct"], color=[ARM_COL[s] for s in ARMS])
    b.set_yticks(y); b.set_yticklabels([ARM_LABEL[s] for s in ARMS])
    b.set_xlabel("Share from the single top country (%)")
    b.set_title("B  Single-country dominance, One Health arms")
    b.set_xlim(0, 100)
    for yi, s in zip(y, ARMS):
        r = dom.loc[s]
        b.text(r["top1_share_pct"] + 1.5, yi,
               f"{r['top1_country']} {r['top1_share_pct']:.0f}%\n(n={int(r['n_genomes']):,}; {int(r['n_countries'])} countries)",
               va="center", fontsize=7)
    fig.suptitle("Geographic concentration of public S. aureus/MRSA genomic visibility",
                 fontsize=10, fontweight="bold", y=1.04)
    fig.tight_layout()
    save(fig, "figS1_geographic_concentration")

    # =======================================================================
    # FIGURE S2 — temporal deposition dynamics (within the primary cohort)
    # =======================================================================
    years = list(range(2015, 2026))
    coll = df["collection_year"].value_counts().reindex(years).fillna(0)
    dep = df["target_creation_year"].value_counts().reindex(years).fillna(0)
    nh_share = (df[df["source_category"].isin(NONHUMAN)].groupby("collection_year").size()
                .reindex(years).fillna(0) / coll.replace(0, np.nan) * 100)

    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.3))
    a.plot(years, coll.values, "-o", color=OI["blue"], lw=1.8, ms=4, label="by collection year")
    a.plot(years, dep.values, "-s", color=OI["orange"], lw=1.8, ms=4, label="by deposition year")
    a.set_xlabel("Year"); a.set_ylabel("Deposited genomes (n)")
    a.set_title("A  Collection vs deposition year")
    a.legend(frameon=False, fontsize=7.5, loc="upper right")
    a.set_xticks(years[::2])

    b.plot(years, nh_share.values, "-o", color=OI["green"], lw=1.8, ms=4)
    b.set_xlabel("Collection year")
    b.set_ylabel("Non-human share of genomes (%)")
    b.set_title("B  One Health share over time")
    b.set_xticks(years[::2])
    b.set_ylim(0, max(6, float(np.nanmax(nh_share.values)) * 1.25))
    fig.suptitle("Temporal deposition dynamics, primary cohort (collection-year 2015–2025)",
                 fontsize=10, fontweight="bold", y=1.04)
    fig.tight_layout()
    save(fig, "figS2_temporal_dynamics")

    # =======================================================================
    # FIGURE S3 — regional visibility
    # =======================================================================
    rr = reg.sort_values("records_per_million")
    y = np.arange(len(rr))
    short = {"Middle East, North Africa, Afghanistan & Pakistan": "MENA + Afg/Pak",
             "Latin America & Caribbean": "Latin America & Carib.",
             "Europe & Central Asia": "Europe & Central Asia"}
    labels = [short.get(r, r) for r in rr["region"]]
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    bars = ax.barh(y, rr["records_per_million"], color=OI["blue"])
    # highlight the two most under-represented regions
    for yi, region in zip(y, rr["region"]):
        if region in ("South Asia", "Sub-Saharan Africa"):
            bars[yi].set_color(OI["vermillion"])
    ax.set_xscale("log")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Genomes per million population (log scale)")
    ax.set_title("Public genomic visibility by World Bank region", fontsize=9.5)
    for yi, r in zip(y, rr.itertuples()):
        ax.text(r.records_per_million * 1.1, yi,
                f"{r.records_per_million:.2f}  (n={int(r.records):,}; {int(r.countries_with_genome)}/{int(r.countries_total)} countries)",
                va="center", fontsize=7)
    ax.set_xlim(rr["records_per_million"].min() * 0.5, rr["records_per_million"].max() * 12)
    fig.tight_layout()
    save(fig, "figS3_regional_visibility")

    # =======================================================================
    # narrative report
    # =======================================================================
    L = []
    L.append("# Technical Appendix — supporting analyses (script 07)\n")
    L.append(f"Primary cohort: collection-year 2015–2025, n={n_total:,}.\n")

    L.append("## Source-classification QA (table 38)\n")
    L.append("| source | n | % total | % host | % isol.source | % any field | n no-field |")
    L.append("|---|---|---|---|---|---|---|")
    for _, r in qa.iterrows():
        L.append(f"| {r['source_category']} | {r['n']:,} | {r['pct_of_total']} | {r['pct_with_host']} | "
                 f"{r['pct_with_isolation_source']} | {r['pct_with_any_informative_field']} | {r['n_with_no_informative_field']:,} |")
    unk = qa[qa.source_category == "unknown"].iloc[0]
    L.append(f"\n**Of {int(unk['n']):,} 'unknown' genomes, {int(unk['n_with_no_informative_field']):,} "
             f"({100*unk['n_with_no_informative_field']/unk['n']:.0f}%) carried no host, isolation-source or "
             "source-type field at all — 'unknown' reflects missing metadata, not classifier failure.**\n")

    L.append("## Geographic concentration (table 39)\n")
    L.append("| stratum | genomes | countries | Gini | top country | top-1 % | top-3 % | top-10 % |")
    L.append("|---|---|---|---|---|---|---|---|")
    for _, r in conc.iterrows():
        L.append(f"| {r['stratum']} | {r['n_genomes']:,} | {r['n_countries']} | {r['gini_across_countries']} | "
                 f"{r['top1_country']} | {r['top1_share_pct']} | {r['top3_share_pct']} | {r['top10_share_pct']} |")

    L.append("\n## Regional visibility (table 40)\n")
    L.append("| region | records | per million | non-human/M | coverage |")
    L.append("|---|---|---|---|---|")
    for _, r in reg.iterrows():
        L.append(f"| {r['region']} | {r['records']:,} | {r['records_per_million']} | {r['non_human_per_million']} | "
                 f"{r['countries_with_genome']}/{r['countries_total']} ({r['coverage_pct']}%) |")

    L.append("\n## Window & date-definition sensitivity (table 41)\n")
    L.append("| window | records | MRSA % | non-human % | countries |")
    L.append("|---|---|---|---|---|")
    for _, r in w.iterrows():
        L.append(f"| {r['window']} | {r['records']:,} | {r['mrsa_pct']} | {r['nonhuman_pct']} | {r['unique_countries']} |")

    L.append("\n## Guardrails\n")
    L.append("- Concentration and single-country dominance are properties of the public repository "
             "(a convenience sample), not of MRSA biology; they quantify how narrow the evidence base is.\n"
             "- Regional per-million rates measure sequencing/deposition intensity, not disease burden.\n"
             "- The collection-year decline after ~2022 reflects submission lag (recent infections not yet "
             "deposited), not a real fall in incidence — hence the deposition-year comparison in Fig S2.\n")

    (REPORTS / "07_appendix.md").write_text("\n".join(L))
    print("Wrote tables 38-42, figures S1-S3 and reports/07_appendix.md")
    print(f"Gini across countries = {gini(cc.values):.3f}; top-10 share = {100*cc.head(10).sum()/total:.1f}%")


if __name__ == "__main__":
    main()
