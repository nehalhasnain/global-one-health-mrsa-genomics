#!/usr/bin/env python3
"""
Feasibility extraction for Global One Health MRSA public-data study.

Inputs:
- NCBI Pathogen Detection Staphylococcus aureus AMR metadata TSV (latest PDG)
- World Bank country metadata API

Outputs:
- processed feasibility CSV summaries under tables/
- markdown report under reports/

This is a metadata-first feasibility test; it does not download assemblies.
"""
from __future__ import annotations

import csv
import datetime as dt
import io
import json
import math
import os
import re
import sys
import urllib.request
import ssl
from collections import Counter, defaultdict
from pathlib import Path

import certifi
import pandas as pd

# macOS Python sometimes fails NCBI/World Bank HTTPS verification with a stale or
# intercepted certificate chain. Use certifi explicitly for all urllib/pandas HTTPS calls.
SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
ssl._create_default_https_context = lambda: SSL_CONTEXT

PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data" / "raw"
PROCESSED = PROJECT / "data" / "processed"
TABLES = PROJECT / "tables"
REPORTS = PROJECT / "reports"
LOGS = PROJECT / "logs"
for d in [RAW, PROCESSED, TABLES, REPORTS, LOGS]:
    d.mkdir(parents=True, exist_ok=True)

NCBI_BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Staphylococcus_aureus/"
WINDOW_START = 2015
WINDOW_END = 2025
RUN_DATE = dt.datetime.now().astimezone().isoformat(timespec="seconds")

USECOLS = [
    "#label", "asm_acc", "biosample_acc", "collection_date", "geo_loc_name", "host",
    "isolation_source", "source_type", "target_acc", "target_creation_date",
    "computed_types", "number_amr_genes", "AMR_genotypes", "AMR_genotypes_core",
    "number_virulence_genes", "virulence_genotypes", "amrfinder_version", "refgene_db_version",
]


def fetch_text(url: str, timeout: int = 60) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")


def latest_pdg_dir() -> str:
    html = fetch_text(NCBI_BASE)
    hrefs = re.findall(r'href="(PDG000000073\.\d+/)"', html)
    if not hrefs:
        raise RuntimeError("No PDG directories found at NCBI Staphylococcus_aureus path")
    return max(hrefs, key=lambda h: int(h.rstrip("/").split(".")[-1]))


def url_exists(url: str) -> bool:
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=30):
            return True
    except Exception:
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                r.read(1)
            return True
        except Exception:
            return False


def read_ncbi_metadata(url: str) -> pd.DataFrame:
    # pandas can stream URL directly; keep only feasibility columns.
    df = pd.read_csv(
        url,
        sep="\t",
        usecols=lambda c: c in USECOLS,
        dtype=str,
        na_values=["NULL", "", "NA", "N/A", "nan", "NaN"],
        keep_default_na=True,
        low_memory=False,
    )
    missing = [c for c in USECOLS if c not in df.columns]
    if missing:
        print(f"WARNING: missing expected columns: {missing}", file=sys.stderr)
    return df


def extract_year(x) -> float:
    if pd.isna(x):
        return math.nan
    s = str(x)
    m = re.search(r"(19\d{2}|20\d{2})", s)
    if not m:
        return math.nan
    y = int(m.group(1))
    if y < 1900 or y > 2035:
        return math.nan
    return y


def parse_country(raw) -> str | None:
    if pd.isna(raw):
        return None
    s = str(raw).strip()
    if not s or s.upper() == "NULL":
        return None
    # NCBI often uses "Country: state/city".
    country = s.split(":", 1)[0].strip()
    country = re.sub(r"\s+", " ", country)
    if not country:
        return None
    aliases = {
        "USA": "United States", "U.S.A.": "United States", "US": "United States",
        "United States of America": "United States",
        "UK": "United Kingdom", "England": "United Kingdom", "Scotland": "United Kingdom",
        "Wales": "United Kingdom", "Northern Ireland": "United Kingdom",
        "Republic of Korea": "Korea, Rep.", "South Korea": "Korea, Rep.",
        "Korea": "Korea, Rep.", "Viet Nam": "Vietnam", "Russia": "Russian Federation",
        "Iran": "Iran, Islamic Rep.", "Venezuela": "Venezuela, RB",
        "Czech Republic": "Czechia", "Democratic Republic of the Congo": "Congo, Dem. Rep.",
        "DRC": "Congo, Dem. Rep.", "Congo": "Congo, Rep.",
        "Egypt": "Egypt, Arab Rep.", "Syria": "Syrian Arab Republic",
        "Slovak Republic": "Slovakia", "Turkey": "Turkiye",
        "Türkiye": "Turkiye", "Laos": "Lao PDR", "Tanzania": "Tanzania",
        "Bolivia": "Bolivia", "Moldova": "Moldova", "Brunei": "Brunei Darussalam",
        "Cape Verde": "Cabo Verde", "Micronesia": "Micronesia, Fed. Sts.",
        "Gambia": "Gambia, The", "Bahamas": "Bahamas, The", "Yemen": "Yemen, Rep.",
    }
    return aliases.get(country, country)


def normalize_text(*vals) -> str:
    parts = []
    for v in vals:
        if pd.notna(v):
            parts.append(str(v).lower())
    return " | ".join(parts)


def token_set(s: str) -> set[str]:
    return set(re.findall(r"[a-z]+", s.lower()))


def classify_source(row) -> str:
    text = normalize_text(row.get("host"), row.get("isolation_source"), row.get("source_type"))
    toks = token_set(text)

    # Specific scientific/common names first. Avoid raw substring matching for cat/cattle.
    livestock_terms = {
        "bos taurus", "sus scrofa", "gallus gallus", "ovis aries", "capra hircus", "bubalus bubalis",
        "cattle", "cow", "bovine", "calf", "swine", "pig", "porcine", "pork",
        "chicken", "poultry", "broiler", "layer", "turkey", "duck", "goose",
        "sheep", "ovine", "goat", "caprine", "buffalo", "camel", "llama", "alpaca",
    }
    companion_terms = {
        "felis catus", "canis lupus familiaris", "dog", "dogs", "cat", "cats", "feline", "canine",
        "horse", "equine", "pet", "companion",
    }
    environment_terms = {
        "wastewater", "sewage", "effluent", "environment", "environmental", "water", "soil",
        "surface", "air", "dust", "sediment", "farm environment", "sludge",
    }
    food_terms = {
        "food", "meat", "retail", "milk", "dairy", "cheese", "beef", "pork", "chicken meat",
        "turkey meat", "carcass", "seafood", "fish", "egg", "eggs", "ready to eat",
    }
    human_terms = {
        "homo sapiens", "human", "patient", "clinical", "hospital", "blood", "wound", "nasal",
        "urine", "sputum", "respiratory", "abscess", "skin", "infection", "screening",
    }

    def has_phrase(phrases):
        return any(p in text for p in phrases if " " in p)

    # Token-sensitive single-word checks.
    if has_phrase(livestock_terms) or toks.intersection({t for t in livestock_terms if " " not in t}):
        return "livestock"
    if has_phrase(companion_terms) or toks.intersection({t for t in companion_terms if " " not in t}):
        return "companion_animal"
    if has_phrase(environment_terms) or toks.intersection({t for t in environment_terms if " " not in t}):
        return "environment"
    if has_phrase(food_terms) or toks.intersection({t for t in food_terms if " " not in t}):
        return "food"
    if has_phrase(human_terms) or toks.intersection({t for t in human_terms if " " not in t}):
        return "human_or_clinical"
    return "unknown"


def parse_gene_set(s) -> set[str]:
    if pd.isna(s):
        return set()
    out = set()
    for part in str(s).replace('"', '').split(','):
        g = part.strip()
        if not g:
            continue
        g = g.split('=', 1)[0].strip()
        out.add(g)
    return out


def has_gene(row, genes: set[str]) -> bool:
    gs = set()
    gs.update(parse_gene_set(row.get("AMR_genotypes")))
    gs.update(parse_gene_set(row.get("AMR_genotypes_core")))
    return any(g.lower() in {x.lower() for x in gs} for g in genes)


def load_worldbank_countries() -> pd.DataFrame:
    # Offline-reproducible: prefer the frozen cache (00_freeze_raw_sources.py).
    cache = RAW / "worldbank_countries.csv"
    if cache.exists():
        c = pd.read_csv(cache, dtype=str)
        cols = ["country_name_wb", "iso2", "iso3", "region", "income_id", "income_group"]
        return c[[col for col in cols if col in c.columns]]
    url = "https://api.worldbank.org/v2/country/all?format=json&per_page=400"
    data = json.loads(fetch_text(url, timeout=60))
    rows = []
    for r in data[1]:
        if r.get("region", {}).get("id") == "NA":
            continue
        rows.append({
            "country_name_wb": r.get("name"),
            "iso2": r.get("iso2Code"),
            "iso3": r.get("id"),
            "region": r.get("region", {}).get("value"),
            "income_id": r.get("incomeLevel", {}).get("id"),
            "income_group": r.get("incomeLevel", {}).get("value"),
        })
    return pd.DataFrame(rows)


def load_worldbank_indicator(indicator: str, value_name: str, years=(2025, 2024, 2023, 2022, 2021)) -> pd.DataFrame:
    # Pull a recent range; use latest non-null value per country.
    start, end = min(years), max(years)
    # Offline-reproducible: prefer the frozen cache (00_freeze_raw_sources.py).
    cache = RAW / "worldbank_indicators.csv"
    if cache.exists():
        c = pd.read_csv(cache)
        c = c[(c["indicator"] == indicator) & (c["year"] >= start) & (c["year"] <= end)]
        if c.empty:
            return pd.DataFrame(columns=["iso3", value_name, f"{value_name}_year"])
        c = c.sort_values(["iso3", "year"]).groupby("iso3", as_index=False).tail(1)
        out = c.rename(columns={"value": value_name, "year": f"{value_name}_year"})[["iso3", value_name, f"{value_name}_year"]]
        if value_name.startswith("population"):
            out[value_name] = out[value_name].astype("int64")
        return out
    url = f"https://api.worldbank.org/v2/country/all/indicator/{indicator}?format=json&per_page=20000&date={start}:{end}"
    try:
        data = json.loads(fetch_text(url, timeout=90))
        rows = []
        for r in data[1]:
            if r.get("value") is None:
                continue
            rows.append({"iso3": r.get("countryiso3code"), "year": int(r.get("date")), value_name: r.get("value")})
        df = pd.DataFrame(rows)
        if df.empty:
            return pd.DataFrame(columns=["iso3", value_name, f"{value_name}_year"])
        df = df.sort_values(["iso3", "year"]).groupby("iso3", as_index=False).tail(1)
        return df.rename(columns={"year": f"{value_name}_year"})[["iso3", value_name, f"{value_name}_year"]]
    except Exception as e:
        print(f"WARNING: failed World Bank indicator {indicator}: {e}", file=sys.stderr)
        return pd.DataFrame(columns=["iso3", value_name, f"{value_name}_year"])


def write_csv(df: pd.DataFrame, name: str) -> Path:
    p = TABLES / name
    df.to_csv(p, index=False)
    return p


def pct(n, d):
    return round(100 * n / d, 2) if d else 0.0


def find_frozen_metadata():
    """Return (Path, pdg) of a frozen metadata snapshot in data/raw, or (None, None).

    Prefers the snapshot pinned in latest_ncbi_metadata_url.txt so the freeze
    reproduces the published results rather than whatever is newest on-disk.
    """
    files = sorted(RAW.glob("ncbi_s_aureus_amr_metadata.*.tsv.gz"))
    if not files:
        return None, None
    pin = RAW / "latest_ncbi_metadata_url.txt"
    if pin.exists():
        pinned_pdg = pin.read_text().split("/AMR/", 1)[0].rstrip("/").split("/")[-1]
        for f in files:
            if pinned_pdg in f.name:
                return f, pinned_pdg
    f = files[-1]
    return f, f.name[len("ncbi_s_aureus_amr_metadata."):-len(".tsv.gz")]


def main():
    # Offline-reproducible: if a frozen metadata snapshot exists (00_freeze_raw_sources.py),
    # pin to it instead of the newest snapshot on the NCBI server.
    frozen_meta, frozen_pdg = find_frozen_metadata()
    if frozen_meta is not None:
        pdg_clean = frozen_pdg
        pdg = pdg_clean + "/"
        metadata_url = f"{NCBI_BASE}{pdg}AMR/{pdg_clean}.amr.metadata.tsv"
        clusters_url = f"{NCBI_BASE}{pdg}Clusters/{pdg_clean}.reference_target.all_isolates.tsv"
        print(f"Frozen PDG (offline): {pdg_clean}")
        print(f"Reading frozen metadata: {frozen_meta.name}")
        print(f"Run date: {RUN_DATE}")
        df = read_ncbi_metadata(str(frozen_meta))
    else:
        pdg = latest_pdg_dir()
        pdg_clean = pdg.rstrip("/")
        metadata_url = f"{NCBI_BASE}{pdg}AMR/{pdg_clean}.amr.metadata.tsv"
        clusters_url = f"{NCBI_BASE}{pdg}Clusters/{pdg_clean}.reference_target.all_isolates.tsv"
        print(f"Latest PDG: {pdg_clean}")
        print(f"Metadata URL: {metadata_url}")
        print(f"Run date: {RUN_DATE}")
        df = read_ncbi_metadata(metadata_url)
    df.columns = [c.lstrip("#") if c == "#label" else c for c in df.columns]
    if "label" not in df.columns and "#label" in df.columns:
        df = df.rename(columns={"#label": "label"})

    # Save columns manifest and a small raw sample for provenance.
    (RAW / "latest_ncbi_metadata_url.txt").write_text(metadata_url + "\n", encoding="utf-8")
    df.head(1000).to_csv(RAW / "ncbi_s_aureus_metadata_sample_1000.csv", index=False)

    df["collection_year"] = df.get("collection_date", pd.Series(index=df.index, dtype=str)).map(extract_year)
    df["target_creation_year"] = df.get("target_creation_date", pd.Series(index=df.index, dtype=str)).map(extract_year)
    df["analysis_year"] = df["collection_year"].fillna(df["target_creation_year"])
    df["country_parsed"] = df.get("geo_loc_name", pd.Series(index=df.index, dtype=str)).map(parse_country)
    df["source_category"] = df.apply(classify_source, axis=1)
    df["mecA_positive"] = df.apply(lambda r: has_gene(r, {"mecA"}), axis=1)
    df["mecC_positive"] = df.apply(lambda r: has_gene(r, {"mecC"}), axis=1)
    df["mrsa_marker_positive"] = df["mecA_positive"] | df["mecC_positive"]
    df["has_country"] = df["country_parsed"].notna()
    df["has_host"] = df.get("host", pd.Series(index=df.index, dtype=str)).notna()
    df["has_source"] = df.get("isolation_source", pd.Series(index=df.index, dtype=str)).notna() | df.get("source_type", pd.Series(index=df.index, dtype=str)).notna()
    df["has_collection_year"] = df["collection_year"].notna()
    df["has_target_year"] = df["target_creation_year"].notna()

    primary = df[(df["collection_year"] >= WINDOW_START) & (df["collection_year"] <= WINDOW_END)].copy()
    sensitivity = df[(df["analysis_year"] >= WINDOW_START) & (df["analysis_year"] <= WINDOW_END)].copy()

    # Join World Bank country metadata.
    wb = load_worldbank_countries()
    pop = load_worldbank_indicator("SP.POP.TOTL", "population_latest")
    wb = wb.merge(pop, on="iso3", how="left")
    wb_name_map = wb.copy()
    primary = primary.merge(wb_name_map, left_on="country_parsed", right_on="country_name_wb", how="left")
    sensitivity = sensitivity.merge(wb_name_map, left_on="country_parsed", right_on="country_name_wb", how="left")

    # Overall counts.
    overall = pd.DataFrame([
        {"metric": "latest_pdg", "value": pdg_clean},
        {"metric": "run_date", "value": RUN_DATE},
        {"metric": "total_ncbi_rows_all_dates", "value": len(df)},
        {"metric": "primary_rows_collection_year_2015_2025", "value": len(primary)},
        {"metric": "sensitivity_rows_collection_or_target_year_2015_2025", "value": len(sensitivity)},
        {"metric": "primary_mrsa_marker_positive", "value": int(primary["mrsa_marker_positive"].sum())},
        {"metric": "primary_mecA_positive", "value": int(primary["mecA_positive"].sum())},
        {"metric": "primary_mecC_positive", "value": int(primary["mecC_positive"].sum())},
        {"metric": "primary_unique_countries_parsed", "value": int(primary["country_parsed"].nunique())},
        {"metric": "primary_unique_countries_matched_worldbank", "value": int(primary["iso3"].nunique())},
        {"metric": "clusters_all_isolates_url_exists", "value": str(url_exists(clusters_url))},
        {"metric": "metadata_url", "value": metadata_url},
        {"metric": "clusters_url", "value": clusters_url},
    ])
    write_csv(overall, "01_overall_counts.csv")

    # Metadata completeness.
    completeness_rows = []
    for name, subset in [("all_dates", df), ("primary_collection_2015_2025", primary), ("sensitivity_analysis_year_2015_2025", sensitivity)]:
        d = len(subset)
        for col, label in [
            ("has_country", "country"), ("has_host", "host"), ("has_source", "source_or_source_type"),
            ("has_collection_year", "collection_year"), ("has_target_year", "target_creation_year"),
        ]:
            n = int(subset[col].sum()) if col in subset.columns else 0
            completeness_rows.append({"dataset": name, "field": label, "non_missing_n": n, "total_n": d, "non_missing_pct": pct(n, d)})
    write_csv(pd.DataFrame(completeness_rows), "02_metadata_completeness.csv")

    # Source summary.
    src = primary.groupby("source_category", dropna=False).agg(
        total_records=("label", "size"),
        mrsa_marker_positive=("mrsa_marker_positive", "sum"),
        mecA_positive=("mecA_positive", "sum"),
        mecC_positive=("mecC_positive", "sum"),
        unique_countries=("country_parsed", "nunique"),
    ).reset_index()
    src["mrsa_marker_pct"] = (src["mrsa_marker_positive"] / src["total_records"] * 100).round(2)
    src = src.sort_values("total_records", ascending=False)
    write_csv(src, "03_source_mrsa_summary.csv")

    # Year x source.
    year_src = primary.groupby(["collection_year", "source_category"], dropna=False).agg(
        records=("label", "size"), mrsa_marker_positive=("mrsa_marker_positive", "sum")
    ).reset_index().sort_values(["collection_year", "source_category"])
    write_csv(year_src, "04_year_source_counts.csv")

    # Income summary.
    inc = primary.groupby("income_group", dropna=False).agg(
        records=("label", "size"),
        mrsa_marker_positive=("mrsa_marker_positive", "sum"),
        nonhuman_records=("source_category", lambda s: int(s.isin(["livestock", "companion_animal", "food", "environment"]).sum())),
        unique_countries=("country_parsed", "nunique"),
    ).reset_index().sort_values("records", ascending=False)
    inc["mrsa_marker_pct"] = (inc["mrsa_marker_positive"] / inc["records"] * 100).round(2)
    write_csv(inc, "05_income_group_summary.csv")

    # Country summary with population-normalized rate for matched countries.
    country = primary.groupby(["country_parsed", "iso3", "income_group", "region"], dropna=False).agg(
        records=("label", "size"),
        mrsa_marker_positive=("mrsa_marker_positive", "sum"),
        nonhuman_records=("source_category", lambda s: int(s.isin(["livestock", "companion_animal", "food", "environment"]).sum())),
        human_or_clinical_records=("source_category", lambda s: int((s == "human_or_clinical").sum())),
        unknown_source_records=("source_category", lambda s: int((s == "unknown").sum())),
    ).reset_index()
    pop_map = wb[["iso3", "population_latest", "population_latest_year"]].drop_duplicates("iso3")
    country = country.drop(columns=[c for c in ["population_latest", "population_latest_year"] if c in country.columns], errors="ignore").merge(pop_map, on="iso3", how="left")
    country["records_per_million_pop"] = (country["records"] / country["population_latest"] * 1_000_000).round(4)
    country["mrsa_per_million_pop"] = (country["mrsa_marker_positive"] / country["population_latest"] * 1_000_000).round(4)
    country = country.sort_values("records", ascending=False)
    write_csv(country.head(100), "06_top100_country_visibility.csv")
    write_csv(country.sort_values("nonhuman_records", ascending=False).head(100), "07_top100_country_nonhuman_visibility.csv")

    # Top raw metadata values for classifier QA.
    for col, name in [("host", "08_top_hosts.csv"), ("isolation_source", "09_top_isolation_sources.csv"), ("source_type", "10_top_source_types.csv"), ("geo_loc_name", "11_top_geo_loc_name.csv")]:
        if col in primary.columns:
            top = primary[col].fillna("<missing>").value_counts().head(100).rename_axis(col).reset_index(name="records")
            write_csv(top, name)

    # Gene-by-source feasibility for top AMR genes.
    gene_counts = []
    for _, row in primary.iterrows():
        genes = parse_gene_set(row.get("AMR_genotypes")) | parse_gene_set(row.get("AMR_genotypes_core"))
        for g in genes:
            gene_counts.append((row["source_category"], g))
    if gene_counts:
        gdf = pd.DataFrame(gene_counts, columns=["source_category", "gene"])
        gtab = gdf.groupby(["source_category", "gene"]).size().reset_index(name="records")
        gtab = gtab.sort_values(["source_category", "records"], ascending=[True, False])
        top_gene_table = gtab.groupby("source_category").head(30).reset_index(drop=True)
        write_csv(top_gene_table, "12_top30_amr_genes_by_source.csv")

    # Save a compressed feasibility-level row table for next steps.
    keep = [
        "label", "asm_acc", "biosample_acc", "collection_date", "collection_year", "target_creation_date",
        "target_creation_year", "geo_loc_name", "country_parsed", "iso3", "income_group", "region",
        "host", "isolation_source", "source_type", "source_category", "target_acc", "computed_types",
        "number_amr_genes", "AMR_genotypes", "AMR_genotypes_core", "number_virulence_genes",
        "virulence_genotypes", "mecA_positive", "mecC_positive", "mrsa_marker_positive",
    ]
    keep = [c for c in keep if c in primary.columns]
    primary[keep].to_csv(PROCESSED / "s_aureus_primary_collection_year_2015_2025_feasibility_rows.csv.gz", index=False, compression="gzip")

    # Report.
    total = len(primary)
    mrsa = int(primary["mrsa_marker_positive"].sum())
    nonhuman = int(primary["source_category"].isin(["livestock", "companion_animal", "food", "environment"]).sum())
    unknown = int((primary["source_category"] == "unknown").sum())
    matched = int(primary["iso3"].notna().sum())
    report = []
    report.append("# Feasibility Test 1 — NCBI Pathogen Detection Staphylococcus aureus / MRSA, 2015–2025\n")
    report.append(f"Run date: {RUN_DATE}\n")
    report.append(f"Latest NCBI PDG snapshot: `{pdg_clean}`\n")
    report.append(f"Metadata URL: {metadata_url}\n")
    report.append(f"Cluster all-isolates URL available: `{url_exists(clusters_url)}`\n")
    report.append("## Primary window definition\n")
    report.append(f"Primary filter: collection year {WINDOW_START}–{WINDOW_END}. Records with missing collection year are not counted in the primary window; they are tracked in the sensitivity count.\n")
    report.append("## Headline feasibility counts\n")
    report.append(f"- Total NCBI rows all dates: {len(df):,}\n")
    report.append(f"- Primary rows with collection year 2015–2025: {total:,}\n")
    report.append(f"- Sensitivity rows using collection year or target-creation year 2015–2025: {len(sensitivity):,}\n")
    report.append(f"- MRSA-marker-positive rows in primary window (mecA or mecC): {mrsa:,} ({pct(mrsa,total)}%)\n")
    report.append(f"- Non-human One Health rows in primary window (livestock + companion + food + environment): {nonhuman:,} ({pct(nonhuman,total)}%)\n")
    report.append(f"- Unknown source rows in primary window: {unknown:,} ({pct(unknown,total)}%)\n")
    report.append(f"- Rows matched to World Bank country metadata: {matched:,} ({pct(matched,total)}%)\n")
    report.append("## Source-category summary\n\n")
    report.append(src.to_markdown(index=False))
    report.append("\n\n## Income-group summary\n\n")
    report.append(inc.to_markdown(index=False))
    report.append("\n\n## Top 20 countries by public S. aureus/MRSA sequence visibility\n\n")
    report.append(country.head(20).to_markdown(index=False))
    report.append("\n\n## Preliminary go/no-go interpretation\n")
    if nonhuman >= 1000 and mrsa >= 1000:
        report.append("- GO for the full One Health interface paper: the primary window has enough non-human records and enough mecA/mecC-positive records for source-stratified analysis.\n")
    elif mrsa >= 1000:
        report.append("- GO for global MRSA genomic visibility and inequity paper; One Health interface can remain a secondary analysis unless non-human classification improves.\n")
    else:
        report.append("- CAUTION: MRSA-marker-positive denominator may be low under current definitions; inspect AMR genotype parsing and consider widening the window or using additional GenBank/BioSample metadata.\n")
    report.append("- Interpret all counts as public genomic visibility, not incidence/prevalence/burden.\n")
    report.append("- Next test: inspect classifier QA tables and then join NCBI SNP cluster all-isolates metadata to find mixed human–animal–environment clusters.\n")
    report.append("\n## Output files\n")
    for p in sorted(TABLES.glob("*.csv")):
        report.append(f"- `{p.relative_to(PROJECT)}`\n")
    report.append(f"- `data/processed/{(PROCESSED / 's_aureus_primary_collection_year_2015_2025_feasibility_rows.csv.gz').name}`\n")
    report_path = REPORTS / "feasibility_test_1_ncbi_s_aureus_mrsa_2015_2025.md"
    report_path.write_text("".join(report), encoding="utf-8")
    print(f"Report written: {report_path}")
    print(f"Primary rows: {total:,}; MRSA-marker-positive: {mrsa:,}; non-human: {nonhuman:,}")


if __name__ == "__main__":
    main()
