#!/usr/bin/env python3
"""Compare contemporary 2020–2025 vs 11-year 2015–2025 windows for the MRSA public-data study."""
from __future__ import annotations

import re
import ssl
import datetime as dt
from pathlib import Path

import certifi
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
TABLES = PROJECT / "tables"
REPORTS = PROJECT / "reports"
for d in [TABLES, REPORTS]:
    d.mkdir(parents=True, exist_ok=True)

SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
ssl._create_default_https_context = lambda: SSL_CONTEXT

URL = (PROJECT / "data/raw/latest_ncbi_metadata_url.txt").read_text().strip()
RUN_DATE = dt.datetime.now().astimezone().isoformat(timespec="seconds")
USECOLS = [
    "#label", "collection_date", "target_creation_date", "geo_loc_name", "host",
    "isolation_source", "source_type", "AMR_genotypes", "AMR_genotypes_core",
]
NONHUMAN = ["livestock", "companion_animal", "food", "environment"]


def extract_year_series(s: pd.Series) -> pd.Series:
    return s.astype("string").str.extract(r"(19\d{2}|20\d{2})", expand=False).astype("Int64")


def parse_country_series(s: pd.Series) -> pd.Series:
    out = s.astype("string").str.split(":", n=1).str[0].str.strip()
    out = out.mask(out.str.upper().isin(["NULL", "", "<NA>"]))
    return out


def classify_sources(df: pd.DataFrame) -> pd.Series:
    text = (
        df.get("host", pd.Series(index=df.index, dtype="string")).fillna("").astype(str).str.lower()
        + " | " + df.get("isolation_source", pd.Series(index=df.index, dtype="string")).fillna("").astype(str).str.lower()
        + " | " + df.get("source_type", pd.Series(index=df.index, dtype="string")).fillna("").astype(str).str.lower()
    )
    source = pd.Series("unknown", index=df.index, dtype="string")

    # Priority order mirrors prior feasibility script. Token-safe cat/dog terms avoid cat/cattle errors.
    livestock_re = re.compile(
        r"bos taurus|sus scrofa|gallus gallus|ovis aries|capra hircus|bubalus bubalis|"
        r"\b(cattle|cow|bovine|calf|swine|pig|porcine|pork|chicken|poultry|broiler|layer|turkey|duck|goose|sheep|ovine|goat|caprine|buffalo|camel)\b"
    )
    companion_re = re.compile(
        r"felis catus|canis lupus familiaris|\b(dog|dogs|cat|cats|feline|canine|horse|equine|pet|companion)\b"
    )
    env_re = re.compile(
        r"farm environment|\b(wastewater|sewage|effluent|environment|environmental|water|soil|surface|air|dust|sediment|sludge)\b"
    )
    food_re = re.compile(
        r"chicken meat|turkey meat|ready to eat|\b(food|meat|retail|milk|dairy|cheese|beef|pork|carcass|seafood|fish|egg|eggs)\b"
    )
    human_re = re.compile(
        r"homo sapiens|\b(human|patient|clinical|hospital|blood|wound|nasal|urine|sputum|respiratory|abscess|skin|infection|screening)\b"
    )

    masks = [
        ("livestock", text.str.contains(livestock_re, na=False)),
        ("companion_animal", text.str.contains(companion_re, na=False)),
        ("environment", text.str.contains(env_re, na=False)),
        ("food", text.str.contains(food_re, na=False)),
        ("human_or_clinical", text.str.contains(human_re, na=False)),
    ]
    assigned = pd.Series(False, index=df.index)
    for label, mask in masks:
        mask = mask & ~assigned
        source.loc[mask] = label
        assigned |= mask
    return source


def summarize_window(df: pd.DataFrame, label: str, year_col: str, start: int, end: int) -> dict:
    sub = df[(df[year_col] >= start) & (df[year_col] <= end)].copy()
    is_nonhuman = sub["source_category"].isin(NONHUMAN)
    return {
        "window": label,
        "year_definition": year_col,
        "start": start,
        "end": end,
        "records": len(sub),
        "mrsa_marker_positive": int(sub["mrsa_marker_positive"].sum()),
        "mrsa_pct": round(float(sub["mrsa_marker_positive"].mean() * 100), 2) if len(sub) else 0,
        "nonhuman": int(is_nonhuman.sum()),
        "nonhuman_pct": round(float(is_nonhuman.mean() * 100), 2) if len(sub) else 0,
        "livestock": int((sub["source_category"] == "livestock").sum()),
        "companion_animal": int((sub["source_category"] == "companion_animal").sum()),
        "food": int((sub["source_category"] == "food").sum()),
        "environment": int((sub["source_category"] == "environment").sum()),
        "human_or_clinical": int((sub["source_category"] == "human_or_clinical").sum()),
        "unknown": int((sub["source_category"] == "unknown").sum()),
        "unique_countries": int(sub["country_parsed"].nunique()),
    }


def trend_table(df: pd.DataFrame, year_col: str, start=2015, end=2025) -> pd.DataFrame:
    sub = df[(df[year_col] >= start) & (df[year_col] <= end)].copy()
    g = sub.groupby(year_col, dropna=False).agg(
        records=("label", "size"),
        mrsa_marker_positive=("mrsa_marker_positive", "sum"),
        nonhuman=("source_category", lambda s: int(s.isin(NONHUMAN).sum())),
        livestock=("source_category", lambda s: int((s == "livestock").sum())),
        companion_animal=("source_category", lambda s: int((s == "companion_animal").sum())),
        food=("source_category", lambda s: int((s == "food").sum())),
        environment=("source_category", lambda s: int((s == "environment").sum())),
        human_or_clinical=("source_category", lambda s: int((s == "human_or_clinical").sum())),
        unknown=("source_category", lambda s: int((s == "unknown").sum())),
        unique_countries=("country_parsed", "nunique"),
    ).reset_index().rename(columns={year_col: "year"})
    g["mrsa_pct"] = (g["mrsa_marker_positive"] / g["records"] * 100).round(2)
    g["nonhuman_pct"] = (g["nonhuman"] / g["records"] * 100).round(2)
    return g


def source_summary(df: pd.DataFrame, year_col: str, start: int, end: int) -> pd.DataFrame:
    sub = df[(df[year_col] >= start) & (df[year_col] <= end)].copy()
    out = sub.groupby("source_category", dropna=False).agg(
        records=("label", "size"),
        mrsa_marker_positive=("mrsa_marker_positive", "sum"),
        unique_countries=("country_parsed", "nunique"),
    ).reset_index().sort_values("records", ascending=False)
    out["mrsa_pct"] = (out["mrsa_marker_positive"] / out["records"] * 100).round(2)
    return out


def main():
    print(f"Reading NCBI metadata: {URL}")
    df = pd.read_csv(
        URL,
        sep="\t",
        usecols=lambda c: c in USECOLS,
        dtype=str,
        na_values=["NULL", "", "NA", "N/A"],
        keep_default_na=True,
        low_memory=False,
    ).rename(columns={"#label": "label"})
    df["collection_year"] = extract_year_series(df["collection_date"])
    df["deposition_year"] = extract_year_series(df["target_creation_date"])
    amr = (df["AMR_genotypes"].fillna("") + "," + df["AMR_genotypes_core"].fillna("")).str.lower()
    df["mrsa_marker_positive"] = amr.str.contains(r"\bmeca\b|\bmecc\b", regex=True, na=False)
    df["source_category"] = classify_sources(df)
    df["country_parsed"] = parse_country_series(df["geo_loc_name"])

    comparisons = pd.DataFrame([
        summarize_window(df, "deposition_2020_2025", "deposition_year", 2020, 2025),
        summarize_window(df, "collection_2020_2025", "collection_year", 2020, 2025),
        summarize_window(df, "deposition_2015_2025", "deposition_year", 2015, 2025),
        summarize_window(df, "collection_2015_2025", "collection_year", 2015, 2025),
    ])
    comparisons.to_csv(TABLES / "18_window_comparison_2020_2025_vs_2015_2025.csv", index=False)

    dep_trend = trend_table(df, "deposition_year")
    col_trend = trend_table(df, "collection_year")
    dep_trend.to_csv(TABLES / "19_deposition_year_trend.csv", index=False)
    col_trend.to_csv(TABLES / "20_collection_year_trend.csv", index=False)

    source_summary(df, "deposition_year", 2020, 2025).to_csv(TABLES / "21_source_summary_deposition_2020_2025.csv", index=False)
    source_summary(df, "collection_year", 2020, 2025).to_csv(TABLES / "22_source_summary_collection_2020_2025.csv", index=False)
    source_summary(df, "deposition_year", 2015, 2025).to_csv(TABLES / "23_source_summary_deposition_2015_2025.csv", index=False)
    source_summary(df, "collection_year", 2015, 2025).to_csv(TABLES / "24_source_summary_collection_2015_2025.csv", index=False)

    exact = df[df["deposition_year"].isin([2020, 2025])].groupby("deposition_year").agg(
        records=("label", "size"),
        mrsa_marker_positive=("mrsa_marker_positive", "sum"),
        nonhuman=("source_category", lambda s: int(s.isin(NONHUMAN).sum())),
        unique_countries=("country_parsed", "nunique"),
    ).reset_index().rename(columns={"deposition_year": "year"})
    exact.to_csv(TABLES / "25_exact_deposition_2020_and_2025.csv", index=False)

    report = REPORTS / "window_choice_2020_2025_vs_2015_2025.md"
    report.write_text(
        "# Window choice check: 2020–2025 vs 2015–2025\n\n"
        f"Run date: {RUN_DATE}\n\n"
        f"NCBI metadata URL: {URL}\n\n"
        "## Window summary\n\n"
        + comparisons.to_markdown(index=False)
        + "\n\n## Deposition-year trend, 2015–2025\n\n"
        + dep_trend.to_markdown(index=False)
        + "\n\n## Collection-year trend, 2015–2025\n\n"
        + col_trend.to_markdown(index=False)
        + "\n\n## Exact deposition years 2020 and 2025\n\n"
        + exact.to_markdown(index=False)
        + "\n\n## Interpretation\n\n"
        "The 2020–2025 window is feasible and cleaner for contemporary trends. "
        "The 2015–2025 window is stronger for historical context, trend stability, and source/income comparisons. "
        "A good manuscript design is to use 2015–2025 as the main denominator and explicitly analyze 2020–2025 as a contemporary post-2020 subperiod; "
        "if journal scope demands a sharper contemporary framing, use 2020–2025 as primary and keep 2015–2019 as historical baseline.\n",
        encoding="utf-8",
    )

    print("WINDOW SUMMARY")
    print(comparisons.to_string(index=False))
    print("\nDEPOSITION TREND")
    print(dep_trend.to_string(index=False))
    print("\nCOLLECTION TREND")
    print(col_trend.to_string(index=False))
    print(f"\nReport written: {report}")


if __name__ == "__main__":
    main()
