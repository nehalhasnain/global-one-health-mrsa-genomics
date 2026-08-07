#!/usr/bin/env python3
"""Feasibility Test 2: NCBI SNP cluster overlap across One Health source categories."""
from __future__ import annotations

import datetime as dt
import ssl
import urllib.request
from pathlib import Path

import certifi
import pandas as pd

PROJECT = Path(__file__).resolve().parents[1]
TABLES = PROJECT / "tables"
REPORTS = PROJECT / "reports"
RAW = PROJECT / "data" / "raw"
PROCESSED = PROJECT / "data" / "processed"
for d in [TABLES, REPORTS, RAW, PROCESSED]:
    d.mkdir(parents=True, exist_ok=True)

SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())
ssl._create_default_https_context = lambda: SSL_CONTEXT

NCBI_BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Staphylococcus_aureus/"
PDG = (RAW / "latest_ncbi_metadata_url.txt").read_text().split("/AMR/", 1)[0].rstrip("/").split("/")[-1]
CLUSTER_URL = f"{NCBI_BASE}{PDG}/Clusters/{PDG}.reference_target.all_isolates.tsv"
RUN_DATE = dt.datetime.now().astimezone().isoformat(timespec="seconds")


def main():
    rows_path = PROCESSED / "s_aureus_primary_collection_year_2015_2025_feasibility_rows.csv.gz"
    rows = pd.read_csv(rows_path, dtype=str, low_memory=False)
    for b in ["mrsa_marker_positive", "mecA_positive", "mecC_positive"]:
        if b in rows.columns:
            rows[b] = rows[b].map({"True": True, "False": False, True: True, False: False}).fillna(False).astype(bool)

    # Offline-reproducible: prefer the frozen clusters snapshot (00_freeze_raw_sources.py).
    frozen_clusters = RAW / f"ncbi_s_aureus_clusters_all_isolates.{PDG}.tsv.gz"
    cluster_src = frozen_clusters if frozen_clusters.exists() else CLUSTER_URL
    print(f"Reading cluster map: {cluster_src}")
    clusters = pd.read_csv(cluster_src, sep="\t", dtype=str, na_values=["NULL", "", "NA"], keep_default_na=True)
    clusters.to_csv(RAW / "ncbi_s_aureus_cluster_map_sample_1000.csv", index=False)
    merged = rows.merge(clusters[["target_acc", "PDS_acc"]], on="target_acc", how="left")
    merged.to_csv(PROCESSED / "s_aureus_primary_2015_2025_with_pds_clusters.csv.gz", index=False, compression="gzip")

    clustered = merged[merged["PDS_acc"].notna()].copy()
    informative = clustered[~clustered["source_category"].isin(["unknown"])]

    # Cluster-level summary.
    grp = informative.groupby("PDS_acc").agg(
        records=("target_acc", "nunique"),
        source_categories=("source_category", lambda s: ";".join(sorted(set(s.dropna())))),
        n_source_categories=("source_category", lambda s: len(set(s.dropna()))),
        countries=("country_parsed", lambda s: ";".join(sorted(set([x for x in s.dropna().astype(str) if x and x != 'nan']))[:12])),
        n_countries=("country_parsed", lambda s: len(set([x for x in s.dropna().astype(str) if x and x != 'nan']))),
        mrsa_marker_positive=("mrsa_marker_positive", "sum"),
        years=("collection_year", lambda s: ";".join(sorted(set(s.dropna().astype(str)))[:12])),
    ).reset_index()
    mixed = grp[grp["n_source_categories"] >= 2].sort_values(["n_source_categories", "records"], ascending=False)
    mixed.to_csv(TABLES / "13_mixed_source_pds_clusters.csv", index=False)

    # Count specific cross-interface patterns.
    patterns = []
    source_sets = grp["source_categories"].fillna("").map(lambda x: set(x.split(";")) if x else set())
    for label, needed in [
        ("human_livestock", {"human_or_clinical", "livestock"}),
        ("human_companion", {"human_or_clinical", "companion_animal"}),
        ("human_food", {"human_or_clinical", "food"}),
        ("human_environment", {"human_or_clinical", "environment"}),
        ("livestock_food", {"livestock", "food"}),
        ("livestock_environment", {"livestock", "environment"}),
        ("three_or_more_sources", set()),
    ]:
        if label == "three_or_more_sources":
            mask = grp["n_source_categories"] >= 3
        else:
            mask = source_sets.map(lambda s: needed.issubset(s))
        patterns.append({"pattern": label, "clusters": int(mask.sum()), "records_in_clusters": int(grp.loc[mask, "records"].sum())})
    pd.DataFrame(patterns).to_csv(TABLES / "14_cross_interface_cluster_patterns.csv", index=False)

    # Representative clusters for manual inspection.
    reps = []
    for pds in mixed.head(25)["PDS_acc"]:
        sub = merged[merged["PDS_acc"] == pds].copy()
        cols = [c for c in ["PDS_acc", "target_acc", "collection_year", "country_parsed", "source_category", "host", "isolation_source", "source_type", "mrsa_marker_positive", "AMR_genotypes_core"] if c in sub.columns]
        reps.append(sub[cols].head(50))
    if reps:
        pd.concat(reps).to_csv(TABLES / "15_representative_mixed_cluster_records.csv", index=False)

    source_cluster_summary = informative.groupby("source_category").agg(
        records_with_cluster=("target_acc", "nunique"),
        unique_pds_clusters=("PDS_acc", "nunique"),
        mrsa_marker_positive=("mrsa_marker_positive", "sum"),
    ).reset_index().sort_values("records_with_cluster", ascending=False)
    source_cluster_summary.to_csv(TABLES / "16_source_cluster_availability.csv", index=False)

    report = []
    report.append("# Feasibility Test 2 — Mixed-source NCBI SNP cluster overlap\n")
    report.append(f"Run date: {RUN_DATE}\n")
    report.append(f"Cluster map URL: {CLUSTER_URL}\n")
    report.append(f"Rows in primary 2015–2025 dataset: {len(rows):,}\n")
    report.append(f"Rows with PDS cluster ID: {len(clustered):,} ({len(clustered)/len(rows)*100:.2f}%)\n")
    report.append(f"Informative non-unknown source rows with PDS cluster ID: {len(informative):,}\n")
    report.append(f"Total informative PDS clusters: {grp['PDS_acc'].nunique():,}\n")
    report.append(f"Mixed-source PDS clusters (>=2 source categories): {len(mixed):,}\n")
    report.append("\n## Cross-interface pattern counts\n\n")
    patterns_df = pd.read_csv(TABLES / "14_cross_interface_cluster_patterns.csv")
    report.append(patterns_df.to_markdown(index=False))
    report.append("\n\n## Source cluster availability\n\n")
    report.append(source_cluster_summary.to_markdown(index=False))
    report.append("\n\n## Top 20 mixed-source clusters\n\n")
    report.append(mixed.head(20).to_markdown(index=False))
    report.append("\n\n## Interpretation guardrail\n")
    report.append("Mixed PDS clusters are genomic-overlap / transmission-compatible signals only. They do not prove direct human–animal–environment transmission without epidemiologic linkage, sampling design, and household/farm metadata.\n")
    report.append("\n## Outputs\n")
    for name in ["13_mixed_source_pds_clusters.csv", "14_cross_interface_cluster_patterns.csv", "15_representative_mixed_cluster_records.csv", "16_source_cluster_availability.csv"]:
        report.append(f"- `tables/{name}`\n")
    report.append("- `data/processed/s_aureus_primary_2015_2025_with_pds_clusters.csv.gz`\n")
    out = REPORTS / "feasibility_test_2_mixed_source_pds_clusters.md"
    out.write_text("".join(report), encoding="utf-8")
    print(f"Report written: {out}")
    print(f"Mixed-source clusters: {len(mixed):,}")


if __name__ == "__main__":
    main()
