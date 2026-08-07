#!/usr/bin/env python3
"""
00_freeze_raw_sources.py — freeze every live-fetched raw source into data/raw/
so the whole pipeline (01, 02, 05, 07) becomes fully offline-reproducible.

Two remote dependencies are normally streamed at run time and never stored:

  1. NCBI Pathogen Detection Staphylococcus aureus AMR metadata TSV + the
     reference_target.all_isolates cluster TSV (a specific PDG snapshot).
  2. World Bank country classifications + a handful of indicators, fetched live
     from the World Bank REST API.

Because both are pulled live, re-running the pipeline later can silently produce
different numbers (NCBI rolls to a new PDG snapshot; the World Bank revises
indicator values). This script pins them to disk.

Snapshot is taken from data/raw/latest_ncbi_metadata_url.txt (the exact PDG the
published results were built on) -- NOT the newest snapshot on the server -- so
the freeze reproduces the existing outputs rather than a moving target.

Writes to data/raw/:
  * ncbi_s_aureus_amr_metadata.<PDG>.tsv.gz          (full metadata, gzipped)
  * ncbi_s_aureus_clusters_all_isolates.<PDG>.tsv.gz (full clusters, gzipped)
  * worldbank_countries.csv                          (iso3/income/region universe)
  * worldbank_indicators.csv                         (tidy: indicator,iso3,year,value)
  * frozen_sources_manifest.json                     (URLs, sha256, sizes, fetch time)

Usage:
  python3 scripts/00_freeze_raw_sources.py            # download + freeze
  python3 scripts/00_freeze_raw_sources.py --verify   # re-hash local files vs manifest
"""
from __future__ import annotations

import argparse
import datetime as dt
import gzip
import hashlib
import json
import ssl
import sys
import urllib.request
from pathlib import Path

import certifi
import pandas as pd

SSL_CONTEXT = ssl.create_default_context(cafile=certifi.where())

PROJECT = Path(__file__).resolve().parents[1]
RAW = PROJECT / "data" / "raw"
RAW.mkdir(parents=True, exist_ok=True)

NCBI_BASE = "https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Staphylococcus_aureus/"
URL_PIN = RAW / "latest_ncbi_metadata_url.txt"
MANIFEST = RAW / "frozen_sources_manifest.json"

# World Bank indicators to freeze (superset of what 01/05/07 consume).
WB_INDICATORS = {
    "SP.POP.TOTL": "Population, total",
    "NY.GDP.PCAP.CD": "GDP per capita (current US$)",
    "SH.XPD.CHEX.PC.CD": "Current health expenditure per capita (current US$)",
    "SH.MED.PHYS.ZS": "Physicians (per 1,000 people)",
    "SH.MED.BEDS.ZS": "Hospital beds (per 1,000 people)",
}
WB_YEAR_START, WB_YEAR_END = 2015, 2025


def pinned_snapshot() -> tuple[str, str, str]:
    """Return (pdg, metadata_url, clusters_url) from the pinned URL file."""
    if not URL_PIN.exists():
        sys.exit(f"ERROR: {URL_PIN} not found; cannot determine the pinned PDG snapshot.")
    meta_url = URL_PIN.read_text().strip()
    pdg = meta_url.split("/AMR/", 1)[0].rstrip("/").split("/")[-1]
    clusters_url = f"{NCBI_BASE}{pdg}/Clusters/{pdg}.reference_target.all_isolates.tsv"
    return pdg, meta_url, clusters_url


def download_gzip(url: str, dest: Path) -> dict:
    """Stream a (large) text URL to a gzip file, hashing the *uncompressed* bytes.

    gzip mtime is fixed to 0 so the compressed file is reproducible byte-for-byte.
    Returns manifest fields for this source.
    """
    print(f"  downloading {url}")
    sha = hashlib.sha256()
    n_raw = 0
    req = urllib.request.Request(url, headers={"User-Agent": "mrsa-onehealth-freeze/1.0"})
    with urllib.request.urlopen(req, timeout=120, context=SSL_CONTEXT) as resp:
        declared = resp.headers.get("Content-Length")
        last_mod = resp.headers.get("Last-Modified")
        tmp = dest.with_suffix(dest.suffix + ".part")
        with gzip.GzipFile(filename="", mode="wb", fileobj=open(tmp, "wb"), mtime=0) as gz:
            while True:
                chunk = resp.read(1 << 20)  # 1 MiB
                if not chunk:
                    break
                n_raw += len(chunk)
                sha.update(chunk)
                gz.write(chunk)
                if n_raw % (32 << 20) < (1 << 20):
                    print(f"    {n_raw/1e6:8.1f} MB", end="\r", flush=True)
    if declared is not None and int(declared) != n_raw:
        tmp.unlink(missing_ok=True)
        sys.exit(f"ERROR: size mismatch for {url}: declared {declared}, got {n_raw}")
    tmp.replace(dest)
    print(f"    saved {dest.name}  ({n_raw:,} raw bytes -> {dest.stat().st_size:,} gz bytes)")
    return {
        "url": url,
        "http_last_modified": last_mod,
        "uncompressed_bytes": n_raw,
        "uncompressed_sha256": sha.hexdigest(),
        "local_file": dest.name,
        "local_gz_bytes": dest.stat().st_size,
        "local_gz_sha256": _sha256_file(dest),
    }


def _sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def wb_fetch(url: str, timeout: int = 90):
    req = urllib.request.Request(url, headers={"User-Agent": "mrsa-onehealth-freeze/1.0"})
    with urllib.request.urlopen(req, timeout=timeout, context=SSL_CONTEXT) as r:
        return json.loads(r.read())


def freeze_worldbank() -> dict:
    """Freeze WB country universe + indicator series to CSV. Returns manifest fields."""
    print("  fetching World Bank country classifications ...")
    data = wb_fetch("https://api.worldbank.org/v2/country/all?format=json&per_page=400")
    crows = []
    for r in data[1]:
        if r.get("region", {}).get("id") == "NA":  # drop aggregates, keep real countries
            continue
        crows.append({
            "iso3": r.get("id"),
            "iso2": r.get("iso2Code"),
            "country_name_wb": r.get("name"),
            "region": r.get("region", {}).get("value"),
            "income_id": r.get("incomeLevel", {}).get("id"),
            "income_group": r.get("incomeLevel", {}).get("value"),
        })
    countries = pd.DataFrame(crows).sort_values("iso3").reset_index(drop=True)
    countries.to_csv(RAW / "worldbank_countries.csv", index=False)
    print(f"    saved worldbank_countries.csv  ({len(countries)} countries)")

    print(f"  fetching {len(WB_INDICATORS)} World Bank indicators, {WB_YEAR_START}-{WB_YEAR_END} ...")
    irows = []
    for ind, name in WB_INDICATORS.items():
        url = (f"https://api.worldbank.org/v2/country/all/indicator/{ind}"
               f"?format=json&per_page=20000&date={WB_YEAR_START}:{WB_YEAR_END}")
        data = wb_fetch(url)
        got = 0
        for r in data[1] if len(data) > 1 and data[1] else []:
            if r.get("value") is None:
                continue
            iso3 = r.get("countryiso3code")
            if not iso3:
                continue
            irows.append({"indicator": ind, "iso3": iso3, "year": int(r["date"]), "value": r["value"]})
            got += 1
        print(f"    {ind:20s} {got:6d} non-null obs")
    indicators = pd.DataFrame(irows).sort_values(["indicator", "iso3", "year"]).reset_index(drop=True)
    indicators.to_csv(RAW / "worldbank_indicators.csv", index=False)
    print(f"    saved worldbank_indicators.csv  ({len(indicators)} rows)")
    return {
        "api": "https://api.worldbank.org/v2",
        "indicators": WB_INDICATORS,
        "year_range": [WB_YEAR_START, WB_YEAR_END],
        "countries_file": "worldbank_countries.csv",
        "countries_sha256": _sha256_file(RAW / "worldbank_countries.csv"),
        "n_countries": int(len(countries)),
        "indicators_file": "worldbank_indicators.csv",
        "indicators_sha256": _sha256_file(RAW / "worldbank_indicators.csv"),
        "n_indicator_rows": int(len(indicators)),
    }


def do_freeze() -> None:
    pdg, meta_url, clusters_url = pinned_snapshot()
    print(f"Pinned snapshot: {pdg}")
    meta_dest = RAW / f"ncbi_s_aureus_amr_metadata.{pdg}.tsv.gz"
    clust_dest = RAW / f"ncbi_s_aureus_clusters_all_isolates.{pdg}.tsv.gz"

    manifest = {
        "frozen_at": dt.datetime.now().astimezone().isoformat(timespec="seconds"),
        "pdg_snapshot": pdg,
        "ncbi": {},
        "worldbank": {},
    }
    print("NCBI metadata:")
    manifest["ncbi"]["metadata"] = download_gzip(meta_url, meta_dest)
    print("NCBI clusters:")
    manifest["ncbi"]["clusters"] = download_gzip(clusters_url, clust_dest)
    print("World Bank:")
    manifest["worldbank"] = freeze_worldbank()

    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(f"\nWrote manifest: {MANIFEST.relative_to(PROJECT)}")
    print("Freeze complete. The pipeline (01/02/05/07) will now read these local files first.")


def do_verify() -> None:
    if not MANIFEST.exists():
        sys.exit(f"ERROR: {MANIFEST} not found; run the freeze first.")
    m = json.loads(MANIFEST.read_text())
    ok = True

    def check(path: Path, expect_sha: str, label: str) -> None:
        nonlocal ok
        if not path.exists():
            print(f"  MISSING  {label}: {path.name}")
            ok = False
            return
        actual = _sha256_file(path)
        status = "OK  " if actual == expect_sha else "FAIL"
        if actual != expect_sha:
            ok = False
        print(f"  {status} {label}: {path.name}")

    print(f"Verifying against manifest frozen_at {m.get('frozen_at')} (snapshot {m.get('pdg_snapshot')})")
    check(RAW / m["ncbi"]["metadata"]["local_file"], m["ncbi"]["metadata"]["local_gz_sha256"], "NCBI metadata")
    check(RAW / m["ncbi"]["clusters"]["local_file"], m["ncbi"]["clusters"]["local_gz_sha256"], "NCBI clusters")
    check(RAW / m["worldbank"]["countries_file"], m["worldbank"]["countries_sha256"], "WB countries")
    check(RAW / m["worldbank"]["indicators_file"], m["worldbank"]["indicators_sha256"], "WB indicators")
    print("All frozen sources verified." if ok else "VERIFICATION FAILED.")
    sys.exit(0 if ok else 1)


def do_wb_only() -> None:
    """Re-freeze only the World Bank CSVs (no 158 MB NCBI re-download) and update
    the manifest's worldbank section. Use this if a freeze captured a truncated-
    precision WB API response (the WB API's load-balanced nodes serialize floats
    inconsistently) and the covariate columns no longer reproduce the tables."""
    print("World Bank (re-freeze only):")
    wb = freeze_worldbank()
    if MANIFEST.exists():
        m = json.loads(MANIFEST.read_text())
    else:
        m = {"pdg_snapshot": pinned_snapshot()[0], "ncbi": {}}
    m["worldbank"] = wb
    m["worldbank_refrozen_at"] = dt.datetime.now().astimezone().isoformat(timespec="seconds")
    MANIFEST.write_text(json.dumps(m, indent=2) + "\n", encoding="utf-8")
    print(f"Updated manifest worldbank section: {MANIFEST.relative_to(PROJECT)}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify", action="store_true", help="re-hash local frozen files against the manifest")
    ap.add_argument("--wb-only", action="store_true", help="re-freeze only the World Bank CSVs (skip NCBI download)")
    args = ap.parse_args()
    if args.verify:
        do_verify()
    elif args.wb_only:
        do_wb_only()
    else:
        do_freeze()
