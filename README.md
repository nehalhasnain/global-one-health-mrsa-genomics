# Global One Health MRSA genomic sequence visibility at the human–animal–environment interface, 2015–2025

Zenodo-ready research archive accompanying the manuscript by Hasnain N, Mannan MA, and Hossain S.

## What this archive contains

This archive contains the manuscript source, a readable Technical Appendix, complete machine-readable supplementary tables S1–S11, all main and supplementary figures, the Python analysis/figure scripts, frozen public input files, processed analysis files, reports, checksums, and an analysis manifest.

The source of truth for the manuscript text is `manuscript/MANUSCRIPT_v2.md`. The source of truth for the supplementary methodology and archive crosswalk is `manuscript/TECHNICAL_APPENDIX.md`. The complete numerical tables are the CSV files under `tables/`; they are intentionally not duplicated as hundreds of wide Markdown rows inside the appendix.

## Study snapshot

- NCBI Pathogen Detection snapshot: `PDG000000073.1234`
- Primary cohort: collection years 2015–2025, n = 86,947
- MRSA marker: `mecA` or `mecC`
- One Health categories: human/clinical, livestock, companion animal, food, environment, unknown
- World Bank country and indicator data: frozen locally for this analysis version

All results are measures of public genomic visibility and deposition intensity. They must not be interpreted as MRSA incidence, prevalence, disease burden, or proof of direct transmission.

## Directory structure

```text
manuscript/       Manuscript and readable Technical Appendix
 tables/          Complete S1–S11 CSVs plus all upstream derived tables
figures/          Figures 1–5 and S1–S3 as PNG and PDF
scripts/          Python scripts 00–07
reports/          Script-generated analysis reports
data/raw/         Frozen public NCBI and World Bank inputs
data/processed/   Frozen compressed processed datasets used downstream
```

## Reproducibility

Use Python 3.10+ in an environment containing the packages listed in `requirements.txt`.

```bash
python -m pip install -r requirements.txt
python scripts/00_freeze_raw_sources.py --verify
python scripts/01_feasibility_ncbi_s_aureus.py
python scripts/02_feasibility_mixed_source_clusters.py
python scripts/03_window_choice_2020_2025_vs_2015_2025.py
python scripts/04_source_stats.py
python scripts/05_inequity.py
python scripts/06_figures.py
python scripts/07_appendix.py
```

The `--verify` mode is offline. Scripts prefer the frozen inputs under `data/raw/` and the frozen processed file under `data/processed/`. Script `00` without `--verify` is the only intentional network-refresh operation and can produce a different analysis if providers revise their data.

## Provenance

- `data/raw/frozen_sources_manifest.json` records the pinned NCBI URLs, World Bank source, sizes, and hashes.

## Citation

Hasnain N, Mannan MA, Hossain S. Global One Health MRSA genomic sequence visibility at the human–animal–environment interface, 2015–2025. Zenodo archive, version 1.0.0. DOI: to be assigned by Zenodo.

Please cite the Zenodo DOI after deposit and retain the NCBI Pathogen Detection snapshot identifier in any reuse.

## License

The Zenodo license should be selected at deposit after author agreement. A recommended split is CC BY 4.0 for manuscript/data/figures and an OSI-approved source-code license for the Python scripts. No license should be inferred from this working directory alone.
