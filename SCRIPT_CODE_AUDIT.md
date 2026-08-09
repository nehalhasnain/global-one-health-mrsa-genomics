# Script and archive audit

Generated: 2026-07-21T11:11:21+00:00

## Checks performed

- Python syntax compilation for 8 scripts: return code 0.
- `scripts/00_freeze_raw_sources.py --verify`: return code 0.
- Offline smoke test: `scripts/07_appendix.py`: return code 0.

## Classification

- `00_freeze_raw_sources.py`: public-source freezer and checksum verifier; network only when run without `--verify`.
- `01_feasibility_ncbi_s_aureus.py` through `05_inequity.py`: primary extraction and statistical pipeline; frozen local inputs are preferred when present.
- `06_figures.py`: main-text figure generation.
- `07_appendix.py`: supplementary tables 38–42, supplementary Figures S1–S3, and appendix report generation; tested offline in this archive.

## Important interpretation boundary

This is a public repository visibility analysis. The scripts do not estimate MRSA incidence or prevalence in the underlying population and do not establish direct transmission from mixed-source clusters.

## Output status

The archive contains the complete CSV outputs used by the manuscript, figures in PNG/PDF formats, frozen public inputs, processed compressed data, code, reports, and checksums. A Zenodo DOI must be inserted into the manuscript and `CITATION.cff` after Zenodo creates the record.

## Captured command output

### Python compile

```text

```

### Frozen-input verification

```text
Verifying against manifest frozen_at 2026-07-09T22:13:14+06:00 (snapshot PDG000000073.1234)
  OK   NCBI metadata: ncbi_s_aureus_amr_metadata.PDG000000073.1234.tsv.gz
  OK   NCBI clusters: ncbi_s_aureus_clusters_all_isolates.PDG000000073.1234.tsv.gz
  OK   WB countries: worldbank_countries.csv
  OK   WB indicators: worldbank_indicators.csv
All frozen sources verified.
/opt/anaconda3/lib/python3.12/site-packages/pandas/core/computation/expressions.py:23: UserWarning: Pandas requires version '2.10.2' or newer of 'numexpr' (version '2.8.7' currently installed).
  from pandas.core.computation.check import NUMEXPR_INSTALLED
/opt/anaconda3/lib/python3.12/site-packages/pandas/core/arrays/masked.py:56: UserWarning: Pandas requires version '1.4.2' or newer of 'bottleneck' (version '1.3.7' currently installed).
  from pandas.core import (

```

### Offline appendix smoke test

```text
Fetching World Bank region + population (certifi SSL)...
  wrote figures/figS1_geographic_concentration.png / .pdf
  wrote figures/figS2_temporal_dynamics.png / .pdf
  wrote figures/figS3_regional_visibility.png / .pdf
Wrote tables 38-42, figures S1-S3 and reports/07_appendix.md
Gini across countries = 0.885; top-10 share = 81.3%
/opt/anaconda3/lib/python3.12/site-packages/pandas/core/computation/expressions.py:23: UserWarning: Pandas requires version '2.10.2' or newer of 'numexpr' (version '2.8.7' currently installed).
  from pandas.core.computation.check import NUMEXPR_INSTALLED
/opt/anaconda3/lib/python3.12/site-packages/pandas/core/arrays/masked.py:56: UserWarning: Pandas requires version '1.4.2' or newer of 'bottleneck' (version '1.3.7' currently installed).
  from pandas.core import (

```
