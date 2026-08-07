# Feasibility Test 1 — NCBI Pathogen Detection Staphylococcus aureus / MRSA, 2015–2025
Run date: 2026-07-06T23:39:23+06:00
Latest NCBI PDG snapshot: `PDG000000073.1234`
Metadata URL: https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Staphylococcus_aureus/PDG000000073.1234/AMR/PDG000000073.1234.amr.metadata.tsv
Cluster all-isolates URL available: `True`
## Primary window definition
Primary filter: collection year 2015–2025. Records with missing collection year are not counted in the primary window; they are tracked in the sensitivity count.
## Headline feasibility counts
- Total NCBI rows all dates: 169,116
- Primary rows with collection year 2015–2025: 86,947
- Sensitivity rows using collection year or target-creation year 2015–2025: 124,904
- MRSA-marker-positive rows in primary window (mecA or mecC): 48,161 (55.39%)
- Non-human One Health rows in primary window (livestock + companion + food + environment): 7,183 (8.26%)
- Unknown source rows in primary window: 13,914 (16.0%)
- Rows matched to World Bank country metadata: 85,000 (97.76%)
## Source-category summary

| source_category   |   total_records |   mrsa_marker_positive |   mecA_positive |   mecC_positive |   unique_countries |   mrsa_marker_pct |
|:------------------|----------------:|-----------------------:|----------------:|----------------:|-------------------:|------------------:|
| human_or_clinical |           65850 |                  41705 |           41575 |             130 |                108 |             63.33 |
| unknown           |           13914 |                   3784 |            3675 |             109 |                 40 |             27.2  |
| livestock         |            4005 |                   1339 |            1322 |              17 |                 49 |             33.43 |
| food              |            1815 |                    595 |             595 |               0 |                 39 |             32.78 |
| companion_animal  |             930 |                    440 |             439 |               1 |                 20 |             47.31 |
| environment       |             433 |                    298 |             298 |               0 |                 24 |             68.82 |

## Income-group summary

| income_group        |   records |   mrsa_marker_positive |   nonhuman_records |   unique_countries |   mrsa_marker_pct |
|:--------------------|----------:|-----------------------:|-------------------:|-------------------:|------------------:|
| High income         |     70700 |                  39695 |               4155 |                 42 |             56.15 |
| Upper middle income |     12403 |                   6154 |               2602 |                 29 |             49.62 |
| nan                 |      1947 |                   1563 |                124 |                 18 |             80.28 |
| Lower middle income |      1535 |                    720 |                246 |                 17 |             46.91 |
| Low income          |       362 |                     29 |                 56 |                  9 |              8.01 |

## Top 20 countries by public S. aureus/MRSA sequence visibility

| country_parsed   | iso3   | income_group        | region                                            |   records |   mrsa_marker_positive |   nonhuman_records |   human_or_clinical_records |   unknown_source_records |   population_latest |   population_latest_year |   records_per_million_pop |   mrsa_per_million_pop |
|:-----------------|:-------|:--------------------|:--------------------------------------------------|----------:|-----------------------:|-------------------:|----------------------------:|-------------------------:|--------------------:|-------------------------:|--------------------------:|-----------------------:|
| United States    | USA    | High income         | North America                                     |     25672 |                  15442 |                591 |                       24952 |                      129 |         3.41785e+08 |                     2025 |                   75.1116 |                45.1805 |
| United Kingdom   | GBR    | High income         | Europe & Central Asia                             |      9955 |                   1302 |                 39 |                         572 |                     9344 |         6.9487e+07  |                     2025 |                  143.264  |                18.7373 |
| China            | CHN    | Upper middle income | East Asia & Pacific                               |      8684 |                   5121 |               2188 |                        6066 |                      430 |         1.40658e+09 |                     2025 |                    6.1738 |                 3.6407 |
| Netherlands      | NLD    | High income         | Europe & Central Asia                             |      7259 |                   7185 |                176 |                        6953 |                      130 |         1.80876e+07 |                     2025 |                  401.324  |               397.233  |
| Denmark          | DNK    | High income         | Europe & Central Asia                             |      5446 |                   3891 |                190 |                        5226 |                       30 |         6.00917e+06 |                     2025 |                  906.282  |               647.51   |
| Australia        | AUS    | High income         | East Asia & Pacific                               |      3678 |                   2539 |                 33 |                        3643 |                        2 |         2.76144e+07 |                     2025 |                  133.191  |                91.9447 |
| Japan            | JPN    | High income         | East Asia & Pacific                               |      3361 |                   1779 |                574 |                        2365 |                      422 |         1.23367e+08 |                     2025 |                   27.244  |                14.4204 |
| France           | FRA    | High income         | Europe & Central Asia                             |      2229 |                    411 |                225 |                        1725 |                      279 |         6.87203e+07 |                     2025 |                   32.4358 |                 5.9808 |
| Germany          | DEU    | High income         | Europe & Central Asia                             |      2064 |                   1042 |                335 |                        1372 |                      357 |         8.34912e+07 |                     2025 |                   24.7212 |                12.4803 |
| Saudi Arabia     | SAU    | High income         | Middle East, North Africa, Afghanistan & Pakistan |      1720 |                   1487 |                 59 |                         843 |                      818 |         3.69736e+07 |                     2025 |                   46.5197 |                40.2179 |
| Switzerland      | CHE    | High income         | Europe & Central Asia                             |      1678 |                   1296 |                 86 |                        1265 |                      327 |         9.09244e+06 |                     2025 |                  184.549  |               142.536  |
| New Zealand      | NZL    | High income         | East Asia & Pacific                               |      1627 |                    407 |                926 |                         701 |                        0 |         5.3247e+06  |                     2025 |                  305.557  |                76.4362 |
| Mexico           | MEX    | Upper middle income | Latin America & Caribbean                         |      1503 |                     64 |                 80 |                        1423 |                        0 |         1.31947e+08 |                     2025 |                   11.3909 |                 0.485  |
| Spain            | ESP    | High income         | Europe & Central Asia                             |      1098 |                    316 |                126 |                         694 |                      278 |         4.93551e+07 |                     2025 |                   22.2469 |                 6.4026 |
| Norway           | NOR    | High income         | Europe & Central Asia                             |       956 |                    357 |                 32 |                         745 |                      179 |         5.61087e+06 |                     2025 |                  170.384  |                63.6265 |
| nan              | nan    | nan                 | nan                                               |       728 |                    608 |                 96 |                         242 |                      390 |       nan           |                      nan |                  nan      |               nan      |
| Italy            | ITA    | High income         | Europe & Central Asia                             |       699 |                    323 |                172 |                         509 |                       18 |         5.89157e+07 |                     2025 |                   11.8644 |                 5.4824 |
| not provided     | nan    | nan                 | nan                                               |       585 |                    583 |                  0 |                         585 |                        0 |       nan           |                      nan |                  nan      |               nan      |
| Thailand         | THA    | Upper middle income | East Asia & Pacific                               |       457 |                    118 |                 60 |                         397 |                        0 |         7.16199e+07 |                     2025 |                    6.3809 |                 1.6476 |
| Brazil           | BRA    | Upper middle income | Latin America & Caribbean                         |       354 |                    279 |                 56 |                         293 |                        5 |         2.12812e+08 |                     2025 |                    1.6634 |                 1.311  |

## Preliminary go/no-go interpretation
- GO for the full One Health interface paper: the primary window has enough non-human records and enough mecA/mecC-positive records for source-stratified analysis.
- Interpret all counts as public genomic visibility, not incidence/prevalence/burden.
- Next test: inspect classifier QA tables and then join NCBI SNP cluster all-isolates metadata to find mixed human–animal–environment clusters.

## Output files
- `tables/01_overall_counts.csv`
- `tables/02_metadata_completeness.csv`
- `tables/03_source_mrsa_summary.csv`
- `tables/04_year_source_counts.csv`
- `tables/05_income_group_summary.csv`
- `tables/06_top100_country_visibility.csv`
- `tables/07_top100_country_nonhuman_visibility.csv`
- `tables/08_top_hosts.csv`
- `tables/09_top_isolation_sources.csv`
- `tables/10_top_source_types.csv`
- `tables/11_top_geo_loc_name.csv`
- `tables/12_top30_amr_genes_by_source.csv`
- `data/processed/s_aureus_primary_collection_year_2015_2025_feasibility_rows.csv.gz`
