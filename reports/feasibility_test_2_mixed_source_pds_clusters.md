# Feasibility Test 2 — Mixed-source NCBI SNP cluster overlap
Run date: 2026-07-06T23:45:59+06:00
Cluster map URL: https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Staphylococcus_aureus/PDG000000073.1234/Clusters/PDG000000073.1234.reference_target.all_isolates.tsv
Rows in primary 2015–2025 dataset: 86,947
Rows with PDS cluster ID: 47,725 (54.89%)
Informative non-unknown source rows with PDS cluster ID: 43,035
Total informative PDS clusters: 9,408
Mixed-source PDS clusters (>=2 source categories): 224

## Cross-interface pattern counts

| pattern               |   clusters |   records_in_clusters |
|:----------------------|-----------:|----------------------:|
| human_livestock       |        115 |                  2113 |
| human_companion       |         19 |                   983 |
| human_food            |         26 |                   850 |
| human_environment     |         47 |                  1477 |
| livestock_food        |         29 |                   698 |
| livestock_environment |         26 |                   932 |
| three_or_more_sources |         21 |                  1479 |

## Source cluster availability

| source_category   |   records_with_cluster |   unique_pds_clusters |   mrsa_marker_positive |
|:------------------|-----------------------:|----------------------:|-----------------------:|
| human_or_clinical |                  38502 |                  8453 |                  25269 |
| livestock         |                   2895 |                   715 |                   1038 |
| food              |                    949 |                   293 |                    393 |
| companion_animal  |                    356 |                    90 |                    228 |
| environment       |                    333 |                   107 |                    250 |

## Top 20 mixed-source clusters

| PDS_acc         |   records | source_categories                                        |   n_source_categories | countries                                                                                      |   n_countries |   mrsa_marker_positive | years                                                                 |
|:----------------|----------:|:---------------------------------------------------------|----------------------:|:-----------------------------------------------------------------------------------------------|--------------:|-----------------------:|:----------------------------------------------------------------------|
| PDS000175068.20 |       298 | environment;food;human_or_clinical;livestock             |                     4 | Argentina;China;Czechia;Denmark;Finland;Germany;Hungary;Italy;Japan;Netherlands;Slovenia;Spain |            15 |                    296 | 2015.0;2016.0;2017.0;2018.0;2019.0;2020.0;2021.0;2022.0;2023.0        |
| PDS000238066.5  |       232 | companion_animal;environment;human_or_clinical;livestock |                     4 | Belgium;China;Denmark;Finland;Germany;Hungary;Italy;Netherlands;Slovenia                       |             9 |                    232 | 2015.0;2016.0;2017.0;2018.0;2019.0;2020.0;2021.0;2022.0;2023.0;2025.0 |
| PDS000093745.10 |        69 | companion_animal;environment;human_or_clinical;livestock |                     4 | Belgium;Canada;France;Netherlands                                                              |             4 |                     69 | 2017.0;2018.0;2019.0;2020.0;2021.0;2022.0;2023.0                      |
| PDS000069649.42 |        19 | environment;food;human_or_clinical;livestock             |                     4 | China;Colombia;Denmark;Korea, Rep.;Saudi Arabia;South Africa;United States                     |             7 |                      0 | 2015.0;2016.0;2017.0;2018.0;2019.0;2020.0;2021.0;2024.0               |
| PDS000069053.33 |        14 | environment;food;human_or_clinical;livestock             |                     4 | China;France;South Africa;Sri Lanka;Tanzania;United States                                     |             6 |                      0 | 2016.0;2017.0;2020.0;2021.0;2022.0;2023.0                             |
| PDS000144312.8  |       241 | companion_animal;human_or_clinical;livestock             |                     3 | Australia;Fiji;Netherlands;Samoa;not provided                                                  |             5 |                    241 | 2015.0;2016.0;2017.0;2018.0;2019.0;2020.0;2022.0                      |
| PDS000067777.38 |        98 | food;human_or_clinical;livestock                         |                     3 | Denmark;Germany;Hungary;Italy;Netherlands;Spain;Switzerland;Thailand                           |             8 |                     97 | 2015.0;2016.0;2017.0;2018.0;2019.0;2020.0;2021.0;2022.0               |
| PDS000226212.9  |        79 | companion_animal;human_or_clinical;livestock             |                     3 | Denmark;France;Netherlands;Switzerland                                                         |             4 |                     79 | 2015.0;2016.0;2017.0;2018.0;2019.0;2020.0;2021.0;2022.0;2023.0        |
| PDS000239414.3  |        76 | environment;human_or_clinical;livestock                  |                     3 | China                                                                                          |             1 |                     76 | 2016.0;2017.0;2018.0;2019.0;2020.0;2021.0;2022.0;2023.0               |
| PDS000226510.6  |        75 | food;human_or_clinical;livestock                         |                     3 | Denmark;Germany;Netherlands;Switzerland                                                        |             4 |                     75 | 2015.0;2016.0;2017.0;2018.0;2019.0;2020.0;2021.0;2022.0               |
| PDS000090674.17 |        59 | environment;food;human_or_clinical                       |                     3 | Australia;China;Denmark                                                                        |             3 |                     59 | 2015.0;2016.0;2017.0;2018.0;2019.0;2020.0;2021.0;2022.0               |
| PDS000068767.4  |        58 | environment;human_or_clinical;livestock                  |                     3 | China                                                                                          |             1 |                     58 | 2016.0                                                                |
| PDS000068755.17 |        43 | companion_animal;human_or_clinical;livestock             |                     3 | Czechia;Denmark;France;Germany;Netherlands;Switzerland                                         |             6 |                     43 | 2015.0;2016.0;2017.0;2018.0;2019.0;2020.0;2021.0;2022.0;2023.0        |
| PDS000170130.2  |        42 | food;human_or_clinical;livestock                         |                     3 | China                                                                                          |             1 |                      0 | 2016.0;2019.0                                                         |
| PDS000221219.5  |        17 | food;human_or_clinical;livestock                         |                     3 | Belgium;China;Czechia;Hungary;Japan;Spain                                                      |             6 |                     17 | 2016.0;2017.0;2018.0;2019.0;2020.0;2021.0                             |
| PDS000112001.4  |        15 | food;human_or_clinical;livestock                         |                     3 | China                                                                                          |             1 |                     15 | 2017.0;2019.0;2025.0                                                  |
| PDS000115307.14 |        14 | environment;human_or_clinical;livestock                  |                     3 | China                                                                                          |             1 |                     12 | 2019.0;2021.0;2022.0                                                  |
| PDS000181740.1  |         9 | environment;food;livestock                               |                     3 | China                                                                                          |             1 |                      0 | 2020.0                                                                |
| PDS000068724.24 |         8 | food;human_or_clinical;livestock                         |                     3 | Brazil;China;Russian Federation                                                                |             3 |                      0 | 2015.0;2017.0;2019.0;2021.0                                           |
| PDS000094049.1  |         7 | environment;human_or_clinical;livestock                  |                     3 | China                                                                                          |             1 |                      0 | 2020.0                                                                |

## Interpretation guardrail
Mixed PDS clusters are genomic-overlap / transmission-compatible signals only. They do not prove direct human–animal–environment transmission without epidemiologic linkage, sampling design, and household/farm metadata.

## Outputs
- `tables/13_mixed_source_pds_clusters.csv`
- `tables/14_cross_interface_cluster_patterns.csv`
- `tables/15_representative_mixed_cluster_records.csv`
- `tables/16_source_cluster_availability.csv`
- `data/processed/s_aureus_primary_2015_2025_with_pds_clusters.csv.gz`
