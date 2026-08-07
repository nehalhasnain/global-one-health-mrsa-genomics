# Technical Appendix — supporting analyses (script 07)

Primary cohort: collection-year 2015–2025, n=86,947.

## Source-classification QA (table 38)

| source | n | % total | % host | % isol.source | % any field | n no-field |
|---|---|---|---|---|---|---|
| human_or_clinical | 65,850 | 75.74 | 99.5 | 83.9 | 100.0 | 0 |
| livestock | 4,005 | 4.61 | 88.7 | 75.7 | 100.0 | 0 |
| companion_animal | 930 | 1.07 | 95.9 | 68.5 | 100.0 | 0 |
| food | 1,815 | 2.09 | 9.1 | 97.3 | 100.0 | 0 |
| environment | 433 | 0.5 | 24.7 | 99.1 | 100.0 | 0 |
| unknown | 13,914 | 16.0 | 8.4 | 7.2 | 12.1 | 12,230 |

**Of 13,914 'unknown' genomes, 12,230 (88%) carried no host, isolation-source or source-type field at all — 'unknown' reflects missing metadata, not classifier failure.**

## Geographic concentration (table 39)

| stratum | genomes | countries | Gini | top country | top-1 % | top-3 % | top-10 % |
|---|---|---|---|---|---|---|---|
| all_sources | 86,219 | 115 | 0.885 | United States | 29.8 | 51.4 | 81.3 |
| mrsa_marker_positive | 47,553 | 96 | 0.886 | United States | 32.5 | 58.4 | 86.4 |
| livestock | 4,005 | 49 | 0.756 | China | 31.3 | 58.5 | 78.6 |
| companion_animal | 930 | 20 | 0.748 | Japan | 43.0 | 72.7 | 95.2 |
| food | 1,815 | 39 | 0.793 | China | 37.9 | 63.9 | 85.2 |
| environment | 433 | 24 | 0.781 | China | 52.4 | 74.1 | 91.5 |

## Regional visibility (table 40)

| region | records | per million | non-human/M | coverage |
|---|---|---|---|---|
| Europe & Central Asia | 33,921 | 36.679 | 2.0664 | 36/58 (62.1%) |
| North America | 25,980 | 68.936 | 1.6398 | 2/3 (66.7%) |
| East Asia & Pacific | 18,699 | 7.92 | 1.6205 | 12/37 (32.4%) |
| Middle East, North Africa, Afghanistan & Pakistan | 2,213 | 2.774 | 0.2131 | 15/23 (65.2%) |
| Latin America & Caribbean | 2,195 | 3.338 | 0.2752 | 10/42 (23.8%) |
| Sub-Saharan Africa | 1,535 | 1.219 | 0.2026 | 18/48 (37.5%) |
| South Asia | 457 | 0.275 | 0.0589 | 4/6 (66.7%) |

## Window & date-definition sensitivity (table 41)

| window | records | MRSA % | non-human % | countries |
|---|---|---|---|---|
| deposition_2020_2025 | 146,304 | 54.45 | 7.42 | 141 |
| collection_2020_2025 | 39,840 | 47.71 | 8.4 | 91 |
| deposition_2015_2025 | 155,137 | 55.31 | 7.5 | 149 |
| collection_2015_2025 | 86,947 | 55.39 | 8.26 | 117 |

## Guardrails

- Concentration and single-country dominance are properties of the public repository (a convenience sample), not of MRSA biology; they quantify how narrow the evidence base is.
- Regional per-million rates measure sequencing/deposition intensity, not disease burden.
- The collection-year decline after ~2022 reflects submission lag (recent infections not yet deposited), not a real fall in incidence — hence the deposition-year comparison in Fig S2.
