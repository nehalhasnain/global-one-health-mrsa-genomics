# Phase 4 — Global inequity in public MRSA/S. aureus genomic sequence visibility

Primary window: collection-year 2015-2025 (n=86,947). World Bank covariates fetched live (population, GDP pc, health exp pc, physicians/1k, beds/1k).

## Visibility by income group (per-million normalised)

| income group | records | records/million | non-human/million | unique countries |
|---|---|---|---|---|
| High income | 70,700 | 51.056 | 3.0005 | 42 |
| Upper middle income | 12,403 | 4.056 | 0.8509 | 29 |
| Lower middle income | 1,535 | 0.535 | 0.0857 | 17 |
| Low income | 362 | 0.498 | 0.077 | 9 |

**High-income vs low-income records/million ratio ≈ 102.5:1.**

## Country coverage (share of countries with ≥1 deposited genome)

| income group | with genome | total countries | coverage % |
|---|---|---|---|
| High income | 42 | 86 | 48.8 |
| Upper middle income | 29 | 59 | 49.2 |
| Lower middle income | 17 | 47 | 36.2 |
| Low income | 9 | 25 | 36.0 |

## Metadata completeness by income group (%)

| income group | n | host | source | collection date |
|---|---|---|---|---|
| High income | 70,700 | 81.6 | 71.3 | 100.0 |
| Upper middle income | 12,403 | 84.4 | 70.7 | 100.0 |
| Lower middle income | 1,535 | 92.9 | 73.9 | 100.0 |
| Low income | 362 | 99.7 | 91.7 | 100.0 |

## Submission lag (deposition − collection year) by income

Kruskal-Wallis p = 1.64e-101.

| income group | n | median lag (yr) | IQR |
|---|---|---|---|
| High income | 70,700 | 4.0 | 3.0-5.0 |
| Upper middle income | 12,402 | 4.0 | 3.0-5.0 |
| Lower middle income | 1,535 | 3.0 | 2.0-5.0 |
| Low income | 362 | 2.0 | 0.0-6.0 |

## Negative-binomial regression — genome count ~ income, offset log(population)

Reference = High income. IRR<1 means fewer deposited genomes per capita.

| term | IRR | 95% CI | p |
|---|---|---|---|
| Upper middle income | 0.715 | 0.508-1.006 | 5.44e-02 |
| Lower middle income | 0.068 | 0.045-0.103 | 2.28e-37 |
| Low income | 0.08 | 0.047-0.135 | 2.66e-21 |

## AMR/virulence gene profile by income (human isolates only)

| income group | n | mecA | acq.tet | erm | msr(A) | FQ-QRDR | PVL | scn |
|---|---|---|---|---|---|---|---|---|
| High income | 53,752 | 64.8 | 16.2 | 31.3 | 22.5 | 36.4 | 23.1 | 87.5 |
| Upper middle income | 9,323 | 53.1 | 25.4 | 51.1 | 3.1 | 29.9 | 14.4 | 86.7 |
| Lower middle income | 1,154 | 52.3 | 32.7 | 24.1 | 16.1 | 46.5 | 46.7 | 95.2 |
| Low income | 306 | 9.5 | 31.4 | 14.4 | 5.2 | 12.1 | 28.4 | 75.2 |

## Contemporary sensitivity (collection-year 2020-2025)

| income group | records | records/million | non-human/million |
|---|---|---|---|
| High income | 33,693 | 24.331 | 1.5353 |
| Upper middle income | 4,271 | 1.397 | 0.3463 |
| Lower middle income | 689 | 0.24 | 0.0289 |
| Low income | 231 | 0.318 | 0.0413 |

## Guardrails

- **Visibility ≠ burden.** Per-million rates measure public *sequencing/deposition* intensity, not MRSA incidence or prevalence.
- The gene-by-income profile reflects **which genomes were deposited** (lineage mix + assembly quality, both income-correlated), not AMR burden by healthcare strength; interpret descriptively.
- The regression is over countries **with ≥1 genome**; the ~half of low-income countries with zero genomes (coverage table) are the strongest inequity signal and sit outside the model.
