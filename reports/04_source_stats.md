# Phase 3 — AMR / virulence / MRSA comparison across One Health source categories

Snapshot input: `s_aureus_primary_collection_year_2015_2025_feasibility_rows.csv.gz` (collection-year 2015-2025, n=86,947).

Statistics: per-gene omnibus Pearson chi-square across the five informative source categories (human, livestock, companion animal, food, environment), Benjamini-Hochberg FDR across genes; pairwise odds ratios vs human reference (Fisher exact p, Wald 95% CI, Haldane 0.5 correction for zero cells). 'unknown' source excluded from tests.

## Source denominators

| source | n |
|---|---|
| human_or_clinical | 65,850 |
| livestock | 4,005 |
| companion_animal | 930 |
| food | 1,815 |
| environment | 433 |
| unknown | 13,914 |

**32/33 genes differ significantly by source (q<0.05).**

## Prevalence by source (%), selected markers

| gene | human_or_clinical | livestock | companion_animal | food | environment |
|---|---|---|---|---|---|
| mecA (methicillin) | 63.14 | 33.01 | 47.2 | 32.78 | 68.82 |
| mecC (methicillin) | 0.2 | 0.42 | 0.11 | 0.0 | 0.0 |
| PVL [lukS-PV/lukF-PV] | 22.85 | 1.6 | 2.69 | 6.83 | 9.47 |
| scn (IEC / staph. complement inhibitor) | 87.52 | 14.21 | 64.3 | 44.63 | 71.59 |
| sak (IEC / staphylokinase) | 79.52 | 12.68 | 50.97 | 38.07 | 63.28 |
| tst (TSST-1) | 8.53 | 4.84 | 15.38 | 8.04 | 10.16 |
| blaZ (penicillinase) | 73.7 | 51.66 | 74.09 | 57.74 | 71.13 |
| any acquired tetracycline [tet(K/L/M/O)] | 17.9 | 37.85 | 34.62 | 31.18 | 44.57 |
| any MLSb [erm(*)] | 33.93 | 30.19 | 38.92 | 26.01 | 50.12 |
| aac(6')-Ie/aph(2'')-Ia (gentamicin) | 11.25 | 19.25 | 34.62 | 11.57 | 35.33 |

## Strongest source differences (top 15 by q-value)

| gene | chi2 | p | q(BH) | min exp |
|---|---|---|---|---|
| mecA (methicillin) | 2123.909 | 0.00e+00 | 0.00e+00 | 170.77 |
| tet(L) | 5559.049 | 0.00e+00 | 0.00e+00 | 9.08 |
| sak (IEC / staphylokinase) | 10462.613 | 0.00e+00 | 0.00e+00 | 110.98 |
| scn (IEC / staph. complement inhibitor) | 15808.094 | 0.00e+00 | 0.00e+00 | 77.74 |
| PVL [lukS-PV/lukF-PV] | 1489.888 | 0.00e+00 | 0.00e+00 | 90.73 |
| tet(M) | 1422.669 | 8.40e-307 | 4.62e-306 | 38.92 |
| any acquired tetracycline [tet(K/L/M/O)] | 1420.727 | 2.21e-306 | 1.04e-305 | 85.29 |
| mph(C) | 1372.697 | 5.75e-296 | 2.37e-295 | 74.09 |
| msr(A) | 1325.777 | 8.58e-286 | 3.15e-285 | 76.76 |
| blaZ (penicillinase) | 1103.356 | 1.42e-237 | 4.68e-237 | 120.87 |
| dfrG (trimethoprim) | 957.987 | 4.54e-206 | 1.36e-205 | 41.23 |
| aac(6')-Ie/aph(2'')-Ia (gentamicin) | 899.21 | 2.47e-193 | 6.79e-193 | 52.54 |
| erm(A) | 839.975 | 1.68e-180 | 4.27e-180 | 64.42 |
| erm(B) | 703.695 | 5.52e-151 | 1.30e-150 | 13.87 |
| any aminoglycoside [aac/ant/aph] | 604.247 | 1.87e-129 | 4.11e-129 | 180.22 |

## Gene burden (median [IQR]) by source

| metric | human_or_clinical | livestock | companion_animal | food | environment | Kruskal p |
|---|---|---|---|---|---|---|
| number_amr_genes | 11 [8-14] | 6 [3-15] | 10 [7-15] | 8 [5-12] | 13 [8-19] | 1.06e-269 |
| number_virulence_genes | 16 [14-19] | 13 [10-14] | 14 [10-18] | 13 [11-15] | 16 [13-19] | 0.00e+00 |

## Guardrails

- Prevalence reflects **deposited genomes**, not population prevalence; convenience/deposition bias applies, especially to small non-human arms.
- Odds ratios contrast source categories within the public repository; they are associations in deposited data, not transmission or causal effects.
- Small denominators (environment n<500, companion animal n<1000) widen CIs; interpret those ORs cautiously.
