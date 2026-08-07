# Technical Appendix

### Global One Health MRSA genomic sequence visibility at the human–animal–environment interface, 2015–2025

Hasnain N, Mannan MA, Hossain S

This appendix supports the Methods and Results of the main manuscript. It (i) documents the analysis pipeline and every classification and statistical decision in reproducible detail; (ii) reports supporting analyses that substantiate the study's stated limitations — geographic concentration, single-country dominance of the small non-human arms, regional visibility, temporal deposition dynamics, and source-classifier validation; and (iii) documents the supplementary tables and figures, which are supplied as separate machine-readable files in this Zenodo archive. All results derive from a single frozen data snapshot and regenerate deterministically from the archived code (Section A). Throughout, we interpret every quantity as a measure of public genomic **visibility** — the intensity of sequencing and deposition — and never as MRSA incidence, prevalence, disease burden, or transmission.

---

## A. Reproducibility and data provenance

**Data snapshot.** NCBI Pathogen Detection *Staphylococcus aureus* AMR metadata, snapshot **PDG000000073.1234**. The metadata table (`.amr.metadata.tsv`) and the SNP-cluster membership table (`.reference_target.all_isolates.tsv`) were retrieved once from the NCBI Pathogen Detection FTP tree and frozen; the exact source URL is stored in `data/raw/latest_ncbi_metadata_url.txt`. World Bank Open Data (income classification, total population, GDP per capita, current health expenditure per capita, physicians per 1,000, hospital beds per 1,000; and, for this appendix, region assignment) were retrieved from the World Bank API v2 and frozen locally using the most recent non-missing value per country over 2015–2025; the frozen files and hashes are included in this archive.

**Primary cohort.** Records with a parsed collection year in 2015–2025 (**n = 86,947**). A deposition-year (target-creation) window and a contemporary collection-year 2020–2025 subset were retained for sensitivity analyses (Section C4).

**Pipeline.** Analysis is a deterministic seven-stage Python pipeline; each stage writes versioned tables and a narrative report:

| Script | Role | Key outputs |
|---|---|---|
| `01_feasibility_ncbi_s_aureus.py` | Ingest snapshot; parse dates, country, source; call MRSA markers; join World Bank; build frozen row table | `data/processed/…_feasibility_rows.csv.gz`; tables 01–12 |
| `02_feasibility_mixed_source_clusters.py` | Join SNP-cluster (PDS) membership; enumerate mixed-source clusters | tables 13–16 |
| `03_window_choice_2020_2025_vs_2015_2025.py` | Collection- vs deposition-year window comparison | tables 18–25 |
| `04_source_stats.py` | Per-source gene prevalence, omnibus χ²/FDR, OR vs human, gene burden | tables 26–29 |
| `05_inequity.py` | Visibility, coverage, metadata, lag, negative-binomial model, gene-by-income | tables 30–37 |
| `06_figures.py` | Main display items (Figures 1–5) | `figures/fig1–5` |
| `07_appendix.py` | **This appendix**: classifier QA, concentration, regional visibility, window sensitivity, master gene table, Figures S1–S3 | tables 38–42; `figures/figS1–S3` |

**Software.** Python 3; pandas, NumPy, SciPy (`chi2_contingency`, `fisher_exact`, `kruskal`), statsmodels (`GLM` Poisson/Negative-Binomial; `proportion_confint`), matplotlib (Agg backend). The pipeline is single-threaded and free of random seeds; every table and figure is a deterministic function of the frozen snapshot and the World Bank pull. Figures are rendered at 300 dpi in a colour-blind-safe Okabe-Ito palette and exported as both `.png` and `.pdf`.

**Data availability.** All code and derived tables/figures are included in the accompanying Zenodo deposit; the DOI will be inserted after Zenodo creates the record. No individual-level or restricted data are used; the study reanalyses openly deposited genome metadata and public World Bank indicators.

---

## B. Supplementary Methods

### B1. Source-classification algorithm

Each isolate was assigned to exactly one of six categories — **human/clinical, livestock, companion animal, food, environment, unknown** — by matching a normalised concatenation of the `host`, `isolation_source` and `source_type` fields against curated term dictionaries. Two design choices protect against classic false matches:

1. **Tokenisation.** Single-word terms are matched on whole word tokens (`[a-z]+`), not substrings, so "cat" is not matched inside "cattle" and "pig" is not matched inside "pigeon"; multi-word scientific and product terms (e.g. *bos taurus*, "chicken meat", "farm environment") are matched as phrases.
2. **Fixed priority order.** When a record matches more than one dictionary, categories are resolved in the order **livestock → companion animal → environment → food → human/clinical**. This deliberately keeps, for example, a "swine farm" record in livestock rather than environment, and a "retail chicken meat" record in food only if no animal-host token is present.

Dictionaries (verbatim):

- **Livestock:** *bos taurus, sus scrofa, gallus gallus, ovis aries, capra hircus, bubalus bubalis*; cattle, cow, bovine, calf, swine, pig, porcine, pork, chicken, poultry, broiler, layer, turkey, duck, goose, sheep, ovine, goat, caprine, buffalo, camel, llama, alpaca.
- **Companion animal:** *felis catus, canis lupus familiaris*; dog(s), cat(s), feline, canine, horse, equine, pet, companion.
- **Environment:** wastewater, sewage, effluent, environment(al), water, soil, surface, air, dust, sediment, sludge, "farm environment".
- **Food:** food, meat, retail, milk, dairy, cheese, beef, pork, carcass, seafood, fish, egg(s), "chicken meat", "turkey meat", "ready to eat".
- **Human/clinical:** *homo sapiens*; human, patient, clinical, hospital, blood, wound, nasal, urine, sputum, respiratory, abscess, skin, infection, screening.

Records with no informative host/source/type field, or with text matching no dictionary, were retained as **unknown** and excluded from all source-comparison hypothesis tests. Classifier validation is reported in Section C1 (Table S4).

### B2. AMR and virulence gene panel

Genotype strings were taken from the AMRFinderPlus fields already computed and curated by NCBI Pathogen Detection (`AMR_genotypes`, `AMR_genotypes_core`, `virulence_genotypes`). Tokens were normalised to lower case and stripped of AMRFinderPlus qualifiers (`=PARTIAL`, `=MISTRANSLATION`, `=POINT`, `=PARTIAL_END_OF_CONTIG`) before membership testing. The locked 33-marker panel:

- **Methicillin:** *mecA*, and the divergent homologue *mecC* (pre-computed boolean columns).
- **β-lactam / other AMR:** *blaZ*; acquired tetracycline *tet(K)/tet(L)/tet(M)/tet(O)* — the intrinsic chromosomal *tet(38)* (present in ≈99.6% of genomes) is **deliberately excluded** from "acquired tetracycline"; MLSB *erm(A)/erm(B)/erm(C)*, *msr(A)*, *mph(C)*; aminoglycoside *aac(6′)-Ie/aph(2″)-Ia* and any *aac/ant/aph*; trimethoprim *dfrG*, *dfrS1*; fusidic acid *fusB*, *fusC*; fluoroquinolone-associated *gyrA/parC* QRDR point mutations.
- **Virulence:** Panton-Valentine leukocidin (*lukS-PV* + *lukF-PV*), immune-evasion cluster *scn* and *sak*, toxic-shock toxin *tst*, enterotoxins *sea–see*, exfoliative toxins *eta/etb*, collagen adhesin *cna*.

The immune-evasion chemotaxis inhibitor *chp* returned zero tokens in this snapshot (not emitted by AMRFinderPlus here) and is therefore not measurable; it is excluded rather than reported as absent.

### B3. Country covariates

`geo_loc_name` was parsed to a country string (taking the text before the first colon) and mapped to World Bank names/ISO3 through an alias table (e.g. USA/UK constituent nations, Korea Rep., Türkiye, Czechia, Congo Dem./Rep.). Country records were joined to World Bank income group, total population and the four health-system indicators, using the most recent non-missing value per country.

### B4. Statistical procedures

- **Visibility** = deposited genomes per million population, by income group, World Bank region and country; and **country coverage** = share of countries in a stratum with ≥1 deposited genome, denominator = all World Bank countries in that stratum.
- **Proportions** are reported with **Wilson** 95% CIs.
- **Per-source gene prevalence** was compared by **Pearson χ²** across the five informative source categories, with **Benjamini-Hochberg** FDR correction across the 33 genes; the minimum expected cell count is reported per test.
- **Pairwise odds ratios versus the human reference** used **Fisher exact** p-values with **Wald** 95% CIs and a **0.5 Haldane–Anscombe** correction for zero cells.
- **Per-country deposition counts** were modelled by **negative-binomial** regression on income group with a **log-population offset**; the dispersion parameter α was estimated by the **Cameron–Trivedi** auxiliary regression and fixed for the final fit to yield proper Wald CIs.
- **Gene burden and submission lag** across groups were compared by **Kruskal–Wallis**.
- **Concentration** across countries was summarised by the **Gini coefficient** and the **Lorenz curve** (this appendix, Section C2), and by cumulative top-k country shares.

---

## C. Supplementary results (supporting analyses)

### C1. Source-classification quality control and validation (Table S4)

The classifier's behaviour is transparent and its residual is small. Every human, livestock, companion-animal, food and environmental isolate carried at least one informative host/source field (100% in each), so no positive assignment rests on empty metadata. Crucially, of the 13,914 **unknown** genomes, **12,230 (88%)** carried **no** host, isolation-source or source-type field at all — "unknown" therefore reflects **missing metadata, not classifier failure**; only 1,684 genomes (1.9% of the cohort) had descriptive text that matched no dictionary.

**Table S4. Source-classification basis and residual.**

| Source category | n | % of cohort | % with host | % with isolation source | % with any informative field | n with no field |
|---|---|---|---|---|---|---|
| Human/clinical | 65,850 | 75.7 | 99.5 | 83.9 | 100.0 | 0 |
| Livestock | 4,005 | 4.6 | 88.7 | 75.7 | 100.0 | 0 |
| Companion animal | 930 | 1.1 | 95.9 | 68.5 | 100.0 | 0 |
| Food | 1,815 | 2.1 | 9.1 | 97.3 | 100.0 | 0 |
| Environment | 433 | 0.5 | 24.7 | 99.1 | 100.0 | 0 |
| Unknown | 13,914 | 16.0 | 8.4 | 7.2 | 12.1 | 12,230 |

The classification is further validated **biologically** by the results themselves: the human-specific immune-evasion-cluster genes *scn* and *sak* and Panton-Valentine leukocidin form a steep, monotonic human→companion→food→livestock gradient (main-text Figure 3; Section C5), exactly as expected from φ3-prophage host adaptation. A source labelling that were noise could not reproduce this textbook signal. Exemplar raw `host`/`isolation_source` values per category are tabulated in `tables/38b_source_classification_exemplars.csv`.

### C2. Geographic concentration of genomic visibility (Figure S1, Table S5)

Deposited genomes are extraordinarily concentrated: across the 115 countries with ≥1 genome the **Gini coefficient is 0.885**, the single top country (United States) contributes **29.8%**, and the **top 10 countries account for 81.3%** of all genomes (Figure S1A; Table S5). The MRSA-marker-positive subset is at least as concentrated (Gini 0.886; top 10 = 86.4%).

This concentration is most consequential inside the small non-human arms, each of which is dominated by a single country: **China contributes 52.4% of all environmental genomes, 37.9% of food and 31.3% of livestock; Japan contributes 43.0% of companion-animal genomes** (Figure S1B; Table S5). The top three countries supply 58–74% of each non-human arm. This is the quantitative basis for the manuscript's caution that the small non-human and low-income arms are susceptible to single-project/single-country dominance, and that their gene frequencies partly reflect one national sampling programme rather than a global reservoir.

![Figure S1. Geographic concentration of public *S. aureus*/MRSA genomic visibility.](figures/figS1_geographic_concentration.png){width=6.5in}

**Table S5. Geographic concentration and single-country dominance.**

| Stratum | Genomes | Countries | Gini | Top country | Top-1 % | Top-3 % | Top-10 % |
|---|---|---|---|---|---|---|---|
| All sources | 86,219 | 115 | 0.885 | United States | 29.8 | 51.4 | 81.3 |
| MRSA-marker + | 47,553 | 96 | 0.886 | United States | 32.5 | 58.4 | 86.4 |
| Livestock | 4,005 | 49 | 0.756 | China | 31.3 | 58.5 | 78.6 |
| Companion animal | 930 | 20 | 0.748 | Japan | 43.0 | 72.7 | 95.2 |
| Food | 1,815 | 39 | 0.793 | China | 37.9 | 63.9 | 85.2 |
| Environment | 433 | 24 | 0.781 | China | 52.4 | 74.1 | 91.5 |

*(Concentration is computed over the 86,219 genomes with a parsed country; 728 records lacked a resolvable `geo_loc_name`.)*

### C3. Regional visibility (Figure S3, Table S6)

Normalising to population by World Bank region reproduces and sharpens the income gradient. Visibility spans **~250-fold**, from **68.9 genomes per million in North America** and 36.7 in Europe & Central Asia down to **1.22 in Sub-Saharan Africa and 0.28 in South Asia** (Figure S3; Table S6). South Asia — home to roughly a quarter of the world's population — contributes just **457 genomes (0.5% of the cohort)** from 4 of 6 countries, and Sub-Saharan Africa 1,535 (18 of 48 countries). Non-human One Health visibility per million is an order of magnitude lower again in these regions. Regional coverage is as informative as the rate: only 10 of 42 Latin-American & Caribbean countries and 18 of 48 Sub-Saharan African countries have deposited any genome.

![Figure S3. Public genomic visibility by World Bank region.](figures/figS3_regional_visibility.png){width=6.5in}

**Table S6. Public genomic visibility by World Bank region.**

| Region | Records | % of cohort | Genomes / million | Non-human / million | Countries with genome / total (coverage) |
|---|---|---|---|---|---|
| North America | 25,980 | 29.9 | 68.94 | 1.640 | 2 / 3 (66.7%) |
| Europe & Central Asia | 33,921 | 39.0 | 36.68 | 2.066 | 36 / 58 (62.1%) |
| East Asia & Pacific | 18,699 | 21.5 | 7.92 | 1.621 | 12 / 37 (32.4%) |
| Latin America & Caribbean | 2,195 | 2.5 | 3.34 | 0.275 | 10 / 42 (23.8%) |
| Middle East, N. Africa, Afg. & Pak. | 2,213 | 2.5 | 2.77 | 0.213 | 15 / 23 (65.2%) |
| Sub-Saharan Africa | 1,535 | 1.8 | 1.22 | 0.203 | 18 / 48 (37.5%) |
| South Asia | 457 | 0.5 | 0.28 | 0.059 | 4 / 6 (66.7%) |

### C4. Temporal deposition dynamics and window choice (Figure S2, Table S7)

Two date fields — collection year and target-creation (deposition) year — tell different stories, and the difference justifies the study design. Within the primary cohort, genomes by **collection** year rise to a 2019 peak and then **fall** after ~2022, whereas genomes by **deposition** year keep rising into 2023 (Figure S2A). The apparent recent decline in collection-year counts is a **submission-lag artefact** — infections from 2023–2025 are largely not yet deposited — not a real fall in occurrence, which is why deposition lag is reported explicitly (main-text Table 1) and why disease trends are not inferred. The non-human share of genomes is broadly stable at 5–13% across collection years (Figure S2B; the 2025 uptick reflects the near-absence of recent human deposits and small denominators).

![Figure S2. Temporal deposition dynamics in the primary cohort.](figures/figS2_temporal_dynamics.png){width=6.5in}

The headline composition is robust to the window and date definition chosen (Table S7): whether one uses collection or deposition year, and 2015–2025 or 2020–2025, the MRSA-marker fraction stays 48–55% and the non-human share 7.4–8.4%. We adopt collection-year 2015–2025 as the primary denominator for epidemiological interpretability, and report the contemporary collection-year 2020–2025 subset as a sensitivity analysis (Table S2).

**Table S7. Window and date-definition sensitivity.**

| Window | Date field | Records | MRSA-marker % | Non-human % | Countries |
|---|---|---|---|---|---|
| 2020–2025 | Deposition | 146,304 | 54.5 | 7.4 | 141 |
| 2020–2025 | Collection | 39,840 | 47.7 | 8.4 | 91 |
| 2015–2025 | Deposition | 155,137 | 55.3 | 7.5 | 149 |
| **2015–2025** | **Collection (primary)** | **86,947** | **55.4** | **8.3** | **117** |

*The "Countries" column counts distinct parsed `geo_loc_name` strings in each window; for the primary cohort this is 117, of which 115 resolve to a World Bank–matched country and enter the geographic-concentration analysis (Section C2, Table S5) and 97 additionally carry a World Bank income-group classification (main-text Table 1). Records with no resolvable country (n = 728 in the primary cohort) are excluded from all country-level analyses.*

### C5. Full gene prevalence and host differentiation (Tables S1, S8)

Gene prevalence differed significantly by source for **32 of 33 markers** (omnibus χ², FDR q < 0.05); only enterotoxin *see* (present in <0.02% of genomes) showed no source difference. The complete per-source prevalence table with Wilson CIs, omnibus χ² and BH q-values is archived as Table S1 in `tables/42_gene_master_by_source.csv` (with source-detail files `tables/26_gene_prevalence_by_source.csv` and `tables/27_gene_omnibus_tests.csv`); all pairwise odds ratios versus human are in `tables/28_gene_or_vs_human.csv`. The two opposing axes are captured by the livestock-versus-human odds ratios (main-text Figure 4):

| Gene | Livestock vs human OR (95% CI) | Direction |
|---|---|---|
| *scn* (immune-evasion) | 0.024 (0.022–0.026) | depleted in livestock |
| *sak* (immune-evasion) | 0.037 (0.034–0.041) | depleted in livestock |
| PVL (*lukS/F-PV*) | 0.055 (0.043–0.070) | depleted in livestock |
| *dfrG* (trimethoprim) | 2.94 (2.72–3.19) | enriched in livestock |
| *tet(M)* | 2.91 (2.68–3.16) | enriched in livestock |
| Acquired tetracycline | 2.79 (2.61–2.99) | enriched in livestock |
| *cna* (collagen adhesin) | 2.01 (1.86–2.17) | enriched in livestock |

**Gene burden (Table S8).** Acquired-AMR gene burden was lowest in livestock (median 6, IQR 3–15) and highest in environmental isolates (13, 8–19); virulence-gene burden was highest in human isolates (16, 14–19) (Kruskal–Wallis p ≪ 10⁻²⁰⁰ for both; `tables/29_gene_burden_by_source.csv`). Because gene detection depends on assembly quality, which correlates with sequencing capacity, these burdens are interpreted as properties of deposited genomes, not of the underlying bacterial populations.

### C6. Income-stratified robustness (Tables S2, S3, S9–S11)

The inequity is stable across specifications. The negative-binomial model (`tables/34`) gives per-country deposition incidence-rate ratios versus high-income of **0.068 (0.045–0.103)** for lower-middle and **0.080 (0.047–0.135)** for low-income countries, with a modest, non-significant upper-middle contrast (0.715; 0.508–1.006) — i.e. the pooled upper-middle gap is driven by large under-depositing populations, not typical-country behaviour. The contemporary 2020–2025 subset preserves the gradient (**Table S2**, high 24.3 vs low 0.32 genomes/million; `tables/37`). The metadata paradox is quantified in Tables S10–S11: host metadata completeness *rises* as income falls (81.6% high → 99.7% low) and median submission lag *shortens* (4 → 2 years), consistent with a shift from bulk automated surveillance deposition toward smaller, richly curated research submissions. Income-stratified human gene profiles (**Table S3**; `tables/35`) show acquired tetracycline rising as income falls (16.2% → 32.7%) and PVL peaking in lower-middle-income isolates (46.7%), reflecting lineage composition and interpreted descriptively.

---

## D. Supplementary figure legends

**Figure S1. Geographic concentration of public *S. aureus*/MRSA genomic visibility.** (A) Lorenz curve of deposited genomes across the 115 countries with ≥1 genome; the dashed diagonal is perfect equality, and the shaded gap indexes concentration (Gini = 0.885; top 10 countries = 81%). (B) Share of genomes contributed by the single top country within each One Health arm, showing single-country dominance (China of environment/food/livestock; Japan of companion animals); arm size and country count annotated.

**Figure S2. Temporal deposition dynamics, primary cohort (collection-year 2015–2025).** (A) Deposited genomes by collection year (blue) versus deposition/target-creation year (orange); the divergence after ~2022 is submission lag, not declining occurrence. (B) Non-human One Health share of genomes by collection year.

**Figure S3. Public genomic visibility by World Bank region.** Deposited genomes per million population (log scale) by region; the two most under-represented regions (Sub-Saharan Africa, South Asia) are highlighted. Genome count and country coverage annotated.

---
---
---

## E. Supplementary data files and Zenodo archive

The Technical Appendix is intentionally kept as a readable methods-and-interpretation document. The complete supplementary tables are submitted as separate, machine-readable CSV files in the accompanying Zenodo archive under `tables/`. This avoids an unreadably wide appendix while preserving every row, denominator, confidence interval, test statistic, model estimate, country record, and classifier audit exemplar.

### E1. Supplementary table crosswalk

| Table | Contents | Complete file(s) in this archive |
|---|---|---|
| Table S1 | Full 33-gene prevalence by source, Wilson CI, omnibus χ²/BH q, and pairwise ORs versus human | `tables/42_gene_master_by_source.csv`; `tables/26_gene_prevalence_by_source.csv`; `tables/27_gene_omnibus_tests.csv`; `tables/28_gene_or_vs_human.csv` |
| Table S2 | Contemporary collection-year 2020–2025 visibility by income group | `tables/37_visibility_by_income_2020_2025.csv` |
| Table S3 | AMR/virulence gene prevalence by income group among human isolates | `tables/35_amr_by_income_human.csv` |
| Table S4 | Source-classification QC and raw metadata exemplars | `tables/38_source_classification_qa.csv`; `tables/38b_source_classification_exemplars.csv` |
| Table S5 | Geographic concentration and single-country dominance | `tables/39_geographic_concentration.csv` |
| Table S6 | Regional genomic visibility and country coverage | `tables/40_regional_visibility.csv` |
| Table S7 | Collection/deposition window sensitivity | `tables/41_window_sensitivity.csv`; `tables/18_window_comparison_2020_2025_vs_2015_2025.csv` |
| Table S8 | AMR and virulence gene burden by source | `tables/29_gene_burden_by_source.csv` |
| Table S9 | Country-level counts, population, covariates, and visibility rates | `tables/33_country_level_dataset.csv` |
| Table S10 | Negative-binomial visibility model | `tables/34_negbin_visibility_by_income.csv` |
| Table S11 | Metadata completeness and submission lag by income group | `tables/31_metadata_completeness_by_income.csv`; `tables/32_submission_lag_by_income.csv` |

The remaining numbered files in `tables/` are upstream derived outputs that document the main-text tables, cluster analysis, temporal trends, source summaries, and exploratory metadata counts. They are retained for provenance and regeneration rather than discarded.

### E2. Supplementary figures

Supplementary Figures S1–S3 are supplied in both 300-dpi PNG and vector PDF formats under `figures/`. Main-text Figures 1–5 are also supplied in both formats. Figure captions and interpretation guardrails are given above and in the main manuscript.

### E3. Reproducibility route

The reproducible route is:

1. `python scripts/00_freeze_raw_sources.py --verify` — verify the pinned public NCBI and World Bank inputs against `data/raw/frozen_sources_manifest.json`.
2. Run scripts `01` through `05` in numerical order to regenerate the primary derived data and statistical tables. The scripts prefer the frozen local inputs when present; they fall back to live public endpoints only when the frozen inputs are absent.
3. Run `python scripts/06_figures.py` to regenerate main-text figures.
4. Run `python scripts/07_appendix.py` to regenerate Tables 38–42, the supplementary figures, and the supporting report.

All quantities describe public genomic visibility and deposition intensity. They are not estimates of MRSA incidence, prevalence, disease burden, or direct transmission.

### E4. Archive boundary

The NCBI Pathogen Detection and World Bank inputs included here are public sources frozen for reproducibility. The archive contains no patient-level identifiers, restricted clinical data, or private laboratory records. The NCBI snapshot and World Bank values may be revised by their providers after this archive is created; the frozen files and checksums define the exact analysis version used for the manuscript.

---

## F. Interpretation guardrails

1. **Visibility is not burden.** Every per-million rate measures public sequencing/deposition intensity, not MRSA incidence, prevalence or mortality. Income and regional gradients in visibility must not be read as gradients in disease.
2. **The repository is a convenience sample.** Geographic concentration (Gini 0.885) and single-country dominance of the non-human arms mean small strata can be driven by one national programme; their gene frequencies and cluster memberships carry wide, partly unquantifiable uncertainty.
3. **Gene content is confounded by lineage and assembly quality.** Cross-source and cross-income gene contrasts reflect which clones were sequenced and how well they assembled; MLST/clonal-complex assignments were absent from this snapshot (`computed_types` empty), so lineage was not directly adjusted. These contrasts are descriptive.
4. **Mixed-source clusters are genomic-overlap signals, not transmission.** Absent sampling frames, epidemiological linkage and directionality, shared SNP clusters indicate compatibility with cross-interface circulation only.
5. **Recent years are lag-dominated.** Collection-year counts after ~2022 undercount true activity because deposition lags collection; contemporary trends are read on deposition year or the 2020–2025 subset with this caveat.
