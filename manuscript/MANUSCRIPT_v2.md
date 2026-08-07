# Global One Health MRSA genomic sequence visibility at the human–animal–environment interface, 2015–2025

**Running title:** Global One Health MRSA genomic visibility

**Authors:** Nehal Hasnain,¹ Muhammad Abdul Mannan,¹ Sajjad Hossain²

¹ Department of Microbiology and Parasitology, Faculty of Animal Science & Veterinary Medicine, Sher-e-Bangla Agricultural University (SAU), Dhaka, Bangladesh

² Faculty of Animal Science & Veterinary Medicine, Sher-e-Bangla Agricultural University (SAU), Dhaka, Bangladesh

**Corresponding author:** Nehal Hasnain, Department of Microbiology and Parasitology, Faculty of Animal Science & Veterinary Medicine, Sher-e-Bangla Agricultural University (SAU), Sher-e-Bangla Nagar, Dhaka 1207, Bangladesh; email nehal.has9@gmail.com

**ORCID:** Nehal Hasnain 0009-0009-2648-1475; Muhammad Abdul Mannan 0000-0002-9626-2560; Sajjad Hossain (pending)

**Article type:** Research

**Display items:** 5 figures, 1 table (+ Additional file 1, Technical Appendix: supplementary figures S1–S3 and tables S1–S11)

---

## Abstract

**Background:** Public pathogen genomics is increasingly used for antimicrobial resistance (AMR) surveillance. However, deposited sequences may reflect sequencing capacity more than disease occurrence. The balance of human, animal, food, and environmental genomes at the One Health interface is also poorly understood. We measured the global visibility of publicly deposited *Staphylococcus aureus* and meticillin-resistant *S. aureus* (MRSA) genomes across One Health source categories and country income groups from 2015 to 2025. We also compared AMR and virulence gene content and examined genomic overlap between source categories.

**Methods:** We analysed *S. aureus* metadata from NCBI Pathogen Detection (snapshot PDG000000073.1234; collection years 2015–2025; n = 86,947). We classified isolates as human, livestock, companion animal, food, environmental, or unknown using tokenised host and source matching. We obtained MRSA markers (*mecA* and *mecC*), AMR genes, and virulence genes from AMRFinderPlus fields. We linked country records to World Bank income, population, and health-system indicators. We calculated genomes per million population, country coverage, metadata completeness, and submission lag. We also fitted a negative-binomial visibility model, compared gene prevalence by source using χ² tests with Benjamini-Hochberg correction, calculated odds ratios against human isolates, and assessed mixed-source SNP-cluster overlap.

**Results:** High-income countries contributed 81.3% of deposited genomes, and human sources contributed 75.7%. Non-human One Health sources contributed 8.3%. Genome visibility was 51.1 per million population in high-income countries and 0.50 in low-income countries, a ratio of approximately 102:1. Only 36% of low-income countries had at least one deposited genome. Human immune-evasion genes were strongly associated with human isolates. For example, *scn* occurred in 87.5% of human isolates and 14.2% of livestock isolates. The corresponding values for *sak* were 79.5% and 12.7%, and for Panton-Valentine leukocidin they were 22.9% and 1.6%. Acquired tetracycline resistance showed the opposite pattern, occurring in 17.9% of human isolates and 44.6% of environmental isolates. We identified 224 mixed-source SNP clusters.

**Conclusions:** Public MRSA genomic surveillance is highly unequal and mainly represents human sources. The data identify surveillance gaps and show genomic overlap across One Health sources. They should not be interpreted as measures of disease prevalence or proof of direct transmission.

**Keywords:** *Staphylococcus aureus*; MRSA; One Health; genomic surveillance; antimicrobial resistance; health equity; NCBI Pathogen Detection

---

## Background

Whole-genome sequencing is now used for antimicrobial resistance (AMR) surveillance, outbreak investigation, and tracking high-risk lineages <sup><a href="#ref-1">1</a>, <a href="#ref-2">2</a></sup>. Meticillin-resistant *Staphylococcus aureus* (MRSA) causes drug-resistant infections and remains an important pathogen with substantial morbidity and mortality <sup><a href="#ref-3">3</a>–<a href="#ref-4">4</a>–<a href="#ref-5">5</a></sup>. It is also a One Health pathogen. It can circulate among humans, livestock, companion animals, food, and the environment <sup><a href="#ref-5">5</a>, <a href="#ref-6">6</a></sup>. Large public repositories, including NCBI Pathogen Detection, now contain hundreds of thousands of *S. aureus* genomes. NCBI records include AMR and virulence calls generated with AMRFinderPlus <sup><a href="#ref-7">7</a>, <a href="#ref-8">8</a></sup>. Researchers widely use these records for secondary analyses.

Many analyses assume that deposited genomes represent the underlying epidemiology. In practice, sequencing and data sharing depend on laboratory capacity, funding, infrastructure, and national data-sharing policies <sup><a href="#ref-9">9</a>, <a href="#ref-10">10</a></sup>. If genomes come mainly from a small group of countries or sources, researchers may mistake visibility for disease burden. One Health analyses may also rely on reservoirs that have little representation in public databases <sup><a href="#ref-11">11</a>, <a href="#ref-12">12</a></sup>. Previous genomic studies showed that AMR and virulence genes occur in *S. aureus* from several animal hosts <sup><a href="#ref-13">13</a></sup>. However, they did not measure how the public data themselves are distributed across countries, income groups, and animal, food, or environmental reservoirs. They also did not assess companion animals separately from livestock or examine virulence content and metadata quality together.

We measured the global visibility of public *S. aureus* and MRSA genomes at the human–animal–environment interface from 2015 to 2025 <sup><a href="#ref-14">14</a>, <a href="#ref-15">15</a></sup>. We described the distribution of deposited genomes across five One Health source categories and World Bank income groups. We adjusted visibility measures for population and health-system indicators. We also assessed metadata completeness and submission lag, compared AMR and virulence gene content between sources, and examined mixed-source SNP-cluster overlap. We treated the repository as an object of study. We therefore interpreted the results as measures of surveillance visibility, not as measures of MRSA prevalence or transmission.

## Methods

### Data source and study population

We downloaded *S. aureus* metadata from NCBI Pathogen Detection (snapshot PDG000000073.1234). We included records with collection years from 2015 to 2025 (n = 86,947). We retained a deposition-year window and a contemporary 2020–2025 subset for sensitivity analyses. Each record represented one isolate. The metadata included country, host, isolation source, collection date, target-creation date as a proxy for deposition date, and AMRFinderPlus AMR and virulence genotype strings. The records also included gene counts and *mecA* and *mecC* status.

### Source classification

We classified isolates as human or clinical, livestock, companion animal, food, environmental, or unknown. The livestock group included cattle, swine, poultry, sheep, goat, and buffalo. The companion-animal group included dog and cat. The food group included meat, milk, and retail food. The environmental group included wastewater, farm environment, water, soil, and surfaces.

We used tokenised matching of host and isolation-source fields. The matching dictionary included scientific and common names. Tokenisation prevented substring errors. For example, the token "cat" was not matched within "cattle". We retained records without an informative host or source as unknown. We excluded unknown records from source-comparison tests.

### AMR and virulence gene panel

We normalised gene tokens by removing AMRFinderPlus qualifiers (`=PARTIAL`, `=MISTRANSLATION`, `=POINT`, `=PARTIAL_END_OF_CONTIG`). We analysed a prespecified panel of AMR and virulence determinants.

The AMR panel included *mecA*, the divergent *mecC* <sup><a href="#ref-16">16</a></sup>, *blaZ*, and acquired tetracycline genes (*tet(K/L/M/O)*). We excluded the intrinsic chromosomal *tet(38)*. The panel also included macrolide-lincosamide-streptogramin B genes [*erm(A/B/C)*, *msr(A)*, and *mph(C)*], aminoglycoside genes [*aac(6′)-Ie/aph(2″)-Ia* and other *aac/ant/aph* genes], trimethoprim genes (*dfrG* and *dfrS1*), fusidic acid genes (*fusB* and *fusC*), and the fluoroquinolone-associated *gyrA/parC* quinolone resistance-determining region mutations.

The virulence panel included Panton-Valentine leukocidin [*lukS-PV* and *lukF-PV*], immune-evasion cluster genes (*scn* and *sak*), toxic-shock toxin *tst*, enterotoxins *sea*–*see*, exfoliative toxins *eta* and *etb*, and collagen adhesin *cna*. AMRFinderPlus did not emit the chemotaxis-inhibitor gene *chp* in this snapshot. We therefore did not analyse *chp*.

### Country covariates

We parsed country names to ISO3 codes. We then linked the codes to World Bank income classification, total population, and health-system indicators. The health-system indicators were GDP per capita, current health expenditure per capita, physicians per 1,000 population, and hospital beds per 1,000 population. For each country, we used the most recent non-missing value.

### Statistical analysis

We expressed genome visibility as deposited genomes per million population. We calculated this measure for each country and income group. We also calculated country coverage as the proportion of countries in each income group with at least one deposited genome. The denominator included all World Bank countries.

We summarised metadata completeness for host, source, and collection date by income group. We calculated submission lag as the target-creation year minus the collection year. We compared submission lag between income groups with the Kruskal-Wallis test.

We used a negative-binomial regression model to estimate per-country genome counts as a function of income group. The model included log population as an offset. We estimated dispersion using the method of Cameron and Trivedi. We report incidence-rate ratios (IRR) with 95% confidence intervals (CI).

We compared per-source gene prevalence across the five informative source categories with Pearson χ² tests. We applied the Benjamini-Hochberg false-discovery-rate (FDR) correction across genes <sup><a href="#ref-17">17</a></sup>. We calculated pairwise odds ratios (OR) against human isolates using Fisher exact tests. We used Wald 95% CI and a 0.5 Haldane-Anscombe correction for zero cells. We report proportions with Wilson 95% CI.

We assessed mixed-source overlap by counting NCBI SNP clusters, or PDS accessions, that contained isolates from more than one source category. We performed the analyses with Python using pandas, SciPy, and statsmodels. We archived all code and derived tables openly, as described in the Data availability section.

## Results

### Dataset composition

Among 86,947 *S. aureus* genomes collected from 2015 to 2025, 48,161 (55.4%) had a *mecA* or *mecC* MRSA marker. Human or clinical isolates dominated the dataset. They accounted for 65,850 genomes (75.7%). Non-human One Health sources accounted for 8.3% of all genomes. These included livestock (4,005; 4.6%), food (1,815; 2.1%), companion animals (930; 1.1%), and the environment (433; 0.5%). The source was unknown for 13,914 genomes (16.0%) (Figure 2A; Table 1).

![Figure 2. One Health source composition of deposited genomes.](figures/fig2_source_composition.png){width=6.5in}

Every isolate assigned to a non-unknown source category had informative host or isolation-source metadata. This represented 100% of the assigned non-unknown records. Among the unknown records, 88% lacked a host, isolation-source, or source-type field. The unknown group therefore reflected missing metadata more than classifier failure (Supplementary Table S4).

Genomes came from 97 countries with an assigned World Bank income classification. In total, 115 countries had a resolvable location (Supplementary Table S5). High-income countries contributed 81.3% of all genomes. Low-income countries contributed 0.4% (362 genomes).

### Global inequity in genomic visibility

Genome visibility differed by two orders of magnitude between income groups. High-income countries had 51.1 genomes per million population. The values were 4.06 in upper-middle-income countries, 0.54 in lower-middle-income countries, and 0.50 in low-income countries. The high-income to low-income ratio was approximately 102:1 (Figure 1A; Table 1).

![Figure 1. Global inequity in public MRSA and *S. aureus* genomic sequence visibility.](figures/fig1_visibility_inequity.png){width=6.5in}

Non-human One Health visibility showed a similar difference. The rates were 3.00 per million in high-income countries and 0.077 per million in low-income countries, a ratio of approximately 39:1.

Country coverage was incomplete in all income groups. Forty-two of 86 high-income countries (48.8%) and 9 of 25 low-income countries (36.0%) had deposited at least one *S. aureus* genome (Figure 1B). Most low-income countries therefore had no genome in the repository.

The deposited genomes were also concentrated in a small number of countries. The Gini coefficient across the 115 contributing countries was 0.885. The top 10 countries contributed 81.3% of all genomes. Per-capita visibility varied by approximately 250-fold across World Bank regions. It ranged from 68.9 genomes per million population in North America to 0.28 in South Asia (Supplementary Figures S1 and S3; Supplementary Tables S5 and S6).

The negative-binomial model showed lower per-country deposition rates in lower-middle-income countries than in high-income countries (IRR 0.068; 95% CI 0.045–0.103). The rate was also lower in low-income countries (IRR 0.080; 95% CI 0.047–0.135). The upper-middle-income contrast was smaller and was not statistically significant at the country level (IRR 0.715; 95% CI 0.508–1.006). The pooled upper-middle-income gap therefore appears to reflect large populations with low deposition rather than the behaviour of a typical country. The findings were similar in the contemporary 2020–2025 subset (Supplementary Table S2).

### Metadata quality and submission lag

Metadata completeness was highest in the least-represented settings. Host information was available for 99.7% of genomes from low-income countries and 81.6% of genomes from high-income countries. Isolation source was available for 91.7% and 71.3%, respectively.

Median submission lag also differed by income group. It was 4 years in high- and upper-middle-income countries, 3 years in lower-middle-income countries, and 2 years in low-income countries. The Kruskal-Wallis test gave p < 10⁻¹⁰⁰.

These patterns are consistent with a difference in deposition practice. High-income settings may contribute more bulk, automated surveillance deposits. Other settings may contribute smaller research submissions with more detailed curation (Table 1).

The dataset composition remained similar across analysis windows and date definitions. The MRSA-marker fraction ranged from 48% to 55%. The non-human share ranged from 7.4% to 8.4% when we compared collection year and deposition year and used either 2015–2025 or 2020–2025 (Supplementary Table S7). Genome counts by collection year declined after 2022. This pattern reflects submission lag rather than a decline in occurrence, based on comparison of collection-year and deposition-year distributions (Supplementary Figure S2).

### AMR and virulence content across One Health sources

Gene prevalence differed by source for 32 of 33 markers (χ², FDR q < 0.05; Figure 3; Supplementary Table S1).

![Figure 3. AMR and virulence gene prevalence across One Health source categories.](figures/fig3_gene_heatmap.png){width=6.5in}

Two patterns were clear. The first pattern involved human-associated immune-evasion genes. The immune-evasion cluster genes *scn* and *sak* were most common in human isolates. The prevalence of *scn* was 87.5% in human isolates, 64.3% in companion animals, 44.6% in food isolates, 71.6% in environmental isolates, and 14.2% in livestock isolates. The corresponding values for *sak* were 79.5%, 51.0%, 38.1%, 63.3%, and 12.7% <sup><a href="#ref-18">18</a>–<a href="#ref-19">19</a>–<a href="#ref-20">20</a></sup>. Panton-Valentine leukocidin was also more common in human isolates than in livestock isolates (22.9% versus 1.6%). The human-clinical resistance markers *msr(A)* and *mph(C)* showed a similar pattern, with prevalence of approximately 19% in human isolates and approximately 1% in livestock isolates. The same pattern was observed for *fusC*.

Companion-animal isolates had intermediate values for several human-associated genes. This finding is consistent with sharing between humans and pets and with a possible reservoir role for companion animals <sup><a href="#ref-21">21</a>–<a href="#ref-22">22</a>–<a href="#ref-23">23</a></sup>.

The second pattern involved agricultural resistance markers. Acquired tetracycline resistance occurred in 17.9% of human isolates, 37.9% of livestock isolates, and 44.6% of environmental isolates. The genes *tet(L)* and *dfrG* and the collagen adhesin gene *cna* were also more common in animal and environmental isolates.

Odds ratios against human isolates confirmed these differences (Figure 4). Livestock isolates had lower prevalence of *scn*, *sak*, and PVL. They had higher prevalence of acquired tetracycline resistance, *tet(M)*, *dfrG*, and *cna*. Only the enterotoxin gene *see*, which occurred in fewer than 0.02% of genomes, showed no difference between sources.

![Figure 4. Host differentiation of gene carriage, livestock versus human reservoirs.](figures/fig4_host_forest.png){width=6.5in}

AMR gene burden was lowest in livestock isolates, with a median of 6 genes (IQR 3–15). It was highest in environmental isolates, with a median of 13 genes (IQR 8–19). Virulence gene burden was highest in human isolates, with a median of 16 genes (IQR 14–19).

Among human isolates, gene profiles also differed by income setting (Supplementary Table S3). Acquired tetracycline resistance was more common in lower-middle-income isolates than in high-income isolates (32.7% versus 16.2%). PVL was also more common in lower-middle-income isolates (46.7%). This pattern is consistent with the known association between PVL and community-associated *S. aureus* and with reports of high PVL prevalence in The Gambia <sup><a href="#ref-24">24</a>–<a href="#ref-25">25</a>–<a href="#ref-26">26</a></sup>. These income contrasts describe the lineages represented in the deposited genomes. They do not show that resistance or virulence differs because of country income.

### Mixed-source genomic overlap

Among informative clusters, 224 NCBI SNP clusters contained isolates from more than one source category (Figure 5). Human–livestock overlap was most common, with 115 clusters and 2,113 genomes. Human–environment overlap occurred in 47 clusters and included 1,477 genomes. Livestock–food overlap occurred in 29 clusters and included 698 genomes. Livestock–environment overlap occurred in 26 clusters and included 932 genomes. Human–food overlap occurred in 26 clusters and included 850 genomes. Human–companion overlap occurred in 19 clusters and included 983 genomes. Twenty-one clusters included three or more source categories and contained 1,479 genomes.

![Figure 5. Mixed-source SNP-cluster overlap at the One Health interface.](figures/fig5_cluster_overlap.png){width=6.5in}

These clusters show genomic overlap in the deposited data. They are compatible with transmission across the One Health interface, but they are not evidence of direct transmission because the dataset lacks epidemiological links and defined sampling frames.

## Discussion

### Principal findings

We analysed 86,947 public *S. aureus* genomes collected from 2015 to 2025. Three findings were central. First, genomic visibility was highly unequal and focused on human sources. High-income countries had about 100 times more genomes per capita than low-income countries. Most low-income countries had no deposited genome. Animals, food, and environmental sources together made up less than one tenth of the dataset.

Second, the deposited genomes showed consistent differences in gene content between host reservoirs. Third, the repository reflected sequencing capacity more than the underlying epidemiology of MRSA. This limits the conclusions that secondary analyses can support. Public repositories should therefore be studied as data sources with their own sampling structure. They should not be treated as a complete representation of disease.

### Visibility does not measure disease burden

The high-income to low-income difference was approximately 100-fold for genomes per million population. The difference across World Bank regions was approximately 250-fold. The rates ranged from 68.9 genomes per million population in North America to 1.22 in sub-Saharan Africa and 0.28 in South Asia.

This distribution contrasts with the known geographic burden of AMR. Sub-Saharan Africa and South Asia have some of the highest bacterial-AMR mortality rates worldwide <sup><a href="#ref-3">3</a>, <a href="#ref-4">4</a>, <a href="#ref-27">27</a></sup>. *S. aureus* is among the leading pathogens contributing to this burden <sup><a href="#ref-4">4</a></sup>. Genome deposition per million population measures sequencing and data-sharing activity. It does not measure MRSA incidence or prevalence. The income gradient must therefore not be interpreted as a disease gradient.

The comparison still has an important implication. Countries with high AMR burden and limited surveillance systems <sup><a href="#ref-11">11</a>, <a href="#ref-12">12</a></sup> are also the least visible in public genomic databases. Global AMR priorities increasingly use genomic data <sup><a href="#ref-1">1</a>, <a href="#ref-2">2</a></sup>. These data are least complete in the regions that may need the greatest attention. Analyses of public repositories should model visibility bias. Otherwise, pooled analyses may give too little weight to these regions.

### A structural inequity with an actionable diagnosis

This pattern is unlikely to be unique to *S. aureus*. Similar differences affected SARS-CoV-2 genomic surveillance, even during the pandemic when sequencing received exceptional funding and attention <sup><a href="#ref-9">9</a>, <a href="#ref-10">10</a></sup>. A chronic bacterial pathogen without comparable global mobilisation may face an even wider gap.

Our findings also help identify the main constraint. Low-income countries had the highest completeness for host and isolation-source metadata. They also had the shortest submission lag. Host completeness was 99.7% in low-income countries and 81.6% in high-income countries. The median lag was 2 years in low-income countries and 4 years in high-income countries.

These findings suggest that under-representation is not mainly a problem of poor curation or unwillingness to share. When low-income countries contributed genomes, the submissions were smaller but often well annotated. The larger problem appears to be access to sequencing capacity and funding. Investment in sequencing throughput could therefore increase visibility while using existing data-sharing practices.

### The One Health blind spot

The second gap was sectoral. One Health AMR surveillance aims to connect human, animal, food, and environmental data <sup><a href="#ref-1">1</a>, <a href="#ref-14">14</a></sup>. In this dataset, non-human sources accounted for only 8.3% of deposited genomes. Environmental genomes accounted for 0.5%.

A single country also dominated each non-human source group. China provided 31% to 52% of livestock, food, and environmental genomes. Japan provided 43% of companion-animal genomes. The available One Health data therefore reflect a small number of national programmes rather than a global sample.

This creates two risks. Researchers may treat poorly sampled reservoirs as unimportant. They may also draw One Health conclusions from compartments represented by only a few thousand genomes. More sequencing beyond clinical settings is needed. Metagenomic and targeted surveillance of urban wastewater and farm environments can provide population-level AMR signals and complement clinical isolate sequencing <sup><a href="#ref-15">15</a>, <a href="#ref-28">28</a></sup>.

### A defensible biological core at the interface

The deposited genomes showed a clear source-related pattern in gene content. This pattern is the strongest biological result in the study, although it must still be interpreted within the limits of biased sampling.

The first pattern involved human-associated genes. The immune-evasion cluster genes *scn* and *sak* occur on the human-specific β-haemolysin-converting φ3 prophage <sup><a href="#ref-18">18</a>–<a href="#ref-19">19</a>–<a href="#ref-20">20</a></sup>. Both genes were much more common in human isolates than in livestock isolates. The livestock-versus-human odds ratios were 0.024 for *scn*, 0.037 for *sak*, and 0.055 for PVL. Companion-animal isolates had intermediate values, including 64.3% for *scn* and 51.0% for *sak*.

The observed pattern from human to companion animal, food, and livestock is compatible with host-associated differences <sup><a href="#ref-29">29</a></sup>. It is also consistent with sharing of human-associated lineages between people and pets <sup><a href="#ref-21">21</a>–<a href="#ref-22">22</a>–<a href="#ref-23">23</a></sup>. The pattern supports the source classification used in this study. However, it does not prove host adaptation because the dataset lacks complete sampling frames and lineage adjustment.

The second pattern involved agricultural resistance markers. Acquired tetracycline resistance occurred in 17.9% of human isolates, 37.9% of livestock isolates, and 44.6% of environmental isolates. The same source pattern was seen for *tet(M)*, *dfrG*, and *cna*. This pattern is consistent with antimicrobial selection in agricultural settings and with the biology of livestock-associated CC398 <sup><a href="#ref-6">6</a>, <a href="#ref-30">30</a>–<a href="#ref-31">31</a>–<a href="#ref-32">32</a>–<a href="#ref-33">33</a></sup>. MRSA has also been recovered from retail meat and environmental water <sup><a href="#ref-28">28</a>, <a href="#ref-34">34</a>, <a href="#ref-35">35</a>, <a href="#ref-36">36</a></sup>.

The 224 mixed-source SNP clusters were most often shared between human and livestock sources. This finding is compatible with circulation across the interface. It does not establish transmission direction or direct transmission. The dataset does not include the epidemiological links and sampling frames required for that conclusion.

### Implications for reuse of public repositories

These findings have a direct methodological implication for studies that reuse public genome repositories. Previous analyses have compared gene content and pathogenicity between human and animal *S. aureus* isolates <sup><a href="#ref-13">13</a>, <a href="#ref-37">37</a></sup>. Such analyses can be misleading if they assume that deposition reflects the underlying epidemiology.

Visibility bias, dominance by single countries, and lineage confounding are major features of these data. Deposited gene frequencies reflect which clones were sequenced. They do not necessarily reflect which clones circulate in the population. These factors can create false geographic or cross-sector differences if researchers do not account for them.

Our income-stratified results illustrate this issue. Acquired tetracycline resistance and PVL were more common among human isolates from lower-middle-income countries than among those from high-income countries. We report these results descriptively because they reflect the lineages represented in the repository <sup><a href="#ref-24">24</a>–<a href="#ref-25">25</a>–<a href="#ref-26">26</a></sup>.

Within these limits, NCBI Pathogen Detection remains useful. It can show where genomic surveillance exists and where it is absent. It can help prioritise capacity building. It can also generate cautious and testable hypotheses about the human–animal–environment interface. It cannot estimate MRSA prevalence, rank disease burden between countries, or establish transmission.

### Toward equitable One Health genomic surveillance

The global AMR response requires better coverage of these gaps. The findings support three priorities.

First, low- and lower-middle-income countries need greater access to sequencing and bioinformatics capacity. Their metadata were often complete. The main missing resource appears to be sequencing throughput.

Second, national surveillance programmes should integrate human, veterinary, food-chain, and environmental sequencing. Shared metadata standards and timely deposition would make the data more useful <sup><a href="#ref-14">14</a></sup>.

Third, programmes should use scalable methods such as wastewater genomics <sup><a href="#ref-15">15</a>, <a href="#ref-28">28</a></sup>. These approaches can provide useful signals without requiring high-volume clinical isolate sequencing.

Linking these activities with global systems such as GLASS <sup><a href="#ref-38">38</a></sup>, together with fair data-sharing and benefit frameworks, could improve the global record of AMR. Countries that generate surveillance data should also benefit from their use. Future public genomic data should reflect the distribution of AMR more closely, rather than the distribution of sequencing capacity <sup><a href="#ref-27">27</a></sup>.

### Limitations

This study has several limitations. First, visibility is not disease burden. Genome rates per million population measure sequencing and deposition activity. They do not measure MRSA incidence or prevalence. The income gradient in visibility must not be interpreted as a disease gradient.

Second, the repository is a convenience sample. Small non-human and low-income groups may be strongly affected by individual projects. One country contributed 31% to 52% of each non-human One Health group. China dominated the livestock, food, and environmental groups, while Japan dominated the companion-animal group (Supplementary Table S5).

Third, gene detection depends on assembly quality and the analytical pipeline. Some gene differences may therefore reflect technical variation rather than biology. We report these comparisons descriptively <sup><a href="#ref-39">39</a></sup>.

Fourth, lineage confounds income and source comparisons. Country and source gene frequencies depend on which clones were sequenced. Fifth, this snapshot did not include MLST or clonal-complex assignments. We could not adjust directly for lineage. MLST would provide a complementary framework for future analyses <sup><a href="#ref-40">40</a></sup>.

Finally, mixed-source clusters show genomic similarity among deposited isolates. Without defined sampling frames and epidemiological data, they cannot establish the direction or occurrence of transmission.

## Conclusions

Public MRSA and *S. aureus* genomic surveillance is dominated by high-income countries and human clinical genomes. Most low-income countries and most non-human One Health reservoirs are poorly represented or absent. Within these limits, public repositories can help map surveillance gaps, guide capacity building, and generate cautious hypotheses about genomic overlap at the human–animal–environment interface. They cannot be used to infer MRSA prevalence or transmission. Equitable, metadata-complete, and One Health-inclusive genomic surveillance is needed for the global AMR response <sup><a href="#ref-27">27</a>, <a href="#ref-36">36</a>, <a href="#ref-38">38</a></sup>.

---

## List of abbreviations

AMR: antimicrobial resistance; MRSA: meticillin-resistant *Staphylococcus aureus*; WGS: whole-genome sequencing; NCBI: National Center for Biotechnology Information; IEC: immune-evasion cluster; PVL: Panton-Valentine leukocidin; MLS_B: macrolide–lincosamide–streptogramin B; QRDR: quinolone resistance-determining region; MLST: multilocus sequence typing; CC: clonal complex; SNP: single-nucleotide polymorphism; FDR: false discovery rate; OR: odds ratio; IRR: incidence-rate ratio; CI: confidence interval; IQR: interquartile range; GDP: gross domestic product; GLASS: Global Antimicrobial Resistance and Use Surveillance System; LA-MRSA: livestock-associated MRSA; CA-MRSA: community-associated MRSA.

## Declarations

### Ethics approval and consent to participate

Not applicable. This study analysed openly deposited and de-identified bacterial genome metadata from NCBI Pathogen Detection and public country-level indicators from World Bank Open Data. It involved no human participants, human tissue, identifiable personal data, or animal experiments. Ethical approval and consent were therefore not required.

### Consent for publication

Not applicable.

### Availability of data and materials

All data analysed in this study are publicly available. We obtained *Staphylococcus aureus* genome metadata from the [NCBI Pathogen Detection snapshot](https://ftp.ncbi.nlm.nih.gov/pathogen/Results/Staphylococcus_aureus/) PDG000000073.1234. The snapshot aggregates and uniformly re-annotates GenBank/INSDC submissions. We obtained country-level population, income, and health-system indicators from [World Bank Open Data](https://data.worldbank.org). The processing and analysis code, the frozen record of the source data URL, the derived per-country and per-source tables, and the code used to regenerate each figure are archived in the accompanying Zenodo deposit (DOI to be assigned at deposit). These materials are documented in the Technical Appendix (Additional file 1).

### Competing interests

All authors declare no competing interests.

### Funding

No external funding was received.

### Authors' contributions

N.H. designed the study, performed the analyses, and drafted the manuscript. M.A.M. provided critical input on the analysis and supervised the work. S.H. contributed to the One Health framing, bioinformatics support, and data interpretation. All authors read and approved the final manuscript.

### Use of artificial intelligence

The authors used a generative artificial-intelligence assistant, a large language model, to refine language and edit parts of the manuscript. The assistant also helped write and debug the data-analysis and figure-generation code. The authors critically reviewed, validated, and revised all AI-assisted text and code. They take full responsibility for the accuracy and integrity of the manuscript. The tool was not used as a source of scientific data or interpretation. It is not listed as an author, in accordance with ICMJE recommendations and Springer Nature editorial policy.

### Acknowledgements

We thank the researchers and the originating and submitting laboratories worldwide that generated and deposited the *Staphylococcus aureus* genome sequences reanalysed in this study. These sequences were submitted to GenBank/INSDC and aggregated and uniformly processed through NCBI Pathogen Detection. We also acknowledge the World Bank for the open population, income, and health-system indicators used to normalise genomic visibility.

---

## Table 1. Dataset composition and genomic visibility by World Bank income group, 2015–2025

| Income group | Genomes (n) | % of total | MRSA-marker + (%) | Non-human (n) | Countries with genome / total | Genomes per million | Non-human per million | Host metadata (%) | Median submission lag (yr) |
|---|---|---|---|---|---|---|---|---|---|
| High income | 70,700 | 81.3 | 56.1 | 4,155 | 42 / 86 (48.8%) | 51.06 | 3.00 | 81.6 | 4 |
| Upper-middle income | 12,403 | 14.3 | 49.6 | 2,602 | 29 / 59 (49.2%) | 4.06 | 0.85 | 84.4 | 4 |
| Lower-middle income | 1,535 | 1.8 | 46.9 | 246 | 17 / 47 (36.2%) | 0.54 | 0.086 | 92.9 | 3 |
| Low income | 362 | 0.4 | 8.0 | 56 | 9 / 25 (36.0%) | 0.50 | 0.077 | 99.7 | 2 |

*Non-human = livestock, companion-animal, food and environmental sources combined. Records of unknown income (n=1,947) are excluded from this table.*

---

## Figure legends

**Figure 1. Global inequity in public MRSA/*S. aureus* genomic sequence visibility, 2015–2025.** (A) Deposited genomes per million population by World Bank income group, shown on a log scale. (B) Country-level coverage, defined as the percentage of countries in each income group with at least one deposited genome. Counts are shown.

**Figure 2. One Health source composition of deposited genomes.** (A) Number of deposited genomes by source category, shown on a log scale. (B) Non-human One Health genomes per million population by income group, shown on a log scale.

**Figure 3. AMR and virulence gene prevalence across One Health source categories.** The heatmap shows gene prevalence (%) in deposited genomes. AMR determinants are shown in the upper block and virulence determinants in the lower block. The source categories are human, livestock, companion animal, food, and environmental.

**Figure 4. Host differentiation of gene carriage, livestock versus human reservoirs.** Odds ratios on a log scale are shown as points. Bars show 95% confidence intervals for selected genes in livestock relative to human isolates. Blue indicates lower prevalence in livestock. Orange indicates higher prevalence in livestock. The dashed line marks an odds ratio of 1.

**Figure 5. Mixed-source SNP-cluster overlap at the One Health interface.** The figure shows the number of genomes in NCBI SNP clusters that span more than one source category. Overlap indicates genomic similarity or transmission-compatible signals. It does not prove direct transmission.

---

## Supplementary Information

**Additional file 1: Technical Appendix.** The appendix provides analysis provenance, supplementary methods, source-classification dictionaries, the gene panel, and statistical procedures. It also contains supplementary figures S1–S3 and tables S1–S11.

**Supplementary figures**
- **Figure S1.** Geographic concentration of genomic visibility. The figure shows a Lorenz curve across contributing countries (Gini = 0.885) and single-country dominance within each One Health arm.
- **Figure S2.** Temporal deposition dynamics. The figure shows deposited genomes by collection year and deposition year and the non-human share over time.
- **Figure S3.** Public genomic visibility by World Bank region, measured as genomes per million population.

**Supplementary tables**
- **Table S1.** Full 33-gene prevalence by source with Wilson 95% CI, omnibus χ² and Benjamini-Hochberg q-values, and pairwise odds ratios versus human.
- **Table S2.** Contemporary sensitivity analysis for collection years 2020–2025, showing visibility by income group.
- **Table S3.** AMR and virulence gene prevalence by income group among human isolates.
- **Table S4.** Source-classification basis and residual quality-control results.
- **Table S5.** Geographic concentration and single-country dominance by source.
- **Table S6.** Public genomic visibility by World Bank region.
- **Table S7.** Sensitivity analysis for analysis windows and date definitions.
- **Table S8.** AMR and virulence gene burden by source.
- **Table S9.** Country-level dataset with genome counts, population, income, health-system covariates, and visibility rates.
- **Table S10.** Full output from the negative-binomial visibility model.
- **Table S11.** Metadata completeness and submission lag by income group.

## References

<a id="ref-1"></a>1. Djordjevic SP, Jarocki VM, Seemann T, Cummins ML, Watt AE, Drigo B, et al. Genomic surveillance for antimicrobial resistance: a One Health perspective. Nat Rev Genet. 2024;25:142–157. https://doi.org/10.1038/s41576-023-00649-y

<a id="ref-2"></a>2. Wheeler NE, Price V, Cunningham-Oakes E, Tsang KK, Nunn JG, Midega JT, et al. Innovations in genomic antimicrobial resistance surveillance. Lancet Microbe. 2023;4:e1063–e1070. https://doi.org/10.1016/S2666-5247(23)00285-9

<a id="ref-3"></a>3. Murray CJL, Ikuta KS, Sharara F, Swetschinski L, Robles Aguilar G, Gray A, et al. Global burden of bacterial antimicrobial resistance in 2019: a systematic analysis. Lancet. 2022;399:629–655. https://doi.org/10.1016/S0140-6736(21)02724-0

<a id="ref-4"></a>4. Ikuta KS, Swetschinski LR, Robles Aguilar G, Sharara F, Mestrovic T, Gray AP, et al. Global mortality associated with 33 bacterial pathogens in 2019: a systematic analysis for the Global Burden of Disease Study 2019. Lancet. 2022;400:2221–2248. https://doi.org/10.1016/S0140-6736(22)02185-7

<a id="ref-5"></a>5. Turner NA, Sharma-Kuinkel BK, Maskarinec SA, Eichenberger EM, Shah PP, Carugati M, et al. Methicillin-resistant Staphylococcus aureus: an overview of basic and clinical research. Nat Rev Microbiol. 2019;17:203–218. https://doi.org/10.1038/s41579-018-0147-4

<a id="ref-6"></a>6. Park S, Ronholm J. Staphylococcus aureus in Agriculture: Lessons in Evolution from a Multispecies Pathogen. Clin Microbiol Rev. 2021;34. https://doi.org/10.1128/CMR.00182-20

<a id="ref-7"></a>7. Feldgarden M, Brover V, Gonzalez-Escalona N, Frye JG, Haendiges J, Haft DH, et al. AMRFinderPlus and the Reference Gene Catalog facilitate examination of the genomic links among antimicrobial resistance, stress response, and virulence. Sci Rep. 2021;11. https://doi.org/10.1038/s41598-021-91456-0

<a id="ref-8"></a>8. Feldgarden M, Brover V, Fedorov B, Haft DH, Prasad AB, Klimke W. Curation of the AMRFinderPlus databases: applications, functionality and impact. Microb Genom. 2022;8. https://doi.org/10.1099/mgen.0.000832

<a id="ref-9"></a>9. Brito AF, Semenova E, Dudas G, Hassler GW, Kalinich CC, Kraemer MUG, et al. Global disparities in SARS-CoV-2 genomic surveillance. Nat Commun. 2022;13. https://doi.org/10.1038/s41467-022-33713-y

<a id="ref-10"></a>10. Smith EA, Fleming DF, Lackritz EM, Ulrich AK. Inequities and global declines in SARS-CoV-2 genomic data availability hinder response to emerging variants. Npj Viruses. 2026;4. https://doi.org/10.1038/s44298-026-00176-7

<a id="ref-11"></a>11. Gandra S, Alvarez-Uria G, Turner P, Joshi J, Limmathurotsakul D, van Doorn HR. Antimicrobial Resistance Surveillance in Low- and Middle-Income Countries: Progress and Challenges in Eight South Asian and Southeast Asian Countries. Clin Microbiol Rev. 2020;33. https://doi.org/10.1128/CMR.00048-19

<a id="ref-12"></a>12. Iskandar K, Molinier L, Hallit S, Sartelli M, Hardcastle TC, Haque M, et al. Surveillance of antimicrobial resistance in low- and middle-income countries: a scattered picture. Antimicrob Resist Infect Control. 2021;10. https://doi.org/10.1186/s13756-021-00931-w

<a id="ref-13"></a>13. Bruce SA, Smith JT, Mydosh JL, Ball J, Needle DB, Gibson R, et al. Shared antibiotic resistance and virulence genes in Staphylococcus aureus from diverse animal hosts. Sci Rep. 2022;12. https://doi.org/10.1038/s41598-022-08230-z

<a id="ref-14"></a>14. Queenan K, Häsler B, Rushton J. A One Health approach to antimicrobial resistance surveillance: is there a business case for it? Int J Antimicrob Agents. 2016;48:422–427. https://doi.org/10.1016/j.ijantimicag.2016.06.014

<a id="ref-15"></a>15. Hendriksen RS, Munk P, Njage P, van Bunnik B, McNally L, Lukjancenko O, et al. Global monitoring of antimicrobial resistance based on metagenomics analyses of urban sewage. Nat Commun. 2019;10. https://doi.org/10.1038/s41467-019-08853-3

<a id="ref-16"></a>16. García-Álvarez L, Holden MT, Lindsay H, Webb CR, Brown DF, Curran MD, et al. Meticillin-resistant Staphylococcus aureus with a novel mecA homologue in human and bovine populations in the UK and Denmark: a descriptive study. Lancet Infect Dis. 2011;11:595–603. https://doi.org/10.1016/S1473-3099(11)70126-8

<a id="ref-17"></a>17. Benjamini Y, Hochberg Y. Controlling the false discovery rate: a practical and powerful approach to multiple testing. J R Stat Soc Series B Stat Methodol. 1995;57:289–300. https://doi.org/10.1111/j.2517-6161.1995.tb02031.x

<a id="ref-18"></a>18. van Wamel WJB, Rooijakkers SHM, Ruyken M, van Kessel KPM, van Strijp JAG. The innate immune modulators staphylococcal complement inhibitor and chemotaxis inhibitory protein of Staphylococcus aureus are located on beta-hemolysin-converting bacteriophages. J Bacteriol. 2006;188:1310–1315. https://doi.org/10.1128/JB.188.4.1310-1315.2006

<a id="ref-19"></a>19. McCarthy AJ, Lindsay JA. Staphylococcus aureus innate immune evasion is lineage-specific: A bioinfomatics study. Infect Genet Evol. 2013;19:7–14. https://doi.org/10.1016/j.meegid.2013.06.012

<a id="ref-20"></a>20. Cuny C, Abdelbary M, Layer F, Werner G, Witte W. Prevalence of the immune evasion gene cluster in Staphylococcus aureus CC398. Vet Microbiol. 2015;177:219–223. https://doi.org/10.1016/j.vetmic.2015.02.031

<a id="ref-21"></a>21. Bramble M, Morris D, Tolomeo P, Lautenbach E. Potential Role of Pet Animals in Household Transmission of Methicillin-Resistant Staphylococcus aureus: A Narrative Review. Vector Borne Zoonotic Dis. 2011;11:617–620. https://doi.org/10.1089/vbz.2010.0025

<a id="ref-22"></a>22. Scott Weese J. Antimicrobial resistance in companion animals. Anim Health Res Rev. 2008;9:169–176. https://doi.org/10.1017/S1466252308001485

<a id="ref-23"></a>23. Khairullah A, Sudjarwo S, Effendi M, Ramandinianto S, Gelolodo M, Widodo A, et al. Pet animals as reservoirs for spreading methicillin-resistant Staphylococcus aureus to human health. J Adv Vet Anim Res. 2023;10:1. https://doi.org/10.5455/javar.2023.j641

<a id="ref-24"></a>24. DeLeo FR, Otto M, Kreiswirth BN, Chambers HF. Community-associated meticillin-resistant Staphylococcus aureus. Lancet. 2010;375:1557–1568. https://doi.org/10.1016/S0140-6736(09)61999-1

<a id="ref-25"></a>25. Lo WT, Wang CC. Panton-Valentine Leukocidin in the Pathogenesis of Community-associated Methicillin-resistant Staphylococcus aureus Infection. Pediatr Neonatol. 2011;52:59–65. https://doi.org/10.1016/j.pedneo.2011.02.008

<a id="ref-26"></a>26. Darboe S, Dobreniecki S, Jarju S, Jallow M, Mohammed NI, Wathuo M, et al. Prevalence of Panton-Valentine Leukocidin (PVL) and Antimicrobial Resistance in Community-Acquired Clinical Staphylococcus aureus in an Urban Gambian Hospital: A 11-Year Period Retrospective Pilot Study. Front Cell Infect Microbiol. 2019;9. https://doi.org/10.3389/fcimb.2019.00170

<a id="ref-27"></a>27. Sartorius B, Gray AP, Davis Weaver N, Robles Aguilar G, Swetschinski LR, Ikuta KS, et al. The burden of bacterial antimicrobial resistance in the WHO African region in 2019: a cross-country systematic analysis. Lancet Glob Health. 2024;12:e201–e216. https://doi.org/10.1016/S2214-109X(23)00539-9

<a id="ref-28"></a>28. Knight ME, Webster G, Perry WB, Baldwin A, Rushton L, Pass DA, et al. National-scale antimicrobial resistance surveillance in wastewater: A comparative analysis of HT qPCR and metagenomic approaches. Water Res. 2024;262:121989. https://doi.org/10.1016/j.watres.2024.121989

<a id="ref-29"></a>29. Matuszewska M, Murray GGR, Harrison EM, Holmes MA, Weinert LA. The Evolutionary Genomics of Host Specificity in Staphylococcus aureus. Trends Microbiol. 2020;28:465–477. https://doi.org/10.1016/j.tim.2019.12.007

<a id="ref-30"></a>30. Price LB, Stegger M, Hasman H, Aziz M, Larsen J, Andersen PS, et al. Staphylococcus aureus CC398: Host Adaptation and Emergence of Methicillin Resistance in Livestock. mBio. 2012;3. https://doi.org/10.1128/MBIO.00305-11

<a id="ref-31"></a>31. Cuny C, Wieler L, Witte W. Livestock-Associated MRSA: The Impact on Humans. Antibiotics (Basel). 2015;4:521–543. https://doi.org/10.3390/antibiotics4040521

<a id="ref-32"></a>32. Sieber RN, Skov RL, Nielsen J, Schulz J, Price LB, Aarestrup FM, et al. Drivers and Dynamics of Methicillin-Resistant Livestock-Associated Staphylococcus aureus CC398 in Pigs and Humans in Denmark. mBio. 2018;9. https://doi.org/10.1128/MBIO.02142-18

<a id="ref-33"></a>33. Kadlec K, Feßler AT, Hauschild T, Schwarz S. Novel and uncommon antimicrobial resistance genes in livestock-associated methicillin-resistant Staphylococcus aureus. Clin Microbiol Infect. 2012;18:745–755. https://doi.org/10.1111/j.1469-0691.2012.03842.x

<a id="ref-34"></a>34. Tang Y, Larsen J, Kjeldgaard J, Andersen PS, Skov R, Ingmer H. Methicillin-resistant and -susceptible Staphylococcus aureus from retail meat in Denmark. Int J Food Microbiol. 2017;249:72–76. https://doi.org/10.1016/j.ijfoodmicro.2017.03.001

<a id="ref-35"></a>35. González-Machado C, Alonso-Calleja C, Capita R. Prevalence and types of methicillin-resistant Staphylococcus aureus (MRSA) in meat and meat products from retail outlets and in samples of animal origin collected in farms, slaughterhouses and meat processing facilities. A review. Food Microbiol. 2024;123:104580. https://doi.org/10.1016/j.fm.2024.104580

<a id="ref-36"></a>36. Hadjirin NF, Lay EM, Paterson GK, Harrison EM, Peacock SJ, Parkhill J, et al. Detection of livestock-associated meticillin-resistant Staphylococcus aureus CC398 in retail pork, United Kingdom, February 2015. Euro Surveill. 2015;20. https://doi.org/10.2807/1560-7917.ES2015.20.24.21156

<a id="ref-37"></a>37. Randad PR, Dillen CA, Ortines RV, Mohr D, Aziz M, Price LB, et al. Comparison of livestock-associated and community-associated Staphylococcus aureus pathogenicity in a mouse model of skin and soft tissue infection. Sci Rep. 2019;9. https://doi.org/10.1038/s41598-019-42919-y

<a id="ref-38"></a>38. Ajulo S, Awosile B. Global antimicrobial resistance and use surveillance system (GLASS 2022): Investigating the relationship between antimicrobial resistance and antimicrobial consumption data across the participating countries. PLoS One. 2024;19:e0297921. https://doi.org/10.1371/journal.pone.0297921

<a id="ref-39"></a>39. Spies R, Crook DW, Peto TEA, Fowler PW, Turner R, Thai H, et al. Evaluating 12 automated, whole-genome sequencing analysis pipelines for Mycobacterium tuberculosis complex: a comparative study. Lancet Microbe. 2025;6:101210. https://doi.org/10.1016/j.lanmic.2025.101210

<a id="ref-40"></a>40. Liu Y, Ji Y. Multilocus Sequence Typing of Staphylococcus aureus. Methods Mol Biol. 2020:95–102. https://doi.org/10.1007/978-1-4939-9849-4_7
