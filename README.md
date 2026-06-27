# KOA-WM genetic overlap analyses

This repository contains the analysis scripts used for the manuscript:

**Dissecting pleiotropic genetic overlap between knee osteoarthritis liability and white matter microstructure**

The study evaluates the genetic overlap between knee osteoarthritis (KOA) liability and diffusion MRI-derived white matter (WM) microstructure phenotypes. The final manuscript uses genome-wide genetic correlation, local genetic correlation, conditional and conjunctional false discovery rate enrichment, locus-level pruning, functional annotation, ENIGMA external consistency, FinnGen pain triangulation, GWAS Catalog context, exploratory disease-state miRNA analysis, MHC/non-MHC sensitivity, local colocalization, CAUSE, weighted multivariable Mendelian randomization, GSMR2/HEIDI network analysis, OA subtype lookup, BMI/pain/MHC evidence stratification, and WM genetic factor/domain-level convergence.

This repository follows the style of a manuscript-oriented genetics code release. Each script corresponds to one analysis module in the final manuscript. Scripts used only for exploratory testing, manuscript editing, Word formatting, reference-order checks, or supplementary-document rebuilding are not included.

## Repository structure

```text
KOA-WM-Genetics/
├── README.md
├── code/
│   ├── 01. KOA_WM_Global_Genetic_Correlation.R
│   ├── 02. KOA_WM_Local_Genetic_Correlation_LAVA.R
│   ├── 03. KOA_WM_condFDR_conjFDR_Enrichment.py
│   ├── 04. KOA_WM_Locus_Pruning_Merging_and_Directionality.R
│   ├── 05. KOA_WM_Functional_Annotation_and_Enrichment.R
│   ├── 06. KOA_WM_ENIGMA_External_Consistency.R
│   ├── 07. KOA_WM_FinnGen_Pain_Triangulation_and_GWAS_Catalog_Context.R
│   ├── 08. KOA_WM_Clinical_Multiomics_miRNA_Context.R
│   ├── 09. KOA_WM_MHC_nonMHC_Colocalization_CAUSE_Sensitivity.R
│   ├── 10. KOA_WM_MVMR_GSMR2_HEIDI_Network.R
│   ├── 11. KOA_WM_OA_Subtype_BMI_Pain_MHC_Evidence_Integration.R
│   ├── 12. KOA_WM_WM_Genetic_Factor_Domain_Convergence.R
│   └── 13. KOA_WM_Publication_Tables_and_Figures.R
├── data_manifest/
│   ├── input_data_sources.tsv
│   └── wm_trait_metadata.tsv
└── results/
    ├── tables/
    ├── figures/
    └── logs/
```

## Required input data

The repository does not redistribute third-party GWAS, imaging genetics, clinical cohort, LD reference, FUMA, g:Profiler, or GWAS Catalog files. Users should download the required files from the original resources and update `data_manifest/input_data_sources.tsv` before running the scripts.

Required inputs include:

* KOA GWAS summary statistics from the Musculoskeletal Knowledge Portal osteoarthritis genetics resource.
* UK Biobank Brain Imaging Genetics diffusion MRI-derived FA and MD white matter summary statistics.
* LDSC-compatible summary statistics for KOA, WM phenotypes, BMI, depression, FinnGen pain, hip OA, hand OA, and all-OA.
* 1000 Genomes Phase 3 European reference files and HapMap3 SNP list.
* European LD-score reference files.
* LAVA reference files.
* ENIGMA DTI GWAS summary statistics.
* FinnGen R11 PAIN endpoint summary statistics.
* BMI and depression GWAS summary statistics.
* Hip OA, hand OA, and all-OA GWAS summary statistics.
* GWAS Catalog association files.
* Published KOA clinical multi-omics cohort data.
* FUMA SNP2GENE outputs generated from candidate SNPs.
* g:Profiler enrichment outputs generated from candidate gene sets and pleiotropic locus context.
* Pairwise WM-WM LDSC genetic-correlation files used for the Stage6 WM genetic factor layer.

## Local configuration

Scripts use environment variables and do not contain user-specific local paths.

```r
Sys.setenv(
  KOA_WM_PROJECT_DIR = "path/to/KOA-WM-Genetics",
  KOA_WM_DATA_DIR = "path/to/local/input_data",
  KOA_WM_RESULTS_DIR = "path/to/local/results",
  LDSC_REF_DIR = "path/to/eur_w_ld_chr",
  LAVA_REF_DIR = "path/to/LAVA/reference",
  PLINK_EXE = "path/to/plink",
  GCTA_EXE = "path/to/gcta64"
)
```

The scripts write cleaned tables to `results/tables/`, figure source files to `results/figures/`, and logs or software-derived intermediate summaries to `results/logs/`.

## R packages

Core data processing and plotting packages:

```text
data.table, dplyr, tidyr, tibble, readr, stringr, purrr, magrittr,
ggplot2, ggrepel, patchwork, cowplot, scales, RColorBrewer, viridis
```

Genomic interval and annotation packages:

```text
GenomicRanges, IRanges, S4Vectors, GenomeInfoDb, rtracklayer, biomaRt,
org.Hs.eg.db, AnnotationDbi
```

Statistical genetics and causal-sensitivity packages:

```text
Matrix, coloc, cause, TwoSampleMR, MendelianRandomization, ieugwasr,
MRPRESSO, gwasglue
```

Clinical multi-omics and miRNA analysis packages:

```text
limma, edgeR, lme4, lmerTest, broom, broom.mixed, cluster, factoextra
```

Functional enrichment packages:

```text
gprofiler2, clusterProfiler, enrichplot, DOSE
```

Workbook/table export packages used for source-table preparation:

```text
openxlsx, writexl
```

Not every package is called in every script. The package list reflects the complete final analysis workflow and source-table/figure generation used for the manuscript.

## Python packages

```text
python >= 3.8, numpy, scipy, pandas, matplotlib, seaborn
```

## External software and web resources

```text
LDSC, HDL, LAVA, pleioFDR/condFDR/conjFDR, PLINK 1.9, MTAG,
GCTA-GSMR2, FUMA SNP2GENE, g:Profiler, miRDB, GWAS Catalog
```

## Analysis modules

### 01. KOA_WM_Global_Genetic_Correlation.R

**Purpose**

Summarize genome-wide genetic correlation between KOA liability and diffusion MRI-derived WM phenotypes using LDSC and HDL.

**Main analyses**

* Parse LDSC genetic-correlation logs for KOA-WM trait pairs.
* Parse HDL genetic-correlation outputs.
* Combine LDSC and HDL results across FA and MD traits.
* Apply Benjamini-Hochberg FDR correction within method.
* Summarize positive and negative rg distributions and direction-balance indices.
* Generate source data for the global rg panel.

**Main outputs**

* `Supplementary_Table_1_WM_trait_metadata.tsv`
* `Supplementary_Table_2_Global_LDSC_HDL_rg.tsv`
* `Figure_1A_global_rg_source.png`
* `Figure_1A_direction_balance_source.tsv`

**Manuscript link**

Genome-wide genetic correlation; Figure 1A; Supplementary Tables 1-2.

### 02. KOA_WM_Local_Genetic_Correlation_LAVA.R

**Purpose**

Summarize local genetic correlation between KOA and WM phenotypes using LAVA.

**Main analyses**

* Read LAVA local bivariate genetic-correlation output.
* Identify nominally significant local KOA-WM associations.
* Apply FDR correction across local tests.
* Separate FA and MD FDR-corrected local signals.
* Generate local genetic correlation scatter and volcano source plots.

**Main outputs**

* `Supplementary_Table_3_LAVA_nominal_local_rg.tsv`
* `Supplementary_Table_4_LAVA_FA_FDR_significant.tsv`
* `Supplementary_Table_5_LAVA_MD_FDR_significant.tsv`
* `Figure_1B_LAVA_local_rg_scatter_source.png`
* `Figure_1C_LAVA_volcano_source.png`

**Manuscript link**

Local genetic correlation; Figure 1B-C; Supplementary Tables 3-5.

### 03. KOA_WM_condFDR_conjFDR_Enrichment.py

**Purpose**

Run conditional Q-Q, conditional FDR, and conjunctional FDR analyses for KOA-WM phenotype pairs.

**Main analyses**

* Harmonize KOA and WM summary statistics by SNP.
* Exclude the extended MHC interval and 8p23.1 region in the primary enrichment analysis.
* Generate conditional Q-Q plots in both trait directions.
* Estimate condFDR for KOA conditioned on WM and WM conditioned on KOA.
* Estimate conjFDR as the maximum of the two reciprocal condFDR values.
* Export SNP-level enrichment records for downstream pruning and locus consolidation.

**Main outputs**

* `Supplementary_Table_6_7_condFDR_outputs.tsv`
* `conjFDR_SNP_level_enrichment_records.tsv`
* `Figure_2A_F_condFDR_conjFDR_count_source.tsv`
* Conditional Q-Q source figures corresponding to Supplementary Figures 1-13.

**Manuscript link**

Variant-level cross-trait enrichment; Figure 2A-F; Supplementary Figures 1-13; Supplementary Tables 6-7.

### 04. KOA_WM_Locus_Pruning_Merging_and_Directionality.R

**Purpose**

Convert SNP-level condFDR/conjFDR enrichment records into locus-level and tract-level summaries used in the manuscript.

**Main analyses**

* Define independent lead regions from conjFDR-enriched SNP records.
* Merge nearby lead regions into candidate pleiotropic loci.
* Summarize tract-domain counts and hemispheric specificity.
* Evaluate effect-direction concordance for KOA-FA and KOA-MD lead signals.
* Generate source summaries for tract-level bar charts and locus-level Manhattan/LocusZoom follow-up.

**Main outputs**

* `Supplementary_Table_8_tract_category_counts.tsv`
* `Supplementary_Table_9_hemispheric_specificity.tsv`
* `Supplementary_Table_10_top_lead_SNPs.tsv`
* `Supplementary_Table_11_KOA_FA_effect_direction.tsv`
* `Supplementary_Table_12_KOA_MD_effect_direction.tsv`
* `Supplementary_Table_13_distinct_cross_trait_loci.tsv`
* `Figure_2G_I_tract_hemisphere_source.png`
* `Figure_3A_distinct_loci_source.png`

**Manuscript link**

Variant-level cross-trait enrichment and locus consolidation; Figure 2G-I; Figure 3; Supplementary Figures 14-25; Supplementary Tables 8-13.

### 05. KOA_WM_Functional_Annotation_and_Enrichment.R

**Purpose**

Summarize functional annotation and enrichment results for candidate KOA-WM enrichment regions.

**Main analyses**

* Parse FUMA SNP2GENE outputs for candidate and LD-linked SNPs.
* Summarize functional categories, CADD, RegulomeDB, and chromatin-state annotations.
* Summarize positional, eQTL, and chromatin-interaction gene mapping.
* Parse g:Profiler enrichment outputs.
* Generate source data for functional annotation panels.

**Main outputs**

* `Supplementary_Table_14_candidate_shared_regions.tsv`
* `Supplementary_Table_15_candidate_and_LD_linked_SNPs.tsv`
* `Supplementary_Table_16_functional_annotation_summary.tsv`
* `Figure_4A_B_functional_annotation_source.png`
* `Figure_4C_gene_mapping_source.tsv`
* `Figure_4D_gprofiler_enrichment_source.tsv`

**Manuscript link**

Functional annotation and enrichment analysis; Figure 4; Supplementary Tables 14-16.

### 06. KOA_WM_ENIGMA_External_Consistency.R

**Purpose**

Evaluate external imaging consistency of retained KOA-WM locus signals using ENIGMA DTI summary statistics.

**Main analyses**

* Compare retained KOA-WM tract-specific loci with ENIGMA DTI signals.
* Apply a +/-500 kb locus-overlap window.
* Estimate overlap rates for retained tract-specific signals.
* Classify external imaging consistency as substantial, intermediate, moderate, or limited.
* Export detailed ENIGMA overlap summaries for signals retained in downstream triangulation.

**Main outputs**

* `Supplementary_Table_17_ENIGMA_overlap_summary.tsv`
* `Supplementary_Table_18_ENIGMA_substantial_overlap.tsv`
* `Supplementary_Table_19_ENIGMA_CCG_R_FA.tsv`
* `Supplementary_Table_20_ENIGMA_SCP_L_FA.tsv`
* `Supplementary_Table_21_ENIGMA_SLF_L_MD.tsv`
* `Supplementary_Table_22_ENIGMA_SLF_R_MD.tsv`

**Manuscript link**

External imaging consistency; Supplementary Tables 17-22.

### 07. KOA_WM_FinnGen_Pain_Triangulation_and_GWAS_Catalog_Context.R

**Purpose**

Evaluate broad pain-liability triangulation and GWAS Catalog pleiotropic context.

**Main analyses**

* Compare KOA-WM, KOA-pain, and pain-WM conjFDR locus sets.
* Quantify pairwise and three-way locus overlap using a +/-500 kb window.
* Summarize GWAS Catalog context around loci.
* Classify GWAS Catalog overlap density.
* Exclude OA, pain, and WM/brain traits to summarize broader pleiotropic traits.
* Summarize g:Profiler pathway context for pleiotropic loci.

**Main outputs**

* `Supplementary_Table_23_FinnGen_pain_triangulation.tsv`
* `Supplementary_Table_24_GWAS_Catalog_locus_context.tsv`
* `Supplementary_Table_25_GWAS_Catalog_density_classification.tsv`
* `Supplementary_Table_26_pleiotropic_traits.tsv`
* `Supplementary_Table_27_gProfiler_pleiotropy_pathways.tsv`
* `Supplementary_Table_28_locus_pleiotropy_distribution.tsv`
* `Supplementary_Figure_26_GWAS_Catalog_context_source.png`

**Manuscript link**

Pain triangulation and GWAS Catalog context; Supplementary Figure 26; Supplementary Tables 23-28.

### 08. KOA_WM_Clinical_Multiomics_miRNA_Context.R

**Purpose**

Reanalyze a published KOA clinical multi-omics cohort as exploratory disease-state miRNA context.

**Main analyses**

* Define clinical endotypes using K-means clustering of standardized baseline clinical variables.
* Summarize endotype-specific differentially expressed miRNAs across plasma, synovial fluid, and urine.
* Predict high-confidence miRNA target genes using miRDB score filtering.
* Test overlap or depletion of GWAS-mapped genes in miRNA target networks using Fisher exact tests.
* Fit linear mixed-effects models linking miRNA expression to baseline pain intensity.
* Prioritize miRNAs by joint differential expression and pain-association evidence.

**Main outputs**

* `Supplementary_Table_29_differentially_expressed_miRNAs.tsv`
* `Supplementary_Table_30_miRDB_high_confidence_targets.tsv`
* `Supplementary_Table_31_miRNA_GWAS_gene_overlap_tests.tsv`
* `Supplementary_Table_32_prioritized_miRNAs.tsv`
* `Supplementary_Table_33_miRNA_pain_mixed_models.tsv`
* `Supplementary_Figure_27_miRNA_context_source.png`

**Manuscript link**

Exploratory clinical molecular context; Supplementary Figures 27-29; Supplementary Tables 29-33.

### 09. KOA_WM_MHC_nonMHC_Colocalization_CAUSE_Sensitivity.R

**Purpose**

Evaluate MHC dependence, local colocalization, CAUSE sharing-versus-causal models, and post hoc evidence tiers.

**Main analyses**

* Annotate KOA-WM index loci by MHC/HLA location.
* Parse no-MHC LDSC sensitivity summaries.
* Parse MTAG-derived sensitivity summaries used only for colocalization and CAUSE sensitivity layers.
* Run or summarize coloc.abf within +/-1 Mb regions around KOA-WM index loci.
* Classify PP.H4 support as strong, moderate, same-locus distinct-signal, or no colocalization support.
* Summarize external-trait locus sensitivity for BMI, FinnGen pain, and depression.
* Summarize CAUSE model comparisons between sharing and causal models.
* Create post hoc evidence tiers.

**Main outputs**

* `Supplementary_Table_34_MHC_region_stratification_summary.tsv`
* `Supplementary_Table_35_no_MHC_LDSC_rg_summary.tsv`
* `Supplementary_Table_36_MTAG_MR_sensitivity.tsv`
* `Supplementary_Table_37_local_coloc_summary.tsv`
* `Supplementary_Table_38_external_trait_locus_sensitivity.tsv`
* `Supplementary_Table_39_CAUSE_model_comparison.tsv`
* `Supplementary_Table_40_final_posthoc_evidence_tiers.tsv`

**Manuscript link**

Post hoc causal-sensitivity analyses; Supplementary Figures 30-31; Supplementary Tables 34-40.

### 10. KOA_WM_MVMR_GSMR2_HEIDI_Network.R

**Purpose**

Evaluate conditional sensitivity using weighted multivariable MR and GSMR2/HEIDI network edges.

**Main analyses**

* Build union no-MHC instrument sets for KOA, BMI, and FinnGen pain.
* Fit weighted multivariable MR models for WM outcomes conditioned on KOA, BMI, and FinnGen pain effects.
* Compute instrument-strength summaries and compare univariable versus conditional KOA-WM estimates.
* Summarize GSMR2 network quality control and HEIDI-filtered effect estimates.
* Classify network edges as null, nominal, or FDR-significant.

**Main outputs**

* `Supplementary_Table_41_MHC_annotation_index_loci.tsv`
* `Supplementary_Table_42_MHC_nonMHC_counts_by_trait.tsv`
* `Supplementary_Table_43_nonMHC_coloc_summary.tsv`
* `Supplementary_Table_44_nonMHC_external_trait_sensitivity.tsv`
* `Supplementary_Table_45_nonMHC_final_evidence_tiers.tsv`
* `Supplementary_Table_46_noMHC_MVMR_instruments.tsv`
* `Supplementary_Table_47_weighted_MVMR_estimates.tsv`
* `Supplementary_Table_48_univariable_vs_MVMR_comparison.tsv`
* `Supplementary_Table_49_GSMR2_HEIDI_QC_summary.tsv`
* `Supplementary_Table_50_GSMR2_network_effect_estimates.tsv`
* `Supplementary_Table_51_integrated_network_interpretation.tsv`
* `Supplementary_Figure_36_GSMR2_network_mean_F_statistics.png`

**Manuscript link**

Weighted MVMR, GSMR2/HEIDI, and network sensitivity; Supplementary Figures 32-36; Supplementary Tables 41-51.

### 11. KOA_WM_OA_Subtype_BMI_Pain_MHC_Evidence_Integration.R

**Purpose**

Integrate OA subtype specificity with BMI, pain, MHC, colocalization, and GSMR2 evidence.

**Main analyses**

* Summarize OA subtype source-file manifest and LDSC heritability QC for hip OA, hand OA, and all-OA files.
* Perform +/-500 kb lookup of KOA-WM index loci in hip OA, hand OA, and all-OA summary statistics.
* Classify loci into KOA-enriched, lower-limb OA, hand/general OA, general OA only, pan-OA, or weak/unclassified groups.
* Integrate BMI, FinnGen pain, depression, MHC/HLA, colocalization, and GSMR2 evidence.
* Generate final locus-level and trait-level evidence strata.
* Summarize evidence-stratum counts used in the Stage5 manuscript interpretation.

**Main outputs**

* `Supplementary_Table_52_OA_subtype_source_manifest.tsv`
* `Supplementary_Table_53_OA_subtype_LDSC_h2_QC.tsv`
* `Supplementary_Table_54_OA_subtype_locus_lookup.tsv`
* `Supplementary_Table_55_OA_subtype_specificity_classification.tsv`
* `Supplementary_Table_56_integrated_locus_evidence_strata.tsv`
* `Supplementary_Table_57_trait_level_evidence_summary.tsv`
* `Supplementary_Table_58_locus_stratum_counts.tsv`
* `Supplementary_Table_59_trait_stratum_counts.tsv`
* `Supplementary_Table_60_key_GSMR2_edges_for_Stage5.tsv`
* `Supplementary_Figure_37_Stage5_integrated_evidence_heatmap.png`
* `Supplementary_Figure_38_Stage5_locus_evidence_stratum_counts.png`

**Manuscript link**

OA subtype-aware evidence integration; Supplementary Figures 37-38; Supplementary Tables 52-60.

### 12. KOA_WM_WM_Genetic_Factor_Domain_Convergence.R

**Purpose**

Summarize WM genetic factor structure, external-trait alignment, FA/MD directionality, domain heterogeneity, and immune-metabolic-pain mechanism-axis convergence.

**Main analyses**

* Build the retained WM-WM LDSC genetic-correlation matrix.
* Regularize the matrix to the nearest positive definite correlation matrix.
* Construct descriptive general WM, FA-enriched, and MD-enriched genetic factor summaries.
* Project KOA, BMI, FinnGen pain, and depression genetic correlations onto WM factor loadings.
* Summarize tract-domain evidence scores.
* Test FA/MD directionality of KOA-WM rg estimates.
* Classify locus-trait rows into BMI/metabolic-sensitive, pain-sensitive, MHC/immune-sensitive, BMI-and-pain-sensitive, or weak/unclassified mechanism-axis groups.

**Main outputs**

* `Supplementary_Table_61_WM_factor_trait_manifest.tsv`
* `Supplementary_Table_62_pairwise_WM_LDSC_rg.tsv`
* `Supplementary_Table_63_external_trait_to_WM_rg.tsv`
* `Supplementary_Table_64_WM_genetic_factor_loadings.tsv`
* `Supplementary_Table_65_external_trait_factor_alignment.tsv`
* `Supplementary_Table_66_factor_interpretation_summary.tsv`
* `Supplementary_Table_67_domain_level_Fisher_tests.tsv`
* `Supplementary_Table_68_metric_level_Fisher_tests.tsv`
* `Supplementary_Table_69_permutation_and_FA_MD_directionality.tsv`
* `Supplementary_Table_70_mechanism_axis_classification.tsv`
* `Supplementary_Figure_39_Stage6_factor_external_rg.png`
* `Supplementary_Figure_40_Stage6_domain_evidence_score.png`
* `Supplementary_Figure_41_Stage6_mechanism_axis_by_domain.png`
* `Supplementary_Figure_42_WM_rg_heatmap_matrix.tsv`
* `Supplementary_Figure_43_WM_factor_loadings.png`

**Manuscript link**

WM genetic factor/domain convergence; Supplementary Figures 39-43; Supplementary Tables 61-70.

### 13. KOA_WM_Publication_Tables_and_Figures.R

**Purpose**

Assemble publication-ready summary tables and source figures from cleaned outputs.

**Main analyses**

* Check the presence of Supplementary Tables 1-70.
* Generate the evidence-tiered Table 1 source file.
* Generate compact Stage5 and Stage6 summary plots from the final evidence tables.
* Export publication-level source files without rerunning statistical analyses.

**Main outputs**

* `publication_table_inventory.tsv`
* `Table_1_evidence_tiered_interpretation.tsv`
* `Publication_summary_Stage5_counts.png`
* `Publication_summary_Stage6_KOA_factor_alignment.png`

**Manuscript link**

Table 1 and final publication source-file assembly.

## Run order

Run the scripts in numerical order after setting paths and preparing input files:

```bash
Rscript "code/01. KOA_WM_Global_Genetic_Correlation.R"
Rscript "code/02. KOA_WM_Local_Genetic_Correlation_LAVA.R"
python  "code/03. KOA_WM_condFDR_conjFDR_Enrichment.py"
Rscript "code/04. KOA_WM_Locus_Pruning_Merging_and_Directionality.R"
Rscript "code/05. KOA_WM_Functional_Annotation_and_Enrichment.R"
Rscript "code/06. KOA_WM_ENIGMA_External_Consistency.R"
Rscript "code/07. KOA_WM_FinnGen_Pain_Triangulation_and_GWAS_Catalog_Context.R"
Rscript "code/08. KOA_WM_Clinical_Multiomics_miRNA_Context.R"
Rscript "code/09. KOA_WM_MHC_nonMHC_Colocalization_CAUSE_Sensitivity.R"
Rscript "code/10. KOA_WM_MVMR_GSMR2_HEIDI_Network.R"
Rscript "code/11. KOA_WM_OA_Subtype_BMI_Pain_MHC_Evidence_Integration.R"
Rscript "code/12. KOA_WM_WM_Genetic_Factor_Domain_Convergence.R"
Rscript "code/13. KOA_WM_Publication_Tables_and_Figures.R"
```

## Interpretation boundary

The scripts follow the interpretation used in the final manuscript. condFDR and conjFDR identify SNP-level P-value enrichment records and do not establish shared causal variants. MHC/HLA signals are excluded from primary discovery and used as immune-genetic context. ENIGMA provides external imaging consistency, not full tract-level replication. FinnGen PAIN is used as a broad pain-liability comparator, not a KOA-specific pain phenotype. The clinical miRNA layer is interpreted as disease-state molecular context and does not provide genotype-linked molecular validation. Stage5 and Stage6 organize evidence across OA subtype, BMI, pain, MHC, WM metric, and tract-domain layers.

## Citation

If using this code, please cite the associated manuscript and the original data/software resources described in the manuscript.
