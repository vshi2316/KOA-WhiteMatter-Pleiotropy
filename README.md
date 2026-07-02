# KOA White Matter Pleiotropy

Repository for the analysis code accompanying the manuscript:

**Shared genetic architecture between knee osteoarthritis and white matter microstructure implicates immune metabolic and pain sensitive multimorbidity in older adults**

The workflow evaluates genetic overlap between knee osteoarthritis (KOA) liability and diffusion magnetic resonance imaging derived white matter (WM) microstructure phenotypes. The analysis combines genome wide genetic correlation, local genetic correlation, conditional and conjunctional false discovery rate enrichment, locus consolidation, functional annotation, external imaging consistency, FinnGen pain triangulation, GWAS Catalog context, older adult phenotype context, MHC sensitivity, local colocalization, CAUSE, weighted multivariable Mendelian randomization, GSMR2 with HEIDI filtering, osteoarthritis subtype lookup, and WM genetic factor analysis.

The code release contains the manuscript facing analysis scripts. Input data from third party cohorts, genome wide association studies, LD reference panels, FUMA, g:Profiler, FinnGen, GWAS Catalog, ENIGMA, and aging cohort resources must be obtained from the original sources under their access conditions.

## Repository Structure

```text
README.md
01. KOA_WM_Global_Genetic_Correlation.R
02. KOA_WM_Local_Genetic_Correlation_LAVA.R
03. KOA_WM_condFDR_conjFDR_Enrichment.py
04. KOA_WM_Locus_Pruning_Merging_and_Directionality.R
05. KOA_WM_Functional_Annotation_and_Enrichment.R
06. KOA_WM_ENIGMA_External_Consistency.R
07. KOA_WM_FinnGen_Pain_Triangulation_and_GWAS_Catalog_Context.R
08. KOA_WM_OlderAdult_Phenotype_Context.R
09. KOA_WM_MHC_nonMHC_Colocalization_CAUSE_Sensitivity.R
10. KOA_WM_MVMR_GSMR2_HEIDI_Network.R
11. KOA_WM_OA_Subtype_BMI_Pain_MHC_Evidence_Integration.R
12. KOA_WM_WM_Genetic_Factor_Domain_Convergence.R
```

## Required Input Data

The scripts expect local paths to prepared summary statistics, intermediate analysis outputs, annotation files, and harmonized cohort files. The repository provides code only.

1. KOA genome wide association summary statistics from the Musculoskeletal Knowledge Portal osteoarthritis genetics resource.

2. UK Biobank Brain Imaging Genetics diffusion magnetic resonance imaging summary statistics for fractional anisotropy and mean diffusivity WM phenotypes.

3. LDSC compatible summary statistics for KOA, WM phenotypes, body mass index, depression, FinnGen pain, hip osteoarthritis, hand osteoarthritis, and all osteoarthritis.

4. HDL output files for KOA and WM global genetic correlation analyses.

5. LAVA reference files and LAVA local genetic correlation output.

6. 1000 Genomes Phase 3 European reference files, HapMap3 SNP list, European LD score reference files, and PLINK compatible LD reference files.

7. Manifest of KOA and WM summary statistic pairs for conditional and conjunctional false discovery rate analyses.

8. FUMA SNP2GENE output for candidate SNP annotation.

9. g:Profiler enrichment output for candidate gene sets and GWAS Catalog pleiotropy context.

10. ENIGMA DTI summary statistics and significant loci used for external imaging consistency checks.

11. FinnGen R11 PAIN endpoint summary statistics.

12. GWAS Catalog association files.

13. MTAG derived summary statistics used for local colocalization and CAUSE sensitivity analyses.

14. Local colocalization input manifest and harmonized exposure and outcome files for KOA and WM locus windows.

15. CAUSE model comparison summaries for KOA, WM, and external trait pairs.

16. MVMR union instrument file, harmonized MVMR effect file, univariable KOA to WM Mendelian randomization estimates, GSMR2 quality control summary, and GSMR2 network effect estimates.

17. Hip osteoarthritis, hand osteoarthritis, and all osteoarthritis source manifest, LDSC heritability quality control table, and locus lookup table.

18. Pairwise WM to WM LDSC genetic correlation matrix and external trait to WM LDSC genetic correlation table.

19. Harmonized older adult cohort CSV files for the phenotype context module:

```text
HRS_2018_stage7.csv
ELSA_stage7_long.csv
MHAS_stage7_long.csv
CHARLS_stage7_long.csv
SHARE_stage7_long.csv
```

## Local Configuration

Environment variables keep local paths outside the code.

```r
Sys.setenv(
  KOA_WM_PROJECT_DIR = "path/to/KOA-WhiteMatter-Pleiotropy",
  KOA_WM_DATA_DIR = "path/to/local/input_data",
  KOA_WM_RESULTS_DIR = "path/to/local/results",
  KOA_WM_STAGE7_INPUT_DIR = "path/to/stage7_aging_cohorts",
  KOA_WM_STAGE7_RESULTS_DIR = "path/to/stage7_results",
  LDSC_REF_DIR = "path/to/eur_w_ld_chr",
  LAVA_REF_DIR = "path/to/LAVA/reference",
  PLINK_EXE = "path/to/plink",
  GCTA_EXE = "path/to/gcta64"
)
```

Default output folders are created under `results/tables`, `results/figures`, `results/logs`, and `results/stage7_older_adult_phenotype_context` unless environment variables specify other locations.

## Software

R packages used across the workflow:

```text
data.table, dplyr, tidyr, tibble, readr, stringr, purrr, magrittr,
ggplot2, ggrepel, patchwork, cowplot, scales, RColorBrewer, viridis,
GenomicRanges, IRanges, S4Vectors, GenomeInfoDb, rtracklayer, biomaRt,
org.Hs.eg.db, AnnotationDbi, Matrix, coloc, cause, TwoSampleMR,
MendelianRandomization, ieugwasr, MRPRESSO, gwasglue, gprofiler2,
clusterProfiler, enrichplot, DOSE, broom, broom.mixed, metafor,
openxlsx, writexl
```

Python packages used by the conditional and conjunctional false discovery rate script:

```text
python >= 3.8, numpy, scipy, pandas, matplotlib, seaborn
```

External software and web resources used by the workflow:

```text
LDSC, HDL, LAVA, pleioFDR, PLINK 1.9, MTAG, GCTA GSMR2,
FUMA SNP2GENE, g:Profiler, GWAS Catalog
```

## Analysis Modules

### 01. KOA_WM_Global_Genetic_Correlation.R

Purpose:

Summarize genome wide genetic correlation between KOA liability and diffusion magnetic resonance imaging derived WM phenotypes using LDSC and HDL.

Main outputs:

```text
Supplementary_Table_1_WM_trait_metadata.tsv
Supplementary_Table_2_Global_LDSC_HDL_rg.tsv
Figure_1A_global_rg_source.png
Figure_1A_direction_balance_source.tsv
```

Manuscript link:

Genome wide genetic correlation, Figure 1A, Supplementary Tables 1 and 2.

### 02. KOA_WM_Local_Genetic_Correlation_LAVA.R

Purpose:

Summarize local genetic correlation between KOA and WM phenotypes using LAVA.

Main outputs:

```text
Supplementary_Table_3_LAVA_nominal_local_rg.tsv
Supplementary_Table_4_LAVA_FA_FDR_significant.tsv
Supplementary_Table_5_LAVA_MD_FDR_significant.tsv
Figure_1B_LAVA_local_rg_scatter_source.png
Figure_1C_LAVA_volcano_source.png
```

Manuscript link:

Local genetic correlation, Figure 1B and Figure 1C, Supplementary Tables 3 to 5.

### 03. KOA_WM_condFDR_conjFDR_Enrichment.py

Purpose:

Run conditional false discovery rate, conjunctional false discovery rate, and conditional quantile quantile analyses for KOA and WM phenotype pairs.

Main outputs:

```text
Supplementary_Table_6_7_condFDR_outputs.tsv
conjFDR_SNP_level_enrichment_records.tsv
Figure_2A_F_condFDR_conjFDR_count_source.tsv
Conditional quantile quantile plots for Supplementary Figures 1 to 13
```

Manuscript link:

Variant level cross trait enrichment, Figure 2A to Figure 2F, Supplementary Figures 1 to 13, Supplementary Tables 6 and 7.

### 04. KOA_WM_Locus_Pruning_Merging_and_Directionality.R

Purpose:

Convert SNP level enrichment records into locus level and tract level summaries.

Main outputs:

```text
Supplementary_Table_8_tract_category_counts.tsv
Supplementary_Table_9_hemispheric_specificity.tsv
Supplementary_Table_10_top_lead_SNPs.tsv
Supplementary_Table_11_KOA_FA_effect_direction.tsv
Supplementary_Table_12_KOA_MD_effect_direction.tsv
Supplementary_Table_13_distinct_cross_trait_loci.tsv
Figure_2G_I_tract_hemisphere_source.png
Figure_3A_distinct_loci_source.png
```

Manuscript link:

Locus consolidation and directionality, Figure 2G to Figure 2I, Figure 3, Supplementary Figures 14 to 25, Supplementary Tables 8 to 13.

### 05. KOA_WM_Functional_Annotation_and_Enrichment.R

Purpose:

Summarize functional annotation and enrichment results for candidate KOA and WM enrichment regions.

Main outputs:

```text
Supplementary_Table_14_candidate_shared_regions.tsv
Supplementary_Table_15_candidate_and_LD_linked_SNPs.tsv
Supplementary_Table_16_functional_annotation_summary.tsv
Figure_4A_B_functional_annotation_source.png
Figure_4C_gene_mapping_source.tsv
Figure_4D_gprofiler_enrichment_source.tsv
```

Manuscript link:

Functional annotation and enrichment analysis, Figure 4, Supplementary Tables 14 to 16.

### 06. KOA_WM_ENIGMA_External_Consistency.R

Purpose:

Evaluate external imaging consistency of retained KOA and WM locus signals using ENIGMA DTI summary statistics.

Main outputs:

```text
Supplementary_Table_17_ENIGMA_overlap_summary.tsv
Supplementary_Table_18_ENIGMA_substantial_overlap.tsv
Supplementary_Table_19_ENIGMA_CCG_R_FA.tsv
Supplementary_Table_20_ENIGMA_SCP_L_FA.tsv
Supplementary_Table_21_ENIGMA_SLF_L_MD.tsv
Supplementary_Table_22_ENIGMA_SLF_R_MD.tsv
```

Manuscript link:

External imaging consistency, Supplementary Tables 17 to 22.

### 07. KOA_WM_FinnGen_Pain_Triangulation_and_GWAS_Catalog_Context.R

Purpose:

Evaluate FinnGen pain triangulation and GWAS Catalog pleiotropy context.

Main outputs:

```text
Supplementary_Table_23_FinnGen_pain_triangulation.tsv
Supplementary_Table_24_GWAS_Catalog_locus_context.tsv
Supplementary_Table_25_GWAS_Catalog_density_classification.tsv
Supplementary_Table_26_pleiotropic_traits.tsv
Supplementary_Table_27_gProfiler_pleiotropy_pathways.tsv
Supplementary_Table_28_locus_pleiotropy_distribution.tsv
Supplementary_Figure_26_GWAS_Catalog_context_source.png
```

Manuscript link:

FinnGen pain triangulation and GWAS Catalog context, Supplementary Figure 26, Supplementary Tables 23 to 28.

### 08. KOA_WM_OlderAdult_Phenotype_Context.R

Purpose:

Generate the older adult phenotype context analyses for arthritis, pain, walking difficulty, falls, activities of daily living, instrumental activities of daily living, functional burden, and cognition sensitivity outcomes.

Main outputs:

```text
Stage7_HRS_variable_coverage.csv
Stage7_auxiliary_variable_coverage.csv
Stage7_HRS_primary_models.csv
Stage7_HRS_pain_attenuation_gradient.csv
Stage7_auxiliary_cross_sectional_models.csv
Stage7_auxiliary_lagged_models.csv
Stage7_cross_cohort_random_effects_meta.csv
Stage7_leave_one_cohort_out_meta.csv
Stage7_outcome_specificity_summary.csv
Figure_5_older_adult_phenotype_context.pdf
Figure_5_older_adult_phenotype_context.png
Stage7_sessionInfo.txt
```

Manuscript link:

Older adult phenotype context, Figure 5, Supplementary Figures 27 to 29, Supplementary Tables 29 to 33.

### 09. KOA_WM_MHC_nonMHC_Colocalization_CAUSE_Sensitivity.R

Purpose:

Evaluate MHC dependence, local colocalization, CAUSE sharing and causal model comparisons, and post hoc evidence tiers.

Main outputs:

```text
Supplementary_Table_34_MHC_region_stratification_summary.tsv
Supplementary_Table_35_no_MHC_LDSC_rg_summary.tsv
Supplementary_Table_36_MTAG_MR_sensitivity.tsv
Supplementary_Table_37_local_coloc_summary.tsv
Supplementary_Table_38_external_trait_locus_sensitivity.tsv
Supplementary_Table_39_CAUSE_model_comparison.tsv
Supplementary_Table_40_final_posthoc_evidence_tiers.tsv
```

Manuscript link:

Post hoc genetic sensitivity, Supplementary Figures 30 and 31, Supplementary Tables 34 to 40.

### 10. KOA_WM_MVMR_GSMR2_HEIDI_Network.R

Purpose:

Evaluate conditional sensitivity using weighted multivariable Mendelian randomization and GSMR2 with HEIDI filtering.

Main outputs:

```text
Supplementary_Table_41_MHC_annotation_index_loci.tsv
Supplementary_Table_42_MHC_nonMHC_counts_by_trait.tsv
Supplementary_Table_43_nonMHC_coloc_summary.tsv
Supplementary_Table_44_nonMHC_external_trait_sensitivity.tsv
Supplementary_Table_45_nonMHC_final_evidence_tiers.tsv
Supplementary_Table_46_noMHC_MVMR_instruments.tsv
Supplementary_Table_47_weighted_MVMR_estimates.tsv
Supplementary_Table_48_univariable_vs_MVMR_comparison.tsv
Supplementary_Table_49_GSMR2_HEIDI_QC_summary.tsv
Supplementary_Table_50_GSMR2_network_effect_estimates.tsv
Supplementary_Table_51_integrated_network_interpretation.tsv
Supplementary_Figure_36_GSMR2_network_mean_F_statistics.png
```

Manuscript link:

Weighted multivariable Mendelian randomization and GSMR2 with HEIDI filtering, Supplementary Figures 32 to 36, Supplementary Tables 41 to 51.

### 11. KOA_WM_OA_Subtype_BMI_Pain_MHC_Evidence_Integration.R

Purpose:

Integrate osteoarthritis subtype specificity with body mass index, pain, MHC, colocalization, and GSMR2 evidence.

Main outputs:

```text
Supplementary_Table_52_OA_subtype_source_manifest.tsv
Supplementary_Table_53_OA_subtype_LDSC_h2_QC.tsv
Supplementary_Table_54_OA_subtype_locus_lookup.tsv
Supplementary_Table_55_OA_subtype_specificity_classification.tsv
Supplementary_Table_56_integrated_locus_evidence_strata.tsv
Supplementary_Table_57_trait_level_evidence_summary.tsv
Supplementary_Table_58_locus_stratum_counts.tsv
Supplementary_Table_59_trait_stratum_counts.tsv
Supplementary_Table_60_key_GSMR2_edges_for_Stage5.tsv
Supplementary_Figure_37_Stage5_integrated_evidence_heatmap.png
Supplementary_Figure_38_Stage5_locus_evidence_stratum_counts.png
```

Manuscript link:

Osteoarthritis subtype aware evidence integration, Supplementary Figures 37 and 38, Supplementary Tables 52 to 60.

### 12. KOA_WM_WM_Genetic_Factor_Domain_Convergence.R

Purpose:

Summarize WM genetic factor structure, external trait alignment, FA and MD directionality, tract domain heterogeneity, and immune metabolic pain mechanism axis convergence.

Main outputs:

```text
Supplementary_Table_61_WM_factor_trait_manifest.tsv
Supplementary_Table_62_pairwise_WM_LDSC_rg.tsv
Supplementary_Table_63_external_trait_to_WM_rg.tsv
Supplementary_Table_64_WM_genetic_factor_loadings.tsv
Supplementary_Table_65_external_trait_factor_alignment.tsv
Supplementary_Table_66_factor_interpretation_summary.tsv
Supplementary_Table_67_domain_level_Fisher_tests.tsv
Supplementary_Table_68_metric_level_Fisher_tests.tsv
Supplementary_Table_69_permutation_and_FA_MD_directionality.tsv
Supplementary_Table_70_mechanism_axis_classification.tsv
Supplementary_Figure_39_Stage6_factor_external_rg.png
Supplementary_Figure_40_Stage6_domain_evidence_score.png
Supplementary_Figure_41_Stage6_mechanism_axis_by_domain.png
Supplementary_Figure_42_WM_rg_heatmap_matrix.tsv
Supplementary_Figure_43_WM_factor_loadings.png
```

Manuscript link:

WM genetic factor and domain convergence, Supplementary Figures 39 to 43, Supplementary Tables 61 to 70.

## Run Order

Run scripts from the repository root after configuring paths and preparing input files.

```bash
Rscript "01. KOA_WM_Global_Genetic_Correlation.R"
Rscript "02. KOA_WM_Local_Genetic_Correlation_LAVA.R"
python  "03. KOA_WM_condFDR_conjFDR_Enrichment.py"
Rscript "04. KOA_WM_Locus_Pruning_Merging_and_Directionality.R"
Rscript "05. KOA_WM_Functional_Annotation_and_Enrichment.R"
Rscript "06. KOA_WM_ENIGMA_External_Consistency.R"
Rscript "07. KOA_WM_FinnGen_Pain_Triangulation_and_GWAS_Catalog_Context.R"
Rscript "08. KOA_WM_OlderAdult_Phenotype_Context.R"
Rscript "09. KOA_WM_MHC_nonMHC_Colocalization_CAUSE_Sensitivity.R"
Rscript "10. KOA_WM_MVMR_GSMR2_HEIDI_Network.R"
Rscript "11. KOA_WM_OA_Subtype_BMI_Pain_MHC_Evidence_Integration.R"
Rscript "12. KOA_WM_WM_Genetic_Factor_Domain_Convergence.R"
```

## Interpretation Boundary

The conditional and conjunctional false discovery rate modules identify cross trait P value enrichment and require locus level sensitivity checks before biological interpretation. MHC signals are handled as immune genetic context because complex LD can affect enrichment based discovery. ENIGMA provides external imaging consistency across coordinated DTI resources. FinnGen PAIN is used as a broad pain liability comparator. The older adult cohort module quantifies functional burden associated with arthritis and pain in HRS, ELSA, MHAS, CHARLS, and SHARE. Direct translation of KOA and WM genetic signals requires genotype matched brain imaging cohorts.

## Citation

Please cite the associated manuscript and the original data and software resources listed in the manuscript when using this code.
