#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(coloc)
  library(cause)
})

project_dir <- Sys.getenv("KOA_WM_PROJECT_DIR", ".")
data_dir <- Sys.getenv("KOA_WM_DATA_DIR", file.path(project_dir, "data"))
results_dir <- Sys.getenv("KOA_WM_RESULTS_DIR", file.path(project_dir, "results"))
tables_dir <- file.path(results_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) data.table::fread(path, data.table = FALSE)
write_tsv <- function(x, path) data.table::fwrite(as.data.frame(x), path, sep = "\t", na = "NA")

index_loci <- read_tsv(file.path(tables_dir, "Supplementary_Table_13_distinct_cross_trait_loci.tsv")) |>
  mutate(CHR = as.character(CHR), BP = as.numeric(BP),
         is_mhc = CHR == "6" & BP >= 28477797 & BP <= 33448354)

coloc_manifest <- read_tsv(Sys.getenv("KOA_WM_COLOC_INPUT_MANIFEST", file.path(data_dir, "coloc_input_manifest.tsv")))
coloc_summary <- lapply(seq_len(nrow(coloc_manifest)), function(i) {
  row <- coloc_manifest[i, ]
  d1 <- read_tsv(row$exposure_file)
  d2 <- read_tsv(row$outcome_file)
  dat <- inner_join(d1, d2, by = "SNP", suffix = c("_KOA", "_WM"))
  fit <- coloc.abf(dataset1 = list(beta = dat$BETA_KOA, varbeta = dat$SE_KOA^2, snp = dat$SNP, type = "cc"),
                   dataset2 = list(beta = dat$BETA_WM, varbeta = dat$SE_WM^2, snp = dat$SNP, type = "quant"))
  cbind(locus_id = row$locus_id, as.data.frame(t(fit$summary)))
}) |> bind_rows()

coloc_summary <- coloc_summary |>
  mutate(PP.H4.abf = as.numeric(PP.H4.abf), PP.H3.abf = as.numeric(PP.H3.abf),
         coloc_class = case_when(PP.H4.abf >= 0.80 ~ "strong_colocalization",
                                 PP.H4.abf >= 0.50 ~ "moderate_colocalization",
                                 PP.H3.abf > PP.H4.abf ~ "same_locus_distinct_signal",
                                 TRUE ~ "no_colocalization_support"))

external <- read_tsv(Sys.getenv("KOA_WM_EXTERNAL_TRAIT_SENSITIVITY", file.path(data_dir, "external_trait_locus_sensitivity.tsv")))
cause_summary <- read_tsv(Sys.getenv("KOA_WM_CAUSE_SUMMARY", file.path(data_dir, "cause_model_comparison_summary.tsv"))) |>
  mutate(cause_class = case_when(sharing_vs_causal_p < 0.05 & delta_elpd_sharing_minus_causal > 0 ~ "sharing_favored",
                                 causal_vs_sharing_p < 0.05 & delta_elpd_causal_minus_sharing > 0 ~ "causal_favored",
                                 TRUE ~ "no_clear_preference"))

tiers <- index_loci |>
  left_join(coloc_summary |> select(locus_id, PP.H4.abf, coloc_class), by = "locus_id") |>
  left_join(external, by = "locus_id") |>
  mutate(evidence_tier = case_when(coloc_class %in% c("strong_colocalization", "moderate_colocalization") & !is_mhc ~ "shared_locus_support",
                                   bmi_signal | pain_signal | is_mhc ~ "external_trait_or_MHC_sensitive",
                                   TRUE ~ "not_robust"))

write_tsv(index_loci, file.path(tables_dir, "Supplementary_Table_34_MHC_region_stratification_summary.tsv"))
write_tsv(read_tsv(file.path(tables_dir, "Supplementary_Table_2_Global_LDSC_HDL_rg.tsv")) |> filter(method == "LDSC"), file.path(tables_dir, "Supplementary_Table_35_no_MHC_LDSC_rg_summary.tsv"))
write_tsv(read_tsv(Sys.getenv("KOA_WM_MTAG_MR_SUMMARY", file.path(data_dir, "mtag_mr_summary.tsv"))), file.path(tables_dir, "Supplementary_Table_36_MTAG_MR_sensitivity.tsv"))
write_tsv(coloc_summary, file.path(tables_dir, "Supplementary_Table_37_local_coloc_summary.tsv"))
write_tsv(external, file.path(tables_dir, "Supplementary_Table_38_external_trait_locus_sensitivity.tsv"))
write_tsv(cause_summary, file.path(tables_dir, "Supplementary_Table_39_CAUSE_model_comparison.tsv"))
write_tsv(tiers, file.path(tables_dir, "Supplementary_Table_40_final_posthoc_evidence_tiers.tsv"))

message("Completed MHC/non-MHC, colocalization, CAUSE, and post hoc evidence-tier analysis.")
