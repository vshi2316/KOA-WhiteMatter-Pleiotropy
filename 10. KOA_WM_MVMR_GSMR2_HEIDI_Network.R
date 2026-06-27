#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

project_dir <- Sys.getenv("KOA_WM_PROJECT_DIR", ".")
data_dir <- Sys.getenv("KOA_WM_DATA_DIR", file.path(project_dir, "data"))
results_dir <- Sys.getenv("KOA_WM_RESULTS_DIR", file.path(project_dir, "results"))
tables_dir <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) data.table::fread(path, data.table = FALSE)
write_tsv <- function(x, path) data.table::fwrite(as.data.frame(x), path, sep = "\t", na = "NA")

instruments <- read_tsv(Sys.getenv("KOA_WM_MVMR_INSTRUMENTS", file.path(data_dir, "mvmr_union_instruments.tsv")))
mvmr_dat <- read_tsv(Sys.getenv("KOA_WM_MVMR_HARMONIZED", file.path(data_dir, "mvmr_harmonized_effects.tsv")))

fit_mvmr <- function(dat) {
  fit <- lm(beta_wm ~ beta_koa + beta_bmi + beta_pain - 1, weights = 1 / se_wm^2, data = dat)
  co <- summary(fit)$coefficients
  data.frame(exposure = rownames(co), beta = co[, "Estimate"], se = co[, "Std. Error"], p = co[, "Pr(>|t|)"], row.names = NULL)
}

mvmr <- mvmr_dat |> group_by(wm_label) |> group_modify(~fit_mvmr(.x)) |> ungroup() |> group_by(exposure) |> mutate(q = p.adjust(p, "BH")) |> ungroup()
univ <- read_tsv(Sys.getenv("KOA_WM_UNIVARIABLE_MR", file.path(data_dir, "univariable_koa_wm_estimates.tsv")))
comparison <- univ |>
  left_join(mvmr |> filter(exposure == "beta_koa") |> select(wm_label, beta_mvmr = beta, se_mvmr = se, p_mvmr = p), by = "wm_label") |>
  mutate(attenuation_percent = (1 - abs(beta_mvmr) / abs(beta_univariable)) * 100,
         direction_change = sign(beta_univariable) != sign(beta_mvmr))

gsmr_qc <- read_tsv(Sys.getenv("KOA_WM_GSMR_QC", file.path(data_dir, "gsmr2_network_qc.tsv")))
gsmr <- read_tsv(Sys.getenv("KOA_WM_GSMR_EFFECTS", file.path(data_dir, "gsmr2_network_effects.tsv"))) |> mutate(q = p.adjust(p, "BH"))
network <- gsmr |> mutate(edge_class = case_when(q < 0.05 ~ "FDR_significant", p < 0.05 ~ "nominal", TRUE ~ "null"))

write_tsv(read_tsv(file.path(tables_dir, "Supplementary_Table_34_MHC_region_stratification_summary.tsv")), file.path(tables_dir, "Supplementary_Table_41_MHC_annotation_index_loci.tsv"))
write_tsv(instruments |> count(wm_label, is_mhc), file.path(tables_dir, "Supplementary_Table_42_MHC_nonMHC_counts_by_trait.tsv"))
write_tsv(read_tsv(file.path(tables_dir, "Supplementary_Table_37_local_coloc_summary.tsv")), file.path(tables_dir, "Supplementary_Table_43_nonMHC_coloc_summary.tsv"))
write_tsv(read_tsv(file.path(tables_dir, "Supplementary_Table_38_external_trait_locus_sensitivity.tsv")), file.path(tables_dir, "Supplementary_Table_44_nonMHC_external_trait_sensitivity.tsv"))
write_tsv(read_tsv(file.path(tables_dir, "Supplementary_Table_40_final_posthoc_evidence_tiers.tsv")), file.path(tables_dir, "Supplementary_Table_45_nonMHC_final_evidence_tiers.tsv"))
write_tsv(instruments, file.path(tables_dir, "Supplementary_Table_46_noMHC_MVMR_instruments.tsv"))
write_tsv(mvmr, file.path(tables_dir, "Supplementary_Table_47_weighted_MVMR_estimates.tsv"))
write_tsv(comparison, file.path(tables_dir, "Supplementary_Table_48_univariable_vs_MVMR_comparison.tsv"))
write_tsv(gsmr_qc, file.path(tables_dir, "Supplementary_Table_49_GSMR2_HEIDI_QC_summary.tsv"))
write_tsv(gsmr, file.path(tables_dir, "Supplementary_Table_50_GSMR2_network_effect_estimates.tsv"))
write_tsv(network, file.path(tables_dir, "Supplementary_Table_51_integrated_network_interpretation.tsv"))

f_stats <- instruments |> group_by(exposure) |> summarise(mean_F = mean(F_stat, na.rm = TRUE), .groups = "drop")
p <- f_stats |> ggplot(aes(x = reorder(exposure, mean_F), y = mean_F)) +
  geom_col(fill = "#2B8CBE") + geom_hline(yintercept = 10, linetype = 2, color = "red") +
  coord_flip() + theme_bw(base_size = 9) + labs(x = NULL, y = "Mean F statistic")
ggsave(file.path(figures_dir, "Supplementary_Figure_36_GSMR2_network_mean_F_statistics.png"), p, width = 6, height = 4, dpi = 300)

message("Completed weighted MVMR and GSMR2/HEIDI network analysis.")
