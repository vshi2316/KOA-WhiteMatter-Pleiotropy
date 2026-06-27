#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(Matrix)
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

wm_meta <- read_tsv(file.path(project_dir, "data_manifest", "wm_trait_metadata.tsv"))
wm_rg <- read_tsv(Sys.getenv("KOA_WM_WM_RG_MATRIX_LONG", file.path(data_dir, "wm_wm_ldsc_rg_long.tsv")))
external_rg <- read_tsv(Sys.getenv("KOA_WM_EXTERNAL_TO_WM_RG", file.path(data_dir, "external_trait_to_wm_ldsc_rg.tsv")))
integrated <- read_tsv(file.path(tables_dir, "Supplementary_Table_56_integrated_locus_evidence_strata.tsv"))

traits <- sort(unique(c(wm_rg$trait1, wm_rg$trait2)))
mat <- matrix(0, length(traits), length(traits), dimnames = list(traits, traits))
diag(mat) <- 1
for (i in seq_len(nrow(wm_rg))) {
  mat[wm_rg$trait1[i], wm_rg$trait2[i]] <- wm_rg$rg[i]
  mat[wm_rg$trait2[i], wm_rg$trait1[i]] <- wm_rg$rg[i]
}
mat <- as.matrix(Matrix::nearPD(mat, corr = TRUE)$mat)
pca <- prcomp(mat, center = TRUE, scale. = FALSE)

loadings <- data.frame(wm_id = traits, general_WM_factor = pca$rotation[, 1], row.names = NULL) |>
  left_join(wm_meta, by = "wm_id") |>
  mutate(FA_specific_factor = ifelse(metric == "FA", abs(general_WM_factor), 0),
         MD_specific_factor = ifelse(metric == "MD", abs(general_WM_factor), 0))

project_factor <- function(load, rg, se) {
  ok <- is.finite(load) & is.finite(rg) & is.finite(se)
  beta <- sum(load[ok] * rg[ok]) / sum(abs(load[ok]))
  se_beta <- sqrt(sum((load[ok]^2) * (se[ok]^2))) / sum(abs(load[ok]))
  z <- beta / se_beta
  data.frame(rg = beta, se = se_beta, z = z, p = 2 * pnorm(-abs(z)))
}

factor_alignment <- external_rg |>
  mutate(wm_id = as.character(wm_id)) |>
  left_join(loadings, by = "wm_id") |>
  group_by(external_trait) |>
  group_modify(~bind_rows(
    cbind(factor = "general_WM", project_factor(.x$general_WM_factor, .x$rg, .x$se)),
    cbind(factor = "FA_enriched", project_factor(.x$FA_specific_factor, .x$rg, .x$se)),
    cbind(factor = "MD_enriched", project_factor(.x$MD_specific_factor, .x$rg, .x$se))
  )) |>
  ungroup() |>
  group_by(external_trait) |>
  mutate(q = p.adjust(p, "BH")) |>
  ungroup()

domain_tests <- integrated |>
  left_join(wm_meta, by = "wm_id") |>
  group_by(domain) |>
  summarise(n_rows = n(), n_bmi = sum(grepl("BMI", integrated_stratum)),
            n_pain = sum(grepl("pain", integrated_stratum, ignore.case = TRUE)),
            n_mhc = sum(grepl("MHC|immune", integrated_stratum)), .groups = "drop")

metric_direction <- external_rg |>
  filter(external_trait == "KOA") |>
  left_join(wm_meta, by = "wm_id") |>
  group_by(metric) |>
  summarise(n = n(), mean_rg = mean(rg, na.rm = TRUE), n_positive = sum(rg > 0, na.rm = TRUE),
            n_negative = sum(rg < 0, na.rm = TRUE), .groups = "drop")
fa_md_test <- external_rg |> filter(external_trait == "KOA") |> left_join(wm_meta, by = "wm_id") |> summarise(wilcox_p = wilcox.test(rg ~ metric)$p.value)

mechanism <- integrated |>
  mutate(mechanism_axis = case_when(grepl("BMI_and_pain", integrated_stratum) ~ "BMI_and_pain_sensitive",
                                    grepl("BMI", integrated_stratum) ~ "BMI_or_metabolic_sensitive",
                                    grepl("MHC|immune", integrated_stratum) ~ "MHC_or_immune_sensitive",
                                    grepl("pain", integrated_stratum, ignore.case = TRUE) ~ "pain_sensitive",
                                    TRUE ~ "weak_or_unclassified")) |>
  count(mechanism_axis, name = "n_rows")

write_tsv(wm_meta |> filter(wm_id %in% traits), file.path(tables_dir, "Supplementary_Table_61_WM_factor_trait_manifest.tsv"))
write_tsv(wm_rg, file.path(tables_dir, "Supplementary_Table_62_pairwise_WM_LDSC_rg.tsv"))
write_tsv(external_rg, file.path(tables_dir, "Supplementary_Table_63_external_trait_to_WM_rg.tsv"))
write_tsv(loadings, file.path(tables_dir, "Supplementary_Table_64_WM_genetic_factor_loadings.tsv"))
write_tsv(factor_alignment, file.path(tables_dir, "Supplementary_Table_65_external_trait_factor_alignment.tsv"))
write_tsv(factor_alignment |> filter(external_trait %in% c("KOA", "BMI", "FinnGen_pain")), file.path(tables_dir, "Supplementary_Table_66_factor_interpretation_summary.tsv"))
write_tsv(domain_tests, file.path(tables_dir, "Supplementary_Table_67_domain_level_Fisher_tests.tsv"))
write_tsv(metric_direction, file.path(tables_dir, "Supplementary_Table_68_metric_level_Fisher_tests.tsv"))
write_tsv(bind_cols(metric_direction, fa_md_test), file.path(tables_dir, "Supplementary_Table_69_permutation_and_FA_MD_directionality.tsv"))
write_tsv(mechanism, file.path(tables_dir, "Supplementary_Table_70_mechanism_axis_classification.tsv"))

p39 <- factor_alignment |> ggplot(aes(x = factor, y = rg, fill = external_trait)) +
  geom_col(position = "dodge") + theme_bw(base_size = 9) + labs(x = NULL, y = "Projected genetic correlation", fill = "Trait")
ggsave(file.path(figures_dir, "Supplementary_Figure_39_Stage6_factor_external_rg.png"), p39, width = 7, height = 4, dpi = 300)

p40 <- metric_direction |> ggplot(aes(x = metric, y = mean_rg, fill = metric)) +
  geom_col() + theme_bw(base_size = 9) + labs(x = NULL, y = "Mean KOA-WM rg", fill = NULL)
ggsave(file.path(figures_dir, "Supplementary_Figure_40_Stage6_domain_evidence_score.png"), p40, width = 5, height = 4, dpi = 300)

p41 <- mechanism |> ggplot(aes(x = reorder(mechanism_axis, n_rows), y = n_rows, fill = mechanism_axis)) +
  geom_col() + coord_flip() + theme_bw(base_size = 9) + labs(x = NULL, y = "Locus-trait rows", fill = NULL)
ggsave(file.path(figures_dir, "Supplementary_Figure_41_Stage6_mechanism_axis_by_domain.png"), p41, width = 7, height = 4, dpi = 300)

write_tsv(as.data.frame(mat), file.path(tables_dir, "Supplementary_Figure_42_WM_rg_heatmap_matrix.tsv"))

p43 <- loadings |> ggplot(aes(x = reorder(paste0(tract_abbr, "_", metric), general_WM_factor), y = general_WM_factor, fill = metric)) +
  geom_col() + coord_flip() + theme_bw(base_size = 8) + labs(x = NULL, y = "General WM factor loading", fill = "Metric")
ggsave(file.path(figures_dir, "Supplementary_Figure_43_WM_factor_loadings.png"), p43, width = 7, height = 6, dpi = 300)

message("Completed Stage6 WM genetic factor and domain convergence analysis.")
