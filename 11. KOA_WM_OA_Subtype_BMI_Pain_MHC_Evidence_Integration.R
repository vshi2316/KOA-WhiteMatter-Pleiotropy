#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
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

subtype_manifest <- read_tsv(Sys.getenv("KOA_WM_OA_SUBTYPE_MANIFEST", file.path(data_dir, "oa_subtype_source_file_manifest.tsv")))
subtype_h2 <- read_tsv(Sys.getenv("KOA_WM_OA_SUBTYPE_H2", file.path(data_dir, "oa_subtype_ldsc_h2_qc.tsv")))
subtype_lookup <- read_tsv(Sys.getenv("KOA_WM_OA_SUBTYPE_LOOKUP", file.path(data_dir, "oa_subtype_locus_lookup.tsv")))
index_loci <- read_tsv(file.path(tables_dir, "Supplementary_Table_13_distinct_cross_trait_loci.tsv"))
external <- read_tsv(file.path(tables_dir, "Supplementary_Table_38_external_trait_locus_sensitivity.tsv"))
coloc <- read_tsv(file.path(tables_dir, "Supplementary_Table_37_local_coloc_summary.tsv"))
gsmr <- read_tsv(file.path(tables_dir, "Supplementary_Table_50_GSMR2_network_effect_estimates.tsv"))

subtype_class <- subtype_lookup |>
  mutate(hip_hit = hip_oa_p < 0.05, hand_hit = hand_oa_p < 0.05, all_hit = all_oa_p < 0.05,
         subtype_class = case_when(koa_p < 0.05 & !hip_hit & !hand_hit & !all_hit ~ "KOA_enriched",
                                   koa_p < 0.05 & hip_hit & !hand_hit ~ "lower_limb_OA",
                                   hand_hit & all_hit ~ "hand_or_general_OA",
                                   all_hit & !hip_hit & !hand_hit ~ "general_OA_only",
                                   koa_p < 0.05 & hip_hit & hand_hit & all_hit ~ "pan_OA",
                                   TRUE ~ "weak_or_unclassified"))

integrated <- index_loci |>
  left_join(subtype_class |> select(locus_id, subtype_class), by = "locus_id") |>
  left_join(external, by = "locus_id") |>
  left_join(coloc |> select(locus_id, PP.H4.abf, coloc_class), by = "locus_id") |>
  mutate(bmi_sensitive = bmi_signal %in% TRUE,
         pain_sensitive = pain_signal %in% TRUE,
         depression_sensitive = depression_signal %in% TRUE,
         mhc_sensitive = is_mhc %in% TRUE,
         integrated_stratum = case_when(mhc_sensitive ~ "MHC_or_immune_region_pleiotropy",
                                        bmi_sensitive & pain_sensitive ~ "BMI_and_pain_sensitive_pleiotropy",
                                        bmi_sensitive ~ "BMI_or_metabolic_sensitive_pleiotropy",
                                        pain_sensitive ~ "pain_sensitive_pleiotropy",
                                        subtype_class == "KOA_enriched" ~ "residual_KOA_enriched_noncolocalized",
                                        TRUE ~ "weak_or_unclassified"))

trait_summary <- integrated |>
  group_by(wm_label, metric, domain, integrated_stratum) |>
  summarise(n_loci = n(), .groups = "drop") |>
  group_by(wm_label) |>
  mutate(primary_trait_stratum = integrated_stratum[which.max(n_loci)]) |>
  ungroup()

locus_counts <- integrated |> count(integrated_stratum, name = "n_loci") |> arrange(desc(n_loci))
trait_counts <- trait_summary |> count(primary_trait_stratum, name = "n_traits") |> arrange(desc(n_traits))
gsmr_key <- gsmr |> filter(grepl("BMI|pain|KOA", exposure, ignore.case = TRUE) | grepl("BMI|pain|KOA", outcome, ignore.case = TRUE))

write_tsv(subtype_manifest, file.path(tables_dir, "Supplementary_Table_52_OA_subtype_source_manifest.tsv"))
write_tsv(subtype_h2, file.path(tables_dir, "Supplementary_Table_53_OA_subtype_LDSC_h2_QC.tsv"))
write_tsv(subtype_lookup, file.path(tables_dir, "Supplementary_Table_54_OA_subtype_locus_lookup.tsv"))
write_tsv(subtype_class, file.path(tables_dir, "Supplementary_Table_55_OA_subtype_specificity_classification.tsv"))
write_tsv(integrated, file.path(tables_dir, "Supplementary_Table_56_integrated_locus_evidence_strata.tsv"))
write_tsv(trait_summary, file.path(tables_dir, "Supplementary_Table_57_trait_level_evidence_summary.tsv"))
write_tsv(locus_counts, file.path(tables_dir, "Supplementary_Table_58_locus_stratum_counts.tsv"))
write_tsv(trait_counts, file.path(tables_dir, "Supplementary_Table_59_trait_stratum_counts.tsv"))
write_tsv(gsmr_key, file.path(tables_dir, "Supplementary_Table_60_key_GSMR2_edges_for_Stage5.tsv"))

p1 <- integrated |> count(wm_label, integrated_stratum) |> pivot_wider(names_from = integrated_stratum, values_from = n, values_fill = 0) |>
  pivot_longer(-wm_label, names_to = "stratum", values_to = "n") |>
  ggplot(aes(x = stratum, y = wm_label, fill = n)) + geom_tile(color = "white") +
  theme_bw(base_size = 8) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, fill = "Loci")
ggsave(file.path(figures_dir, "Supplementary_Figure_37_Stage5_integrated_evidence_heatmap.png"), p1, width = 8, height = 6, dpi = 300)

p2 <- locus_counts |> ggplot(aes(x = reorder(integrated_stratum, n_loci), y = n_loci, fill = integrated_stratum)) +
  geom_col() + coord_flip() + theme_bw(base_size = 9) + labs(x = NULL, y = "Number of loci", fill = NULL)
ggsave(file.path(figures_dir, "Supplementary_Figure_38_Stage5_locus_evidence_stratum_counts.png"), p2, width = 7, height = 4, dpi = 300)

message("Completed Stage5 OA subtype, BMI, pain, and MHC evidence integration.")
