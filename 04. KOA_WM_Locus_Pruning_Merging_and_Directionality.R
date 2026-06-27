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
pick_col <- function(x, candidates, default = NA) {
  hit <- intersect(candidates, names(x))
  if (length(hit)) x[[hit[1]]] else default
}

wm_meta <- read_tsv(file.path(project_dir, "data_manifest", "wm_trait_metadata.tsv"))
conj <- read_tsv(file.path(tables_dir, "conjFDR_SNP_level_enrichment_records.tsv"))
conj$wm_id <- as.character(pick_col(conj, "wm_id", sub(".*_(\\d+)$", "\\1", conj$pair_id)))
conj <- conj |>
  mutate(CHR = as.character(CHR), BP = as.numeric(BP)) |>
  left_join(wm_meta, by = "wm_id")

lead_by_region <- conj |>
  arrange(CHR, wm_label, BP, conjFDR) |>
  group_by(CHR, wm_label) |>
  mutate(region_index = cumsum(c(TRUE, diff(BP) > 500000))) |>
  group_by(CHR, wm_label, region_index) |>
  slice_min(conjFDR, n = 1, with_ties = FALSE) |>
  ungroup()

merged_loci <- lead_by_region |>
  arrange(CHR, BP) |>
  mutate(locus_id = paste0("L", row_number()), locus_start = pmax(1, BP - 500000), locus_end = BP + 500000))

beta_koa <- pick_col(lead_by_region, c("BETA_KOA", "Z_KOA"), NA_real_)
beta_wm <- pick_col(lead_by_region, c("BETA_WM", "Z_WM"), NA_real_)
directionality <- lead_by_region |>
  mutate(beta_koa = as.numeric(beta_koa), beta_wm = as.numeric(beta_wm),
         direction_class = case_when(
           is.na(beta_koa) | is.na(beta_wm) ~ "unknown",
           sign(beta_koa) == sign(beta_wm) ~ "concordant",
           TRUE ~ "opposite"
         )) |>
  count(metric, direction_class, name = "n_loci")

tract_counts <- conj |> count(metric, domain, tract_abbr, name = "n_records") |> arrange(metric, desc(n_records))
hemisphere <- conj |> count(metric, laterality, name = "n_records") |> group_by(metric) |> mutate(percent = n_records / sum(n_records) * 100) |> ungroup()

write_tsv(tract_counts, file.path(tables_dir, "Supplementary_Table_8_tract_category_counts.tsv"))
write_tsv(hemisphere, file.path(tables_dir, "Supplementary_Table_9_hemispheric_specificity.tsv"))
write_tsv(lead_by_region, file.path(tables_dir, "Supplementary_Table_10_top_lead_SNPs.tsv"))
write_tsv(directionality |> filter(metric == "FA"), file.path(tables_dir, "Supplementary_Table_11_KOA_FA_effect_direction.tsv"))
write_tsv(directionality |> filter(metric == "MD"), file.path(tables_dir, "Supplementary_Table_12_KOA_MD_effect_direction.tsv"))
write_tsv(merged_loci, file.path(tables_dir, "Supplementary_Table_13_distinct_cross_trait_loci.tsv"))

p1 <- tract_counts |>
  ggplot(aes(x = reorder(tract_abbr, n_records), y = n_records, fill = domain)) +
  geom_col() + coord_flip() + facet_wrap(~metric, scales = "free") +
  theme_bw(base_size = 9) + labs(x = NULL, y = "Number of enrichment records", fill = "Domain")
ggsave(file.path(figures_dir, "Figure_2G_I_tract_hemisphere_source.png"), p1, width = 7, height = 6, dpi = 300)

p2 <- merged_loci |>
  ggplot(aes(x = BP, y = -log10(conjFDR), color = wm_label)) +
  geom_point(size = 2) + facet_wrap(~CHR, scales = "free_x") +
  theme_bw(base_size = 8) + labs(x = "Genomic position", y = "-log10(conjFDR)", color = "WM trait")
ggsave(file.path(figures_dir, "Figure_3A_distinct_loci_source.png"), p2, width = 10, height = 5, dpi = 300)

message("Completed locus pruning, merging, tract summaries, and directionality analysis.")
