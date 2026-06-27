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

wm_meta <- read_tsv(file.path(project_dir, "data_manifest", "wm_trait_metadata.tsv"))
lava_results <- read_tsv(Sys.getenv("KOA_WM_LAVA_RESULTS", file.path(data_dir, "lava_koa_wm_results.tsv"))) |>
  mutate(wm_id = as.character(wm_id), p = as.numeric(p), rho = as.numeric(rho)) |>
  left_join(wm_meta, by = "wm_id") |>
  mutate(q = p.adjust(p, method = "BH"), nominal = p < 0.05, fdr_significant = q < 0.05)

write_tsv(lava_results |> filter(nominal), file.path(tables_dir, "Supplementary_Table_3_LAVA_nominal_local_rg.tsv"))
write_tsv(lava_results |> filter(metric == "FA", fdr_significant), file.path(tables_dir, "Supplementary_Table_4_LAVA_FA_FDR_significant.tsv"))
write_tsv(lava_results |> filter(metric == "MD", fdr_significant), file.path(tables_dir, "Supplementary_Table_5_LAVA_MD_FDR_significant.tsv"))

p1 <- lava_results |>
  ggplot(aes(x = genomic_midpoint, y = -log10(p), color = rho)) +
  geom_point(size = 1.1, alpha = 0.8) + facet_wrap(~metric, ncol = 1) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "red") +
  scale_color_gradient2(low = "#2166AC", mid = "grey85", high = "#B2182B", midpoint = 0) +
  theme_bw(base_size = 9) + labs(x = "Genomic position", y = "-log10(P)", color = "rho")
ggsave(file.path(figures_dir, "Figure_1B_LAVA_local_rg_scatter_source.png"), p1, width = 8, height = 5, dpi = 300)

p2 <- lava_results |>
  ggplot(aes(x = rho, y = -log10(p), color = fdr_significant)) +
  geom_point(size = 1.3, alpha = 0.85) + facet_wrap(~metric) +
  scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "#B03A2E")) +
  theme_bw(base_size = 9) + labs(x = "Local genetic correlation (rho)", y = "-log10(P)", color = "FDR q < 0.05")
ggsave(file.path(figures_dir, "Figure_1C_LAVA_volcano_source.png"), p2, width = 8, height = 4, dpi = 300)

message("Completed KOA-WM LAVA local genetic correlation analysis.")
