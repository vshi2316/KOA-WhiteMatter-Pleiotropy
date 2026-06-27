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

parse_ldsc_log <- function(path) {
  txt <- if (file.exists(path)) readLines(path, warn = FALSE) else character()
  grab <- function(pattern) {
    z <- grep(pattern, txt, value = TRUE)
    if (!length(z)) return(NA_real_)
    suppressWarnings(as.numeric(sub(".*?([-+]?[0-9.]+(?:e[-+]?[0-9]+)?).*", "\\1", z[1], ignore.case = TRUE)))
  }
  data.frame(rg = grab("Genetic Correlation"), se = grab("Standard Error"), z = grab("^Z-score"), p = grab("^P:"), stringsAsFactors = FALSE)
}

wm_meta <- read_tsv(file.path(project_dir, "data_manifest", "wm_trait_metadata.tsv"))
ldsc_log_dir <- Sys.getenv("KOA_WM_LDSC_LOG_DIR", file.path(data_dir, "ldsc_koa_wm_logs"))
hdl_file <- Sys.getenv("KOA_WM_HDL_RESULTS", file.path(data_dir, "hdl_koa_wm_results.tsv"))

ldsc <- lapply(seq_len(nrow(wm_meta)), function(i) {
  wm_id <- wm_meta$wm_id[i]
  log_file <- file.path(ldsc_log_dir, paste0("KOA_vs_WM_", wm_id, ".log"))
  cbind(wm_meta[i, ], parse_ldsc_log(log_file), method = "LDSC")
}) |> bind_rows()

hdl <- if (file.exists(hdl_file)) {
  read_tsv(hdl_file) |> mutate(wm_id = as.character(wm_id), method = "HDL") |> left_join(wm_meta, by = "wm_id")
} else {
  data.frame()
}

global_rg <- bind_rows(ldsc, hdl) |>
  mutate(p = as.numeric(p), rg = as.numeric(rg),
         q = ave(p, method, FUN = function(x) p.adjust(x, method = "BH")),
         nominal = is.finite(p) & p < 0.05,
         fdr_significant = is.finite(q) & q < 0.05,
         direction = ifelse(rg > 0, "positive", ifelse(rg < 0, "negative", "zero"))) |>
  arrange(method, metric, tract_abbr)

direction_balance <- global_rg |>
  group_by(method) |>
  summarise(n_traits = n(), n_positive = sum(direction == "positive", na.rm = TRUE),
            n_negative = sum(direction == "negative", na.rm = TRUE),
            mean_positive_rg = mean(rg[direction == "positive"], na.rm = TRUE),
            mean_negative_rg = mean(rg[direction == "negative"], na.rm = TRUE),
            direction_balance_index = abs(n_positive - n_negative) / n_traits * 100,
            .groups = "drop")

write_tsv(wm_meta, file.path(tables_dir, "Supplementary_Table_1_WM_trait_metadata.tsv"))
write_tsv(global_rg, file.path(tables_dir, "Supplementary_Table_2_Global_LDSC_HDL_rg.tsv"))
write_tsv(direction_balance, file.path(tables_dir, "Figure_1A_direction_balance_source.tsv"))

p <- global_rg |> mutate(label = paste0(tract_abbr, "_", metric)) |>
  ggplot(aes(x = reorder(label, rg), y = rg, fill = nominal)) +
  geom_col(width = 0.75) + coord_flip() + facet_wrap(~method, scales = "free_y") +
  scale_fill_manual(values = c("FALSE" = "grey75", "TRUE" = "#B03A2E")) +
  theme_bw(base_size = 9) + labs(x = NULL, y = "Genetic correlation (rg)", fill = "P < 0.05")
ggsave(file.path(figures_dir, "Figure_1A_global_rg_source.png"), p, width = 8, height = 8, dpi = 300)

message("Completed KOA-WM global LDSC/HDL genetic correlation analysis.")
