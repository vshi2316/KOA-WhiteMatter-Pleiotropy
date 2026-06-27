#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

project_dir <- Sys.getenv("KOA_WM_PROJECT_DIR", ".")
data_dir <- Sys.getenv("KOA_WM_DATA_DIR", file.path(project_dir, "data"))
results_dir <- Sys.getenv("KOA_WM_RESULTS_DIR", file.path(project_dir, "results"))
tables_dir <- file.path(results_dir, "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) data.table::fread(path, data.table = FALSE)
write_tsv <- function(x, path) data.table::fwrite(as.data.frame(x), path, sep = "\t", na = "NA")

primary <- read_tsv(Sys.getenv("KOA_WM_RETAINED_LOCI", file.path(data_dir, "retained_koa_wm_loci_for_enigma.tsv")))
enigma <- read_tsv(Sys.getenv("KOA_WM_ENIGMA_HITS", file.path(data_dir, "enigma_dti_significant_loci.tsv")))
window_bp <- 500000

overlap <- primary
overlap$n_enigma_overlap <- vapply(seq_len(nrow(primary)), function(i) {
  nrow(enigma |>
         filter(as.character(CHR) == as.character(primary$CHR[i]),
                BP >= primary$BP[i] - window_bp,
                BP <= primary$BP[i] + window_bp,
                wm_label == primary$wm_label[i]))
}, integer(1))
overlap <- overlap |>
  mutate(overlap_rate = n_enigma_overlap / pmax(n_primary_loci, 1) * 100,
         consistency_class = case_when(overlap_rate >= 70 ~ "substantial",
                                       overlap_rate >= 45 ~ "intermediate",
                                       overlap_rate >= 20 ~ "moderate",
                                       TRUE ~ "limited"))

summary <- overlap |> summarise(n_loci = n(), mean_overlap_rate = mean(overlap_rate, na.rm = TRUE),
                                n_substantial = sum(consistency_class == "substantial"),
                                n_retained_for_triangulation = sum(overlap_rate >= 20))

write_tsv(summary, file.path(tables_dir, "Supplementary_Table_17_ENIGMA_overlap_summary.tsv"))
write_tsv(overlap |> filter(consistency_class == "substantial"), file.path(tables_dir, "Supplementary_Table_18_ENIGMA_substantial_overlap.tsv"))
write_tsv(overlap |> filter(wm_label == "CCG-R_FA"), file.path(tables_dir, "Supplementary_Table_19_ENIGMA_CCG_R_FA.tsv"))
write_tsv(overlap |> filter(wm_label == "SCP-L_FA"), file.path(tables_dir, "Supplementary_Table_20_ENIGMA_SCP_L_FA.tsv"))
write_tsv(overlap |> filter(wm_label == "SLF-L_MD"), file.path(tables_dir, "Supplementary_Table_21_ENIGMA_SLF_L_MD.tsv"))
write_tsv(overlap |> filter(wm_label == "SLF-R_MD"), file.path(tables_dir, "Supplementary_Table_22_ENIGMA_SLF_R_MD.tsv"))

message("Completed ENIGMA external consistency analysis.")
