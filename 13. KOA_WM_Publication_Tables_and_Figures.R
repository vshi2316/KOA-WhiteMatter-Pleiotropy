#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

project_dir <- Sys.getenv("KOA_WM_PROJECT_DIR", ".")
results_dir <- Sys.getenv("KOA_WM_RESULTS_DIR", file.path(project_dir, "results"))
tables_dir <- file.path(results_dir, "tables")
figures_dir <- file.path(results_dir, "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) data.table::fread(path, data.table = FALSE)
write_tsv <- function(x, path) data.table::fwrite(as.data.frame(x), path, sep = "\t", na = "NA")

inventory <- data.frame(expected_table = sprintf("Supplementary_Table_%d", 1:70))
available <- list.files(tables_dir, pattern = "^Supplementary_Table_.*\\.tsv$", full.names = FALSE)
inventory$present <- vapply(inventory$expected_table, function(x) any(startsWith(available, paste0(x, "_"))), logical(1))
write_tsv(inventory, file.path(tables_dir, "publication_table_inventory.tsv"))

table1 <- data.frame(
  analysis_layer = c("Genome-wide genetic correlation", "Local genetic correlation", "Variant-level enrichment",
                     "External imaging consistency", "Pain triangulation", "Functional annotation",
                     "Causal-sensitivity analysis", "OA subtype and external-trait stratification",
                     "WM factor and domain-level convergence"),
  evidence_constraint = c("No FDR-significant global KOA-WM rg",
                          "Sparse regional signals with tract-level heterogeneity",
                          "condFDR/conjFDR records represent P-value enrichment",
                          "Limited ENIGMA consistency for most tract signals",
                          "No three-way KOA-WM-pain locus overlap",
                          "Few mapped genes and hypothesis-generating enrichment",
                          "Low PP.H4, sharing-favored CAUSE patterns, and weak MVMR KOA instruments",
                          "Most apparent KOA-enriched loci reassigned to BMI, pain, or immune-sensitive contexts",
                          "Descriptive general and MD-enriched WM factor alignment with preselection caveat")
)
write_tsv(table1, file.path(tables_dir, "Table_1_evidence_tiered_interpretation.tsv"))

if (file.exists(file.path(tables_dir, "Supplementary_Table_56_integrated_locus_evidence_strata.tsv"))) {
  stage5 <- read_tsv(file.path(tables_dir, "Supplementary_Table_56_integrated_locus_evidence_strata.tsv"))
  p <- stage5 |> count(integrated_stratum, name = "n_locus_trait_rows") |>
    ggplot(aes(x = reorder(integrated_stratum, n_locus_trait_rows), y = n_locus_trait_rows, fill = integrated_stratum)) +
    geom_col() + coord_flip() + theme_bw(base_size = 9) + labs(x = NULL, y = "Locus-trait rows", fill = NULL)
  ggsave(file.path(figures_dir, "Publication_summary_Stage5_counts.png"), p, width = 7, height = 4, dpi = 300)
}

if (file.exists(file.path(tables_dir, "Supplementary_Table_65_external_trait_factor_alignment.tsv"))) {
  stage6 <- read_tsv(file.path(tables_dir, "Supplementary_Table_65_external_trait_factor_alignment.tsv"))
  p <- stage6 |> filter(external_trait == "KOA") |>
    ggplot(aes(x = reorder(factor, rg), y = rg, fill = q < 0.05)) +
    geom_col() + coord_flip() + theme_bw(base_size = 9) + labs(x = NULL, y = "KOA projected rg", fill = "FDR q < 0.05")
  ggsave(file.path(figures_dir, "Publication_summary_Stage6_KOA_factor_alignment.png"), p, width = 6, height = 4, dpi = 300)
}

message("Completed publication-ready table and figure assembly.")
