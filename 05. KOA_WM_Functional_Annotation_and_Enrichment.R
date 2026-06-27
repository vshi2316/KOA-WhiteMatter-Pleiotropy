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

fuma_dir <- Sys.getenv("KOA_WM_FUMA_DIR", file.path(data_dir, "fuma_outputs"))
gprofiler_dir <- Sys.getenv("KOA_WM_GPROFILER_DIR", file.path(data_dir, "gprofiler_outputs"))

regions <- read_tsv(file.path(tables_dir, "Supplementary_Table_13_distinct_cross_trait_loci.tsv"))
snps <- read_tsv(file.path(fuma_dir, "candidate_and_ld_proxy_snps.tsv"))
genes <- read_tsv(file.path(fuma_dir, "mapped_genes.tsv"))
enrichment <- read_tsv(file.path(gprofiler_dir, "gprofiler_enrichment.tsv"))

annotation <- snps |>
  mutate(metric = ifelse(grepl("FA", phenotype, ignore.case = TRUE), "FA", "MD"),
         cadd_high = as.numeric(CADD) > 12.37,
         regulomedb_high = suppressWarnings(as.numeric(RegulomeDB)) < 2,
         open_chromatin = suppressWarnings(as.numeric(chromatin_state)) <= 7) |>
  group_by(metric, functional_category) |>
  summarise(n_snps = n(), n_cadd_high = sum(cadd_high, na.rm = TRUE),
            n_regulomedb_high = sum(regulomedb_high, na.rm = TRUE),
            n_open_chromatin = sum(open_chromatin, na.rm = TRUE), .groups = "drop")

gene_summary <- genes |> count(metric, gene, mapping_strategy, name = "n_snps")
enrichment_clean <- enrichment |>
  mutate(q = as.numeric(if ("adjusted_p_value" %in% names(.)) adjusted_p_value else p_value)) |>
  arrange(q)

write_tsv(regions, file.path(tables_dir, "Supplementary_Table_14_candidate_shared_regions.tsv"))
write_tsv(snps, file.path(tables_dir, "Supplementary_Table_15_candidate_and_LD_linked_SNPs.tsv"))
write_tsv(annotation, file.path(tables_dir, "Supplementary_Table_16_functional_annotation_summary.tsv"))
write_tsv(gene_summary, file.path(tables_dir, "Figure_4C_gene_mapping_source.tsv"))
write_tsv(enrichment_clean, file.path(tables_dir, "Figure_4D_gprofiler_enrichment_source.tsv"))

p <- annotation |>
  ggplot(aes(x = functional_category, y = n_snps, fill = metric)) +
  geom_col(position = "dodge") + coord_flip() +
  theme_bw(base_size = 9) + labs(x = NULL, y = "SNP count", fill = "Metric")
ggsave(file.path(figures_dir, "Figure_4A_B_functional_annotation_source.png"), p, width = 7, height = 5, dpi = 300)

message("Completed functional annotation and enrichment analysis.")
