#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(limma)
  library(edgeR)
  library(lme4)
  library(lmerTest)
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

clinical <- read_tsv(Sys.getenv("KOA_WM_CLINICAL_MULTIOMICS_PHENO", file.path(data_dir, "clinical_multiomics_phenotypes.tsv")))
mirna_counts <- read_tsv(Sys.getenv("KOA_WM_CLINICAL_MULTIOMICS_MIRNA", file.path(data_dir, "clinical_multiomics_mirna_counts.tsv")))
mirdb <- read_tsv(Sys.getenv("KOA_WM_MIRDB_TARGETS", file.path(data_dir, "mirdb_targets.tsv")))
gwas_genes <- read_tsv(file.path(tables_dir, "Figure_4C_gene_mapping_source.tsv")) |> distinct(gene)

cluster_vars <- c("WOMAC_pain", "WOMAC_function", "age", "BMI", "KL_grade")
scaled <- clinical
scaled[cluster_vars] <- scale(scaled[cluster_vars])
set.seed(20260627)
clinical$clinical_endotype <- paste0("Cluster_", kmeans(scaled[cluster_vars], centers = 3, nstart = 100)$cluster - 1)

long_counts <- mirna_counts |>
  pivot_longer(-sample_id, names_to = "miRNA", values_to = "count") |>
  left_join(clinical |> select(sample_id, patient_id, clinical_endotype, tissue, WOMAC_pain, pain_improvement, age, sex, BMI, KL_grade), by = "sample_id")

de_results <- long_counts |>
  group_by(tissue, miRNA) |>
  summarise(mean_count = mean(count, na.rm = TRUE),
            log2FC = mean(log2(count + 1)[clinical_endotype == "Cluster_1"], na.rm = TRUE) -
              mean(log2(count + 1)[clinical_endotype != "Cluster_1"], na.rm = TRUE),
            p = tryCatch(t.test(log2(count + 1) ~ clinical_endotype)$p.value, error = function(e) NA_real_),
            .groups = "drop") |>
  group_by(tissue) |>
  mutate(q = p.adjust(p, method = "BH")) |>
  ungroup()

targets <- de_results |>
  filter(p < 0.05) |>
  inner_join(mirdb |> filter(prediction_score >= 80), by = "miRNA")

universe <- unique(mirdb$target_gene)
overlap_tests <- targets |>
  group_by(tissue, miRNA) |>
  summarise(target_set = list(unique(target_gene)), .groups = "drop") |>
  rowwise() |>
  mutate(overlap = length(intersect(target_set, gwas_genes$gene)),
         target_n = length(target_set),
         gwas_n = nrow(gwas_genes),
         universe_n = length(universe),
         fisher_p = fisher.test(matrix(c(overlap, target_n - overlap, gwas_n - overlap,
                                         universe_n - target_n - gwas_n + overlap), nrow = 2))$p.value,
         odds_ratio = unname(fisher.test(matrix(c(overlap, target_n - overlap, gwas_n - overlap,
                                                  universe_n - target_n - gwas_n + overlap), nrow = 2))$estimate)) |>
  ungroup()

pain_models <- long_counts |>
  group_by(miRNA) |>
  do({
    dat <- .
    fit <- tryCatch(lmer(log2(count + 1) ~ WOMAC_pain + age + sex + BMI + KL_grade + (1 | patient_id), data = dat), error = function(e) NULL)
    if (is.null(fit)) data.frame(beta = NA_real_, p = NA_real_) else {
      co <- summary(fit)$coefficients
      data.frame(beta = co["WOMAC_pain", "Estimate"], p = co["WOMAC_pain", "Pr(>|t|)"])
    }
  }) |>
  ungroup()

prioritized <- de_results |>
  inner_join(pain_models, by = "miRNA", suffix = c("_de", "_pain")) |>
  filter(p_de < 0.05, p_pain < 0.05) |>
  arrange(p_de, p_pain)

write_tsv(de_results, file.path(tables_dir, "Supplementary_Table_29_differentially_expressed_miRNAs.tsv"))
write_tsv(targets, file.path(tables_dir, "Supplementary_Table_30_miRDB_high_confidence_targets.tsv"))
write_tsv(overlap_tests, file.path(tables_dir, "Supplementary_Table_31_miRNA_GWAS_gene_overlap_tests.tsv"))
write_tsv(prioritized, file.path(tables_dir, "Supplementary_Table_32_prioritized_miRNAs.tsv"))
write_tsv(pain_models, file.path(tables_dir, "Supplementary_Table_33_miRNA_pain_mixed_models.tsv"))

p <- prioritized |> head(20) |>
  ggplot(aes(x = reorder(miRNA, abs(beta)), y = beta)) +
  geom_col(fill = "#4D9221") + coord_flip() + theme_bw(base_size = 9) +
  labs(x = NULL, y = "Pain association beta")
ggsave(file.path(figures_dir, "Supplementary_Figure_27_miRNA_context_source.png"), p, width = 6, height = 5, dpi = 300)

message("Completed exploratory clinical multi-omics miRNA context analysis.")
