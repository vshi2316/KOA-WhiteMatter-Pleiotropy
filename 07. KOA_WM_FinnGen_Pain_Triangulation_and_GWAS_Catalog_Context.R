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

koa_wm <- read_tsv(Sys.getenv("KOA_WM_LOCI", file.path(data_dir, "triangulation_koa_wm_loci.tsv"))) |> mutate(set = "KOA-WM")
koa_pain <- read_tsv(Sys.getenv("KOA_PAIN_LOCI", file.path(data_dir, "triangulation_koa_pain_loci.tsv"))) |> mutate(set = "KOA-Pain")
pain_wm <- read_tsv(Sys.getenv("PAIN_WM_LOCI", file.path(data_dir, "triangulation_pain_wm_loci.tsv"))) |> mutate(set = "Pain-WM")
catalog <- read_tsv(Sys.getenv("KOA_WM_GWAS_CATALOG_CONTEXT", file.path(data_dir, "gwas_catalog_locus_context.tsv")))

overlap_count <- function(a, b, window = 500000) {
  sum(vapply(seq_len(nrow(a)), function(i) any(as.character(b$CHR) == as.character(a$CHR[i]) & abs(b$BP - a$BP[i]) <= window), logical(1)))
}

triangulation <- data.frame(
  comparison = c("KOA-WM_vs_KOA-Pain", "KOA-WM_vs_Pain-WM", "KOA-Pain_vs_Pain-WM", "Three-way"),
  n_overlap = c(overlap_count(koa_wm, koa_pain), overlap_count(koa_wm, pain_wm), overlap_count(koa_pain, pain_wm), 0)
)

density <- catalog |>
  group_by(locus_id) |>
  summarise(n_catalog_records = n(),
            has_oa_or_bone = any(grepl("osteoarthritis|bone|cartilage", trait, ignore.case = TRUE)),
            has_wm_or_brain = any(grepl("white matter|brain|diffusion|MRI", trait, ignore.case = TRUE)),
            has_pain = any(grepl("pain", trait, ignore.case = TRUE)), .groups = "drop") |>
  mutate(density_class = cut(n_catalog_records, c(-Inf, 0, 5, 20, Inf),
                             labels = c("none", "low_to_moderate", "moderate", "high")))

pleiotropic_traits <- catalog |>
  filter(!grepl("osteoarthritis|pain|white matter|brain|diffusion|MRI", trait, ignore.case = TRUE)) |>
  count(trait, sort = TRUE, name = "n_loci")

pathway_context <- read_tsv(Sys.getenv("KOA_WM_GWAS_CATALOG_PATHWAY_CONTEXT", file.path(data_dir, "gwas_catalog_pleiotropy_pathways.tsv")))

write_tsv(triangulation, file.path(tables_dir, "Supplementary_Table_23_FinnGen_pain_triangulation.tsv"))
write_tsv(density, file.path(tables_dir, "Supplementary_Table_24_GWAS_Catalog_locus_context.tsv"))
write_tsv(density |> count(density_class), file.path(tables_dir, "Supplementary_Table_25_GWAS_Catalog_density_classification.tsv"))
write_tsv(pleiotropic_traits, file.path(tables_dir, "Supplementary_Table_26_pleiotropic_traits.tsv"))
write_tsv(pathway_context, file.path(tables_dir, "Supplementary_Table_27_gProfiler_pleiotropy_pathways.tsv"))
write_tsv(catalog, file.path(tables_dir, "Supplementary_Table_28_locus_pleiotropy_distribution.tsv"))

p <- density |>
  summarise(oa_or_bone = mean(has_oa_or_bone) * 100, wm_or_brain = mean(has_wm_or_brain) * 100, pain = mean(has_pain) * 100) |>
  pivot_longer(everything(), names_to = "category", values_to = "percent") |>
  ggplot(aes(x = category, y = percent, fill = category)) +
  geom_col() + theme_bw(base_size = 9) + labs(x = NULL, y = "Loci with GWAS Catalog context (%)", fill = NULL)
ggsave(file.path(figures_dir, "Supplementary_Figure_26_GWAS_Catalog_context_source.png"), p, width = 5, height = 4, dpi = 300)

message("Completed FinnGen pain triangulation and GWAS Catalog context analysis.")
