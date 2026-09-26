# ============================================================
# Stage 10: Automated Cell-Type Annotation with SingleR
# Dataset: GSE144735 - KUL3 colorectal cancer scRNA-seq
# Input: completed DEG-annotated Seurat object
# ============================================================

# ------------------------------------------------------------
# 1. Install required packages (only needed the first time)
# ------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "SummarizedExperiment"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

install.packages("rlang")
packageVersion("rlang")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install("celldex", ask = FALSE, update = FALSE)


BiocManager::install("celldex", ask = FALSE, update = FALSE)

BiocManager::install("celldex", ask = FALSE, update = FALSE)

# ------------------------------------------------------------
# 2. Load packages
# ------------------------------------------------------------

library(Seurat)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(dplyr)
library(ggplot2)

# ------------------------------------------------------------
# 3. Project folders
# ------------------------------------------------------------

# Main project folder
project_dir <- "D:/github/CRC_stemness_invasive_border_scrnaseq"

# Save all CSV and RDS result files here
results_dir <- file.path(project_dir, "results")

# Save figures here
figures_dir <- file.path(project_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------
# 4. Load the completed DEG-annotated object
# ------------------------------------------------------------

seurat_obj <- readRDS(
  "D:/github/CRC_stemness_invasive_border_scrnaseq/results/09_DEG_annotated_seurat_obj.rds"
)

DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "seurat_clusters"

print(table(seurat_obj$seurat_clusters))

# ------------------------------------------------------------
# 5. Convert Seurat object to SingleCellExperiment
# ------------------------------------------------------------

sce_obj <- as.SingleCellExperiment(
  seurat_obj,
  assay = "RNA"
)

# Ensure normalized log-expression data are available
if (!"logcounts" %in% assayNames(sce_obj)) {
  logcounts(sce_obj) <- log1p(counts(sce_obj))
}

# ------------------------------------------------------------
# 6. Load Human Primary Cell Atlas reference
# Downloads automatically the first time only
# ------------------------------------------------------------

hpca_ref <- HumanPrimaryCellAtlasData()

# Keep only genes shared by your data and the reference
common_genes <- intersect(
  rownames(sce_obj),
  rownames(hpca_ref)
)

sce_obj <- sce_obj[common_genes, ]
hpca_ref <- hpca_ref[common_genes, ]

message("Shared genes used for SingleR: ", length(common_genes))

# ------------------------------------------------------------
# 7. Run automated SingleR annotation by cluster
# ------------------------------------------------------------

singleR_result <- SingleR(
  test = sce_obj,
  ref = hpca_ref,
  labels = hpca_ref$label.main,
  clusters = as.character(seurat_obj$seurat_clusters)
)

# Save the complete SingleR result

saveRDS(
  singleR_result,
  file.path(results_dir, "10_SingleR_cluster_results.rds")
)

# ------------------------------------------------------------
# 8. Make a table of SingleR labels for each cluster
# ------------------------------------------------------------

singleR_cluster_labels <- data.frame(
  cluster = rownames(singleR_result),
  singleR_raw_label = singleR_result$labels,
  singleR_pruned_label = singleR_result$pruned.labels,
  stringsAsFactors = FALSE
)

# Use the pruned label if confident.
# Otherwise, mark the cluster as uncertain.
singleR_cluster_labels$singleR_annotation <-
  singleR_cluster_labels$singleR_pruned_label

singleR_cluster_labels$singleR_annotation[
  is.na(singleR_cluster_labels$singleR_annotation)
] <- "Ambiguous / review markers"

print(singleR_cluster_labels)

write.csv(
  singleR_cluster_labels,
  file.path(results_dir, "10_SingleR_cluster_labels.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 9. Add automated labels into Seurat metadata
# ------------------------------------------------------------

raw_label_lookup <- setNames(
  singleR_cluster_labels$singleR_raw_label,
  singleR_cluster_labels$cluster
)

pruned_label_lookup <- setNames(
  singleR_cluster_labels$singleR_pruned_label,
  singleR_cluster_labels$cluster
)

final_label_lookup <- setNames(
  singleR_cluster_labels$singleR_annotation,
  singleR_cluster_labels$cluster
)

cluster_ids <- as.character(seurat_obj$seurat_clusters)

seurat_obj$singleR_raw_label <- unname(
  raw_label_lookup[cluster_ids]
)

seurat_obj$singleR_pruned_label <- unname(
  pruned_label_lookup[cluster_ids]
)

seurat_obj$singleR_annotation <- unname(
  final_label_lookup[cluster_ids]
)

print(table(seurat_obj$singleR_annotation))

# ------------------------------------------------------------
# 10. Automated annotation UMAP
# ------------------------------------------------------------

singleR_umap <- DimPlot(
  object = seurat_obj,
  reduction = "umap_after",
  group.by = "singleR_annotation",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.25
) +
  ggtitle("Automated Cell-Type Annotation: SingleR")

print(singleR_umap)

ggsave(
  filename = file.path(
    figures_dir,
    "10_SingleR_automated_annotation_UMAP.png"
  ),
  plot = singleR_umap,
  width = 14,
  height = 8,
  dpi = 300
)

# ------------------------------------------------------------
# 11. Compare canonical, DEG, and automated annotations
# ------------------------------------------------------------

annotation_comparison <- data.frame(
  cell_barcode = colnames(seurat_obj),
  cluster = as.character(seurat_obj$seurat_clusters),
  canonical_annotation = seurat_obj$canonical_annotation,
  deg_annotation = seurat_obj$deg_annotation,
  singleR_annotation = seurat_obj$singleR_annotation,
  stringsAsFactors = FALSE
)

write.csv(
  annotation_comparison,
  file.path(results_dir, "10_annotation_comparison_per_cell.csv"),
  row.names = FALSE
)

cluster_comparison <- annotation_comparison %>%
  group_by(
    cluster,
    canonical_annotation,
    deg_annotation,
    singleR_annotation
  ) %>%
  summarise(number_of_cells = n(), .groups = "drop")

print(cluster_comparison)

write.csv(
  cluster_comparison,
  file.path(results_dir, "10_annotation_comparison_by_cluster.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 12. Save final fully annotated object
# ------------------------------------------------------------

saveRDS(
  seurat_obj,
  file.path(results_dir, "10_final_annotated_seurat_obj.rds")
)

message("SUCCESS: Final automated SingleR annotation completed and saved.")