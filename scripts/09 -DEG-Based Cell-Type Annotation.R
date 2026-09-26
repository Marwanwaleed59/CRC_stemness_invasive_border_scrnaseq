# ============================================================
# Stage 09: DEG-Based Cell-Type Annotation
# Dataset: GSE144735 (KUL3 colorectal cancer scRNA-seq)
# Input: completed canonical-marker annotation object
# ============================================================

# ---- 1. Packages ----
library(Seurat)
library(dplyr)
library(ggplot2)

# ---- 2. Project folders ----
project_dir <- "D:/github/CRC_stemness_invasive_border_scrnaseq"

results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 3. Load your completed canonical annotation object ----
seurat_obj <- readRDS("D:/github/CRC_stemness_invasive_border_scrnaseq/results/09_harmony_annotated.rds")


DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "seurat_clusters"

# Check the clusters present in your object
print(table(seurat_obj$seurat_clusters))

# ============================================================
# 4. Find Differentially Expressed Genes (DEGs) per cluster
# ============================================================

markers_file <- file.path(results_dir, "08_all_cluster_markers.csv")

if (file.exists(markers_file)) {
  
  message("Existing DEG marker table found. Loading it...")
  markers_all <- read.csv(markers_file)
  
} else {
  
  message("Finding DEGs for all clusters. This may take several minutes...")
  
  markers_all <- FindAllMarkers(
    object = seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  write.csv(
    markers_all,
    file.path(results_dir, "09_all_cluster_DEGs.csv"),
    row.names = FALSE
  )
}

# Save a clean copy for this DEG stage
write.csv(
  markers_all,
  file.path(results_dir, "09_all_cluster_DEGs.csv"),
  row.names = FALSE
)

# ---- 5. Select top 10 DEGs from every cluster ----
# Works with both Seurat versions: avg_log2FC or avg_logFC

fc_column <- if ("avg_log2FC" %in% colnames(markers_all)) {
  "avg_log2FC"
} else {
  "avg_logFC"
}

top10_degs <- markers_all %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(
    order_by = .data[[fc_column]],
    n = 10,
    with_ties = FALSE
  ) %>%
  ungroup()

write.csv(
  top10_degs,
  file.path(results_dir, "09_top10_DEGs_per_cluster.csv"),
  row.names = FALSE
)

print(top10_degs)

# ============================================================
# 6. DEG heatmap
# ============================================================

top_deg_genes <- unique(top10_degs$gene)
top_deg_genes <- top_deg_genes[top_deg_genes %in% rownames(seurat_obj)]

deg_heatmap <- DoHeatmap(
  object = seurat_obj,
  features = top_deg_genes,
  group.by = "seurat_clusters",
  size = 3
) +
  ggtitle("Top Differentially Expressed Genes per Cluster")

print(deg_heatmap)

ggsave(
  filename = file.path(figures_dir, "09_top10_DEG_heatmap.png"),
  plot = deg_heatmap,
  width = 14,
  height = 12,
  dpi = 300
)

# ============================================================
# 7. DEG-based cell-type annotation
# Labels are supported by the cluster DEG results.
# ============================================================

seurat_obj$deg_annotation <- dplyr::recode(
  as.character(seurat_obj$seurat_clusters),
  
  "0"  = "T cells",
  "1"  = "Monocytes",
  "2"  = "T cells",
  "3"  = "Myeloid cells",
  "4"  = "Monocytes",
  "5"  = "B cells",
  "6"  = "T cells",
  "7"  = "NK cells",
  "8"  = "Monocytes",
  "9"  = "Dendritic cells",
  "10" = "T cells",
  "11" = "Monocytes",
  "12" = "Myeloid cells",
  "13" = "T cells",
  "14" = "Monocytes",
  "15" = "T cells",
  "16" = "NK cells",
  "17" = "Monocytes",
  "18" = "Platelets",
  "19" = "B / Plasma-like cells",
  
  .default = "Unassigned - review DEGs"
)

# Check number of cells in each DEG-based label
print(table(seurat_obj$deg_annotation))

# ============================================================
# 8. DEG-based annotation UMAP
# ============================================================

deg_umap <- DimPlot(
  object = seurat_obj,
  reduction = "umap_after",
  group.by = "deg_annotation",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.25
) +
  ggtitle("DEG-Based Annotation")

print(deg_umap)

ggsave(
  filename = file.path(figures_dir, "09_DEG_based_annotation_UMAP.png"),
  plot = deg_umap,
  width = 13,
  height = 8,
  dpi = 300
)

# ============================================================
# 9. Save final DEG-annotated Seurat object
# ============================================================

saveRDS(
  seurat_obj,
  file.path(results_dir, "09_DEG_annotated_seurat_obj.rds")
)

message("SUCCESS: DEG-based annotation completed and saved.")