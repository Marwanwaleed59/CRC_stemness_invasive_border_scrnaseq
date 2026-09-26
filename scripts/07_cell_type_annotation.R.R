# =========================================================
# Cell-Type Annotation of Harmony-Integrated scRNA-seq Data
# Dataset: GSE144735 – KUL3 Colorectal Cancer Cohort
# Author: Mohammed Jasim Mohammed Al-Bushhab
# =========================================================

# ---- 1. Load packages ----
library(Seurat)
library(dplyr)
library(ggplot2)

# ---- 2. Set folders ----
project_dir <- "D:/github/CRC_stemness_invasive_border_scrnaseq"
results_dir <- file.path(project_dir, "results")
figures_dir <- file.path(project_dir, "figures")

dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 3. Load final Harmony object ----
seurat_obj <- readRDS("D:/github/CRC_stemness_invasive_border_scrnaseq/results/07_harmony_integrated_unannotated.rds")

DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "seurat_clusters"

# Confirm your final cluster numbers
table(seurat_obj$seurat_clusters)

# =========================================================
# 4. Literature-based canonical markers for CRC scRNA-seq
# =========================================================

canonical_markers <- list(
  
  "T cells" = c("CD3D", "CD3E", "TRAC", "CD247"),
  
  "NK cells" = c("NKG7", "KLRD1", "GNLY", "XCL1", "XCL2"),
  
  "B cells" = c("MS4A1", "CD79A", "CD74", "CD37", "CD19"),
  
  "Plasma cells" = c("JCHAIN", "MZB1", "IGHG1", "IGKC"),
  
  "Myeloid cells" = c("LST1", "TYROBP", "AIF1", "LILRB1"),
  
  "Monocytes" = c("LYZ", "S100A8", "S100A9", "FCN1", "CD14"),
  
  "Macrophages" = c("C1QA", "C1QB", "C1QC", "APOE", "FOLR2"),
  
  "Dendritic cells" = c("FCER1A", "CD1C", "CLEC10A", "HLA-DRA"),
  
  "Mast cells" = c("TPSAB1", "TPSB2", "KIT", "MS4A2"),
  
  "Platelets" = c("PPBP", "PF4", "NRGN"),
  
  "Epithelial cells" = c("EPCAM", "KRT8", "KRT18", "KRT19"),
  
  "Fibroblasts / CAFs" = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"),
  
  "Endothelial cells" = c("PECAM1", "VWF", "KDR", "EMCN"),
  
  "Pericytes" = c("RGS5", "CSPG4", "MCAM", "NOTCH3")
)

# Keep only genes available in your Seurat object
canonical_markers_available <- lapply(
  canonical_markers,
  function(x) x[x %in% rownames(seurat_obj)]
)

canonical_markers_available <- Filter(
  function(x) length(x) > 0,
  canonical_markers_available
)

# =========================================================
# 5. DotPlot — identifies which cluster has each marker set
# =========================================================

p_dotplot <- DotPlot(
  seurat_obj,
  features = canonical_markers_available
) +
  RotatedAxis() +
  ggtitle("Canonical Marker Genes by Harmony Cluster")

print(p_dotplot)

ggsave(
  filename = file.path(figures_dir, "08_canonical_marker_dotplot.png"),
  plot = p_dotplot,
  width = 16,
  height = 10,
  dpi = 300
)

# =========================================================
# 6. FeaturePlot — confirms marker location on the UMAP
# =========================================================

markers_to_check <- c(
  "CD3D", "TRAC",
  "NKG7", "GNLY",
  "MS4A1", "CD79A",
  "JCHAIN", "MZB1",
  "LYZ", "LST1",
  "S100A8", "S100A9",
  "C1QC", "APOE",
  "FCER1A", "CD1C",
  "TPSAB1", "KIT",
  "PPBP", "PF4",
  "EPCAM", "KRT19",
  "COL1A1", "DCN",
  "PECAM1", "VWF",
  "RGS5", "CSPG4"
)

markers_to_check <- markers_to_check[
  markers_to_check %in% rownames(seurat_obj)
]

p_featureplot <- FeaturePlot(
  seurat_obj,
  features = markers_to_check,
  reduction = "umap_after",
  ncol = 4,
  order = TRUE
)

print(p_featureplot)

ggsave(
  filename = file.path(figures_dir, "08_marker_featureplots.png"),
  plot = p_featureplot,
  width = 16,
  height = 22,
  dpi = 300
)

# =========================================================
# 7. Find top marker genes for every cluster
# =========================================================

markers_all <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers_all,
  file.path(results_dir, "08_all_cluster_markers.csv"),
  row.names = FALSE
)

top5_markers <- markers_all %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5, with_ties = FALSE) %>%
  ungroup()

write.csv(
  top5_markers,
  file.path(results_dir, "08_top5_markers_per_cluster.csv"),
  row.names = FALSE
)

print(top5_markers)

marker_summary <- top5_markers %>%
  group_by(cluster) %>%
  summarise(
    top5_genes = paste(gene, collapse = ", "),
    .groups = "drop"
  )

print(marker_summary, n = 18)

# =========================================================
# 8. Heatmap of top cluster markers
# =========================================================

p_heatmap <- DoHeatmap(
  seurat_obj,
  features = unique(top5_markers$gene)
) +
  ggtitle("Top 5 Differentially Expressed Marker Genes per Cluster")

print(p_heatmap)
# =========================================================

ggsave(
  filename = file.path(figures_dir, "08_top5_marker_heatmap.png"),
  plot = p_heatmap,
  width = 15,
  height = 12,
  dpi = 300
)

# =========================================================
# 9. Save object with annotation evidence
# =========================================================

saveRDS(
  seurat_obj,
  file.path(results_dir, "08_annotation_evidence.rds")
)

cat("\n[SUCCESS] Marker checking and annotation evidence completed.\n")
cat("Review the DotPlot and top-5 marker table before assigning final labels.\n")

# =========================================================
# 9. DEG-based cell-type annotation and final UMAP
# =========================================================

seurat_obj$deg_annotation <- dplyr::recode(
  as.character(seurat_obj$seurat_clusters),

  "0"  = "Activated CD4 T cells",
  "1"  = "Epithelial/tumour-like cells",
  "2"  = "Plasma cells",
  "3"  = "FCN1+ Monocytes",
  "4"  = "NK cells",
  "5"  = "Naive B cells",
  "6"  = "PI16+ Fibroblasts",
  "7"  = "Epithelial/tumour-like cells",
  "8"  = "Blood endothelial cells",
  "9"  = "Cycling cells (lineage unresolved)",
  "10" = "COL10A1+ CAFs",
  "11" = "RGS5+ Pericytes",
  "12" = "CCL13+ Macrophage-like cells",
  "13" = "TAC1+ Neural-like cells",
  "14" = "PLP1+ Schwann-like cells",
  "15" = "Mast cells",
  "16" = "Lymphatic endothelial cells",
  "17" = "B cells"
)

Idents(seurat_obj) <- "deg_annotation"

p_deg <- DimPlot(
  seurat_obj,
  reduction = "umap_after",
  group.by = "deg_annotation",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  ggtitle("DEG-Based Annotation of Harmony-Integrated KUL3 Cells")

print(p_deg)

ggsave(
  filename = file.path(figures_dir, "09_deg_based_annotation_umap.png"),
  plot = p_deg,
  width = 14,
  height = 10,
  dpi = 300
)

saveRDS(
  seurat_obj,
  file.path(results_dir, "09_harmony_annotated.rds")
)