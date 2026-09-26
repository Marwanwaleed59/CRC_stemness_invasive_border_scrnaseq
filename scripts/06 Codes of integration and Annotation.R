(!requireNamespace("harmony", quietly = TRUE)) 
install.packages("harmony")
library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

seurat_obj <- readRDS("D:/github/CRC_stemness_invasive_border_scrnaseq/results/06_baseline_before_harmony.rds")

library(Seurat)
library(ggplot2)
***********************************************************************
p1 <- DimPlot(seurat_obj, reduction = "pca", group.by = "orig.ident") +
  ggtitle("Before Integration")

print(p1)
***********************************************************************
seurat_obj$patient_id <- sub("-.*", "", seurat_obj$orig.ident)
library(harmony)
reduction = "pca"
***********************************************************************
seurat_obj <- RunHarmony(
    seurat_obj,
    group.by.vars = "patient_id",
    reduction.use = "pca",
    plot_convergence = TRUE,
    nclust = 50,
    max_iter = 10,
    early_stop = TRUE
  )
*************************************************************************
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:15, reduction.name = "umap_after")
*************************************************************************
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:15)
**************************************************************************
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
**************************************************************************
p2 <- DimPlot(seurat_obj, reduction = "umap_after", group.by = "orig.ident") + ggtitle("After Integration")
print(p2)
***************************************************************************
p1 + p2
***************************************************************************
seurat_obj <- FindNeighbors(seurat_obj, reduction = "harmony", dims = 1:15)
***************************************************************************
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
***************************************************************************
p_clusters <- DimPlot(
  seurat_obj, 
  reduction = "umap_after", 
  group.by = "seurat_clusters", 
  label = TRUE, 
  repel = TRUE
) + ggtitle("Harmony Integrated Clusters")
print(p_clusters)
***********************************************************************
table(Idents(seurat_obj))
***********************************************************************
canonical_markers <- list(
    "T cells" = c("CD3D", "CD3E", "TRAC"),
    "NK cells" = c("NKG7", "GNLY", "KLRD1"),
    "B cells" = c("MS4A1", "CD79A", "CD74"),
    "Plasma cells" = c("MZB1", "JCHAIN", "IGHG1"),
    "Myeloid cells" = c("LYZ", "LST1", "TYROBP"),
    "Macrophages" = c("C1QC", "C1QA", "APOE"),
    "Dendritic cells" = c("FCER1A", "CD1C", "CLEC10A"),
    "Platelets" = c("PPBP", "PF4"),
    "Epithelial/tumour" = c("EPCAM", "KRT8", "KRT18", "KRT19"),
    "Fibroblasts/CAFs" = c("COL1A1", "COL1A2", "DCN", "LUM"),
    "Endothelial cells" = c("PECAM1", "VWF", "KDR")
  )
*****************************************************************************************
DotPlot(
    seurat_obj,
    features = canonical_markers
  ) +
  RotatedAxis() +
  ggtitle("Canonical Markers by Cluster")
************************************************
  DefaultAssay(seurat_obj) <- "RNA"
Idents(seurat_obj) <- "seurat_clusters"
*************************************
markers_all <- FindAllMarkers(
  seurat_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)
****************************************************************************
FeaturePlot(
    seurat_obj,
    features = unlist(canonical_markers),
    reduction = "umap_after",
    ncol = 4
  )
**********************************************************************************
VlnPlot(seurat_obj, features = "MS4A1")
*********************************************************************************
# Check canonical markers on the Harmony UMAP

markers_to_check <- c(
  "CD3D",    # T cells
  "MS4A1",   # B cells
  "NKG7",    # NK cells
  "LYZ",     # Myeloid cells
  "CD14",    # Classical monocytes
  "FCGR3A",  # Non-classical monocytes
  "CD1C",    # Dendritic cells
  "PPBP"     # Platelets
)

FeaturePlot(
  seurat_obj,
  features = markers_to_check,
  reduction = "umap_after",
  ncol = 4,
  order = TRUE
)  
*********************************************************************************
seurat_obj$canonical_annotation <- dplyr::recode(
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
    "19" = "B cells / Plasma cells"
  )
*********************************************************************
seurat_obj$canonical_annotation <- dplyr::recode(
    as.character(seurat_obj$seurat_clusters),
    
    "0"  = "CD3D+ TRAC+ T cells",
    "1"  = "LYZ+ CD14+ Monocytes",
    "2"  = "CD3D+ TRAC+ T cells",
    "3"  = "LYZ+ LST1+ Myeloid cells",
    "4"  = "LYZ+ CD14+ Monocytes",
    "5"  = "MS4A1+ CD79A+ B cells",
    "6"  = "CD3D+ TRAC+ T cells",
    "7"  = "NKG7+ GNLY+ NK cells",
    "8"  = "LYZ+ CD14+ Monocytes",
    "9"  = "CD1C+ FCER1A+ Dendritic cells",
    "10" = "CD3D+ TRAC+ T cells",
    "11" = "LYZ+ CD14+ Monocytes",
    "12" = "LYZ+ LST1+ Myeloid cells",
    "13" = "CD3D+ TRAC+ T cells",
    "14" = "LYZ+ CD14+ Monocytes",
    "15" = "CD3D+ TRAC+ T cells",
    "16" = "NKG7+ GNLY+ NK cells",
    "17" = "LYZ+ CD14+ Monocytes",
    "18" = "PPBP+ Platelets",
    "19" = "B / plasma-like cells (verify)"
  )

Idents(seurat_obj) <- "canonical_annotation" 

DimPlot(
  seurat_obj,
  reduction = "umap_after",
  group.by = "canonical_annotation",
  label = TRUE,
  repel = TRUE,
  pt.size = 0.3
) +
  NoLegend()
*********************************************************************
  p_canonical <- DimPlot(
  seurat_obj, 
  reduction = "umap_after", 
  group.by = "canonical_annotation", 
  label = TRUE, 
  repel = TRUE
) + ggtitle("1. Canonical Marker-Based Annotation")
print(p_canonical)
**********************************************************************
markers_all <- FindAllMarkers(
  seurat_obj, 
  only.pos = TRUE, 
  min.pct = 0.25, 
  logfc.threshold = 0.25
)
************************************************************************
top5_markers <- markers_all %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC)
************************************************************************
DoHeatmap(seurat_obj, features = top5_markers$gene) + 
  ggtitle("DEG Top Markers Heatmap")
**************************************************************************
p_deg <- DimPlot(
    seurat_obj,
    reduction = "umap",
    group.by = "deg_annotation",
    label = TRUE,
    repel = TRUE
  ) + ggtitle("2. DEG-Based Annotation")

print(p_deg)
***************************************************************************
  # Install SingleR only if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }

BiocManager::install(c("SingleR", "celldex"), ask = FALSE, update = FALSE)
***************************************************************************
library(SingleR)
library(celldex)

saveRDS(
  seurat_obj,
  "D:/github/CRC_stemness_invasive_border_scrnaseq/results/06_baseline_before_harmony.rds"
)



***************************************************************************
hpca_ref <- celldex::HumanPrimaryCellAtlasData()
***************************************************
sce_data <- as.SingleCellExperiment(seurat_obj)
***************************************************
singler_res <- SingleR(
  test = sce_data, 
  ref = hpca_ref, 
  labels = hpca_ref$label.main
)
*****************************************************************
seurat_obj$SingleR_hpca_annotation <- singler_res$labels
p_singler <- DimPlot(
  seurat_obj, 
  reduction = "umap_after", 
  group.by = "SingleR_hpca_annotation", 
  label = TRUE, 
  repel = TRUE
) + ggtitle("3. SingleR Automated Annotation (HPCA)")

print(p_singler)