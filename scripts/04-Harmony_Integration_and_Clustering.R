# ==============================================================================
# Stage 2.1: Harmony_Integration_and_Clustering
# Project: GSE144735 - Single-Cell RNA Sequencing Analysis
# Author: Nada Qenawy
# ==============================================================================

# 1. Load Required Libraries
if (!requireNamespace("harmony", quietly = TRUE)) {
  install.packages("harmony")
}

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)

# 2. Load Pre-processed Seurat Object
seurat_obj <- readRDS("D:/github/CRC_stemness_invasive_border_scrnaseq/results/06_baseline_before_harmony.rds")

# 3. Visualize Pre-Integration UMAP Grouped by Sample
p1 <- DimPlot(
  seurat_obj,
  reduction = "umap",
  group.by = "orig.ident"
) +
  ggtitle("Before Integration")

print(p1)

# 4. Perform Harmony Integration Across Samples (orig.ident)
seurat_obj <- seurat_obj %>%
  RunHarmony(
    group.by.vars = "orig.ident",  
    plot_convergence = TRUE,
    nclust = 50,
    max_iter = 10,
    early_stop = TRUE
  )

set.seed(1234)

seurat_obj <- RunHarmony(
  seurat_obj,
  group.by.vars = "patient_id",
  reduction = "pca"
)

set.seed(1234)

seurat_obj <- RunUMAP(
  seurat_obj,
  reduction = "harmony",
  dims = 1:15,
  reduction.name = "umap_harmony"
)



# 5. Run UMAP on Harmony Embeddings
seurat_obj <- RunUMAP(seurat_obj, reduction = "harmony", dims = 1:15, reduction.name = "umap_after")

# 6. Perform Clustering on Integrated Data
DimPlot(
  seurat_obj,
  reduction = "umap_harmony",
  group.by = "orig.ident"
) + ggtitle("After Harmony Integration")

# 7. Visualize Post-Integration UMAP Grouped by Sample
p2 <- DimPlot(seurat_obj, reduction = "umap_after", group.by = "orig.ident") + 
  ggtitle("After Integration")
print(p2)

# 8. Compare Before and After Integration
p1 + p2

# 9. Visualize Final Integrated Clusters
p_clusters <- DimPlot(
  seurat_obj, 
  reduction = "umap_after", 
  group.by = "seurat_clusters", 
  label = TRUE, 
  repel = TRUE
) + ggtitle("Harmony Integrated Clusters")

print(p_clusters)
##############
saveRDS(
  seurat_obj,
  "D:/github/CRC_stemness_invasive_border_scrnaseq/results/07_harmony_integrated_unannotated.rds"
)

file.exists(
  "D:/github/CRC_stemness_invasive_border_scrnaseq/results/07_harmony_integrated_unannotated.rds"
)
