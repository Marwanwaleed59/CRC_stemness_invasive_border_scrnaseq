# ==============================================================================
# Stage 2.3: Graph-Based Cell Clustering (FindNeighbors & FindClusters)
# Project: GSE144735 - Single-Cell RNA Sequencing Analysis
# Author: Nada Qenawy
# ==============================================================================

# ---- 1. Setup & Environment ----
library(Seurat)

dir.create("results", showWarnings = FALSE)

# ---- 2. Load PCA Checkpoint ----
seurat_obj <- readRDS("results/03_after_pca.rds")

# ---- 3. Find Nearest Neighbors ----
# Construct SNN graph using top 10 principal components
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)

# ---- 4. Find Clusters ----
# Identify cell clusters at resolution 0.5
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# ---- 5. Save Checkpoint ----
saveRDS(seurat_obj, file = "results/04_after_clustering.rds")

cat("\n[SUCCESS] Cell clustering completed successfully.\n")