# ==============================================================================
# Stage 2.2: Principal Component Analysis (PCA) & Dimensionality Assessment
# Project: GSE144735 - Single-Cell RNA Sequencing Analysis
# Author: Nada Qenawy
# ==============================================================================

# ---- 1. Setup & Environment ----
library(Seurat)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---- 2. Load Normalized & Scaled Data ----
seurat_obj <- readRDS("results/02_normalized_scaled.rds")

# ---- 3. Run Principal Component Analysis (PCA) ----
# Run PCA on the top highly variable features
seurat_obj <- RunPCA(
  seurat_obj,
  features = VariableFeatures(object = seurat_obj)
)

# ---- 4. PCA Visualizations ----
# A) PCA Scatter Plot
pca_dimplot <- DimPlot(seurat_obj, reduction = "pca")
ggsave(
  filename = "figures/06_pca_dimplot.png",
  plot = pca_dimplot,
  width = 7,
  height = 6,
  dpi = 300
)

# B) Elbow Plot to determine significant PCs
pca_elbow <- ElbowPlot(seurat_obj)
ggsave(
  filename = "figures/07_pca_elbow_plot.png",
  plot = pca_elbow,
  width = 7,
  height = 5,
  dpi = 300
)

# C) Dimension Heatmap for top 10 PCs
png("figures/08_pca_dim_heatmap.png", width = 10, height = 8, units = "in", res = 300)
DimHeatmap(
  seurat_obj,
  dims = 1:10,
  cells = 500,
  balanced = TRUE
)
dev.off()

# ---- 5. Save Checkpoint ----
saveRDS(seurat_obj, file = "results/03_after_pca.rds")

cat("\n[SUCCESS] PCA analysis and visualization completed successfully.\n")