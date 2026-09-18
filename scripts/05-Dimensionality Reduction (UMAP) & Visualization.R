# ==============================================================================
# Stage 2.4: Non-linear Dimensionality Reduction (UMAP) & Visualization
# Project: GSE144735 - Single-Cell RNA Sequencing Analysis
# Author: Nada Qenawy
# ==============================================================================

# ---- 1. Setup & Environment ----
library(Seurat)
library(ggplot2)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---- 2. Load Clustered Data ----
seurat_obj <- readRDS("results/04_after_clustering.rds")

# ---- 3. Run UMAP ----
# Embed cells in 2D space using UMAP based on top 10 PCs
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# ---- 4. UMAP Visualization ----
umap_plot <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  theme_classic() +
  labs(title = "UMAP Clustering (Res 0.5)")

# Save UMAP Plot
ggsave(
  filename = "figures/09_umap_clusters.png",
  plot = umap_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ---- 5. Save Final Clustered Checkpoint ----
saveRDS(seurat_obj, file = "results/05_seurat_obj_umap_clustered.rds")

cat("\n[SUCCESS] UMAP reduction and plotting completed successfully.\n")