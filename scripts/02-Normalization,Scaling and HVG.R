# ==============================================================================
# Stage 2.1: Data Normalization, Scaling & Feature Selection (HVGs)
# Project: GSE144735 - Single-Cell RNA Sequencing Analysis
# Author: Nada Qenawy
# ==============================================================================

# ---- 1. Setup & Environment ----
library(Seurat)
library(ggplot2)

# Ensure output directories exist
dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---- 2. Load Checkpoint Data ----
# Load object generated from QC stage
seurat_obj <- readRDS("results/01_after_qc.rds")

# ---- 3. Log Normalization ----
# Normalize gene expression counts per cell using standard log-normalization
seurat_obj <- NormalizeData(
  seurat_obj,
  normalization.method = "LogNormalize",
  scale.factor = 10000
)

# ---- 4. Feature Selection (Highly Variable Genes - HVGs) ----
# Identify top 2000 highly variable features for downstream analysis
seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = "vst",
  nfeatures = 2000
)

# Extract top 10 most variable genes for labeling
top10_genes <- head(VariableFeatures(seurat_obj), 10)

# ---- 5. HVG Visualization ----
var_plot <- VariableFeaturePlot(seurat_obj)
var_plot_labeled <- LabelPoints(
  plot = var_plot,
  points = top10_genes,
  repel = TRUE
)

# Save HVG plot to figures directory
ggsave(
  filename = "figures/05_highly_variable_genes.png",
  plot = var_plot_labeled,
  width = 8,
  height = 6,
  dpi = 300
)

# ---- 6. Data Scaling ----
# Scale data so mean expression is 0 and variance is 1 across cells
seurat_obj <- ScaleData(seurat_obj)

# ---- 7. Save Checkpoint ----
saveRDS(seurat_obj, file = "results/02_normalized_scaled.rds")

cat("\n[SUCCESS] Normalization, HVG identification, and Scaling completed successfully.\n")