
# ============================================
# Stage 1-2: Load, Explore, and Quality Control
# GSE144735 - KUL3 Colorectal Cancer Cohort
# Author: Marwan
# ============================================

# ---- Setup ----
library(Seurat)
library(ggplot2)
library(data.table)
library(Matrix)
library(SingleCellExperiment)
library(scDblFinder)

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

# ---- Stage 1: Load Data ----

counts <- fread(
  "data/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"
)

annotation <- fread(
  "data/GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"
)

# Gene and cell names
gene_names <- counts[[1]]
cell_names <- colnames(counts)[-1]

# Create sparse count matrix
count_sparse <- Matrix(
  as.matrix(counts[, -1, with = FALSE]),
  sparse = TRUE
)

rownames(count_sparse) <- gene_names
colnames(count_sparse) <- cell_names

# Confirm that cells match between counts and annotation
stopifnot(
  identical(colnames(count_sparse), annotation$Index)
)

# ---- Create Seurat Object ----

seurat_obj <- CreateSeuratObject(
  counts = count_sparse,
  project = "GSE144735_KUL3",
  min.cells = 0,
  min.features = 0
)

# Add annotation
annotation_df <- as.data.frame(annotation)
rownames(annotation_df) <- annotation_df$Index

seurat_obj <- AddMetaData(
  seurat_obj,
  metadata = annotation_df
)

# ---- Stage 2: Quality Control ----

# Calculate mitochondrial percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(
  seurat_obj,
  pattern = "^MT-"
)

cells_before_qc <- ncol(seurat_obj)

cat(
  "Cells before QC:",
  cells_before_qc,
  "\n"
)

# ---- QC Visualization ----

qc_violin <- VlnPlot(
  seurat_obj,
  features = c(
    "nFeature_RNA",
    "nCount_RNA",
    "percent.mt"
  ),
  ncol = 3
)

ggsave(
  "figures/01_qc_violin_before_filtering.png",
  plot = qc_violin,
  width = 12,
  height = 5,
  dpi = 300
)

qc_scatter <- FeatureScatter(
  seurat_obj,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)

ggsave(
  "figures/02_qc_count_vs_features.png",
  plot = qc_scatter,
  width = 7,
  height = 6,
  dpi = 300
)

# ---- Low-quality Cell Filtering ----

seurat_obj <- subset(
  seurat_obj,
  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 20
)

cells_after_low_quality_filter <- ncol(seurat_obj)

cat(
  "Cells after low-quality filtering:",
  cells_after_low_quality_filter,
  "\n"
)

# ---- Doublet Detection ----

sce_d1 <- as.SingleCellExperiment(seurat_obj)

set.seed(100)

sce_d1 <- scDblFinder(sce_d1)

# Store doublet score and classification
seurat_obj$doublet_score <-
  colData(sce_d1)$scDblFinder.score

seurat_obj$doublet_class <-
  colData(sce_d1)$scDblFinder.class

# Display doublet counts
print(table(seurat_obj$doublet_class))

# ---- Doublet Visualization ----

doublet_plot <- ggplot(
  seurat_obj@meta.data,
  aes(
    x = nCount_RNA,
    y = doublet_score,
    color = doublet_class
  )
) +
  geom_point(
    size = 1,
    alpha = 0.6
  ) +
  labs(
    x = "nCount_RNA",
    y = "Doublet Score",
    color = "Doublet Classification"
  ) +
  theme_classic()

ggsave(
  "figures/03_doublet_score_vs_nCount_RNA.png",
  plot = doublet_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ---- Remove Predicted Doublets ----

seurat_obj <- subset(
  seurat_obj,
  subset = doublet_class == "singlet"
)

cells_after_doublet_removal <- ncol(seurat_obj)

cat(
  "Cells after doublet removal:",
  cells_after_doublet_removal,
  "\n"
)

# ---- Post-doublet Visualization ----

doublet_after_plot <- VlnPlot(
  seurat_obj,
  features = "doublet_score",
  group.by = "doublet_class"
)

ggsave(
  "figures/04_doublet_score_after_removal.png",
  plot = doublet_after_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# ---- QC Summary ----

number_of_doublets <-
  cells_after_low_quality_filter -
  cells_after_doublet_removal

qc_summary <- data.frame(
  Metric = c(
    "Cells before QC",
    "Cells after low-quality filtering",
    "Predicted doublets removed",
    "Final cells after QC"
  ),
  Number_of_cells = c(
    cells_before_qc,
    cells_after_low_quality_filter,
    number_of_doublets,
    cells_after_doublet_removal
  )
)

write.csv(
  qc_summary,
  "results/01_qc_summary.csv",
  row.names = FALSE
)

# ---- Save QC Metadata ----

metadata_qc <- seurat_obj@meta.data

write.csv(
  metadata_qc,
  "results/01_metadata_after_QC.csv",
  row.names = TRUE
)

# ---- Save QC Checkpoint ----

saveRDS(
  seurat_obj,
  "results/01_after_qc.rds"
)

# ---- Save Software Versions ----

writeLines(
  capture.output(sessionInfo()),
  "results/01_sessionInfo.txt"
)

cat("\nQC completed successfully.\n")
cat(
  "Final number of cells:",
  cells_after_doublet_removal,
  "\n"
)

