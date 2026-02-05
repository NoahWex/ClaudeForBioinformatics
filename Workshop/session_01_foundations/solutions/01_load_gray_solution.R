# Solution: 01_load_gray.R
# Session 1, Exercise 2 - Data Loading and Exploration
#
# This is a reference solution. Participant solutions may vary.

library(Seurat)
library(ggplot2)

# ==============================================================================
# Part 1: Load the data
# ==============================================================================

message("Loading gray.rds...")
gray <- readRDS("../data/raw/gray.rds")

# Basic object info
cat("\n=== Basic Object Information ===\n")
cat("Cells:", ncol(gray), "\n")
cat("Genes:", nrow(gray), "\n")
cat("Assays:", paste(Assays(gray), collapse = ", "), "\n")
cat("Default assay:", DefaultAssay(gray), "\n")

# ==============================================================================
# Part 2: Explore metadata
# ==============================================================================

cat("\n=== Metadata Exploration ===\n")

# Get metadata
meta <- gray@meta.data

# List columns with types
cat("\nMetadata columns:\n")
for (col in colnames(meta)) {
  col_class <- class(meta[[col]])[1]
  n_unique <- length(unique(meta[[col]]))
  cat(sprintf("  - %s (%s): %d unique values\n", col, col_class, n_unique))
}

# For categorical columns, show value counts
cat("\n=== Categorical Column Details ===\n")

categorical_cols <- sapply(meta, function(x) {
  is.factor(x) || is.character(x) || (is.numeric(x) && length(unique(x)) < 20)
})

for (col in names(categorical_cols)[categorical_cols]) {
  cat(sprintf("\n%s:\n", col))
  tbl <- sort(table(meta[[col]]), decreasing = TRUE)
  if (length(tbl) <= 10) {
    print(tbl)
  } else {
    cat(sprintf("  (%d unique values, showing top 10)\n", length(tbl)))
    print(head(tbl, 10))
  }
}

# Identify likely cell type column
# Look for columns with names containing "cell", "type", "cluster", "annotation"
cell_type_patterns <- c("cell", "type", "cluster", "annot", "label")
likely_celltype_cols <- colnames(meta)[
  grepl(paste(cell_type_patterns, collapse = "|"), colnames(meta), ignore.case = TRUE)
]
cat("\nLikely cell type columns:", paste(likely_celltype_cols, collapse = ", "), "\n")

# ==============================================================================
# Part 3: QC metrics
# ==============================================================================

cat("\n=== QC Metrics ===\n")

# Calculate percent mitochondrial
gray[["percent.mt"]] <- PercentageFeatureSet(gray, pattern = "^MT-")

# Summary statistics
cat("\nnFeature_RNA summary:\n")
print(summary(gray$nFeature_RNA))

cat("\nnCount_RNA summary:\n")
print(summary(gray$nCount_RNA))

cat("\npercent.mt summary:\n")
print(summary(gray$percent.mt))

# Create QC violin plot
qc_plot <- VlnPlot(
  gray,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0
)

# Save plot
output_dir <- "outputs"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

ggsave(
  file.path(output_dir, "qc_metrics_gray.pdf"),
  qc_plot,
  width = 12,
  height = 4
)
cat("\nSaved QC plot to:", file.path(output_dir, "qc_metrics_gray.pdf"), "\n")

# ==============================================================================
# Summary for manifest
# ==============================================================================

cat("\n=== Summary for Manifest ===\n")
cat("File: gray.rds\n")
cat("n_cells:", ncol(gray), "\n")
cat("cell_type_column: [CHECK OUTPUT ABOVE]\n")
cat("batch_column: [CHECK OUTPUT ABOVE]\n")

message("\nExploration complete!")
