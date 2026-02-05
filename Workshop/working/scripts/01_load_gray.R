# 01_load_gray.R
# Load and explore Gray et al. breast tissue scRNA-seq data

library(Seurat)

# ==============================================================================
# Load the data
# ==============================================================================
data_path <- "/share/crsp/lab/dalawson/share/3_Downloaded_Datasets/iHBCA_Reed_2024/component_studies/gray.rds"

message("Loading gray.rds...")
gray <- readRDS(data_path)
message("Done!")

# ==============================================================================
# Basic object information
# ==============================================================================
cat("\n=== Basic Object Information ===\n")
cat("Cells:", ncol(gray), "\n")
cat("Genes:", nrow(gray), "\n")
cat("Assays:", paste(Assays(gray), collapse = ", "), "\n")
cat("Default assay:", DefaultAssay(gray), "\n")

# ==============================================================================
# Metadata exploration
# ==============================================================================
cat("\n=== Metadata Columns ===\n")

meta <- gray@meta.data

# List all columns with types
for (col in colnames(meta)) {
  col_class <- class(meta[[col]])[1]
  n_unique <- length(unique(meta[[col]]))
  cat(sprintf("  %s (%s): %d unique\n", col, col_class, n_unique))
}

# Key categorical columns
cat("\n=== Cell Type Distribution ===\n")
print(table(meta$cell_type))

cat("\n=== Cell Subtype Distribution ===\n")
print(table(meta$cellSubtype))

cat("\n=== Donor Distribution ===\n")
print(table(meta$donor_id))

cat("\n=== Genotype Distribution ===\n")
print(table(meta$genotype_1))
