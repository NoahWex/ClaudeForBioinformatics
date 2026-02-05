#!/usr/bin/env Rscript
#===============================================================================
# 05_qc_normalize.R - QC and Normalization Pipeline
# Workshop: Claude for Bioinformatics, Session 3
#
# Usage:
#   Rscript 05_qc_normalize.R --study gray --manifest manifest.yaml
#   Rscript 05_qc_normalize.R --study gray --manifest manifest.yaml --test
#
# Inputs:
#   - Harmonized Seurat object: outputs/processed/{study}_harmonized.rds
#   - Manifest with QC parameters
#
# Outputs:
#   - Processed object: outputs/processed/{study}_processed.rds
#   - QC report: outputs/figures/{study}_qc_report.pdf
#
# Required packages: Seurat, argparse, ggplot2, patchwork
#===============================================================================

library(Seurat)
library(argparse)
library(ggplot2)
library(patchwork)

#-------------------------------------------------------------------------------
# Logging helper
#-------------------------------------------------------------------------------
log_message <- function(...) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(timestamp, " ", ...)
}

#-------------------------------------------------------------------------------
# Argument parsing
#-------------------------------------------------------------------------------
parser <- ArgumentParser(description = "QC and normalization for scRNA-seq data")
parser$add_argument("--study", required = TRUE, help = "Study name (e.g., 'gray')")
parser$add_argument("--manifest", required = TRUE, help = "Path to manifest.yaml")
parser$add_argument("--test", action = "store_true",
                    help = "Test mode: subset to 1000 cells")

args <- parser$parse_args()

log_message("Starting QC pipeline")
log_message("Study: ", args$study)
log_message("Manifest: ", args$manifest)
log_message("Test mode: ", args$test)

#-------------------------------------------------------------------------------
# Load configuration
#-------------------------------------------------------------------------------
log_message("Loading manifest...")
source("../shared/manifest_utils.R")
manifest <- load_manifest(args$manifest)

qc_params <- get_qc_params(manifest)
log_message("QC parameters loaded:")
log_message("  min_genes: ", qc_params$min_genes)
log_message("  max_genes: ", qc_params$max_genes)
log_message("  max_mito_pct: ", qc_params$max_mito_pct)

#-------------------------------------------------------------------------------
# Load data
#-------------------------------------------------------------------------------
input_dir <- get_output_dir(manifest, "processed", create = FALSE)
input_file <- file.path(input_dir, paste0(args$study, "_harmonized.rds"))

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

log_message("Loading data from: ", input_file)
seurat_obj <- readRDS(input_file)
log_message("Loaded ", ncol(seurat_obj), " cells, ", nrow(seurat_obj), " genes")

# Test mode: subset
if (args$test) {
  log_message("TEST MODE: Subsetting to 1000 cells")
  set.seed(42)
  cells_keep <- sample(colnames(seurat_obj), min(1000, ncol(seurat_obj)))
  seurat_obj <- subset(seurat_obj, cells = cells_keep)
  log_message("Subsetted to ", ncol(seurat_obj), " cells")
}

#-------------------------------------------------------------------------------
# Calculate QC metrics
#-------------------------------------------------------------------------------
log_message("Calculating QC metrics...")

# Mitochondrial genes
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# Ribosomal genes
seurat_obj[["percent.ribo"]] <- PercentageFeatureSet(seurat_obj, pattern = "^RP[SL]")

# Store pre-filter counts
n_cells_pre <- ncol(seurat_obj)

# QC summary
log_message("QC metrics summary:")
log_message("  nFeature_RNA: ", round(mean(seurat_obj$nFeature_RNA)), " (mean)")
log_message("  nCount_RNA: ", round(mean(seurat_obj$nCount_RNA)), " (mean)")
log_message("  percent.mt: ", round(mean(seurat_obj$percent.mt), 2), "% (mean)")

#-------------------------------------------------------------------------------
# QC Filtering
#-------------------------------------------------------------------------------
log_message("Applying QC filters...")

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA >= qc_params$min_genes &
           nFeature_RNA <= qc_params$max_genes &
           percent.mt <= qc_params$max_mito_pct
)

n_cells_post <- ncol(seurat_obj)
n_filtered <- n_cells_pre - n_cells_post

log_message("Filtered ", n_filtered, " cells (", round(100 * n_filtered / n_cells_pre, 1), "%)")
log_message("Remaining: ", n_cells_post, " cells")

if (n_cells_post == 0) {
  stop("No cells passed QC filters! Check your thresholds.")
}

gc()  # Free memory

#-------------------------------------------------------------------------------
# Normalization
#-------------------------------------------------------------------------------
log_message("Normalizing data...")

norm_method <- qc_params$normalization_method

if (tolower(norm_method) == "sctransform") {
  log_message("Using SCTransform (regressing percent.mt)...")
  seurat_obj <- SCTransform(
    seurat_obj,
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )
} else {
  log_message("Using LogNormalize...")
  seurat_obj <- NormalizeData(
    seurat_obj,
    normalization.method = "LogNormalize",
    scale.factor = qc_params$scale_factor,
    verbose = FALSE
  )
}

gc()

#-------------------------------------------------------------------------------
# Variable features
#-------------------------------------------------------------------------------
log_message("Finding variable features...")

seurat_obj <- FindVariableFeatures(
  seurat_obj,
  selection.method = qc_params$selection_method,
  nfeatures = qc_params$n_variable_features,
  verbose = FALSE
)

log_message("Selected ", length(VariableFeatures(seurat_obj)), " variable features")

#-------------------------------------------------------------------------------
# Scaling (if not SCTransform)
#-------------------------------------------------------------------------------
if (tolower(norm_method) != "sctransform") {
  log_message("Scaling data...")
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
}

gc()

#-------------------------------------------------------------------------------
# PCA
#-------------------------------------------------------------------------------
log_message("Running PCA...")
seurat_obj <- RunPCA(seurat_obj, npcs = 50, verbose = FALSE)
log_message("Computed 50 PCs")

gc()

#-------------------------------------------------------------------------------
# Generate QC Report Figure
#-------------------------------------------------------------------------------
log_message("Generating QC report figure...")

# Create output directory
fig_dir <- get_output_dir(manifest, "figures", create = TRUE)
if (args$test) {
  fig_dir <- file.path(fig_dir, "test")
  dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
}

# Pre vs post filtering (store in object for plotting)
filter_summary <- data.frame(
  Stage = c("Pre-filter", "Post-filter"),
  Cells = c(n_cells_pre, n_cells_post)
)

p1 <- ggplot(filter_summary, aes(x = Stage, y = Cells, fill = Stage)) +
  geom_col() +
  geom_text(aes(label = Cells), vjust = -0.5) +
  scale_fill_manual(values = c("Pre-filter" = "gray70", "Post-filter" = "steelblue")) +
  labs(title = paste(args$study, "- Cell Filtering"), y = "Number of Cells") +
  theme_minimal() +
  theme(legend.position = "none")

# QC distributions
p2 <- VlnPlot(seurat_obj, features = "nFeature_RNA", pt.size = 0) +
  labs(title = "Genes per Cell") +
  theme(legend.position = "none")

p3 <- VlnPlot(seurat_obj, features = "nCount_RNA", pt.size = 0) +
  labs(title = "UMIs per Cell") +
  theme(legend.position = "none")

p4 <- VlnPlot(seurat_obj, features = "percent.mt", pt.size = 0) +
  labs(title = "% Mitochondrial") +
  theme(legend.position = "none")

# Combine
report_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = paste("QC Report:", args$study),
    subtitle = paste("Processed:", Sys.Date())
  )

# Save
output_fig <- file.path(fig_dir, paste0(args$study, "_qc_report.pdf"))
ggsave(output_fig, report_plot, width = 10, height = 8)
log_message("Saved QC report: ", output_fig)

#-------------------------------------------------------------------------------
# Save processed object
#-------------------------------------------------------------------------------
output_dir <- get_output_dir(manifest, "processed", create = TRUE)
if (args$test) {
  output_dir <- file.path(output_dir, "test")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
}

output_file <- file.path(output_dir, paste0(args$study, "_processed.rds"))
log_message("Saving processed object to: ", output_file)
saveRDS(seurat_obj, output_file)

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
log_message("======================================")
log_message("QC Pipeline Complete")
log_message("Study: ", args$study)
log_message("Input cells: ", n_cells_pre)
log_message("Output cells: ", n_cells_post)
log_message("Cells filtered: ", n_filtered, " (", round(100 * n_filtered / n_cells_pre, 1), "%)")
log_message("Output: ", output_file)
log_message("QC Report: ", output_fig)
log_message("======================================")
