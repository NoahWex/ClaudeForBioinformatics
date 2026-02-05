#!/usr/bin/env Rscript
#===============================================================================
# 06_integrate.R - Dataset Integration with Harmony
# Workshop: Claude for Bioinformatics, Session 4
#
# Usage:
#   Rscript 06_integrate.R --manifest manifest.yaml
#
# Inputs:
#   - Processed Seurat objects from outputs/processed/
#
# Outputs:
#   - Integrated object: outputs/integrated/integrated.rds
#   - Validation figure: outputs/figures/integration_validation.pdf
#===============================================================================

library(Seurat)
library(harmony)
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
parser <- ArgumentParser(description = "Integrate scRNA-seq datasets with Harmony")
parser$add_argument("--manifest", required = TRUE, help = "Path to manifest.yaml")

args <- parser$parse_args()

log_message("Starting integration pipeline")
log_message("Manifest: ", args$manifest)

#-------------------------------------------------------------------------------
# Load configuration
#-------------------------------------------------------------------------------
log_message("Loading manifest...")
source("../shared/manifest_utils.R")
manifest <- load_manifest(args$manifest)

integration_params <- get_integration_params(manifest)
log_message("Integration parameters:")
log_message("  Method: ", integration_params$method)
log_message("  Batch key: ", integration_params$batch_key)
log_message("  Dimensions: ", integration_params$dims)
log_message("  Resolution: ", integration_params$resolution)

#-------------------------------------------------------------------------------
# Load processed objects
#-------------------------------------------------------------------------------
log_message("Loading processed objects...")

processed_dir <- get_output_dir(manifest, "processed", create = FALSE)
study_names <- names(manifest$data$studies)

seurat_list <- list()
for (study_name in study_names) {
  file_path <- file.path(processed_dir, paste0(study_name, "_processed.rds"))
  if (!file.exists(file_path)) {
    stop("Processed file not found: ", file_path)
  }
  log_message("  Loading ", study_name, "...")
  seurat_list[[study_name]] <- readRDS(file_path)
  log_message("    Cells: ", ncol(seurat_list[[study_name]]))
}

#-------------------------------------------------------------------------------
# Merge objects
#-------------------------------------------------------------------------------
log_message("Merging objects...")

# Merge all objects
merged <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = study_names
)

log_message("Merged object: ", ncol(merged), " cells")

# Verify study_source
log_message("Cells per study:")
print(table(merged$study_source))

gc()

#-------------------------------------------------------------------------------
# Re-normalize merged data
#-------------------------------------------------------------------------------
log_message("Normalizing merged data...")

merged <- NormalizeData(merged, verbose = FALSE)
merged <- FindVariableFeatures(merged, nfeatures = 2000, verbose = FALSE)
merged <- ScaleData(merged, verbose = FALSE)

log_message("Running PCA...")
merged <- RunPCA(merged, npcs = 50, verbose = FALSE)

gc()

#-------------------------------------------------------------------------------
# Run Harmony integration
#-------------------------------------------------------------------------------
log_message("Running Harmony integration...")

merged <- RunHarmony(
  merged,
  group.by.vars = integration_params$batch_key,
  reduction = "pca",
  dims.use = 1:integration_params$dims,
  assay.use = "RNA",
  verbose = TRUE
)

log_message("Harmony complete")
gc()

#-------------------------------------------------------------------------------
# Find neighbors and clusters using Harmony
#-------------------------------------------------------------------------------
log_message("Finding neighbors on Harmony embeddings...")

merged <- FindNeighbors(
  merged,
  reduction = "harmony",
  dims = 1:integration_params$dims,
  verbose = FALSE
)

log_message("Clustering at resolution ", integration_params$resolution, "...")

merged <- FindClusters(
  merged,
  resolution = integration_params$resolution,
  verbose = FALSE
)

log_message("Found ", length(unique(merged$seurat_clusters)), " clusters")

#-------------------------------------------------------------------------------
# Run UMAP
#-------------------------------------------------------------------------------
log_message("Running UMAP on Harmony embeddings...")

merged <- RunUMAP(
  merged,
  reduction = "harmony",
  dims = 1:integration_params$dims,
  verbose = FALSE
)

gc()

#-------------------------------------------------------------------------------
# Validation
#-------------------------------------------------------------------------------
log_message("Generating validation plots...")

# Create figure output directory
fig_dir <- get_output_dir(manifest, "figures", create = TRUE)

# UMAP by study
p1 <- DimPlot(merged, group.by = "study_source", pt.size = 0.1) +
  ggtitle("Study Source") +
  theme_classic()

# UMAP by cell type
p2 <- DimPlot(merged, group.by = "cell_type_harmonized", pt.size = 0.1,
              label = TRUE, label.size = 2.5, repel = TRUE) +
  ggtitle("Cell Type") +
  theme_classic() +
  NoLegend()

# UMAP by cluster
p3 <- DimPlot(merged, group.by = "seurat_clusters", pt.size = 0.1,
              label = TRUE, label.size = 3) +
  ggtitle("Clusters") +
  theme_classic()

# Combine
validation_plot <- (p1 | p2 | p3) +
  plot_annotation(
    title = "Integration Validation",
    subtitle = paste("Method:", integration_params$method,
                     "| Cells:", ncol(merged))
  )

# Save
validation_fig <- file.path(fig_dir, "integration_validation.pdf")
ggsave(validation_fig, validation_plot, width = 15, height = 5)
log_message("Saved validation plot: ", validation_fig)

# Print mixing statistics
log_message("\n=== Integration Quality Check ===")
log_message("Cells per study per cluster:")
print(table(merged$study_source, merged$seurat_clusters))

# Check cell type preservation
log_message("\nCell types per cluster:")
ct_cluster <- table(merged$cell_type_harmonized, merged$seurat_clusters)
# Show dominant cell type per cluster
dominant_ct <- apply(ct_cluster, 2, function(x) names(which.max(x)))
log_message("Dominant cell type per cluster:")
print(dominant_ct)

#-------------------------------------------------------------------------------
# Save integrated object
#-------------------------------------------------------------------------------
output_dir <- get_output_dir(manifest, "integrated", create = TRUE)
output_file <- file.path(output_dir, "integrated.rds")

log_message("Saving integrated object to: ", output_file)
saveRDS(merged, output_file)

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
log_message("\n======================================")
log_message("Integration Complete")
log_message("Method: ", integration_params$method)
log_message("Total cells: ", ncol(merged))
log_message("Studies: ", paste(study_names, collapse = ", "))
log_message("Clusters: ", length(unique(merged$seurat_clusters)))
log_message("Output: ", output_file)
log_message("Validation: ", validation_fig)
log_message("======================================")
