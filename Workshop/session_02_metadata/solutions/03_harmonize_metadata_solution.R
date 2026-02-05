# Solution: 03_harmonize_metadata.R
# Session 2, Exercise 2 - Metadata Harmonization
#
# This is a reference solution. Participant solutions may vary.

library(Seurat)
library(yaml)

# Source manifest utilities
source("../shared/manifest_utils.R")

# ==============================================================================
# Load Configuration
# ==============================================================================

message("Loading manifest...")
manifest <- load_manifest("manifest.yaml")
print_manifest_summary(manifest)

# Load cell type mapping
mapping_file <- file.path(
  dirname("manifest.yaml"),
  manifest$harmonization$mapping_file
)
cell_type_mapping <- yaml::read_yaml(mapping_file)

# ==============================================================================
# Load Data
# ==============================================================================

message("\nLoading datasets...")

study_paths <- get_study_paths(manifest)
seurat_objects <- list()

for (study_name in names(study_paths)) {
  message("  Loading ", study_name, "...")
  seurat_objects[[study_name]] <- readRDS(study_paths[[study_name]])
  message("    Cells: ", ncol(seurat_objects[[study_name]]))
}

# ==============================================================================
# Apply Harmonization
# ==============================================================================

harmonize_study <- function(obj, study_name, manifest, mapping) {
  message("\nHarmonizing ", study_name, "...")

  # Get study-specific config
  study_config <- manifest$harmonization$studies[[study_name]]
  cell_type_col <- study_config$cell_type_column
  batch_col <- study_config$batch_column

  # Add study source
  obj$study_source <- study_name

  # Harmonize cell types (L2 level)
  original_labels <- obj[[cell_type_col, drop = TRUE]]
  l2_map <- unlist(mapping$L2_mapping)

  harmonized_labels <- l2_map[as.character(original_labels)]

  # Report unmapped

  unmapped_mask <- is.na(harmonized_labels)
  if (any(unmapped_mask)) {
    unmapped_types <- unique(original_labels[unmapped_mask])
    warning(
      study_name, ": ", sum(unmapped_mask), " cells with unmapped labels: ",
      paste(unmapped_types, collapse = ", ")
    )
    # Keep original for unmapped
    harmonized_labels[unmapped_mask] <- as.character(original_labels[unmapped_mask])
  }

  # Store harmonized labels
  obj$cell_type_harmonized <- harmonized_labels

  # Rename batch column to unified name
  unified_cols <- manifest$harmonization$unified_columns
  if (!is.null(batch_col) && batch_col %in% colnames(obj@meta.data)) {
    obj[[unified_cols$sample]] <- obj[[batch_col, drop = TRUE]]
  }

  # Store original cell type for reference
  obj$cell_type_original <- original_labels

  message("  Harmonized ", ncol(obj), " cells")
  message("  Unique harmonized labels: ", length(unique(harmonized_labels)))

  return(obj)
}

# Apply to all studies
for (study_name in names(seurat_objects)) {
  seurat_objects[[study_name]] <- harmonize_study(
    seurat_objects[[study_name]],
    study_name,
    manifest,
    cell_type_mapping
  )
}

# ==============================================================================
# Validation
# ==============================================================================

message("\n=== Validation ===")

validate_harmonization <- function(objects) {
  all_valid <- TRUE

  # Check 1: All cells have harmonized labels
  for (name in names(objects)) {
    n_na <- sum(is.na(objects[[name]]$cell_type_harmonized))
    if (n_na > 0) {
      warning(name, ": ", n_na, " cells missing harmonized labels")
      all_valid <- FALSE
    }
  }

  # Check 2: Compare cell types across datasets
  all_types <- lapply(objects, function(x) unique(x$cell_type_harmonized))
  shared_types <- Reduce(intersect, all_types)
  message("Shared cell types across all datasets: ", length(shared_types))

  for (name in names(objects)) {
    unique_to_study <- setdiff(all_types[[name]], shared_types)
    if (length(unique_to_study) > 0) {
      message("  ", name, " unique types: ", paste(unique_to_study, collapse = ", "))
    }
  }

  # Check 3: Cell counts preserved
  message("\nCell counts:")
  for (name in names(objects)) {
    message("  ", name, ": ", ncol(objects[[name]]))
  }

  # Check 4: Study source correctly assigned
  for (name in names(objects)) {
    if (!all(objects[[name]]$study_source == name)) {
      warning(name, ": study_source incorrectly assigned")
      all_valid <- FALSE
    }
  }

  return(all_valid)
}

validation_passed <- validate_harmonization(seurat_objects)

if (validation_passed) {
  message("\nAll validation checks passed!")
} else {
  message("\nSome validation checks failed - review warnings above")
}

# ==============================================================================
# Save Harmonized Objects
# ==============================================================================

output_dir <- get_output_dir(manifest, "processed", create = TRUE)

message("\nSaving harmonized objects...")
for (study_name in names(seurat_objects)) {
  output_path <- file.path(output_dir, paste0(study_name, "_harmonized.rds"))
  saveRDS(seurat_objects[[study_name]], output_path)
  message("  Saved: ", output_path)
}

# ==============================================================================
# Summary
# ==============================================================================

message("\n=== Harmonization Summary ===")
message("Studies processed: ", paste(names(seurat_objects), collapse = ", "))
message("Total cells: ", sum(sapply(seurat_objects, ncol)))
message("Output directory: ", output_dir)
message("\nHarmonization complete!")
