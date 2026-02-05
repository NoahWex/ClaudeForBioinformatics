# manifest_utils.R - Utility functions for manifest-driven workflows
#
# This module provides functions for reading manifest files and
# constructing paths based on manifest configuration.
#
# Usage:
#   source("shared/manifest_utils.R")
#   manifest <- load_manifest("manifest.yaml")
#   paths <- get_study_paths(manifest)

library(yaml)

#' Load and validate a manifest file
#'
#' @param manifest_path Path to the manifest YAML file
#' @return Parsed manifest as a list
#' @export
load_manifest <- function(manifest_path) {
  if (!file.exists(manifest_path)) {
    stop("Manifest file not found: ", manifest_path)
  }

  manifest <- yaml::read_yaml(manifest_path)

  # Validate required fields
  required <- c("project_name", "data")
  missing <- setdiff(required, names(manifest))
  if (length(missing) > 0) {
    stop("Manifest missing required fields: ", paste(missing, collapse = ", "))
  }

  # Store manifest directory for relative path resolution
  manifest$manifest_dir <- dirname(normalizePath(manifest_path))

  message("Loaded manifest: ", manifest$project_name)
  return(manifest)
}


#' Get absolute paths for all studies in manifest
#'
#' Automatically detects HPC vs local environment and uses appropriate base path.
#'
#' @param manifest Loaded manifest object
#' @param force_env Optional: force "hpc" or "local" environment
#' @return Named list of study paths
#' @export
get_study_paths <- function(manifest, force_env = NULL) {
  base_path <- resolve_data_path(manifest, force_env)

  paths <- lapply(manifest$data$studies, function(study) {
    file.path(base_path, study$file)
  })

  names(paths) <- names(manifest$data$studies)
  return(paths)
}

#' Resolve the appropriate data base path for current environment
#'
#' @param manifest Loaded manifest object
#' @param force_env Optional: force "hpc" or "local" environment
#' @return Base path string
#' @export
resolve_data_path <- function(manifest, force_env = NULL) {
  # Check for dual-path manifest structure
  if (!is.null(manifest$data$base_path_hpc) && !is.null(manifest$data$base_path_local)) {
    env <- force_env
    if (is.null(env)) {
      # Auto-detect: HPC paths start with /share or /dfs
      env <- if (file.exists("/share/crsp")) "hpc" else "local"
    }

    base_path <- if (env == "hpc") {
      manifest$data$base_path_hpc
    } else {
      manifest$data$base_path_local
    }
    message("Using ", env, " data path: ", base_path)
  } else if (!is.null(manifest$data$base_path)) {
    # Legacy single-path manifest
    base_path <- file.path(manifest$manifest_dir, manifest$data$base_path)
  } else {
    stop("No data base_path configured in manifest")
  }

  return(base_path)
}


#' Get path for a specific study
#'
#' @param manifest Loaded manifest object
#' @param study_name Name of the study (e.g., "gray", "pal")
#' @param force_env Optional: force "hpc" or "local" environment
#' @return Absolute path to the study file
#' @export
get_study_path <- function(manifest, study_name, force_env = NULL) {
  if (!study_name %in% names(manifest$data$studies)) {
    stop("Study not found in manifest: ", study_name,
         "\nAvailable studies: ", paste(names(manifest$data$studies), collapse = ", "))
  }

  base_path <- resolve_data_path(manifest, force_env)
  return(file.path(base_path, manifest$data$studies[[study_name]]$file))
}


#' Get output directory path
#'
#' @param manifest Loaded manifest object
#' @param output_type Type of output: "processed", "integrated", "reports", "figures"
#' @param create If TRUE, create directory if it doesn't exist
#' @return Absolute path to the output directory
#' @export
get_output_dir <- function(manifest, output_type = "processed", create = TRUE) {
  dir_key <- paste0(output_type, "_dir")

  if (!dir_key %in% names(manifest$outputs)) {
    stop("Output type not configured: ", output_type)
  }

  dir_path <- file.path(manifest$manifest_dir, manifest$outputs[[dir_key]])

  if (create && !dir.exists(dir_path)) {
    dir.create(dir_path, recursive = TRUE)
    message("Created output directory: ", dir_path)
  }

  return(dir_path)
}


#' Get QC parameters from manifest
#'
#' @param manifest Loaded manifest object
#' @return List of QC parameters
#' @export
get_qc_params <- function(manifest) {
  if (!"qc" %in% names(manifest)) {
    message("No QC parameters in manifest, using defaults")
    return(list(
      min_genes = 200,
      max_genes = 5000,
      max_mito_pct = 20,
      normalization_method = "LogNormalize",
      scale_factor = 10000,
      n_variable_features = 2000,
      selection_method = "vst"
    ))
  }
  return(manifest$qc)
}


#' Get integration parameters from manifest
#'
#' @param manifest Loaded manifest object
#' @return List of integration parameters
#' @export
get_integration_params <- function(manifest) {
  if (!"integration" %in% names(manifest)) {
    message("No integration parameters in manifest, using defaults")
    return(list(
      method = "harmony",
      batch_key = "study_source",
      dims = 30,
      resolution = 0.5
    ))
  }
  return(manifest$integration)
}


#' Get cell type mapping from manifest
#'
#' @param manifest Loaded manifest object
#' @return Named vector for cell type harmonization
#' @export
get_cell_type_map <- function(manifest) {
  if (!"harmonization" %in% names(manifest) ||
      !"cell_type_map" %in% names(manifest$harmonization)) {
    warning("No cell type mapping in manifest")
    return(NULL)
  }

  map <- unlist(manifest$harmonization$cell_type_map)
  return(map)
}


#' Apply cell type harmonization to a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param manifest Loaded manifest object
#' @param source_col Column containing original cell types
#' @param target_col Column to create with harmonized cell types
#' @return Seurat object with new harmonized column
#' @export
harmonize_cell_types <- function(seurat_obj, manifest, source_col, target_col = "cell_type_harmonized") {
  cell_type_map <- get_cell_type_map(manifest)

  if (is.null(cell_type_map)) {
    stop("No cell type mapping defined in manifest")
  }

  original_types <- seurat_obj[[source_col, drop = TRUE]]
  harmonized <- cell_type_map[as.character(original_types)]

  # Check for unmapped types
  unmapped <- is.na(harmonized)
  if (any(unmapped)) {
    unmapped_types <- unique(original_types[unmapped])
    warning("Unmapped cell types found: ", paste(unmapped_types, collapse = ", "))
    harmonized[unmapped] <- as.character(original_types[unmapped])
  }

  seurat_obj[[target_col]] <- harmonized
  return(seurat_obj)
}


#' List all studies in manifest with metadata
#'
#' @param manifest Loaded manifest object
#' @return Data frame with study information
#' @export
list_studies <- function(manifest) {
  studies_info <- lapply(names(manifest$data$studies), function(name) {
    study <- manifest$data$studies[[name]]
    data.frame(
      study = name,
      file = study$file,
      description = ifelse(is.null(study$description), NA, study$description),
      n_cells = ifelse(is.null(study$n_cells), NA, study$n_cells),
      cell_type_col = ifelse(is.null(study$cell_type_column), NA, study$cell_type_column),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, studies_info)
}


#' Print manifest summary
#'
#' @param manifest Loaded manifest object
#' @export
print_manifest_summary <- function(manifest) {
  cat("=====================================\n")
  cat("Manifest Summary\n")
  cat("=====================================\n")
  cat("Project:", manifest$project_name, "\n")
  cat("Author:", ifelse(is.null(manifest$author), "Not specified", manifest$author), "\n")
  cat("\n")

  cat("Studies:\n")
  studies <- list_studies(manifest)
  print(studies, row.names = FALSE)
  cat("\n")

  if (!is.null(manifest$qc)) {
    cat("QC Parameters:\n")
    cat("  Min genes:", manifest$qc$min_genes, "\n")
    cat("  Max genes:", manifest$qc$max_genes, "\n")
    cat("  Max mito %:", manifest$qc$max_mito_pct, "\n")
  }
  cat("=====================================\n")
}
