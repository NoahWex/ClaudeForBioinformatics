# Exercise 1: QC Pipeline Script

**Objective**: Write an R script that performs quality control and normalization on Seurat objects.

**Time**: ~40 minutes

---

## Background

Before integration, each dataset needs:
1. Quality control filtering (removing low-quality cells)
2. Normalization (making expression values comparable)
3. Variable feature selection (identifying informative genes)
4. Scaling (centering and standardizing expression)

This script will be run on HPC in the next exercises.

## Your Task

### Part 1: Basic QC filtering (15 min)

Ask Claude to write the QC script:

```
Please create working/scripts/05_qc_normalize.R that:

1. Takes command-line arguments:
   - --study: study name (e.g., "gray")
   - --manifest: path to manifest.yaml

2. Loads the harmonized Seurat object from outputs/processed/

3. Calculates QC metrics:
   - percent.mt (mitochondrial genes)
   - percent.ribo (ribosomal genes)

4. Filters cells based on manifest QC parameters:
   - min_genes, max_genes
   - max_mito_pct

5. Reports how many cells were filtered

Use the argparse package for command-line arguments.
Read QC thresholds from the manifest using manifest_utils.R.
```

**Review**: Does the script handle edge cases (no cells passing filter)?

### Part 2: Normalization and feature selection (15 min)

Extend the script:

```
Please extend the script to add:

1. Normalization using SCTransform or LogNormalize
   - Use the method specified in manifest
   - If SCTransform, regress out percent.mt

2. Variable feature selection
   - Use n_variable_features from manifest

3. PCA computation
   - Compute 50 PCs

4. Save the processed object to outputs/processed/{study}_processed.rds

5. Generate a QC report figure showing:
   - Pre vs post filtering cell counts
   - Distribution of nFeature, nCount, percent.mt

Save figure to outputs/figures/{study}_qc_report.pdf
```

### Part 3: Make it robust (10 min)

Add error handling and logging:

```
Please update the script with:

1. Logging: Print timestamped messages for each major step
2. Error handling: Catch errors and report meaningful messages
3. Validation: Check that input file exists before processing
4. Memory: Add gc() calls after major steps to free memory

Also add a header comment block documenting:
- Purpose
- Usage example
- Required packages
- Input/output files
```

---

## Expected Script Structure

```r
#!/usr/bin/env Rscript
# 05_qc_normalize.R
# QC filtering and normalization for scRNA-seq data
#
# Usage: Rscript 05_qc_normalize.R --study gray --manifest manifest.yaml
#
# Inputs:
#   - Harmonized Seurat object from outputs/processed/
#   - Manifest with QC parameters
#
# Outputs:
#   - Processed Seurat object: outputs/processed/{study}_processed.rds
#   - QC report: outputs/figures/{study}_qc_report.pdf

library(Seurat)
library(argparse)
library(ggplot2)

# ... argument parsing ...

# ... load data ...

# ... QC filtering ...

# ... normalization ...

# ... save outputs ...
```

## Testing Locally

Before running on HPC, test with a small subset:

```
Please add a --test flag that:
- Subsets to 1000 random cells
- Runs faster for debugging
- Outputs to a test/ subdirectory
```

Test command:
```bash
Rscript working/scripts/05_qc_normalize.R --study gray --manifest working/manifest.yaml --test
```

## Discussion Questions

1. Why use command-line arguments instead of hardcoding values?
2. What's the advantage of reading parameters from the manifest?
3. How does SCTransform differ from LogNormalize?

## Checkpoint Verification

Before moving on, verify:
- [ ] Script accepts command-line arguments
- [ ] Script reads QC parameters from manifest
- [ ] Script handles errors gracefully
- [ ] Test run completes without errors
- [ ] Output files are created in correct locations

## Common Issues

### argparse not installed
```r
install.packages("argparse")
```
Or use the container which has it pre-installed.

### Memory issues locally
Use the `--test` flag or run on HPC with more memory.

### File paths
Remember that paths are relative to where you run the script, not where the script is located.

## Next Steps

In [Exercise 2](02_slurm_jobs.md), you'll wrap this script in a SLURM job for HPC execution.
