# Exercise 2: Data Loading and Exploration

**Objective**: Use Claude to write R code for loading and exploring a Seurat object.

**Time**: ~45 minutes

---

## Background

Before integrating datasets, you need to understand what's in them. In this exercise, you'll use Claude to load a Seurat object and explore its structure, metadata, and quality metrics.

## Your Task

### Part 1: Load the data (10 min)

Ask Claude to help you write R code to load the Gray dataset:

```
I have a Seurat object at:
/share/crsp/lab/dalawson/share/3_Downloaded_Datasets/iHBCA_Reed_2024/component_studies/gray.rds

Please write an R script that:
1. Loads the required libraries (Seurat)
2. Reads the RDS file
3. Prints basic information about the object (number of cells, genes, assays)

Save the script to working/scripts/01_load_gray.R
```

**Watch for**: Claude using `Write` to create the R script.

### Part 2: Explore metadata (15 min)

Now ask Claude to help you explore the metadata:

```
Please extend the script to explore the metadata of the Seurat object:
1. List all metadata columns with their data types
2. For categorical columns (factors/characters), show the unique values and counts
3. Identify which column likely contains cell type annotations
4. Print a summary of cells per sample/batch if such a column exists

Add this to the existing script.
```

**Discussion**: Look at the output - what columns do you see? Which might be useful for integration?

### Part 3: QC metrics visualization (15 min)

Ask Claude to add QC metric exploration:

```
Please add code to calculate and visualize QC metrics:
1. Calculate percent mitochondrial genes (genes starting with "MT-")
2. Create a violin plot showing nFeature_RNA, nCount_RNA, and percent.mt
3. Save the plot to working/outputs/qc_metrics_gray.pdf

You can use Seurat's built-in plotting functions.
```

**Watch for**: How Claude handles the visualization code and file paths.

### Part 4: Document findings in manifest (5 min)

Based on your exploration, update the manifest:

```
Based on the exploration, please update working/manifest.yaml with:
- The actual number of cells in gray.rds
- The column name that contains cell type annotations
- Any batch/sample column identified

Use the Edit tool to update just the gray study section.
```

---

## Expected Outcome

You should have:
- `working/scripts/01_load_gray.R` - Complete exploration script
- `working/outputs/qc_metrics_gray.pdf` - QC visualization
- Updated `manifest.yaml` with gray dataset metadata

## Sample Exploration Output

Your exploration should reveal something like:

```
Seurat Object: 52,681 cells × 20,052 genes
Assays: RNA

Key metadata columns:
  - cell_type (character): 6 unique values
    (basal cell, fibroblast, lactocyte, leukocyte, progenitor..., vascular...)
  - cellSubtype (character): 19 unique values
  - donor_id (character): 16 unique values
  - genotype_1 (character): 4 values (BRCA1, BRCA2, WT, RAD51C)
  - percent.mt (numeric): mitochondrial percentage
```

## Tips for Working with Claude

1. **Start simple, then add complexity** - Get basic loading working before adding visualizations
2. **Share errors** - If code doesn't work, paste the error message to Claude
3. **Be specific about file paths** - Relative vs absolute paths matter
4. **Request comments** - Ask Claude to add comments explaining the code

## Discussion Questions

1. What cell types are present in the Gray dataset?
2. How many samples/batches are there?
3. What QC filtering might be needed based on the metrics?
4. How does the metadata structure compare to what you expected?

## Optional Extension

If you finish early, explore the clinical metadata:

```
Please analyze the clinical metadata in the Gray dataset:
1. How many donors have each genotype (BRCA1, BRCA2, WT, RAD51C)?
2. What is the distribution of cell types per genotype?
3. Create a table showing cells per donor, colored by genotype.
```

## Checkpoint Verification

Before moving on, verify you have:
- [ ] R script that loads and explores gray.rds
- [ ] QC visualization saved to outputs/
- [ ] Manifest updated with discovered metadata fields
- [ ] Notes on cell type column names and sample structure

## Next Steps

In Session 2, you'll dive deeper into the metadata structure and establish QC thresholds for filtering.
