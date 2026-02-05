# Breast Tissue scRNA-seq Analysis

Analysis of Gray et al. breast tissue single-cell RNA-seq data.

## Project Structure

```
working/
├── scripts/    # R analysis scripts
├── data/       # Input data (symlinks or references)
├── outputs/    # Processed data, organized by date
├── reports/    # Quarto/Rmd analysis reports
└── figures/    # Publication-quality figures
```

## Data Source

Gray et al. breast tissue scRNA-seq from iHBCA (Reed et al. 2024)
- ~53,000 cells
- 16 donors
- Cell types: basal, fibroblast, lactocyte, leukocyte, progenitor, vascular

## Usage

Scripts are numbered in execution order:
- `01_*.R` - Data loading and exploration
- `02_*.R` - QC and filtering
- `03_*.R` - Processing and visualization
