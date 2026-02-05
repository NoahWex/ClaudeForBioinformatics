# Exercise 2: QC Assessment & Filtering Strategy

**Objective**: Evaluate quality metrics and define filtering thresholds for the dataset.

**Time**: ~35 minutes

**Claude Features**: Edit tool, Bash (R execution), iterative refinement

---

## Background

Quality control is critical for single-cell analysis. Low-quality cells (damaged, doublets, empty droplets) can confound downstream analyses. In this exercise, you'll:
1. Visualize QC metrics
2. Make threshold decisions
3. Update your manifest with chosen parameters
4. Estimate filtering impact

## Your Task

### Part 1: Visualize QC distributions (10 min)

Ask Claude to create QC visualizations:

```
Please write an R script that creates QC visualizations for the Gray dataset:

1. Violin plot of nFeature_RNA by cell_type
2. Violin plot of percent.mt by cell_type
3. Scatter plot of nCount_RNA vs nFeature_RNA, colored by percent.mt
4. Histogram of percent.mt with proposed threshold line at 10%

Save plots to outputs/figures/qc_plots.pdf
Save the script to scripts/02_qc_visualization.R
```

**Observe**: How do the metrics differ by cell type? Are there obvious outliers?

### Part 2: Explore threshold options (15 min)

Now explore how different thresholds affect the data:

```
Please analyze how different QC thresholds would affect the dataset:

Test these threshold combinations:
1. Lenient: min_genes=200, max_genes=6000, max_mito=20
2. Standard: min_genes=500, max_genes=5000, max_mito=15
3. Stringent: min_genes=800, max_genes=4000, max_mito=10

For each, show:
- Total cells retained
- Cells retained per cell_type
- Cells retained per donor_id
- Any cell types or donors that lose >50% of cells

Create a comparison table.
```

**Discussion questions**:
- Which threshold balances quality vs. data retention?
- Are any cell types disproportionately affected?
- Do different donors have different quality profiles?

### Part 3: Update manifest with chosen thresholds (5 min)

Based on your analysis, update the manifest:

```
Please use the Edit tool to update my manifest.yaml with QC thresholds.

I want to use:
- min_genes: [your choice]
- max_genes: [your choice]
- max_mito_pct: [your choice]

Add a comment explaining why I chose these values.
```

**Note**: This is where you practice using the Edit tool rather than rewriting the whole file.

### Part 4: Document your decisions (5 min)

Create a QC decision document:

```
Please create outputs/qc_decisions.md documenting:
1. The thresholds I chose and why
2. Expected cell count after filtering (total and per cell type)
3. Any caveats or cell types to watch
4. Whether I should apply additional filtering later
```

---

## Claude Tools in Action

This exercise emphasizes the **Edit** tool:

```
# Instead of asking Claude to rewrite the whole manifest:
"Please use the Edit tool to change the max_mito_pct value from 20 to 15"

# Claude will use targeted edits rather than rewriting
```

Benefits of Edit:
- Preserves other parts of the file
- Shows exactly what changed
- Less error-prone than full rewrites

## Expected Findings

You might discover:

**QC by Cell Type**:
| Cell Type | Median Genes | Median Mito% | Notes |
|-----------|--------------|--------------|-------|
| basal cell | ~2500 | ~3% | High quality |
| fibroblast | ~2000 | ~4% | High quality |
| leukocyte | ~1500 | ~5% | Lower gene counts expected |
| lactocyte | ~3000 | ~8% | Higher mito (metabolically active) |

**Threshold Impact** (example):
| Threshold | Cells Kept | % Retained |
|-----------|------------|------------|
| Lenient | 51,500 | 97.8% |
| Standard | 48,200 | 91.5% |
| Stringent | 42,100 | 79.9% |

## Tips

1. **Consider biology** - Lactocytes may have higher mito% due to metabolic activity
2. **Check balance** - Ensure all cell types and donors retain cells
3. **Start lenient** - You can always filter more later
4. **Document everything** - Future you will thank present you

## Threshold Guidelines

General starting points (adjust based on your exploration):

| Metric | Guideline | Rationale |
|--------|-----------|-----------|
| min_genes | 200-500 | Remove empty droplets |
| max_genes | 2500-5000 | Remove potential doublets |
| max_mito% | 10-20% | Remove damaged cells |

Tissue-specific considerations:
- Breast epithelium may have higher baseline mito%
- Immune cells often have lower gene counts
- Stromal cells vary widely

## Checkpoint Verification

Before moving on, you should have:
- [ ] QC visualization plots saved
- [ ] Understanding of how thresholds affect different cell types
- [ ] Manifest updated with chosen QC parameters
- [ ] QC decisions documented with rationale
- [ ] Proficiency with the Edit tool

## Common Issues

### Removing too many cells of one type
- Solution: Consider cell type-specific thresholds or looser global thresholds

### Keeping obvious low-quality cells
- Solution: Visualize after filtering to verify improvement

### Unclear where to draw the line
- Solution: Look for natural breaks in the distribution; when in doubt, keep cells

## Next Steps

In Session 3, you'll apply these QC filters on HPC and run the full processing pipeline.
