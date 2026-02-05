# Exercise 2: Visualization

**Objective**: Create publication-quality visualizations of the integrated data.

**Time**: ~25 minutes

---

## Background

Good visualizations communicate your results effectively. In this exercise, you'll create UMAPs, feature plots, and summary figures suitable for publication or presentation.

## Your Task

### Part 1: Basic UMAP plots (10 min)

Ask Claude to create a visualization script:

```
Please create working/scripts/07_visualize.R that:

1. Loads the integrated object from outputs/integrated/integrated.rds
2. Creates the following UMAP plots:
   - Colored by study_source (to show mixing)
   - Colored by cell_type_harmonized (to show biology)
   - Colored by seurat_clusters (to show clustering)

3. Uses a consistent color palette
4. Adds proper titles and legends
5. Saves as outputs/figures/umaps_overview.pdf (multi-page)
```

### Part 2: Publication-quality formatting (10 min)

Refine the plots:

```
Please update the script to make publication-quality figures:

1. Use a clean theme (theme_classic or similar)
2. Set appropriate point size for the number of cells
3. Use colorblind-friendly palettes where possible
4. Remove redundant axis labels on UMAPs
5. Add cell count annotations to legends or subtitles
6. Set figure dimensions appropriate for publication (e.g., 8x6 inches)

Also create a combined panel figure with all three UMAPs side by side.
Save as outputs/figures/figure1_umaps.pdf
```

### Part 3: Additional visualizations (5 min)

Add supplementary figures:

```
Please add code to generate:

1. Feature plots for key marker genes:
   - EPCAM (epithelial)
   - PTPRC (immune)
   - COL1A1 (stromal)

2. A stacked bar plot showing cell type proportions per study

3. A table/heatmap of cell counts per study and cell type

Save all to outputs/figures/
```

---

## Color Palette Recommendations

### For cell types (categorical)
```r
library(RColorBrewer)
# For up to 12 categories
pal <- brewer.pal(12, "Paired")

# For more categories, use viridis discrete
library(viridisLite)
pal <- viridis(n_categories, option = "turbo")
```

### For studies (2-3 categories)
```r
study_colors <- c("gray" = "#E69F00", "pal" = "#56B4E9")
```

### For expression (continuous)
```r
# Seurat default or viridis
FeaturePlot(obj, features = "gene", cols = viridis(100))
```

## Example Code Patterns

### Clean UMAP
```r
DimPlot(obj, group.by = "cell_type", label = TRUE, label.size = 3) +
  theme_classic() +
  theme(
    legend.position = "right",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank()
  ) +
  ggtitle("Cell Types") +
  NoAxes()
```

### Multi-panel with patchwork
```r
library(patchwork)

p1 <- DimPlot(obj, group.by = "study_source") + ggtitle("Study")
p2 <- DimPlot(obj, group.by = "cell_type") + ggtitle("Cell Type")
p3 <- DimPlot(obj, group.by = "seurat_clusters") + ggtitle("Clusters")

combined <- p1 | p2 | p3
ggsave("combined.pdf", combined, width = 15, height = 5)
```

### Proportion bar plot
```r
library(ggplot2)

props <- obj@meta.data %>%
  group_by(study_source, cell_type_harmonized) %>%
  summarise(n = n()) %>%
  mutate(prop = n / sum(n))

ggplot(props, aes(x = study_source, y = prop, fill = cell_type_harmonized)) +
  geom_col() +
  theme_classic() +
  labs(y = "Proportion", x = "Study", fill = "Cell Type")
```

## Iterative Refinement with Claude

This is where Claude excels! Try prompts like:

```
The legend is overlapping the plot. Can you move it below?
```

```
The colors for T_cell and B_cell are too similar. Can you use more distinct colors?
```

```
Can you make the point size smaller? There are 130k cells and it's too crowded.
```

## Discussion Questions

1. What makes a UMAP "publication quality"?
2. How do you choose point size for large datasets?
3. When would you use a split UMAP vs overlay?

## Checkpoint Verification

Before moving on, verify:
- [ ] Overview UMAP PDF created
- [ ] Combined panel figure created
- [ ] Feature plots for marker genes
- [ ] Proportion plot created
- [ ] Figures are readable and well-formatted

## Next Steps

In [Exercise 3](03_report_generation.md), you'll assemble these figures into a reproducible report.
