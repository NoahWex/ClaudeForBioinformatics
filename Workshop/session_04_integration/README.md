# Session 4: Visualization, Reporting & Claude Customization

**Duration**: 3 hours
**Theme**: Analysis visualization, reproducible reports, customizing Claude

## Learning Objectives

By the end of this session, participants will be able to:

1. **Create publication-quality visualizations** - Generate UMAPs, violin plots, heatmaps
2. **Iterate on plots with Claude** - Refine aesthetics through dialogue
3. **Build reproducible reports** - Use Quarto for documentation
4. **Customize Claude behavior** - Write rules and understand hooks
5. **Apply workshop skills independently** - Design analyses for your own data

## Schedule

| Time | Segment | Type | Notes |
|------|---------|------|-------|
| 0:00 | Visualization strategies | Demo | Seurat plotting, ggplot2 |
| 0:15 | Creating visualizations | Hands-on | Exercise 1 |
| 0:55 | Break | - | 10 minutes |
| 1:05 | Report generation | Demo + Hands-on | Exercise 2 |
| 1:40 | Claude customization | Demo + Hands-on | Exercise 3 |
| 2:20 | Putting it together | Discussion | Full workflow review |
| 2:45 | Wrap-up & next steps | Discussion | Resources, Q&A |

## Prerequisites

Completed Session 3, or reset to checkpoint:
```bash
../recovery/sync_checkpoint.sh cp3_processed
```

## Key Claude Features Introduced

| Feature | Usage |
|---------|-------|
| `WebFetch` tool | Access documentation |
| Iterative refinement | Improve plots through dialogue |
| Rules files | Customize Claude behavior |
| Hooks | Intercept and modify commands |

## Exercises

1. **[Visualization](exercises/01_visualization.md)** - Create analysis plots
2. **[Report Generation](exercises/02_report_generation.md)** - Build Quarto report
3. **[Claude Customization](exercises/03_claude_customization.md)** - Write custom rules

## Templates

- `templates/report_template.qmd` - Quarto report template

## Checkpoint

At the end of this session, you should have:
- UMAP visualizations colored by cell type, donor, genotype
- Violin plots of marker genes
- Rendered HTML report
- Custom Claude rule file
- Understanding of full workshop workflow

## Visualization Strategies

### Dimensionality Reduction
```r
# UMAP colored by cell type
DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type")

# Split by genotype
DimPlot(seurat_obj, reduction = "umap", split.by = "genotype_1")
```

### Expression Plots
```r
# Violin plot of marker genes
VlnPlot(seurat_obj, features = c("KRT14", "KRT18", "CD3D"))

# Feature plot on UMAP
FeaturePlot(seurat_obj, features = c("EPCAM", "VIM"))
```

### Working with Claude

Ask Claude to help refine visualizations:
```
"Make the UMAP plot larger, use a colorblind-friendly palette,
and add a title showing the number of cells per cell type"

"Can you create a heatmap of the top 5 marker genes per cell type?"
```

## Claude Customization

### Rules Files

Create rules to change Claude's default behavior:

```markdown
# ~/.claude/rules/my_lab.md

## Data Conventions
- Always use ISO dates (YYYY-MM-DD) in filenames
- Save figures as both PDF and PNG
- Use `dalawson_lab` for HPC account
```

### Hooks

Hooks intercept commands before/after execution:
- PreToolUse: Validate before running
- PostToolUse: Log or modify results

See existing examples in `~/.claude/hooks/`

### Why Customize?

- **Consistency**: Same conventions across all sessions
- **Safety**: Prevent accidental destructive commands
- **Efficiency**: Skip repetitive instructions

## Resources

- [Seurat Visualization Vignette](https://satijalab.org/seurat/articles/visualization_vignette.html)
- [Quarto Documentation](https://quarto.org/)
- [Claude Code Customization](https://docs.anthropic.com/claude-code)
