# Exercise 3: Report Generation

**Objective**: Create a reproducible Quarto report documenting your analysis.

**Time**: ~25 minutes

---

## Background

A good analysis report:
- Documents your methods and parameters
- Shows key results with figures
- Is reproducible (can be re-rendered)
- Serves as a record for collaborators and future you

Quarto (successor to R Markdown) lets you combine narrative, code, and outputs.

## Your Task

### Part 1: Create report structure (10 min)

Ask Claude to create a Quarto document:

```
Please create working/reports/integration_report.qmd that:

1. Has YAML header with:
   - Title: "Breast Tissue scRNA-seq Integration"
   - Author: (from manifest or placeholder)
   - Date: auto-generated
   - Format: html with table of contents

2. Includes sections for:
   - Introduction (project background)
   - Methods (brief description of pipeline)
   - Data Overview (cell counts, studies)
   - QC Summary (filtering statistics)
   - Integration Results (UMAPs)
   - Cell Type Composition
   - Conclusions

3. Uses R code chunks that reference our generated figures
   rather than re-running all analysis
```

### Part 2: Add content and figures (10 min)

Fill in the report:

```
Please update the report to:

1. Add actual content to each section (2-3 sentences each)

2. Include figures we generated:
   - QC report from Session 3
   - Integration UMAPs from Exercise 2
   - Cell type proportions

3. Add a table showing:
   - Study name
   - Original cell count
   - Post-QC cell count
   - Cell types present

4. Include the key parameters from manifest.yaml:
   - QC thresholds
   - Integration method and parameters

Use knitr::include_graphics() for existing figures.
```

### Part 3: Render the report (5 min)

```
Please help me render the report:

1. Show me the command to render the Quarto document
2. How do I render it on HPC using the container?
3. Where will the HTML output be saved?
```

---

## Quarto Basics

### YAML Header
```yaml
---
title: "My Analysis Report"
author: "Your Name"
date: today
format:
  html:
    toc: true
    toc-depth: 3
    code-fold: true
    theme: cosmo
---
```

### Code Chunks
````markdown
```{r}
#| label: load-data
#| echo: false
#| message: false

library(Seurat)
obj <- readRDS("../outputs/integrated/integrated.rds")
```
````

### Including Figures
````markdown
```{r}
#| label: fig-umap
#| fig-cap: "UMAP visualization of integrated data"
#| echo: false

knitr::include_graphics("../outputs/figures/figure1_umaps.pdf")
```
````

### Tables with kable
````markdown
```{r}
#| label: tbl-summary
#| tbl-cap: "Cell counts by study"

summary_table <- data.frame(
  Study = c("Gray", "Pal"),
  Cells = c(52700, 80000)
)
knitr::kable(summary_table)
```
````

## Report Template

See `templates/report_template.qmd` for a starting point.

## Rendering Commands

### Local (with Quarto installed)
```bash
quarto render reports/integration_report.qmd
```

### On HPC
```bash
singularity exec $CONTAINER quarto render reports/integration_report.qmd
```

### From R
```r
quarto::quarto_render("reports/integration_report.qmd")
```

## Tips for Good Reports

1. **Don't re-run expensive computations** - Load saved objects
2. **Use caching** for chunks that take time
3. **Fold code** (`code-fold: true`) for cleaner reading
4. **Include sessionInfo()** for reproducibility
5. **Date your report** to track versions

## Discussion Questions

1. When would you use Quarto vs plain R scripts?
2. How do you handle reports for analyses still in progress?
3. What's the minimum a report should include for reproducibility?

## Checkpoint Verification

Before moving on, verify:
- [ ] Quarto document created with all sections
- [ ] Figures successfully included
- [ ] Report renders without errors
- [ ] HTML output is readable and complete

## Common Issues

### Quarto not found
```bash
# Check if installed
quarto --version

# On HPC, may need to load module or use container
module load quarto
```

### Figure paths
Paths are relative to the .qmd file location, not the project root.

### Missing packages
Add a setup chunk that loads all required packages.

## Next Steps

In [Exercise 4](04_claude_customization.md), you'll learn to customize Claude's behavior for your lab.
