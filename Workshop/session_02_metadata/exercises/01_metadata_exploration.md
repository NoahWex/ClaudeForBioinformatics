# Exercise 1: Metadata Exploration

**Objective**: Systematically explore the Gray dataset metadata using Claude's tools.

**Time**: ~40 minutes

**Claude Features**: Read, Bash (R execution), Task/Explore agent

---

## Background

Before analyzing single-cell data, you need to understand what information is available. The Gray dataset contains rich metadata including cell type annotations, donor information, and clinical variables.

## Your Task

### Part 1: Get a metadata overview (10 min)

Ask Claude to explore the metadata structure:

```
I have a Seurat object at:
/share/crsp/lab/dalawson/share/3_Downloaded_Datasets/iHBCA_Reed_2024/component_studies/gray.rds

Please write and run R code that:
1. Loads the object
2. Shows all metadata column names
3. Reports the data type of each column
4. Shows dimensions (cells x genes)

Save the script to scripts/01_explore_metadata.R
```

**Observe**: How does Claude handle running R code? What tools does it use?

### Part 2: Explore categorical columns (15 min)

Now dig into the categorical metadata:

```
Please extend the analysis to:
1. For each categorical column (factors or characters), show:
   - Number of unique values
   - The actual values (if < 20 unique)
   - Count of cells per value
2. Focus especially on:
   - cell_type
   - cellSubtype
   - donor_id
   - genotype_1

Create a summary table.
```

**Discussion**:
- How many cell types are there?
- How many donors?
- What genotypes are represented?

### Part 3: Explore numeric columns (10 min)

Understand the QC metrics:

```
Now analyze the numeric metadata columns:
1. For nCount_RNA, nFeature_RNA, and percent.mt:
   - Show summary statistics (min, max, median, mean)
   - Show distribution by cell_type
2. Identify any other numeric columns
3. Flag any potential outliers
```

**Key questions**:
- What's the range of genes detected per cell?
- What's the mitochondrial percentage distribution?
- Do different cell types have different QC profiles?

### Part 4: Explore sample structure (5 min)

Understand the experimental design:

```
Please analyze the sample structure:
1. How many cells per donor?
2. What is the distribution of cell types per donor?
3. Are all genotypes represented across multiple donors?
4. Any donors with unusually few/many cells?
```

---

## Claude Tools in Action

During this exercise, notice how Claude uses different tools:

| Task | Tool Used | Why |
|------|-----------|-----|
| Writing R script | Write | Creates the script file |
| Running R code | Bash | Executes Rscript command |
| Reading results | Read | Views output files |
| Finding patterns | Grep | Searches for specific strings |

## Expected Findings

Your exploration should reveal something like:

**Dataset Overview**:
- 52,681 cells × 20,052 genes
- ~35 metadata columns

**Cell Types**:
| cell_type | Count |
|-----------|-------|
| basal cell | ~X,XXX |
| fibroblast | ~X,XXX |
| lactocyte | ~X,XXX |
| leukocyte | ~X,XXX |
| progenitor cell... | ~X,XXX |
| vascular lymphangioblast | ~X,XXX |

**Donors**: 16 unique (donor_id)

**Genotypes**: BRCA1, BRCA2, WT, RAD51C

## Tips

1. **Ask iteratively** - Start broad, then zoom in on interesting findings
2. **Request visualizations** - "Show me a histogram of percent.mt"
3. **Cross-tabulate** - "How do cell types distribute across genotypes?"
4. **Save your work** - Ask Claude to save outputs to files

## Document Your Findings

Create a summary document:

```
Please create outputs/metadata_summary.md containing:
1. Dataset dimensions
2. All metadata columns with descriptions
3. Cell type breakdown
4. Donor/genotype structure
5. QC metric ranges
6. Any interesting observations
```

## Checkpoint Verification

Before moving on, you should have:
- [ ] Script that explores metadata
- [ ] Understanding of all available columns
- [ ] Cell type and subtype counts
- [ ] Donor/genotype structure documented
- [ ] QC metric ranges identified

## Next Steps

In [Exercise 2](02_qc_assessment.md), you'll use this understanding to set QC thresholds.
