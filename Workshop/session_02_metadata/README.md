# Session 2: Metadata Exploration & QC Assessment

**Duration**: 2.5 hours
**Theme**: Deep-dive into dataset structure, quality metrics, and filtering decisions

## Learning Objectives

By the end of this session, participants will be able to:

1. **Explore metadata systematically** - Use Claude to examine all columns and their distributions
2. **Assess data quality** - Identify QC metrics and set filtering thresholds
3. **Use pattern matching** - Leverage Claude's Grep tool to search code and documentation
4. **Use Task agents** - Delegate exploratory work to specialized agents
5. **Make data-driven decisions** - Choose appropriate parameters based on exploration

## Schedule

| Time | Segment | Type | Notes |
|------|---------|------|-------|
| 0:00 | Metadata landscape | Demo | Review Session 1 outputs |
| 0:20 | Deep metadata exploration | Hands-on | Exercise 1 |
| 1:00 | Break | - | 10 minutes |
| 1:10 | QC metrics assessment | Demo | Understanding quality thresholds |
| 1:35 | Filtering strategy | Hands-on | Exercise 2 |
| 2:10 | Checkpoint review | Discussion | Verify understanding |
| 2:30 | Wrap-up | - | Preview Day 3 |

## Prerequisites

Completed Session 1, or reset to checkpoint:
```bash
../recovery/sync_checkpoint.sh cp1_data_loaded
```

## Key Claude Features Introduced

| Feature | Usage |
|---------|-------|
| `Grep` tool | Search for patterns in files |
| `Edit` tool | Modify existing files |
| `Task` agent (Explore) | Delegated codebase exploration |
| Iterative questioning | Refine understanding through dialogue |

## Exercises

1. **[Metadata Exploration](exercises/01_metadata_exploration.md)** - Systematically explore all metadata
2. **[QC Assessment](exercises/02_qc_assessment.md)** - Evaluate quality metrics and set thresholds

## Checkpoint

At the end of this session, you should have:
- Complete understanding of all metadata columns
- QC threshold decisions documented
- Updated manifest with QC parameters
- Understanding of Grep, Edit, and Task tools

Use `../recovery/sync_checkpoint.sh cp2_explored` if you need to catch up.

## Key Concepts

### Gray Dataset Metadata Overview

The Gray et al. dataset contains rich metadata:

**Cell annotations:**
- `cell_type`: 6 broad categories (basal cell, fibroblast, lactocyte, leukocyte, progenitor, vascular lymphangioblast)
- `cellSubtype`: 19 finer-grained subtypes

**Sample information:**
- `donor_id`: 16 unique donors
- `genotype_1`: BRCA1, BRCA2, WT, RAD51C

**Clinical metadata:**
- `surgeryType`, `menopause`, `numberOfPregnancies`, etc.

**QC metrics:**
- `nCount_RNA`: Total UMI counts per cell
- `nFeature_RNA`: Number of detected genes per cell
- `percent.mt`: Mitochondrial read percentage

### QC Decision Framework

1. **Low-quality cells** - Few genes, high mitochondrial %
2. **Doublets** - Unusually high gene/UMI counts
3. **Empty droplets** - Very few counts

Typical thresholds (adjust based on your data):
- Min genes: 200-500
- Max genes: 2500-5000 (depends on sequencing depth)
- Max mito %: 5-20% (tissue-dependent)

### Exploration Strategies

**Systematic approach:**
1. List all columns and their types
2. Check unique values for categorical columns
3. Summarize distributions for numeric columns
4. Identify potential issues or anomalies

**Claude-assisted exploration:**
```
"Show me the distribution of nFeature_RNA grouped by cell_type"
"What percentage of cells have percent.mt > 10?"
"Are there any donors with unusually low cell counts?"
```

## Common Issues

### Understanding column relationships
- Some columns are redundant (ontology IDs vs labels)
- Focus on columns you'll actually use

### Setting appropriate thresholds
- Thresholds depend on tissue type and experimental design
- Use visualization (histograms) to inform decisions
- Start lenient, tighten if needed

### Missing or NA values
- Check for completeness in key columns
- Decide how to handle missing clinical metadata

## Resources

- Seurat QC tutorial: https://satijalab.org/seurat/articles/pbmc3k_tutorial
- Workshop manifest: `../data/workshop_data.yaml`
- Manifest utilities: `../shared/manifest_utils.R`
