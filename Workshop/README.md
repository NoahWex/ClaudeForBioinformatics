# Claude Code for Bioinformatics Workshop

A hands-on workshop teaching Claude Code through a scRNA-seq analysis vignette.

## What This Is

**This is a Claude Code workshop**, not a single-cell analysis course. The bioinformatics serves as a realistic context to learn:
- How to work with Claude's tools (Read, Write, Edit, Bash, Grep)
- How to run commands on HPC through Claude
- How to iterate with Claude to solve problems
- How to customize Claude's behavior for lab workflows

## Overview

**Duration**: 4 days, 2.5-3 hours per session
**Audience**: Graduate students with R proficiency and basic CLI experience
**Infrastructure**: UCI HPC3 cluster
**Dataset**: Gray et al. breast tissue scRNA-seq (~53K cells)

### What You'll Learn

| Day | Topic | Claude Features |
|-----|-------|-----------------|
| 1 | Foundations | Write, Read, Bash, Glob |
| 2 | Metadata Exploration | Grep, Edit, Task/Explore |
| 3 | HPC Processing | SLURM jobs, /hpc skill |
| 4 | Visualization & Customization | WebFetch, Rules, Skills |

## Prerequisites

### Required
- [ ] UCI HPC3 account with SSH access
- [ ] Claude Code CLI installed and configured
- [ ] Familiarity with R (loading packages, basic syntax)
- [ ] Basic command line experience (cd, ls, pwd)

### Helpful
- [ ] Seurat package exposure (not required - we'll learn as we go)

## Getting Started

### 1. Verify Your Setup

```bash
cd Workshop/recovery
./verify_setup.sh
```

### 2. Prepare Your Workspace

```bash
cp -r Workshop/participant_workspace Workshop/working
cd Workshop/working
```

### 3. Start Session 1

Open the session directory and follow the README:
```
session_01_foundations/README.md
```

## Directory Structure

```
Workshop/
├── README.md                    # This file
│
├── session_01_foundations/      # Day 1: Claude basics
├── session_02_metadata/         # Day 2: Exploration & QC
├── session_03_hpc/             # Day 3: HPC processing
├── session_04_integration/     # Day 4: Viz & customization
│
├── shared/                     # Common resources
├── data/                       # Dataset info
├── participant_workspace/      # Template for your work
└── recovery/                   # Checkpoint recovery
```

## Checkpoint System

Each session ends with a checkpoint. If you fall behind:

```bash
./recovery/sync_checkpoint.sh cp1_data_loaded ./working
```

| Checkpoint | End of | Description |
|------------|--------|-------------|
| `cp1_data_loaded` | Day 1 | Project + loaded data |
| `cp2_explored` | Day 2 | Metadata explored, QC set |
| `cp3_processed` | Day 3 | Processed on HPC |
| `cp4_complete` | Day 4 | Visualizations + report |

## Dataset

**Gray et al.** breast tissue scRNA-seq from iHBCA:
- 52,681 cells
- 6 cell types (basal, fibroblast, lactocyte, leukocyte, progenitor, vascular)
- 16 donors with BRCA1/BRCA2/WT genotypes
- Rich clinical metadata

Location: `/share/crsp/lab/dalawson/share/3_Downloaded_Datasets/iHBCA_Reed_2024/component_studies/gray.rds`

## Resources

- [Claude Code Docs](https://docs.anthropic.com/claude-code)
- [UCI HPC3 Guide](https://rcic.uci.edu/hpc3/)
- [Seurat Tutorials](https://satijalab.org/seurat/)

---

*Developed for the Dawson Lab, UCI*
