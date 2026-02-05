# CLAUDE.md

## Project Overview

**[Project Name]** - [Brief 1-2 sentence description of what this project does and what data it works with.]

**Dataset**: [Key dataset characteristics - samples, cells, conditions, etc.]

## Current Phase

| Phase | Location | Status |
|-------|----------|--------|
| **[Phase 1]** | `project/[location]/` | [Status] |
| **[Phase 2]** | `project/[location]/` | [Status] |

## Orientation: Understanding Current State

**To understand what's being worked on:**
1. Check `project/[analyses_path]/config/experiment_catalog.yaml` for active experiments
2. Look at experiments with `status: development` for ongoing work
3. Each experiment has `experiment_metadata.yaml` describing its purpose

**To understand available data:**
1. `project/config/raw_data_manifest.yaml` - source of truth (raw data + processing outputs)
2. `project/config/reference_datasets.yaml` - external references and models
3. `project/[analyses_path]/config/analysis_outputs.yaml` - contributions from completed experiments

**To understand project decisions:**
1. [Location of decision documentation]
2. Experiment README.md files document rationale and findings

## Behavioral Rules

Rules in `.claude/rules/` provide context-specific guidance:
- `[rule-name].md` - [brief description]
- `[folder]/[rule-name].md` - [brief description]

## Quick Reference

### Manifests

| Manifest | Purpose |
|----------|---------|
| `raw_data_manifest.yaml` | Source of truth: raw data + processing outputs |
| `reference_datasets.yaml` | External references and models |
| `container_registry.yaml` | HPC containers and environments |
| `experiment_catalog.yaml` | Experiment registry |
| `analysis_outputs.yaml` | Contributions from completed experiments |

### Key Data

- **Primary dataset**: [description and location]
- **Key statistics**: [e.g., N cells, N samples]
- **Primary objects**: [key data object paths]

### Analysis Workflow

```r
# Start new experiment
source("project/[path]/scripts/init_experiment.R")
init_experiment("name", "ResearchArea")

# Complete experiment
source("project/[path]/scripts/complete_experiment.R")
complete_experiment("name_YYYYMMDD")

# Deprecate if unproductive
source("project/[path]/scripts/deprecate_experiment.R")
deprecate_experiment("name_YYYYMMDD", reason = "...")
```

Research areas: `[Area1]`, `[Area2]`, `[Area3]`, `[Area4]`

## HPC Environment

**Cluster**: [Cluster name] | **Account**: `[account]` | **Partition**: `[partition]`

**Containers**:
- **R_env**: `[/path/to/R/container.sif]`
- **Python_env**: `[/path/to/Python/container.sif]`

**Bind mounts**:
```bash
-B [/path/to/data] -B [/path/to/scratch]
```

## External References

- **[Reference 1]**: `[path or URL]` - [description]
- **[Reference 2]**: `[path or URL]` - [description]

## Scientific Context

**Goal**: [Brief description of scientific goals]

**Key concepts**:
- **[Concept 1]**: [description]
- **[Concept 2]**: [description]

## Skills

This project has custom skills available - use `Skill` tool to invoke:
- `project-status` - Check project health and progress
- `project-discovery` - Understand project structure
