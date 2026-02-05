# Claude Code Best Practices for Scientific Computing

**Extracted from**: Spatial HBCA Resource Project
**Version**: 1.1
**Date**: 2026-02-02

This guide documents infrastructure patterns and frameworks developed for working with Claude Code in scientific computing environments, particularly for single-cell genomics and HPC-based analysis pipelines.

## Overview

These best practices cover:

1. **Project Configuration** - How to structure `CLAUDE.md` and Claude Code settings
2. **Contextual Rules** - Domain-specific guidance that activates based on file paths
3. **Anti-Sycophancy Framework** - Maintaining objectivity and scope discipline
4. **Manifest-Driven Development** - YAML-based configuration management
5. **Experiment Lifecycle** - Managing analysis experiments with version control
6. **HPC Infrastructure** - SLURM patterns and Singularity container gotchas
7. **Custom Skills and Agents** - Extending Claude Code's capabilities
8. **Deprecation Policy** - Preserving context when retiring code
9. **Graph Memory System** - Plan and knowledge tracking across sessions

## Quick Start

1. Copy the `templates/` directory to your project
2. Adapt `CLAUDE.md` to your project's specifics
3. Configure `.claude/settings.local.json` for your environment
4. Add contextual rules to `.claude/rules/` as needed

## File Organization

```
ClaudeForBioinformatics_NW/
├── README.md                        # This file
├── GETTING_STARTED.md               # Quick start guide
├── 01_CLAUDE_MD_STRUCTURE.md        # Project instruction file patterns
├── 02_SETTINGS_AND_HOOKS.md         # Claude Code configuration
├── 03_RULES_FRAMEWORK.md            # Contextual rules system
├── 04_TRACE_FRAMEWORK.md            # Anti-sycophancy and objectivity
├── 05_MANIFEST_PATTERNS.md          # YAML manifest design
├── 06_EXPERIMENT_LIFECYCLE.md       # Analysis experiment management
├── 07_HPC_INFRASTRUCTURE.md         # SLURM and container patterns
├── 08_SKILLS_AND_AGENTS.md          # Custom extensions
├── 09_DEPRECATION_POLICY.md         # Handling retired code
├── 10_HPC_PLUGIN.md                 # HPC submission plugin (NEW)
├── 11_INTERFACE_AWARE_PLANNING.md   # Multi-component planning (NEW)
├── 12_SANDBOX_ENVIRONMENT.md        # Learning sandbox guide
├── 13_GRAPH_MEMORY_SYSTEM.md        # Plan and knowledge tracking (NEW)
├── graph-memory/                    # Portable toolkit (NEW)
│   ├── plan                         # CLI entry point
│   ├── tools/                       # Python modules
│   ├── templates/                   # Plan/knowledge templates
│   └── tests/                       # Verification tests
├── templates/
│   ├── CLAUDE.md                    # Project instruction template
│   ├── settings.local.json          # Settings template
│   ├── trace-framework.md           # Objectivity framework
│   ├── experiment_metadata.yaml     # Experiment config template
│   ├── slurm_template_R.sh          # SLURM R job template
│   └── slurm_template_python.sh     # SLURM Python job template
└── sandbox/                         # Hands-on learning environment
    ├── README.md                    # Sandbox quickstart
    ├── config/                      # Configurable settings
    ├── setup-guides/                # Step-by-step installation
    ├── global-setup/                # Global Claude Code config + production hpc-toolkit
    │   ├── hpc-toolkit/             # Unified HPC CLI (bin/hpc, gate.sh, config.sh.template)
    │   ├── skills/hpc/              # HPC skill
    │   └── rules/                   # Global rules
    ├── project-template/            # Project structure template
    └── exercises/                   # Guided learning exercises
```

## Core Principles

### 1. Source of Truth Philosophy

Every data path, configuration parameter, and external reference should have exactly ONE authoritative source:

- **Raw data paths** → `raw_data_manifest.yaml`
- **Reference datasets** → `reference_datasets.yaml`
- **Container configs** → `container_registry.yaml`
- **Experiment tracking** → `experiment_catalog.yaml`

Scripts load from manifests, never hardcode paths.

### 2. Contextual Rules Over Generic Instructions

Instead of putting all guidance in `CLAUDE.md`, use path-triggered rules:

```yaml
---
paths: project/02_Analyses/**/*.R
---
# This rule only activates when working on R analysis scripts
```

### 3. Explicit Scope Discipline

Before modifying anything:
- State what's in-scope vs out-of-scope
- Ask permission before extending scope
- Implement only what's explicitly requested

### 4. Anti-Sycophancy

Configure Claude to:
- Prioritize factual accuracy over validation
- Avoid unwarranted praise
- Disagree when evidence supports it
- Use neutral acknowledgments ("Understood", "Noted")

### 5. Deprecation Over Deletion

When retiring code:
- Move to `.deprecated/` directory
- Always create `DEPRECATED.md` explaining why
- Include de-deprecation instructions
- Preserve context for future reference

## Adaptation Guide

### For Different Scientific Domains

Replace domain-specific content:
- Cell type annotations → your domain's categories
- Anatomical positions → your experimental conditions
- Compartments → your biological groupings

### For Different HPC Environments

Update in templates:
- Container paths
- Partition/account names
- Bind mount paths
- Module load commands

### For Different Languages

The patterns work for:
- R + Python (primary examples)
- Jupyter notebooks
- Shell scripts
- Any script-based analysis

## Dependencies

Claude Code features used:
- Custom instructions (`CLAUDE.md`)
- Settings and hooks (`.claude/settings.local.json`)
- Rules system (`.claude/rules/`)
- Skills (`.claude/skills/`)
- Agents (`.claude/agents/`)

## Changelog

### v1.2 (2026-02-04)
- **Graph Memory System**: Added plan and knowledge tracking toolkit (`graph-memory/`)
- **Documentation**: New `13_GRAPH_MEMORY_SYSTEM.md` comprehensive guide
- **Templates**: Plan and knowledge entry templates
- **Verification**: Automated verification tests for installation

### v1.1 (2026-02-02)
- **HPC Toolkit**: Sandbox now ships the production `hpc-toolkit` (unified CLI with 8 subcommands, gate hook, config.sh)
- **Deprecated**: Old `hpc-sandbox` plugin moved to `.deprecated/`
- **Setup Guides**: Updated `03_plugin_install.md` and `04_verify_setup.md` for toolkit installation
- **Configurable**: Added `config.sh.template` so lab members only edit one file for their environment
- **GETTING_STARTED**: Added global setup section distinguishing project-level vs user-level config

## License

These patterns are shared for community benefit. Adapt freely for your projects.
