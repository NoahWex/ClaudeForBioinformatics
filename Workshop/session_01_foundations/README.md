# Session 1: Foundations & Project Setup

**Duration**: 2.5 hours
**Theme**: Claude basics, reproducible project structure, data loading

## Learning Objectives

By the end of this session, participants will be able to:

1. **Navigate Claude Code interface** - Use basic commands, get help, understand tool outputs
2. **Create project structure with Claude** - Use Write and Bash tools to scaffold directories
3. **Configure a manifest file** - Understand manifest-driven workflows and YAML syntax
4. **Load and explore Seurat objects** - Use Claude to write R code for data exploration
5. **Iterate with Claude** - Practice the prompt-response-refine cycle

## Schedule

| Time | Segment | Type | Notes |
|------|---------|------|-------|
| 0:00 | Introduction to Claude Code | Demo | Interface, /help, tool concepts |
| 0:25 | Project scaffolding | Hands-on | Exercise 1 |
| 1:00 | Break | - | 10 minutes |
| 1:10 | Loading & exploring data | Demo + Hands-on | Exercise 2 |
| 1:55 | Checkpoint review | Discussion | Verify outputs, Q&A |
| 2:15 | Wrap-up | - | Preview Day 2 |

## Prerequisites

Before this session, participants should have:
- [ ] Claude Code CLI installed and authenticated
- [ ] HPC3 account with SSH access
- [ ] Basic R familiarity (loading packages, reading data)
- [ ] Run `verify_setup.sh` with no critical failures

## Key Claude Features Introduced

| Feature | Usage |
|---------|-------|
| `/help` | Get help on Claude commands |
| `Write` tool | Create new files |
| `Bash` tool | Run shell commands (mkdir, ls, etc.) |
| `Read` tool | View file contents |
| `Glob` tool | Find files by pattern |

## Exercises

1. **[Project Scaffolding](exercises/01_project_scaffold.md)** - Create the workshop project structure
2. **[Data Loading](exercises/02_data_loading.md)** - Load and explore a Seurat object

## Checkpoint

At the end of this session, you should have:
- `working/` directory with standard project structure
- `manifest.yaml` with basic configuration
- `gray.rds` loaded and initial exploration complete
- Understanding of basic Claude Code workflow

Use `../recovery/sync_checkpoint.sh cp1_data_loaded` if you need to catch up.

## Common Issues

### Claude not responding
- Check API key configuration: `claude config`
- Verify internet connection

### File not found errors
- Use absolute paths when possible
- Check current working directory with `pwd`

### R package errors
- Most packages available in HPC container
- Local R may have different versions

## Resources

- [Claude Code Documentation](https://docs.anthropic.com/claude-code)
- [Seurat Vignettes](https://satijalab.org/seurat/articles/)
- Workshop shared utilities: `../shared/`
