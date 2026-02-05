# Graph Memory System

A plan and knowledge management toolkit for Claude Code projects.

## Quick Start

```bash
# From your project root
./graph-memory/plan --help

# Query available plans
./graph-memory/plan query ready

# Validate a plan file
./graph-memory/plan validate path/to/plan.md

# Rebuild the index
./graph-memory/plan index rebuild
```

## Installation

### Option A: Project-Embedded (Recommended)

Copy the `graph-memory/` directory into your project and create the data directories:

```bash
cp -r graph-memory/ your-project/
mkdir -p your-project/project/.plans/{active_plans,.archive}
mkdir -p your-project/project/.knowledge/{lessons,research,patterns}
```

### Option B: Global Installation

Copy to `~/.claude/` for use across projects:

```bash
cp -r graph-memory/ ~/.claude/graph-memory/
```

Then create `.plans/` and `.knowledge/` in each project as needed.

## Environment Variables

- `PLAN_DIR` - Override plans directory location
- `KNOWLEDGE_DIR` - Override knowledge directory location

## Commands

See `./plan --help` for full command reference.

## Documentation

See `13_GRAPH_MEMORY_SYSTEM.md` in the package root for comprehensive documentation.
