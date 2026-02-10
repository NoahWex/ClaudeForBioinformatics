# Graph Memory System

A plan and knowledge management toolkit for tracking work across Claude Code sessions.

## Overview

The graph memory system provides:

- **Plans**: Structured task documents with dependency tracking and graph queries
- **Knowledge Base**: Persistent storage for lessons, research, and patterns
- **Graph Queries**: Find related plans, blocked work, file references, and clusters

This enables Claude to understand project context across sessions without re-exploring the codebase each time.

## Installation

### Option A: Project-Embedded (Recommended)

Copy tools into your project and create data directories:

```bash
# Copy the toolkit
cp -r graph-memory/ your-project/

# Create plan directories
mkdir -p your-project/project/.plans/{active_plans,.archive}

# Create knowledge directories
mkdir -p your-project/project/.knowledge/{lessons,research,patterns}
```

The CLI discovers `.plans/` by walking up from the current directory.

### Option B: Global Installation

Install once in `~/.claude/` for use across projects:

```bash
# Copy to global location
cp -r graph-memory/ ~/.claude/graph-memory/

# Add to PATH (optional)
echo 'export PATH="$HOME/.claude/graph-memory:$PATH"' >> ~/.bashrc
```

Then create `.plans/` and `.knowledge/` in each project as needed.

### Environment Variables

Override discovery with environment variables:

```bash
export PLAN_DIR=/path/to/project/.plans
export KNOWLEDGE_DIR=/path/to/project/.knowledge
```

## Quick Start

### Initialize a Project

```bash
# Create directories
mkdir -p project/.plans/{active_plans,.archive}
mkdir -p project/.knowledge/{lessons,research,patterns}

# Copy templates
cp graph-memory/templates/plan_template.md project/.plans/active_plans/my_area/first_plan.md

# Edit the plan with your details
```

### Create Your First Plan

Plans live in `.plans/active_plans/` organized by area:

```
.plans/
├── active_plans/
│   ├── integration/
│   │   └── api_refactor.md
│   └── analysis/
│       └── performance_audit.md
├── .archive/
└── .index.json
```

### Query Plans

```bash
# Plans ready to work on (not blocked)
./plan query ready

# Blocked plans and their blockers
./plan query blocked

# Plans referencing a file
./plan query touches config/settings.yaml

# All pending tasks across plans
./plan query tasks
```

## Commands Reference

### Plan Queries

| Command | Description |
|---------|-------------|
| `plan query ready` | Plans with no unmet blockers |
| `plan query blocked` | Plans and what blocks them |
| `plan query context <id>` | Cluster context for a plan |
| `plan query downstream <id>` | Plans depending on this one |
| `plan query touches <path>` | Plans referencing a file pattern |
| `plan query anchors` | MNN anchors (cross-modal links) |
| `plan query communities` | All clusters and members |
| `plan query history <id>` | Succession chain for a plan |
| `plan query tasks` | Pending tasks across all plans |
| `plan query tasks --all` | Include completed tasks |

### Graph Visualization

```bash
plan graph              # Table format (default)
plan graph --ascii      # ASCII DAG view
plan graph --mermaid    # Mermaid markdown
plan graph --dot        # Graphviz DOT
```

### Index Operations

```bash
plan validate <file>    # Check plan structure
plan index rebuild      # Regenerate .index.json
```

### Sync Operations

```bash
plan sync [path]        # Sync draft to active_plans + sequence YAML
plan nudge [path]       # Write sync reminder for next session
plan reconcile [seq]    # Cross-reference plan files vs sequence YAML
plan reconcile --fix <seq>  # Fix status mismatches (YAML authoritative)
```

`plan sync` is the primary way agents share plan state. It copies the draft to
`active_plans/`, rebuilds the index, reconciles with sequence YAML, and removes
the matching breadcrumb from `.pending_sync.d/`.

### Knowledge Operations

```bash
plan knowledge list              # All entries by type
plan knowledge search <query>    # Search by tag/content
plan knowledge add <path>        # Index a new entry
plan knowledge show <id>         # Display entry details
plan knowledge rebuild           # Rebuild .index.yaml
```

### Plan Completion

```bash
plan complete <id>      # Debrief and archive
```

This prompts for lessons learned before archiving.

### Sequence Operations

| Command | Description |
|---------|-------------|
| `plan sequence list` | All sequences with progress |
| `plan sequence status <name>` | Full status with vision and plans |
| `plan sequence vision <name>` | Vision section only |
| `plan sequence next <name>` | Next ready plan (JSON with `--json`) |
| `plan sequence validate <name>` | Check for errors and cycles |
| `plan sequence scope <name>` | Vision with decisions |

### Session Management

| Command | Description |
|---------|-------------|
| `plan enter <seq/plan>` | Claim plan, show vision, create .md |
| `plan pause "notes"` | Save notes, release ownership |
| `plan seq-complete "summary"` | Mark complete, unblock downstream |

## Plan Document Format

### Frontmatter (YAML)

```yaml
---
area: integration           # Directory grouping
name: api_refactor          # Plan identifier (area/name is unique)
version: 1                  # Increment on major revisions
updated: 2026-02-04         # Last modified
supersedes: null            # Previous plan ID if replacing
blocked_by: []              # Must complete first
blocks: []                  # Waiting on this
status: in_progress         # queued | in_progress | blocked | complete
draft_path: /path/to/draft  # Absolute path to this draft file (set by agent)
jobs: []                    # Optional SLURM job IDs
---
```

### Required Sections

| Section | Content Type | Purpose |
|---------|--------------|---------|
| `## Inputs` | Pointer list | Files consumed |
| `## Outputs` | Pointer list | Files produced |
| `## Config` | Pointer list | Configuration files |
| `## Scripts` | Pointer list | Implementation code |
| `## Tasks` | Checklist | Work items |
| `## Current State` | Freeform | Status description |
| `## Discussion` | Dated entries | Context and decisions |
| `## Next Actions` | Freeform | Immediate steps |

### Pointer Syntax

```markdown
- `path/to/file.ext` — description
- `path/to/file.ext:42` — with line reference
- `path/{study}/outputs/*.rds` — with wildcards
```

### Task Syntax

```markdown
- [ ] Pending task
- [x] Completed task
```

### Dated Discussion

```markdown
## Discussion

### 2026-02-04: Initial planning

Context and decisions made today.

### 2026-02-05: Progress update

What changed and why.
```

## Knowledge Entry Types

### Lessons

Hard-won debugging knowledge. For problems you never want to debug again.

```yaml
---
id: crsp_stale_handle
type: lesson
source: debugging
tags: [crsp, filesystem, hdf5]
severity: high
---

# CRSP Stale File Handle

## Problem
...
```

### Research

Literature summaries and method comparisons.

```yaml
---
id: milo_differential_abundance
type: research
source: literature
tags: [milo, differential_abundance, single_cell]
zotero_key: ABCD1234
---
```

### Patterns

Validated reusable code patterns.

```yaml
---
id: slurm_array_pattern
type: pattern
source: experiment
tags: [slurm, hpc, array_jobs]
---
```

## Sequence Coordination

Sequences coordinate multiple plans toward a shared end goal. While individual plans track tasks within a single work unit, sequences manage the dependencies and flow between them.

See **[14_SEQUENCE_COORDINATION.md](14_SEQUENCE_COORDINATION.md)** for the full standalone guide.

### Quick Reference

```bash
# Overview
plan sequence list              # All sequences with progress
plan sequence status <name>     # Full status with vision + plans
plan sequence vision <name>     # Vision only (quick re-orient)

# Session management
plan enter <sequence>/<plan>    # Claim plan, show vision, create .md
plan pause "notes"              # Save progress notes, release claim
plan seq-complete "summary"     # Mark complete, unblock downstream

# Validation
plan sequence validate <name>   # Check for errors, cycles
plan sequence next <name>       # Next ready plan (for automation)
```

### Sequence YAML

Sequences live in `.plans/sequences/<name>.yaml`:

```yaml
name: my_pipeline
status: active
auto_progress: prompt     # prompt | silent | disabled
vision:
  end_goal: "What success looks like"
  why_it_matters: "Connection to project goals"
  current_state: "Where things stand"
  decisions:
    - date: "2026-02-05"
      fork: "What was decided"
      chose: "What was chosen"
      rationale: "Why"
  exit_criteria:
    - "Verifiable criterion"
plans:
  - id: area/plan_name
    summary: "What this plan does"
    status: queued
    depends_on: [upstream_plan]
```

See `graph-memory/templates/sequence_template.yaml` for a full template.

### When to Use Sequences vs Individual Plans

| Situation | Use |
|-----------|-----|
| Single task, no dependencies | Individual plan |
| 2-3 related tasks, simple order | Individual plans with `blocked_by` |
| 4+ plans with branching dependencies | Sequence |
| Shared end goal across work areas | Sequence |
| Parallel sessions on related work | Sequence (ownership tracking) |

## Architecture

### Graph Construction

The system builds a dependency graph from plan metadata:

1. **Explicit edges**: `blocked_by` and `blocks` frontmatter
2. **Dataflow edges**: Output paths matching input paths (inferred)
3. **Resource edges**: Shared file references (lower weight)
4. **Succession edges**: `supersedes` relationships

### SNN Clustering

Plans are grouped by Shared Nearest Neighbor (SNN) clustering:

- Plans in the same area form workstreams
- Clusters pool inputs, outputs, tasks, and discussion
- Use `plan query context <id>` to see cluster context

### Index Structure

`.index.json` contains:

```json
{
  "active_plans": { "area/name": { "version": 1, "status": "..." } },
  "archive": { "2026-02-04__area_name_v1": { ... } },
  "clusters": { "area_single": { "members": [...], "pooled_inputs": [...] } },
  "graphs": [{ "from": "...", "to": "...", "type": "dependency" }]
}
```

## Verification

After installation, run the verification tests:

```bash
cd graph-memory
./tests/verify_install.sh
```

Expected output:

```
Graph Memory System Verification
=================================
[OK] Python 3 available
[OK] CLI script exists
[OK] CLI invocation
...
All checks passed!
```

## Integration with CLAUDE.md

Add graph memory commands to your project's CLAUDE.md:

```markdown
## Graph Memory

Track plans and knowledge with the graph memory system:

\`\`\`bash
./plan query ready      # Plans to work on
./plan query tasks      # Pending tasks
./plan knowledge list   # Lessons and patterns
\`\`\`
```

## Tips

1. **Prefer updating existing knowledge** over creating new entries
2. Use `plan complete <id>` for structured debriefs
3. Run `plan index rebuild` after manual edits
4. Use `plan query touches <path>` to find related plans before major changes
5. Check `plan query blocked` regularly to unblock work
