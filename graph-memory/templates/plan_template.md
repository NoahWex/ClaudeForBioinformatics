---
area: AREA_NAME               # Directory grouping (e.g., "integration", "analysis")
name: PLAN_NAME               # Plan identifier (area/name is unique ID)
version: 1                    # Increment on major revisions
updated: YYYY-MM-DD           # Last modified date
supersedes: null              # Previous plan ID if this replaces another
blocked_by: []                # Plan IDs that must complete first
blocks: []                    # Plan IDs waiting on this
status: queued                # queued | in_progress | blocked | complete
jobs: []                      # Optional: SLURM job IDs if applicable
---

# Plan Title

Brief description of what this plan accomplishes.

## Inputs

Files and data this plan consumes:

- `path/to/input/file.ext` — description of input
- `path/to/config.yaml` — configuration file

## Outputs

Files and artifacts this plan produces:

- `path/to/output/file.ext` — description of output

## Config

Configuration files and settings:

- `path/to/config.yaml:50` — specific line reference if needed

## Scripts

Implementation code:

- `path/to/script.py` — main implementation
- `path/to/helper.R` — supporting code

## Tasks

- [ ] First task to complete
- [ ] Second task to complete
- [ ] Final verification step

## Current State

Brief description of where things stand. Update this as work progresses.

## Discussion

### YYYY-MM-DD: Initial planning

Context, decisions made, and rationale for this plan.

## Next Actions

1. Immediate next step
2. Following step

## Completion

When all tasks are complete, run:

```bash
./plan complete AREA/NAME
```

This will:
1. Prompt for knowledge debrief (lessons learned)
2. Archive the plan with timestamp
3. Update the graph index
