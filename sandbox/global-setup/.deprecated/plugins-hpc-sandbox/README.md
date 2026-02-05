# HPC Sandbox Plugin

A Claude Code plugin for safe, documented HPC job management on UCI HPC3.

## Overview

This plugin provides:
- **Form-based job submission** via `hpc-submit` wrapper
- **Permission hooks** that block raw sbatch and prompt for sensitive operations
- **Audit logging** for all submissions
- **Safety guidelines** loaded into Claude's context

## Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    Claude attempts HPC command               │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│              PreToolUse Hook (validate-hpc-command.sh)       │
│                                                              │
│  • Read-only (squeue, sacct) → Pass through                 │
│  • Raw sbatch → BLOCK, redirect to hpc-submit               │
│  • hpc-submit → Allow (wrapper handles validation)          │
│  • scancel, srun → Prompt user for confirmation             │
│  • Dangerous patterns → Prompt with warning                 │
└─────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────┐
│                    hpc-submit wrapper                        │
│                                                              │
│  Required:    --script, --purpose, --outputs                │
│  Conditional: --rerun-reason, --risk-note, --overwrites     │
│  Optional:    --trivial, --depends-on, --inline, --dry-run  │
│                                                              │
│  Validates → Displays review → Submits → Logs to audit.jsonl│
└─────────────────────────────────────────────────────────────┘
```

## Components

### bin/hpc-submit
Form-based wrapper that requires explicit documentation before job submission.

```bash
# Standard submission
hpc-submit \
  --script /path/to/script.sh \
  --purpose "Build Milo object for Gray study" \
  --outputs "/path/to/outputs/gray_milo.rds"

# Inline heredoc
hpc-submit --inline \
  --purpose "Quick test" \
  --outputs "stdout" \
  --trivial << 'EOF'
#!/bin/bash
#SBATCH --time=00:05:00
echo "Hello"
EOF
```

### hooks/validate-hpc-command.sh
PreToolUse hook that intercepts Bash commands and:
- Allows read-only HPC commands (squeue, sacct, sinfo)
- Blocks raw sbatch with redirect to hpc-submit
- Prompts for confirmation on scancel and srun
- Warns on dangerous patterns (rm -rf, etc.)

### skills/hpc-operations/SKILL.md
Comprehensive documentation loaded when HPC skill is triggered.

### commands/hpc.md, hpc-status.md
Slash command definitions for `/hpc` and `/hpc status`.

## Required Fields by Scenario

| Scenario | Required Fields |
|----------|-----------------|
| New job | `--script`, `--purpose`, `--outputs` |
| Output exists | + `--rerun-reason`, `--overwrites` |
| Heavy resources (≥64GB or ≥24h) | + `--risk-note` (unless `--trivial`) |
| Chained job | + `--depends-on JOB_ID` |
| Inline content | `--inline` instead of `--script` |

## Audit Log

All submissions are logged to `~/.hpc-submit/audit.jsonl`:

```json
{"timestamp":"2026-01-21T10:30:00Z","job_id":"12345678","script":"/path/to/script.sh","purpose":"Build Milo","outputs":"/path/to/output.rds","resources":{"mem":"64G","cpus":"4","time":"2:00:00"},"trivial":false,"inline":false}
```

## Installation

The plugin is installed at `~/.claude/plugins/hpc-sandbox/`.

Ensure the hook is registered in Claude Code settings:
```json
{
  "hooks": {
    "PreToolUse": [
      {
        "matcher": "Bash",
        "command": "~/.claude/plugins/hpc-sandbox/hooks/validate-hpc-command.sh"
      }
    ]
  }
}
```

## Safety Philosophy

The goal is to make Claude **think before submitting** by requiring explicit justification:

1. **Purpose** - Why is this job being run?
2. **Outputs** - What will it produce?
3. **Risk awareness** - What could go wrong? (for heavy resources)
4. **Rerun justification** - Why run again? (if outputs exist)

Arguments become auditable documentation. Shell history shows submissions with full rationale.
