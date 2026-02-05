---
name: hpc
description: Use this skill when the user asks to run jobs on HPC, submit SLURM jobs, check job status, work on their analysis pipeline, debug HPC issues, manage git workspaces, or anything involving UCI HPC3 cluster. Trigger phrases include "run on HPC", "submit job", "sbatch", "check queue", "squeue", "HPC3", "cluster", "SLURM", "workspace", "git mirror", "DFS mirror".
---

# HPC Operations

Use `~/.claude/hpc-toolkit/bin/hpc` for all HPC operations. All paths must be HPC paths (starting with `/`).

## CLI Reference (hpc --help)

```
hpc - Unified CLI for UCI HPC3 operations

USAGE:
  hpc <command> [args]

COMMANDS:
  submit      Submit a batch job
  shell       Start interactive session
  status      Check jobs (squeue)
  logs        View job output
  cancel      Cancel a job
  file        File operations (ls, cat, rm, cp, mv, write)
  git         Git operations (login node, no srun)
  workspace   Manage fast git mirrors on DFS

EXAMPLES:
  hpc submit /path/script.sh --purpose "Run analysis" --outputs "/path/results/"
  hpc shell --purpose "Debug issue" --time 1:00:00
  hpc status
  hpc logs 12345678
  hpc cancel 12345678 --reason "Wrong parameters"
  hpc file ls /share/crsp/lab/YOUR_LAB/YOUR_USER/
  hpc git status
  hpc workspace init /share/crsp/lab/YOUR_LAB/YOUR_USER/YourProject

Run 'hpc <command> --help' for command-specific help.
```

## Submit Command (hpc submit --help)

```
hpc submit - Submit a batch job

USAGE:
  hpc submit <script_path> --purpose "..." --outputs "..." [options]

REQUIRED:
  <script_path>     Path to script ON HPC (must start with /)
  --purpose "..."   Why running this job
  --outputs "..."   Where outputs go

OPTIONS:
  --sbatch "..."    Pass-through sbatch args (--array, --mem, etc.)
  --dry-run         Show what would happen without submitting

EXAMPLE:
  hpc submit /share/crsp/lab/YOUR_LAB/YOUR_USER/run.sh \
    --purpose "Process batch 1" \
    --outputs "/share/crsp/lab/YOUR_LAB/YOUR_USER/results/" \
    --sbatch "--array=1-10 --mem=64G --time=08:00:00"
```

## Commands Summary

| Task | Command |
|------|---------|
| Submit job | `hpc submit /path/script.sh --purpose "..." --outputs "..." [--sbatch "..."]` |
| Run command | `hpc shell --cmd "..." --purpose "..." --time 1:00:00 [--module MOD]` |
| Status | `hpc status [job_id]` |
| Logs | `hpc logs <job_id>` |
| Cancel | `hpc cancel <job_id>` |
| Files | `hpc file ls/cat/write/rm/cp/mv /path` |
| Git | `hpc git <args>` (**all git writes must use this**) |
| Workspace | `hpc workspace <init\|list\|sync\|status\|destroy>` |

## Quick Reference

**Partitions:**
| Partition | Time | Notes |
|-----------|------|-------|
| `free` | 3 days | No charge, preemptible |
| `standard` | 14 days | Default |
| `free-gpu` | 3 days | GPU, preemptible |
| `gpu` | 14 days | GPU, needs `_gpu` account |

**Modules:** Use `ssh SSH_ALIAS 'module avail <name>'` to find versions. (Replace SSH_ALIAS with your config value.)

**Paths:** (set in `~/.claude/hpc-toolkit/config.sh`)
- Lab: `$CRSP_HPC_PREFIX/YOUR_USER/`
- Scratch: `/pub/$USER/` or cluster-specific

## Submit Examples

Basic submission:
```bash
hpc submit /path/to/script.sh --purpose "Run analysis" --outputs "/path/to/outputs/"
```

With extra sbatch arguments (memory, time, array jobs):
```bash
hpc submit /path/to/script.sh \
  --purpose "High-memory analysis" \
  --outputs "/path/to/outputs/" \
  --sbatch "--mem=180G --time=08:00:00"
```

Array job submission:
```bash
hpc submit /path/to/script.sh \
  --purpose "Process samples 1,3,8" \
  --outputs "/path/to/outputs/" \
  --sbatch "--array=1,3,8 --mem=32G --time=04:00:00"
```

Dry run (preview without submitting):
```bash
hpc submit /path/to/script.sh --purpose "..." --outputs "..." --dry-run
```

## SLURM Template

```bash
#!/bin/bash
#SBATCH --job-name=name
#SBATCH --account=YOUR_ACCOUNT
#SBATCH --partition=standard
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

module load R/4.2.2
Rscript analysis.R
```

## Git Operations

All git **write** operations must go through `hpc git` — CRSP's macOS mount caches `.git/objects/` aggressively, causing desync between local and HPC views.

| Task | Command |
|------|---------|
| Status | `hpc git status` or local `git status` |
| Stage files | `hpc git add <paths>` |
| Commit | `hpc git commit -m "message"` |
| Log | `hpc git log --oneline -5` or local |
| Diff | `hpc git diff --stat` or local |

Read-only git commands work locally. Write commands are blocked by a PreToolUse hook.

### Git Workspaces

Fast git mirrors on DFS. After `hpc workspace init <path>`, `hpc git` auto-uses the mirror.

| Task | Command |
|------|---------|
| Init | `hpc workspace init /share/crsp/.../RepoName` |
| List | `hpc workspace list` |
| Status | `hpc workspace status <name>` |
| Sync | `hpc workspace sync <name>` |
| Destroy | `hpc workspace destroy <name>` |

## Allowed Direct SSH (read-only)

```bash
ssh SSH_ALIAS 'squeue -u $USER'
ssh SSH_ALIAS 'sacct -j <id>'
ssh SSH_ALIAS 'sinfo'
ssh SSH_ALIAS 'module avail <name>'
```

## UCI RCIC Documentation

| Topic | URL |
|-------|-----|
| SLURM | https://rcic.uci.edu/slurm/slurm.html |
| Examples | https://rcic.uci.edu/slurm/examples.html |
| Modules | https://rcic.uci.edu/software/modules.html |
| Software | https://rcic.uci.edu/software/software.html |
