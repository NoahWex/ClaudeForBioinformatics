# HPC Toolkit Guide

Use `~/.claude/hpc-toolkit/bin/hpc` for HPC operations. All paths must be HPC paths (starting with `/`).

## Commands

| Task | Command |
|------|---------|
| Submit job | `~/.claude/hpc-toolkit/bin/hpc submit /path/script.sh --purpose "..." --outputs "..." [--sbatch "extra args"]` |
| Run command | `~/.claude/hpc-toolkit/bin/hpc shell --cmd "..." --purpose "..." --time 1:00:00 [--module MOD]` |
| Status | `~/.claude/hpc-toolkit/bin/hpc status [job_id]` |
| Logs | `~/.claude/hpc-toolkit/bin/hpc logs <job_id>` |
| Cancel | `~/.claude/hpc-toolkit/bin/hpc cancel <job_id>` |
| Files | `~/.claude/hpc-toolkit/bin/hpc file ls/cat/write/rm/cp/mv /path` |
| Git | `~/.claude/hpc-toolkit/bin/hpc git <args>` (login node, **writes must use this**) |
| Workspace | `~/.claude/hpc-toolkit/bin/hpc workspace <init\|list\|sync\|status\|destroy>` |

## Git Operations

All git writes go through the login node to avoid CRSP cache issues.

| Task | Command |
|------|---------|
| Status | `hpc git status` |
| Stage | `hpc git add <paths>` |
| Commit | `hpc git commit -m "message"` |
| Log | `hpc git log --oneline -5` |
| Diff | `hpc git diff --stat` |

Read-only commands (status, log, diff) also work locally.
Write commands are blocked locally by hook.

## Quick Reference

**Partitions:**
| Partition | Time | Notes |
|-----------|------|-------|
| `free` | 3 days | No charge, preemptible |
| `standard` | 14 days | Default |
| `free-gpu` | 3 days | GPU, preemptible |
| `gpu` | 14 days | GPU, needs `_gpu` account |

**Modules:** Use `module avail <name>` to find versions.

**Paths:** (configure in `~/.claude/hpc-toolkit/config.sh`)
- Lab: `$CRSP_HPC_PREFIX/YOUR_USER/`
- Scratch: `/pub/$USER/` or `/dfs7/YOUR_ACCOUNT/`

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

All git **write** operations must go through `hpc git` due to CRSP filesystem caching.
Local read-only git commands (status, log, diff) work fine.

| Task | Command |
|------|---------|
| Status | `hpc git status` or local `git status` |
| Stage files | `hpc git add <paths>` |
| Commit | `hpc git commit -m "message"` |
| Log | `hpc git log --oneline -5` or local |
| Diff | `hpc git diff --stat` or local |
| Reset | `hpc git reset <args>` |

**Why**: CRSP's macOS SMB mount caches `.git/objects/` directory listings. Writes from the local mount create objects invisible to HPC (and vice versa), causing `fatal: bad object HEAD`, orphaned commits, and index.lock conflicts.

## Git Workspaces (DFS Acceleration)

Workspaces keep a bare git mirror on fast DFS (`/pub/$USER/git-mirrors/`) while the working tree stays on CRSP. `hpc git` auto-detects workspaces — no workflow change needed.

| Task | Command |
|------|---------|
| Create workspace | `hpc workspace init /share/crsp/lab/YOUR_LAB/YOUR_USER/YOUR_PROJECT` |
| List workspaces | `hpc workspace list` |
| Check sync state | `hpc workspace status YOUR_PROJECT` |
| Manual sync | `hpc workspace sync YOUR_PROJECT` |
| Remove workspace | `hpc workspace destroy YOUR_PROJECT` |

**How it works**: After `workspace init`, `hpc git` sets `GIT_DIR` to the DFS mirror and `GIT_WORK_TREE` to the CRSP path. Write operations (commit, add, merge, etc.) auto-sync back to CRSP `.git/` so local macOS reads stay current.

**Repo detection**: `hpc git` resolves the repo by: (1) `--repo` flag, (2) macOS CWD mapped to HPC path, (3) hardcoded default.

## Allowed Direct SSH

```bash
ssh SSH_ALIAS 'squeue -u $USER'
ssh SSH_ALIAS 'sacct -j <id>'
ssh SSH_ALIAS 'sinfo'
ssh SSH_ALIAS 'module avail <name>'
```

## UCI RCIC Docs

| Topic | URL |
|-------|-----|
| SLURM | https://rcic.uci.edu/slurm/slurm.html |
| Examples | https://rcic.uci.edu/slurm/examples.html |
| Modules | https://rcic.uci.edu/software/modules.html |
| Software | https://rcic.uci.edu/software/software.html |
