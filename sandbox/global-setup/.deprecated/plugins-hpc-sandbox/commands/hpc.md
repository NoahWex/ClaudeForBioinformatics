---
name: hpc
description: Execute HPC operations on UCI HPC3
usage: /hpc [command]
examples:
  - /hpc status
  - /hpc submit job.sh
  - /hpc jobs
  - /hpc cancel 12345
---

# /hpc Command

Interact with UCI HPC3 cluster.

**Requires**: `hpc3` alias in `~/.ssh/config` with key-based authentication.

## Usage

`/hpc [subcommand] [args]`

## Subcommands

| Command | Description |
|---------|-------------|
| `status` | Check your job queue (squeue) |
| `jobs` | Show recent job history (sacct) |
| `info` | Show cluster partition info |
| `submit <script>` | Submit a job script |
| `cancel <job_id>` | Cancel a running job |
| `logs <job_id>` | View output logs for a job |

## Examples

```
/hpc status
→ Shows your current jobs in queue

/hpc jobs
→ Shows recent job history with resource usage

/hpc submit /path/to/script.sh
→ Submits job script to SLURM

/hpc cancel 12345678
→ Cancels the specified job

/hpc logs 12345678
→ Shows output from job
```

## Permission Model

- **Read-only commands** (status, jobs, info, logs): Can be auto-approved via settings
- **Write commands** (submit, cancel): Prompt for approval by default

**If a command is blocked or denied**: Do not attempt workarounds. Ask the user how to proceed.
