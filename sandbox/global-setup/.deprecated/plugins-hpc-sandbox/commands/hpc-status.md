---
name: hpc-status
description: Quick check of HPC job queue
usage: /hpc-status
---

# /hpc-status

Quick shortcut to check your HPC job queue.

Equivalent to: `ssh hpc3 "squeue -u nwechter"`

**Note**: Requires `hpc3` alias in `~/.ssh/config` with key-based auth.

This is a read-only operation and executes immediately without review.
