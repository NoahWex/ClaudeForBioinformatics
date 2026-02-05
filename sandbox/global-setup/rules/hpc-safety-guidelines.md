# HPC Safety Guidelines

These guidelines apply when working with HPC cluster operations.

## Core Principle

HPC jobs consume shared resources and may run for hours/days. Mistakes waste compute allocation, block other users, and can corrupt data. **Always verify before acting.**

## Before Job Submission (sbatch)

1. **Verify the script exists and is readable**
2. **Confirm resource requests are reasonable** - flag unusually large requests (>24h, >100GB, >32 cores)
3. **Validate input files exist** - check that data paths resolve
4. **Check for obvious errors** - missing shebang, hardcoded paths, missing module loads
5. **Verify output directories exist**

Present to user:
```
Submitting: /path/to/script.sh
Resources: [CPUs], [memory], [time]
Account: your_account
```

## Before Job Cancellation (scancel)

1. **Show job details first** - use `squeue -j <jobid>`
2. **Confirm ownership** - never cancel jobs belonging to other users
3. **Check job state**:
   - PENDING: Safe to cancel
   - RUNNING: Work will be lost - confirm with user
   - COMPLETING: May corrupt output - warn strongly

## Before Destructive Commands

Commands requiring extra caution:
- `rm -rf` - Show exact paths, count files affected
- `scancel --user` - Kills ALL jobs, list first
- `mv` to overwrite - Check if destination exists
- Writing to shared paths - Verify ownership

**Protocol:**
1. State the action clearly
2. Show what will be affected
3. Explain irreversibility if applicable
4. Let user execute destructive commands themselves when appropriate

## Severity Levels

| Level | Examples | Action |
|-------|----------|--------|
| **LOW** | Quick test job, status checks | Proceed with brief summary |
| **MEDIUM** | Long-running job, overwriting cached results | Show impact, get confirmation |
| **HIGH** | Cancelling running job, deleting processed outputs | Explicit warning, require confirmation |
| **CRITICAL** | `rm -rf` on data, `scancel --user` | Provide command only, user executes |

## Balance: Safety vs Usability

**Do:**
- Batch safe operations (multiple squeue checks)
- Provide copy-paste commands for verification
- Trust explicit user requests

**Don't:**
- Block routine workflows (squeue, sacct are always safe)
- Over-warn on reversible actions
- Execute `rm -rf` directly - provide command for user
