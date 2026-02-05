# Exercise 1: First HPC Command

Learn how read-only HPC commands work with Claude Code hooks.

## Objective

Run SLURM status commands and observe how the hook system:
1. Injects quick reference documentation on first use
2. Auto-approves read-only commands
3. Shows results without prompts

## Prerequisites

- Completed setup (see setup-guides/)
- SSH connection working (`ssh hpc3 hostname`)
- Claude Code installed and plugin registered

## Steps

### Step 1: Start Claude Code

Open a terminal and navigate to your project:

```bash
cd ~/your-project  # or any directory
claude
```

### Step 2: Check Job Queue

In Claude, ask:

> "Check my HPC job queue"

**Expected behavior:**
1. First HPC command triggers quick reference injection
2. You see the HPC QUICK REFERENCE message
3. Claude runs `ssh hpc3 "squeue -u $USER"`
4. Results displayed (empty queue or your jobs)

**Observe:**
- No permission prompt appeared
- Command executed immediately
- Hook provided documentation context

### Step 3: Check Job History

Ask:

> "Show my recent job history"

**Expected behavior:**
1. Claude runs `ssh hpc3 "sacct -u $USER --starttime=today"`
2. Recent jobs displayed
3. No permission prompt

### Step 4: Check Cluster Status

Ask:

> "What's the cluster status?"

**Expected behavior:**
1. Claude runs `ssh hpc3 "sinfo -p standard"`
2. Partition info displayed
3. No prompt needed

### Step 5: List Remote Files

Ask:

> "List the contents of /share/crsp/lab on HPC"

**Expected behavior:**
1. Claude runs `ssh hpc3 "ls /share/crsp/lab"`
2. Lab directories listed
3. Read-only, no prompt

## Verification Checklist

- [ ] Quick reference appeared on first HPC command
- [ ] squeue executed without prompting
- [ ] sacct executed without prompting
- [ ] sinfo executed without prompting
- [ ] ls on HPC executed without prompting
- [ ] All results displayed correctly

## What You Learned

1. **Read-only commands are auto-approved** - squeue, sacct, sinfo, ls, cat, tail
2. **First-use documentation injection** - Hook provides context automatically
3. **SSH alias works** - `hpc3` resolves to full hostname
4. **Socket persistence** - Repeated commands don't trigger DUO

## Common Issues

### "Permission denied"
- SSH key not configured correctly
- Run `ssh -v hpc3 hostname` to debug

### Quick reference not appearing
- Hook may not be registered
- Check `~/.claude/settings.json` has hooks configured

### Commands prompting for permission
- Check settings.json allows these patterns
- Verify hook script is executable

## Next Exercise

→ [02_submit_test_job.md](02_submit_test_job.md) - Submit a job using hpc submit
