# 04: Verify Setup

Run through these tests to confirm everything works.

## 1. SSH Connection

```bash
ssh hpc3 hostname
```

**Expected:** `hpc3-xx-xx.rcic.uci.edu` (or your cluster's hostname)

**If fails:** Check `~/.ssh/config`, VPN status. See [troubleshooting.md](troubleshooting.md).

## 2. SLURM Access

```bash
ssh hpc3 "squeue -u \$USER"
```

**Expected:** Empty queue or your jobs listed.

## 3. CRSP Access

```bash
ssh hpc3 "ls /share/crsp/lab"
```

**Expected:** List of lab directories.

## 4. Claude Code

```bash
claude --version
```

**Expected:** Version number displayed.

## 5. Toolkit Installed

```bash
~/.claude/hpc-toolkit/bin/hpc --help
```

**Expected:** Usage information with commands: submit, shell, status, logs, cancel, file, git, workspace.

## 6. Config Loaded

```bash
source ~/.claude/hpc-toolkit/config.sh && echo "SSH_ALIAS=$SSH_ALIAS CRSP_HPC_PREFIX=$CRSP_HPC_PREFIX"
```

**Expected:** Your configured values printed.

## 7. Hook Test — Read-Only Commands

Start Claude Code and test:

```bash
claude
```

Ask: "Check my HPC job queue with squeue"

**Expected:**
1. First command injects quick reference docs (CLAUDE_GUIDE.md)
2. Command executes without being blocked
3. Queue results displayed

## 8. Hook Test — Blocked sbatch

In the same Claude session, ask: "Submit a test job with sbatch test.sh"

**Expected:**
1. Command **BLOCKED** by gate.sh
2. Message suggests using `hpc submit` instead
3. No job submitted

## 9. Hook Test — Git Write Block

In a Claude session inside a CRSP-mounted project, ask: "Run git commit -m test"

**Expected:**
1. Command **BLOCKED** by gate.sh
2. Message suggests using `hpc git commit -m test` instead

## 10. Toolkit Dry Run

```bash
~/.claude/hpc-toolkit/bin/hpc status
```

**Expected:** Your job queue (or empty).

```bash
~/.claude/hpc-toolkit/bin/hpc submit /path/to/test.sh \
  --purpose "Verify setup" \
  --outputs "/tmp/" \
  --dry-run
```

**Expected:** `[DRY RUN] Would execute: ssh ... sbatch ...`

## Complete Checklist

- [ ] SSH connection works
- [ ] SLURM commands work
- [ ] CRSP accessible
- [ ] Claude Code runs
- [ ] Toolkit installed (`hpc --help` works)
- [ ] Config loaded (SSH_ALIAS, CRSP paths)
- [ ] Read-only commands auto-allowed
- [ ] Raw sbatch blocked
- [ ] Local git writes blocked on CRSP mount
- [ ] `hpc status` works
- [ ] `hpc submit --dry-run` works

## Setup Complete

You're ready to start the exercises:

→ [../exercises/01_first_hpc_command.md](../exercises/01_first_hpc_command.md)

## If Something Fails

See [troubleshooting.md](troubleshooting.md) for common issues.
