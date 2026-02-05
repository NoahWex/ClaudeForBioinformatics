# Exercise 2: Submit Test Job

Learn how the hpc-toolkit enforces documented job submissions.

## Objective

Understand how:
1. Raw sbatch is blocked with redirect message
2. `hpc submit` requires explicit justification
3. --dry-run allows safe testing
4. Submissions are audited

## Prerequisites

- Completed Exercise 1
- Toolkit installed with hooks registered

## Steps

### Step 1: Try Raw sbatch (See It Blocked)

In Claude, ask:

> "Submit a test job with sbatch"

**Expected behavior:**
1. Claude attempts `ssh hpc3 "sbatch test.sh"`
2. Gate hook BLOCKS the command
3. Message redirects to `hpc submit`
4. No job submitted

**Observe the redirect message:**
```
BLOCKED: Use hpc submit instead of raw sbatch.
Suggested: ~/.claude/hpc-toolkit/bin/hpc submit /path/to/script.sh --purpose "..." --outputs "..."
```

### Step 2: Use hpc submit with Dry Run

Ask Claude:

> "Use hpc submit to submit a test job with dry-run"

Or run directly:

```bash
~/.claude/hpc-toolkit/bin/hpc submit /path/to/your/test.sh \
  --purpose "Exercise 2 test job" \
  --outputs "stdout" \
  --dry-run
```

**Expected behavior:**
1. `[DRY RUN]` message shows what would execute
2. No actual job submitted

**Note:** You need a script file on HPC. If you don't have one, create one first:

```bash
~/.claude/hpc-toolkit/bin/hpc file write /pub/$USER/exercise2_test.sh --content '#!/bin/bash
#SBATCH --job-name=exercise2_test
#SBATCH --time=00:01:00
#SBATCH --mem=1G
#SBATCH --partition=free
#SBATCH --account=YOUR_ACCOUNT

echo "Hello from Exercise 2!"
hostname
date'
```

Then:

```bash
~/.claude/hpc-toolkit/bin/hpc submit /pub/$USER/exercise2_test.sh \
  --purpose "Exercise 2 test job" \
  --outputs "stdout" \
  --dry-run
```

### Step 3: Check Status Commands

```bash
~/.claude/hpc-toolkit/bin/hpc status
```

**Expected:** Your job queue (or empty).

### Step 4: Actual Submission (Optional)

Remove `--dry-run` to submit for real:

```bash
~/.claude/hpc-toolkit/bin/hpc submit /pub/$USER/exercise2_test.sh \
  --purpose "Exercise 2 verification" \
  --outputs "stdout"
```

**Expected behavior:**
1. Job submitted
2. Job ID returned

**Monitor and cleanup:**
```bash
~/.claude/hpc-toolkit/bin/hpc status
~/.claude/hpc-toolkit/bin/hpc cancel JOB_ID --reason "Exercise complete"
```

### Step 5: Check Audit Log

View what was recorded:

```bash
cat ~/.claude/hpc-toolkit/logs/audit.jsonl | tail -5
```

**Expected:** JSON entries for your submissions with timestamps, purpose, and job IDs.

## Verification Checklist

- [ ] Raw sbatch was blocked with redirect
- [ ] `hpc submit --dry-run` showed what would execute without submitting
- [ ] `hpc status` returned job queue
- [ ] (Optional) Real submission worked
- [ ] Audit log contains submission record

## What You Learned

1. **Raw sbatch is blocked** — must use `hpc submit`
2. **Required fields** — purpose and outputs must be specified
3. **--dry-run is safe** — preview without submitting
4. **Audit trail** — all submissions logged to `audit.jsonl`

## Next Exercise

→ [03_git_via_ssh.md](03_git_via_ssh.md) - Git operations on CRSP repos
