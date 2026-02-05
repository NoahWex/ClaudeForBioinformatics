---
name: HPC Operations
description: >
  Use this skill when the user asks to run jobs on HPC, submit SLURM jobs,
  check job status, work on their analysis pipeline, debug HPC issues,
  or anything involving UCI HPC3 cluster. Trigger phrases include:
  "run on HPC", "submit job", "sbatch", "check queue", "squeue",
  "HPC3", "cluster", "SLURM", "debug on HPC", "work on pipeline".
version: 3.0.0
---

# HPC Operations for UCI HPC3

You are managing HPC operations for dalawson_lab on UCI HPC3.

## Important: Job Submission Protocol

**Raw `sbatch` is disabled.** You must use the `hpc-submit` wrapper for all job submissions.

The wrapper enforces explicit documentation via required CLI arguments. This ensures every submission is justified and auditable.

---

## Quick Reference

| Item | Value |
|------|-------|
| **Connection** | `ssh hpc3` (requires ~/.ssh/config alias) |
| **User** | nwechter |
| **Account** | dalawson_lab |
| **Partition** | standard (or free for preemptible) |
| **Wrapper** | `~/.claude/plugins/hpc-sandbox/bin/hpc-submit` |

---

## Job Submission with hpc-submit

### Required Arguments

Every submission requires:

| Argument | Description |
|----------|-------------|
| `--script PATH` | SLURM script to submit (or use `--inline`) |
| `--purpose "TEXT"` | Why this job exists and why it's being run now |
| `--outputs "PATH"` | Expected output files/directories |

### Conditional Arguments

These are required based on context:

| Argument | When Required |
|----------|---------------|
| `--rerun-reason "TEXT"` | When output files already exist |
| `--risk-note "TEXT"` | When memory ≥ 64GB or time ≥ 24h (unless `--trivial`) |
| `--overwrites` | Explicit flag to acknowledge overwriting existing outputs |

### Optional Arguments

| Argument | Description |
|----------|-------------|
| `--trivial` | Skip risk-note requirement for quick/test jobs |
| `--depends-on JOB_ID` | Create job dependency (SLURM afterok) |
| `--hypothesis "TEXT"` | For debug jobs: what you expect to learn |
| `--dry-run` | Show what would happen without submitting |
| `--sbatch-arg "ARG"` | Pass additional args to sbatch (repeatable) |
| `--inline` | Read script content from stdin (for heredocs) |

---

## Usage Examples

### Standard Script Submission

```bash
~/.claude/plugins/hpc-sandbox/bin/hpc-submit \
  --script /share/crsp/lab/dalawson/nwechter/project/run/build_milo.sh \
  --purpose "Build canonical Milo object for Gray study" \
  --outputs "/share/crsp/lab/dalawson/nwechter/project/outputs/gray_milo.rds"
```

### Re-run After Failure (Output Exists)

```bash
~/.claude/plugins/hpc-sandbox/bin/hpc-submit \
  --script /path/to/run/build_milo.sh \
  --purpose "Rebuild gray Milo after fixing cell ID parsing" \
  --outputs "/path/to/outputs/gray_milo.rds" \
  --rerun-reason "Previous run failed: barcode format mismatch in cell IDs" \
  --overwrites
```

### Heavy Resources (Memory ≥ 64GB or Time ≥ 24h)

```bash
~/.claude/plugins/hpc-sandbox/bin/hpc-submit \
  --script /path/to/run/integrate_scvi.sh \
  --purpose "Full cohort scVI integration (258K cells)" \
  --outputs "/path/to/outputs/integrated.h5ad" \
  --risk-note "258K cells requires 128GB; tested successfully on 50K subset"
```

### Job Chaining (Dependencies)

```bash
~/.claude/plugins/hpc-sandbox/bin/hpc-submit \
  --script /path/to/run/run_da.sh \
  --purpose "Run DA analysis on gray Milo (depends on build completing)" \
  --outputs "/path/to/outputs/gray_da.rds" \
  --depends-on 12345678
```

### Quick Test (Escape Hatch)

```bash
~/.claude/plugins/hpc-sandbox/bin/hpc-submit \
  --script /path/to/run/test_metrics.sh \
  --purpose "Validate metric calculation on 100-cell test set" \
  --outputs "/dev/null" \
  --trivial
```

### Inline Heredoc

```bash
~/.claude/plugins/hpc-sandbox/bin/hpc-submit --inline \
  --purpose "Quick connectivity test" \
  --outputs "stdout" \
  --trivial << 'EOF'
#!/bin/bash
#SBATCH --job-name=connectivity_test
#SBATCH --time=00:05:00
#SBATCH --mem=1G
#SBATCH --partition=free
#SBATCH --account=dalawson_lab

echo "Hello from HPC"
hostname
date
EOF
```

### Dry Run (Preview Without Submitting)

Add `--dry-run` to any command to see what would happen:

```bash
~/.claude/plugins/hpc-sandbox/bin/hpc-submit \
  --script /path/to/script.sh \
  --purpose "Test submission" \
  --outputs "/path/to/outputs" \
  --dry-run
```

---

## Read-Only Commands (Auto-Approved)

These commands work directly via SSH without the wrapper:

```bash
# Check job queue
ssh hpc3 "squeue -u nwechter"

# Recent job history
ssh hpc3 "sacct -u nwechter --starttime=today -o JobID,JobName,State,Elapsed,MaxRSS"

# Cluster status
ssh hpc3 "sinfo -p standard"

# Storage quota
ssh hpc3 "quota"

# List files
ssh hpc3 "ls -la /path/to/directory"

# View file contents
ssh hpc3 "cat /path/to/file"

# Tail job output
ssh hpc3 "tail -100 /path/to/logs/job_12345.out"
```

---

## Job Management

### Cancel a Job (Prompts for Confirmation)

```bash
ssh hpc3 "scancel JOB_ID"
```

### Check Job Efficiency After Completion

```bash
ssh hpc3 "sacct -j JOB_ID -o JobID,MaxRSS,MaxVMSize,Elapsed,State"
```

---

## Permission Architecture

### SLURM Commands

| Command | Behavior |
|---------|----------|
| `ssh hpc3 "squeue ..."` | Auto-approved (read-only) |
| `ssh hpc3 "sacct ..."` | Auto-approved (read-only) |
| `ssh hpc3 "sbatch ..."` | **BLOCKED** - use `hpc-submit` |
| `hpc-submit ...` | Validates arguments, then submits |
| `ssh hpc3 "scancel ..."` | Prompts for confirmation |
| `ssh hpc3 "srun ..."` | **Discouraged** - use batch jobs instead |

**If a command is blocked:** Do NOT attempt workarounds. Explain what you need and ask the user how to proceed.

**Interactive sessions (`srun`):** Avoid using interactive sessions. Batch jobs via `hpc-submit` are preferred because they are documented, auditable, and reproducible. If debugging requires interactivity, ask the user to run the session manually.

### Git Commands (via SSH)

Git operations on CRSP-hosted repos are managed with the following tiers:

| Category | Commands | Behavior |
|----------|----------|----------|
| **Safe** | `status`, `log`, `diff`, `show`, `add`, `stash`, `tag`, `fetch` | Auto-approved |
| **Supervised** | `commit`, `checkout`, `switch`, `branch`, `merge`, `rebase`, `push`, `pull`, `reset` | Prompts for confirmation |
| **Dangerous** | `push --force`, `reset --hard`, `clean -f`, `branch -D` | Prompts with strong warning |

**Examples:**
```bash
# Auto-approved (read-only and staging)
ssh hpc3 "cd /path/to/repo && git status"
ssh hpc3 "cd /path/to/repo && git add ."

# Prompts for confirmation (commits and branch ops)
ssh hpc3 "cd /path/to/repo && git commit -m 'message'"
ssh hpc3 "cd /path/to/repo && git push origin main"
ssh hpc3 "cd /path/to/repo && git checkout feature-branch"

# Prompts with strong warning (dangerous)
ssh hpc3 "cd /path/to/repo && git push --force"  # Avoid!
```

---

## SLURM Script Template

```bash
#!/bin/bash
#SBATCH --job-name=my_job
#SBATCH --partition=standard
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/logs/%x_%j.out
#SBATCH --error=/share/crsp/lab/dalawson/nwechter/logs/%x_%j.err

set -euo pipefail

echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: ${SLURMD_NODENAME}"
echo "Start: $(date)"

module load singularity

singularity exec \
    --bind /share/crsp/lab/dalawson/share:/share/crsp/lab/dalawson/share:ro \
    --bind /share/crsp/lab/dalawson/nwechter:/share/crsp/lab/dalawson/nwechter:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /dfs8:/dfs8:ro \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    /path/to/container.sif \
    python /path/to/script.py

echo "Done: $(date)"
```

---

## Common Gotchas

### Always Load Singularity First
```bash
module load singularity  # Before any singularity exec
```

### Don't Bind /scratch
Not all nodes have `/scratch`. Only bind CRSP paths.

### Use Absolute Paths for Logs
```bash
#SBATCH --output=/full/absolute/path/logs/%x_%j.out  # Correct
#SBATCH --output=logs/%x_%j.out                       # Wrong - goes to sbatch cwd
```

### CRSP Stale File Handle
Write large files to `/tmp` first, then copy to CRSP.

### R Packages in Containers
Container `/home/jovyan` is read-only. Bind user R library:
```bash
--bind "/data/homezvol0/nwechter/.../R/library:/home/jovyan/R/library:ro"
--env "R_LIBS_USER=/home/jovyan/R/library"
```

---

## Documentation Reference

**Proactively consult UCI RCIC documentation** when:
- Encountering unfamiliar SLURM options or errors
- Debugging module/software issues
- Optimizing resource requests
- Working with containers or storage

| Topic | URL |
|-------|-----|
| SLURM | https://rcic.uci.edu/slurm/ |
| Storage | https://rcic.uci.edu/storage/ |
| Software/Modules | https://rcic.uci.edu/software/ |
| Containers | https://rcic.uci.edu/singularity/ |

**Use WebFetch proactively** to look up documentation before guessing at solutions. RCIC docs are authoritative for UCI HPC3.
