# SLURM Job Practices

Best practices for HPC job scripts, derived from dalawson_lab pipeline patterns.

## Script Organization

**Two-script pattern:**
- `submit_*.sh` - Validates prerequisites, submits jobs, handles dependencies
- `run_*.sh` - Worker script with SBATCH directives, executes in container

## Worker Script Structure (run_*.sh)

### Header
```bash
#!/bin/bash
#SBATCH --job-name=experiment_stepN
#SBATCH --partition=standard      # or free-gpu
#SBATCH --account=dalawson_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=/full/absolute/path/to/logs/step_N_%j.out
#SBATCH --error=/full/absolute/path/to/logs/step_N_%j.err

set -euo pipefail
```

### Initialization
```bash
PROJECT_ROOT="/share/crsp/lab/dalawson/nwechter/ProjectName"
cd "$PROJECT_ROOT"
module load singularity/3.11.3
```

### Banner Output
```bash
echo "============================================================================"
echo "Step N: Description"
echo "============================================================================"
echo "Job ID:           ${SLURM_JOB_ID}"
echo "Node:             ${SLURMD_NODENAME}"
echo "Memory:           ${SLURM_MEM_PER_NODE:-N/A}"
echo "CPUs:             ${SLURM_CPUS_PER_TASK:-N/A}"
echo "Working Dir:      $(pwd)"
echo "Start Time:       $(date)"
echo "============================================================================"
```

### Container Execution
```bash
singularity exec \
    --bind /share/crsp/lab/dalawson/share:/share/crsp/lab/dalawson/share:ro \
    --bind /share/crsp/lab/dalawson/nwechter:/share/crsp/lab/dalawson/nwechter:rw \
    --bind /dfs7:/dfs7:ro \
    --bind /dfs8:/dfs8:ro \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    /path/to/container.sif \
    Rscript /path/to/script.R
```

### Completion Banner
```bash
echo "============================================================================"
echo "Completed: $(date)"
echo "============================================================================"
```

## Log Directory Convention

| Directory | Use |
|-----------|-----|
| `logs/` | Production pipeline - canonical scripts that produce tracked outputs |
| `.dev/logs/` | Development/experimental - debugging, testing, exploration |

Production logs are part of the project record. Dev logs are gitignored ephemeral artifacts.

## Job Submission Pattern

When presenting SLURM job submissions to the user:

```bash
# Submit
sbatch /full/path/to/run/script.sh

# Submit and tail output (glob waits for file to exist)
JOB=$(sbatch --parsable /full/path/to/run/script.sh) && \
  while [ ! -f /full/path/to/logs/*${JOB}* ]; do sleep 1; done && \
  tail -f /full/path/to/logs/*${JOB}*

# Check status
squeue -u $USER
```

**Preferences:**
- Always use absolute paths to .sh scripts
- Use `--parsable` to capture job ID
- Glob `*${JOB}*` catches both .out and .err
- Show exact sbatch command first
- Include how to check job status
