# Exercise 2: Writing SLURM Jobs

**Objective**: Create and submit a SLURM batch script to run your QC pipeline on HPC.

**Time**: ~35 minutes

---

## Background

SLURM (Simple Linux Utility for Resource Management) is the job scheduler on UCI HPC3. To run your analysis:

1. Write a batch script with resource requests
2. Submit with `sbatch`
3. Monitor with `squeue`
4. Check results when complete

## Your Task

### Part 1: Basic SLURM script (15 min)

Ask Claude to create a SLURM job script:

```
Please create working/jobs/run_qc_gray.sh, a SLURM batch script that:

1. Requests resources:
   - Account: dalawson_lab
   - Partition: standard
   - Time: 2 hours
   - Memory: 32GB
   - CPUs: 4
   - Job name: qc_gray

2. Sets up the environment:
   - Loads singularity module
   - Sets working directory

3. Runs the QC script in a container:
   - Container path from ../shared/container_info.md
   - Runs: Rscript scripts/05_qc_normalize.R --study gray --manifest manifest.yaml

4. Includes:
   - Email notification on completion
   - Output and error log files
```

**Review**: Check the `#SBATCH` directives are correct.

### Part 2: Understand the components (10 min)

Ask Claude to explain the script:

```
Please explain each part of the SLURM script:
1. What does each #SBATCH directive do?
2. Why do we load singularity?
3. What does the singularity exec command do?
4. Where will output logs be written?
```

### Part 3: Submit and monitor (10 min)

Using the HPC skill or SSH:

```
Please help me:
1. Transfer the job script to HPC3 (or use the shared CRSP path)
2. Submit the job with sbatch
3. Check the job status with squeue
4. Show me how to check the output logs
```

**Note**: If using Claude's `/hpc` skill, it can handle submission directly.

---

## SLURM Script Template

See `templates/single_job.sh` for a ready-to-use template:

```bash
#!/bin/bash
#SBATCH --job-name=qc_gray
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=logs/qc_gray_%j.out
#SBATCH --error=logs/qc_gray_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=your_email@uci.edu

# Load modules
module load singularity

# Set paths
CONTAINER="/share/crsp/lab/dalawson/share/0_Resources/containers/r_seurat_4.4.2.sif"
WORKDIR="/share/crsp/lab/dalawson/share/0_Resources/ClaudeForBioinformatics_NW/Workshop/working"

# Change to working directory
cd $WORKDIR

# Create logs directory if needed
mkdir -p logs

# Run analysis
echo "Starting QC pipeline for gray..."
echo "Time: $(date)"

singularity exec $CONTAINER \
    Rscript scripts/05_qc_normalize.R \
    --study gray \
    --manifest manifest.yaml

echo "Completed at: $(date)"
```

## Key SLURM Commands

```bash
# Submit job
sbatch jobs/run_qc_gray.sh

# Check job status
squeue -u $USER

# Cancel job
scancel <jobid>

# View job details
scontrol show job <jobid>

# Check completed job info
sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed,MaxRSS

# View output logs
cat logs/qc_gray_<jobid>.out
cat logs/qc_gray_<jobid>.err
```

## Understanding Resource Requests

| Directive | Purpose | How to Choose |
|-----------|---------|---------------|
| `--time` | Max runtime | Estimate + 50% buffer |
| `--mem` | RAM | Check data size, add buffer |
| `--cpus-per-task` | CPU cores | More for parallelizable tasks |
| `--partition` | Queue type | `standard` for most, `free` for testing |

## Discussion Questions

1. What happens if your job exceeds the requested time?
2. How would you request more memory if the job fails?
3. Why use a container instead of module-loaded R?

## Troubleshooting

### Job fails immediately
Check the error log:
```bash
cat logs/qc_gray_<jobid>.err
```

### Out of memory (OOM)
```bash
sacct -j <jobid> --format=JobID,MaxRSS,ReqMem
```
Increase `--mem` if MaxRSS is near ReqMem.

### Module not found
```bash
module avail singularity
```

### Container not found
Verify the path is accessible from compute nodes:
```bash
ls -la /share/crsp/lab/dalawson/share/0_Resources/containers/
```

## Checkpoint Verification

Before moving on, verify:
- [ ] SLURM script created with correct directives
- [ ] Job submitted successfully (got job ID)
- [ ] Can monitor job with squeue
- [ ] Know where to find output logs
- [ ] Understand basic troubleshooting steps

## Next Steps

In [Exercise 3](03_array_jobs.md), you'll learn to process multiple samples in parallel using array jobs.
