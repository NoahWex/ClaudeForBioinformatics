# Exercise 3: Array Jobs for Parallel Processing

**Objective**: Use SLURM array jobs to process multiple datasets in parallel.

**Time**: ~30 minutes

---

## Background

Array jobs let you run the same script multiple times with different inputs. Instead of writing separate jobs for each study, you write one script that uses `$SLURM_ARRAY_TASK_ID` to determine which input to process.

Benefits:
- Single script manages multiple jobs
- Jobs run in parallel (when resources available)
- Easier to manage than many individual jobs

## Your Task

### Part 1: Create array job script (15 min)

Ask Claude to create a manifest-driven array job:

```
Please create working/jobs/run_qc_array.sh that:

1. Uses SLURM array to process both "gray" and "pal" studies
2. Maps array indices to study names:
   - Index 0 -> gray
   - Index 1 -> pal

3. Uses the same resources as before but with array directive:
   #SBATCH --array=0-1

4. Logs each array task to separate files:
   - logs/qc_${STUDY}_%A_%a.out
   - %A = job ID, %a = array index

5. Uses the study name to call the QC script

Include comments explaining how the array indexing works.
```

### Part 2: Understand array variables (5 min)

Ask Claude to explain:

```
Please explain:
1. What is SLURM_ARRAY_TASK_ID?
2. What's the difference between %A and %a in log filenames?
3. How would you add more studies to the array?
4. What if I want to process a subset of studies?
```

### Part 3: Submit and monitor array job (10 min)

```
Please help me:
1. Submit the array job
2. Show how to monitor all array tasks
3. Explain how to cancel specific array tasks vs the whole job
4. Check logs for a specific task
```

---

## Array Job Template

See `templates/array_job.sh`:

```bash
#!/bin/bash
#SBATCH --job-name=qc_array
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --array=0-1
#SBATCH --output=logs/qc_%A_%a.out
#SBATCH --error=logs/qc_%A_%a.err

# Load modules
module load singularity

# Define studies - ORDER MATTERS (matches array indices)
STUDIES=(gray pal)

# Get study for this task
STUDY=${STUDIES[$SLURM_ARRAY_TASK_ID]}

echo "================================================"
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing study: $STUDY"
echo "Time: $(date)"
echo "================================================"

# Paths
CONTAINER="/share/crsp/lab/dalawson/share/0_Resources/containers/r_seurat_4.4.2.sif"
WORKDIR="/share/crsp/lab/dalawson/share/0_Resources/ClaudeForBioinformatics_NW/Workshop/working"

cd $WORKDIR
mkdir -p logs

# Run QC for this study
singularity exec $CONTAINER \
    Rscript scripts/05_qc_normalize.R \
    --study $STUDY \
    --manifest manifest.yaml

EXIT_CODE=$?

echo "================================================"
echo "Completed: $STUDY"
echo "Exit code: $EXIT_CODE"
echo "Time: $(date)"
echo "================================================"

exit $EXIT_CODE
```

## Array Job Commands

```bash
# Submit array job
sbatch jobs/run_qc_array.sh

# Check all tasks in array
squeue -u $USER -o "%.10i %.9P %.20j %.8u %.2t %.10M %.6D %R"

# Cancel entire array
scancel <jobid>

# Cancel specific task (e.g., task 1)
scancel <jobid>_1

# Cancel range of tasks
scancel <jobid>_[0-1]

# Check status of all tasks
sacct -j <jobid> --format=JobID,State,ExitCode,Elapsed
```

## Advanced: Manifest-Driven Arrays

For more flexibility, read study names from the manifest:

```bash
#!/bin/bash
#SBATCH --array=0-1  # Update based on number of studies

# Read studies from manifest
MANIFEST="manifest.yaml"
STUDIES=($(grep -A 100 "studies:" $MANIFEST | grep "^    [a-z]" | cut -d: -f1 | tr -d ' '))

# Or hardcode for simplicity
STUDIES=(gray pal)

STUDY=${STUDIES[$SLURM_ARRAY_TASK_ID]}
```

## Controlling Parallel Execution

```bash
# Limit concurrent tasks (e.g., max 2 at a time)
#SBATCH --array=0-9%2

# Skip certain indices
#SBATCH --array=0,2,4,6  # Only even indices

# Range with step
#SBATCH --array=0-10:2   # 0,2,4,6,8,10
```

## Discussion Questions

1. When would you use array jobs vs separate job scripts?
2. How do you handle dependencies between array tasks?
3. What if different studies need different resources?

## Common Patterns

### Reading from a file list
```bash
# studies.txt contains one study name per line
STUDY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" studies.txt)
```

### Using manifest with yq
```bash
# If yq is available
STUDY=$(yq e ".data.studies | keys | .[$SLURM_ARRAY_TASK_ID]" manifest.yaml)
```

## Checkpoint Verification

Before moving on, verify:
- [ ] Array job script created
- [ ] Understand array index to study mapping
- [ ] Successfully submitted array job
- [ ] Can monitor individual tasks
- [ ] Both studies processed successfully
- [ ] Processed objects saved to outputs/

## Next Steps

In Session 4, you'll integrate these processed objects and generate reports.
