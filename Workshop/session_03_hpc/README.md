# Session 3: Processing Pipeline & HPC

**Duration**: 3 hours
**Theme**: Containerized workflows, SLURM jobs, parallel processing

## Learning Objectives

By the end of this session, participants will be able to:

1. **Understand containers** - Know when and why to use Singularity containers
2. **Write SLURM job scripts** - Create properly formatted batch scripts
3. **Use the HPC toolkit** - Leverage Claude's HPC-aware tools
4. **Design array jobs** - Process multiple samples in parallel
5. **Monitor and troubleshoot jobs** - Use squeue, sacct, and logs

## Schedule

| Time | Segment | Type | Notes |
|------|---------|------|-------|
| 0:00 | Container concepts | Demo | Why containers, Singularity basics |
| 0:15 | QC pipeline script | Hands-on | Exercise 1 |
| 0:55 | Break | - | 10 minutes |
| 1:05 | SLURM job writing | Demo + Hands-on | Exercise 2 |
| 1:40 | Array jobs for parallelization | Demo + Hands-on | Exercise 3 |
| 2:10 | Job submission & monitoring | Hands-on | Live job submission |
| 2:30 | Checkpoint review | Discussion | Verify outputs |
| 2:50 | Wrap-up | - | Preview Day 4 |

## Prerequisites

Completed Session 2, or reset to checkpoint:
```bash
../recovery/sync_checkpoint.sh cp2_explored
```

Also need:
- SSH access to HPC3
- Familiarity with basic shell commands

## Key Claude Features Introduced

| Feature | Usage |
|---------|-------|
| `/hpc` skill | HPC-aware operations |
| HPC toolkit | Job submission, monitoring |
| Bash tool | squeue, sacct, module commands |

## HPC Quick Reference

### Connecting to HPC3
```bash
ssh hpc3.rcic.uci.edu
```

### Key SLURM Commands
```bash
sbatch script.sh     # Submit job
squeue -u $USER      # Check your jobs
scancel <jobid>      # Cancel job
sacct -j <jobid>     # Job history/resources
```

### Module System
```bash
module load singularity  # Load Singularity
module load R/4.4.2      # Alternative to container
module list              # Show loaded modules
```

## Exercises

1. **[QC Pipeline](exercises/01_qc_pipeline.md)** - Write R script for QC and normalization
2. **[SLURM Jobs](exercises/02_slurm_jobs.md)** - Create and submit a batch job
3. **[Array Jobs](exercises/03_array_jobs.md)** - Process multiple samples in parallel

## Templates

Pre-built templates in `templates/`:
- `single_job.sh` - Basic single-task job
- `array_job.sh` - Manifest-driven array job

## Checkpoint

At the end of this session, you should have:
- QC + normalization R script
- Tested SLURM job scripts
- Processed (normalized) Seurat objects
- Understanding of job monitoring

Use `../recovery/sync_checkpoint.sh cp3_processed` if you need to catch up.

## Container Information

See `../shared/container_info.md` for details on:
- Available containers
- How to use them in jobs
- Common issues and solutions

## SLURM Resource Guidelines

| Job Type | Time | Memory | CPUs |
|----------|------|--------|------|
| Quick test | 30m | 8G | 1 |
| Single sample QC | 2h | 32G | 4 |
| Large integration | 8h | 128G | 16 |

## Common Issues

### Job stuck in PENDING
```bash
squeue -j <jobid> -o "%R"  # Show reason
```
Common reasons: Resources unavailable, account limits

### Out of memory
Increase `--mem` in SLURM script. Check actual usage:
```bash
sacct -j <jobid> --format=JobID,MaxRSS,MaxVMSize
```

### Module not found
Use `module avail <name>` to find correct module name

### Container not found
Verify path exists: `ls -la /share/crsp/lab/dalawson/share/0_Resources/containers/`

## Resources

- [UCI HPC3 Documentation](https://rcic.uci.edu/hpc3/)
- [SLURM Documentation](https://slurm.schedmd.com/documentation.html)
- [Singularity User Guide](https://docs.sylabs.io/guides/latest/user-guide/)
- Container info: `../shared/container_info.md`
