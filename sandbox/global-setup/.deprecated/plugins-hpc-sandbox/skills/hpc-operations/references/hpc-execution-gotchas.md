# HPC Execution Gotchas

Hard-won lessons from container and SLURM debugging on UCI HPC3.

## Pipeline Correctness Principle

**Deterministic, linear pipelines must produce correct results on first execution.**

- No downstream patches to fix incorrect outputs
- Verify logic before submitting HPC jobs
- Trace through code paths to confirm expected behavior
- If outputs are wrong, fix the source and re-run - never post-process

**Before running any pipeline:**
1. Verify input → output mapping is correct
2. Check that lookup tables have matching keys
3. Confirm file paths resolve correctly on HPC
4. Test with a small subset if possible

## Module Loading

### Always load singularity before exec

```bash
# WRONG - singularity command not found
singularity exec ...

# CORRECT - always load module first
module load singularity
singularity exec ...
```

This is especially easy to miss in SLURM scripts where the login node has singularity but compute nodes may not.

## Container Mount Failures

### `/scratch` does not exist on all nodes

```
FATAL: container creation failed: mount /scratch->/scratch error: mount source /scratch doesn't exist
```

**Fix**: Remove `/scratch` from bind mounts. Use only:
```bash
--bind /share/crsp/lab/dalawson/share:/share/crsp/lab/dalawson/share:ro \
--bind /share/crsp/lab/dalawson/nwechter:/share/crsp/lab/dalawson/nwechter:rw \
--bind /dfs7:/dfs7:ro \
--bind /dfs8:/dfs8:ro
```

### Missing R packages (miloR, argparse, etc.)

```
Error in library(miloR) : there is no package called 'miloR'
```

**Fix**: Bind R_LIBS_USER and set environment:
```bash
R_LIBS_USER="/data/homezvol0/nwechter/biojhub4_dir/Rocky8_jupyter_base_R4.3.3_Spatial.sif/R/library"

singularity exec \
    --bind "$R_LIBS_USER:/home/jovyan/R/library:ro" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    ...
```

## SLURM Log Paths

### Relative paths go to sbatch working directory

```bash
# WRONG - goes wherever you run sbatch from
#SBATCH --output=logs/job_%j.out

# CORRECT - always goes to the right place
#SBATCH --output=/share/crsp/lab/dalawson/nwechter/Project/logs/job_%j.out
```

## Python Cache Directories

Container paths are read-only. Redirect caches:
```bash
--env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
--env "MPLCONFIGDIR=/tmp/matplotlib_config"
```

## CRSP Filesystem

### shutil.copy2 fails

```python
# WRONG - xattr not supported on CRSP
shutil.copy2(src, dst)

# CORRECT
shutil.copy(src, dst)
```

### Stale file handle with gzip/HDF5

Writing large files then immediately reading them (e.g., gzip after writeMM) causes "Stale file handle" errors on CRSP.

```r
# WRONG - stale file handle
writeMM(counts, "/crsp/path/matrix.mtx")
system2("gzip", c("-f", "/crsp/path/matrix.mtx"))  # fails

# CORRECT - write to /tmp, then copy
writeMM(counts, "/tmp/matrix.mtx")
system2("gzip", c("-f", "/tmp/matrix.mtx"))
file.copy("/tmp/matrix.mtx.gz", "/crsp/path/matrix.mtx.gz")
```

Same pattern applies to HDF5 files with rhdf5 - write to /tmp first.

## R source() Gotchas

### CLI blocks execute when source()'d

Scripts with `if (!interactive()) { ... }` CLI sections will execute when source()'d in a container (non-interactive). Command line args from the parent script leak through.

```r
# WRONG - triggers on any --all arg
if (!interactive()) {
  if ("--all" %in% args) harmonize_all_studies()
}

# CORRECT - check if this script is the main script
.is_main_script <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  script_arg <- args[grep("--file=", args)]
  if (length(script_arg) == 0) return(FALSE)
  basename(sub("--file=", "", script_arg)) == "this_script.R"
}

if (!interactive() && .is_main_script()) { ... }
```

## Seurat Compatibility

### v3 objects fail in v5

```
Error: no slot of name "images" for this object of class "Seurat"
```

**Fix**:
```r
srt <- readRDS("old_seurat.rds")
srt <- UpdateSeuratObject(srt)
```

## Cell ID Format Gotchas

### 10x barcode `-1` suffix

Barcodes have `-N` suffix: `AAACCTGAGTATGACA-1`

```r
# WRONG - validates to 18 chars, fails
nchar(barcode_full) == 16

# CORRECT - strip suffix first
barcode <- sub("-[0-9]+$", "", barcode_full)
nchar(barcode) == 16
```

### Inconsistent formats across datasets

Always inspect actual data before assuming format:
```r
head(colnames(seurat_obj))
```
