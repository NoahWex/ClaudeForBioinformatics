# Container Information for UCI HPC3

This document describes the Singularity containers available for workshop exercises.

## Primary R Container

**Path**: `/share/crsp/lab/dalawson/share/0_Resources/containers/r_seurat_4.4.2.sif`

### Contents
- R 4.4.2
- Seurat 5.x
- Harmony
- SeuratWrappers
- Common Bioconductor packages

### Usage in SLURM Scripts

```bash
#!/bin/bash
#SBATCH --job-name=seurat_analysis
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

module load singularity

CONTAINER="/share/crsp/lab/dalawson/share/0_Resources/containers/r_seurat_4.4.2.sif"

singularity exec $CONTAINER Rscript my_script.R
```

### Interactive Usage

For testing scripts interactively:

```bash
# Request interactive node
srun --pty --account=dalawson_lab --partition=standard --time=1:00:00 --mem=16G /bin/bash

# Load singularity
module load singularity

# Start R in container
singularity exec /share/crsp/lab/dalawson/share/0_Resources/containers/r_seurat_4.4.2.sif R
```

### Binding Paths

The container automatically binds common paths. For custom paths:

```bash
singularity exec --bind /custom/path:/custom/path $CONTAINER Rscript script.R
```

## Module System Alternative

If you prefer modules over containers:

```bash
module load R/4.4.2
module load hdf5/1.14.3  # For h5 file support
```

Note: Module-based R may have different package versions than the container.

## Checking Available Packages

Inside the container:

```r
# List all installed packages
installed.packages()[, c("Package", "Version")]

# Check specific package
packageVersion("Seurat")
```

## Common Issues

### Memory Errors
If you get memory errors, increase `--mem` in your SLURM script:
```bash
#SBATCH --mem=64G
```

### Missing Packages
If a package is missing from the container, you have options:
1. Install to a local library: `R_LIBS_USER=~/R/libs Rscript -e "install.packages('pkg')"`
2. Request it be added to the lab container
3. Use a different container

### Permission Denied
Ensure your paths are accessible from HPC3:
- CRSP paths: `/share/crsp/lab/dalawson/...`
- DFS paths: `/dfs3b/dalawson-lab/...`

## Resources

- [HPC3 Singularity Documentation](https://rcic.uci.edu/hpc3/singularity.html)
- [Seurat Documentation](https://satijalab.org/seurat/)
- [Lab Container Recipes](https://github.com/dalawson-lab/containers)
