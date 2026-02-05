# Project-Specific Rules

Rules loaded for this project.

## HPC Execution

Follow these practices for HPC jobs:

1. **Use absolute paths** in SLURM scripts
2. **Always load singularity** before container exec
3. **Don't bind /scratch** - not available on all nodes
4. **Use --cleanenv** for R containers
5. **Set temp directories** for Python (NUMBA_CACHE_DIR, etc.)

## Data Management

- Register all inputs in `project/config/data_manifest.yaml`
- Register all outputs in `project/config/outputs_manifest.yaml`
- Never hardcode paths in scripts

## Code Style

- Include header comments explaining script purpose
- Use set -euo pipefail in bash scripts
- Include input validation before heavy computation
