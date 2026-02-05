#!/bin/bash
#SBATCH --job-name=experiment_analysis
#SBATCH --partition=standard
#SBATCH --account=lab_account
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/full/path/to/experiment/.dev/logs/analysis_%j.out
#SBATCH --error=/full/path/to/experiment/.dev/logs/analysis_%j.err

# For GPU jobs, uncomment:
# #SBATCH --partition=free-gpu
# #SBATCH --gres=gpu:1

# Exit on error
set -euo pipefail

# =============================================================================
# CONFIGURATION
# =============================================================================

PROJECT_ROOT="/full/path/to/project"
cd "$PROJECT_ROOT"

# Load singularity module (CRITICAL - not available by default on compute nodes)
module load singularity/3.11.3

# Container configuration
CONTAINER="/path/to/python/container.sif"

# Script to run
SCRIPT="$PROJECT_ROOT/project/02_Analyses/analyses/ResearchArea/experiment_YYYYMMDD/scripts/analysis.py"

# =============================================================================
# BANNER
# =============================================================================

echo "============================================================================"
echo "Experiment Analysis (Python)"
echo "============================================================================"
echo "Job ID:           ${SLURM_JOB_ID}"
echo "Node:             ${SLURMD_NODENAME}"
echo "CPUs:             ${SLURM_CPUS_PER_TASK:-1}"
echo "Memory:           ${SLURM_MEM_PER_NODE:-N/A}"
echo "Project Root:     $PROJECT_ROOT"
echo "Container:        $CONTAINER"
echo "Script:           $SCRIPT"
echo "============================================================================"

# =============================================================================
# VALIDATION
# =============================================================================

if [ ! -f "$SCRIPT" ]; then
    echo "ERROR: Script not found: $SCRIPT"
    exit 1
fi

if [ ! -f "$CONTAINER" ]; then
    echo "ERROR: Container not found: $CONTAINER"
    exit 1
fi

# =============================================================================
# EXECUTION
# =============================================================================

echo ""
echo "Starting analysis at $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

# For GPU jobs, add --nv flag to singularity exec

singularity exec \
    --pwd "$PROJECT_ROOT" \
    --bind /dfs7:/dfs7:ro \
    --bind /dfs8:/dfs8:ro \
    --bind "/share/crsp/lab/labname:/share/crsp/lab/labname:rw" \
    --env "PROJECT_ROOT=$PROJECT_ROOT" \
    --env "NUMBA_CACHE_DIR=/tmp/numba_cache" \
    --env "MPLCONFIGDIR=/tmp/matplotlib_config" \
    "$CONTAINER" python3 "$SCRIPT" || {
        EXIT_CODE=$?
        echo ""
        echo "ERROR: Script failed with exit code $EXIT_CODE"
        exit $EXIT_CODE
    }

# =============================================================================
# COMPLETION
# =============================================================================

echo ""
echo "============================================================================"
echo "Analysis completed successfully at $(date '+%Y-%m-%d %H:%M:%S')"
echo "============================================================================"

# Why do Python programmers prefer dark mode?
# Because light attracts bugs.
