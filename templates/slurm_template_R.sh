#!/bin/bash
#SBATCH --job-name=experiment_analysis
#SBATCH --partition=standard
#SBATCH --account=lab_account
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=/full/path/to/experiment/.dev/logs/analysis_%j.out
#SBATCH --error=/full/path/to/experiment/.dev/logs/analysis_%j.err

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
CONTAINER="/path/to/R/container.sif"
R_LIBS_USER="/path/to/user/R/library"

# Script to run
SCRIPT="$PROJECT_ROOT/project/02_Analyses/analyses/ResearchArea/experiment_YYYYMMDD/scripts/analysis.R"

# =============================================================================
# BANNER
# =============================================================================

echo "============================================================================"
echo "Experiment Analysis"
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
# TEMP DIRECTORY SETUP (CRITICAL for Rmd rendering)
# =============================================================================

PIPELINE_TEMP_ROOT="$PROJECT_ROOT/.pipeline_temp"
mkdir -p "$PIPELINE_TEMP_ROOT"
JOB_TEMP_DIR=$(mktemp -d "$PIPELINE_TEMP_ROOT/job_${SLURM_JOB_ID}_XXXXXX")

# Cleanup trap
trap 'echo "CLEANUP: Removing $JOB_TEMP_DIR"; rm -rf "$JOB_TEMP_DIR"' EXIT

export TMPDIR="$JOB_TEMP_DIR"
export TMP="$JOB_TEMP_DIR"
export TEMP="$JOB_TEMP_DIR"

echo "Temp directory:   $JOB_TEMP_DIR"
echo "============================================================================"

# =============================================================================
# EXECUTION
# =============================================================================

echo ""
echo "Starting analysis at $(date '+%Y-%m-%d %H:%M:%S')"
echo ""

singularity exec \
    --cleanenv --containall --pwd "$PROJECT_ROOT" \
    --bind "$R_LIBS_USER:/home/jovyan/R/library:ro" \
    --bind "$PROJECT_ROOT:$PROJECT_ROOT:rw" \
    --bind "/share/crsp/lab/labname/share:/share/crsp/lab/labname/share:ro" \
    --bind "$TMPDIR:$TMPDIR:rw" \
    --env "TMPDIR=$TMPDIR" \
    --env "TMP=$TMPDIR" \
    --env "TEMP=$TMPDIR" \
    --env "HOME=/home/jovyan" \
    --env "R_LIBS_USER=/home/jovyan/R/library" \
    --env "RSTUDIO_PANDOC=/usr/lib/rstudio-server/bin/quarto/bin/tools" \
    --env "LANG=en_US.UTF-8" \
    --env "LC_ALL=en_US.UTF-8" \
    --env "PROJECT_ROOT=$PROJECT_ROOT" \
    "$CONTAINER" Rscript "$SCRIPT" || {
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

# Why do R programmers make terrible DJs?
# Because they can't stop dropping their library() calls.
