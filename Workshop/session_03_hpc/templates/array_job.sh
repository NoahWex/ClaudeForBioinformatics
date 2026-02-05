#!/bin/bash
#===============================================================================
# SLURM Array Job Template
# Workshop: Claude for Bioinformatics
#
# Usage: sbatch array_job.sh
#
# This template processes multiple studies in parallel using SLURM arrays.
# Each array task processes one study based on its index.
#===============================================================================

#-------------------------------------------------------------------------------
# SLURM Resource Requests
#-------------------------------------------------------------------------------
#SBATCH --job-name=qc_array
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

#-------------------------------------------------------------------------------
# Array Configuration
#-------------------------------------------------------------------------------
# Array indices 0 and 1 for two studies
# Add %N to limit concurrent tasks: --array=0-1%1 runs one at a time
#SBATCH --array=0-1

#-------------------------------------------------------------------------------
# Output Logging
#-------------------------------------------------------------------------------
# %A = array job ID, %a = array task ID
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

#-------------------------------------------------------------------------------
# Notifications (optional)
#-------------------------------------------------------------------------------
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=your_email@uci.edu

#===============================================================================
# Environment Setup
#===============================================================================

set -e
module load singularity

#-------------------------------------------------------------------------------
# Configuration
#-------------------------------------------------------------------------------
CONTAINER="/share/crsp/lab/dalawson/share/0_Resources/containers/r_seurat_4.4.2.sif"
WORKDIR="/share/crsp/lab/dalawson/share/0_Resources/ClaudeForBioinformatics_NW/Workshop/working"
SCRIPT="scripts/05_qc_normalize.R"

#-------------------------------------------------------------------------------
# Study Mapping
# IMPORTANT: Order must match array indices!
# Index 0 -> gray
# Index 1 -> pal
#-------------------------------------------------------------------------------
STUDIES=(gray pal)

# Get study name for this array task
STUDY=${STUDIES[$SLURM_ARRAY_TASK_ID]}

# Validate study exists
if [ -z "$STUDY" ]; then
    echo "ERROR: No study mapped to array index $SLURM_ARRAY_TASK_ID"
    exit 1
fi

#===============================================================================
# Job Execution
#===============================================================================

echo "========================================"
echo "Array Job: $SLURM_JOB_NAME"
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "Processing Study: $STUDY"
echo "Node: $SLURMD_NODENAME"
echo "Started: $(date)"
echo "========================================"

cd $WORKDIR
mkdir -p logs

echo "Container: $CONTAINER"
echo "Script: $SCRIPT"
echo "========================================"

# Run the analysis for this study
singularity exec $CONTAINER \
    Rscript $SCRIPT \
    --study $STUDY \
    --manifest manifest.yaml

EXIT_CODE=$?

echo "========================================"
echo "Finished: $STUDY"
echo "Time: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================"

exit $EXIT_CODE
