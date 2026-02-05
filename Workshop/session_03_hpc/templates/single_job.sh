#!/bin/bash
#===============================================================================
# SLURM Single Job Template
# Workshop: Claude for Bioinformatics
#
# Usage: sbatch single_job.sh
# Customize the variables below for your analysis.
#===============================================================================

#-------------------------------------------------------------------------------
# SLURM Resource Requests
#-------------------------------------------------------------------------------
#SBATCH --job-name=my_analysis
#SBATCH --account=dalawson_lab
#SBATCH --partition=standard
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4

#-------------------------------------------------------------------------------
# Output Logging
#-------------------------------------------------------------------------------
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
# %x = job name, %j = job ID

#-------------------------------------------------------------------------------
# Notifications (optional)
#-------------------------------------------------------------------------------
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=your_email@uci.edu

#===============================================================================
# Environment Setup
#===============================================================================

# Exit on error
set -e

# Load required modules
module load singularity

#-------------------------------------------------------------------------------
# Configuration - CUSTOMIZE THESE
#-------------------------------------------------------------------------------
CONTAINER="/share/crsp/lab/dalawson/share/0_Resources/containers/r_seurat_4.4.2.sif"
WORKDIR="/share/crsp/lab/dalawson/share/0_Resources/ClaudeForBioinformatics_NW/Workshop/working"
SCRIPT="scripts/05_qc_normalize.R"
STUDY="gray"

#===============================================================================
# Job Execution
#===============================================================================

echo "========================================"
echo "Job: $SLURM_JOB_NAME"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Started: $(date)"
echo "========================================"

# Change to working directory
cd $WORKDIR

# Create logs directory if needed
mkdir -p logs

# Print configuration
echo "Container: $CONTAINER"
echo "Script: $SCRIPT"
echo "Study: $STUDY"
echo "========================================"

# Run the analysis
singularity exec $CONTAINER \
    Rscript $SCRIPT \
    --study $STUDY \
    --manifest manifest.yaml

# Capture exit code
EXIT_CODE=$?

echo "========================================"
echo "Finished: $(date)"
echo "Exit code: $EXIT_CODE"
echo "========================================"

exit $EXIT_CODE
