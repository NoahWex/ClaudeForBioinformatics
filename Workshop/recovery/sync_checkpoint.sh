#!/bin/bash
# sync_checkpoint.sh - Reset participant workspace to a checkpoint
#
# Usage: ./sync_checkpoint.sh <checkpoint_name>
# Example: ./sync_checkpoint.sh cp1_data_loaded
#
# Available checkpoints:
#   cp1_data_loaded  - Day 1 complete: Project scaffold + loaded data
#   cp2_explored   - Day 2 complete: Harmonized metadata
#   cp3_processed    - Day 3 complete: Processed Seurat objects
#   cp4_complete   - Day 4 complete: Integrated object + report

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSHOP_DIR="$(dirname "$SCRIPT_DIR")"

# Checkpoint mapping to session directories
declare -A CHECKPOINT_PATHS
CHECKPOINT_PATHS["cp1_data_loaded"]="session_01_foundations/checkpoint/cp1_data_loaded"
CHECKPOINT_PATHS["cp2_explored"]="session_02_metadata/checkpoint/cp2_explored"
CHECKPOINT_PATHS["cp3_processed"]="session_03_hpc/checkpoint/cp3_processed"
CHECKPOINT_PATHS["cp4_complete"]="session_04_integration/checkpoint/cp4_complete"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

usage() {
    echo "Usage: $0 <checkpoint_name> [target_directory]"
    echo ""
    echo "Available checkpoints:"
    for cp in "${!CHECKPOINT_PATHS[@]}"; do
        echo "  - $cp"
    done
    echo ""
    echo "If target_directory is not specified, uses ./working/"
    exit 1
}

if [ -z "$1" ]; then
    usage
fi

CHECKPOINT="$1"
TARGET_DIR="${2:-./working}"

# Validate checkpoint name
if [ -z "${CHECKPOINT_PATHS[$CHECKPOINT]}" ]; then
    echo -e "${RED}Error: Unknown checkpoint '$CHECKPOINT'${NC}"
    usage
fi

CHECKPOINT_PATH="$WORKSHOP_DIR/${CHECKPOINT_PATHS[$CHECKPOINT]}"

# Check checkpoint exists
if [ ! -d "$CHECKPOINT_PATH" ]; then
    echo -e "${RED}Error: Checkpoint directory not found: $CHECKPOINT_PATH${NC}"
    echo "The checkpoint may not have been created yet."
    exit 1
fi

# Warn if target exists
if [ -d "$TARGET_DIR" ]; then
    echo -e "${YELLOW}Warning: Target directory exists: $TARGET_DIR${NC}"
    read -p "This will overwrite existing files. Continue? [y/N] " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Aborted."
        exit 1
    fi
fi

# Create target if needed
mkdir -p "$TARGET_DIR"

# Sync checkpoint to target
echo -e "${GREEN}Syncing checkpoint '$CHECKPOINT' to '$TARGET_DIR'...${NC}"
rsync -av --delete "$CHECKPOINT_PATH/" "$TARGET_DIR/"

# Verify sync
if [ $? -eq 0 ]; then
    echo ""
    echo -e "${GREEN}Success! Reset to checkpoint: $CHECKPOINT${NC}"
    echo "Working directory: $TARGET_DIR"
    echo ""
    echo "Next steps:"
    case $CHECKPOINT in
        cp1_data_loaded)
            echo "  - Continue with Session 2: Metadata exploration and harmonization"
            ;;
        cp2_explored)
            echo "  - Continue with Session 3: QC pipeline and HPC processing"
            ;;
        cp3_processed)
            echo "  - Continue with Session 4: Integration and visualization"
            ;;
        cp4_complete)
            echo "  - Workshop complete! Review final outputs."
            ;;
    esac
else
    echo -e "${RED}Error during sync. Please check the output above.${NC}"
    exit 1
fi
