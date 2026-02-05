#!/bin/bash
# Gate 2: Post-ExitPlanMode sync hook
# Archives stale plans, writes new plan, rebuilds index
#
# Project-aware: Detects current project from git root and uses its .plans/

set -e

DRAFT_DIR="$HOME/.claude/plans"
FORCE_FILE="$HOME/.claude/plans/.force"

# Force bypass: if .force file exists, skip sync entirely
# This keeps the plan ephemeral (not archived or indexed)
if [[ -f "$FORCE_FILE" ]]; then
    rm -f "$FORCE_FILE"
    echo "Force bypass active - skipping sync (plan is ephemeral)"
    exit 0
fi

# Detect project root from git or find .plans directory
detect_project_plans_dir() {
    # Try git root first
    local git_root
    git_root=$(git rev-parse --show-toplevel 2>/dev/null || echo "")

    if [[ -n "$git_root" ]]; then
        # Check for .plans in git root or project/ subdirectory
        if [[ -d "$git_root/.plans" ]]; then
            echo "$git_root/.plans"
            return
        fi
        if [[ -d "$git_root/project/.plans" ]]; then
            echo "$git_root/project/.plans"
            return
        fi
    fi

    # Walk up from CWD looking for .plans
    local current="$PWD"
    while [[ "$current" != "/" ]]; do
        if [[ -d "$current/.plans" ]]; then
            echo "$current/.plans"
            return
        fi
        if [[ -d "$current/project/.plans" ]]; then
            echo "$current/project/.plans"
            return
        fi
        current=$(dirname "$current")
    done

    # No .plans found
    echo ""
}

# Detect project .plans directory
PLANS_DIR=$(detect_project_plans_dir)

if [[ -z "$PLANS_DIR" ]]; then
    echo "No .plans directory found in current project - skipping sync"
    exit 0
fi

ACTIVE_DIR="$PLANS_DIR/active_plans"
ARCHIVE_DIR="$PLANS_DIR/.archive"
TOOLS_DIR="$PLANS_DIR/tools"

# Find the most recently modified plan draft
LATEST_PLAN=$(find "$DRAFT_DIR" -name "*.md" -type f -mmin -5 2>/dev/null | head -1)

if [[ -z "$LATEST_PLAN" ]]; then
    echo "No recent plan drafts to sync"
    exit 0
fi

# Extract frontmatter fields
extract_field() {
    local field=$1
    local file=$2
    grep "^${field}:" "$file" | head -1 | sed "s/^${field}:[[:space:]]*//"
}

AREA=$(extract_field "area" "$LATEST_PLAN")
NAME=$(extract_field "name" "$LATEST_PLAN")
VERSION=$(extract_field "version" "$LATEST_PLAN")
UPDATED=$(extract_field "updated" "$LATEST_PLAN")

if [[ -z "$AREA" || -z "$NAME" ]]; then
    echo "ERROR: Could not extract area/name from plan"
    exit 1
fi

# Determine target path
TARGET_DIR="$ACTIVE_DIR/$AREA"
TARGET_PATH="$TARGET_DIR/${NAME}.md"

# Check if there's an existing plan to archive
if [[ -f "$TARGET_PATH" ]]; then
    OLD_UPDATED=$(extract_field "updated" "$TARGET_PATH")
    OLD_VERSION=$(extract_field "version" "$TARGET_PATH")

    # Only archive if the dates differ
    if [[ "$OLD_UPDATED" != "$UPDATED" ]]; then
        ARCHIVE_NAME="${OLD_UPDATED}_${AREA}_${NAME}.md"
        echo "Archiving previous version: $ARCHIVE_NAME"
        mkdir -p "$ARCHIVE_DIR"
        cp "$TARGET_PATH" "$ARCHIVE_DIR/$ARCHIVE_NAME"
    fi
fi

# Write new plan
echo "Writing plan to: $TARGET_PATH"
mkdir -p "$TARGET_DIR"
cp "$LATEST_PLAN" "$TARGET_PATH"

# Rebuild index
echo "Rebuilding index..."
if command -v python3 &> /dev/null && [[ -f "$TOOLS_DIR/graph.py" ]]; then
    python3 "$TOOLS_DIR/graph.py" "$PLANS_DIR" --rebuild
else
    echo "WARNING: Python3 or graph.py not available, skipping index rebuild"
fi

echo "Plan synced successfully: $AREA/$NAME (v$VERSION)"
