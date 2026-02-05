#!/bin/bash --noprofile
# Gate 1: Pre-ExitPlanMode validation hook
# Validates plan structure before user review
#
# Project-aware: Only validates if current project has .plans/ directory
#
# Uses JSON output format for Claude Code hooks:
# - permissionDecision: "deny" blocks with explanation
# - permissionDecision: "allow" proceeds
# - Always exit 0; use JSON structure for control

set -euo pipefail

# JSON output helpers
json_deny() {
    local reason="$1"
    local context="${2:-$1}"
    cat << EOF
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "deny",
    "permissionDecisionReason": "$reason",
    "additionalContext": "$context"
  }
}
EOF
}

json_allow() {
    cat << EOF
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "allow"
  }
}
EOF
}

# Escape strings for JSON embedding (strips trailing newline)
json_escape() {
    python3 -c 'import sys,json; print(json.dumps(sys.stdin.read().rstrip("\n"))[1:-1])'
}

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

# Consume stdin (Claude Code passes tool params as JSON)
INPUT=$(cat)

# Check if current project has a .plans directory
PROJECT_PLANS=$(detect_project_plans_dir)
if [[ -z "$PROJECT_PLANS" ]]; then
    # No .plans in current project - skip validation
    json_allow
    exit 0
fi

# Force bypass: if .force file exists, skip validation
# Note: Don't delete here - post_exit_plan.sh will delete after skipping sync
FORCE_FILE="$HOME/.claude/plans/.force"
if [[ -f "$FORCE_FILE" ]]; then
    json_allow
    exit 0
fi

# Find the most recently modified plan in ~/.claude/plans/
DRAFT_DIR="$HOME/.claude/plans"

if [[ ! -d "$DRAFT_DIR" ]]; then
    json_allow
    exit 0
fi

# Get the most recently modified .md file (within last 5 minutes)
LATEST_PLAN=$(find "$DRAFT_DIR" -name "*.md" -type f -mmin -5 2>/dev/null | head -1)

if [[ -z "$LATEST_PLAN" ]]; then
    json_allow
    exit 0
fi

# Validation function - returns error message or empty string
# Fixed: Now extracts actual YAML frontmatter instead of grepping entire file
# (prevents false positives from code block examples)
validate_plan() {
    local plan_file="$1"

    # Check file starts with frontmatter delimiter
    if ! head -1 "$plan_file" | grep -q "^---$"; then
        echo "Missing YAML frontmatter (file must start with ---)"
        return
    fi

    # Extract frontmatter (content between first two --- lines)
    local frontmatter
    frontmatter=$(awk '/^---$/{if(++c==1)next; if(c==2)exit} c==1' "$plan_file")

    if [[ -z "$frontmatter" ]]; then
        echo "Empty or malformed YAML frontmatter"
        return
    fi

    # Check required frontmatter fields within extracted block only
    local fields=("area" "name" "version" "updated")
    for field in "${fields[@]}"; do
        if ! echo "$frontmatter" | grep -q "^${field}:"; then
            echo "Missing required frontmatter field: $field"
            return
        fi
    done

    # Check for required sections in body (grep still works here - sections are unique headers)
    local sections=("## Inputs" "## Outputs" "## Config" "## Scripts" "## Tasks" "## Current State" "## Discussion" "## Next Actions")

    for section in "${sections[@]}"; do
        if ! grep -q "^$section" "$plan_file"; then
            echo "Missing required section: $section"
            return
        fi
    done

    # All checks passed
    echo ""
}

# Run validation
ERROR=$(validate_plan "$LATEST_PLAN")

if [[ -n "$ERROR" ]]; then
    REASON=$(echo "Plan validation: $ERROR" | json_escape)
    CONTEXT=$(echo "Fix in $LATEST_PLAN then retry ExitPlanMode" | json_escape)
    json_deny "$REASON" "$CONTEXT"
    exit 0
fi

# Validation passed
json_allow
exit 0
