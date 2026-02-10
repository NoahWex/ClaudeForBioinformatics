#!/bin/bash --noprofile
# Pre-ExitPlanMode validation hook
# Validates plan structure before user review.
#
# Finds the draft plan via draft_path frontmatter (deterministic) with
# a -mmin -15 recency fallback for plans written without draft_path.
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
    echo "BLOCKED: $reason" >&2
    echo "  $context" >&2
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

json_escape() {
    python3 -c 'import sys,json; print(json.dumps(sys.stdin.read().rstrip("\n"))[1:-1])'
}

# Detect project plans directory (walks up from CWD)
detect_project_plans_dir() {
    local git_root
    git_root=$(git rev-parse --show-toplevel 2>/dev/null || echo "")
    if [[ -n "$git_root" ]]; then
        [[ -d "$git_root/.plans" ]] && echo "$git_root/.plans" && return
        [[ -d "$git_root/project/.plans" ]] && echo "$git_root/project/.plans" && return
    fi
    local current="$PWD"
    while [[ "$current" != "/" ]]; do
        [[ -d "$current/.plans" ]] && echo "$current/.plans" && return
        [[ -d "$current/project/.plans" ]] && echo "$current/project/.plans" && return
        current=$(dirname "$current")
    done
    echo ""
}

# Consume stdin
INPUT=$(cat)

# Only validate if current project uses the plan system
PROJECT_PLANS=$(detect_project_plans_dir)
[[ -z "$PROJECT_PLANS" ]] && { json_allow; exit 0; }

# Force bypass
FORCE_FILE="$HOME/.claude/plans/.force"
[[ -f "$FORCE_FILE" ]] && { json_allow; exit 0; }

# Find the draft plan to validate
DRAFT_DIR="$HOME/.claude/plans"
[[ ! -d "$DRAFT_DIR" ]] && { json_allow; exit 0; }

LATEST_PLAN=""

# Primary: find file where draft_path == own path (deterministic)
for candidate in "$DRAFT_DIR"/*.md; do
    [[ -f "$candidate" ]] || continue
    dp=$(awk '/^---$/{if(++c==1)next; if(c==2)exit} c==1 && /^draft_path:/{sub(/^draft_path:[[:space:]]*/, ""); print}' "$candidate" 2>/dev/null)
    [[ -z "$dp" ]] && continue
    real_dp=$(realpath "$dp" 2>/dev/null || echo "$dp")
    real_candidate=$(realpath "$candidate" 2>/dev/null || echo "$candidate")
    if [[ "$real_dp" == "$real_candidate" ]]; then
        LATEST_PLAN="$candidate"
        break
    fi
done

# Fallback: most recently modified .md within 15 minutes
if [[ -z "$LATEST_PLAN" ]]; then
    LATEST_PLAN=$(find "$DRAFT_DIR" -maxdepth 1 -name "*.md" -type f -mmin -15 2>/dev/null \
        | xargs ls -t 2>/dev/null | head -1)
fi

[[ -z "$LATEST_PLAN" ]] && { json_allow; exit 0; }

# Validate plan structure
validate_plan() {
    local plan_file="$1"

    if ! head -1 "$plan_file" | grep -q "^---$"; then
        echo "Missing YAML frontmatter (file must start with ---)"
        return
    fi

    local frontmatter
    frontmatter=$(awk '/^---$/{if(++c==1)next; if(c==2)exit} c==1' "$plan_file")
    [[ -z "$frontmatter" ]] && { echo "Empty or malformed YAML frontmatter"; return; }

    # Required fields with non-empty values
    local required_fields=("area" "name" "version" "updated" "status")
    for field in "${required_fields[@]}"; do
        if ! echo "$frontmatter" | grep -q "^${field}:"; then
            echo "Missing required frontmatter field: $field"
            return
        fi
    done

    # Validate status enum
    local status_val
    status_val=$(echo "$frontmatter" | grep "^status:" | sed 's/^status:[[:space:]]*//')
    case "$status_val" in
        queued|in_progress|blocked|complete) ;;
        *)
            echo "Invalid status: '$status_val' (expected: queued, in_progress, blocked, complete)"
            return
            ;;
    esac

    # Fields that must be present (any value including null/[])
    local present_fields=("supersedes" "blocked_by" "blocks" "jobs")
    for field in "${present_fields[@]}"; do
        if ! echo "$frontmatter" | grep -q "^${field}:"; then
            echo "Missing frontmatter field: $field"
            return
        fi
    done

    # Required sections
    local sections=("## Inputs" "## Outputs" "## Config" "## Scripts" "## Tasks" "## Current State" "## Discussion" "## Next Actions" "## Completion")
    for section in "${sections[@]}"; do
        if ! grep -q "^$section" "$plan_file"; then
            echo "Missing required section: $section"
            return
        fi
    done

    echo ""
}

# Run validation
ERROR=$(validate_plan "$LATEST_PLAN")

if [[ -n "$ERROR" ]]; then
    echo "BLOCKED: Plan validation: $ERROR" >&2
    echo "  Expected frontmatter fields:" >&2
    echo "    REQUIRED: area, name, version, updated, status" >&2
    echo "    PRESENT:  supersedes, blocked_by, blocks, jobs" >&2
    echo "    OPTIONAL: draft_path, sequence, plan_target" >&2
    echo "  Fix in: $LATEST_PLAN" >&2
    REASON=$(echo "Plan validation: $ERROR" | json_escape)
    CONTEXT=$(echo "Fix in $LATEST_PLAN then retry ExitPlanMode" | json_escape)
    json_deny "$REASON" "$CONTEXT"
    exit 0
fi

json_allow
exit 0
