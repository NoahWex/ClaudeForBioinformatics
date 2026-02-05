#!/bin/bash --noprofile
# validate-hpc-command.sh - PreToolUse security gate for HPC commands
#
# This hook uses JSON permissionDecision to:
# - BLOCK dangerous chained commands (exit 0 + deny)
# - ASK for confirmation on sensitive operations (exit 0 + ask)
# - PASS safe read-only commands (exit 0, no JSON)
#
# Exit codes:
#   0 = continue (with optional JSON decision)
#   2 = block with error message

set -euo pipefail

# ============================================================================
# Session state tracking (for first-use documentation injection)
# ============================================================================
SESSION_STATE_DIR="/tmp/.claude-hpc-hooks"
mkdir -p "$SESSION_STATE_DIR" 2>/dev/null || true
# Use PPID (parent process = Claude Code) to track session
# Falls back to daily file if PPID unavailable
SESSION_ID="${PPID:-$(date +%Y%m%d)}"
SESSION_STATE_FILE="$SESSION_STATE_DIR/docs-shown-$SESSION_ID"

is_first_hpc_command() {
    if [ ! -f "$SESSION_STATE_FILE" ]; then
        return 0  # First command
    fi
    return 1
}

mark_docs_shown() {
    touch "$SESSION_STATE_FILE" 2>/dev/null || true
}

# Quick reference docs (injected on first HPC command)
# Note: Python json.dumps handles escaping, keep this readable
read -r -d '' HPC_QUICK_REFERENCE << 'DOCS' || true
=== HPC QUICK REFERENCE (first command this session) ===

JOB SUBMISSION: Raw sbatch is BLOCKED. Use the wrapper with FULL PATH:
  ~/.claude/plugins/hpc-sandbox/bin/hpc-submit \
    --script /path/to/script.sh \
    --purpose "Why this job" \
    --outputs "/path/to/outputs"

CONDITIONAL ARGS:
  --rerun-reason "..." (if outputs exist)
  --overwrites (acknowledge overwrite)
  --risk-note "..." (if mem>=64GB or time>=24h)
  --sbatch-arg "--array=1-10" (for array jobs)

READ-ONLY (auto-allowed): squeue, sacct, sinfo, ls, cat, tail
SUPERVISED (ask): scancel, git commit/push/branch ops
BLOCKED: raw sbatch, git push --force, git reset --hard

Full docs: ~/.claude/plugins/hpc-sandbox/skills/hpc-operations/SKILL.md
=== END QUICK REFERENCE ===
DOCS

# Read JSON from stdin
INPUT=$(cat)
COMMAND=$(echo "$INPUT" | python3 -c "import sys, json; d=json.load(sys.stdin); print((d.get('tool_input') or {}).get('command', ''))" 2>/dev/null || echo "")

# Empty command? Pass through
[ -z "$COMMAND" ] && exit 0

# ============================================================================
# Check if HPC-related (only validate HPC commands)
# ============================================================================
is_hpc_command() {
    echo "$1" | grep -qE "(ssh[[:space:]]+(hpc3|hpc3\.rcic\.uci\.edu))|^s(batch|run|queue|cancel|info|acct)[[:space:]]"
}

# Not HPC? Pass through
if ! is_hpc_command "$COMMAND"; then
    exit 0
fi

# ============================================================================
# Extract remote command from ssh hpc3 "..."
# ============================================================================
extract_remote_cmd() {
    local cmd="$1"
    # Use [[:space:]] instead of \s for POSIX compatibility (macOS sed)
    if echo "$cmd" | grep -qE '^ssh[[:space:]]+(hpc3|hpc3\.rcic\.uci\.edu)[[:space:]]+'; then
        echo "$cmd" | sed -E 's/^ssh[[:space:]]+(hpc3|hpc3\.rcic\.uci\.edu)[[:space:]]+["'"'"']?//' | sed -E 's/["'"'"']?$//'
    else
        echo "$cmd"
    fi
}

REMOTE_CMD=$(extract_remote_cmd "$COMMAND")

# ============================================================================
# Parse into chain segments
# ============================================================================
parse_chain() {
    local cmd="$1"
    echo "$cmd" | python3 -c '
import sys
import re
cmd = sys.stdin.read().strip()
segments = re.split(r"\s*(?:&&|\|\||;)\s*", cmd)
for seg in segments:
    seg = seg.strip()
    if seg:
        print(seg)
'
}

# ============================================================================
# Pattern detection functions
# ============================================================================
is_dangerous() {
    local cmd="$1"
    echo "$cmd" | grep -qiE 'rm\s+(-[rf]+\s+)*[/~]' && return 0
    echo "$cmd" | grep -qiE 'rm\s+(-[rf]+\s+)*\$HOME' && return 0
    echo "$cmd" | grep -qiE 'mkfs\.' && return 0
    echo "$cmd" | grep -qiE 'dd\s+.*of=/dev/' && return 0
    echo "$cmd" | grep -qE ':\(\)\s*\{' && return 0
    echo "$cmd" | grep -qiE 'chmod\s+(-[rR]+\s+)*777\s+/' && return 0
    echo "$cmd" | grep -qE '>\s*/etc/' && return 0
    # Catch writes to device files, but exclude /dev/null (safe stderr redirect)
    echo "$cmd" | grep -qE '>\s*/dev/' && ! echo "$cmd" | grep -qE '>\s*/dev/null' && return 0
    return 1
}

is_job_submission() {
    local cmd="$1"
    echo "$cmd" | grep -qE '^sbatch\s' && return 0
    return 1
}

is_job_cancellation() {
    local cmd="$1"
    echo "$cmd" | grep -qE '^scancel\s' && return 0
    return 1
}

is_interactive_session() {
    local cmd="$1"
    echo "$cmd" | grep -qE '^srun\s' && return 0
    return 1
}

is_destructive() {
    local cmd="$1"
    echo "$cmd" | grep -qiE '^rm\s' && return 0
    echo "$cmd" | grep -qiE '^mv\s.*-f' && return 0
    return 1
}

# ============================================================================
# Git command detection
# ============================================================================
is_git_command() {
    local cmd="$1"
    echo "$cmd" | grep -qE '^git\s' && return 0
    return 1
}

# Git: safe read-only and staging operations (auto-allow)
is_git_safe() {
    local cmd="$1"
    echo "$cmd" | grep -qE '^git\s+(status|log|diff|show|ls-files|remote|config|fetch)\b' && return 0
    echo "$cmd" | grep -qE '^git\s+(add|stash|tag)\b' && return 0
    return 1
}

# Git: commits and branch management (supervised - ask)
is_git_supervised() {
    local cmd="$1"
    # Commits require supervision
    echo "$cmd" | grep -qE '^git\s+commit\b' && return 0
    # Branch operations
    echo "$cmd" | grep -qE '^git\s+(checkout|switch|branch|merge|rebase|cherry-pick)\b' && return 0
    echo "$cmd" | grep -qE '^git\s+(push|pull)\b' && return 0
    echo "$cmd" | grep -qE '^git\s+reset\b' && return 0
    return 1
}

# Git: dangerous operations (block/strong warning)
is_git_dangerous() {
    local cmd="$1"
    echo "$cmd" | grep -qE '^git\s+push\s+.*(\-f|\-\-force)' && return 0
    echo "$cmd" | grep -qE '^git\s+reset\s+.*\-\-hard' && return 0
    echo "$cmd" | grep -qE '^git\s+clean\s+.*\-f' && return 0
    echo "$cmd" | grep -qE '^git\s+branch\s+.*\-D' && return 0
    return 1
}

# ============================================================================
# JSON response helpers
# ============================================================================
json_deny() {
    local reason="$1"
    local context="$2"
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

json_ask() {
    local reason="$1"
    local context="$2"
    cat << EOF
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "ask",
    "permissionDecisionReason": "$reason",
    "additionalContext": "$context"
  }
}
EOF
}

json_allow() {
    local context="$1"
    cat << EOF
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "allow",
    "additionalContext": "$context"
  }
}
EOF
}

# ============================================================================
# Main validation
# ============================================================================

SEGMENTS=$(parse_chain "$REMOTE_CMD")
SEGMENT_COUNT=$(echo "$SEGMENTS" | grep -c . || echo 0)

# Check for dangerous patterns in any segment
DANGEROUS_SEGMENT=""
while IFS= read -r segment; do
    [ -z "$segment" ] && continue
    if is_dangerous "$segment"; then
        DANGEROUS_SEGMENT="$segment"
        break
    fi
done <<< "$SEGMENTS"

# ASK: Chained command with dangerous segment
if [ -n "$DANGEROUS_SEGMENT" ] && [ "$SEGMENT_COUNT" -gt 1 ]; then
    json_ask \
        "Dangerous pattern in chained command" \
        "WARNING: Chained command contains dangerous pattern: $DANGEROUS_SEGMENT. Recommend running dangerous commands separately. Proceed only if intentional."
    exit 0
fi

# ASK: Single dangerous command
if [ -n "$DANGEROUS_SEGMENT" ]; then
    json_ask \
        "Destructive operation: $DANGEROUS_SEGMENT" \
        "This command may cause data loss. Verify: 1) Target paths are correct 2) This is intentional 3) Backups exist if needed."
    exit 0
fi

# ============================================================================
# hpc-submit wrapper - always allow (wrapper handles validation)
# ============================================================================
if echo "$COMMAND" | grep -qE '^(hpc-submit|claude-sbatch)[[:space:]]|^~/.claude/plugins/hpc-sandbox/bin/(hpc-submit|claude-sbatch)'; then
    # Wrapper runs locally, handles its own validation flow
    exit 0
fi

# ============================================================================
# Raw sbatch via SSH - redirect to hpc-submit wrapper
# ============================================================================
if is_job_submission "$REMOTE_CMD"; then
    # Extract script path (POSIX-compatible, no grep -P)
    SCRIPT_PATH=$(echo "$REMOTE_CMD" | sed -n 's/.*sbatch[[:space:]]\{1,\}\([^[:space:]]\{1,\}\).*/\1/p')
    [ -z "$SCRIPT_PATH" ] && SCRIPT_PATH="<script>"

    json_deny \
        "Use hpc-submit wrapper instead" \
        "Raw sbatch is disabled. Use: ~/.claude/plugins/hpc-sandbox/bin/hpc-submit --script $SCRIPT_PATH --purpose \\\"...\\\" --outputs \\\"...\\\""
    exit 0
fi

# ASK: Job cancellation (scancel)
if is_job_cancellation "$REMOTE_CMD"; then
    # Extract job ID (POSIX-compatible, no grep -P)
    JOB_ID=$(echo "$REMOTE_CMD" | sed -n 's/.*scancel[[:space:]]\{1,\}\([^[:space:]]\{1,\}\).*/\1/p')
    [ -z "$JOB_ID" ] && JOB_ID="unknown"
    json_ask \
        "Job cancellation: $JOB_ID" \
        "Before cancelling: 1) Verify this is the correct job 2) Check job state (RUNNING jobs will lose progress) 3) Confirm you want to cancel."
    exit 0
fi

# ASK: Interactive session (srun)
if is_interactive_session "$REMOTE_CMD"; then
    json_ask \
        "Interactive session request" \
        "Before starting: 1) Review resource request (time, memory, partition) 2) Confirm this is needed."
    exit 0
fi

# ASK: Other destructive commands (rm, mv -f)
if is_destructive "$REMOTE_CMD"; then
    json_ask \
        "Potentially destructive: $REMOTE_CMD" \
        "This command may modify or delete files. Verify the paths and confirm."
    exit 0
fi

# ============================================================================
# Git operations via SSH
# ============================================================================

# Check each segment for git commands
while IFS= read -r segment; do
    [ -z "$segment" ] && continue

    if is_git_command "$segment"; then
        # ALLOW: Safe git operations (no JSON = pass through)
        if is_git_safe "$segment"; then
            # Safe read-only/staging operation, continue checking other segments
            continue
        fi

        # ASK: Dangerous git operations (with strong warning)
        if is_git_dangerous "$segment"; then
            json_ask \
                "Dangerous git operation" \
                "WARNING: $segment - This can cause permanent data loss. Force push overwrites remote history. Hard reset discards uncommitted changes. Proceed only if certain."
            exit 0
        fi

        # ASK: Commits and branch management (supervised)
        if is_git_supervised "$segment"; then
            json_ask \
                "Git operation requires confirmation" \
                "Supervised: $segment - Verify: 1) Correct branch/remote 2) Changes are intentional 3) Commit message is descriptive."
            exit 0
        fi

        # Unknown git command - ask to be safe
        json_ask \
            "Unrecognized git operation" \
            "Git command not in allowlist: $segment - Please verify this is safe."
        exit 0
    fi
done <<< "$SEGMENTS"

# ============================================================================
# All checks passed - inject docs on first command, then allow
# ============================================================================
if is_first_hpc_command; then
    mark_docs_shown
    # Escape the quick reference for JSON
    ESCAPED_DOCS=$(echo "$HPC_QUICK_REFERENCE" | python3 -c 'import sys,json; print(json.dumps(sys.stdin.read())[1:-1])')
    json_allow "$ESCAPED_DOCS"
    exit 0
fi

# Subsequent commands: no JSON = use permission system default
exit 0
