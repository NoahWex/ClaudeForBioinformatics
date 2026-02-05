#!/bin/bash --noprofile
# gate.sh - Validation gate for HPC commands
#
# - First HPC command: inject docs
# - hpc toolkit commands: allow
# - Read-only queries (squeue, sacct, sinfo, module avail): allow
# - Everything else: redirect to toolkit

set -euo pipefail

TOOLKIT_DIR="$HOME/.claude/hpc-toolkit"
SESSION_STATE_DIR="/tmp/.claude-hpc-session"
SESSION_ID="${PPID:-$(date +%Y%m%d)}"
SESSION_STATE_FILE="$SESSION_STATE_DIR/docs-shown-$SESSION_ID"

# Load user config for SSH_ALIAS
CONFIG_FILE="$TOOLKIT_DIR/config.sh"
if [[ -f "$CONFIG_FILE" ]]; then
    # shellcheck source=/dev/null
    source "$CONFIG_FILE"
fi
HPC_ALIAS="${SSH_ALIAS:-hpc3}"

mkdir -p "$SESSION_STATE_DIR" 2>/dev/null || true

get_docs() {
    cat "$TOOLKIT_DIR/docs/CLAUDE_GUIDE.md" 2>/dev/null || echo "HPC toolkit docs not found"
}

json_escape() {
    python3 -c 'import sys,json; print(json.dumps(sys.stdin.read())[1:-1])'
}

json_deny() {
    cat << EOF
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "deny",
    "permissionDecisionReason": "$1",
    "additionalContext": "$2"
  }
}
EOF
}

json_allow() {
    cat << EOF
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "allow",
    "additionalContext": "$1"
  }
}
EOF
}

is_hpc_ssh() {
    echo "$1" | grep -qE "^ssh[[:space:]]+${HPC_ALIAS}[[:space:]]+"
}

is_hpc_command() {
    echo "$1" | grep -qE "(ssh[[:space:]]+${HPC_ALIAS})|^hpc[[:space:]]"
}

is_toolkit_command() {
    echo "$1" | grep -qE '^~/.claude/hpc-toolkit/bin/hpc[[:space:]]|/hpc-toolkit/bin/hpc[[:space:]]'
}

is_allowed_query() {
    local cmd="$1"
    # Read-only job queries
    echo "$cmd" | grep -qE '^(squeue|sacct|sinfo|scontrol[[:space:]]+show|module[[:space:]]+avail)[[:space:]]?' && return 0
    return 1
}

extract_remote_cmd() {
    echo "$1" | sed -E "s/^ssh[[:space:]]+${HPC_ALIAS}[[:space:]]+[\"']?//" | sed -E "s/[\"']?\$//"
}

# Main
INPUT=$(cat)
COMMAND=$(echo "$INPUT" | python3 -c "import sys, json; d=json.load(sys.stdin); print((d.get('tool_input') or {}).get('command', ''))" 2>/dev/null || echo "")

[ -z "$COMMAND" ] && exit 0

# Block local git write commands on CRSP-mounted repos
# CRSP's macOS SMB mount caches .git/objects/ directory listings aggressively.
# Git writes from the local mount create objects invisible to HPC, causing
# "bad object HEAD", stale refs, and index.lock conflicts.
is_local_git_write() {
    local cmd="$1"
    # Match: git <subcmd> or git -C <path> <subcmd>
    echo "$cmd" | grep -qE '^git[[:space:]]' || return 1

    # Extract the git subcommand (skip flags like -C, --no-pager)
    local subcmd
    subcmd=$(echo "$cmd" | sed -E 's/^git[[:space:]]+(-C[[:space:]]+"[^"]*"[[:space:]]+|-C[[:space:]]+[^[:space:]]+[[:space:]]+|--no-pager[[:space:]]+|--git-dir=[^[:space:]]+[[:space:]]+)*//' | awk '{print $1}')

    # Read-only commands are safe locally
    case "$subcmd" in
        status|log|diff|show|cat-file|rev-parse|for-each-ref|ls-files|ls-tree|check-ignore|count-objects|describe|name-rev|shortlog|whatchanged|blame|grep)
            return 1 ;;
    esac

    # Everything else is a write — block it
    return 0
}

if is_local_git_write "$COMMAND"; then
    # Extract git args for the suggestion
    GIT_ARGS=$(echo "$COMMAND" | sed -E 's/^git[[:space:]]+//')
    json_deny \
        "Git writes must use hpc git (CRSP cache coherency)" \
        "CRSP's macOS mount caches .git/objects/ aggressively. Local git writes cause desync with HPC.\n\nTry: ~/.claude/hpc-toolkit/bin/hpc git $GIT_ARGS"
    exit 0
fi

is_hpc_command "$COMMAND" || is_toolkit_command "$COMMAND" || exit 0

# All toolkit commands pass through to Claude Code's permission system
# - Read-only ops (status, logs, file cat, git log, etc.) are in settings.json allow list → no prompt
# - Destructive ops (submit, git push, file rm, etc.) are in settings.json ask list → yes/no prompt
is_toolkit_command "$COMMAND" && exit 0

# Helper: inject docs on first HPC command if not yet shown
maybe_inject_docs() {
    if [ ! -f "$SESSION_STATE_FILE" ]; then
        touch "$SESSION_STATE_FILE" 2>/dev/null || true
        DOCS=$(get_docs | json_escape)
        echo "=== HPC TOOLKIT GUIDE ===\n$DOCS\n=== END GUIDE ==="
    fi
}

# SSH commands to HPC
if is_hpc_ssh "$COMMAND"; then
    REMOTE=$(extract_remote_cmd "$COMMAND")

    # Allow read-only queries (with docs on first use)
    if is_allowed_query "$REMOTE"; then
        DOCS_MSG=$(maybe_inject_docs)
        if [ -n "$DOCS_MSG" ]; then
            json_allow "$DOCS_MSG"
        fi
        exit 0
    fi

    # Everything else is denied — redirect to toolkit
    # Mark session and inject docs with the deny message
    touch "$SESSION_STATE_FILE" 2>/dev/null || true
    DOCS=$(get_docs | json_escape)

    HPC="~/.claude/hpc-toolkit/bin/hpc"
    SUGGESTION="Use $HPC"
    if echo "$REMOTE" | grep -qE '^sbatch[[:space:]]'; then
        SCRIPT=$(echo "$REMOTE" | sed -E 's/^sbatch[[:space:]]+//' | tr ' ' '\n' | grep -v '^-' | head -1)
        SUGGESTION="$HPC submit $SCRIPT --purpose \\\"...\\\" --outputs \\\"...\\\""
    elif echo "$REMOTE" | grep -qE '^srun[[:space:]]'; then
        SUGGESTION="$HPC shell --cmd \\\"...\\\" --purpose \\\"...\\\" --time 1:00:00"
    elif echo "$REMOTE" | grep -qE '^(python|Rscript|bash)[[:space:]]'; then
        SUGGESTION="$HPC shell --cmd \\\"$REMOTE\\\" --purpose \\\"...\\\" --time 1:00:00"
    elif echo "$REMOTE" | grep -qE '^(cd[[:space:]]+.*&&[[:space:]]*)?git[[:space:]]'; then
        GIT_ARGS=$(echo "$REMOTE" | sed -E 's/^(cd[[:space:]]+.*&&[[:space:]]*)?git[[:space:]]+//')
        SUGGESTION="$HPC git $GIT_ARGS"
    elif echo "$REMOTE" | grep -qE '^(ls|cat|head|tail|rm|cp|mv|mkdir)[[:space:]]'; then
        SUGGESTION="$HPC file ${REMOTE%% *} ..."
    elif echo "$REMOTE" | grep -qE '^scancel[[:space:]]'; then
        JOB_ID=$(echo "$REMOTE" | sed -E 's/^scancel[[:space:]]+//' | head -1)
        SUGGESTION="$HPC cancel $JOB_ID"
    fi

    json_deny \
        "Use hpc toolkit" \
        "Try: $SUGGESTION\n\n=== HPC TOOLKIT GUIDE ===\n$DOCS\n=== END GUIDE ==="
    exit 0
fi

exit 0
