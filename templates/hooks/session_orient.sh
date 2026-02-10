#!/bin/bash
# SessionStart hook: inject sequence state and pending sync warnings
#
# On every session start:
# 1. Checks CRSP/project directory accessibility
# 2. Iterates .pending_sync.d/ for per-plan sync reminders
# 3. Checks for hook event divergence (global vs project settings)
# 4. Injects active sequence status
#
# Install: Add to .claude/settings.local.json under hooks.SessionStart

cd "$CLAUDE_PROJECT_DIR" 2>/dev/null || {
    echo '{"hookSpecificOutput":{"hookEventName":"SessionStart","additionalContext":"WARNING: Could not cd to project directory. Network mount may be stale."},"systemMessage":"Project directory unavailable"}'
    exit 0
}

# Check if plan CLI exists in this project
[ -f "./plan" ] || exit 0

# Check for pending sync breadcrumbs
PENDING_DIR="$HOME/.claude/plans/.pending_sync.d"
PENDING_MSGS=""

if [ -d "$PENDING_DIR" ]; then
    for breadcrumb in "$PENDING_DIR"/*; do
        [ -f "$breadcrumb" ] || continue
        PENDING_PATH=$(cat "$breadcrumb" 2>/dev/null)
        PLAN_NAME=$(basename "$breadcrumb")
        if [ -n "$PENDING_PATH" ]; then
            PENDING_MSGS="${PENDING_MSGS}PENDING SYNC [${PLAN_NAME}]: Run ./plan sync ${PENDING_PATH} BEFORE any other work.\n"
        fi
    done
fi

# Backward compat: legacy single-file breadcrumb
LEGACY_PENDING="$HOME/.claude/plans/.pending_sync"
if [ -f "$LEGACY_PENDING" ]; then
    LEGACY_PATH=$(cat "$LEGACY_PENDING" 2>/dev/null)
    if [ -n "$LEGACY_PATH" ]; then
        PENDING_MSGS="${PENDING_MSGS}PENDING SYNC [legacy]: Run ./plan sync ${LEGACY_PATH} BEFORE any other work.\n"
    fi
fi

# Check for hook event divergence between global and project settings
HOOK_WARN=""
GLOBAL_SETTINGS="$HOME/.claude/settings.json"
PROJECT_SETTINGS=".claude/settings.local.json"
if [ -f "$GLOBAL_SETTINGS" ] && [ -f "$PROJECT_SETTINGS" ]; then
    DIVERGENCE=$(python3 -c "
import json
try:
    with open('$GLOBAL_SETTINGS') as f:
        g = set(json.load(f).get('hooks', {}).keys())
    with open('$PROJECT_SETTINGS') as f:
        p = set(json.load(f).get('hooks', {}).keys())
    missing = g - p
    if missing:
        print('Hook events in global but NOT in project: ' + ', '.join(sorted(missing)) + '. Project hooks override global per-event.')
except:
    pass
" 2>/dev/null)
    [ -n "$DIVERGENCE" ] && HOOK_WARN="$DIVERGENCE"
fi

# Get sequence status
RESULT=$(python3 ./project/.plans/tools/cli.py sequence list 2>/dev/null)

if [ -n "$RESULT" ] || [ -n "$PENDING_MSGS" ] || [ -n "$HOOK_WARN" ]; then
    PENDING_MSGS="$PENDING_MSGS" HOOK_WARN="$HOOK_WARN" python3 -c "
import sys, json, os

seq_output = sys.stdin.read().strip()
pending_msgs = os.environ.get('PENDING_MSGS', '')
hook_warn = os.environ.get('HOOK_WARN', '')

parts = []
if hook_warn:
    parts.append('HOOK WARNING: ' + hook_warn)
if pending_msgs:
    parts.append('ACTION REQUIRED:\n' + pending_msgs.strip())
if seq_output:
    parts.append('ACTIVE SEQUENCES:\n' + seq_output)
parts.append('To orient: ./plan sequence status <name>')
parts.append('To work a plan: ./plan enter <sequence>/<plan>')

ctx = '\n\n'.join(parts)
sys_msg = 'Plan sync pending!' if pending_msgs else 'Sequence status loaded.'

print(json.dumps({
    'hookSpecificOutput': {
        'hookEventName': 'SessionStart',
        'additionalContext': ctx
    },
    'systemMessage': sys_msg
}))
" <<< "${RESULT:-}" 2>/dev/null
fi
exit 0
