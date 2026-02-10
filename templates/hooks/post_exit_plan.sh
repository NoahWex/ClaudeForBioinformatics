#!/bin/bash
# Post-ExitPlanMode hook
#
# After plan approval, writes a per-plan breadcrumb to .pending_sync.d/
# so sync instructions survive context reset across sessions.
#
# Breadcrumbs are named <area>__<name> and contain the draft file path.
# The session_orient hook reads these on next session start.
# ./plan sync removes the matching breadcrumb after successful sync.

# Consume stdin (Claude Code passes tool result as JSON)
INPUT=$(cat)

# Extract file path from ExitPlanMode tool response
FILE_PATH=$(echo "$INPUT" | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    print(data.get('tool_response', {}).get('filePath', ''))
except:
    print('')
" 2>/dev/null)

[[ -z "$FILE_PATH" ]] && exit 0

# Extract area, name, and draft_path from plan frontmatter
read -r AREA NAME DRAFT_PATH <<< "$(awk '
/^---$/ { if(++c==1) next; if(c==2) exit }
c==1 && /^area:/ { sub(/^area:[[:space:]]*/, ""); area=$0 }
c==1 && /^name:/ { sub(/^name:[[:space:]]*/, ""); name=$0 }
c==1 && /^draft_path:/ { sub(/^draft_path:[[:space:]]*/, ""); dp=$0 }
END { print area, name, dp }
' "$FILE_PATH" 2>/dev/null)"

# Use draft_path if available, otherwise FILE_PATH
SYNC_PATH="${DRAFT_PATH:-$FILE_PATH}"

# Write per-plan breadcrumb
PENDING_DIR="$HOME/.claude/plans/.pending_sync.d"
mkdir -p "$PENDING_DIR"

if [[ -n "$AREA" && -n "$NAME" ]]; then
    BREADCRUMB="$PENDING_DIR/${AREA}__${NAME}"
else
    BREADCRUMB="$PENDING_DIR/$(basename "$FILE_PATH" .md)"
fi

echo "$SYNC_PATH" > "$BREADCRUMB"

# Clean up legacy single-file breadcrumb
rm -f "$HOME/.claude/plans/.pending_sync" 2>/dev/null

cat << EOF
{
  "hookSpecificOutput": {
    "hookEventName": "PostToolUse",
    "additionalContext": "REQUIRED NEXT STEP: Run ./plan sync $SYNC_PATH to publish this plan. The plan will NOT be visible to other sessions until synced."
  },
  "systemMessage": "Plan approved. Agent should run: ./plan sync $SYNC_PATH"
}
EOF
