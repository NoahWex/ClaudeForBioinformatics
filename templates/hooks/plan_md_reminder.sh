#!/bin/bash
# PreToolUse hook: remind about plan structure when editing .md in plan mode
#
# Matcher: Edit|Write (only fires for file editing tools)
# Only activates when permission_mode == "plan" and target is .md
#
# Install: Add to .claude/settings.local.json under hooks.PreToolUse
# with matcher: "Edit|Write"

input=$(cat)

mode=$(echo "$input" | jq -r '.permission_mode // empty')
file_path=$(echo "$input" | jq -r '.tool_input.file_path // empty')

if [[ "$mode" == "plan" ]] && [[ "$file_path" == *.md ]]; then
    # Point to your project's interface-aware planning rules
    # Customize this path for your project:
    RULES_PATH="$HOME/.claude/rules/planning/interface-aware-planning.md"
    if [[ -f "$RULES_PATH" ]]; then
        cat "$RULES_PATH"
    fi
fi
