#!/bin/bash
# PreCompact hook: inject plan state before context compaction
#
# When context is compacted, the agent loses detailed state.
# This hook re-injects a summary so compacted context retains orientation.
#
# Install: Add to .claude/settings.local.json under hooks.PreCompact

# Consume stdin
INPUT=$(cat)

# Find project directory
PROJECT_DIR="${CLAUDE_PROJECT_DIR:-$(pwd)}"
cd "$PROJECT_DIR" 2>/dev/null || exit 0

# Check if plan CLI exists
[ -f "./plan" ] || exit 0

# Gather plan state
STATE=$(python3 -c "
import sys, json
sys.path.insert(0, 'project/.plans/tools')

try:
    from sequence import SequenceManager
    from pathlib import Path

    plans_dir = Path('project/.plans')
    mgr = SequenceManager(plans_dir)
    mgr.load_all()

    parts = []

    # Current session
    session_file = plans_dir / '.current_session'
    if session_file.exists():
        lines = session_file.read_text().strip().split('\n')
        if len(lines) >= 3:
            parts.append(f'ACTIVE SESSION: {lines[0]}/{lines[1]} (session {lines[2]})')

    # Sequence progress
    for name, seq in sorted(mgr.sequences.items()):
        if seq.status != 'active':
            continue
        completed = sum(1 for p in seq.plans if p.status == 'complete')
        total = len(seq.plans)
        ready = [p.short_name for p in seq.plans if p.status == 'queued']
        in_prog = [p.short_name for p in seq.plans if p.status == 'in_progress']
        line = f'{name}: {completed}/{total}'
        if in_prog:
            line += f' (active: {\", \".join(in_prog)})'
        if ready:
            line += f' (ready: {\", \".join(ready[:3])})'
        parts.append(line)

    # Pending syncs
    import os
    pending_dir = os.path.expanduser('~/.claude/plans/.pending_sync.d')
    if os.path.isdir(pending_dir):
        for f in os.listdir(pending_dir):
            path = open(os.path.join(pending_dir, f)).read().strip()
            parts.append(f'PENDING SYNC: {f} -> {path}')

    if parts:
        print(json.dumps({
            'hookSpecificOutput': {
                'hookEventName': 'PreCompact',
                'additionalContext': 'PLAN STATE (pre-compaction summary):\n' + '\n'.join(parts)
            }
        }))
except Exception:
    pass
" 2>/dev/null)

if [ -n "$STATE" ]; then
    echo "$STATE"
fi
