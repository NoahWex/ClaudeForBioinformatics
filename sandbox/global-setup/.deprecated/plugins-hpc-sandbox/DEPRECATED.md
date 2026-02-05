# DEPRECATED: hpc-sandbox Plugin

**Superseded by:** `hpc-toolkit/` (the unified CLI at `~/.claude/hpc-toolkit/`)

**Date:** 2026-02-02

## What Changed

The original `hpc-sandbox` plugin used:
- `bin/hpc-submit` — single-purpose job submission wrapper
- `hooks/validate-hpc-command.sh` — PreToolUse hook for command validation
- `skills/hpc-operations/` — skill with embedded gotchas and practices
- `commands/hpc.md`, `hpc-status.md` — slash command definitions

The replacement `hpc-toolkit` provides:
- `bin/hpc` — unified CLI with 8 subcommands (submit, shell, status, logs, cancel, file, git, workspace)
- `hooks/gate.sh` — streamlined PreToolUse gate with session-level doc injection
- Git workspace acceleration via DFS mirrors
- Configurable via `config.sh` (no hardcoded paths)

## De-deprecation

If you need to restore this plugin, copy contents back to `~/.claude/plugins/hpc-sandbox/` and update `~/.claude/plugins/installed_plugins.json`. Note that the plugin API and the toolkit approach are architecturally different — they should not be used simultaneously.
