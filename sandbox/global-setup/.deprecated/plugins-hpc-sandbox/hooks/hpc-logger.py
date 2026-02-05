#!/usr/bin/env python3
"""
hpc-logger.py - PostToolUse logging for HPC commands
Logs all HPC-related commands to ~/.claude/hpc/logs/commands.jsonl
"""

import json
import sys
import os
import re
from datetime import datetime
from pathlib import Path

LOG_DIR = Path.home() / ".claude" / "hpc" / "logs"
COMMANDS_LOG = LOG_DIR / "commands.jsonl"

def is_hpc_command(command: str) -> bool:
    """Check if command is HPC-related."""
    patterns = [
        r"ssh.*hpc",
        r"^s(batch|run|queue|cancel|info|acct)\b",
        r"^module\s",
        r"^conda\s+activate",
    ]
    return any(re.search(p, command) for p in patterns)

def log_command(entry: dict):
    """Append entry to commands log."""
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    with open(COMMANDS_LOG, "a") as f:
        f.write(json.dumps(entry) + "\n")

def main():
    try:
        input_data = json.loads(sys.stdin.read())
    except json.JSONDecodeError:
        sys.exit(0)

    tool_name = input_data.get("tool_name", "")
    tool_input = input_data.get("tool_input", {})
    tool_result = input_data.get("tool_result", {})

    # Only log Bash commands
    if tool_name != "Bash":
        sys.exit(0)

    command = tool_input.get("command", "")

    # Only log HPC-related commands
    if not is_hpc_command(command):
        sys.exit(0)

    # Create log entry
    entry = {
        "timestamp": datetime.now().isoformat(),
        "command": command,
        "exit_code": tool_result.get("exit_code"),
        "success": tool_result.get("exit_code", 1) == 0,
    }

    log_command(entry)
    sys.exit(0)

if __name__ == "__main__":
    main()
