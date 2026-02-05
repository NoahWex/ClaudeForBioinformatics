# 03: Install HPC Toolkit

Install the `hpc-toolkit` — a unified CLI for HPC operations with safety guardrails, git workspace acceleration, and audit logging.

## Prerequisites

- SSH access to HPC (completed in [01_ssh_config.md](01_ssh_config.md))
- Claude Code installed (completed in [02_claude_code_install.md](02_claude_code_install.md))

## What Gets Installed

```
~/.claude/
├── hpc-toolkit/
│   ├── bin/hpc              # Unified CLI (submit, shell, status, logs, cancel, file, git, workspace)
│   ├── hooks/gate.sh        # PreToolUse hook — blocks unsafe commands, injects docs
│   ├── docs/CLAUDE_GUIDE.md # Quick reference injected on first HPC command per session
│   ├── config.sh            # Your environment configuration
│   └── logs/                # Audit trail (auto-created)
├── skills/hpc/SKILL.md      # HPC skill for auto-invocation
└── rules/
    ├── trace-framework.md    # Objectivity and scope discipline
    └── hpc-safety-guidelines.md  # HPC operation safety rules
```

## Steps

### 1. Copy Toolkit Files

From the sandbox directory:

```bash
# Create target directories
mkdir -p ~/.claude/hpc-toolkit/bin
mkdir -p ~/.claude/hpc-toolkit/hooks
mkdir -p ~/.claude/hpc-toolkit/docs
mkdir -p ~/.claude/hpc-toolkit/logs
mkdir -p ~/.claude/skills/hpc
mkdir -p ~/.claude/rules

# Copy toolkit
cp global-setup/hpc-toolkit/bin/hpc       ~/.claude/hpc-toolkit/bin/
cp global-setup/hpc-toolkit/hooks/gate.sh ~/.claude/hpc-toolkit/hooks/
cp global-setup/hpc-toolkit/docs/CLAUDE_GUIDE.md ~/.claude/hpc-toolkit/docs/

# Copy skill and rules
cp global-setup/skills/hpc/SKILL.md       ~/.claude/skills/hpc/
cp global-setup/rules/trace-framework.md  ~/.claude/rules/
cp global-setup/rules/hpc-safety-guidelines.md ~/.claude/rules/
```

### 2. Make Scripts Executable

```bash
chmod +x ~/.claude/hpc-toolkit/bin/hpc
chmod +x ~/.claude/hpc-toolkit/hooks/gate.sh
```

### 3. Create Your Configuration

```bash
cp global-setup/hpc-toolkit/config.sh.template ~/.claude/hpc-toolkit/config.sh
```

Edit `~/.claude/hpc-toolkit/config.sh`:

```bash
# Required: Your SSH alias (must match ~/.ssh/config Host entry)
SSH_ALIAS="hpc3"

# Required: CRSP mount paths
# Find your macOS mount path:
#   ls ~/Library/Application\ Support/CRSP\ Desktop/Volumes.noindex/
CRSP_MAC_PREFIX="$HOME/Library/Application Support/CRSP Desktop/Volumes.noindex/CRSP Lab - dalawson.localized"
CRSP_HPC_PREFIX="/share/crsp/lab/dalawson"

# Optional: Your HPC username (blank = auto-detect)
HPC_USER=""

# Optional: Default repo path for hpc git fallback
FALLBACK_REPO=""
```

### 4. Configure Claude Code Settings

Create or update `~/.claude/settings.json`:

```json
{
  "permissions": {
    "allow": [
      "Bash(~/.claude/hpc-toolkit/bin/hpc *)",
      "Bash(ssh hpc3 *)"
    ]
  },
  "hooks": {
    "SessionStart": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "cat ~/.claude/rules/trace-framework.md 2>/dev/null || true"
          }
        ]
      }
    ],
    "PreToolUse": [
      {
        "matcher": "Bash",
        "hooks": [
          {
            "type": "command",
            "command": "~/.claude/hpc-toolkit/hooks/gate.sh",
            "timeout": 5
          }
        ]
      }
    ]
  }
}
```

**Note:** Replace `hpc3` in the `ssh hpc3 *` permission with your `SSH_ALIAS` value if different.

If you already have a `settings.json`, merge these entries into your existing config.

### 5. Verify Installation

Continue to [04_verify_setup.md](04_verify_setup.md).

## How the Toolkit Works

### Gate Hook (gate.sh)

Runs before every Bash command Claude attempts:
- **First HPC command** in a session → injects CLAUDE_GUIDE.md as context
- **Toolkit commands** (`hpc *`) → allowed through
- **Read-only queries** (`squeue`, `sacct`, `sinfo`) → allowed through
- **Raw sbatch, python, git writes** → **blocked** with suggested toolkit command
- **Local git writes on CRSP mount** → **blocked** with redirect to `hpc git`

### Audit Log

All `hpc` commands are logged to `~/.claude/hpc-toolkit/logs/audit.jsonl` with timestamps and purpose.

### Git Workspace Acceleration (Optional)

`hpc workspace init <crsp-repo-path>` creates a bare git mirror on fast DFS storage. After init, `hpc git` auto-uses the mirror. See `10_HPC_PLUGIN.md` for details.

## Troubleshooting

### Hook not triggering

```bash
# Verify settings.json syntax
python3 -c "import json; json.load(open('$HOME/.claude/settings.json'))"

# Check hook script is executable
ls -la ~/.claude/hpc-toolkit/hooks/gate.sh

# Start a new Claude Code session (hooks load at session start)
```

### "Command not found" errors

```bash
# Check bin/hpc exists and is executable
ls -la ~/.claude/hpc-toolkit/bin/hpc

# Test directly
~/.claude/hpc-toolkit/bin/hpc --help
```

### Config not loading

```bash
# Check config.sh exists
cat ~/.claude/hpc-toolkit/config.sh

# Test sourcing it
source ~/.claude/hpc-toolkit/config.sh && echo "SSH_ALIAS=$SSH_ALIAS"
```

## Next Step

→ [04_verify_setup.md](04_verify_setup.md) - Verify complete setup
