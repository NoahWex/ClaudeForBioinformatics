# Claude Code Installation

Install and configure the Claude Code CLI.

## Installation Options

### Via npm (Recommended)

```bash
npm install -g @anthropic-ai/claude-code
```

### Via Homebrew (macOS)

```bash
brew install claude-code
```

### Via Direct Download

Visit https://claude.ai/ and follow download instructions.

## Verify Installation

```bash
# Check version
claude --version

# Show help
claude --help
```

## Initial Authentication

First run will prompt for authentication:

```bash
claude

# Opens browser for Anthropic login
# Follow prompts to authenticate
```

## Directory Structure

Claude Code creates:

```
~/.claude/
├── settings.json          # Global settings
├── rules/                 # Global rules (optional)
├── plugins/               # Installed plugins
└── projects/              # Project session data
```

## Basic Configuration

Create global settings at `~/.claude/settings.json`:

```json
{
  "permissions": {
    "allow": [
      "Read",
      "Glob",
      "Grep"
    ],
    "deny": []
  }
}
```

This allows:
- Reading files
- Searching with glob patterns
- Searching with grep

Other tools prompt for permission by default.

## Test Basic Operation

```bash
# Start in a test directory
mkdir -p ~/claude-test
cd ~/claude-test

# Start Claude Code
claude

# Try: "What files are in this directory?"
# Try: "Create a simple Python hello world"
```

## Common Settings

### Allow Shell Commands

```json
{
  "permissions": {
    "allow": [
      "Read",
      "Glob",
      "Grep",
      "Bash(git status)",
      "Bash(git log *)",
      "Bash(ls *)"
    ]
  }
}
```

### Deny Patterns

```json
{
  "permissions": {
    "deny": [
      "Bash(rm -rf *)",
      "Bash(sudo *)"
    ]
  }
}
```

## Verify Complete

```bash
# 1. Version check
claude --version
# Expected: claude-code X.Y.Z

# 2. Help available
claude --help
# Expected: Usage information

# 3. Can start session
cd ~/claude-test
claude
# Expected: Claude prompt appears
# Type: "exit" to leave
```

## Next Step

→ [03_plugin_install.md](03_plugin_install.md) - Install HPC sandbox plugin
