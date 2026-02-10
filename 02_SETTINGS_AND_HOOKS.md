# Claude Code Settings and Hooks

Configure Claude Code's behavior with settings files and session hooks.

## File Locations

| File | Scope | Purpose |
|------|-------|---------|
| `.claude/settings.json` | Project (committed) | Shared project settings |
| `.claude/settings.local.json` | Local (gitignored) | User-specific settings |
| `~/.claude/settings.json` | Global | User defaults |

```mermaid
graph TD
    G["~/.claude/settings.json<br/>(Global defaults)"] --> P[".claude/settings.json<br/>(Project — committed)"]
    P --> L[".claude/settings.local.json<br/>(Local — gitignored)"]
    G -->|"overridden by"| P
    P -->|"overridden by"| L
    style G fill:#d0e8f9,stroke:#007
    style P fill:#f9ecd0,stroke:#960
    style L fill:#d0f9d4,stroke:#070
```

## Settings Structure

### Basic Template (`.claude/settings.json`)

```json
{
  "permissions": {
    "allow": []
  }
}
```

This minimal version commits no permissions, leaving them for local configuration.

### Full Local Settings (`.claude/settings.local.json`)

```json
{
  "permissions": {
    "allow": [
      "Bash(ls:*)",
      "Bash(mkdir:*)",
      "Bash(cp:*)",
      "Bash(mv:*)",
      "Bash(chmod:*)",
      "Bash(find:*)",
      "Bash(grep:*)",
      "Bash(python3:*)",
      "Bash(Rscript:*)",
      "Bash(git:*)",
      "WebSearch",
      "WebFetch(domain:github.com)",
      "WebFetch(domain:docs.example.com)",
      "Skill(project-status)",
      "Skill(project-discovery)"
    ],
    "deny": [],
    "additionalDirectories": [
      "/path/to/reference/codebase/",
      "/path/to/shared/templates/"
    ]
  },
  "hooks": {
    "SessionStart": [
      {
        "hooks": [
          {
            "type": "command",
            "command": "cat .claude/trace-framework.md"
          }
        ]
      }
    ]
  },
  "enabledPlugins": {
    "r-lsp": true,
    "pyright-lsp@claude-plugins-official": true
  }
}
```

## Permission Patterns

### Permission Types

| Type | Behavior |
|------|----------|
| `allow` | Auto-approve, no prompt |
| `ask` | Yes/No prompt every time (no "always allow" option) |
| `deny` | Block entirely |

```json
{
  "permissions": {
    "allow": ["Bash(ls:*)"],     // Auto-approved
    "ask": ["Bash(rm -rf:*)"],   // Always prompts yes/no
    "deny": ["Bash(shutdown:*)"] // Blocked
  }
}
```

**Use `ask` for destructive operations** you want to review but not block entirely.

### Command Permissions

```json
"Bash(command:*)"           // Allow command with any args
"Bash(git check-ignore:*)"  // Specific subcommand
```

### Domain-Specific Web Access

```json
"WebFetch(domain:docs.scvi-tools.org)",
"WebFetch(domain:satijalab.github.io)",
"WebFetch(domain:www.nature.com)"
```

### Skill Permissions

```json
"Skill(project-status)",
"Skill(marker-debate)"
```

## Additional Directories

Grant Claude access to directories outside the project:

```json
"additionalDirectories": [
  "/path/to/reference/analysis/code/",
  "/path/to/shared/templates/"
]
```

Use cases:
- Reference codebases for patterns
- Shared lab templates
- External data directories

## Session Hooks

```mermaid
sequenceDiagram
    participant CC as Claude Code
    participant H as gate.sh hook
    participant CMD as Command
    CC->>H: PreToolUse: Bash("ssh hpc3 sbatch ...")
    H->>H: Pattern match: raw sbatch
    H-->>CC: ❌ BLOCK — use hpc submit instead
    CC->>H: PreToolUse: Bash("hpc submit ...")
    H->>H: Pattern match: toolkit command
    H-->>CMD: ✅ PASS — execute
    CMD-->>CC: Result
```

### SessionStart Hook

Execute commands when a session begins:

```json
"hooks": {
  "SessionStart": [
    {
      "hooks": [
        {
          "type": "command",
          "command": "cat .claude/trace-framework.md"
        }
      ]
    }
  ]
}
```

This injects the trace-framework content at session start, ensuring objectivity rules are loaded before any interaction.

### Hook Types

| Type | Description |
|------|-------------|
| `command` | Run shell command, output injected into context |

### Hook Events

| Event | When | Can Block? |
|-------|------|------------|
| `SessionStart` | New session begins | No |
| `PreToolUse` | Before a tool executes | Yes (`permissionDecision: deny`) |
| `PostToolUse` | After a tool executes | No |
| `PreCompact` | Before context compaction | No (observe only) |

### PreToolUse Matchers

Restrict when a hook fires using the `matcher` field:

```json
{
  "matcher": "ExitPlanMode",
  "hooks": [{"type": "command", "command": "~/.claude/hooks/pre_exit_plan.sh"}]
}
```

Matchers support pipe-delimited patterns:

```json
{
  "matcher": "Edit|Write",
  "hooks": [{"type": "command", "command": "~/.claude/hooks/plan_md_reminder.sh"}]
}
```

### Hook Output Format

Hooks communicate via JSON on stdout:

```json
{
  "hookSpecificOutput": {
    "hookEventName": "PreToolUse",
    "permissionDecision": "allow",
    "additionalContext": "Message visible to agent (not user)"
  },
  "systemMessage": "Message visible to user terminal (not agent)"
}
```

- `additionalContext` → injected into agent context (the agent sees this)
- `systemMessage` → shown in terminal (the user sees this)
- `permissionDecision` → `"allow"` or `"deny"` (PreToolUse only)

### Hook Event Override Behavior

**Project-level hooks override global hooks per event.** If `.claude/settings.local.json`
defines `SessionStart`, global `~/.claude/settings.json` hooks for `SessionStart` do NOT fire.
Other events (e.g., `PreToolUse`) from global settings still fire unless also overridden.

This means: if you add a project-level `SessionStart` hook, you must duplicate any global
`SessionStart` hooks you still want to run. The session orient hook (see below) warns
when global events are missing from the project config.

### Plan Coordination Hooks

For projects using the sequence coordination system, five hooks form the
draft→validate→sync pipeline. See [14_SEQUENCE_COORDINATION.md](14_SEQUENCE_COORDINATION.md)
for the full pipeline design.

Template scripts are provided in `templates/hooks/`:

| Hook | Event | Matcher | Purpose |
|------|-------|---------|---------|
| `pre_exit_plan.sh` | PreToolUse | ExitPlanMode | Validate plan structure before approval |
| `post_exit_plan.sh` | PostToolUse | ExitPlanMode | Write sync breadcrumb for next session |
| `session_orient.sh` | SessionStart | — | Inject sequence status + pending sync warnings |
| `pre_compact_plan.sh` | PreCompact | — | Preserve plan state through context compaction |
| `plan_md_reminder.sh` | PreToolUse | Edit\|Write | Remind about plan structure in plan mode |

Full settings configuration example:

```json
{
  "hooks": {
    "PreToolUse": [
      {
        "matcher": "ExitPlanMode",
        "hooks": [{"type": "command", "command": "~/.claude/hooks/pre_exit_plan.sh", "timeout": 10}]
      },
      {
        "matcher": "Edit|Write",
        "hooks": [{"type": "command", "command": "~/.claude/hooks/plan_md_reminder.sh", "timeout": 3}]
      }
    ],
    "PostToolUse": [
      {
        "matcher": "ExitPlanMode",
        "hooks": [{"type": "command", "command": "~/.claude/hooks/post_exit_plan.sh", "timeout": 5}]
      }
    ],
    "SessionStart": [
      {
        "hooks": [{"type": "command", "command": "~/.claude/hooks/session_orient.sh", "timeout": 10}]
      }
    ],
    "PreCompact": [
      {
        "hooks": [{"type": "command", "command": "~/.claude/hooks/pre_compact_plan.sh", "timeout": 5}]
      }
    ]
  }
}
```

## Language Server Plugins

Enable language-specific intelligence:

```json
"enabledPlugins": {
  "r-lsp": true,
  "pyright-lsp@claude-plugins-official": true,
  "typescript-lsp@claude-plugins-official": true
}
```

## Best Practices

### 1. Separate Project vs Local Settings

**Project settings** (committed):
- Minimal permissions
- Project-wide configurations
- Plugin enablement

**Local settings** (gitignored):
- User-specific permissions
- Additional directories
- Session hooks

### 2. Progressive Permission Granting

Start restrictive, add permissions as needed:

```json
// Start with nothing
"allow": []

// Add as you use them
"allow": [
  "Bash(ls:*)",  // Added when needed listing
  "Bash(git:*)" // Added when needed git operations
]
```

### 3. Domain-Specific Web Access

Grant only documentation sites you actually use:

```json
// Good: Specific domains
"WebFetch(domain:docs.scvi-tools.org)"

// Avoid: Overly broad
"WebFetch(*)"
```

### 4. Reference Directories

Include mature codebases as references:

```json
"additionalDirectories": [
  "../reference_project/"  // Working patterns to copy
]
```

## Common Permission Patterns

### Scientific Computing

```json
"allow": [
  "Bash(python3:*)",
  "Bash(Rscript:*)",
  "Bash(R --slave -e:*)",
  "Bash(conda:*)",
  "Bash(mamba:*)"
]
```

### HPC Operations (with Toolkit)

Split read-only commands (auto-allow) from destructive operations (always prompt):

```json
"allow": [
  "Bash(~/.claude/hpc-toolkit/bin/hpc status*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc logs*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc file ls*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc file cat*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc git status*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc git log*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc git diff*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc workspace list*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc --help*)"
],
"ask": [
  "Bash(~/.claude/hpc-toolkit/bin/hpc submit*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc cancel*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc git push*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc git commit*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc file rm*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc file mv*)",
  "Bash(~/.claude/hpc-toolkit/bin/hpc workspace sync*)"
]
```

This pattern ensures destructive HPC operations always prompt yes/no without the option to auto-approve permanently.

### Git Operations

```json
"allow": [
  "Bash(git status:*)",
  "Bash(git diff:*)",
  "Bash(git log:*)",
  "Bash(git check-ignore:*)"
]
```

### File Operations

```json
"allow": [
  "Bash(ls:*)",
  "Bash(mkdir:*)",
  "Bash(cp:*)",
  "Bash(mv:*)",
  "Bash(chmod:*)",
  "Bash(touch:*)",
  "Bash(wc:*)",
  "Bash(unzip:*)"
]
```

## Troubleshooting

### Permission Denied Errors

If a command is blocked:
1. Check exact command syntax in error
2. Add specific permission pattern
3. Use wildcards cautiously

### Hook Not Executing

Verify:
1. File exists at specified path
2. Command syntax is valid
3. JSON is properly formatted
4. **Check for event override**: If your project defines `hooks` in `settings.local.json`,
   global hooks for the same event names are silently suppressed. Add missing events to
   the project config.

### Plugin Not Working

Check:
1. Plugin is installed
2. Correct plugin identifier used
3. `enabledPlugins` key spelled correctly
