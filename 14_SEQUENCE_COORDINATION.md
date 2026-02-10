# Sequence Coordination System

Coordinate multi-plan workflows with dependency tracking, vision persistence, and session management.

## Problem

Complex work areas (pipelines, multi-study analyses) span many plans. Individual plans track tasks, but nothing tracks the overall goal, plan ordering, or session state across plans.

Without sequences:
- Vision gets lost as you move between plans
- Dependencies between plans are implicit
- No way to know which plans are ready to work on
- Session handoffs lose context

## Concepts

### Sequence

A YAML file that defines an ordered set of plans with:
- **Vision**: The end goal and why it matters
- **Plans**: What needs to happen, with dependencies
- **State**: What's done, what's in progress, what's blocked

### Vision

The "why" that stays visible throughout work:

```yaml
vision:
  end_goal: "What success looks like"
  why_it_matters: "Connection to larger project"
  current_state: "Where things stand right now"
  decisions:
    - date: "2026-02-05"
      fork: "What was being decided"
      chose: "What was chosen"
      rationale: "Why, with tradeoffs"
  exit_criteria:
    - "Verifiable criterion 1"
    - "Verifiable criterion 2"
```

`end_goal` and `exit_criteria` are required. Everything else is optional but recommended.

### Session Ownership

Plans can be claimed by sessions to prevent conflicts:
- `plan enter` claims a plan and records a session ID
- `plan pause` saves progress notes and releases the claim
- `plan seq-complete` marks done and unblocks downstream plans

## Setup

### Directory Structure

```
.plans/
├── sequences/
│   └── my_pipeline.yaml      # Sequence definitions
├── active_plans/
│   └── area/
│       └── plan_name.md      # Plan files (created on demand)
├── .current_session           # Active session tracking (gitignored)
└── .gitignore                 # Should include .current_session
```

### Creating a Sequence

1. Copy the template:
   ```bash
   cp graph-memory/templates/sequence_template.yaml .plans/sequences/my_pipeline.yaml
   ```

2. Edit the YAML:
   - Define the vision (end goal, exit criteria)
   - List plans in dependency order
   - Set `depends_on` for each plan (using short names)

3. Validate:
   ```bash
   ./plan sequence validate my_pipeline
   ```

## Workflow

### 1. Check Status

```bash
# Quick overview of all sequences
./plan sequence list

# Full status with vision for one sequence
./plan sequence status my_pipeline

# Vision only (re-orient quickly)
./plan sequence vision my_pipeline
```

### 2. Enter a Plan

```bash
./plan enter my_pipeline/plan_name
```

This:
- Claims the plan (sets ownership)
- Creates a `.md` file if one doesn't exist
- Shows the vision reminder and prior notes
- Saves session info for `pause`/`seq-complete`

### 3. Work the Plan

Work normally within the plan. The plan `.md` file tracks tasks, discussion, and state.

### 4. Pause or Complete

**Pause** (save state for later):
```bash
./plan pause "Completed steps 1-3, step 4 needs input data from upstream"
```

**Complete** (mark done, unblock downstream):
```bash
./plan seq-complete "Built all Milo objects, validated cell counts"
```

On completion:
- Plan status changes to `complete`
- Downstream plans with all deps met become `queued`
- If `auto_progress: prompt`, you're asked to enter the next ready plan

## Sequence YAML Schema

```yaml
# Required fields
name: string              # Unique identifier
description: string       # Brief description
status: active|paused|complete
auto_progress: prompt|silent|disabled
created: "YYYY-MM-DD"
updated: "YYYY-MM-DD"

# Vision (required)
vision:
  end_goal: string        # REQUIRED
  why_it_matters: string  # Recommended
  current_state: string   # Updated as work progresses
  decisions: []           # Dated decision records
  exit_criteria: []       # REQUIRED, verifiable items

# Plans (required, at least one)
plans:
  - id: area/name         # Full plan ID
    summary: string       # One-line description
    status: queued|in_progress|complete|blocked
    depends_on: []        # Short names of upstream plans
    # Set during work:
    owner: string         # Session ID (set by enter)
    notes: string         # Progress notes (set by pause)
    blocked_reason: string
    completed_summary: string
```

### Plan Dependencies

Dependencies use **short names** (the part after the last `/`):

```yaml
plans:
  - id: analysis/data_prep
    depends_on: []

  - id: analysis/model_fit
    depends_on:
      - data_prep          # Short name of analysis/data_prep

  - id: analysis/validation
    depends_on:
      - model_fit
      - data_prep          # Can depend on multiple plans
```

Circular dependencies are detected by `plan sequence validate`.

### Auto-Progress

Controls what happens after `plan seq-complete`:

| Value | Behavior |
|-------|----------|
| `prompt` | Shows next ready plan, asks to enter |
| `silent` | Automatically enters next ready plan |
| `disabled` | Just reports what's ready |

## Vision Section Design

The vision answers five questions that stay visible throughout work:

1. **End goal** — What does success look like?
2. **Why it matters** — How does this connect to the larger project?
3. **Current state** — Where do things stand right now?
4. **Decisions** — What forks were taken, with rationale?
5. **Exit criteria** — How do we know when we're done?

### Decisions as Fork Records

Record decisions when the work could have gone multiple ways:

```yaml
decisions:
  - date: "2026-02-05"
    fork: "Test region selection approach"
    chose: "Data-driven from cluster markers"
    rationale: >
      Manual selection biases toward expected results.
      Data-driven approach uses marker gene expression
      to find regions with clear cell type boundaries.
```

### Updating Current State

Update `current_state` as work progresses. This is the first thing you see when re-entering a sequence:

```yaml
current_state: >
  Phase 1 complete: all segmentation methods tested on Pat1_P1.
  Phase 2 in progress: ProSeg and Baysor running on full cohort.
  Blocked: Cellpose needs GPU allocation approval.
```

## Plan Lifecycle: Draft → Validate → Sync

Plans follow a pipeline from draft creation through validation to the shared index.

### Overview

```
1. Agent writes plan       2. Pre-exit validates     3. Post-exit breadcrumb    4. Next session sees it
~/.claude/plans/foo.md  →  hook finds via            hook writes to             session_orient.sh iterates
  includes draft_path:     draft_path match           .pending_sync.d/           .pending_sync.d/*
  in frontmatter           validates structure        <area>__<name>             shows ACTION REQUIRED

5. Agent syncs
./plan sync <path>      →  copies to active_plans/ → rebuilds index → reconciles YAML → removes breadcrumb
```

### `draft_path` Field

Every plan should include `draft_path:` in its YAML frontmatter — the absolute path to
the plan file itself. This is how the pre-exit hook finds the correct file to validate
when multiple sessions are active.

```yaml
---
area: analysis
name: performance_audit
draft_path: /Users/you/.claude/plans/random-session-name.md
status: in_progress
# ... other fields
---
```

- **Written by**: the agent (when creating the plan file)
- **Read by**: pre-exit hook (find the right file), post-exit hook (write breadcrumb)
- **Validation**: recognized but not required (fallback: `-mmin -15` recency scan)
- **After sync**: ignored in `active_plans/` copies

Without `draft_path`, the hook falls back to a time-based scan (`-mmin -15`), which
can pick the wrong file in multi-session scenarios.

### Breadcrumbs (`.pending_sync.d/`)

After `ExitPlanMode`, context may reset. The post-exit hook writes a breadcrumb
file to `~/.claude/plans/.pending_sync.d/` named `<area>__<name>`, containing the
draft path. This is:

- **Multi-session safe**: each plan gets its own breadcrumb file
- **Read by** `session_orient.sh` at session start → shows per-plan ACTION REQUIRED
- **Cleaned up by** `./plan sync` → deletes the matching breadcrumb after successful sync

### Hooks

| Hook | Event | Purpose |
|------|-------|---------|
| `pre_exit_plan.sh` | PreToolUse (ExitPlanMode) | Finds plan via `draft_path`, validates frontmatter/sections. Blocks if invalid. |
| `post_exit_plan.sh` | PostToolUse (ExitPlanMode) | Writes breadcrumb to `.pending_sync.d/<area>__<name>`. |
| `session_orient.sh` | SessionStart | Iterates `.pending_sync.d/*`, injects sequence status, warns on hook divergence. |
| `pre_compact_plan.sh` | PreCompact | Injects plan state before context compaction so agent retains orientation. |
| `plan_md_reminder.sh` | PreToolUse (Edit\|Write) | Reminds about plan structure when editing `.md` in plan mode. |

Template scripts for all hooks are in `templates/hooks/`.

### Settings Configuration

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

**Important**: Project-level hooks override global hooks **per event**. If the project
defines `SessionStart`, global `SessionStart` hooks do not fire. The session orient
hook warns when global events are missing from the project config.

## Index Integration

Sequences are included in `.index.json` when the index is rebuilt:

```json
{
  "sequences": {
    "my_pipeline": {
      "status": "active",
      "progress": {"completed": 3, "total": 7},
      "ready_plans": ["analysis/next_step"],
      "current_state": "Phase 2 in progress..."
    }
  }
}
```

## Session Handoff

Session prompts are self-contained documents for handing work to new Claude Code sessions. They sit between the sequence (what needs to happen) and execution (the actual work).

### Architecture

```
Sequence YAML          Session Prompts         Plan Files          Memory
(source of truth)      (handoff layer)         (work tracking)     (cross-session)

xenium_l1_l2.yaml  →   sessions/06_Xenium/  →  active_plans/    →  .claude/memory/
  - vision               L1_loading/             06_Xenium/
  - plan ordering           load.md                 L1_load.md
  - dependencies         L2a_pipeline_setup/
  - state                   proseg.md
```

### Directory Structure

```
sessions/
├── {area}/
│   ├── {working_unit}/
│   │   ├── {session_a}.md
│   │   └── {session_b}.md
│   └── {working_unit}/
│       └── {session_c}.md
└── 00_README.md
```

### Session Prompt Format

Session prompts have YAML frontmatter linking them to the sequence:

```yaml
---
sequence: my_pipeline           # Sequence this session belongs to
plan_id: area/plan_name         # Plan ID in the sequence YAML
mode: plan                      # plan | execute | verify
---
```

Sections: Sequence Context, Objective, What Exists, Reference Implementations, Scope, Constraints, Exit Criteria, Next Plans.

See `graph-memory/templates/session_prompt_template.md` for the full template.

### Session Modes

| Mode | When to use | Agent behavior |
|------|-------------|----------------|
| `plan` | Design sessions, complex implementations | Explore codebase, write plan file, get approval, then work |
| `execute` | Well-specified jobs, script submissions | Start implementing immediately |
| `verify` | Ports with reference implementations | Quick check, then proceed |

### Workflow

1. **Check readiness**: `./plan session list` shows which sessions are ready based on sequence dependencies
2. **Copy prompt**: `cat sessions/area/unit/session.md | pbcopy`
3. **Paste into new Claude Code session**
4. **Agent enters mode** based on frontmatter:
   - `plan` mode: explores, writes plan file in context of vision + sequence, gets approval, works
   - `execute` mode: starts implementing immediately
   - `verify` mode: quick verification pass, then proceeds
5. **On completion**: `./plan seq-complete "summary"` — updates sequence, unblocks downstream

### Readiness

`plan session list` cross-references session frontmatter with sequence state:

```
Session Prompts:

  06_Xenium/
    ● 06_Xenium/L1_loading/load  [execute]
    ○ 06_Xenium/L2a_pipeline_setup/proseg  [verify]
    ○ 06_Xenium/L2b_benchmarking/acgs  [plan]

Legend: ● ready  ○ blocked  ✓ complete  ? no sequence
```

### Design Properties

- **Self-contained**: Each prompt includes everything a new session needs. No assumption of prior conversation history.
- **Layered**: Prompts point to `.knowledge/` docs for method details rather than inlining. Focus on WHAT to do, not HOW things work.
- **Sequence-aware**: Frontmatter links prompt to sequence plan, enabling automated readiness checks.
- **Mode-driven**: The `mode` field controls initial posture, reducing unnecessary planning for straightforward work.

## Commands Reference

| Command | Description |
|---------|-------------|
| `plan sequence list` | All sequences with rollup progress |
| `plan sequence status <n>` | Full status: vision + all plans |
| `plan sequence vision <n>` | Vision section only |
| `plan sequence next <n>` | Next ready plan ID |
| `plan sequence next <n> --json` | Next ready plan as JSON |
| `plan sequence validate <n>` | Check for errors and cycles |
| `plan sequence scope <n>` | Vision with full decision history |
| `plan enter <seq/plan>` | Claim plan, create .md, show context |
| `plan pause "notes"` | Save notes, release ownership |
| `plan seq-complete "summary"` | Complete plan, unblock downstream |
| `plan sync [path]` | Sync draft to active_plans + sequence YAML |
| `plan nudge [path]` | Write sync reminder for next session |
| `plan reconcile [seq]` | Cross-reference plan files vs sequence YAML |
| `plan reconcile --fix <seq>` | Fix status mismatches (YAML authoritative) |
| `plan session list` | List session prompts with readiness |
| `plan session show <path>` | Display a session prompt |

## Tips

1. **Start with vision** — Define `end_goal` and `exit_criteria` before listing plans
2. **Keep plans coarse** — Each plan should represent 1-3 sessions of work
3. **Update `current_state`** — This is the cheapest way to maintain context
4. **Record decisions** — Future sessions will thank you for the rationale
5. **Use `plan sequence vision`** — Quick re-orient before starting work
6. **Validate after edits** — `plan sequence validate` catches cycles and typos
