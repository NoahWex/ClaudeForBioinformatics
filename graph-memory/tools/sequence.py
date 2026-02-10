#!/usr/bin/env python3
"""
Sequence manager - Orchestrate multi-plan workflows with dependency tracking.

Features:
- YAML sequence definitions with dependency graphs
- Vision section for end-goal visibility
- Session ownership tracking for coordination
- Ready-plan detection based on completion status
- Circular dependency validation
"""

import json
import re
import yaml
import os
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
from datetime import date


# Sequences are like playlists but for plans. Skip track? Your project falls apart.


@dataclass
class SequencePlan:
    """A plan within a sequence with coordination state."""
    id: str                             # Full plan ID: area/name
    short_name: str                     # Just the name part
    depends_on: list[str]               # Short names of dependencies
    summary: str = ""                   # One-line description
    status: str = "queued"              # queued | in_progress | complete | blocked
    owner: Optional[str] = None         # Session ID when claimed
    notes: Optional[str] = None         # Latest progress/pause notes
    blocked_reason: Optional[str] = None  # If blocked, why
    completed_summary: Optional[str] = None  # Filled on completion
    claimed_at: Optional[str] = None         # ISO date when claimed

    @classmethod
    def from_dict(cls, data: dict) -> "SequencePlan":
        """Parse from YAML dict entry."""
        plan_id = data.get("id", "")
        # Extract short name from full ID
        short_name = plan_id.split("/")[-1] if "/" in plan_id else plan_id
        return cls(
            id=plan_id,
            short_name=short_name,
            depends_on=data.get("depends_on", []),
            summary=data.get("summary", ""),
            status=data.get("status", "queued"),
            owner=data.get("owner"),
            notes=data.get("notes"),
            blocked_reason=data.get("blocked_reason"),
            completed_summary=data.get("completed_summary"),
            claimed_at=data.get("claimed_at")
        )

    def to_dict(self) -> dict:
        """Convert to dict for YAML serialization."""
        d = {
            "id": self.id,
            "summary": self.summary,
            "status": self.status,
            "depends_on": self.depends_on
        }
        if self.owner:
            d["owner"] = self.owner
        if self.notes:
            d["notes"] = self.notes
        if self.blocked_reason:
            d["blocked_reason"] = self.blocked_reason
        if self.completed_summary:
            d["completed_summary"] = self.completed_summary
        if self.claimed_at:
            d["claimed_at"] = self.claimed_at
        return d


@dataclass
class VisionDecision:
    """A recorded decision at a fork point."""
    date: str
    fork: str                   # What was being decided
    chose: str                  # What was chosen
    rationale: str = ""         # Why (with tradeoffs)

    @classmethod
    def from_dict(cls, data: dict) -> "VisionDecision":
        return cls(
            date=str(data.get("date", "")),
            fork=data.get("fork", data.get("decision", "")),  # Support both formats
            chose=data.get("chose", ""),
            rationale=data.get("rationale", "")
        )

    def to_dict(self) -> dict:
        d = {"date": self.date, "fork": self.fork, "chose": self.chose}
        if self.rationale:
            d["rationale"] = self.rationale
        return d


@dataclass
class Vision:
    """Vision section - the 'why' that stays visible."""
    end_goal: str               # What success looks like
    why_it_matters: str = ""    # Connection to larger project
    current_state: str = ""     # Plain-language where things stand
    decisions: list[VisionDecision] = field(default_factory=list)
    exit_criteria: list[str] = field(default_factory=list)

    @classmethod
    def from_dict(cls, data: dict) -> "Vision":
        decisions = []
        for item in data.get("decisions", []):
            decisions.append(VisionDecision.from_dict(item))

        return cls(
            end_goal=data.get("end_goal", data.get("goal", "")),  # Support both
            why_it_matters=data.get("why_it_matters", ""),
            current_state=data.get("current_state", ""),
            decisions=decisions,
            exit_criteria=data.get("exit_criteria", [])
        )

    def to_dict(self) -> dict:
        return {
            "end_goal": self.end_goal,
            "why_it_matters": self.why_it_matters,
            "current_state": self.current_state,
            "decisions": [d.to_dict() for d in self.decisions],
            "exit_criteria": self.exit_criteria
        }

    def validate(self) -> list[str]:
        """Return validation errors for vision section."""
        errors = []
        if not self.end_goal:
            errors.append("Vision missing required field: end_goal")
        if not self.exit_criteria:
            errors.append("Vision missing required field: exit_criteria")
        return errors

    def format_display(self, width: int = 60) -> str:
        """Format vision for terminal display."""
        lines = []
        sep = "═" * width

        lines.append(sep)
        lines.append("END GOAL:")
        for line in self.end_goal.strip().split('\n'):
            lines.append(f"  {line}")
        lines.append("")

        if self.why_it_matters:
            lines.append("WHY IT MATTERS:")
            for line in self.why_it_matters.strip().split('\n'):
                lines.append(f"  {line}")
            lines.append("")

        if self.current_state:
            lines.append("CURRENT STATE:")
            for line in self.current_state.strip().split('\n'):
                lines.append(f"  {line}")
            lines.append("")

        if self.decisions:
            lines.append("KEY DECISIONS:")
            for decision in self.decisions:
                lines.append(f"  [{decision.date}] {decision.fork}")
                if decision.chose:
                    lines.append(f"    → {decision.chose}")
            lines.append("")

        if self.exit_criteria:
            lines.append("EXIT CRITERIA:")
            for i, criterion in enumerate(self.exit_criteria, 1):
                lines.append(f"  {i}. {criterion}")

        lines.append(sep)
        return '\n'.join(lines)


@dataclass
class Sequence:
    """A sequence of ordered plans with vision and coordination."""
    name: str
    description: str
    status: str                     # active | paused | complete
    auto_progress: str              # prompt | silent | disabled
    vision: Vision                  # The 'why' that stays visible
    plans: list[SequencePlan]
    created: str
    updated: str
    sessions_dir: str = ""              # Relative path to session prompts directory
    path: Optional[Path] = None

    @classmethod
    def from_yaml(cls, content: str, path: Optional[Path] = None) -> "Sequence":
        """Parse sequence from YAML content."""
        data = yaml.safe_load(content) or {}

        # Parse plans list
        plans = []
        for plan_data in data.get("plans", []):
            plans.append(SequencePlan.from_dict(plan_data))

        # Parse vision section
        vision_data = data.get("vision")
        scope_data = data.get("scope")

        if scope_data:
            raise ValueError(
                f"Sequence uses deprecated 'scope:' key. "
                f"Migrate to 'vision:' format. See /sequence-format for schema."
            )

        if vision_data:
            vision = Vision.from_dict(vision_data)
        else:
            vision = Vision(end_goal="", exit_criteria=[])

        # Handle date fields
        created = data.get("created", "")
        updated = data.get("updated", "")
        if hasattr(created, "isoformat"):
            created = created.isoformat()
        if hasattr(updated, "isoformat"):
            updated = updated.isoformat()

        # Derive sessions_dir: explicit field or convention
        seq_name = data.get("name", "")
        sessions_dir = data.get("sessions_dir", "")
        if not sessions_dir and seq_name:
            sessions_dir = f"sessions/sequences/{seq_name}/"

        return cls(
            name=seq_name,
            description=data.get("description", ""),
            status=data.get("status", "active"),
            auto_progress=data.get("auto_progress", "prompt"),
            vision=vision,
            plans=plans,
            created=str(created),
            updated=str(updated),
            sessions_dir=sessions_dir,
            path=path
        )

    def to_yaml(self) -> str:
        """Serialize sequence to YAML."""
        data = {
            "name": self.name,
            "description": self.description,
            "status": self.status,
            "auto_progress": self.auto_progress,
            "sessions_dir": self.sessions_dir,
            "created": self.created,
            "updated": date.today().isoformat(),
            "vision": self.vision.to_dict(),
            "plans": [p.to_dict() for p in self.plans]
        }
        return yaml.dump(data, default_flow_style=False, sort_keys=False, allow_unicode=True)

    def save(self):
        """Save sequence back to its YAML file."""
        if not self.path:
            raise ValueError("Sequence has no path set")
        self.updated = date.today().isoformat()
        self.path.write_text(self.to_yaml())

    def get_plan(self, plan_id: str) -> Optional[SequencePlan]:
        """Get a plan by ID or short name."""
        for plan in self.plans:
            if plan.id == plan_id or plan.short_name == plan_id:
                return plan
        return None

    def claim_plan(self, plan_id: str, session_id: str) -> bool:
        """Claim a plan for a session. Returns True if successful."""
        plan = self.get_plan(plan_id)
        if not plan:
            return False
        if plan.owner and plan.owner != session_id:
            return False  # Already owned by another session
        plan.owner = session_id
        plan.status = "in_progress"
        plan.claimed_at = date.today().isoformat()
        return True

    def release_plan(self, plan_id: str) -> bool:
        """Release ownership of a plan."""
        plan = self.get_plan(plan_id)
        if not plan:
            return False
        plan.owner = None
        plan.claimed_at = None
        # Keep status as in_progress if work was done
        return True

    def complete_plan(self, plan_id: str, summary: str) -> list[str]:
        """
        Mark a plan complete and return newly unblocked plans.

        Returns list of plan IDs that are now ready.
        """
        plan = self.get_plan(plan_id)
        if not plan:
            return []

        plan.status = "complete"
        plan.owner = None
        plan.claimed_at = None
        plan.completed_summary = summary

        # Find plans that might be unblocked
        short_to_full = {p.short_name: p.id for p in self.plans}
        unblocked = []

        for p in self.plans:
            if p.status not in ("queued", "blocked"):
                continue
            # Check if all deps are now complete
            all_complete = all(
                self.get_plan(short_to_full.get(dep, "")).status == "complete"
                for dep in p.depends_on
                if self.get_plan(short_to_full.get(dep, ""))
            )
            if all_complete and p.depends_on:  # Had deps, now all complete
                if p.status == "blocked":
                    p.status = "queued"
                    p.blocked_reason = None
                unblocked.append(p.id)

        return unblocked

    def update_notes(self, plan_id: str, notes: str):
        """Update notes for a plan."""
        plan = self.get_plan(plan_id)
        if plan:
            plan.notes = notes

    def update_current_state(self, state: str):
        """Update the vision's current state."""
        self.vision.current_state = state

    def compute_downstream(self, plan_id: str) -> list[SequencePlan]:
        """Plans whose depends_on includes this plan's short_name."""
        plan = self.get_plan(plan_id)
        if not plan:
            return []
        return [p for p in self.plans if plan.short_name in p.depends_on]

    def get_ready_plans(self, plan_statuses: Optional[dict[str, str]] = None) -> list[str]:
        """
        Get plans whose dependencies are all complete.

        Uses internal plan statuses if plan_statuses not provided.
        """
        ready = []
        short_to_full = {p.short_name: p.id for p in self.plans}

        for plan in self.plans:
            # Get status from sequence's own tracking or external
            if plan_statuses:
                plan_status = plan_statuses.get(plan.id, plan.status)
            else:
                plan_status = plan.status

            # Skip if already complete or in progress
            if plan_status in ("complete", "in_progress"):
                continue

            # Check all dependencies
            deps_complete = True
            for dep in plan.depends_on:
                dep_full_id = short_to_full.get(dep)
                if dep_full_id:
                    dep_plan = self.get_plan(dep_full_id)
                    if dep_plan:
                        dep_status = plan_statuses.get(dep_full_id, dep_plan.status) if plan_statuses else dep_plan.status
                        if dep_status != "complete":
                            deps_complete = False
                            break

            if deps_complete:
                ready.append(plan.id)

        return ready

    def compute_progress(self, plan_statuses: Optional[dict[str, str]] = None) -> tuple[int, int]:
        """Compute (completed, total) counts."""
        completed = 0
        for p in self.plans:
            status = plan_statuses.get(p.id, p.status) if plan_statuses else p.status
            if status == "complete":
                completed += 1
        return completed, len(self.plans)

    def validate(self, available_plans: Optional[set[str]] = None) -> list[str]:
        """Validate sequence for errors."""
        errors = []

        if not self.name:
            errors.append("Missing required field: name")

        # Validate vision
        errors.extend(self.vision.validate())

        # Check plans exist (if available_plans provided)
        if available_plans:
            for plan in self.plans:
                if plan.id not in available_plans:
                    errors.append(f"Plan not found: {plan.id}")

        # Check for circular dependencies
        short_to_deps = {p.short_name: p.depends_on for p in self.plans}
        valid_names = set(short_to_deps.keys())

        for plan in self.plans:
            for dep in plan.depends_on:
                if dep not in valid_names:
                    errors.append(f"Unknown dependency '{dep}' in plan {plan.id}")

        # Detect cycles
        visited = set()
        rec_stack = set()

        def has_cycle(node: str, path: list[str]) -> Optional[list[str]]:
            visited.add(node)
            rec_stack.add(node)
            path = path + [node]

            for dep in short_to_deps.get(node, []):
                if dep not in valid_names:
                    continue
                if dep not in visited:
                    cycle_path = has_cycle(dep, path)
                    if cycle_path:
                        return cycle_path
                elif dep in rec_stack:
                    cycle_start = path.index(dep) if dep in path else 0
                    return path[cycle_start:] + [dep]

            rec_stack.remove(node)
            return None

        for plan in self.plans:
            if plan.short_name not in visited:
                cycle = has_cycle(plan.short_name, [])
                if cycle:
                    errors.append(f"Circular dependency: {' -> '.join(cycle)}")
                    break

        if self.auto_progress not in ("prompt", "silent", "disabled"):
            errors.append(f"Invalid auto_progress value: {self.auto_progress}")

        if self.status not in ("active", "paused", "complete"):
            errors.append(f"Invalid status value: {self.status}")

        return errors

    def format_status(self, width: int = 60) -> str:
        """Format full status display for terminal."""
        lines = []
        sep = "═" * width
        thin_sep = "─" * width

        lines.append(sep)
        lines.append(f"SEQUENCE: {self.name}")
        lines.append(sep)
        lines.append("")

        # Vision
        lines.append("END GOAL:")
        for line in self.vision.end_goal.strip().split('\n'):
            lines.append(f"  {line.strip()}")
        lines.append("")

        if self.vision.why_it_matters:
            lines.append("WHY IT MATTERS:")
            for line in self.vision.why_it_matters.strip().split('\n'):
                lines.append(f"  {line.strip()}")
            lines.append("")

        if self.vision.current_state:
            lines.append("CURRENT STATE:")
            for line in self.vision.current_state.strip().split('\n'):
                lines.append(f"  {line.strip()}")
            lines.append("")

        # Plans
        lines.append(thin_sep)
        lines.append(f"{'PLANS':<40} {'STATUS':>18}")
        lines.append(thin_sep)

        status_icons = {
            "complete": "✓",
            "in_progress": "●",
            "queued": "○",
            "blocked": "⊘"
        }

        for plan in self.plans:
            icon = status_icons.get(plan.status, "?")
            status_str = plan.status
            if plan.owner:
                status_str += " (owned)"

            lines.append(f"{icon} {plan.short_name:<38} {status_str:>18}")
            if plan.summary:
                lines.append(f"  └─ {plan.summary[:55]}")
            if plan.notes:
                lines.append(f"  └─ Notes: {plan.notes[:50]}")
            if plan.depends_on:
                deps_str = ", ".join(plan.depends_on)
                lines.append(f"  └─ Depends on: {deps_str}")
            if plan.claimed_at and plan.status == "in_progress":
                try:
                    claimed_date = date.fromisoformat(plan.claimed_at)
                    days = (date.today() - claimed_date).days
                    if days > 1:
                        lines.append(f"  └─ ⚠ STALE CLAIM: claimed {days} days ago ({plan.claimed_at})")
                except (ValueError, TypeError):
                    pass

        lines.append(thin_sep)

        # Ready plans
        ready = self.get_ready_plans()
        if ready:
            lines.append(f"READY TO WORK: {', '.join(ready)}")
        else:
            completed, total = self.compute_progress()
            if completed == total:
                lines.append("SEQUENCE COMPLETE")
            else:
                lines.append("NO PLANS READY (check blockers)")

        lines.append(thin_sep)

        return '\n'.join(lines)


class SequenceManager:
    """Manages all sequences in a project."""

    def __init__(self, plans_dir: Path):
        self.plans_dir = plans_dir
        self.sequences_dir = plans_dir / "sequences"
        self.active_dir = plans_dir / "active_plans"
        self.sequences: dict[str, Sequence] = {}

    def load_sequence(self, name: str) -> Optional[Sequence]:
        """Load a single sequence by name."""
        path = self.sequences_dir / f"{name}.yaml"
        if not path.exists():
            return None

        content = path.read_text()
        seq = Sequence.from_yaml(content, path)
        self.sequences[seq.name] = seq
        return seq

    def load_all(self) -> dict[str, Sequence]:
        """Load all sequences from the sequences directory."""
        self.sequences = {}

        if not self.sequences_dir.exists():
            return self.sequences

        for yaml_file in self.sequences_dir.glob("*.yaml"):
            if yaml_file.name.startswith("._"):
                continue
            try:
                content = yaml_file.read_text()
                seq = Sequence.from_yaml(content, yaml_file)
                if seq.name:
                    self.sequences[seq.name] = seq
            except Exception as e:
                print(f"Warning: Failed to parse {yaml_file}: {e}")

        return self.sequences

    def get_sequence(self, name: str) -> Optional[Sequence]:
        """Get a sequence by name, loading if necessary."""
        if name in self.sequences:
            return self.sequences[name]
        return self.load_sequence(name)

    def find_sequences_containing(self, plan_id: str) -> list[str]:
        """Find all sequences that contain a given plan."""
        self.load_all()
        containing = []

        for seq_name, seq in self.sequences.items():
            for plan in seq.plans:
                if plan.id == plan_id:
                    containing.append(seq_name)
                    break

        return containing

    def get_setting(self, seq_name: str, setting: str) -> Optional[str]:
        """Get a setting value from a sequence."""
        seq = self.get_sequence(seq_name)
        if not seq:
            return None
        return getattr(seq, setting, None)

    def create_plan_file(self, seq_name: str, plan_id: str) -> Optional[Path]:
        """Create a plan .md file from template if it doesn't exist."""
        seq = self.get_sequence(seq_name)
        if not seq:
            return None

        plan = seq.get_plan(plan_id)
        if not plan:
            return None

        # Determine file path
        area, name = plan.id.split("/", 1) if "/" in plan.id else ("", plan.id)
        plan_dir = self.active_dir / area
        plan_path = plan_dir / f"{name}.md"

        if plan_path.exists():
            return plan_path

        # Create from template
        plan_dir.mkdir(parents=True, exist_ok=True)

        today = date.today().isoformat()
        vision_line = seq.vision.end_goal.strip().split('\n')[0][:60]

        content = f'''---
area: {area}
name: {name}
version: 1
updated: {today}
supersedes: null
sequence: {seq.name}
status: in_progress
blocked_by: []
blocks: []
jobs: []
---

# {plan.summary or name}

Part of sequence: {seq.name}
Vision: {vision_line}

## Inputs

- (to be filled during work)

## Outputs

- (to be filled during work)

## Config

- (to be filled during work)

## Scripts

- (to be filled during work)

## Tasks

- [ ] (to be filled during work)

## Current State

{plan.notes or "Starting work."}

## Discussion

### {today}: Initial entry

Work beginning on this plan.

## Next Actions

- Review inputs and define approach

## Completion

When tasks complete, run:
```bash
./plan complete {plan.id}
```
'''
        plan_path.write_text(content)
        return plan_path

    def generate_session_id(self) -> str:
        """Generate a unique session ID."""
        import uuid
        return f"session_{uuid.uuid4().hex[:8]}"

    def _scan_existing_sessions(self, sessions_dir: Path) -> dict[str, Path]:
        """Scan sessions/ recursively, return {plan_id: path} mapping.

        Matches by plan_id frontmatter field, handling reorganized files.
        """
        mapping: dict[str, Path] = {}
        if not sessions_dir.exists():
            return mapping

        for md_file in sessions_dir.rglob("*.md"):
            if md_file.name.startswith("._") or md_file.name.startswith("00_") or md_file.name.upper().startswith("README"):
                continue
            fm = parse_session_frontmatter(md_file)
            plan_id = fm.get("plan_id", "")
            if plan_id:
                mapping[plan_id] = md_file

        return mapping

    def _render_plan_stub(self, seq: Sequence, plan: SequencePlan) -> str:
        """Render a plan .md with standard frontmatter + empty sections."""
        area, name = plan.id.split("/", 1) if "/" in plan.id else ("", plan.id)
        today = date.today().isoformat()

        # Compute blocked_by from depends_on (full IDs)
        short_to_full = {p.short_name: p.id for p in seq.plans}
        blocked_by = [short_to_full[d] for d in plan.depends_on if d in short_to_full]

        # Compute blocks (downstream plans)
        downstream = seq.compute_downstream(plan.id)
        blocks = [p.id for p in downstream]

        blocked_by_yaml = json.dumps(blocked_by)
        blocks_yaml = json.dumps(blocks)

        return f'''---
area: {area}
name: {name}
version: 1
updated: {today}
supersedes: null
sequence: {seq.name}
status: {plan.status}
blocked_by: {blocked_by_yaml}
blocks: {blocks_yaml}
jobs: []
---

# {plan.summary or name}

Part of sequence: {seq.name}

## Inputs

- (to be filled during planning)

## Outputs

- (to be filled during planning)

## Config

- (to be filled during planning)

## Scripts

- (to be filled during planning)

## Tasks

- [ ] (to be filled during planning)

## Current State

Stub created by `plan session generate`. Awaiting session work.

## Discussion

### {today}: Plan created

Auto-generated from sequence `{seq.name}`.

## Next Actions

1. Review inputs and define approach

## Completion

When all tasks are complete, run:

```bash
./plan complete {plan.id}
```
'''

    def _render_session_skeleton(self, seq: Sequence, plan: SequencePlan) -> str:
        """Render a session .md with frontmatter + TODO markers."""
        area, name = plan.id.split("/", 1) if "/" in plan.id else ("", plan.id)

        # Compute downstream plans
        downstream = seq.compute_downstream(plan.id)
        next_plans_lines = ""
        if downstream:
            for dp in downstream:
                next_plans_lines += f"- `{dp.id}` — {dp.summary}\n"
        else:
            next_plans_lines = "- (none — this is a terminal plan)\n"

        return f'''---
sequence: {seq.name}
plan_id: {plan.id}
mode: plan
---

# Session: {name} — {plan.summary or name}

## Sequence Context

You are working on plan `{plan.id}` in the `{seq.name}` sequence.
Sequence YAML: `project/.plans/sequences/{seq.name}.yaml`

## Workflow

### Planning
1. Run `/plan-format` to see the plan file template and required fields
2. Enter plan mode (`EnterPlanMode`)
3. Include `plan_target: {plan.id}` in your plan's frontmatter — the pre-exit hook
   will auto-inject the full frontmatter from the existing stub
4. Write the required sections (Inputs, Outputs, Config, Scripts, Tasks, etc.)
5. Exit plan mode — the pre-exit hook validates your plan format

### After Plan Approval (FIRST action before any implementation)
6. Run `./plan sync <path>` where `<path>` is your plan file
   (e.g. `./plan sync ~/.claude/plans/some-session-name.md`)
   This publishes to `active_plans/` and updates the sequence YAML.
   **Do not start implementation until sync completes.**

## Objective

{plan.summary or "(TODO: describe what this session achieves)"}

## What Exists

- (TODO: list existing files, scripts, configs relevant to this plan)

## Scope

1. (TODO: concrete numbered steps)

## Constraints

- Git writes must use `hpc git` (CRSP mount)
- (TODO: add container, environment, other constraints)

## Exit Criteria

- (TODO: verifiable outcomes)
- Update sequence plan status to complete

## Next Plans (unblocked by completion)

{next_plans_lines}
## Completion

When all work is done:
```bash
./plan seq-complete "summary of what was done"
```
This marks the plan complete in the sequence YAML and unblocks downstream plans.
'''

    def generate_session_skeletons(
        self,
        seq_name: str,
        dry_run: bool = False
    ) -> dict:
        """Generate plan stubs + session prompts for a sequence.

        Additive-only: never overwrites existing files.

        Returns dict with:
            generated: [(plan_id, stub_path, session_path), ...]
            skipped: [(plan_id, reason), ...]
            errors: [(plan_id, error), ...]
        """
        result = {"generated": [], "skipped": [], "errors": []}

        seq = self.get_sequence(seq_name)
        if not seq:
            result["errors"].append(("", f"Sequence not found: {seq_name}"))
            return result

        sessions_dir = self.plans_dir / seq.sessions_dir
        existing_sessions = self._scan_existing_sessions(sessions_dir)

        for plan in seq.plans:
            area, name = plan.id.split("/", 1) if "/" in plan.id else ("", plan.id)

            # Check if plan stub exists
            stub_path = self.active_dir / area / f"{name}.md"
            stub_exists = stub_path.exists()

            # Check if session prompt exists (by plan_id in frontmatter)
            session_exists_path = existing_sessions.get(plan.id)

            if stub_exists and session_exists_path:
                rel_session = session_exists_path.relative_to(sessions_dir) if session_exists_path.is_relative_to(sessions_dir) else session_exists_path
                result["skipped"].append((plan.id, f"stub + session exist (session at {rel_session})"))
                continue

            if dry_run:
                parts = []
                if not stub_exists:
                    parts.append(f"active_plans/{area}/{name}.md")
                if not session_exists_path:
                    parts.append(f"{seq.sessions_dir}{name}.md")
                result["generated"].append((plan.id, parts[0] if parts else "", parts[-1] if len(parts) > 1 else ""))
                continue

            try:
                # Generate stub if missing
                actual_stub_path = None
                if not stub_exists:
                    stub_content = self._render_plan_stub(seq, plan)
                    stub_path.parent.mkdir(parents=True, exist_ok=True)
                    stub_path.write_text(stub_content)
                    actual_stub_path = stub_path

                # Generate session if missing
                # Sessions go under sequences/{seq_name}/ (flat within sequence)
                actual_session_path = None
                if not session_exists_path:
                    session_content = self._render_session_skeleton(seq, plan)
                    session_path = sessions_dir / f"{name}.md"
                    session_path.parent.mkdir(parents=True, exist_ok=True)
                    session_path.write_text(session_content)
                    actual_session_path = session_path

                result["generated"].append((
                    plan.id,
                    str(actual_stub_path) if actual_stub_path else "(exists)",
                    str(actual_session_path) if actual_session_path else "(exists)"
                ))
            except Exception as e:
                result["errors"].append((plan.id, str(e)))

        return result

    @staticmethod
    def update_plan_frontmatter(plan_path: Path, updates: dict[str, str]) -> bool:
        """Update scalar fields in a plan .md file's YAML frontmatter.

        Uses regex replacement to preserve exact formatting. Only handles
        simple scalar values (status, updated, sequence, etc.).

        Args:
            plan_path: Path to the plan .md file
            updates: Dict of field_name → new_value to replace

        Returns:
            True if any field was updated, False otherwise.
        """
        try:
            text = plan_path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            return False

        if not text.startswith("---"):
            return False

        end = text.find("---", 3)
        if end == -1:
            return False

        header = text[3:end]
        changed = False

        for field_name, new_value in updates.items():
            # Match: field_name: <anything> (but not multiline values)
            pattern = re.compile(
                rf'^({re.escape(field_name)}:\s*)(.*)$',
                re.MULTILINE
            )
            match = pattern.search(header)
            if match:
                old_line = match.group(0)
                new_line = f"{match.group(1)}{new_value}"
                if old_line != new_line:
                    header = header[:match.start()] + new_line + header[match.end():]
                    changed = True

        if changed:
            plan_path.write_text(text[:3] + header + text[end:], encoding="utf-8")

        return changed

    def sync_plan_to_yaml(self, plan_id: str, plan_status: str,
                          plan_blocks: list[str] = None,
                          plan_blocked_by: list[str] = None) -> dict:
        """Sync plan file state to sequence YAML using ruamel for round-trip.

        Auto-syncs: status
        Auto-adds: new depends_on entries from plan blocks
        Flags: dependency removals, mismatches

        Returns dict with: synced (list), conflicts (list), errors (list)
        """
        result = {"synced": [], "conflicts": [], "errors": [], "warnings": []}

        # Find containing sequences
        self.load_all()
        containing = self.find_sequences_containing(plan_id)
        if not containing:
            result["warnings"].append(f"Plan {plan_id} not found in any sequence (standalone plan)")
            return result

        try:
            from ruamel.yaml import YAML
            ryaml = YAML()
            ryaml.preserve_quotes = True
        except ImportError:
            result["errors"].append("ruamel.yaml not available; cannot round-trip YAML")
            return result

        for seq_name in containing:
            seq = self.sequences[seq_name]
            if not seq.path:
                continue

            # Round-trip load
            with open(seq.path) as f:
                data = ryaml.load(f)

            plans_list = data.get("plans", [])
            plan_entry = None
            for entry in plans_list:
                if entry.get("id") == plan_id:
                    plan_entry = entry
                    break

            if not plan_entry:
                continue

            changed = False

            # --- Status sync (plan is authoritative) ---
            yaml_status = plan_entry.get("status", "queued")
            if plan_status and plan_status != yaml_status:
                result["synced"].append(
                    f"[{seq_name}] {plan_id}: status {yaml_status} → {plan_status}")
                plan_entry["status"] = plan_status
                changed = True

            # --- Dependency sync (bidirectional conflict detection) ---
            short_name = plan_id.split("/")[-1] if "/" in plan_id else plan_id
            short_to_full = {p.short_name: p.id for p in seq.plans}
            full_to_short = {p.id: p.short_name for p in seq.plans}

            # Check blocks: plan says it blocks X → X should depend_on this plan
            if plan_blocks is not None:
                for blocked_id in plan_blocks:
                    blocked_short = full_to_short.get(blocked_id, blocked_id.split("/")[-1])
                    # Find the downstream plan in YAML
                    for other in plans_list:
                        if other.get("id") == blocked_id:
                            yaml_deps = other.get("depends_on", [])
                            if short_name not in yaml_deps:
                                result["synced"].append(
                                    f"[{seq_name}] {blocked_id}: added depends_on '{short_name}' (from {plan_id} blocks)")
                                yaml_deps.append(short_name)
                                other["depends_on"] = yaml_deps
                                changed = True
                            break

            # Check blocked_by: plan says it's blocked by X → this should depend_on X
            if plan_blocked_by is not None:
                yaml_deps = plan_entry.get("depends_on", [])
                for blocker_id in plan_blocked_by:
                    blocker_short = full_to_short.get(blocker_id, blocker_id.split("/")[-1])
                    if blocker_short not in yaml_deps:
                        result["synced"].append(
                            f"[{seq_name}] {plan_id}: added depends_on '{blocker_short}' (from plan blocked_by)")
                        yaml_deps.append(blocker_short)
                        plan_entry["depends_on"] = yaml_deps
                        changed = True

            # Detect removals (YAML has dep that plan doesn't claim)
            if plan_blocked_by is not None:
                plan_blocker_shorts = set()
                for bid in plan_blocked_by:
                    plan_blocker_shorts.add(full_to_short.get(bid, bid.split("/")[-1]))
                yaml_deps = set(plan_entry.get("depends_on", []))
                removed = yaml_deps - plan_blocker_shorts
                for dep in removed:
                    result["conflicts"].append(
                        f"[{seq_name}] {plan_id}: YAML depends_on '{dep}' not in plan blocked_by (removal?)")

            if plan_blocks is not None:
                plan_block_shorts = set()
                for bid in plan_blocks:
                    plan_block_shorts.add(full_to_short.get(bid, bid.split("/")[-1]))
                # Check downstream: YAML plans that depend on this but aren't in plan blocks
                for other in plans_list:
                    other_deps = other.get("depends_on", [])
                    if short_name in other_deps:
                        other_short = other.get("id", "").split("/")[-1]
                        other_full = other.get("id", "")
                        if other_full not in plan_blocks and other_full != plan_id:
                            result["conflicts"].append(
                                f"[{seq_name}] {other_full}: depends_on '{short_name}' but {plan_id} doesn't list it in blocks (stale?)")

            if changed:
                with open(seq.path, 'w') as f:
                    ryaml.dump(data, f)

        return result

    def detect_all_conflicts(self) -> list[str]:
        """Compare all plan files against sequence YAMLs. Return conflict descriptions."""
        conflicts = []
        self.load_all()

        for seq_name, seq in self.sequences.items():
            for plan in seq.plans:
                area, name = plan.id.split("/", 1) if "/" in plan.id else ("", plan.id)
                plan_path = self.active_dir / area / f"{name}.md"
                if not plan_path.exists():
                    continue

                # Read plan frontmatter
                try:
                    text = plan_path.read_text(encoding="utf-8")
                except (OSError, UnicodeDecodeError):
                    continue

                if not text.startswith("---"):
                    conflicts.append(f"[{seq_name}] {plan.id}: plan file has no frontmatter")
                    continue

                end = text.find("---", 3)
                if end == -1:
                    continue

                try:
                    fm = yaml.safe_load(text[3:end]) or {}
                except Exception:
                    continue

                plan_status = fm.get("status", "")
                plan_blocks = fm.get("blocks", [])
                plan_blocked_by = fm.get("blocked_by", [])

                # Normalize to lists of strings
                if isinstance(plan_blocks, str):
                    plan_blocks = [plan_blocks]
                if isinstance(plan_blocked_by, str):
                    plan_blocked_by = [plan_blocked_by]

                short_to_full = {p.short_name: p.id for p in seq.plans}
                full_to_short = {p.id: p.short_name for p in seq.plans}
                short_name = plan.id.split("/")[-1]

                # Status mismatch
                if plan_status and plan_status != plan.status:
                    conflicts.append(
                        f"[{seq_name}] {plan.id}: status mismatch — plan={plan_status}, yaml={plan.status}")

                # blocked_by vs depends_on
                plan_dep_shorts = set()
                for bid in plan_blocked_by:
                    plan_dep_shorts.add(full_to_short.get(bid, bid.split("/")[-1]))
                yaml_deps = set(plan.depends_on)

                for dep in plan_dep_shorts - yaml_deps:
                    conflicts.append(
                        f"[{seq_name}] {plan.id}: plan blocked_by '{dep}' not in YAML depends_on")
                for dep in yaml_deps - plan_dep_shorts:
                    conflicts.append(
                        f"[{seq_name}] {plan.id}: YAML depends_on '{dep}' not in plan blocked_by")

                # blocks vs downstream depends_on
                plan_block_ids = set(plan_blocks)
                yaml_downstream = set()
                for p in seq.plans:
                    if short_name in p.depends_on:
                        yaml_downstream.add(p.id)

                for bid in plan_block_ids - yaml_downstream:
                    conflicts.append(
                        f"[{seq_name}] {plan.id}: plan blocks '{bid}' but that plan doesn't depend_on '{short_name}' in YAML")
                for did in yaml_downstream - plan_block_ids:
                    conflicts.append(
                        f"[{seq_name}] {plan.id}: YAML has '{did}' depending on '{short_name}' but plan doesn't list it in blocks")

        return conflicts

    def build_index_section(self, plan_statuses: dict[str, str]) -> dict:
        """Build the sequences section for .index.json."""
        self.load_all()
        index = {}

        for name, seq in self.sequences.items():
            completed, total = seq.compute_progress(plan_statuses)
            ready = seq.get_ready_plans(plan_statuses)

            index[name] = {
                "status": seq.status,
                "auto_progress": seq.auto_progress,
                "progress": {"completed": completed, "total": total},
                "ready_plans": ready,
                "plans": [p.id for p in seq.plans],
                "current_state": seq.vision.current_state[:100] if seq.vision.current_state else ""
            }

        return index


def parse_session_frontmatter(path: Path) -> dict:
    """Parse YAML frontmatter from a session prompt file."""
    try:
        text = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return {}

    if not text.startswith("---"):
        return {}

    end = text.find("---", 3)
    if end == -1:
        return {}

    try:
        fm = yaml.safe_load(text[3:end]) or {}
        return fm
    except Exception:
        return {}


# CLI interface for direct invocation
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: sequence.py <command> <args>")
        print("Commands:")
        print("  find-containing <plan_id>  - Find sequences containing a plan")
        print("  get-setting <seq> <key>    - Get a setting value")
        print("  status <seq>               - Show sequence status")
        print("  vision <seq>               - Show vision only")
        sys.exit(1)

    command = sys.argv[1]

    # Find plans_dir
    current = Path(__file__).parent.parent
    if current.name != ".plans":
        for parent in [current] + list(current.parents):
            if (parent / "project" / ".plans").exists():
                current = parent / "project" / ".plans"
                break
            if (parent / ".plans").exists():
                current = parent / ".plans"
                break

    manager = SequenceManager(current)

    if command == "find-containing":
        if len(sys.argv) < 3:
            print("Usage: sequence.py find-containing <plan_id>")
            sys.exit(1)
        plan_id = sys.argv[2]
        containing = manager.find_sequences_containing(plan_id)
        for seq_name in containing:
            print(seq_name)

    elif command == "get-setting":
        if len(sys.argv) < 4:
            print("Usage: sequence.py get-setting <seq_name> <setting>")
            sys.exit(1)
        seq_name = sys.argv[2]
        setting = sys.argv[3]
        value = manager.get_setting(seq_name, setting)
        if value:
            print(value)

    elif command == "status":
        if len(sys.argv) < 3:
            print("Usage: sequence.py status <seq_name>")
            sys.exit(1)
        seq_name = sys.argv[2]
        seq = manager.get_sequence(seq_name)
        if seq:
            print(seq.format_status())
        else:
            print(f"Sequence not found: {seq_name}")
            sys.exit(1)

    elif command == "vision":
        if len(sys.argv) < 3:
            print("Usage: sequence.py vision <seq_name>")
            sys.exit(1)
        seq_name = sys.argv[2]
        seq = manager.get_sequence(seq_name)
        if seq:
            print(seq.vision.format_display())
        else:
            print(f"Sequence not found: {seq_name}")
            sys.exit(1)

    else:
        print(f"Unknown command: {command}")
        sys.exit(1)
