#!/usr/bin/env python3
"""
Plan CLI - Query interface for the planning system.

Commands:
    plan query ready          - Plans ready to work on
    plan query blocked        - Plans that are blocked
    plan query context <id>   - Cluster context for a plan
    plan query downstream <id> - Plans depending on this one
    plan query touches <path> - Plans referencing a file
    plan query anchors        - MNN anchors (cross-modal links)
    plan query communities    - Show all clusters
    plan query history <id>   - Succession history
    plan query tasks          - All pending tasks across plans
    plan query tasks --all    - Include completed tasks
    plan graph                - Visualize plan graph (default: table)
    plan graph --ascii        - ASCII DAG view
    plan graph --mermaid      - Mermaid markdown syntax
    plan graph --dot          - Graphviz DOT syntax
    plan validate <file>      - Validate a plan file
    plan index rebuild        - Rebuild the index
    plan knowledge list       - List all knowledge entries by type
    plan knowledge search <q> - Search entries by tag/content
    plan knowledge add <path> - Index a new entry
    plan knowledge show <id>  - Display entry details
    plan complete <id>        - Run completion debrief and archive
    plan sync <path>          - Sync draft plan to active_plans + sequence YAML
    plan reconcile [seq]      - Cross-reference plan files against sequence YAML
"""

import sys
import json
import argparse
from pathlib import Path
from typing import Optional

# Add tools directory to path for imports
TOOLS_DIR = Path(__file__).parent
sys.path.insert(0, str(TOOLS_DIR))

import yaml
from parse import parse_plan_file, parse_plan
from graph import PlanGraph
from visualize import visualize
from knowledge import KnowledgeIndex, get_knowledge_dir
from sequence import SequenceManager, parse_session_frontmatter


def _is_interactive() -> bool:
    """Check if stdin is a TTY (interactive terminal)."""
    import os
    return os.isatty(0)


def _prompt(message: str, default: str = "") -> str:
    """Prompt user, falling back to default if no stdin available."""
    try:
        return input(message).strip()
    except EOFError:
        return default


# Yet another CLI wrapper. At least this one queries something useful.


def get_plans_dir() -> Path:
    """Find the .plans directory."""
    # Start from tools dir and go up
    current = TOOLS_DIR.parent
    if current.name == ".plans":
        return current

    # Try to find it in project root
    for parent in [current] + list(current.parents):
        plans_dir = parent / "project" / ".plans"
        if plans_dir.exists():
            return plans_dir
        plans_dir = parent / ".plans"
        if plans_dir.exists():
            return plans_dir

    raise FileNotFoundError("Could not find .plans directory")


def cmd_query_ready(graph: PlanGraph) -> None:
    """Show plans ready to work on."""
    ready = graph.get_ready_plans()

    if not ready:
        print("No plans ready to work on.")
        return

    print("Plans ready to work on:")
    for plan_id in ready:
        plan = graph.plans[plan_id]
        completed, total = plan.task_progress()
        print(f"  {plan_id} ({plan.frontmatter.status}) - {completed}/{total} tasks")


def cmd_query_blocked(graph: PlanGraph) -> None:
    """Show blocked plans and their blockers."""
    blocked = graph.get_blocked_plans()

    if not blocked:
        print("No plans are blocked.")
        return

    print("Blocked plans:")
    for plan_id, blockers in blocked:
        print(f"  {plan_id}")
        print(f"    blocked by: {', '.join(blockers)}")


def cmd_query_context(graph: PlanGraph, plan_id: str) -> None:
    """Show cluster context for a plan."""
    # Find the plan's cluster
    cluster = None
    for c_name, c in graph.clusters.items():
        if plan_id in c.members:
            cluster = c
            break

    if not cluster:
        print(f"Plan {plan_id} not found in any cluster")
        return

    print(f"Cluster: {cluster.name}")
    print(f"Members: {', '.join(cluster.members)}")
    print()

    print("Pooled Inputs:")
    for inp in cluster.pooled_inputs[:10]:
        print(f"  - {inp}")
    if len(cluster.pooled_inputs) > 10:
        print(f"  ... and {len(cluster.pooled_inputs) - 10} more")
    print()

    print("Pooled Outputs:")
    for out in cluster.pooled_outputs[:10]:
        print(f"  - {out}")
    if len(cluster.pooled_outputs) > 10:
        print(f"  ... and {len(cluster.pooled_outputs) - 10} more")
    print()

    print("Pooled Tasks:")
    for task in cluster.pooled_tasks[:10]:
        print(f"  {task}")
    if len(cluster.pooled_tasks) > 10:
        print(f"  ... and {len(cluster.pooled_tasks) - 10} more")
    print()

    print("Recent Discussion:")
    for entry in cluster.pooled_discussion[:5]:
        print(f"  [{entry['date']}] {entry['plan']}: {entry['title']}")


def cmd_query_downstream(graph: PlanGraph, plan_id: str) -> None:
    """Show plans that depend on this one."""
    downstream = graph.get_downstream(plan_id)

    if not downstream:
        print(f"No plans depend on {plan_id}")
        return

    print(f"Plans depending on {plan_id}:")
    for dep_id in downstream:
        if dep_id in graph.plans:
            plan = graph.plans[dep_id]
            print(f"  {dep_id} ({plan.frontmatter.status})")
        else:
            print(f"  {dep_id} (not found)")


def cmd_query_touches(graph: PlanGraph, file_pattern: str) -> None:
    """Show plans that reference a file."""
    matching = graph.get_plans_touching_file(file_pattern)

    if not matching:
        print(f"No plans reference files matching '{file_pattern}'")
        return

    print(f"Plans referencing '{file_pattern}':")
    for plan_id in matching:
        print(f"  {plan_id}")


def cmd_query_anchors(graph: PlanGraph) -> None:
    """Show MNN anchors (plans connected via 2+ edge types)."""
    anchors = graph.get_mnn_anchors()

    if not anchors:
        print("No MNN anchors found")
        return

    print("MNN Anchors (cross-modal connections):")
    for p1, p2, types in anchors:
        print(f"  {p1} <-> {p2}")
        print(f"    edge types: {', '.join(types)}")


def cmd_query_communities(graph: PlanGraph) -> None:
    """Show all clusters."""
    if not graph.clusters:
        print("No clusters defined")
        return

    print("Plan clusters:")
    for name, cluster in graph.clusters.items():
        print(f"\n{name}:")
        for member in cluster.members:
            if member in graph.plans:
                plan = graph.plans[member]
                completed, total = plan.task_progress()
                print(f"  - {member} ({plan.frontmatter.status}, {completed}/{total} tasks)")
            else:
                print(f"  - {member} (not found)")


def cmd_query_tasks(graph: PlanGraph, show_complete: bool = False) -> None:
    """Show all tasks across active plans."""
    all_tasks = []

    for plan_id, plan in graph.plans.items():
        tasks = plan.get_all_tasks()
        for task in tasks:
            if show_complete or not task.complete:
                all_tasks.append((plan_id, task))

    if not all_tasks:
        if show_complete:
            print("No tasks found in any active plan.")
        else:
            print("No pending tasks. All done!")
        return

    # Group by plan
    current_plan = None
    pending_count = 0
    complete_count = 0

    for plan_id, task in sorted(all_tasks, key=lambda x: (x[0], x[1].complete)):
        if plan_id != current_plan:
            if current_plan is not None:
                print()
            print(f"{plan_id}:")
            current_plan = plan_id

        marker = "[x]" if task.complete else "[ ]"
        print(f"  {marker} {task.text}")

        if task.complete:
            complete_count += 1
        else:
            pending_count += 1

    print()
    print(f"Total: {pending_count} pending, {complete_count} complete")


def cmd_query_history(graph: PlanGraph, plan_id: str) -> None:
    """Show succession history for a plan."""
    history = graph.get_succession_history(plan_id)

    print(f"Succession history for {plan_id}:")
    for i, version in enumerate(history):
        prefix = "  " if i < len(history) - 1 else "->"

        # Get version info
        plan = graph.plans.get(version) or graph.archive.get(version)
        if plan:
            info = f"v{plan.frontmatter.version} ({plan.frontmatter.updated})"
        else:
            info = "(not found)"

        print(f"{prefix} {version} {info}")


def cmd_graph(graph: PlanGraph, format: str = "table") -> None:
    """Visualize the plan graph."""
    # Extract edges and plan info for visualize module
    edges = [e.to_dict() for e in graph.edges]
    plans = {}
    for plan_id, plan in graph.plans.items():
        plans[plan_id] = {
            'status': plan.frontmatter.status,
            'version': plan.frontmatter.version
        }

    output = visualize(edges, plans, format)
    print(output)


def cmd_validate(file_path: str) -> int:
    """Validate a plan file."""
    path = Path(file_path)

    if not path.exists():
        print(f"File not found: {path}")
        return 1

    plan = parse_plan_file(path)
    errors = plan.validate()

    if errors:
        print("Validation errors:")
        for e in errors:
            print(f"  - {e}")
        return 1

    print(f"Plan '{plan.id}' is valid")
    completed, total = plan.task_progress()
    print(f"  Status: {plan.frontmatter.status}")
    print(f"  Tasks: {completed}/{total} complete")
    print(f"  Pointers: {len(plan.get_all_pointers())}")
    return 0


def cmd_index_rebuild(plans_dir: Path) -> None:
    """Rebuild the index."""
    graph = PlanGraph(plans_dir)
    graph.build_graph()
    graph.save_index()

    print(f"Index rebuilt:")
    print(f"  Active plans: {len(graph.plans)}")
    print(f"  Archived: {len(graph.archive)}")
    print(f"  Edges: {len(graph.edges)}")
    print(f"  Clusters: {len(graph.clusters)}")


def cmd_knowledge_list(index: KnowledgeIndex) -> None:
    """List all knowledge entries by type."""
    index.scan_entries()
    grouped = index.list_all()

    total = sum(len(entries) for entries in grouped.values())
    if total == 0:
        print("No knowledge entries found.")
        print("Add entries to project/.knowledge/ with proper frontmatter.")
        return

    for entry_type, entries in grouped.items():
        if entries:
            print(f"\n{entry_type.upper()}S ({len(entries)}):")
            for entry in sorted(entries, key=lambda e: e.updated, reverse=True):
                tags = ', '.join(entry.tags[:3])
                severity = f" [{entry.severity}]" if entry.severity else ""
                print(f"  {entry.id}: {entry.title[:50]}{severity}")
                if tags:
                    print(f"    tags: {tags}")


def cmd_knowledge_search(index: KnowledgeIndex, query: str) -> None:
    """Search entries by tag or content."""
    index.scan_entries()
    results = index.search(query)

    if not results:
        print(f"No entries found matching '{query}'")
        return

    print(f"Found {len(results)} entries matching '{query}':")
    for entry in results:
        print(f"  {entry.id}: {entry.title} ({entry.type})")
        if entry.tags:
            print(f"    tags: {', '.join(entry.tags[:5])}")


def cmd_knowledge_add(index: KnowledgeIndex, file_path: str) -> int:
    """Index a new knowledge entry."""
    path = Path(file_path)

    if not path.exists():
        print(f"File not found: {path}")
        return 1

    entry = index.add_entry(path)

    if entry:
        print(f"Entry indexed: {entry.id}")
        print(f"  Type: {entry.type}")
        print(f"  Tags: {', '.join(entry.tags)}")
        print(f"  Path: {entry.path}")
        return 0
    else:
        print(f"Failed to parse entry from {path}")
        print("Ensure file has proper YAML frontmatter with: id, type, created, source")
        return 1


def cmd_knowledge_show(index: KnowledgeIndex, entry_id: str) -> int:
    """Display entry details."""
    index.scan_entries()
    entry = index.get_entry(entry_id)

    if not entry:
        print(f"Entry not found: {entry_id}")
        # Suggest similar entries
        index.scan_entries()
        similar = [e for e in index.entries.keys() if entry_id.lower() in e.lower()]
        if similar:
            print(f"Did you mean: {', '.join(similar[:5])}")
        return 1

    print(f"ID: {entry.id}")
    print(f"Type: {entry.type}")
    print(f"Source: {entry.source}")
    print(f"Path: {entry.path}")
    print(f"Tags: {', '.join(entry.tags)}")
    if entry.severity:
        print(f"Severity: {entry.severity}")
    if entry.zotero_key:
        print(f"Zotero: {entry.zotero_key}")
    print(f"Created: {entry.created}")
    print(f"Updated: {entry.updated}")
    print()
    print("--- Content Preview ---")
    # Show first 500 chars of content
    preview = entry.content[:500]
    if len(entry.content) > 500:
        preview += "\n..."
    print(preview)
    return 0


def cmd_knowledge_rebuild(index: KnowledgeIndex) -> None:
    """Rebuild the knowledge index."""
    count = index.rebuild_index()
    print(f"Knowledge index rebuilt: {count} entries")


# === Sequence Commands ===

def cmd_sequence_list(plans_dir: Path) -> None:
    """List all sequences with progress."""
    graph = PlanGraph(plans_dir)
    graph.build_graph()

    seq_manager = SequenceManager(plans_dir)
    seq_manager.load_all()

    if not seq_manager.sequences:
        print("No sequences found.")
        print(f"Add sequence definitions to {plans_dir / 'sequences'}/*.yaml")
        return

    plan_statuses = {pid: p.frontmatter.status for pid, p in graph.plans.items()}

    print("Sequences:")
    for name, seq in sorted(seq_manager.sequences.items()):
        completed, total = seq.compute_progress(plan_statuses)
        ready = seq.get_ready_plans(plan_statuses)

        status_icon = {"active": "●", "paused": "◌", "complete": "✓"}.get(seq.status, "?")
        print(f"  {status_icon} {name} ({seq.status}): {completed}/{total} plans")
        if ready:
            print(f"    Ready: {', '.join(ready)}")


def cmd_sequence_status(plans_dir: Path, name: str) -> int:
    """Show detailed sequence status with vision."""
    seq_manager = SequenceManager(plans_dir)
    seq = seq_manager.get_sequence(name)

    if not seq:
        print(f"Sequence not found: {name}")
        return 1

    # Use the sequence's own format_status which includes vision
    print(seq.format_status())
    return 0


def cmd_sequence_vision(plans_dir: Path, name: str) -> int:
    """Show vision section only (quick reference)."""
    seq_manager = SequenceManager(plans_dir)
    seq = seq_manager.get_sequence(name)

    if not seq:
        print(f"Sequence not found: {name}")
        return 1

    print(seq.vision.format_display())
    return 0


def cmd_sequence_next(plans_dir: Path, name: str, as_json: bool = False) -> int:
    """Get next ready plan in sequence (for auto-progression)."""
    graph = PlanGraph(plans_dir)
    graph.build_graph()

    seq_manager = SequenceManager(plans_dir)
    seq = seq_manager.get_sequence(name)

    if not seq:
        if as_json:
            print(json.dumps({"error": f"Sequence not found: {name}"}))
        else:
            print(f"Sequence not found: {name}")
        return 1

    plan_statuses = {pid: p.frontmatter.status for pid, p in graph.plans.items()}
    ready = seq.get_ready_plans(plan_statuses)

    if not ready:
        completed, total = seq.compute_progress(plan_statuses)
        if completed == total:
            reason = "sequence_complete"
        else:
            reason = "all_blocked_or_in_progress"

        if as_json:
            print(json.dumps({"next": None, "reason": reason}))
        else:
            if reason == "sequence_complete":
                print(f"Sequence '{name}' is complete.")
            else:
                print(f"No plans ready in sequence '{name}'.")
        return 0

    next_plan = ready[0]
    if as_json:
        print(json.dumps({"next": next_plan, "sequence": name}))
    else:
        print(f"Next: {next_plan}")

    return 0


def cmd_sequence_validate(plans_dir: Path, name: str) -> int:
    """Validate a sequence for errors."""
    graph = PlanGraph(plans_dir)
    graph.build_graph()

    seq_manager = SequenceManager(plans_dir)
    seq = seq_manager.get_sequence(name)

    if not seq:
        print(f"Sequence not found: {name}")
        return 1

    available_plans = set(graph.plans.keys())
    errors = seq.validate(available_plans)

    if errors:
        print(f"Validation errors for '{name}':")
        for e in errors:
            print(f"  - {e}")
        return 1

    print(f"Sequence '{name}' is valid.")
    print(f"  Plans: {len(seq.plans)}")
    print(f"  Status: {seq.status}")
    print(f"  Auto-progress: {seq.auto_progress}")
    return 0


def cmd_sequence_scope(plans_dir: Path, name: str) -> int:
    """Display sequence scope/vision with decisions for review."""
    seq_manager = SequenceManager(plans_dir)
    seq = seq_manager.get_sequence(name)

    if not seq:
        print(f"Sequence not found: {name}")
        return 1

    # Use vision format (works for both new vision and legacy scope)
    print(f"Sequence: {seq.name}")
    print("=" * 60)
    print()

    print("GOAL:")
    for line in seq.vision.end_goal.strip().split('\n'):
        print(f"  {line.strip()}")
    print()

    if seq.vision.why_it_matters:
        print("WHY IT MATTERS:")
        for line in seq.vision.why_it_matters.strip().split('\n'):
            print(f"  {line.strip()}")
        print()

    if seq.vision.current_state:
        print("CURRENT STATE:")
        for line in seq.vision.current_state.strip().split('\n'):
            print(f"  {line.strip()}")
        print()

    if seq.vision.decisions:
        print("DECISIONS:")
        for decision in seq.vision.decisions:
            print(f"  [{decision.date}] {decision.fork}")
            if decision.chose:
                print(f"    → Chose: {decision.chose}")
            if decision.rationale:
                for line in decision.rationale.strip().split('\n'):
                    print(f"    {line.strip()}")
            print()

    print("EXIT CRITERIA:")
    for i, criterion in enumerate(seq.vision.exit_criteria, 1):
        print(f"  {i}. {criterion}")

    return 0


# === Session Management Commands ===

def cmd_enter(plans_dir: Path, spec: str) -> int:
    """
    Enter a plan from a sequence.

    Spec format: sequence_name/plan_id or sequence_name/plan_short_name
    """
    # Parse spec
    if "/" not in spec:
        print(f"Error: Invalid spec '{spec}'. Use format: sequence/plan")
        return 1

    parts = spec.split("/", 1)
    seq_name = parts[0]
    plan_ref = parts[1]

    seq_manager = SequenceManager(plans_dir)
    seq = seq_manager.get_sequence(seq_name)

    if not seq:
        print(f"Sequence not found: {seq_name}")
        print("Available sequences:")
        seq_manager.load_all()
        for name in seq_manager.sequences.keys():
            print(f"  {name}")
        return 1

    # Find the plan
    plan = seq.get_plan(plan_ref)
    if not plan:
        # Try with area prefix
        for p in seq.plans:
            if p.short_name == plan_ref or p.id.endswith("/" + plan_ref):
                plan = p
                break

    if not plan:
        print(f"Plan not found in sequence: {plan_ref}")
        print("Available plans:")
        for p in seq.plans:
            print(f"  {p.short_name} ({p.status})")
        return 1

    # Check ownership
    session_id = seq_manager.generate_session_id()
    if plan.owner:
        print(f"Warning: Plan is owned by session {plan.owner}")
        response = _prompt("Take ownership anyway? [y/N]: ", "y").lower()
        if response != 'y':
            print("Aborted.")
            return 1

    # Claim the plan
    seq.claim_plan(plan.id, session_id)
    seq.save()

    # Create plan file if needed
    plan_path = seq_manager.create_plan_file(seq_name, plan.id)

    print(f"Entering: {seq.name} / {plan.short_name}")
    print()
    print("VISION REMINDER:")
    vision_line = seq.vision.end_goal.strip().split('\n')[0]
    print(f"  Goal: {vision_line}")
    print(f"  This plan: {plan.summary}")
    print()

    if plan.notes:
        print("PRIOR STATE:")
        print(f"  {plan.notes}")
        print()

    print(f"PLAN FILE: {plan_path}")
    print(f"SESSION ID: {session_id}")
    print("-" * 60)
    print()
    print("To pause: ./plan pause \"<notes>\"")
    print("To complete: ./plan seq-complete \"<summary>\"")

    # Save session info for pause/complete commands
    session_file = plans_dir / ".current_session"
    session_file.write_text(f"{seq_name}\n{plan.id}\n{session_id}\n")

    return 0


def cmd_pause(plans_dir: Path, notes: str) -> int:
    """Pause current plan work with notes."""
    session_file = plans_dir / ".current_session"

    if not session_file.exists():
        print("Error: No active session. Use './plan enter' first.")
        return 1

    lines = session_file.read_text().strip().split('\n')
    if len(lines) < 3:
        print("Error: Invalid session file.")
        return 1

    seq_name, plan_id, session_id = lines[0], lines[1], lines[2]

    seq_manager = SequenceManager(plans_dir)
    seq = seq_manager.get_sequence(seq_name)

    if not seq:
        print(f"Error: Sequence not found: {seq_name}")
        return 1

    # Update notes and release ownership
    seq.update_notes(plan_id, notes)
    seq.release_plan(plan_id)
    seq.save()

    # Clean up session file
    session_file.unlink()

    print(f"Paused: {plan_id}")
    print(f"  Notes saved: \"{notes[:50]}{'...' if len(notes) > 50 else ''}\"")
    print()
    print(f"Resume with: ./plan enter {seq_name}/{plan_id.split('/')[-1]}")

    return 0


def cmd_seq_complete(plans_dir: Path, summary: str) -> int:
    """Complete current plan in sequence and update sequence state."""
    session_file = plans_dir / ".current_session"

    if not session_file.exists():
        print("Error: No active session. Use './plan enter' first.")
        return 1

    lines = session_file.read_text().strip().split('\n')
    if len(lines) < 3:
        print("Error: Invalid session file.")
        return 1

    seq_name, plan_id, session_id = lines[0], lines[1], lines[2]

    seq_manager = SequenceManager(plans_dir)
    seq = seq_manager.get_sequence(seq_name)

    if not seq:
        print(f"Error: Sequence not found: {seq_name}")
        return 1

    # Mark plan complete and get unblocked plans
    unblocked = seq.complete_plan(plan_id, summary)
    seq.save()

    # Update plan .md frontmatter to match YAML
    area, name = plan_id.split("/", 1) if "/" in plan_id else ("", plan_id)
    plan_path = plans_dir / "active_plans" / area / f"{name}.md"
    if plan_path.exists():
        from datetime import date as _date
        updated = seq_manager.update_plan_frontmatter(
            plan_path, {"status": "complete", "updated": _date.today().isoformat()}
        )
        if updated:
            print(f"Updated plan frontmatter: {plan_path.name}")

    # Clean up session file
    session_file.unlink()

    print(f"Completed: {plan_id}")
    print(f"  Summary: \"{summary[:60]}{'...' if len(summary) > 60 else ''}\"")
    print()

    if unblocked:
        print("UNBLOCKED:")
        for p_id in unblocked:
            p = seq.get_plan(p_id)
            print(f"  {p_id} - {p.summary if p else ''}")
        print()

    # Check for next ready plan
    ready = seq.get_ready_plans()

    if not ready:
        completed, total = seq.compute_progress()
        if completed == total:
            print("=" * 60)
            print("SEQUENCE COMPLETE!")
            print("=" * 60)
            print()
            print("Exit criteria to verify:")
            for i, criterion in enumerate(seq.vision.exit_criteria, 1):
                print(f"  {i}. {criterion}")
        else:
            print("No plans currently ready (some may be in progress).")
    else:
        next_plan = ready[0]
        next_p = seq.get_plan(next_plan)

        # Check auto-progress setting
        if seq.auto_progress == "disabled":
            print(f"Next ready: {next_plan}")
        elif seq.auto_progress == "silent":
            print(f"Auto-entering next plan: {next_plan}")
            return cmd_enter(plans_dir, f"{seq_name}/{next_p.short_name}")
        else:  # prompt
            print(f"Next ready: {next_plan}")
            if next_p:
                print(f"  {next_p.summary}")
            response = _prompt("Enter this plan? [Y/n]: ", "n").lower()
            if response != 'n':
                return cmd_enter(plans_dir, f"{seq_name}/{next_p.short_name}")

    return 0


def _parse_session_frontmatter(path: Path) -> dict:
    """Parse YAML frontmatter from a session prompt file.

    Delegates to sequence.parse_session_frontmatter (shared utility).
    """
    return parse_session_frontmatter(path)


def _resolve_draft_path(explicit_path: str = None) -> str:
    """Resolve a draft plan path: explicit > .pending_sync.d/ breadcrumbs > legacy breadcrumb > most recent .md."""
    if explicit_path:
        return explicit_path

    plans_dir = Path.home() / ".claude" / "plans"

    # Check directory-based breadcrumbs (.pending_sync.d/)
    pending_dir = plans_dir / ".pending_sync.d"
    if pending_dir.is_dir():
        breadcrumbs = sorted(pending_dir.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
        for bc in breadcrumbs:
            if bc.is_file():
                path = bc.read_text().strip()
                if path and Path(path).expanduser().exists():
                    return path

    # Legacy: single .pending_sync file
    legacy_pending = plans_dir / ".pending_sync"
    if legacy_pending.exists():
        path = legacy_pending.read_text().strip()
        if path and Path(path).expanduser().exists():
            return path

    # Fallback: most recently modified .md in ~/.claude/plans/
    candidates = sorted(plans_dir.glob("*.md"), key=lambda p: p.stat().st_mtime, reverse=True)
    if candidates:
        return str(candidates[0])

    return None


def cmd_nudge(draft_path: str = None) -> int:
    """Write a sync breadcrumb so the next session start sees ACTION REQUIRED."""
    resolved = _resolve_draft_path(draft_path)
    if not resolved:
        print("No plan file found to nudge. Provide a path or have .md files in ~/.claude/plans/")
        return 1

    pending = Path.home() / ".claude" / "plans" / ".pending_sync"
    pending.parent.mkdir(parents=True, exist_ok=True)
    pending.write_text(resolved)
    print(f"Nudge written: {resolved}")
    print(f"Next session start will see: ACTION REQUIRED: sync this plan")
    return 0


def cmd_sync(plans_dir: Path, draft_path: str) -> int:
    """Sync a draft plan from ~/.claude/plans/ to active_plans and sequence YAML.

    1. Read draft, extract area/name from frontmatter
    2. Archive existing plan in active_plans (if any)
    3. Copy draft to active_plans/{area}/{name}.md
    4. Rebuild index
    5. Reconcile with sequence YAML
    """
    import shutil

    draft = Path(draft_path).expanduser()
    if not draft.exists():
        print(f"Draft plan not found: {draft}")
        return 1

    # Read and parse frontmatter
    try:
        text = draft.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError) as e:
        print(f"Cannot read draft: {e}")
        return 1

    if not text.startswith("---"):
        print(f"Draft has no YAML frontmatter (must start with ---)")
        return 1

    end = text.find("---", 3)
    if end == -1:
        print(f"Draft has unclosed frontmatter (missing closing ---)")
        return 1

    fm = yaml.safe_load(text[3:end]) or {}
    area = fm.get("area", "")
    name = fm.get("name", "")
    version = fm.get("version", "?")

    if not area or not name:
        print(f"Draft frontmatter missing 'area' or 'name' field")
        return 1

    plan_id = f"{area}/{name}"
    active_dir = plans_dir / "active_plans"
    target_dir = active_dir / area
    target_path = target_dir / f"{name}.md"
    archive_dir = plans_dir / ".archive"

    # Archive existing plan before overwriting
    if target_path.exists():
        old_text = target_path.read_text(encoding="utf-8")
        old_end = old_text.find("---", 3)
        old_fm = yaml.safe_load(old_text[3:old_end]) or {} if old_end != -1 else {}
        old_updated = old_fm.get("updated", "unknown")
        old_version = old_fm.get("version", 0)
        # Convert date objects to strings
        if hasattr(old_updated, "isoformat"):
            old_updated = old_updated.isoformat()
        archive_name = f"{old_updated}__{area}_{name}_v{old_version}.md"
        archive_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy(str(target_path), str(archive_dir / archive_name))
        print(f"Archived: {archive_name}")

    # Copy draft to active_plans
    target_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy(str(draft), str(target_path))
    print(f"Synced: {plan_id} (v{version})")

    # Rebuild index
    try:
        from .graph import PlanGraph
    except ImportError:
        graph_path = plans_dir / "tools" / "graph.py"
        if graph_path.exists():
            import subprocess
            subprocess.run(
                [sys.executable, str(graph_path), str(plans_dir), "--rebuild"],
                capture_output=True
            )
            print("Index rebuilt")
    else:
        graph = PlanGraph(plans_dir)
        graph.build_index()
        print("Index rebuilt")

    # Reconcile with sequence YAML
    manager = SequenceManager(plans_dir)
    result = manager.sync_plan_to_yaml(
        plan_id=plan_id,
        plan_status=fm.get("status", ""),
        plan_blocks=fm.get("blocks", []),
        plan_blocked_by=fm.get("blocked_by", [])
    )

    for msg in result.get("synced", []):
        print(f"  YAML: {msg}")
    for msg in result.get("conflicts", []):
        print(f"  CONFLICT: {msg}")
    for msg in result.get("warnings", []):
        print(f"  WARNING: {msg}")
    for msg in result.get("errors", []):
        print(f"  ERROR: {msg}")

    if not result.get("synced") and not result.get("conflicts") and not result.get("errors") and not result.get("warnings"):
        print(f"  Sequence YAML: already in sync")

    # Remove matching breadcrumb from .pending_sync.d/ (written by post_exit_plan.sh)
    pending_dir = Path.home() / ".claude" / "plans" / ".pending_sync.d"
    if pending_dir.is_dir():
        breadcrumb_name = f"{area}__{name}"
        breadcrumb_path = pending_dir / breadcrumb_name
        if breadcrumb_path.exists():
            try:
                breadcrumb_path.unlink()
                print(f"Cleared breadcrumb: {breadcrumb_name}")
            except OSError:
                pass
        # Also check for any breadcrumbs pointing to this draft path
        for bc in pending_dir.iterdir():
            if bc.is_file():
                try:
                    if bc.read_text().strip() == str(draft):
                        bc.unlink()
                except OSError:
                    pass

    # Remove legacy single-file breadcrumb if it exists
    legacy_pending = Path.home() / ".claude" / "plans" / ".pending_sync"
    if legacy_pending.exists():
        try:
            legacy_pending.unlink()
        except OSError:
            pass

    # Auto-reconcile: check for drift between plan files and sequence YAML
    try:
        containing = manager.find_sequences_containing(plan_id)
        if containing:
            quiet_conflicts = manager.detect_all_conflicts()
            if quiet_conflicts:
                print(f"\nReconcile found {len(quiet_conflicts)} conflict(s):")
                for c in quiet_conflicts[:5]:
                    print(f"  {c}")
                if len(quiet_conflicts) > 5:
                    print(f"  ... and {len(quiet_conflicts) - 5} more")
                print("Run: ./plan reconcile --fix <sequence> to resolve")
    except Exception:
        pass  # Non-fatal: reconcile is advisory

    return 0


def cmd_reconcile(plans_dir: Path, seq_name: str = None, plan_id: str = None, fix: bool = False) -> int:
    """Cross-reference plan files against sequence YAML, propagate status, flag conflicts.

    If plan_id given, reconciles just that plan.
    If seq_name given, reconciles all plans in that sequence.
    If neither, detects all conflicts across all sequences.

    With --fix: writes YAML-authoritative status back to plan .md frontmatter.
    """
    from datetime import date as _date

    manager = SequenceManager(plans_dir)

    # Full conflict detection mode
    if not seq_name and not plan_id:
        manager.load_all()
        conflicts = manager.detect_all_conflicts()
        if conflicts:
            print(f"Found {len(conflicts)} conflict(s):")
            for c in conflicts:
                print(f"  {c}")
            if fix:
                print("\n--fix requires a sequence name. Run: ./plan reconcile --fix <sequence>")
            return 1
        else:
            print("No conflicts detected.")
            return 0

    # Single plan sync
    if plan_id:
        # Read plan file frontmatter
        area, name = plan_id.split("/", 1) if "/" in plan_id else ("", plan_id)
        plan_path = plans_dir / "active_plans" / area / f"{name}.md"
        if not plan_path.exists():
            print(f"Plan file not found: {plan_path}")
            return 1

        try:
            text = plan_path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            print(f"Cannot read plan file: {plan_path}")
            return 1

        if not text.startswith("---"):
            print(f"Plan file has no frontmatter: {plan_path}")
            return 1

        end = text.find("---", 3)
        fm = yaml.safe_load(text[3:end]) or {} if end != -1 else {}

        result = manager.sync_plan_to_yaml(
            plan_id=plan_id,
            plan_status=fm.get("status", ""),
            plan_blocks=fm.get("blocks", []),
            plan_blocked_by=fm.get("blocked_by", [])
        )

        for msg in result["synced"]:
            print(f"  SYNCED: {msg}")
        for msg in result["conflicts"]:
            print(f"  CONFLICT: {msg}")
        for msg in result["errors"]:
            print(f"  ERROR: {msg}")

        if not result["synced"] and not result["conflicts"]:
            print(f"  {plan_id}: already in sync")

        return 1 if result["errors"] else 0

    # Sequence-wide reconcile (with optional --fix)
    if seq_name:
        seq = manager.get_sequence(seq_name)
        if not seq:
            print(f"Sequence not found: {seq_name}")
            return 1

        total_synced = 0
        total_conflicts = 0
        total_fixed = 0

        for plan in seq.plans:
            area, name = plan.id.split("/", 1) if "/" in plan.id else ("", plan.id)
            plan_path = plans_dir / "active_plans" / area / f"{name}.md"
            if not plan_path.exists():
                continue

            try:
                text = plan_path.read_text(encoding="utf-8")
            except (OSError, UnicodeDecodeError):
                continue

            if not text.startswith("---"):
                continue

            end = text.find("---", 3)
            fm = yaml.safe_load(text[3:end]) or {} if end != -1 else {}

            plan_status = fm.get("status", "")
            yaml_status = plan.status

            # --fix mode: YAML is authoritative for status
            if fix and plan_status and yaml_status and plan_status != yaml_status:
                updated = manager.update_plan_frontmatter(
                    plan_path,
                    {"status": yaml_status, "updated": _date.today().isoformat()}
                )
                if updated:
                    print(f"  FIXED: {plan.id}: status {plan_status} → {yaml_status}")
                    total_fixed += 1
                    # Re-read after fix so sync_plan_to_yaml sees consistent state
                    text = plan_path.read_text(encoding="utf-8")
                    end = text.find("---", 3)
                    fm = yaml.safe_load(text[3:end]) or {} if end != -1 else {}

            result = manager.sync_plan_to_yaml(
                plan_id=plan.id,
                plan_status=fm.get("status", ""),
                plan_blocks=fm.get("blocks", []),
                plan_blocked_by=fm.get("blocked_by", [])
            )

            for msg in result["synced"]:
                print(f"  SYNCED: {msg}")
                total_synced += 1
            for msg in result["conflicts"]:
                print(f"  CONFLICT: {msg}")
                total_conflicts += 1
            for msg in result["errors"]:
                print(f"  ERROR: {msg}")

        if total_synced == 0 and total_conflicts == 0 and total_fixed == 0:
            print(f"All plans in {seq_name} are in sync.")
        else:
            parts = []
            if total_fixed:
                parts.append(f"{total_fixed} fixed")
            if total_synced:
                parts.append(f"{total_synced} synced")
            if total_conflicts:
                parts.append(f"{total_conflicts} conflicts")
            print(f"\n{', '.join(parts)}")

        return 1 if total_conflicts > 0 else 0


def cmd_session_list(plans_dir: Path) -> int:
    """List session prompts with readiness status."""
    sessions_dir = plans_dir / "sessions"

    if not sessions_dir.exists():
        print("No sessions directory found.")
        print(f"  Expected: {sessions_dir}")
        return 0

    session_files = sorted(
        p for p in sessions_dir.rglob("*.md")
        if not p.name.startswith("._") and not p.name.startswith("00_") and not p.name.upper().startswith("README")
    )

    if not session_files:
        print("No session prompts found.")
        return 0

    seq_manager = SequenceManager(plans_dir)
    seq_manager.load_all()

    from collections import defaultdict
    by_area = defaultdict(list)

    for sf in session_files:
        rel = sf.relative_to(sessions_dir)
        area = rel.parts[0] if len(rel.parts) > 1 else "(ungrouped)"
        by_area[area].append((sf, rel))

    print("Session Prompts:")
    print()

    for area, files in by_area.items():
        print(f"  {area}/")
        for sf, rel in files:
            fm = _parse_session_frontmatter(sf)
            seq_name = fm.get("sequence", "")
            plan_id = fm.get("plan_id", "")
            mode = fm.get("mode", "?")

            status_icon = " "
            if seq_name and plan_id:
                seq = seq_manager.get_sequence(seq_name)
                if seq:
                    ready_plans = seq.get_ready_plans()
                    plan_obj = seq.get_plan(plan_id)
                    if plan_obj and plan_obj.status == "complete":
                        status_icon = "\u2713"
                    elif plan_id in ready_plans:
                        status_icon = "\u25cf"
                    else:
                        status_icon = "\u25cb"
                else:
                    status_icon = "?"

            display = str(rel.with_suffix(""))
            print(f"    {status_icon} {display}  [{mode}]")

        print()

    print("Legend: \u25cf ready  \u25cb blocked  \u2713 complete  ? no sequence")
    return 0


def cmd_session_generate(plans_dir: Path, seq_name: str, dry_run: bool = False) -> int:
    """Generate plan stubs + session prompts from a sequence."""
    seq_manager = SequenceManager(plans_dir)
    result = seq_manager.generate_session_skeletons(seq_name, dry_run=dry_run)

    if result["errors"] and not result["generated"] and not result["skipped"]:
        for plan_id, error in result["errors"]:
            print(f"Error: {error}")
        return 1

    prefix = "[DRY RUN] " if dry_run else ""

    if result["generated"]:
        print(f"{prefix}Generated {len(result['generated'])} plan stubs + session prompts:")
        for plan_id, stub, session in result["generated"]:
            parts = []
            if stub and stub != "(exists)":
                parts.append(stub)
            if session and session != "(exists)":
                parts.append(session)
            print(f"  + {plan_id}")
            for p in parts:
                print(f"    → {p}")

    if result["skipped"]:
        print(f"\n{prefix}Skipped {len(result['skipped'])} (already exist):")
        for plan_id, reason in result["skipped"]:
            print(f"  - {plan_id}: {reason}")

    if result["errors"]:
        print(f"\n{prefix}Errors:")
        for plan_id, error in result["errors"]:
            print(f"  ! {plan_id}: {error}")

    return 0


def cmd_session_show(plans_dir: Path, session_path: str) -> int:
    """Display a session prompt."""
    sessions_dir = plans_dir / "sessions"

    candidates = [
        sessions_dir / session_path,
        sessions_dir / f"{session_path}.md",
    ]

    target = None
    for c in candidates:
        if c.exists() and c.is_file():
            target = c
            break

    if not target:
        print(f"Error: Session prompt not found: {session_path}")
        print(f"  Searched: {sessions_dir}")
        all_sessions = list(sessions_dir.rglob("*.md"))
        search_term = Path(session_path).stem.lower()
        close = [
            str(s.relative_to(sessions_dir))
            for s in all_sessions
            if search_term in s.stem.lower()
        ]
        if close:
            print(f"  Did you mean: {', '.join(close[:5])}")
        return 1

    print(target.read_text())
    return 0


def cmd_complete(plans_dir: Path, plan_id: str) -> int:
    """Run completion debrief and archive a plan."""
    import shutil
    from datetime import date

    # Find the plan file
    active_dir = plans_dir / "active_plans"
    plan_file = None

    # Search for plan by area/name or just name
    for md_file in active_dir.rglob("*.md"):
        if md_file.name.startswith("._"):
            continue
        plan = parse_plan_file(md_file)
        plan_id_full = f"{plan.frontmatter.area}/{plan.frontmatter.name}"
        if plan_id in (plan.frontmatter.name, plan_id_full, plan.id):
            plan_file = md_file
            break

    if not plan_file:
        print(f"Plan not found: {plan_id}")
        print("Available plans:")
        for md_file in active_dir.rglob("*.md"):
            if not md_file.name.startswith("._"):
                p = parse_plan_file(md_file)
                print(f"  {p.frontmatter.area}/{p.frontmatter.name}")
        return 1

    plan = parse_plan_file(plan_file)
    completed, total = plan.task_progress()

    print("=" * 60)
    print(f"COMPLETION DEBRIEF: {plan.frontmatter.area}/{plan.frontmatter.name}")
    print("=" * 60)
    print()

    # Task check
    print(f"Tasks: {completed}/{total} complete")
    if completed < total:
        print()
        print("WARNING: Not all tasks are complete!")
        print("Incomplete tasks:")
        for task in plan.sections.get("Tasks", []):
            if isinstance(task, str) and task.strip().startswith("- [ ]"):
                print(f"  {task.strip()}")
        print()
        response = _prompt("Continue anyway? [y/N]: ", "y").lower()
        if response != 'y':
            print("Aborted.")
            return 1

    print()
    print("-" * 60)
    print("KNOWLEDGE DEBRIEF")
    print("-" * 60)
    print()
    print("Answer the following to capture lessons learned.")
    print("Press Enter to skip if not applicable.")
    print("PREFER UPDATING EXISTING DOCS over creating new ones.")
    print()

    # Debrief questions
    questions = [
        ("outdated", "Any knowledge docs with OUTDATED information to UPDATE?\n  (e.g., 'da_exceptions_registry - add new Pal finding')"),
        ("hpc", "Any new HPC conventions or gotchas learned?\n  (e.g., 'container X needs Y bind mount')"),
        ("errors", "Any errors encountered FREQUENTLY worth documenting?\n  (e.g., 'CRSP stale handle on large writes')"),
        ("searches", "Any web searches that yielded useful info to CONSOLIDATE?\n  (e.g., 'MiloR dispersion estimation - add to milo_differential_abundance.md')"),
    ]

    responses = {}
    for key, q in questions:
        print(f"\n{q}")
        response = _prompt("> ", "").strip()
        if response:
            responses[key] = response

    print()
    print("-" * 60)
    print("SUMMARY")
    print("-" * 60)

    if responses:
        print("\nKnowledge actions to take:")
        for key, response in responses.items():
            print(f"  [{key}] {response}")
        print()
        print("Run these updates BEFORE confirming archive:")
        print("  ./plan knowledge list   # Find existing entries")
        print("  Edit the relevant .knowledge/ files directly")
        print()
    else:
        print("\nNo knowledge updates noted.")

    # Confirm archive
    print()
    response = _prompt("Archive this plan as COMPLETE? [y/N]: ", "y").lower()
    if response != 'y':
        print("Aborted. Plan remains active.")
        return 1

    # Archive the plan
    archive_dir = plans_dir / ".archive"
    archive_dir.mkdir(exist_ok=True)

    today = date.today().isoformat()
    archive_name = f"{today}__{plan.frontmatter.area}_{plan.frontmatter.name}_v{plan.frontmatter.version}.md"
    archive_path = archive_dir / archive_name

    # Update status to complete before archiving
    content = plan_file.read_text()
    content = content.replace("status: in_progress", "status: complete")
    content = content.replace("status: queued", "status: complete")
    archive_path.write_text(content)

    # Remove from active
    plan_file.unlink()

    # Clean up empty parent dirs
    parent = plan_file.parent
    if parent != active_dir and not any(parent.iterdir()):
        parent.rmdir()

    print()
    print(f"Archived: {archive_path.name}")
    print("Plan complete.")

    # Archive superseded plan if it exists in active_plans
    if plan.frontmatter.supersedes:
        sup_id = plan.frontmatter.supersedes
        sup_area, sup_name = sup_id.split("/", 1) if "/" in sup_id else ("", sup_id)
        sup_path = active_dir / sup_area / f"{sup_name}.md"
        if sup_path.exists():
            sup_plan = parse_plan_file(sup_path)
            sup_archive_name = f"{today}__{sup_area}_{sup_name}_v{sup_plan.frontmatter.version}_superseded.md"
            sup_archive_path = archive_dir / sup_archive_name
            sup_content = sup_path.read_text()
            sup_content = sup_content.replace("status: in_progress", "status: complete")
            sup_content = sup_content.replace("status: queued", "status: complete")
            sup_content = sup_content.replace("status: blocked", "status: complete")
            sup_archive_path.write_text(sup_content)
            sup_path.unlink()
            sup_parent = sup_path.parent
            if sup_parent != active_dir and not any(sup_parent.iterdir()):
                sup_parent.rmdir()
            print(f"Superseded plan archived: {sup_archive_name}")

    # Rebuild index
    graph = PlanGraph(plans_dir)
    graph.build_graph()
    graph.save_index()

    return 0


def main():
    parser = argparse.ArgumentParser(description="Plan query interface")
    subparsers = parser.add_subparsers(dest="command", help="Command")

    # Query subcommand
    query_parser = subparsers.add_parser("query", help="Query the plan graph")
    query_parser.add_argument("subcommand", choices=[
        "ready", "blocked", "context", "downstream",
        "touches", "anchors", "communities", "history", "tasks"
    ])
    query_parser.add_argument("argument", nargs="?", help="Argument for subcommand")
    query_parser.add_argument("--all", "-a", action="store_true", help="Include completed items")

    # Graph subcommand
    graph_parser = subparsers.add_parser("graph", help="Visualize plan graph")
    graph_parser.add_argument("--ascii", action="store_true", help="ASCII DAG view")
    graph_parser.add_argument("--mermaid", action="store_true", help="Mermaid markdown syntax")
    graph_parser.add_argument("--dot", action="store_true", help="Graphviz DOT syntax")

    # Validate subcommand
    validate_parser = subparsers.add_parser("validate", help="Validate a plan file")
    validate_parser.add_argument("file", help="Plan file to validate")

    # Index subcommand
    index_parser = subparsers.add_parser("index", help="Index operations")
    index_parser.add_argument("subcommand", choices=["rebuild"])

    # Knowledge subcommand
    knowledge_parser = subparsers.add_parser("knowledge", help="Knowledge base operations")
    knowledge_parser.add_argument("subcommand", choices=["list", "search", "add", "show", "rebuild"])
    knowledge_parser.add_argument("argument", nargs="?", help="Query, path, or entry ID")

    # Complete subcommand
    complete_parser = subparsers.add_parser("complete", help="Run completion debrief and archive")
    complete_parser.add_argument("plan_id", help="Plan ID to complete (name or area/name)")

    # Sync subcommand (draft → active_plans → index → YAML)
    sync_parser = subparsers.add_parser("sync", help="Sync a draft plan to active_plans and sequence YAML")
    sync_parser.add_argument("draft_path", nargs="?", default=None, help="Path to draft plan file (omit to auto-detect from .pending_sync or most recent)")

    # Nudge subcommand (write breadcrumb for next session)
    nudge_parser = subparsers.add_parser("nudge", help="Write a sync reminder for the next session start")
    nudge_parser.add_argument("draft_path", nargs="?", default=None, help="Path to plan file (omit to auto-detect most recent)")

    # Reconcile subcommand (cross-reference plan files ↔ YAML)
    reconcile_parser = subparsers.add_parser("reconcile", help="Cross-reference plan files against sequence YAML")
    reconcile_parser.add_argument("argument", nargs="?", help="Sequence name or plan ID (omit for full conflict check)")
    reconcile_parser.add_argument("--plan", help="Reconcile a specific plan ID")
    reconcile_parser.add_argument("--fix", action="store_true", help="Fix status mismatches by writing YAML status back to plan .md frontmatter")

    # Sequence subcommand
    sequence_parser = subparsers.add_parser("sequence", help="Sequence operations")
    sequence_parser.add_argument("subcommand", choices=["list", "status", "next", "validate", "scope", "vision"])
    sequence_parser.add_argument("argument", nargs="?", help="Sequence name")
    sequence_parser.add_argument("--json", action="store_true", help="Output as JSON (for next command)")

    # Enter subcommand (session management)
    enter_parser = subparsers.add_parser("enter", help="Enter a plan from a sequence")
    enter_parser.add_argument("spec", help="Sequence/plan spec (e.g., da_pipeline/milo_rebuild)")

    # Pause subcommand (session management)
    pause_parser = subparsers.add_parser("pause", help="Pause current plan with notes")
    pause_parser.add_argument("notes", nargs="?", default="", help="Progress notes for next session")

    # Seq-complete subcommand (session management)
    seq_complete_parser = subparsers.add_parser("seq-complete", help="Complete current plan in sequence")
    seq_complete_parser.add_argument("summary", help="Completion summary")

    # Session subcommand
    session_parser = subparsers.add_parser("session", help="Session prompt operations")
    session_parser.add_argument("subcommand", choices=["list", "show", "generate"])
    session_parser.add_argument("argument", nargs="?", help="Session path (for show) or sequence name (for generate)")
    session_parser.add_argument("--dry-run", action="store_true", help="Show what would be generated without creating files")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    try:
        plans_dir = get_plans_dir()
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    if args.command == "validate":
        sys.exit(cmd_validate(args.file))

    if args.command == "graph":
        graph = PlanGraph(plans_dir)
        graph.build_graph()
        # Determine format from flags
        if args.ascii:
            fmt = "ascii"
        elif args.mermaid:
            fmt = "mermaid"
        elif args.dot:
            fmt = "dot"
        else:
            fmt = "table"
        cmd_graph(graph, fmt)
        sys.exit(0)

    if args.command == "index":
        if args.subcommand == "rebuild":
            cmd_index_rebuild(plans_dir)
        sys.exit(0)

    if args.command == "knowledge":
        try:
            knowledge_dir = get_knowledge_dir()
        except FileNotFoundError:
            print("Error: Could not find .knowledge directory")
            print("Create project/.knowledge/ with entries in lessons/, research/, patterns/")
            sys.exit(1)

        index = KnowledgeIndex(knowledge_dir)

        if args.subcommand == "list":
            cmd_knowledge_list(index)
        elif args.subcommand == "search":
            if not args.argument:
                print("Error: search requires a query")
                sys.exit(1)
            cmd_knowledge_search(index, args.argument)
        elif args.subcommand == "add":
            if not args.argument:
                print("Error: add requires a file path")
                sys.exit(1)
            sys.exit(cmd_knowledge_add(index, args.argument))
        elif args.subcommand == "show":
            if not args.argument:
                print("Error: show requires an entry ID")
                sys.exit(1)
            sys.exit(cmd_knowledge_show(index, args.argument))
        elif args.subcommand == "rebuild":
            cmd_knowledge_rebuild(index)
        sys.exit(0)

    if args.command == "complete":
        sys.exit(cmd_complete(plans_dir, args.plan_id))

    if args.command == "nudge":
        sys.exit(cmd_nudge(getattr(args, 'draft_path', None)))

    if args.command == "sync":
        resolved = _resolve_draft_path(getattr(args, 'draft_path', None))
        if not resolved:
            print("No plan file found. Provide a path, or have .md files in ~/.claude/plans/")
            sys.exit(1)
        sys.exit(cmd_sync(plans_dir, resolved))

    if args.command == "reconcile":
        plan_id = getattr(args, 'plan', None)
        seq_name = args.argument if not plan_id else None
        fix = getattr(args, 'fix', False)
        sys.exit(cmd_reconcile(plans_dir, seq_name=seq_name, plan_id=plan_id, fix=fix))

    if args.command == "sequence":
        if args.subcommand == "list":
            cmd_sequence_list(plans_dir)
        elif args.subcommand == "status":
            if not args.argument:
                print("Error: status requires a sequence name")
                sys.exit(1)
            sys.exit(cmd_sequence_status(plans_dir, args.argument))
        elif args.subcommand == "next":
            if not args.argument:
                print("Error: next requires a sequence name")
                sys.exit(1)
            sys.exit(cmd_sequence_next(plans_dir, args.argument, as_json=args.json))
        elif args.subcommand == "validate":
            if not args.argument:
                print("Error: validate requires a sequence name")
                sys.exit(1)
            sys.exit(cmd_sequence_validate(plans_dir, args.argument))
        elif args.subcommand == "scope":
            if not args.argument:
                print("Error: scope requires a sequence name")
                sys.exit(1)
            sys.exit(cmd_sequence_scope(plans_dir, args.argument))
        elif args.subcommand == "vision":
            if not args.argument:
                print("Error: vision requires a sequence name")
                sys.exit(1)
            sys.exit(cmd_sequence_vision(plans_dir, args.argument))
        sys.exit(0)

    if args.command == "enter":
        sys.exit(cmd_enter(plans_dir, args.spec))

    if args.command == "pause":
        sys.exit(cmd_pause(plans_dir, args.notes))

    if args.command == "seq-complete":
        sys.exit(cmd_seq_complete(plans_dir, args.summary))

    if args.command == "session":
        if args.subcommand == "list":
            sys.exit(cmd_session_list(plans_dir))
        elif args.subcommand == "show":
            if not args.argument:
                print("Error: show requires a session path")
                print("Usage: plan session show <area/unit/session>")
                sys.exit(1)
            sys.exit(cmd_session_show(plans_dir, args.argument))
        elif args.subcommand == "generate":
            if not args.argument:
                print("Error: generate requires a sequence name")
                print("Usage: plan session generate <sequence> [--dry-run]")
                sys.exit(1)
            sys.exit(cmd_session_generate(plans_dir, args.argument, dry_run=args.dry_run))

    if args.command == "query":
        graph = PlanGraph(plans_dir)
        graph.build_graph()

        if args.subcommand == "ready":
            cmd_query_ready(graph)
        elif args.subcommand == "blocked":
            cmd_query_blocked(graph)
        elif args.subcommand == "context":
            if not args.argument:
                print("Error: context requires a plan ID")
                sys.exit(1)
            cmd_query_context(graph, args.argument)
        elif args.subcommand == "downstream":
            if not args.argument:
                print("Error: downstream requires a plan ID")
                sys.exit(1)
            cmd_query_downstream(graph, args.argument)
        elif args.subcommand == "touches":
            if not args.argument:
                print("Error: touches requires a file pattern")
                sys.exit(1)
            cmd_query_touches(graph, args.argument)
        elif args.subcommand == "anchors":
            cmd_query_anchors(graph)
        elif args.subcommand == "communities":
            cmd_query_communities(graph)
        elif args.subcommand == "history":
            if not args.argument:
                print("Error: history requires a plan ID")
                sys.exit(1)
            cmd_query_history(graph, args.argument)
        elif args.subcommand == "tasks":
            cmd_query_tasks(graph, show_complete=args.all)


if __name__ == "__main__":
    main()
