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
"""

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Optional

# Add tools directory to path for imports
TOOLS_DIR = Path(__file__).parent
sys.path.insert(0, str(TOOLS_DIR))

from parse import parse_plan_file, parse_plan
from graph import PlanGraph
from visualize import visualize
from knowledge import KnowledgeIndex, get_knowledge_dir


# Yet another CLI wrapper. At least this one queries something useful.


def get_plans_dir(explicit_path: Optional[Path] = None) -> Path:
    """
    Find the .plans directory.

    Resolution order:
    1. Explicit path argument (highest priority)
    2. PLAN_DIR environment variable
    3. CWD-based discovery (.plans or project/.plans)
    4. TOOLS_DIR parent traversal (original fallback)
    """
    # 1. Explicit argument
    if explicit_path:
        if explicit_path.exists():
            return explicit_path
        raise FileNotFoundError(f"Specified plans directory not found: {explicit_path}")

    # 2. Environment variable
    env_path = os.environ.get('PLAN_DIR')
    if env_path:
        path = Path(env_path)
        if path.exists():
            return path
        raise FileNotFoundError(f"PLAN_DIR not found: {env_path}")

    # 3. CWD-based discovery
    cwd = Path.cwd()
    for candidate in [cwd / ".plans", cwd / "project" / ".plans"]:
        if candidate.exists():
            return candidate

    # 4. TOOLS_DIR-based discovery (original fallback for project-embedded install)
    current = TOOLS_DIR.parent
    if current.name == ".plans":
        return current

    for parent in [current] + list(current.parents):
        for subpath in ["project/.plans", ".plans"]:
            plans_dir = parent / subpath
            if plans_dir.exists():
                return plans_dir

    raise FileNotFoundError(
        "Could not find .plans directory.\n"
        "Options:\n"
        "  1. Run from project root containing .plans/ or project/.plans/\n"
        "  2. Set PLAN_DIR environment variable\n"
        "  3. Use --plans-dir argument"
    )


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
        print("Add entries to .knowledge/ with proper frontmatter.")
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
        response = input("Continue anyway? [y/N]: ").strip().lower()
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

    # Debrief questions - generic version for any project
    questions = [
        ("outdated", "Any knowledge docs with OUTDATED information to UPDATE?\n  (e.g., 'config_reference - add new option')"),
        ("gotchas", "Any new gotchas or edge cases discovered?\n  (e.g., 'module X requires Y version')"),
        ("errors", "Any errors encountered FREQUENTLY worth documenting?\n  (e.g., 'timeout on large datasets')"),
        ("searches", "Any web searches that yielded useful info to CONSOLIDATE?\n  (e.g., 'API docs for X - add to reference.md')"),
    ]

    responses = {}
    for key, prompt in questions:
        print(f"\n{prompt}")
        response = input("> ").strip()
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
    response = input("Archive this plan as COMPLETE? [y/N]: ").strip().lower()
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

    # Rebuild index
    graph = PlanGraph(plans_dir)
    graph.build_graph()
    graph.save_index()

    return 0


def main():
    parser = argparse.ArgumentParser(description="Plan query interface")
    parser.add_argument('--plans-dir', type=Path,
                        help='Path to .plans directory (default: auto-discover)')
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

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(1)

    # Validate doesn't need plans_dir - handle it early
    if args.command == "validate":
        sys.exit(cmd_validate(args.file))

    # Get plans directory (may use --plans-dir if provided)
    explicit_plans_dir = getattr(args, 'plans_dir', None)

    try:
        plans_dir = get_plans_dir(explicit_plans_dir)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

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
            print("Create .knowledge/ with entries in lessons/, research/, patterns/")
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
