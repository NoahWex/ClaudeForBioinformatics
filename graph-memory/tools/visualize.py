#!/usr/bin/env python3
"""
Plan graph visualization - Multiple output formats for dependency graphs.

Formats:
- Rich table: Terminal quick-view (default)
- ASCII DAG: Hierarchical tree with arrows
- Mermaid: GitHub/browser rendering
- Graphviz DOT: Publication-quality images
"""

from pathlib import Path
from typing import Optional
from collections import defaultdict

# Why did the graph go to therapy?
# It had too many issues with its edges.


def render_rich_table(edges: list[dict], plans: dict) -> str:
    """Render edges as a Rich table."""
    try:
        from rich.console import Console
        from rich.table import Table

        table = Table(title="Plan Dependencies")
        table.add_column("Source", style="cyan")
        table.add_column("Target", style="green")
        table.add_column("Type", style="yellow")
        table.add_column("Pattern", style="dim")

        for edge in sorted(edges, key=lambda e: (e['from'], e['to'])):
            pattern = edge.get('pattern', '')
            if len(pattern) > 40:
                pattern = '...' + pattern[-37:]
            table.add_row(
                edge['from'],
                edge['to'],
                edge['type'],
                pattern
            )

        console = Console(record=True)
        console.print(table)
        return console.export_text()
    except ImportError:
        # Fallback without rich
        return render_ascii_table(edges, plans)


def render_ascii_table(edges: list[dict], plans: dict) -> str:
    """Simple ASCII table fallback."""
    lines = ["Plan Dependencies", "=" * 80]
    lines.append(f"{'Source':<30} {'Target':<25} {'Type':<12}")
    lines.append("-" * 80)

    for edge in sorted(edges, key=lambda e: (e['from'], e['to'])):
        lines.append(f"{edge['from']:<30} {edge['to']:<25} {edge['type']:<12}")

    return "\n".join(lines)


def render_ascii_dag(edges: list[dict], plans: dict) -> str:
    """
    Render as ASCII DAG with topological ordering.
    Shows dependency flow with arrows.
    """
    # Build adjacency for topological sort
    graph = defaultdict(set)
    all_nodes = set()

    for edge in edges:
        if edge['type'] in ('dependency', 'dataflow'):
            graph[edge['from']].add(edge['to'])
            all_nodes.add(edge['from'])
            all_nodes.add(edge['to'])

    # Also add nodes with no edges
    for plan_id in plans:
        all_nodes.add(plan_id)

    # Compute in-degree for topological sort
    in_degree = defaultdict(int)
    for node in all_nodes:
        in_degree[node] = 0
    for source, targets in graph.items():
        for target in targets:
            in_degree[target] += 1

    # Kahn's algorithm for topological sort
    queue = [n for n in all_nodes if in_degree[n] == 0]
    queue.sort()  # Deterministic ordering
    sorted_nodes = []

    while queue:
        node = queue.pop(0)
        sorted_nodes.append(node)
        for target in sorted(graph[node]):
            in_degree[target] -= 1
            if in_degree[target] == 0:
                queue.append(target)
                queue.sort()

    # Handle cycles (nodes not in sorted list)
    remaining = set(all_nodes) - set(sorted_nodes)
    sorted_nodes.extend(sorted(remaining))

    # Build output with indentation based on depth
    depth = {n: 0 for n in all_nodes}
    for node in sorted_nodes:
        for target in graph[node]:
            depth[target] = max(depth[target], depth[node] + 1)

    lines = ["Plan Dependency Graph", "=" * 60, ""]

    for node in sorted_nodes:
        indent = "  " * depth[node]
        targets = graph.get(node, set())
        status = plans.get(node, {}).get('status', 'unknown')
        status_icon = {'complete': '✓', 'in_progress': '●', 'blocked': '⊗', 'queued': '○'}.get(status, '?')

        if targets:
            arrow = " →→ " + ", ".join(sorted(targets))
        else:
            arrow = ""

        lines.append(f"{indent}{status_icon} {node}{arrow}")

    return "\n".join(lines)


def render_mermaid(edges: list[dict], plans: dict) -> str:
    """
    Render as Mermaid flowchart syntax.
    Can be embedded in markdown for GitHub rendering.
    """
    lines = ["```mermaid", "graph LR"]

    # Define nodes with sanitized IDs
    node_ids = {}
    for plan_id in plans:
        safe_id = plan_id.replace('/', '_').replace('-', '_').replace('.', '_')
        node_ids[plan_id] = safe_id
        lines.append(f'    {safe_id}["{plan_id}"]')

    lines.append("")

    # Define edges with labels
    edge_styles = {
        'dependency': '-->',
        'dataflow': '-.->',
        'resource': '-.-',
        'succession': '==>'
    }

    seen_edges = set()
    for edge in edges:
        source = node_ids.get(edge['from'], edge['from'].replace('/', '_'))
        target = node_ids.get(edge['to'], edge['to'].replace('/', '_'))
        edge_type = edge['type']
        arrow = edge_styles.get(edge_type, '-->')

        edge_key = (source, target, edge_type)
        if edge_key in seen_edges:
            continue
        seen_edges.add(edge_key)

        lines.append(f'    {source} {arrow}|{edge_type}| {target}')

    lines.append("")

    # Style nodes by status
    status_colors = {
        'complete': '#c8e6c9',      # Green
        'in_progress': '#bbdefb',   # Blue
        'blocked': '#ffccbc',       # Orange
        'queued': '#f5f5f5'         # Gray
    }

    for plan_id, info in plans.items():
        safe_id = node_ids.get(plan_id, plan_id.replace('/', '_'))
        status = info.get('status', 'queued')
        color = status_colors.get(status, '#f5f5f5')
        lines.append(f'    style {safe_id} fill:{color}')

    lines.append("```")
    return "\n".join(lines)


def render_graphviz_dot(edges: list[dict], plans: dict) -> str:
    """
    Render as Graphviz DOT syntax.
    Use with: dot -Tpng graph.dot > graph.png
    """
    lines = [
        "digraph PlanGraph {",
        "    rankdir=LR;",
        "    node [shape=box, style=filled, fontname=\"Helvetica\"];",
        "    edge [fontname=\"Helvetica\", fontsize=10];",
        ""
    ]

    # Status to color mapping
    status_colors = {
        'complete': 'lightgreen',
        'in_progress': 'lightblue',
        'blocked': 'lightsalmon',
        'queued': 'lightgray'
    }

    # Define nodes
    for plan_id, info in plans.items():
        safe_id = plan_id.replace('/', '_').replace('-', '_').replace('.', '_')
        status = info.get('status', 'queued')
        color = status_colors.get(status, 'lightgray')
        label = plan_id.replace('_', '\\n')  # Multi-line labels
        lines.append(f'    {safe_id} [label="{plan_id}", fillcolor={color}];')

    lines.append("")

    # Define edges
    edge_styles = {
        'dependency': 'solid',
        'dataflow': 'dashed',
        'resource': 'dotted',
        'succession': 'bold'
    }

    seen_edges = set()
    for edge in edges:
        source = edge['from'].replace('/', '_').replace('-', '_').replace('.', '_')
        target = edge['to'].replace('/', '_').replace('-', '_').replace('.', '_')
        edge_type = edge['type']
        style = edge_styles.get(edge_type, 'solid')

        edge_key = (source, target, edge_type)
        if edge_key in seen_edges:
            continue
        seen_edges.add(edge_key)

        lines.append(f'    {source} -> {target} [label="{edge_type}", style={style}];')

    lines.append("}")
    return "\n".join(lines)


def visualize(
    edges: list[dict],
    plans: dict,
    format: str = "table"
) -> str:
    """
    Render plan graph in specified format.

    Args:
        edges: List of edge dicts with 'from', 'to', 'type' keys
        plans: Dict of plan_id -> plan info with 'status' key
        format: One of 'table', 'ascii', 'mermaid', 'dot'

    Returns:
        Formatted string output
    """
    renderers = {
        'table': render_rich_table,
        'ascii': render_ascii_dag,
        'mermaid': render_mermaid,
        'dot': render_graphviz_dot
    }

    renderer = renderers.get(format, render_rich_table)
    return renderer(edges, plans)


if __name__ == "__main__":
    import sys
    import json

    # Example usage with sample data
    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        sample_edges = [
            {"from": "04_DA/milo_rebuild", "to": "06_Xenium/scaffold", "type": "dependency"},
            {"from": "04_DA/milo_rebuild", "to": "04_DA/da_testing", "type": "dependency"},
            {"from": "04_DA/milo_rebuild", "to": "06_Xenium/scaffold", "type": "dataflow", "pattern": "{study}_milo.rds"},
            {"from": "_Root/GOALS", "to": "04_DA/milo_rebuild", "type": "resource", "pattern": "seurat.rds"}
        ]
        sample_plans = {
            "_Root/GOALS": {"status": "in_progress"},
            "04_DA/milo_rebuild": {"status": "in_progress"},
            "04_DA/da_testing": {"status": "blocked"},
            "06_Xenium/scaffold": {"status": "queued"}
        }

        fmt = sys.argv[2] if len(sys.argv) > 2 else "table"
        print(visualize(sample_edges, sample_plans, fmt))
    else:
        print("Usage: visualize.py --demo [table|ascii|mermaid|dot]")
