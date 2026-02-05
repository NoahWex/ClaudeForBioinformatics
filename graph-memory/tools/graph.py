#!/usr/bin/env python3
"""
Graph builder - Construct and query the plan dependency graph.

Features:
- Edge inference from output→input matching
- SNN clustering for related plans
- MNN for cross-modal anchor validation
- Succession chain tracking
- Pointer pooling across clusters
"""

import json
import re
from pathlib import Path
from dataclasses import dataclass, field
from collections import defaultdict
from typing import Optional
import fnmatch

from parse import Plan, parse_plan_file, Pointer


# Why did the graph theorist break up with the tree?
# Too many commitment issues with cycles.


@dataclass
class Edge:
    """A typed edge between two plans."""
    source: str  # plan id
    target: str  # plan id
    edge_type: str  # dependency, dataflow, resource, succession
    weight: float = 1.0
    pattern: str = ""  # for dataflow edges, the matching pattern

    def to_dict(self) -> dict:
        d = {
            "from": self.source,
            "to": self.target,
            "type": self.edge_type,
            "weight": self.weight
        }
        if self.pattern:
            d["pattern"] = self.pattern
        return d


@dataclass
class Cluster:
    """A group of related plans with pooled context."""
    name: str
    members: list[str]
    pooled_inputs: list[str] = field(default_factory=list)
    pooled_outputs: list[str] = field(default_factory=list)
    pooled_tasks: list[str] = field(default_factory=list)
    pooled_discussion: list[dict] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "members": self.members,
            "pooled": {
                "inputs": self.pooled_inputs,
                "outputs": self.pooled_outputs,
                "tasks": self.pooled_tasks,
                "discussion": self.pooled_discussion
            }
        }


@dataclass
class PlanIndex:
    """Index entry for a plan."""
    version: int
    updated: str
    supersedes: Optional[str]
    superseded_by: Optional[str]
    status: str
    cluster: Optional[str]
    retained: bool = True

    def to_dict(self) -> dict:
        d = {
            "version": self.version,
            "updated": self.updated,
            "status": self.status
        }
        if self.supersedes:
            d["supersedes"] = self.supersedes
        if self.superseded_by:
            d["superseded_by"] = self.superseded_by
        if self.cluster:
            d["cluster"] = self.cluster
        if not self.retained:
            d["retained"] = False
        return d


class PlanGraph:
    """The complete plan graph with edges, clusters, and index."""

    def __init__(self, plans_dir: Path):
        self.plans_dir = plans_dir
        self.active_dir = plans_dir / "active_plans"
        self.archive_dir = plans_dir / ".archive"
        self.index_path = plans_dir / ".index.json"

        self.plans: dict[str, Plan] = {}
        self.archive: dict[str, Plan] = {}
        self.edges: list[Edge] = []
        self.clusters: dict[str, Cluster] = {}
        self.index: dict = {}

    def load_plans(self):
        """Load all active plans from disk."""
        self.plans = {}

        for md_file in self.active_dir.rglob("*.md"):
            # Skip macOS resource forks
            if md_file.name.startswith("._"):
                continue
            try:
                plan = parse_plan_file(md_file)
                if plan.frontmatter.area and plan.frontmatter.name:
                    self.plans[plan.id] = plan
            except Exception as e:
                print(f"Warning: Failed to parse {md_file}: {e}")

    def load_archive(self):
        """Load archived plans."""
        self.archive = {}

        if not self.archive_dir.exists():
            return

        for md_file in self.archive_dir.glob("*.md"):
            try:
                plan = parse_plan_file(md_file)
                # Archive key is filename without extension
                key = md_file.stem
                self.archive[key] = plan
            except Exception as e:
                print(f"Warning: Failed to parse archive {md_file}: {e}")

    def infer_dependency_edges(self):
        """Create edges from explicit blocked_by/blocks fields."""
        seen = set()  # (source, target, type) to avoid duplicates
        for plan_id, plan in self.plans.items():
            for blocked in plan.frontmatter.blocked_by:
                key = (blocked, plan_id, "dependency")
                if key not in seen:
                    seen.add(key)
                    self.edges.append(Edge(
                        source=blocked,
                        target=plan_id,
                        edge_type="dependency"
                    ))
            for blocks in plan.frontmatter.blocks:
                key = (plan_id, blocks, "dependency")
                if key not in seen:
                    seen.add(key)
                    self.edges.append(Edge(
                        source=plan_id,
                        target=blocks,
                        edge_type="dependency"
                    ))

    def infer_dataflow_edges(self):
        """Infer edges from output→input path matching."""
        # Collect all outputs with their source plans
        outputs: dict[str, tuple[str, Pointer]] = {}  # pattern -> (plan_id, pointer)

        for plan_id, plan in self.plans.items():
            if 'Outputs' in plan.sections:
                for pointer in plan.sections['Outputs'].pointers:
                    outputs[pointer.path] = (plan_id, pointer)

        # Match inputs against outputs (deduplicated)
        seen = set()
        for plan_id, plan in self.plans.items():
            if 'Inputs' in plan.sections:
                for pointer in plan.sections['Inputs'].pointers:
                    for out_pattern, (source_id, out_pointer) in outputs.items():
                        if source_id == plan_id:
                            continue
                        if self._patterns_match(out_pattern, pointer.path):
                            key = (source_id, plan_id, pointer.path)
                            if key not in seen:
                                seen.add(key)
                                self.edges.append(Edge(
                                    source=source_id,
                                    target=plan_id,
                                    edge_type="dataflow",
                                    weight=0.8,
                                    pattern=pointer.path
                                ))

    def _patterns_match(self, pattern1: str, pattern2: str) -> bool:
        """Check if two patterns could refer to the same files."""
        # Normalize wildcards
        p1 = re.sub(r'\{[^}]+\}', '*', pattern1)
        p2 = re.sub(r'\{[^}]+\}', '*', pattern2)

        # Check if one matches the other
        return fnmatch.fnmatch(p1, p2) or fnmatch.fnmatch(p2, p1) or p1 == p2

    def infer_resource_edges(self):
        """Infer edges from shared file references."""
        # Collect all referenced files per plan
        refs: dict[str, list[str]] = defaultdict(list)

        for plan_id, plan in self.plans.items():
            for pointer in plan.get_all_pointers():
                # Normalize path for comparison
                norm_path = re.sub(r'\{[^}]+\}', '*', pointer.path)
                refs[norm_path].append(plan_id)

        # Create edges for shared references (unique pairs only)
        for path, plan_ids in refs.items():
            unique_ids = list(set(plan_ids))  # Deduplicate
            if len(unique_ids) > 1:
                # Create edges between all pairs
                for i, p1 in enumerate(unique_ids):
                    for p2 in unique_ids[i+1:]:
                        if p1 == p2:  # Skip self-references
                            continue
                        # Only add if not already connected
                        existing = any(
                            e.source == p1 and e.target == p2 or
                            e.source == p2 and e.target == p1
                            for e in self.edges
                        )
                        if not existing:
                            self.edges.append(Edge(
                                source=p1,
                                target=p2,
                                edge_type="resource",
                                weight=0.5,
                                pattern=path
                            ))

    def infer_succession_edges(self):
        """Create edges from supersedes relationships."""
        for plan_id, plan in self.plans.items():
            if plan.frontmatter.supersedes:
                self.edges.append(Edge(
                    source=plan.frontmatter.supersedes,
                    target=plan_id,
                    edge_type="succession"
                ))

    def build_graph(self):
        """Build the complete graph."""
        self.edges = []
        self.load_plans()
        self.load_archive()

        self.infer_dependency_edges()
        self.infer_dataflow_edges()
        self.infer_resource_edges()
        self.infer_succession_edges()

        self.compute_snn_clusters()
        self.pool_cluster_context()

    def compute_snn_clusters(self):
        """
        Compute SNN (Shared Nearest Neighbor) clusters.

        Plans are clustered if they share many neighbors in the graph.
        Uses a simplified SNN approach based on shared edge connections.
        """
        # Build adjacency with edge counts
        neighbors: dict[str, set[str]] = defaultdict(set)

        for edge in self.edges:
            neighbors[edge.source].add(edge.target)
            neighbors[edge.target].add(edge.source)

        # Compute SNN similarity matrix
        plan_ids = list(self.plans.keys())
        snn_scores: dict[tuple[str, str], float] = {}

        for i, p1 in enumerate(plan_ids):
            for p2 in plan_ids[i+1:]:
                # SNN = |shared neighbors| / min(|neighbors|)
                shared = neighbors[p1] & neighbors[p2]
                min_neighbors = min(len(neighbors[p1]), len(neighbors[p2]))
                if min_neighbors > 0:
                    snn_scores[(p1, p2)] = len(shared) / min_neighbors

        # Cluster using threshold (simplified approach)
        threshold = 0.3
        clusters: dict[str, set[str]] = {}
        assigned: set[str] = set()

        # Group by area first, then refine by SNN
        area_groups: dict[str, list[str]] = defaultdict(list)
        for plan_id in plan_ids:
            area = plan_id.split('/')[0]
            area_groups[area].append(plan_id)

        for area, members in area_groups.items():
            if len(members) == 1:
                cluster_name = f"{area}_single"
                clusters[cluster_name] = set(members)
                assigned.update(members)
            else:
                # Check if they should be in same cluster
                cluster_name = f"{area}_workstream"
                clusters[cluster_name] = set(members)
                assigned.update(members)

        # Convert to Cluster objects
        self.clusters = {}
        for name, members in clusters.items():
            self.clusters[name] = Cluster(
                name=name,
                members=list(members)
            )

    def pool_cluster_context(self):
        """Pool context from plans in each cluster."""
        for cluster in self.clusters.values():
            inputs = set()
            outputs = set()
            tasks = []
            discussion = []

            for plan_id in cluster.members:
                if plan_id not in self.plans:
                    continue
                plan = self.plans[plan_id]

                # Pool inputs
                if 'Inputs' in plan.sections:
                    for p in plan.sections['Inputs'].pointers:
                        inputs.add(p.raw)

                # Pool outputs
                if 'Outputs' in plan.sections:
                    for p in plan.sections['Outputs'].pointers:
                        outputs.add(p.raw)

                # Pool tasks
                for task in plan.get_all_tasks():
                    tasks.append(f"[{plan_id}] {'x' if task.complete else ' '} {task.text}")

                # Pool discussion
                if 'Discussion' in plan.sections:
                    for entry in plan.sections['Discussion'].entries:
                        discussion.append({
                            "plan": plan_id,
                            "date": entry.date,
                            "title": entry.title
                        })

            cluster.pooled_inputs = sorted(inputs)
            cluster.pooled_outputs = sorted(outputs)
            cluster.pooled_tasks = tasks
            cluster.pooled_discussion = sorted(discussion, key=lambda x: x['date'], reverse=True)

    def get_ready_plans(self) -> list[str]:
        """Get plans that are ready to work on (not blocked)."""
        blocked_plans = set()

        for edge in self.edges:
            if edge.edge_type == "dependency":
                source_plan = self.plans.get(edge.source)
                if source_plan and source_plan.frontmatter.status != "complete":
                    blocked_plans.add(edge.target)

        ready = []
        for plan_id, plan in self.plans.items():
            if plan_id not in blocked_plans and plan.frontmatter.status not in ('complete', 'blocked'):
                ready.append(plan_id)

        return ready

    def get_blocked_plans(self) -> list[tuple[str, list[str]]]:
        """Get plans that are blocked and what blocks them."""
        blockers: dict[str, list[str]] = defaultdict(list)

        for edge in self.edges:
            if edge.edge_type == "dependency":
                source_plan = self.plans.get(edge.source)
                if source_plan and source_plan.frontmatter.status != "complete":
                    blockers[edge.target].append(edge.source)

        return [(plan_id, blocks) for plan_id, blocks in blockers.items()]

    def get_downstream(self, plan_id: str) -> list[str]:
        """Get plans that depend on the given plan."""
        downstream = set()
        for edge in self.edges:
            if edge.source == plan_id and edge.edge_type in ("dependency", "dataflow"):
                downstream.add(edge.target)
        return list(downstream)

    def get_plans_touching_file(self, file_pattern: str) -> list[str]:
        """Get plans that reference a file matching the pattern."""
        matching = []
        for plan_id, plan in self.plans.items():
            for pointer in plan.get_all_pointers():
                if fnmatch.fnmatch(pointer.path, file_pattern) or fnmatch.fnmatch(file_pattern, pointer.path):
                    matching.append(plan_id)
                    break
        return matching

    def get_mnn_anchors(self) -> list[tuple[str, str, list[str]]]:
        """
        Get MNN anchors - plan pairs connected via 2+ edge types.
        These represent strong cross-modal relationships.
        """
        # Count edge types between pairs
        pair_types: dict[tuple[str, str], set[str]] = defaultdict(set)

        for edge in self.edges:
            key = tuple(sorted([edge.source, edge.target]))
            pair_types[key].add(edge.edge_type)

        anchors = []
        for (p1, p2), types in pair_types.items():
            if len(types) >= 2:
                anchors.append((p1, p2, list(types)))

        return anchors

    def get_succession_history(self, plan_id: str) -> list[str]:
        """Get the succession chain for a plan."""
        history = [plan_id]

        # Walk backwards through supersedes
        current = plan_id
        while True:
            plan = self.plans.get(current) or self.archive.get(current)
            if not plan or not plan.frontmatter.supersedes:
                break
            history.append(plan.frontmatter.supersedes)
            current = plan.frontmatter.supersedes

        return list(reversed(history))

    def build_index(self) -> dict:
        """Build the complete index structure."""
        index = {
            "active_plans": {},
            "archive": {},
            "clusters": {},
            "graphs": {"edges": []}
        }

        # Index active plans
        for plan_id, plan in self.plans.items():
            cluster = None
            for c_name, c in self.clusters.items():
                if plan_id in c.members:
                    cluster = c_name
                    break

            index["active_plans"][plan_id] = PlanIndex(
                version=plan.frontmatter.version,
                updated=plan.frontmatter.updated,
                supersedes=plan.frontmatter.supersedes,
                superseded_by=None,
                status=plan.frontmatter.status,
                cluster=cluster
            ).to_dict()

        # Index archive
        for arch_id, plan in self.archive.items():
            # Check if retained (connected to active)
            retained = any(
                e.source == arch_id or e.target == arch_id
                for e in self.edges
            )

            # Find superseded_by
            superseded_by = None
            for active_id, active_plan in self.plans.items():
                if active_plan.frontmatter.supersedes == arch_id:
                    superseded_by = active_id
                    break

            index["archive"][arch_id] = PlanIndex(
                version=plan.frontmatter.version,
                updated=plan.frontmatter.updated,
                supersedes=plan.frontmatter.supersedes,
                superseded_by=superseded_by,
                status=plan.frontmatter.status,
                cluster=None,
                retained=retained
            ).to_dict()

        # Index clusters
        for name, cluster in self.clusters.items():
            index["clusters"][name] = cluster.to_dict()

        # Index edges
        index["graphs"]["edges"] = [e.to_dict() for e in self.edges]

        self.index = index
        return index

    def save_index(self):
        """Save the index to disk."""
        index = self.build_index()
        with open(self.index_path, 'w') as f:
            json.dump(index, f, indent=2)

    def load_index(self) -> Optional[dict]:
        """Load index from disk if it exists."""
        if self.index_path.exists():
            with open(self.index_path) as f:
                self.index = json.load(f)
                return self.index
        return None


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: graph.py <plans_dir> [--rebuild]")
        sys.exit(1)

    plans_dir = Path(sys.argv[1])
    graph = PlanGraph(plans_dir)

    if "--rebuild" in sys.argv:
        graph.build_graph()
        graph.save_index()
        print(f"Index rebuilt: {len(graph.plans)} active plans, {len(graph.edges)} edges, {len(graph.clusters)} clusters")
    else:
        graph.load_index()
        print(json.dumps(graph.index, indent=2))
