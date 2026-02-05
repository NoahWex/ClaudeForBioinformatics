#!/usr/bin/env python3
"""
Plan parser - Extract sections and pointers from plan documents.

Handles:
- YAML frontmatter extraction
- Section parsing with content type detection
- Pointer syntax with wildcards ({name}, *, **)
- Checklist parsing for task tracking
"""

import re
import yaml
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional
from datetime import date


# Haha, the parser walks into a bar and says "I'll have a regex"
# The bartender says "We don't serve your type here"


@dataclass
class Pointer:
    """A file/config reference with optional wildcards."""
    raw: str
    path: str
    description: str = ""
    wildcards: list[str] = field(default_factory=list)
    line_ref: Optional[int] = None

    @classmethod
    def parse(cls, line: str) -> Optional["Pointer"]:
        """Parse a pointer line like `path/to/file.ext` — description"""
        # Match: - `path` — description OR - `path:line` — description
        match = re.match(r'^-\s*`([^`]+)`(?:\s*[—-]\s*(.*))?$', line.strip())
        if not match:
            return None

        raw_path = match.group(1)
        description = match.group(2) or ""

        # Extract line reference if present
        line_ref = None
        path = raw_path
        if ':' in raw_path:
            parts = raw_path.rsplit(':', 1)
            if parts[1].isdigit():
                path = parts[0]
                line_ref = int(parts[1])

        # Find wildcards
        wildcards = re.findall(r'\{(\w+)\}', path)
        wildcards.extend(['*'] * path.count('*'))

        return cls(
            raw=raw_path,
            path=path,
            description=description,
            wildcards=wildcards,
            line_ref=line_ref
        )


@dataclass
class TaskItem:
    """A checklist item with completion status."""
    text: str
    complete: bool

    @classmethod
    def parse(cls, line: str) -> Optional["TaskItem"]:
        """Parse a task line like - [x] Do something"""
        match = re.match(r'^-\s*\[([ xX])\]\s*(.+)$', line.strip())
        if not match:
            return None
        return cls(
            text=match.group(2),
            complete=match.group(1).lower() == 'x'
        )


@dataclass
class DiscussionEntry:
    """A dated discussion entry."""
    date: str
    title: str
    content: list[str]

    @classmethod
    def parse_all(cls, lines: list[str]) -> list["DiscussionEntry"]:
        """Parse all discussion entries from section content."""
        entries = []
        current = None

        for line in lines:
            # Match: ### 2026-02-04: Title
            match = re.match(r'^###\s*(\d{4}-\d{2}-\d{2}):\s*(.+)$', line)
            if match:
                if current:
                    entries.append(current)
                current = cls(
                    date=match.group(1),
                    title=match.group(2),
                    content=[]
                )
            elif current and line.strip():
                current.content.append(line)

        if current:
            entries.append(current)

        return entries


@dataclass
class Section:
    """A plan section with parsed content."""
    name: str
    content_type: str  # pointer_list, checklist, dated_entries, freeform
    raw_lines: list[str]
    pointers: list[Pointer] = field(default_factory=list)
    tasks: list[TaskItem] = field(default_factory=list)
    entries: list[DiscussionEntry] = field(default_factory=list)

    # Poolable sections can be merged across plans in a cluster
    POOLABLE = {'Inputs', 'Outputs', 'Config', 'Scripts', 'Tasks', 'Discussion'}

    @property
    def is_poolable(self) -> bool:
        return self.name in self.POOLABLE


@dataclass
class Frontmatter:
    """Plan frontmatter (YAML header)."""
    area: str
    name: str
    version: int
    updated: str
    supersedes: Optional[str] = None
    blocked_by: list[str] = field(default_factory=list)
    blocks: list[str] = field(default_factory=list)
    status: str = "queued"
    jobs: list[str] = field(default_factory=list)

    @classmethod
    def parse(cls, yaml_str: str) -> "Frontmatter":
        """Parse frontmatter from YAML string."""
        data = yaml.safe_load(yaml_str) or {}
        # Convert date objects to strings (YAML parses dates automatically)
        updated = data.get('updated', str(date.today()))
        if hasattr(updated, 'isoformat'):
            updated = updated.isoformat()
        elif not isinstance(updated, str):
            updated = str(updated)
        return cls(
            area=data.get('area', ''),
            name=data.get('name', ''),
            version=data.get('version', 1),
            updated=updated,
            supersedes=data.get('supersedes'),
            blocked_by=data.get('blocked_by', []),
            blocks=data.get('blocks', []),
            status=data.get('status', 'queued'),
            jobs=[str(j) for j in data.get('jobs', [])]
        )

    def validate(self) -> list[str]:
        """Return list of validation errors."""
        errors = []
        if not self.area:
            errors.append("Missing required field: area")
        if not self.name:
            errors.append("Missing required field: name")
        if not isinstance(self.version, int) or self.version < 1:
            errors.append("version must be a positive integer")
        if not self.updated:
            errors.append("Missing required field: updated")
        if self.status not in ('queued', 'in_progress', 'blocked', 'complete'):
            errors.append(f"Invalid status: {self.status}")
        return errors


@dataclass
class Plan:
    """A complete parsed plan document."""
    path: Path
    frontmatter: Frontmatter
    title: str
    sections: dict[str, Section]
    raw_content: str

    REQUIRED_SECTIONS = [
        'Inputs', 'Outputs', 'Config', 'Scripts',
        'Tasks', 'Current State', 'Discussion', 'Next Actions',
        'Completion'
    ]

    @property
    def id(self) -> str:
        """Unique plan identifier: area/name"""
        return f"{self.frontmatter.area}/{self.frontmatter.name}"

    def validate(self) -> list[str]:
        """Return list of validation errors."""
        errors = self.frontmatter.validate()

        for section in self.REQUIRED_SECTIONS:
            if section not in self.sections:
                errors.append(f"Missing required section: ## {section}")

        return errors

    def get_all_pointers(self) -> list[Pointer]:
        """Get all pointers across all sections."""
        pointers = []
        for section in self.sections.values():
            pointers.extend(section.pointers)
        return pointers

    def get_all_tasks(self) -> list[TaskItem]:
        """Get all tasks from Tasks section."""
        if 'Tasks' in self.sections:
            return self.sections['Tasks'].tasks
        return []

    def task_progress(self) -> tuple[int, int]:
        """Return (completed, total) task counts."""
        tasks = self.get_all_tasks()
        completed = sum(1 for t in tasks if t.complete)
        return completed, len(tasks)


def parse_plan(content: str, path: Optional[Path] = None) -> Plan:
    """Parse a plan document from string content."""

    # Extract frontmatter
    frontmatter_match = re.match(r'^---\n(.*?)\n---\n', content, re.DOTALL)
    if frontmatter_match:
        frontmatter = Frontmatter.parse(frontmatter_match.group(1))
        body = content[frontmatter_match.end():]
    else:
        frontmatter = Frontmatter(area='', name='', version=1, updated='')
        body = content

    # Extract title (first # heading)
    title_match = re.match(r'^#\s+(.+)$', body.strip(), re.MULTILINE)
    title = title_match.group(1) if title_match else ""

    # Parse sections
    sections = {}
    section_pattern = re.compile(r'^##\s+(.+)$', re.MULTILINE)
    section_starts = [(m.group(1), m.end()) for m in section_pattern.finditer(body)]

    for i, (name, start) in enumerate(section_starts):
        # Find end of section (next section or end of document)
        end = section_starts[i + 1][1] - len(f"## {section_starts[i + 1][0]}") - 1 if i + 1 < len(section_starts) else len(body)
        section_content = body[start:end].strip()
        lines = section_content.split('\n')

        # Detect content type and parse accordingly
        section = Section(name=name, content_type='freeform', raw_lines=lines)

        # Try parsing as pointer list
        pointers = [Pointer.parse(line) for line in lines]
        pointers = [p for p in pointers if p is not None]
        if pointers:
            section.content_type = 'pointer_list'
            section.pointers = pointers

        # Try parsing as checklist
        tasks = [TaskItem.parse(line) for line in lines]
        tasks = [t for t in tasks if t is not None]
        if tasks:
            section.content_type = 'checklist'
            section.tasks = tasks

        # Try parsing as dated entries
        if name == 'Discussion':
            entries = DiscussionEntry.parse_all(lines)
            if entries:
                section.content_type = 'dated_entries'
                section.entries = entries

        sections[name] = section

    return Plan(
        path=path or Path(""),
        frontmatter=frontmatter,
        title=title,
        sections=sections,
        raw_content=content
    )


def parse_plan_file(path: Path) -> Plan:
    """Parse a plan from a file path."""
    content = path.read_text()
    return parse_plan(content, path)


def expand_wildcard_pattern(pattern: str, values: dict[str, list[str]]) -> list[str]:
    """
    Expand a pattern with wildcards into concrete paths.

    Args:
        pattern: Path pattern like "outputs/{study}/data.rds"
        values: Dict mapping wildcard names to possible values
                e.g., {"study": ["gray", "kumar", "reed"]}

    Returns:
        List of expanded paths
    """
    # Find wildcards in pattern
    wildcards = re.findall(r'\{(\w+)\}', pattern)

    if not wildcards:
        return [pattern]

    # Build all combinations
    from itertools import product

    wildcard_values = []
    for wc in wildcards:
        if wc in values:
            wildcard_values.append(values[wc])
        else:
            wildcard_values.append([f'{{{wc}}}'])  # Keep unexpanded

    expanded = []
    for combo in product(*wildcard_values):
        result = pattern
        for wc, val in zip(wildcards, combo):
            result = result.replace(f'{{{wc}}}', val)
        expanded.append(result)

    return expanded


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: parse.py <plan_file> [--validate]")
        sys.exit(1)

    path = Path(sys.argv[1])
    plan = parse_plan_file(path)

    if "--validate" in sys.argv:
        errors = plan.validate()
        if errors:
            print("Validation errors:")
            for e in errors:
                print(f"  - {e}")
            sys.exit(1)
        else:
            print("Plan is valid")
            sys.exit(0)

    # Default: print summary
    print(f"Plan: {plan.id}")
    print(f"Title: {plan.title}")
    print(f"Status: {plan.frontmatter.status}")
    completed, total = plan.task_progress()
    print(f"Tasks: {completed}/{total} complete")
    print(f"Pointers: {len(plan.get_all_pointers())}")
    print(f"Sections: {list(plan.sections.keys())}")
