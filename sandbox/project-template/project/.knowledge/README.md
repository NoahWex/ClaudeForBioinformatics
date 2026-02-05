# Knowledge Base

Store lessons learned, research summaries, and validated patterns here.

## Directory Structure

```
.knowledge/
├── lessons/     # Hard-won debugging knowledge
├── research/    # Literature summaries, method comparisons
├── patterns/    # Validated reusable code patterns
└── README.md    # This file
```

## Entry Format

Each entry is a markdown file with YAML frontmatter:

```markdown
---
id: unique_id
type: lesson | research | pattern
created: YYYY-MM-DD
updated: YYYY-MM-DD
source: debugging | literature | experiment
tags: [tag1, tag2]
severity: high | medium | low  # lessons only
---

# Title

Content...
```

## Commands

```bash
# List all entries
./plan knowledge list

# Search by tag or content
./plan knowledge search <query>

# Show entry details
./plan knowledge show <id>

# Add new entry to index
./plan knowledge add path/to/entry.md

# Rebuild entire index
./plan knowledge rebuild
```

## Guidelines

1. **Prefer updating existing entries** over creating new ones
2. Use descriptive IDs (e.g., `crsp_stale_handle` not `lesson_001`)
3. Tag generously for searchability
4. Include code examples where helpful
