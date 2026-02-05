# Exercise 4: Claude Customization

**Objective**: Learn to customize Claude Code's behavior with rules and understand the configuration system.

**Time**: ~35 minutes

---

## Background

Claude Code can be customized through:
1. **Rules files** (`.claude/rules/`) - Instructions Claude follows
2. **Settings** (`.claude/settings.json`) - Configuration options
3. **Hooks** - Scripts that run before/after Claude actions
4. **Skills** - Custom commands (advanced)

In this exercise, you'll explore existing customizations and create your own rule.

## Part 1: Explore Existing Rules (10 min)

### Review the HPC safety rules

Ask Claude to show you the existing rules:

```
Please show me the Claude rules that are currently active for this project.
I'm particularly interested in:
1. The anti-sycophancy framework
2. The HPC safety guidelines

Explain what each rule does and why it exists.
```

### Understand rule structure

```
Please explain:
1. Where are rules files stored?
2. How does Claude know to follow them?
3. What's the difference between project rules and global rules?
4. Can rules conflict? How are conflicts resolved?
```

## Part 2: Create Your Own Rule (15 min)

### Design a lab-specific rule

Think about conventions in your lab. Common examples:
- File naming conventions
- Required metadata fields
- Preferred coding style
- Safety checks for specific data types

Ask Claude to help you create a rule:

```
Please help me create a Claude rule for my lab. I want to enforce:

1. File naming: Use YYYY-MM-DD prefix for date-stamped outputs
2. Required fields: All Seurat objects should have 'project' and 'analyst' in metadata
3. Documentation: All R scripts should have a header comment block

Create this as working/.claude/rules/lab_conventions.md
```

### Test the rule

After creating the rule, test it:

```
Now that the lab conventions rule is active, please create a new R script
called working/scripts/test_script.R that loads a Seurat object.

Does Claude follow the header comment convention from the rule?
```

## Part 3: Understand the Configuration System (10 min)

### Settings overview

```
Please explain Claude Code's configuration system:
1. What can be configured in settings.json?
2. What's the difference between .claude/settings.json and .claude/settings.local.json?
3. How do I see what settings are currently active?
```

### Hooks (demonstration)

```
Please explain hooks in Claude Code:
1. What are PreToolUse and PostToolUse hooks?
2. Show me an example of the git protection hook used in this project
3. When would I create my own hook?

Note: Don't create any hooks, just explain them.
```

---

## Rules File Format

Rules are Markdown files that Claude reads and follows:

```markdown
# My Lab Conventions

These rules apply to all work in this project.

## File Naming

- Use ISO date format: YYYY-MM-DD
- Include analyst initials in output files
- Example: 2024-01-15_nw_integration_results.rds

## Code Standards

- All scripts must have header comments including:
  - Purpose
  - Author
  - Date
  - Dependencies

## Data Safety

- Never overwrite raw data files
- Always save processed data to outputs/ directory
```

## Rule File Locations

| Location | Scope |
|----------|-------|
| `~/.claude/rules/` | All projects (global) |
| `.claude/rules/` | This project only |
| `CLAUDE.md` | Project root (legacy) |

## Example Rules from This Workshop

### Anti-sycophancy (trace-framework.md)
- Prevents excessive praise
- Requires factual accuracy
- Uses triggers like `CLARIFY:FIRST`

### HPC Safety (hpc-safety-guidelines.md)
- Verification before job submission
- Severity levels for destructive commands
- Git protection for CRSP-mounted repos

## Best Practices for Rules

1. **Be specific** - Vague rules are inconsistently followed
2. **Include examples** - Show what you want, not just what you don't want
3. **Organize by topic** - One rule file per domain
4. **Test your rules** - Verify Claude follows them
5. **Document why** - Future you will thank you

## Discussion Questions

1. What conventions in your lab would benefit from a rule?
2. How do rules differ from just asking Claude to do something?
3. When might a rule be counterproductive?

## Checkpoint Verification

Before finishing, verify:
- [ ] Understand where rules are stored
- [ ] Created a custom rule file
- [ ] Tested that Claude follows the rule
- [ ] Understand the difference between rules, settings, and hooks

## Workshop Wrap-up

Congratulations! You've completed all four sessions. You now know how to:

1. **Day 1**: Use Claude to scaffold projects and load data
2. **Day 2**: Explore and harmonize metadata across datasets
3. **Day 3**: Write and submit HPC jobs with Claude assistance
4. **Day 4**: Integrate data, create visualizations, generate reports, and customize Claude

## Next Steps for Your Own Work

1. **Start small** - Use Claude for one part of your next analysis
2. **Iterate** - Refine your prompts and rules over time
3. **Share** - Tell colleagues about useful patterns you discover
4. **Contribute** - Suggest improvements to lab rules and templates

## Resources

- [Claude Code Documentation](https://docs.anthropic.com/claude-code)
- Workshop materials: Keep this directory as reference
- Lab shared resources: `/share/crsp/lab/dalawson/share/0_Resources/`
