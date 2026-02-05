# Exercise 4: Create Experiment

Learn to initialize a project using Claude Code patterns.

## Objective

Set up a new project structure using:
1. CLAUDE.md for project instructions
2. Manifest files for configuration
3. Proper directory organization
4. Local Claude settings

## Prerequisites

- Completed Exercises 1-3
- Copied project-template from sandbox

## Steps

### Step 1: Copy Project Template

```bash
# Create new project from template
cp -r /path/to/sandbox/project-template ~/my-new-project
cd ~/my-new-project

# Rename template files
mv CLAUDE.md.template CLAUDE.md
mv .claude/settings.local.json.template .claude/settings.local.json
```

### Step 2: Customize CLAUDE.md

Edit CLAUDE.md with your project details:

```bash
# Open in editor
vi CLAUDE.md
```

Update:
- Project name and description
- HPC account name
- Key paths for your project

### Step 3: Start Claude Code

```bash
cd ~/my-new-project
claude
```

Claude now has your project context from CLAUDE.md.

### Step 4: Test Project Understanding

Ask Claude:

> "What is this project about?"

**Expected:** Claude describes your project based on CLAUDE.md.

> "Where should I put input data?"

**Expected:** Claude refers to data_manifest.yaml.

> "How do I submit HPC jobs?"

**Expected:** Claude explains hpc submit with your project's paths.

### Step 5: Create a Script

Ask Claude:

> "Create a simple Python script that prints 'Hello from project' and save it to project/scripts/"

**Expected:**
1. Claude creates `project/scripts/hello.py`
2. Follows any code style rules from `.claude/rules/`

### Step 6: Register in Manifest

Update the outputs manifest manually or ask:

> "Add the hello.py script to the outputs manifest"

Or edit directly:

```yaml
# project/config/outputs_manifest.yaml
outputs:
  hello_script:
    path: project/scripts/hello.py
    description: "Hello world test script"
    generated_by: Claude Code
    date: "2026-01-22"
```

### Step 7: Create SLURM Script

Ask Claude:

> "Create a SLURM script to run hello.py on HPC, put it in project/run/"

**Expected:** Claude creates `project/run/run_hello.sh` with:
- Proper SBATCH directives
- Your account from CLAUDE.md
- Container execution
- Error handling

### Step 8: Dry-Run Submission

```bash
~/.claude/hpc-toolkit/bin/hpc submit /path/to/project/run/run_hello.sh \
  --purpose "Test hello script on HPC" \
  --outputs "stdout" \
  --dry-run
```

**Expected:** Review box shows your script details.

## Verification Checklist

- [ ] Project template copied and customized
- [ ] CLAUDE.md edited with project details
- [ ] Claude understands project context
- [ ] Script created in correct location
- [ ] Manifest updated
- [ ] SLURM script follows project conventions
- [ ] Dry-run submission works

## Final Project Structure

```
my-new-project/
├── CLAUDE.md                    # Project instructions
├── .claude/
│   ├── settings.local.json      # Local settings
│   └── rules/
│       └── project-rules.md     # Project-specific rules
└── project/
    ├── config/
    │   ├── data_manifest.yaml   # Input data registry
    │   └── outputs_manifest.yaml # Output registry
    ├── scripts/
    │   └── hello.py             # Your script
    ├── run/
    │   └── run_hello.sh         # SLURM job script
    └── outputs/                 # Generated outputs
```

## What You Learned

1. **CLAUDE.md provides context** - Claude understands your project
2. **Manifests track data** - Single source of truth
3. **Templates speed setup** - Consistent project structure
4. **Local rules customize behavior** - Project-specific guidance
5. **Integration works** - HPC submission follows patterns

## Next Steps

From here you can:
1. Add real data to data_manifest.yaml
2. Create analysis scripts
3. Build HPC pipelines
4. Track outputs systematically

## Resources

- Main docs: `../` (ClaudeForBioinformatics_NW)
- HPC patterns: `../07_HPC_INFRASTRUCTURE.md`
- Plugin docs: `../10_HPC_PLUGIN.md`
- Interface planning: `../11_INTERFACE_AWARE_PLANNING.md`

## Congratulations!

You've completed the sandbox exercises. You now understand:
- How read-only commands auto-approve
- How hpc submit enforces documentation
- How git operations are tiered
- How to structure projects for Claude Code

Apply these patterns to your real projects!
