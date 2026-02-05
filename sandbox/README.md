# Claude Code Sandbox

A hands-on learning environment for Claude Code with HPC integration.

## Two Ways to Use This

### Option 1: Claude-Guided Setup (Recommended)

Start Claude Code in this directory and say:

> "I want to set up Claude Code for my HPC workflow. Guide me through it."

Claude will read `CLAUDE.md` and walk you through:
- Assessing what you already have
- Configuring for your environment
- Running setup and verification
- Completing the exercises

### Option 2: Self-Directed Setup

Follow the guides yourself:

```bash
# 1. Read prerequisites
cat setup-guides/00_prerequisites.md

# 2. Edit your configuration
vi config/sandbox.yaml

# 3. Run setup
./setup.sh

# 4. Verify
./setup.sh --verify

# 5. Start exercises
cat exercises/01_first_hpc_command.md
```

You can also ask Claude about specific topics:
> "Explain the hpc-toolkit"
> "What's in 07_HPC_INFRASTRUCTURE.md?"
> "Help me troubleshoot SSH connection issues"

### If Working From Another Directory

If you're not in the sandbox directory, tell Claude where it is:

> "I want to set up Claude Code for HPC. The sandbox is at /share/crsp/lab/dalawson/share/0_Resources/ClaudeForBioinformatics_NW/sandbox/"

Claude will read the CLAUDE.md from there and guide you.

## What's Here

```
sandbox/
├── config/                 # Your settings
│   ├── sandbox.yaml        # EDIT THIS - your cluster config
│   └── profiles/           # Profile templates
├── setup-guides/           # Step-by-step installation
├── global-setup/           # Files for ~/.claude/
├── project-template/       # Starting point for new projects
└── exercises/              # Guided learning
```

## Before You Start

1. **UCI HPC3 account** - Request at https://rcic.uci.edu/account.html
2. **SSH key** - Generate with `ssh-keygen -t ed25519`
3. **CRSP access** - Storage allocation for your lab
4. **Claude Code** - Install from Anthropic

## Configuration

Edit `config/sandbox.yaml`:

```yaml
profile: uci-hpc3
cluster:
  ssh_alias: hpc3
  account: your_lab_account
  user: your_ucinetid
```

## Exercises

| # | Exercise | Skills |
|---|----------|--------|
| 1 | First HPC Command | Read-only SLURM, hook behavior |
| 2 | Submit Test Job | hpc submit, --dry-run |
| 3 | Git via SSH | Permission tiers, CRSP repos |
| 4 | Create Experiment | Project patterns, lifecycle |

## Help

- `./setup.sh --help` - Setup options
- `exercises/*.md` - Detailed instructions
- `setup-guides/troubleshooting.md` - Common issues

## External Docs

- UCI HPC: https://rcic.uci.edu/
- Claude Code: https://docs.anthropic.com/en/docs/claude-code
