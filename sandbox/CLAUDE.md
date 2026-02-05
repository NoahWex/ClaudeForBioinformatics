# Claude Code Sandbox - Onboarding Guide

This file instructs Claude on how to guide users through setup. Users can also read it directly or ask Claude about specific sections.

## What This Is

A hands-on learning environment for using Claude Code with HPC (UCI HPC3). Includes:
- Configurable setup for different environments
- Step-by-step guides
- The hpc-toolkit for safe HPC operations
- Guided exercises

---

## How to Guide Setup

When a user says they want to set this up, follow this process:

### Step 1: Assess Current State

Ask or check:

```
1. Do you have a UCI HPC3 account? (or other SLURM cluster)
2. Have you set up SSH keys for HPC access?
3. Is Claude Code already installed?
4. Have you used the hpc-toolkit before?
```

Based on answers, identify which guides they need.

### Step 2: Determine Profile

Ask:
```
Which cluster will you use?
- UCI HPC3 (default) → profile: uci-hpc3
- Other SLURM cluster → profile: generic-slurm
```

### Step 3: Guide Through Prerequisites

If missing prerequisites, walk through:
- **No HPC account** → Point to `setup-guides/00_prerequisites.md`
- **No SSH setup** → Guide through `setup-guides/01_ssh_config.md`
- **No Claude Code** → Guide through `setup-guides/02_claude_code_install.md`

### Step 4: Configure sandbox.yaml

Help them edit `config/sandbox.yaml`:

```yaml
profile: uci-hpc3  # or generic-slurm

cluster:
  ssh_alias: hpc3           # Their SSH config alias
  account: their_lab        # Their SLURM account
  user: their_ucinetid      # Their username
```

### Step 5: Run Setup

```bash
cd /path/to/sandbox
./setup.sh
```

Then verify:
```bash
./setup.sh --verify
```

### Step 6: Test Connection

Have them run:
```bash
ssh hpc3 hostname
```

If DUO prompts, that's expected on first connection.

### Step 7: Start Exercises

Direct them to `exercises/01_first_hpc_command.md` to learn interactively.

---

## Documentation Index

Point users here when they want to learn more or ask about specific topics.

### Core Concepts

| Topic | File | Description |
|-------|------|-------------|
| Project setup | `../01_CLAUDE_MD_STRUCTURE.md` | How to write CLAUDE.md files |
| Settings & hooks | `../02_SETTINGS_AND_HOOKS.md` | Configure Claude Code behavior |
| Rules framework | `../03_RULES_FRAMEWORK.md` | Path-triggered contextual rules |
| Anti-sycophancy | `../04_TRACE_FRAMEWORK.md` | Objectivity and scope discipline |
| Manifest patterns | `../05_MANIFEST_PATTERNS.md` | YAML-based configuration |
| Experiment lifecycle | `../06_EXPERIMENT_LIFECYCLE.md` | Managing analysis experiments |

### HPC & Infrastructure

| Topic | File | Description |
|-------|------|-------------|
| HPC patterns | `../07_HPC_INFRASTRUCTURE.md` | SLURM scripts, containers, gotchas |
| Skills & agents | `../08_SKILLS_AND_AGENTS.md` | Custom Claude extensions |
| Deprecation | `../09_DEPRECATION_POLICY.md` | Retiring code properly |
| **HPC toolkit** | `../10_HPC_PLUGIN.md` | hpc-toolkit CLI, permissions |
| Interface planning | `../11_INTERFACE_AWARE_PLANNING.md` | Multi-component planning |
| Sandbox overview | `../12_SANDBOX_ENVIRONMENT.md` | This environment explained |

### Setup Guides (Step-by-Step)

| Guide | File | When to Use |
|-------|------|-------------|
| Prerequisites | `setup-guides/00_prerequisites.md` | Before anything else |
| SSH config | `setup-guides/01_ssh_config.md` | Set up HPC connection |
| Claude install | `setup-guides/02_claude_code_install.md` | Install Claude Code CLI |
| Toolkit install | `setup-guides/03_plugin_install.md` | Install hpc-toolkit |
| Verify setup | `setup-guides/04_verify_setup.md` | Confirm everything works |
| Troubleshooting | `setup-guides/troubleshooting.md` | Fix common issues |

### Exercises (Hands-On Learning)

| Exercise | File | Skills Learned |
|----------|------|----------------|
| 1. First HPC command | `exercises/01_first_hpc_command.md` | Read-only commands, hook behavior |
| 2. Submit test job | `exercises/02_submit_test_job.md` | hpc submit, --dry-run |
| 3. Git via SSH | `exercises/03_git_via_ssh.md` | Git permission tiers |
| 4. Create experiment | `exercises/04_create_experiment.md` | Project structure, manifests |

### Configuration Files

| File | Purpose |
|------|---------|
| `config/sandbox.yaml` | User's configuration (edit this) |
| `config/profiles/uci-hpc3.yaml` | UCI HPC3 defaults |
| `config/profiles/generic-slurm.yaml` | Template for other clusters |
| `global-setup/settings.json.template` | Global Claude settings |
| `global-setup/rules/*.md` | Global rules (TRACE, safety) |

### Templates

| Template | Purpose |
|----------|---------|
| `project-template/CLAUDE.md.template` | Starting point for new projects |
| `project-template/.claude/` | Project-level Claude config |
| `project-template/project/config/` | Manifest templates |

---

## Common User Questions

### "What is the hpc-toolkit?"

A unified CLI (`~/.claude/hpc-toolkit/bin/hpc`) that replaces raw SSH commands. Subcommands: `submit`, `shell`, `status`, `logs`, `cancel`, `file`, `git`, `workspace`. Every operation requires purpose documentation and is audit-logged.

See `../10_HPC_PLUGIN.md` for full details.

### "Why is sbatch blocked?"

To enforce documentation. Every job submission should use `hpc submit` which requires purpose and expected outputs. The toolkit creates an audit trail at `~/.claude/hpc-toolkit/logs/audit.jsonl`.

### "What's the TRACE framework?"

An anti-sycophancy and objectivity framework. It prevents Claude from giving unwarranted praise and enforces scope discipline. See `../04_TRACE_FRAMEWORK.md`.

### "How do I customize for my cluster?"

1. Copy `config/profiles/generic-slurm.yaml`
2. Edit with your cluster details
3. Set `profile: your-profile` in `sandbox.yaml`
4. Update paths in the profile

### "What if I don't use UCI HPC?"

The toolkit works with any SSH-accessible SLURM cluster. You need to:
1. Create an SSH alias in `~/.ssh/config`
2. Update `config/sandbox.yaml` with your alias
3. Edit `~/.claude/hpc-toolkit/config.sh` with your SSH_ALIAS and paths

---

## Troubleshooting Guidance

### SSH Issues

If `ssh hpc3 hostname` fails:
1. Check VPN if off-campus
2. Verify SSH config: `cat ~/.ssh/config`
3. Debug: `ssh -v hpc3 hostname`
4. See `setup-guides/troubleshooting.md`

### Toolkit Not Working

If commands aren't being intercepted:
1. Check settings.json has PreToolUse hook for gate.sh
2. Verify toolkit path: `ls ~/.claude/hpc-toolkit/bin/hpc`
3. Make executable: `chmod +x ~/.claude/hpc-toolkit/bin/hpc ~/.claude/hpc-toolkit/hooks/gate.sh`
4. Restart Claude Code

### Hook Errors

If hooks produce errors:
1. Test directly: `echo '{"tool_input":{"command":"ssh hpc3 squeue"}}' | ~/.claude/hpc-toolkit/hooks/gate.sh`
2. Check Python 3 available
3. Check hook script syntax

---

## Quick Commands Reference

```bash
# Run setup
./setup.sh

# Verify setup (non-destructive)
./setup.sh --verify

# Test SSH
ssh hpc3 hostname

# Test hpc-toolkit
~/.claude/hpc-toolkit/bin/hpc --help

# Check job queue
ssh hpc3 "squeue -u \$USER"
```

---

## For Claude: Onboarding Checklist

When guiding a new user, ensure they complete:

- [ ] Prerequisites confirmed (account, SSH key, CRSP access)
- [ ] SSH config created with socket multiplexing
- [ ] Claude Code installed and authenticated
- [ ] `sandbox.yaml` configured for their environment
- [ ] `./setup.sh` run successfully
- [ ] `./setup.sh --verify` passes all checks
- [ ] First SSH connection works (DUO completed)
- [ ] Exercise 1 completed (read-only commands work)
- [ ] Exercise 2 completed (understand hpc submit)

After setup, suggest they:
1. Complete remaining exercises
2. Read documentation topics relevant to their work
3. Copy `project-template/` to start their first project
