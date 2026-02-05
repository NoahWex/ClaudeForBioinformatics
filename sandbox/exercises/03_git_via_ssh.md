# Exercise 3: Git via SSH

Learn how git operations are tiered by the permission system.

## Objective

Understand the three permission tiers for git:
1. **Safe** - Auto-approved read-only and staging
2. **Supervised** - Prompts for commits and branch ops
3. **Dangerous** - Strong warnings for force operations

## Prerequisites

- Completed Exercises 1-2
- A git repository on CRSP (or create one)

## Background

Git repos on CRSP are accessed via SSH:
```bash
ssh hpc3 "cd /path/to/repo && git status"
```

The hook system applies different permission levels based on the git command.

## Steps

### Step 1: Create Test Repository (if needed)

```bash
ssh hpc3 "mkdir -p /share/crsp/lab/YOUR_LAB/YOUR_USER/git-exercise && \
  cd /share/crsp/lab/YOUR_LAB/YOUR_USER/git-exercise && \
  git init && \
  echo 'test' > README.md && \
  git add README.md && \
  git commit -m 'Initial commit'"
```

Set repo path for exercises:
```bash
REPO="/share/crsp/lab/YOUR_LAB/YOUR_USER/git-exercise"
```

### Step 2: Safe Operations (Auto-Approved)

In Claude, ask:

> "Check git status of my repository at [REPO_PATH]"

**Expected behavior:**
1. Claude runs `ssh hpc3 "cd $REPO && git status"`
2. Status displayed immediately
3. No permission prompt

Try more safe commands:

> "Show git log for [REPO_PATH]"
> "Show recent changes with git diff"

**Safe commands (auto-approved):**
- `git status`
- `git log`
- `git diff`
- `git show`
- `git add`
- `git stash`
- `git tag`
- `git fetch`

### Step 3: Supervised Operations (Prompts)

Create a test file and try to commit:

```bash
ssh hpc3 "cd $REPO && echo 'exercise 3' > exercise3.txt && git add exercise3.txt"
```

Now in Claude, ask:

> "Commit the changes in [REPO_PATH] with message 'Exercise 3 test'"

**Expected behavior:**
1. Hook detects `git commit`
2. Permission prompt appears
3. Message: "Git operation requires confirmation"
4. You must approve to proceed

**Supervised commands (prompt):**
- `git commit`
- `git push`
- `git pull`
- `git checkout`
- `git switch`
- `git branch`
- `git merge`
- `git rebase`
- `git reset`

### Step 4: Dangerous Operations (Strong Warning)

**DO NOT actually run these** - just observe the warning:

Ask Claude:

> "Force push to origin in [REPO_PATH]"

**Expected behavior:**
1. Hook detects `git push --force`
2. Strong warning prompt appears
3. Message explains this can cause permanent data loss

**Dangerous commands (strong warning):**
- `git push --force`
- `git reset --hard`
- `git clean -f`
- `git branch -D`

### Step 5: Understand the Tiering

Review what happened:

| Operation | Behavior | Why |
|-----------|----------|-----|
| `git status` | Auto-approved | Read-only, no risk |
| `git add` | Auto-approved | Staging is reversible |
| `git commit` | Prompt | Creates permanent history |
| `git push` | Prompt | Affects remote |
| `git push --force` | Strong warning | Can lose others' work |

### Step 6: Cleanup

If you created a test repository:

```bash
ssh hpc3 "rm -rf /share/crsp/lab/YOUR_LAB/YOUR_USER/git-exercise"
```

## Verification Checklist

- [ ] git status auto-approved
- [ ] git log auto-approved
- [ ] git add auto-approved
- [ ] git commit prompted for confirmation
- [ ] git push would prompt (if tested)
- [ ] git push --force showed strong warning

## What You Learned

1. **Read-only git is safe** - status, log, diff auto-approve
2. **Staging is safe** - git add auto-approves
3. **Commits need confirmation** - Creates permanent history
4. **Force operations warn strongly** - Can cause data loss
5. **Tiering makes sense** - Risk-appropriate guardrails

## Why This Matters

Git on shared storage (CRSP) affects:
- Collaborators who may have cloned
- Remote repositories others push to
- History that can't be easily recovered

The tiering prevents accidental:
- Commits to wrong branch
- Pushes to wrong remote
- Force-pushes that overwrite history

## Common Issues

### "Not a git repository"
- Check path is correct
- Ensure `.git` directory exists

### Prompts for safe commands
- Check hook is correctly detecting command type
- May be picking up a different command in chain

## Next Exercise

→ [04_create_experiment.md](04_create_experiment.md) - Initialize an experiment folder
