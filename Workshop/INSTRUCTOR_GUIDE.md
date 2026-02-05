# Instructor Guide

This guide provides timing cues, common issues, and solutions for running the Claude for Bioinformatics workshop.

## Pre-Workshop Checklist

### One Week Before
- [ ] Verify HPC3 access for all participants
- [ ] Test all checkpoint data files
- [ ] Confirm container accessibility on HPC
- [ ] Send pre-workshop email with setup instructions
- [ ] Prepare backup presentation materials (PDF slides)

### Day Before
- [ ] Run `verify_setup.sh` on instructor machine
- [ ] Test Claude Code authentication
- [ ] Verify data files are accessible
- [ ] Prepare backup Claude API keys (if providing)
- [ ] Test screen sharing setup

### Day Of (30 min before)
- [ ] Open workshop materials
- [ ] Have checkpoint recovery ready
- [ ] Test projector/screen share
- [ ] Have backup slides open

## Session Timing

### Day 1: Foundations (2.5 hours)

| Time | Duration | Activity | Notes |
|------|----------|----------|-------|
| 0:00 | 25 min | Intro to Claude Code | Demo: /help, basic tools |
| 0:25 | 35 min | Exercise 1: Scaffolding | Watch for path issues |
| 1:00 | 10 min | Break | |
| 1:10 | 45 min | Exercise 2: Data loading | Main hands-on |
| 1:55 | 15 min | Checkpoint review | Verify everyone has output |
| 2:10 | 20 min | Q&A / Buffer | Use if running behind |

**Critical checkpoint**: Everyone should have `manifest.yaml` and loaded data by end.

### Day 2: Metadata (2.5 hours)

| Time | Duration | Activity | Notes |
|------|----------|----------|-------|
| 0:00 | 20 min | Review Day 1, metadata intro | Quick demo of Grep |
| 0:20 | 40 min | Exercise 1: Exploration | Compare schemas |
| 1:00 | 10 min | Break | |
| 1:10 | 25 min | Demo: Harmonization strategy | Show cell type dictionary |
| 1:35 | 35 min | Exercise 2: Build mapping | Main hands-on |
| 2:10 | 20 min | Checkpoint review | Verify harmonization |

**Critical checkpoint**: Everyone should have harmonized metadata mapping.

### Day 3: HPC (3 hours)

| Time | Duration | Activity | Notes |
|------|----------|----------|-------|
| 0:00 | 15 min | Container concepts | Quick module demo |
| 0:15 | 40 min | Exercise 1: QC script | Test mode first |
| 0:55 | 10 min | Break | |
| 1:05 | 35 min | Exercise 2: SLURM jobs | Walk through carefully |
| 1:40 | 30 min | Exercise 3: Array jobs | May need extra time |
| 2:10 | 20 min | Submit & monitor | Live submission |
| 2:30 | 20 min | Checkpoint review | Jobs may still be running |
| 2:50 | 10 min | Wrap-up | |

**Critical checkpoint**: Jobs submitted, understand monitoring.

**Note**: Jobs may not complete during session - that's OK.

### Day 4: Integration (3 hours)

| Time | Duration | Activity | Notes |
|------|----------|----------|-------|
| 0:00 | 20 min | Integration overview | Harmony vs Seurat |
| 0:20 | 35 min | Exercise 1: Integration | Main compute step |
| 0:55 | 10 min | Break | |
| 1:05 | 25 min | Exercise 2: Visualization | Iterative refinement |
| 1:30 | 25 min | Exercise 3: Reports | Quarto basics |
| 1:55 | 35 min | Exercise 4: Customization | Rules demo |
| 2:30 | 30 min | Wrap-up, Q&A | Resources, next steps |

**Critical checkpoint**: Integrated object + one visualization.

## Common Issues and Solutions

### Day 1

**Issue**: Claude doesn't respond
- **Solution**: Check `claude config`, verify API key
- **Backup**: Provide workshop API key

**Issue**: File paths with spaces cause errors
- **Solution**: Always quote paths or use tab completion
- **Example**: `cd "/path/with spaces/"`

**Issue**: Can't load RDS file
- **Solution**: Check path is absolute, file exists
- **Backup**: Use checkpoint `cp1_data_loaded`

### Day 2

**Issue**: Metadata column names differ unexpectedly
- **Solution**: Use provided cell type dictionary as reference
- **Note**: Real data is messy - this is a teaching moment

**Issue**: Edit tool fails
- **Solution**: Ensure file was Read first
- **Backup**: Rewrite entire file with Write tool

### Day 3

**Issue**: SSH to HPC fails
- **Solution**: Check VPN, SSH key, hostname
- **Backup**: Demo on your machine, provide commands to copy

**Issue**: Container not found
- **Solution**: Verify path: `ls /share/crsp/lab/dalawson/share/0_Resources/containers/`
- **Backup**: Use module system instead

**Issue**: Job stuck in PENDING
- **Solution**: `squeue -j <id> -o "%R"` to see reason
- **Common**: Resource limits, wrong account

**Issue**: Out of memory
- **Solution**: Increase `--mem`, reduce data size
- **Backup**: Use test mode with 1000 cells

### Day 4

**Issue**: Harmony installation
- **Solution**: Container has it; local install: `install.packages("harmony")`

**Issue**: Integration takes too long
- **Solution**: Use subset for demo; full run overnight
- **Backup**: Provide pre-integrated checkpoint

**Issue**: Quarto not installed
- **Solution**: Use container or `module load quarto`
- **Backup**: Show R Markdown alternative

## Teaching Tips

### Pacing
- Start each exercise with a 2-minute live demo
- Give 5-minute warnings before transitions
- Use checkpoints liberally - falling behind is OK

### Engagement
- Ask participants to share their Claude outputs
- Discuss different approaches to same problem
- Encourage "bad" prompts - learning what doesn't work

### Common Participant Struggles
1. **Prompt specificity** - Too vague or too detailed
2. **File paths** - Relative vs absolute confusion
3. **Tool selection** - Not knowing which tool to use
4. **Iterative refinement** - Expecting perfect first answer

### Demonstrating Failures
Intentionally show:
- What happens with a bad prompt
- How to iterate when Claude misunderstands
- Error recovery strategies

## Checkpoint Recovery

If participants fall significantly behind:

```bash
# Quick reset
cd Workshop/recovery
./sync_checkpoint.sh cp2_harmonized ./working
```

Don't hesitate to use checkpoints - workshop value is in learning the workflow, not completing every step manually.

## Materials Checklist

### Per Session
- [ ] README.md printed/available
- [ ] Exercises open in tabs
- [ ] Solutions ready (don't share until needed)
- [ ] Checkpoint data verified

### Backup Materials
- [ ] PDF slides (if screen share fails)
- [ ] Pre-written commands in text file
- [ ] Checkpoint files on multiple locations

## Post-Workshop

### Immediately After
- [ ] Collect feedback (survey link)
- [ ] Note issues for next time
- [ ] Thank participants

### Within One Week
- [ ] Review feedback
- [ ] Update materials based on issues
- [ ] Send follow-up resources email

## Emergency Procedures

### Claude Service Down
- Switch to pre-recorded demos
- Focus on conceptual discussion
- Provide solutions for self-study

### HPC Issues
- Use local R for non-HPC exercises
- Demo SLURM concepts without submission
- Provide job outputs from previous runs

### Data Access Issues
- Have local copies of small subsets
- Use synthetic data for demos
- Checkpoint files include processed data

## Contact

Workshop organizer: [Your contact info]
Technical issues: [Lab Slack/email]
Claude Code issues: https://github.com/anthropics/claude-code/issues
