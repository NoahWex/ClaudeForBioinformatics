# Troubleshooting

Common issues and solutions.

## SSH Issues

### "Connection refused"

**Cause:** Can't reach HPC3 server

**Solutions:**
1. Check VPN if off-campus
   ```bash
   # Connect to UCI VPN first
   ```

2. Verify server reachable
   ```bash
   ping hpc3.rcic.uci.edu
   ```

3. Check for maintenance
   - Visit https://rcic.uci.edu/status/

### "Permission denied (publickey)"

**Cause:** SSH key not accepted

**Solutions:**
1. Verify key exists
   ```bash
   ls -la ~/.ssh/id_ed25519*
   ```

2. Check key permissions
   ```bash
   chmod 600 ~/.ssh/id_ed25519
   chmod 644 ~/.ssh/id_ed25519.pub
   ```

3. Verify key on HPC3
   ```bash
   # Login via Open OnDemand
   # Check ~/.ssh/authorized_keys contains your public key
   ```

4. Debug connection
   ```bash
   ssh -v hpc3.rcic.uci.edu
   ```

### "Stale socket" or connection hangs

**Cause:** Old socket file from previous session

**Solutions:**
```bash
# Remove stale sockets
rm ~/.ssh/sockets/*@hpc3*

# Reconnect
ssh hpc3
```

### DUO not prompting

**Cause:** Phone not registered or connection issue

**Solutions:**
1. Check Duo Mobile installed and registered
2. Try password login first
   ```bash
   ssh YOUR_UCINETID@hpc3.rcic.uci.edu
   # Enter password, then DUO should prompt
   ```

3. Re-register device at https://www.oit.uci.edu/duo/

## Plugin Issues

### "Command not found" for hpc toolkit

**Cause:** Toolkit not installed or not executable

**Solutions:**
1. Verify toolkit location
   ```bash
   ls ~/.claude/hpc-toolkit/bin/hpc
   ```

2. Make executable
   ```bash
   chmod +x ~/.claude/hpc-toolkit/bin/hpc
   chmod +x ~/.claude/hpc-toolkit/hooks/gate.sh
   ```

3. Re-copy from sandbox
   ```bash
   cp global-setup/hpc-toolkit/bin/hpc ~/.claude/hpc-toolkit/bin/
   cp global-setup/hpc-toolkit/hooks/gate.sh ~/.claude/hpc-toolkit/hooks/
   ```

### Hook not triggering

**Cause:** settings.json misconfigured

**Solutions:**
1. Validate JSON syntax
   ```bash
   python3 -c "import json; json.load(open('$HOME/.claude/settings.json'))"
   ```

2. Check hook path exists
   ```bash
   ls ~/.claude/hpc-toolkit/hooks/gate.sh
   ```

3. Make hook executable
   ```bash
   chmod +x ~/.claude/hpc-toolkit/hooks/gate.sh
   ```

4. Verify settings structure
   ```json
   {
     "hooks": {
       "PreToolUse": [{
         "matcher": "Bash",
         "hooks": [{
           "type": "command",
           "command": "~/.claude/hpc-toolkit/hooks/gate.sh",
           "timeout": 5
         }]
       }]
     }
   }
   ```

### Raw sbatch not blocked

**Cause:** Hook not matching command

**Solutions:**
1. Verify hook registered (see above)

2. Test hook directly
   ```bash
   echo '{"tool_input":{"command":"ssh hpc3 sbatch test.sh"}}' | \
     ~/.claude/hpc-toolkit/hooks/gate.sh
   ```

3. Check hook output for errors

### Config not loading

**Cause:** config.sh missing or has syntax errors

**Solutions:**
1. Verify config exists
   ```bash
   cat ~/.claude/hpc-toolkit/config.sh
   ```

2. Test sourcing it
   ```bash
   source ~/.claude/hpc-toolkit/config.sh && echo "SSH_ALIAS=$SSH_ALIAS"
   ```

3. Create from template if missing
   ```bash
   cp global-setup/hpc-toolkit/config.sh.template ~/.claude/hpc-toolkit/config.sh
   ```

## Claude Code Issues

### Claude not starting

**Solutions:**
1. Check installation
   ```bash
   claude --version
   ```

2. Re-authenticate
   ```bash
   claude --login
   ```

3. Clear cache
   ```bash
   rm -rf ~/.claude/cache
   ```

### Settings not applied

**Solutions:**
1. Verify file location
   - Global: `~/.claude/settings.json`
   - Project: `.claude/settings.local.json`

2. Check JSON syntax
   ```bash
   python3 -c "import json; json.load(open('.claude/settings.local.json'))"
   ```

3. Restart Claude Code

## HPC Issues

### "No allocation" errors

**Cause:** SLURM account not set or wrong

**Solutions:**
1. Check available accounts
   ```bash
   ssh hpc3 "sacctmgr show assoc where user=\$USER -P"
   ```

2. Use correct account in scripts
   ```bash
   #SBATCH --account=your_actual_account
   ```

### Jobs stuck in pending

**Cause:** Resource constraints or priority

**Solutions:**
1. Check reason
   ```bash
   ssh hpc3 "squeue -u \$USER -o '%i %T %r'"
   ```

2. Common reasons:
   - `Priority` - Wait for higher priority jobs
   - `Resources` - Requested resources not available
   - `QOSMaxJobsPerUser` - Too many jobs queued

### Container errors

**Cause:** Mount or permission issues

**Solutions:**
1. Always load singularity
   ```bash
   module load singularity
   ```

2. Check bind mounts exist
   ```bash
   ls /share/crsp/lab/your_lab
   ```

3. Don't bind `/scratch` (not on all nodes)

## Getting Help

If issues persist:

1. **UCI RCIC Help**
   - Email: hpc-support@rcic.uci.edu
   - Tickets: https://rcic.uci.edu/support/

2. **Claude Code Issues**
   - GitHub: https://github.com/anthropics/claude-code/issues

3. **This Sandbox**
   - Check examples in exercises/
   - Review setup-guides/ in order
