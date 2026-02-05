# SSH Configuration

Configure SSH for persistent HPC3 connections with socket multiplexing.

## Create Socket Directory

```bash
mkdir -p ~/.ssh/sockets
chmod 700 ~/.ssh/sockets
```

## SSH Config File

Create or edit `~/.ssh/config`:

```bash
# Main HPC3 connection (interactive use)
Host hpc3
    HostName hpc3.rcic.uci.edu
    User YOUR_UCINETID
    IdentityFile ~/.ssh/id_ed25519
    ControlMaster auto
    ControlPath ~/.ssh/sockets/%r@%h-%p
    ControlPersist 24h
    ServerAliveInterval 60

# Headless mode (for scripts, non-interactive)
Host hpc3-headless
    HostName hpc3.rcic.uci.edu
    User YOUR_UCINETID
    IdentityFile ~/.ssh/id_ed25519
    BatchMode yes
    ControlMaster auto
    ControlPath ~/.ssh/sockets/%r@%h-%p
    ControlPersist 24h
```

**Replace `YOUR_UCINETID` with your actual UCINetID.**

## What Each Setting Does

| Setting | Purpose |
|---------|---------|
| `HostName` | Actual server address |
| `User` | Your login username |
| `IdentityFile` | Path to your private key |
| `ControlMaster auto` | Enable socket multiplexing |
| `ControlPath` | Where to store the socket file |
| `ControlPersist 24h` | Keep socket alive for 24 hours |
| `ServerAliveInterval` | Ping server every 60s to keep alive |
| `BatchMode yes` | Fail immediately if auth needed (headless) |

## First Connection

```bash
# This triggers DUO authentication
ssh hpc3

# You should see:
# Duo two-factor login for YOUR_UCINETID
# Enter a passcode or select one of the following options:
#  1. Duo Push to XXX-XXX-XXXX
```

1. Press `1` for Duo Push
2. Approve on your phone
3. You're connected!

## Verify Socket Works

After first connection:

```bash
# Check socket exists
ls -la ~/.ssh/sockets/

# Test quick reconnect (should be instant, no DUO)
ssh hpc3 hostname

# Expected output: hpc3-xx-xx.rcic.uci.edu
```

## Managing Sockets

```bash
# Check active connections
ssh -O check hpc3

# Close connection (removes socket)
ssh -O exit hpc3

# Force new connection
rm ~/.ssh/sockets/*@hpc3*
ssh hpc3
```

## Off-Campus Access

If connecting from off-campus, you may need UCI VPN:

1. Install GlobalProtect: https://www.oit.uci.edu/vpn/
2. Connect to `vpn.uci.edu`
3. Then `ssh hpc3` works normally

Alternatively, use ProxyJump through a campus machine:

```bash
Host hpc3-via-gateway
    HostName hpc3.rcic.uci.edu
    User YOUR_UCINETID
    ProxyJump gateway.uci.edu
    # ... other settings
```

## Troubleshooting

### "Connection refused"
- Check VPN if off-campus
- Verify hostname: `ping hpc3.rcic.uci.edu`

### "Permission denied (publickey)"
- Verify key exists: `ls ~/.ssh/id_ed25519`
- Check key is in authorized_keys on HPC3
- Try verbose mode: `ssh -v hpc3`

### "Stale socket" errors
```bash
rm ~/.ssh/sockets/*@hpc3*
ssh hpc3
```

### DUO not prompting
- Check phone has Duo Mobile
- Verify device registered at https://www.oit.uci.edu/duo/
- Try password login first: `ssh YOUR_UCINETID@hpc3.rcic.uci.edu`

## Verify Complete

Run these tests:

```bash
# 1. Quick connection test
ssh hpc3 hostname
# Expected: hpc3-xx-xx.rcic.uci.edu

# 2. Check SLURM access
ssh hpc3 "squeue -u \$USER"
# Expected: (empty or your jobs)

# 3. Check CRSP access
ssh hpc3 "ls /share/crsp/lab"
# Expected: list of lab directories
```

## Next Step

→ [02_claude_code_install.md](02_claude_code_install.md) - Install Claude Code CLI
