# Prerequisites

Before setting up Claude Code with HPC integration, ensure you have:

## 1. UCI HPC3 Account

**Request at:** https://rcic.uci.edu/account.html

You'll need:
- UCINetID (your UCI login)
- Faculty sponsor (PI) approval
- Lab allocation (for SLURM account)

**Verify:** Try logging in at https://hpc3.rcic.uci.edu via Open OnDemand.

## 2. SSH Key Pair

Generate a key pair if you don't have one:

```bash
# Generate ED25519 key (recommended)
ssh-keygen -t ed25519 -C "your_email@uci.edu"

# Or RSA if ED25519 not supported
ssh-keygen -t rsa -b 4096 -C "your_email@uci.edu"
```

**Default locations:**
- Private key: `~/.ssh/id_ed25519`
- Public key: `~/.ssh/id_ed25519.pub`

**Add to HPC3:**
1. Log in to https://hpc3.rcic.uci.edu (Open OnDemand)
2. Open Terminal
3. `cat >> ~/.ssh/authorized_keys` and paste your public key
4. Or use `ssh-copy-id` if you can password-login first

## 3. CRSP Storage Access

**Request at:** https://rcic.uci.edu/storage/crsp.html

CRSP provides:
- Persistent storage for research data
- Accessible from HPC3 at `/share/crsp/lab/{lab_name}`
- Mountable on local machines via CRSP Desktop app

**Verify:** Check your allocation exists:
```bash
ssh hpc3 "ls -la /share/crsp/lab/your_lab_name"
```

## 4. DUO Two-Factor Authentication

**Setup at:** https://www.oit.uci.edu/duo/

Required for SSH authentication:
- Install DUO Mobile on your phone
- Register your device with UCI

**How it works:**
- First SSH connection triggers DUO push notification
- Approve on phone → connection established
- Socket persists for 24h (no repeated auth)

## 5. Claude Code License

**Obtain from:** https://claude.ai/

Install options:
```bash
# npm
npm install -g @anthropic-ai/claude-code

# homebrew (macOS)
brew install claude-code
```

**Verify:**
```bash
claude --version
claude --help
```

## 6. Local Tools

Ensure you have:
- **Python 3.8+** - For hook scripts
- **Bash** - Shell scripts
- **Git** - Version control
- **Text editor** - VS Code, vim, etc.

## Checklist

Before proceeding to setup:

- [ ] HPC3 account active
- [ ] SSH key generated
- [ ] Public key added to HPC3
- [ ] CRSP allocation confirmed
- [ ] DUO Mobile installed and registered
- [ ] Claude Code installed
- [ ] Python 3 available (`python3 --version`)

## Next Step

→ [01_ssh_config.md](01_ssh_config.md) - Configure SSH for HPC access
