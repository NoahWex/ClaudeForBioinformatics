#!/bin/bash
# setup.sh - Claude Code HPC Toolkit Setup
#
# Reads config/sandbox.yaml and installs the hpc-toolkit, skills, and rules.
# Run from the sandbox directory.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/config/sandbox.yaml"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# ==============================================================================
# Helper Functions
# ==============================================================================

log_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

log_success() {
    echo -e "${GREEN}[OK]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

check_command() {
    command -v "$1" &> /dev/null
}

# Parse YAML (simple grep-based, works for flat config)
get_config() {
    local key="$1"
    grep -E "^[[:space:]]*${key}:" "$CONFIG_FILE" 2>/dev/null | sed 's/.*:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo ""
}

get_nested_config() {
    local section="$1"
    local key="$2"
    awk "/^${section}:$/,/^[a-z]/" "$CONFIG_FILE" 2>/dev/null | grep -E "^[[:space:]]+${key}:" | sed 's/.*:[[:space:]]*//' | tr -d '"' | tr -d "'" || echo ""
}

# ==============================================================================
# Usage
# ==============================================================================

show_help() {
    cat << 'EOF'
Claude Code HPC Toolkit Setup

Usage: ./setup.sh [OPTIONS]

OPTIONS:
  --help, -h      Show this help
  --verify        Only verify setup, don't install
  --skip-ssh      Skip SSH configuration check
  --skip-toolkit  Skip toolkit installation
  --skip-rules    Skip global rules installation

CONFIGURATION:
  Edit config/sandbox.yaml before running.

STEPS:
  1. Edit config/sandbox.yaml with your settings
  2. Run ./setup.sh
  3. Edit ~/.claude/hpc-toolkit/config.sh for your environment
  4. Run ./setup.sh --verify to confirm

EOF
    exit 0
}

# ==============================================================================
# Verification Functions
# ==============================================================================

verify_ssh() {
    local ssh_alias
    ssh_alias=$(get_nested_config "cluster" "ssh_alias")
    [[ -z "$ssh_alias" ]] && ssh_alias="hpc3"

    log_info "Checking SSH connection to $ssh_alias..."

    if ssh -o BatchMode=yes -o ConnectTimeout=5 "$ssh_alias" hostname &>/dev/null; then
        log_success "SSH connection to $ssh_alias works"
        return 0
    else
        log_warn "SSH connection failed - may need DUO authentication"
        echo "    Try: ssh $ssh_alias"
        return 1
    fi
}

verify_claude() {
    log_info "Checking Claude Code installation..."

    if check_command claude; then
        local version
        version=$(claude --version 2>/dev/null || echo "unknown")
        log_success "Claude Code installed: $version"
        return 0
    else
        log_error "Claude Code not found"
        echo "    Install: npm install -g @anthropic-ai/claude-code"
        return 1
    fi
}

verify_toolkit() {
    log_info "Checking HPC toolkit installation..."

    local toolkit_path="$HOME/.claude/hpc-toolkit"

    if [[ -x "$toolkit_path/bin/hpc" ]]; then
        log_success "hpc-toolkit installed and executable"
    else
        log_warn "hpc-toolkit not installed or not executable"
        echo "    Run: ./setup.sh (without --skip-toolkit)"
        return 1
    fi

    if [[ -x "$toolkit_path/hooks/gate.sh" ]]; then
        log_success "gate.sh hook installed"
    else
        log_warn "gate.sh not found or not executable"
        return 1
    fi

    if [[ -f "$toolkit_path/config.sh" ]]; then
        log_success "config.sh exists"
    else
        log_warn "config.sh not found — copy from config.sh.template and edit"
        return 1
    fi

    return 0
}

verify_hooks() {
    log_info "Checking hooks configuration..."

    local settings_file="$HOME/.claude/settings.json"

    if [[ -f "$settings_file" ]]; then
        if grep -q "gate.sh" "$settings_file" 2>/dev/null; then
            log_success "Gate hook registered in settings.json"
            return 0
        else
            log_warn "gate.sh hook not found in settings.json"
            return 1
        fi
    else
        log_warn "settings.json not found"
        echo "    Create: $settings_file"
        return 1
    fi
}

verify_python() {
    log_info "Checking Python 3..."

    if check_command python3; then
        log_success "Python 3 available"
        return 0
    else
        log_error "Python 3 not found (required for gate.sh)"
        return 1
    fi
}

# ==============================================================================
# Installation Functions
# ==============================================================================

install_toolkit() {
    local toolkit_src="$SCRIPT_DIR/global-setup/hpc-toolkit"
    local toolkit_dest="$HOME/.claude/hpc-toolkit"

    log_info "Installing HPC toolkit..."

    # Create directories
    mkdir -p "$toolkit_dest/bin"
    mkdir -p "$toolkit_dest/hooks"
    mkdir -p "$toolkit_dest/docs"
    mkdir -p "$toolkit_dest/logs"

    # Copy toolkit files
    if [[ -d "$toolkit_src" ]]; then
        cp "$toolkit_src/bin/hpc" "$toolkit_dest/bin/"
        cp "$toolkit_src/hooks/gate.sh" "$toolkit_dest/hooks/"
        cp "$toolkit_src/docs/CLAUDE_GUIDE.md" "$toolkit_dest/docs/"
        chmod +x "$toolkit_dest/bin/hpc"
        chmod +x "$toolkit_dest/hooks/gate.sh"
        log_success "Toolkit installed to $toolkit_dest"
    else
        log_error "Toolkit source not found: $toolkit_src"
        return 1
    fi

    # Copy config template if config doesn't exist
    if [[ ! -f "$toolkit_dest/config.sh" ]]; then
        cp "$toolkit_src/config.sh.template" "$toolkit_dest/config.sh"
        log_warn "config.sh created from template — EDIT IT for your environment"
        echo "    File: $toolkit_dest/config.sh"
    else
        log_success "config.sh already exists (not overwritten)"
    fi
}

install_skills() {
    local skills_src="$SCRIPT_DIR/global-setup/skills"
    local skills_dest="$HOME/.claude/skills"

    log_info "Installing HPC skill..."

    mkdir -p "$skills_dest/hpc"

    if [[ -f "$skills_src/hpc/SKILL.md" ]]; then
        cp "$skills_src/hpc/SKILL.md" "$skills_dest/hpc/"
        log_success "HPC skill installed to $skills_dest/hpc/"
    else
        log_warn "HPC skill source not found"
    fi
}

install_rules() {
    local rules_src="$SCRIPT_DIR/global-setup/rules"
    local rules_dest="$HOME/.claude/rules"

    log_info "Installing global rules..."

    mkdir -p "$rules_dest"

    if [[ -d "$rules_src" ]]; then
        cp "$rules_src"/*.md "$rules_dest/" 2>/dev/null || true
        log_success "Rules installed to $rules_dest"
    else
        log_warn "Rules source not found: $rules_src"
    fi
}

configure_settings() {
    local settings_file="$HOME/.claude/settings.json"
    local settings_template="$SCRIPT_DIR/global-setup/settings.json.template"

    log_info "Configuring settings..."

    # Get SSH alias from config
    local ssh_alias
    ssh_alias=$(get_nested_config "cluster" "ssh_alias")
    [[ -z "$ssh_alias" ]] && ssh_alias="hpc3"

    if [[ -f "$settings_file" ]]; then
        log_info "Existing settings.json found — attempting merge..."

        # Use Python to merge settings (preserves existing permissions, adds new hooks)
        if python3 << MERGE_SCRIPT
import json
import sys

# Load existing settings
with open("$settings_file", "r") as f:
    existing = json.load(f)

# Load template and substitute SSH_ALIAS
with open("$settings_template", "r") as f:
    template_str = f.read().replace("SSH_ALIAS", "$ssh_alias")
    template = json.loads(template_str)

# Merge permissions.allow (deduplicate)
existing_allow = existing.get("permissions", {}).get("allow", [])
template_allow = template.get("permissions", {}).get("allow", [])
merged_allow = list(dict.fromkeys(existing_allow + template_allow))  # preserve order, dedupe

if "permissions" not in existing:
    existing["permissions"] = {}
existing["permissions"]["allow"] = merged_allow

# Merge hooks (add new hooks, don't overwrite existing)
if "hooks" not in existing:
    existing["hooks"] = {}

for hook_type, hook_list in template.get("hooks", {}).items():
    if hook_type not in existing["hooks"]:
        existing["hooks"][hook_type] = hook_list
    else:
        # Check if hook already exists (by command path)
        existing_cmds = set()
        for h in existing["hooks"][hook_type]:
            for inner_hook in h.get("hooks", []):
                existing_cmds.add(inner_hook.get("command", ""))

        for new_hook in hook_list:
            for inner_hook in new_hook.get("hooks", []):
                if inner_hook.get("command", "") not in existing_cmds:
                    existing["hooks"][hook_type].append(new_hook)
                    break

# Write merged settings
with open("$settings_file", "w") as f:
    json.dump(existing, f, indent=2)
    f.write("\n")

print("Merged successfully")
MERGE_SCRIPT
        then
            log_success "Settings merged at $settings_file (SSH alias: $ssh_alias)"
        else
            log_warn "Merge failed — manual merge may be needed"
            echo "    Template: $settings_template"
            echo "    Ensure it has PreToolUse hook for gate.sh and SessionStart for trace-framework"
        fi
    else
        if [[ -f "$settings_template" ]]; then
            # Fresh install — just substitute and copy
            sed "s/SSH_ALIAS/$ssh_alias/g" "$settings_template" > "$settings_file"
            log_success "Settings created at $settings_file (SSH alias: $ssh_alias)"
        else
            log_warn "Settings template not found"
        fi
    fi
}

create_ssh_sockets_dir() {
    log_info "Checking SSH sockets directory..."

    if [[ ! -d "$HOME/.ssh/sockets" ]]; then
        mkdir -p "$HOME/.ssh/sockets"
        chmod 700 "$HOME/.ssh/sockets"
        log_success "Created ~/.ssh/sockets"
    else
        log_success "SSH sockets directory exists"
    fi
}

# ==============================================================================
# Main
# ==============================================================================

VERIFY_ONLY=false
SKIP_SSH=false
SKIP_TOOLKIT=false
SKIP_RULES=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --help|-h)
            show_help
            ;;
        --verify)
            VERIFY_ONLY=true
            shift
            ;;
        --skip-ssh)
            SKIP_SSH=true
            shift
            ;;
        --skip-toolkit)
            SKIP_TOOLKIT=true
            shift
            ;;
        --skip-rules)
            SKIP_RULES=true
            shift
            ;;
        *)
            log_error "Unknown option: $1"
            echo "Run: ./setup.sh --help"
            exit 1
            ;;
    esac
done

# Check config exists
if [[ ! -f "$CONFIG_FILE" ]]; then
    log_error "Config file not found: $CONFIG_FILE"
    echo "    Copy and edit: config/sandbox.yaml"
    exit 1
fi

echo ""
echo "=============================================="
echo "  Claude Code HPC Toolkit Setup"
echo "=============================================="
echo ""

if [[ "$VERIFY_ONLY" == "true" ]]; then
    log_info "Verification mode"
    echo ""

    FAILED=0

    verify_python || ((FAILED++))
    verify_claude || ((FAILED++))
    [[ "$SKIP_SSH" == "false" ]] && (verify_ssh || ((FAILED++)))
    verify_toolkit || ((FAILED++))
    verify_hooks || ((FAILED++))

    echo ""
    if [[ $FAILED -eq 0 ]]; then
        log_success "All checks passed!"
    else
        log_warn "$FAILED check(s) failed"
        echo ""
        echo "Run ./setup.sh to install missing components"
    fi

    exit $FAILED
fi

# Installation mode
log_info "Installation mode"
echo ""

# Prerequisites
verify_python || exit 1
verify_claude || exit 1

# Create directories
create_ssh_sockets_dir

# Install components
[[ "$SKIP_TOOLKIT" == "false" ]] && install_toolkit
[[ "$SKIP_TOOLKIT" == "false" ]] && install_skills
[[ "$SKIP_RULES" == "false" ]] && install_rules
configure_settings

echo ""
log_success "Setup complete!"
echo ""
echo "Next steps:"
echo "  1. EDIT ~/.claude/hpc-toolkit/config.sh for your environment"
echo "  2. Configure SSH: See setup-guides/01_ssh_config.md"
echo "  3. Test SSH: ssh hpc3 hostname"
echo "  4. Verify: ./setup.sh --verify"
echo "  5. Start exercises: exercises/01_first_hpc_command.md"
echo ""
