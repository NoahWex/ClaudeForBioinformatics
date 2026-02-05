#!/bin/bash
# verify_setup.sh - Pre-workshop environment verification
#
# Run this script before the workshop to ensure your environment is ready.
# Usage: ./verify_setup.sh

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKSHOP_DIR="$(dirname "$SCRIPT_DIR")"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

PASS=0
WARN=0
FAIL=0

check_pass() {
    echo -e "  ${GREEN}[PASS]${NC} $1"
    ((PASS++))
}

check_warn() {
    echo -e "  ${YELLOW}[WARN]${NC} $1"
    ((WARN++))
}

check_fail() {
    echo -e "  ${RED}[FAIL]${NC} $1"
    ((FAIL++))
}

echo "=============================================="
echo "Workshop Environment Verification"
echo "=============================================="
echo ""

# ----------------------------------------
echo "1. Claude Code CLI"
echo "----------------------------------------"

if command -v claude &> /dev/null; then
    CLAUDE_VERSION=$(claude --version 2>/dev/null || echo "unknown")
    check_pass "Claude CLI installed: $CLAUDE_VERSION"
else
    check_fail "Claude CLI not found. Install with: npm install -g @anthropic-ai/claude-code"
fi

# ----------------------------------------
echo ""
echo "2. R Environment"
echo "----------------------------------------"

if command -v R &> /dev/null; then
    R_VERSION=$(R --version 2>/dev/null | head -1 || echo "unknown")
    check_pass "R installed: $R_VERSION"
else
    check_fail "R not found"
fi

# Check for Seurat package
if R --quiet -e "packageVersion('Seurat')" 2>/dev/null | grep -q "[0-9]"; then
    SEURAT_VERSION=$(R --quiet -e "cat(as.character(packageVersion('Seurat')))" 2>/dev/null)
    check_pass "Seurat package installed: v$SEURAT_VERSION"
else
    check_warn "Seurat package not found locally (available in HPC container)"
fi

# ----------------------------------------
echo ""
echo "3. HPC Access"
echo "----------------------------------------"

if command -v ssh &> /dev/null; then
    check_pass "SSH available"
else
    check_fail "SSH not found"
fi

# Test HPC connection (timeout after 5 seconds)
if timeout 5 ssh -o BatchMode=yes -o ConnectTimeout=3 hpc3.rcic.uci.edu "echo connected" 2>/dev/null | grep -q "connected"; then
    check_pass "HPC3 SSH connection successful"
else
    check_warn "Could not connect to HPC3 (may need VPN or SSH key setup)"
fi

# ----------------------------------------
echo ""
echo "4. Workshop Files"
echo "----------------------------------------"

if [ -d "$WORKSHOP_DIR" ]; then
    check_pass "Workshop directory exists"
else
    check_fail "Workshop directory not found: $WORKSHOP_DIR"
fi

# Check session directories
for session in session_01_foundations session_02_metadata session_03_hpc session_04_integration; do
    if [ -d "$WORKSHOP_DIR/$session" ]; then
        check_pass "Found $session/"
    else
        check_fail "Missing $session/"
    fi
done

# Check data directory
if [ -d "$WORKSHOP_DIR/data" ]; then
    check_pass "Data directory exists"
else
    check_fail "Data directory not found"
fi

# ----------------------------------------
echo ""
echo "5. Data Files (via manifest)"
echo "----------------------------------------"

# Check data manifest exists
DATA_MANIFEST="$WORKSHOP_DIR/data/workshop_data.yaml"
if [ -f "$DATA_MANIFEST" ]; then
    check_pass "Data manifest exists: workshop_data.yaml"
else
    check_fail "Data manifest not found"
fi

# Check source data at HPC or local path
HPC_DATA="/share/crsp/lab/dalawson/share/3_Downloaded_Datasets/iHBCA_Reed_2024/component_studies"
LOCAL_DATA="$HOME/Library/Application Support/CRSP Desktop/Volumes.noindex/CRSP Lab - dalawson.localized/share/3_Downloaded_Datasets/iHBCA_Reed_2024/component_studies"

if [ -d "$HPC_DATA" ]; then
    DATA_DIR="$HPC_DATA"
    check_pass "HPC data path accessible"
elif [ -d "$LOCAL_DATA" ]; then
    DATA_DIR="$LOCAL_DATA"
    check_pass "Local CRSP data path accessible"
else
    check_warn "Data source not found (need HPC access or CRSP mount)"
    DATA_DIR=""
fi

# Check individual files if path found
if [ -n "$DATA_DIR" ]; then
    for file in gray.rds pal.rds; do
        if [ -e "$DATA_DIR/$file" ]; then
            SIZE=$(du -h "$DATA_DIR/$file" 2>/dev/null | cut -f1)
            check_pass "Found $file ($SIZE)"
        else
            check_warn "Data file not found: $file"
        fi
    done
fi

# ----------------------------------------
echo ""
echo "6. Optional Tools"
echo "----------------------------------------"

# Git
if command -v git &> /dev/null; then
    GIT_VERSION=$(git --version)
    check_pass "$GIT_VERSION"
else
    check_warn "Git not found (optional but recommended)"
fi

# Quarto
if command -v quarto &> /dev/null; then
    QUARTO_VERSION=$(quarto --version 2>/dev/null)
    check_pass "Quarto installed: v$QUARTO_VERSION"
else
    check_warn "Quarto not found (needed for Day 4 report generation)"
fi

# ----------------------------------------
echo ""
echo "=============================================="
echo "Summary"
echo "=============================================="
echo -e "  ${GREEN}Passed:${NC}   $PASS"
echo -e "  ${YELLOW}Warnings:${NC} $WARN"
echo -e "  ${RED}Failed:${NC}   $FAIL"
echo ""

if [ $FAIL -eq 0 ]; then
    if [ $WARN -eq 0 ]; then
        echo -e "${GREEN}All checks passed! You're ready for the workshop.${NC}"
    else
        echo -e "${YELLOW}Some warnings (see above). Core requirements met.${NC}"
    fi
    exit 0
else
    echo -e "${RED}Some checks failed. Please resolve issues before the workshop.${NC}"
    exit 1
fi
