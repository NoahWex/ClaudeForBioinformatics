#!/bin/bash
# Graph Memory System - Installation Verification
# Returns 0 if all checks pass, non-zero on failure
#
# Usage:
#   cd graph-memory && ./tests/verify_install.sh
#
# Environment:
#   PLAN_CMD - Override path to plan script (default: ./plan)

set -euo pipefail

FAILED=0
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
GRAPH_MEMORY_DIR="$(dirname "$SCRIPT_DIR")"

# Default to ./plan relative to graph-memory root
PLAN_CMD="${PLAN_CMD:-$GRAPH_MEMORY_DIR/plan}"

check() {
    local name="$1"
    local cmd="$2"

    if eval "$cmd" &>/dev/null; then
        echo "[OK] $name"
    else
        echo "[FAIL] $name"
        ((FAILED++)) || true
    fi
}

check_output() {
    local name="$1"
    local cmd="$2"
    local expected="$3"

    local output
    if output=$(eval "$cmd" 2>&1); then
        if echo "$output" | grep -q "$expected"; then
            echo "[OK] $name"
        else
            echo "[FAIL] $name (output mismatch)"
            echo "  Expected: $expected"
            echo "  Got: ${output:0:100}..."
            ((FAILED++)) || true
        fi
    else
        echo "[FAIL] $name (command failed)"
        ((FAILED++)) || true
    fi
}

echo "Graph Memory System Verification"
echo "================================="
echo "Graph memory dir: $GRAPH_MEMORY_DIR"
echo "Plan command: $PLAN_CMD"
echo ""

# Core checks
check "Python 3 available" "command -v python3"
check "CLI script exists" "test -x '$PLAN_CMD'"
check "CLI invocation" "'$PLAN_CMD' --help"

# Module checks
check "parse.py exists" "test -f '$GRAPH_MEMORY_DIR/tools/parse.py'"
check "graph.py exists" "test -f '$GRAPH_MEMORY_DIR/tools/graph.py'"
check "knowledge.py exists" "test -f '$GRAPH_MEMORY_DIR/tools/knowledge.py'"
check "visualize.py exists" "test -f '$GRAPH_MEMORY_DIR/tools/visualize.py'"
check "cli.py exists" "test -f '$GRAPH_MEMORY_DIR/tools/cli.py'"

# Import checks
check "parse.py imports" "python3 -c 'import sys; sys.path.insert(0, \"$GRAPH_MEMORY_DIR/tools\"); import parse'"
check "graph.py imports" "python3 -c 'import sys; sys.path.insert(0, \"$GRAPH_MEMORY_DIR/tools\"); import graph'"
check "knowledge.py imports" "python3 -c 'import sys; sys.path.insert(0, \"$GRAPH_MEMORY_DIR/tools\"); import knowledge'"
check "visualize.py imports" "python3 -c 'import sys; sys.path.insert(0, \"$GRAPH_MEMORY_DIR/tools\"); import visualize'"

# Template checks
check "plan_template.md exists" "test -f '$GRAPH_MEMORY_DIR/templates/plan_template.md'"
check "lesson_template.md exists" "test -f '$GRAPH_MEMORY_DIR/templates/lesson_template.md'"

# Validation check with fixture
FIXTURE="$SCRIPT_DIR/fixtures/sample_plan.md"
if [[ -f "$FIXTURE" ]]; then
    check_output "Plan validation" "'$PLAN_CMD' validate '$FIXTURE'" "is valid"
else
    echo "[SKIP] Plan validation (no fixture at $FIXTURE)"
fi

echo ""
echo "---------------------------------"
if [[ $FAILED -eq 0 ]]; then
    echo "All checks passed!"
    exit 0
else
    echo "$FAILED check(s) failed"
    exit 1
fi
