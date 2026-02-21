#!/bin/bash
# Launch all 8 GPUs in parallel
# Run setup_once.sh FIRST!

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

# Verify ALL dependencies were installed
echo "Checking dependencies..."
if ! python -c "import openmm" 2>/dev/null; then
    echo "ERROR: OpenMM not installed. Run ./setup_once.sh first!"
    exit 1
fi
if ! python -c "import numpy" 2>/dev/null; then
    echo "ERROR: NumPy not installed. Run ./setup_once.sh first!"
    exit 1
fi
if ! python -c "import openmmtools" 2>/dev/null; then
    echo "ERROR: OpenMMTools not installed. Run ./setup_once.sh first!"
    exit 1
fi
echo "All dependencies OK"

if [ ! -d "data/fep_pmx" ]; then
    echo "ERROR: Data not found at data/fep_pmx"
    echo "Did you unpack and move fep_pmx to data/fep_pmx?"
    exit 1
fi

echo ""
echo "=== Launching 8 GPUs ==="

# Launch each GPU in background
for i in {0..7}; do
    echo "Launching GPU $i..."
    nohup bash deploy.sh $i > gpu_${i}.log 2>&1 &
done

echo ""
echo "All 8 GPUs launched!"
echo ""
echo "Monitor with:"
echo "  tail -f gpu_0.log"
echo "  grep 'DONE' gpu_*.log | wc -l"
echo "  grep 'COMPLETE' gpu_*.log"
