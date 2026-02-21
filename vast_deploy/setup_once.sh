#!/bin/bash
# Run this ONCE before launching GPUs
# Sets up environment for all GPUs to share

set -e

# Use script directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

echo "=== One-time setup ==="
echo "Working dir: $SCRIPT_DIR"

# Use python -m pip to ensure same environment
python -m pip install --upgrade pip
python -m pip install numpy openmm openmmtools

# Verify ALL dependencies
echo "Verifying installations..."
python -c "import openmm; print(f'OpenMM: {openmm.__version__}')"
python -c "import numpy; print(f'NumPy: {numpy.__version__}')"
python -c "import openmmtools; print(f'OpenMMTools: {openmmtools.__version__}')"

# Check CUDA platform
python -c "from openmm import Platform; p = Platform.getPlatformByName('CUDA'); print('CUDA platform: OK')"

echo ""
echo "=== Setup complete ==="
echo "Now run: ./run_all_gpus.sh"
