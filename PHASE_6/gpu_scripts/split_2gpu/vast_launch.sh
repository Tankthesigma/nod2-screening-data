#!/bin/bash
#===============================================================================
# VAST.AI RTX 5090 LAUNCH SCRIPT
# Runs 10 simulations (200ns) then auto-shuts down
#===============================================================================

set -e

echo "============================================================"
echo "VAST.AI RTX 5090 - NOD2-SCOUT MD SIMULATIONS"
echo "============================================================"
echo ""
echo "Started: $(date)"
echo ""

# Get script directory and navigate to PHASE_6
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR/../.."
echo "Working directory: $(pwd)"
echo ""

#===============================================================================
# ENVIRONMENT SETUP
#===============================================================================
echo "Setting up environment..."

# CRITICAL: Source conda FIRST so 'conda' command works
# Try common conda locations
if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
    source "$HOME/anaconda3/etc/profile.d/conda.sh"
elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
    source "/opt/conda/etc/profile.d/conda.sh"
elif command -v conda &> /dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
else
    echo "ERROR: conda not found! Install Miniconda first."
    exit 1
fi

# Create conda environment if it doesn't exist
if ! conda env list | grep -q "nod2md"; then
    echo "Creating conda environment..."
    conda create -n nod2md python=3.10 -y
fi

# Activate environment
echo "Activating nod2md environment..."
conda activate nod2md

# Install dependencies
echo "Installing dependencies..."
conda install -c conda-forge openmm openmmforcefields openff-toolkit -y
pip install mdtraj

#===============================================================================
# PREFLIGHT CHECKS
#===============================================================================
echo ""
echo "Preflight checks..."

# Check OpenMM
python3 -c "import openmm; print(f'  OpenMM: {openmm.__version__}')"

# Check OpenFF (REQUIRED)
python3 -c "from openff.toolkit import Molecule; print('  OpenFF: OK')" || {
    echo "ERROR: OpenFF not installed!"
    exit 1
}

# Check CUDA
python3 -c "from openmm import Platform; p = Platform.getPlatformByName('CUDA'); print('  CUDA: OK')"

# Check GPU
echo ""
echo "GPU Info:"
nvidia-smi --query-gpu=name,memory.total --format=csv

#===============================================================================
# CREATE OUTPUT DIRECTORIES
#===============================================================================
echo ""
echo "Creating output directories..."
mkdir -p trajectories logs checkpoints

#===============================================================================
# RUN SIMULATIONS
#===============================================================================
echo ""
echo "============================================================"
echo "STARTING 10 SIMULATIONS (200ns total)"
echo "Expected runtime: ~5 hours on RTX 5090"
echo "============================================================"
echo ""

# Run the Python script
python3 gpu_scripts/split_2gpu/vast_5090.py

# Note: The Python script handles shutdown automatically
