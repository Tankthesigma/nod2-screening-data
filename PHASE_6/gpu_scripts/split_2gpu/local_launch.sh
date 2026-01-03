#!/bin/bash
#===============================================================================
# LOCAL RTX 4060 Ti LAUNCH SCRIPT
# Runs 4 simulations (80ns) - NO shutdown
#===============================================================================

set -e

echo "============================================================"
echo "LOCAL RTX 4060 Ti - NOD2-SCOUT MD SIMULATIONS"
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
# ENVIRONMENT SETUP (CRITICAL!)
#===============================================================================
echo "Setting up conda environment..."

# Ensure conda is available
if command -v conda &> /dev/null; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
else
    echo "ERROR: conda not found!"
    echo "Install Miniconda: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

# Check if nod2md environment exists
if conda env list | grep -q "nod2md"; then
    echo "Activating nod2md environment..."
    conda activate nod2md
else
    echo "ERROR: 'nod2md' conda environment not found!"
    echo ""
    echo "Create it with:"
    echo "  conda create -n nod2md python=3.10 -y"
    echo "  conda activate nod2md"
    echo "  conda install -c conda-forge openmm openmmforcefields openff-toolkit -y"
    exit 1
fi

#===============================================================================
# PREFLIGHT CHECKS
#===============================================================================
echo ""
echo "Preflight checks..."

# Check OpenMM
python3 -c "import openmm; print(f'  OpenMM: {openmm.__version__}')" || {
    echo "ERROR: OpenMM not installed!"
    echo "Run: conda install -c conda-forge openmm openmmforcefields openff-toolkit -y"
    exit 1
}

# Check OpenFF (REQUIRED - not optional!)
python3 -c "from openff.toolkit import Molecule; print('  OpenFF: OK')" || {
    echo "ERROR: OpenFF not installed! This is REQUIRED."
    echo "Run: conda install -c conda-forge openff-toolkit openmmforcefields -y"
    exit 1
}

# Check CUDA
python3 -c "from openmm import Platform; p = Platform.getPlatformByName('CUDA'); print('  CUDA: OK')" || {
    echo "ERROR: CUDA not available!"
    echo "Make sure you have an NVIDIA GPU with CUDA drivers"
    exit 1
}

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
echo "  Done"

#===============================================================================
# PREVENT SLEEP (reminder)
#===============================================================================
echo ""
echo "!!! IMPORTANT: Disable sleep mode to prevent interruption !!!"
echo "    Windows: Power Settings -> Sleep: Never"
echo "    Linux: sudo systemctl mask sleep.target suspend.target"
echo ""

#===============================================================================
# RUN SIMULATIONS
#===============================================================================
echo "============================================================"
echo "STARTING 4 SIMULATIONS (80ns total)"
echo "Expected runtime: ~8-10 hours"
echo "============================================================"
echo ""

# Run the Python script
python3 gpu_scripts/split_2gpu/local_4060.py

echo ""
echo "============================================================"
echo "ALL DONE! $(date)"
echo "============================================================"
