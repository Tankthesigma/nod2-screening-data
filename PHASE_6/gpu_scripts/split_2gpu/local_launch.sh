#!/bin/bash
#===============================================================================
# LOCAL RTX 4060 Ti LAUNCH SCRIPT
# Runs 4 simulations (80ns) - NO shutdown
#===============================================================================

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
# PREFLIGHT CHECKS
#===============================================================================
echo "Preflight checks..."

# Check OpenMM
python3 -c "import openmm; print(f'OpenMM: {openmm.__version__}')" || {
    echo "ERROR: OpenMM not installed!"
    echo "Install with: conda install -c conda-forge openmm openmmforcefields openff-toolkit"
    exit 1
}

# Check CUDA
python3 -c "from openmm import Platform; p = Platform.getPlatformByName('CUDA'); print('CUDA: OK')" || {
    echo "ERROR: CUDA not available!"
    exit 1
}

# Check GPU
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
echo "STARTING 4 SIMULATIONS (80ns total)"
echo "Expected runtime: ~8-10 hours on RTX 4060 Ti"
echo "============================================================"
echo ""

# Run the Python script
python3 gpu_scripts/split_2gpu/local_4060.py

echo ""
echo "============================================================"
echo "ALL DONE! $(date)"
echo "============================================================"
