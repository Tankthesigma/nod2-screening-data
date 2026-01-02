#!/bin/bash
#===============================================================================
# NOD2-SCOUT PHASE 6: MASTER MD SIMULATION SCRIPT
#===============================================================================
#
# This script runs the complete MD simulation pipeline on a GPU instance.
#
# Usage:
#   chmod +x master_md_pipeline.sh
#   ./master_md_pipeline.sh
#
# Requirements:
#   - NVIDIA GPU with CUDA support
#   - OpenMM, MDTraj, PDBFixer installed
#   - Prepared structure files in structures/
#
# Estimated time: 12-24 hours on RTX 4090
#===============================================================================

set -e  # Exit on error

echo "=============================================="
echo "NOD2-SCOUT PHASE 6: MD SIMULATION PIPELINE"
echo "=============================================="
echo ""

# Check CUDA
echo "Checking GPU..."
nvidia-smi
echo ""

# Check OpenMM
echo "Checking OpenMM..."
python3 -c "from openmm import Platform; p = Platform.getPlatformByName('CUDA'); print(f'OpenMM CUDA platform: OK')"
echo ""

# Create output directories
mkdir -p trajectories analysis logs

# Compounds to simulate
# FDA drugs + Natural product + Controls
COMPOUNDS="febuxostat ursodiol budesonide natural_cid10592 apo"
N_REPLICATES=3

# Run simulations
for compound in $COMPOUNDS; do
    echo ""
    echo "=============================================="
    echo "SIMULATING: $compound"
    echo "=============================================="

    for rep in $(seq 0 $((N_REPLICATES-1))); do
        echo ""
        echo "--- Replicate $((rep+1))/$N_REPLICATES ---"

        # Check if already completed
        if [ -f "trajectories/${compound}_rep$((rep+1)).dcd" ]; then
            echo "Already completed, skipping..."
            continue
        fi

        # Run simulation
        python3 gpu_scripts/03_run_md_${compound}.py             --pdb structures/complex_${compound}.pdb             --replicate $rep

        # Move outputs
        mv ${compound}_rep*.dcd trajectories/ 2>/dev/null || true
        mv ${compound}_rep*.log logs/ 2>/dev/null || true
        mv ${compound}_rep*.chk trajectories/ 2>/dev/null || true

        echo "Replicate $((rep+1)) complete!"
    done
done

echo ""
echo "=============================================="
echo "RUNNING ANALYSIS"
echo "=============================================="

for compound in $COMPOUNDS; do
    echo "Analyzing $compound..."

    python3 gpu_scripts/04_analyze_${compound}.py         --traj trajectories/${compound}_rep1.dcd         --top structures/complex_${compound}.pdb         --output analysis/
done

echo ""
echo "=============================================="
echo "COMPUTING MM-GBSA"
echo "=============================================="

for compound in $COMPOUNDS; do
    if [ "$compound" != "apo" ]; then
        echo "Computing MM-GBSA for $compound..."

        python3 gpu_scripts/05_mmgbsa_${compound}.py             --pdb structures/complex_${compound}.pdb             --traj trajectories/${compound}_rep1.dcd             --frames 100

        mv ${compound}_mmgbsa.json analysis/ 2>/dev/null || true
    fi
done

echo ""
echo "=============================================="
echo "PACKAGING RESULTS"
echo "=============================================="

tar -czvf md_results_$(date +%Y%m%d).tar.gz     trajectories/*.dcd     analysis/*.json     logs/*.log

echo ""
echo "=============================================="
echo "PIPELINE COMPLETE!"
echo "=============================================="
echo ""
echo "Results packaged in: md_results_$(date +%Y%m%d).tar.gz"
echo ""
echo "Download with:"
echo "  scp user@gpu-instance:~/nod2/md_results_*.tar.gz ./"
echo ""
