#!/bin/bash
# Run solvent FEP simulations on Vast.ai
# Usage: bash run_solvent_fep.sh

set -e

echo "========================================"
echo "SOLVENT FEP SIMULATION RUNNER"
echo "========================================"

# Check for GPU
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader
echo ""

# Create workspace
mkdir -p /workspace/fep_data
cd /workspace/fep_data

# Check if data already extracted
if [ ! -f "fep_pmx/run_fep_gpu.py" ]; then
    echo "Extracting deployment package..."
    unzip -o /workspace/solvent_fep_deployment.zip
fi

# Install dependencies if needed
pip install -q openmm openmmtools numpy 2>/dev/null || true

# Check GPU count
N_GPUS=$(nvidia-smi --query-gpu=index --format=csv,noheader | wc -l)
echo "Detected $N_GPUS GPUs"

# Run simulations
# Febuxostat solvent windows 17-18 (on GPUs 0-1)
echo ""
echo "Running Febuxostat solvent windows 17-18..."
python fep_pmx/run_fep_gpu.py /workspace/fep_data/fep_pmx solvent 17 17 0 &
python fep_pmx/run_fep_gpu.py /workspace/fep_data/fep_pmx solvent 18 18 1 &

# CID_10592 solvent windows 0-19 (distributed across GPUs 2-7, then 0-1 after feb)
echo "Running CID_10592 solvent windows 0-19..."
python fep_pmx_natural/run_fep_gpu.py /workspace/fep_data/fep_pmx_natural solvent 0 2 2 &
python fep_pmx_natural/run_fep_gpu.py /workspace/fep_data/fep_pmx_natural solvent 3 5 3 &
python fep_pmx_natural/run_fep_gpu.py /workspace/fep_data/fep_pmx_natural solvent 6 8 4 &
python fep_pmx_natural/run_fep_gpu.py /workspace/fep_data/fep_pmx_natural solvent 9 11 5 &
python fep_pmx_natural/run_fep_gpu.py /workspace/fep_data/fep_pmx_natural solvent 12 14 6 &
python fep_pmx_natural/run_fep_gpu.py /workspace/fep_data/fep_pmx_natural solvent 15 17 7 &

# Wait for first batch
wait

# Run remaining windows on GPUs 0-1 (now free after feb)
python fep_pmx_natural/run_fep_gpu.py /workspace/fep_data/fep_pmx_natural solvent 18 19 0 &

wait

echo ""
echo "========================================"
echo "ALL SIMULATIONS COMPLETE"
echo "========================================"

# Check results
echo ""
echo "Checking results..."
for i in 17 18; do
    f="fep_pmx/solvent/window_$(printf '%02d' $i)/u_nk.npy"
    if [ -f "$f" ]; then
        echo "[OK] $f"
    else
        echo "[MISSING] $f"
    fi
done

for i in $(seq 0 19); do
    f="fep_pmx_natural/solvent/window_$(printf '%02d' $i)/u_nk.npy"
    if [ -f "$f" ]; then
        echo "[OK] $f"
    else
        echo "[MISSING] $f"
    fi
done

echo ""
echo "To download results:"
echo "  scp -r vastai:/workspace/fep_data/fep_pmx/solvent/window_17 ."
echo "  scp -r vastai:/workspace/fep_data/fep_pmx/solvent/window_18 ."
echo "  scp -r vastai:/workspace/fep_data/fep_pmx_natural/solvent/ ."
