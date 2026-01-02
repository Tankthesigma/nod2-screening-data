#!/bin/bash
#===============================================================================
# NOD2-SCOUT: 9x RTX 5090 PARALLEL MD SIMULATION
#===============================================================================
#
# FIXED VERSION with:
# - OpenFF ligand parameterization
# - Proper NVT→NPT equilibration
# - Random seeds for reproducibility
# - Early crash detection
#
# Launches 14 simulations across 9 GPUs:
#   GPU 0: febuxostat_rep1
#   GPU 1: febuxostat_rep2
#   GPU 2: febuxostat_rep3
#   GPU 3: ursodiol_rep1
#   GPU 4: ursodiol_rep2
#   GPU 5: ursodiol_rep3
#   GPU 6: budesonide_rep1
#   GPU 7: budesonide_rep2
#   GPU 8: budesonide_rep3 + natural (3 reps) + apo + decoy (sequential)
#
# HMR enabled: 4fs timestep → 800-1100 ns/day per GPU
# Expected time: ~3 hours (GPU 8 takes longest)
# Expected cost: ~$12 on 9x RTX 5090 @ $4.20/hr
#
#===============================================================================

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR/.."

echo "============================================================"
echo "NOD2-SCOUT: 9x RTX 5090 PARALLEL MD SIMULATION"
echo "============================================================"
echo ""
echo "Started: $(date)"
echo "Working dir: $(pwd)"
echo ""

# Set CUDA_VISIBLE_DEVICES explicitly for predictable GPU mapping
export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7,8

# Check CUDA
echo "Checking GPUs..."
nvidia-smi --query-gpu=index,name,memory.total --format=csv
echo ""

# Verify we have 9 GPUs
GPU_COUNT=$(nvidia-smi --query-gpu=index --format=csv,noheader | wc -l)
if [ "$GPU_COUNT" -lt 9 ]; then
    echo "WARNING: Only $GPU_COUNT GPUs detected. Need 9 for full parallel."
    echo "Continuing anyway..."
fi

# Create output directories
mkdir -p trajectories logs checkpoints

# Function to run simulation with crash detection
run_gpu() {
    local gpu_script=$1
    local log_file=$2
    echo "Launching: $gpu_script → $log_file"
    nohup python3 parallel_9gpu/$gpu_script > logs/$log_file 2>&1 &
    local pid=$!
    echo $pid
}

# Function to check if process is still alive
check_alive() {
    local pid=$1
    if kill -0 $pid 2>/dev/null; then
        return 0
    else
        return 1
    fi
}

echo "============================================================"
echo "LAUNCHING 9 PARALLEL SIMULATIONS"
echo "============================================================"
echo ""

# Launch all GPUs
PIDS=()
SCRIPTS=()

# GPU 0-7: Single simulations
PIDS+=($(run_gpu "gpu0_febuxostat_rep1.py" "gpu0_febuxostat_rep1.log"))
SCRIPTS+=("gpu0_febuxostat_rep1")
PIDS+=($(run_gpu "gpu1_febuxostat_rep2.py" "gpu1_febuxostat_rep2.log"))
SCRIPTS+=("gpu1_febuxostat_rep2")
PIDS+=($(run_gpu "gpu2_febuxostat_rep3.py" "gpu2_febuxostat_rep3.log"))
SCRIPTS+=("gpu2_febuxostat_rep3")
PIDS+=($(run_gpu "gpu3_ursodiol_rep1.py" "gpu3_ursodiol_rep1.log"))
SCRIPTS+=("gpu3_ursodiol_rep1")
PIDS+=($(run_gpu "gpu4_ursodiol_rep2.py" "gpu4_ursodiol_rep2.log"))
SCRIPTS+=("gpu4_ursodiol_rep2")
PIDS+=($(run_gpu "gpu5_ursodiol_rep3.py" "gpu5_ursodiol_rep3.log"))
SCRIPTS+=("gpu5_ursodiol_rep3")
PIDS+=($(run_gpu "gpu6_budesonide_rep1.py" "gpu6_budesonide_rep1.log"))
SCRIPTS+=("gpu6_budesonide_rep1")
PIDS+=($(run_gpu "gpu7_budesonide_rep2.py" "gpu7_budesonide_rep2.log"))
SCRIPTS+=("gpu7_budesonide_rep2")

# GPU 8: Sequential (6 simulations)
PIDS+=($(run_gpu "gpu8_remaining.py" "gpu8_remaining.log"))
SCRIPTS+=("gpu8_remaining")

echo ""
echo "============================================================"
echo "ALL SIMULATIONS LAUNCHED"
echo "============================================================"
echo ""
echo "PIDs: ${PIDS[@]}"
echo ""

# Early crash detection (wait 5 seconds, then check if all processes still alive)
echo "Checking for early crashes (5 second delay)..."
sleep 5

CRASHED=0
for i in "${!PIDS[@]}"; do
    if ! check_alive ${PIDS[$i]}; then
        echo "  CRASHED: ${SCRIPTS[$i]} (PID ${PIDS[$i]})"
        echo "  Check log: logs/${SCRIPTS[$i]}.log"
        CRASHED=$((CRASHED + 1))
    else
        echo "  RUNNING: ${SCRIPTS[$i]} (PID ${PIDS[$i]})"
    fi
done

if [ $CRASHED -gt 0 ]; then
    echo ""
    echo "WARNING: $CRASHED simulation(s) crashed within 5 seconds!"
    echo "Check the log files for errors before continuing."
    echo ""
fi

echo ""
echo "Monitor progress:"
echo "  tail -f logs/gpu0_febuxostat_rep1.log"
echo "  tail -f logs/gpu8_remaining.log"
echo ""
echo "Check GPU usage:"
echo "  watch -n 5 nvidia-smi"
echo ""
echo "Wait for completion:"
echo "  wait ${PIDS[@]}"
echo ""

# Option to wait
if [ "$1" == "--wait" ]; then
    echo "Waiting for all simulations to complete..."
    wait ${PIDS[@]}
    echo ""
    echo "============================================================"
    echo "ALL SIMULATIONS COMPLETE!"
    echo "============================================================"
    echo "Finished: $(date)"

    # Move outputs
    echo ""
    echo "Moving output files..."
    mv *.dcd trajectories/ 2>/dev/null || true
    mv *.chk checkpoints/ 2>/dev/null || true
    mv *_final.xml checkpoints/ 2>/dev/null || true

    echo ""
    echo "Output locations:"
    echo "  Trajectories: trajectories/"
    echo "  Checkpoints: checkpoints/"
    echo "  Logs: logs/"

    # Package results
    echo ""
    echo "Packaging results..."
    tar -czvf md_results_9gpu_$(date +%Y%m%d_%H%M).tar.gz \
        trajectories/*.dcd \
        checkpoints/*.chk \
        logs/*.log

    echo ""
    echo "Download with:"
    echo "  scp user@gpu-instance:$(pwd)/md_results_9gpu_*.tar.gz ./"
fi

echo ""
echo "Script finished at: $(date)"
