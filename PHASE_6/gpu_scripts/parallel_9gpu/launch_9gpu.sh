#!/bin/bash
#===============================================================================
# NOD2-SCOUT: 9x RTX 5090 PARALLEL MD SIMULATION
#===============================================================================
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
# Expected time: ~1.5-2 hours (GPU 8 takes longest: ~3 hours)
# Expected cost: ~$8-12 on 9x RTX 5090 @ $4.20/hr
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

# Check CUDA
echo "Checking GPUs..."
nvidia-smi --query-gpu=index,name,memory.total --format=csv
echo ""

# Create output directories
mkdir -p trajectories logs checkpoints

# Function to run simulation
run_gpu() {
    local gpu_script=$1
    local log_file=$2
    echo "Launching: $gpu_script → $log_file"
    nohup python3 parallel_9gpu/$gpu_script > logs/$log_file 2>&1 &
    echo $!
}

echo "============================================================"
echo "LAUNCHING 9 PARALLEL SIMULATIONS"
echo "============================================================"
echo ""

# Launch all GPUs
PIDS=()

# GPU 0-7: Single simulations
PIDS+=($(run_gpu "gpu0_febuxostat_rep1.py" "gpu0_febuxostat_rep1.log"))
PIDS+=($(run_gpu "gpu1_febuxostat_rep2.py" "gpu1_febuxostat_rep2.log"))
PIDS+=($(run_gpu "gpu2_febuxostat_rep3.py" "gpu2_febuxostat_rep3.log"))
PIDS+=($(run_gpu "gpu3_ursodiol_rep1.py" "gpu3_ursodiol_rep1.log"))
PIDS+=($(run_gpu "gpu4_ursodiol_rep2.py" "gpu4_ursodiol_rep2.log"))
PIDS+=($(run_gpu "gpu5_ursodiol_rep3.py" "gpu5_ursodiol_rep3.log"))
PIDS+=($(run_gpu "gpu6_budesonide_rep1.py" "gpu6_budesonide_rep1.log"))
PIDS+=($(run_gpu "gpu7_budesonide_rep2.py" "gpu7_budesonide_rep2.log"))

# GPU 8: Sequential (6 simulations)
PIDS+=($(run_gpu "gpu8_remaining.py" "gpu8_remaining.log"))

echo ""
echo "============================================================"
echo "ALL SIMULATIONS LAUNCHED"
echo "============================================================"
echo ""
echo "PIDs: ${PIDS[@]}"
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
