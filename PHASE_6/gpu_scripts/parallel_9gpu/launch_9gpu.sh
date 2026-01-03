#!/bin/bash
#===============================================================================
# NOD2-SCOUT: 9x RTX 5090 PARALLEL MD SIMULATION
#===============================================================================
#
# FULLY FIXED VERSION:
# - Robust error handling (no set -e killing everything)
# - Early crash detection
# - Proper output directory handling
# - CUDA 12.x support for RTX 5090
# - nullglob for safe file operations
#
# GPU Assignment:
#   GPU 0: febuxostat_rep1     GPU 3: ursodiol_rep1      GPU 6: budesonide_rep1
#   GPU 1: febuxostat_rep2     GPU 4: ursodiol_rep2      GPU 7: budesonide_rep2
#   GPU 2: febuxostat_rep3     GPU 5: ursodiol_rep3      GPU 8: remaining (6 sims)
#
# Expected: ~3 hours total, ~$12 on 9x RTX 5090
#===============================================================================

# Enable nullglob for safe file operations
shopt -s nullglob

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR/.."

echo "============================================================"
echo "NOD2-SCOUT: 9x RTX 5090 PARALLEL MD SIMULATION"
echo "============================================================"
echo ""
echo "Started: $(date)"
echo "Working dir: $(pwd)"
echo ""

#===============================================================================
# PREFLIGHT CHECKS
#===============================================================================
echo "============================================================"
echo "PREFLIGHT CHECKS"
echo "============================================================"

# Check Python
if ! command -v python3 &> /dev/null; then
    echo "ERROR: python3 not found!"
    exit 1
fi
echo "  Python: $(python3 --version)"

# Check OpenMM
if ! python3 -c "import openmm" 2>/dev/null; then
    echo "ERROR: OpenMM not installed!"
    echo "  Run: conda install -c conda-forge openmm"
    exit 1
fi
echo "  OpenMM: OK"

# Check OpenFF
if ! python3 -c "from openff.toolkit import Molecule" 2>/dev/null; then
    echo "WARNING: OpenFF not installed (will use GAFF fallback)"
else
    echo "  OpenFF: OK"
fi

# Check CUDA
echo ""
echo "Checking GPUs..."
if command -v nvidia-smi &> /dev/null; then
    nvidia-smi --query-gpu=index,name,memory.total --format=csv
    GPU_COUNT=$(nvidia-smi --query-gpu=index --format=csv,noheader | wc -l)
    echo ""
    echo "  GPUs detected: $GPU_COUNT"
    if [ "$GPU_COUNT" -lt 9 ]; then
        echo "  WARNING: Only $GPU_COUNT GPUs. Need 9 for full parallel."
    fi
else
    echo "ERROR: nvidia-smi not found!"
    exit 1
fi

# Set CUDA environment
export CUDA_VISIBLE_DEVICES=0,1,2,3,4,5,6,7,8

#===============================================================================
# CREATE OUTPUT DIRECTORIES
#===============================================================================
echo ""
echo "Creating output directories..."
mkdir -p trajectories logs checkpoints
echo "  trajectories/ logs/ checkpoints/ - OK"

#===============================================================================
# LAUNCH SIMULATIONS
#===============================================================================
echo ""
echo "============================================================"
echo "LAUNCHING 9 PARALLEL SIMULATIONS"
echo "============================================================"
echo ""

PIDS=()
SCRIPTS=()

launch_sim() {
    local script=$1
    local logfile=$2
    echo "  Launching: $script"
    nohup python3 parallel_9gpu/$script > logs/$logfile 2>&1 &
    local pid=$!
    PIDS+=($pid)
    SCRIPTS+=("$script")
}

# Launch GPU 0-7 (single simulations)
launch_sim "gpu0_febuxostat_rep1.py" "gpu0_febuxostat_rep1.log"
launch_sim "gpu1_febuxostat_rep2.py" "gpu1_febuxostat_rep2.log"
launch_sim "gpu2_febuxostat_rep3.py" "gpu2_febuxostat_rep3.log"
launch_sim "gpu3_ursodiol_rep1.py" "gpu3_ursodiol_rep1.log"
launch_sim "gpu4_ursodiol_rep2.py" "gpu4_ursodiol_rep2.log"
launch_sim "gpu5_ursodiol_rep3.py" "gpu5_ursodiol_rep3.log"
launch_sim "gpu6_budesonide_rep1.py" "gpu6_budesonide_rep1.log"
launch_sim "gpu7_budesonide_rep2.py" "gpu7_budesonide_rep2.log"

# Launch GPU 8 (sequential - 6 simulations)
launch_sim "gpu8_remaining.py" "gpu8_remaining.log"

echo ""
echo "All simulations launched!"
echo "PIDs: ${PIDS[@]}"
echo ""

#===============================================================================
# EARLY CRASH DETECTION
#===============================================================================
echo "Checking for early crashes (10 second delay)..."
sleep 10

CRASHED=0
RUNNING=0
for i in "${!PIDS[@]}"; do
    if kill -0 ${PIDS[$i]} 2>/dev/null; then
        echo "  RUNNING: ${SCRIPTS[$i]} (PID ${PIDS[$i]})"
        ((RUNNING++))
    else
        echo "  CRASHED: ${SCRIPTS[$i]} (PID ${PIDS[$i]})"
        echo "           Check: logs/${SCRIPTS[$i]%.py}.log"
        ((CRASHED++))
    fi
done

echo ""
if [ $CRASHED -gt 0 ]; then
    echo "WARNING: $CRASHED simulation(s) crashed within 10 seconds!"
    echo "Check logs for errors. Common issues:"
    echo "  - Missing input files"
    echo "  - Ligand parameterization failure"
    echo "  - CUDA/GPU errors"
    echo ""
fi

echo "Running: $RUNNING / ${#PIDS[@]} simulations"
echo ""

#===============================================================================
# MONITORING COMMANDS
#===============================================================================
echo "============================================================"
echo "MONITORING"
echo "============================================================"
echo ""
echo "Monitor progress:"
echo "  tail -f logs/gpu0_febuxostat_rep1.log"
echo "  tail -f logs/gpu8_remaining.log"
echo ""
echo "Check GPU usage:"
echo "  watch -n 5 nvidia-smi"
echo ""
echo "Wait for all to complete:"
echo "  wait ${PIDS[@]}"
echo ""

#===============================================================================
# WAIT MODE
#===============================================================================
if [ "$1" == "--wait" ]; then
    echo "============================================================"
    echo "WAITING FOR COMPLETION"
    echo "============================================================"
    echo ""

    # Wait for all processes
    FINAL_RESULTS=()
    for i in "${!PIDS[@]}"; do
        wait ${PIDS[$i]}
        EXIT_CODE=$?
        if [ $EXIT_CODE -eq 0 ]; then
            FINAL_RESULTS+=("${SCRIPTS[$i]}: SUCCESS")
        else
            FINAL_RESULTS+=("${SCRIPTS[$i]}: FAILED (exit $EXIT_CODE)")
        fi
    done

    echo ""
    echo "============================================================"
    echo "ALL SIMULATIONS COMPLETE"
    echo "============================================================"
    echo "Finished: $(date)"
    echo ""

    # Print results
    echo "Results:"
    for result in "${FINAL_RESULTS[@]}"; do
        echo "  $result"
    done
    echo ""

    # Safely move output files (nullglob handles missing files)
    echo "Collecting output files..."

    # Count files in each directory
    DCD_COUNT=$(ls trajectories/*.dcd 2>/dev/null | wc -l)
    CHK_COUNT=$(ls checkpoints/*.chk 2>/dev/null | wc -l)
    LOG_COUNT=$(ls logs/*.log 2>/dev/null | wc -l)

    echo "  Trajectories: $DCD_COUNT .dcd files"
    echo "  Checkpoints: $CHK_COUNT .chk files"
    echo "  Logs: $LOG_COUNT .log files"

    # Package results (only if files exist)
    if [ "$DCD_COUNT" -gt 0 ]; then
        ARCHIVE="md_results_$(date +%Y%m%d_%H%M).tar.gz"
        echo ""
        echo "Packaging results to $ARCHIVE..."

        # Build tar command with existing files only
        TAR_FILES=""
        [ "$DCD_COUNT" -gt 0 ] && TAR_FILES="$TAR_FILES trajectories/*.dcd"
        [ "$CHK_COUNT" -gt 0 ] && TAR_FILES="$TAR_FILES checkpoints/*.chk"
        [ "$LOG_COUNT" -gt 0 ] && TAR_FILES="$TAR_FILES logs/*.log"

        tar -czvf "$ARCHIVE" $TAR_FILES 2>/dev/null || echo "WARNING: Some files missing from archive"

        echo ""
        echo "Download with:"
        echo "  scp user@gpu-instance:$(pwd)/$ARCHIVE ./"
    else
        echo ""
        echo "WARNING: No trajectory files found!"
        echo "Check logs for errors."
    fi
fi

echo ""
echo "============================================================"
echo "Script finished: $(date)"
echo "============================================================"
