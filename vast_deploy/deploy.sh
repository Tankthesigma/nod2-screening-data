#!/bin/bash
# Vast.ai FEP Deployment Script
# Runs assigned windows on a single GPU
# Dependencies must be installed via setup_once.sh first!

set -e

echo "=== FEP Worker ==="

# GPU assignment (8 GPUs, indices 0-7)
GPU_ID=$1
if [ -z "$GPU_ID" ]; then
    echo "ERROR: Usage: ./deploy.sh <GPU_ID 0-7>"
    exit 1
fi

# Check if numeric
if ! [[ "$GPU_ID" =~ ^[0-9]+$ ]]; then
    echo "ERROR: GPU_ID must be a number, got: $GPU_ID"
    exit 1
fi

if [ "$GPU_ID" -lt 0 ] || [ "$GPU_ID" -gt 7 ]; then
    echo "ERROR: GPU_ID must be 0-7, got: $GPU_ID"
    exit 1
fi

echo "Running on GPU $GPU_ID"
export CUDA_VISIBLE_DEVICES=$GPU_ID

# Base directory (relative to script location)
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_DIR="${SCRIPT_DIR}/data/fep_pmx"

if [ ! -d "$DATA_DIR" ]; then
    echo "ERROR: Data directory not found: $DATA_DIR"
    exit 1
fi

# Window assignments per GPU (42 windows across 8 GPUs)
case $GPU_ID in
    0) WINDOWS="wt_complex:18 wt_complex:19 mut_complex:00 mut_complex:01 mut_complex:02" ;;
    1) WINDOWS="mut_complex:03 mut_complex:04 mut_complex:05 mut_complex:06 mut_complex:07" ;;
    2) WINDOWS="mut_complex:08 mut_complex:09 mut_complex:10 mut_complex:11 mut_complex:12" ;;
    3) WINDOWS="mut_complex:13 mut_complex:14 mut_complex:15 mut_complex:16 mut_complex:17" ;;
    4) WINDOWS="mut_complex:18 mut_complex:19 solvent:00 solvent:01 solvent:02" ;;
    5) WINDOWS="solvent:03 solvent:04 solvent:05 solvent:06 solvent:07" ;;
    6) WINDOWS="solvent:08 solvent:09 solvent:10 solvent:11 solvent:12" ;;
    7) WINDOWS="solvent:13 solvent:14 solvent:15 solvent:16 solvent:17 solvent:18 solvent:19" ;;
    *)
        echo "ERROR: Invalid GPU_ID: $GPU_ID"
        exit 1
        ;;
esac

echo "Windows to run: $WINDOWS"
COMPLETED=0
FAILED=0

for WIN in $WINDOWS; do
    SYSTEM=$(echo $WIN | cut -d: -f1)
    NUM=$(echo $WIN | cut -d: -f2)
    WINDOW_DIR="${DATA_DIR}/${SYSTEM}/window_${NUM}"

    echo ""
    echo "=========================================="
    echo "[GPU $GPU_ID] Running ${SYSTEM}/window_${NUM}"
    echo "=========================================="

    if [ ! -d "$WINDOW_DIR" ]; then
        echo "ERROR: Window directory not found: $WINDOW_DIR"
        FAILED=$((FAILED + 1))
        continue
    fi

    if ! cd "$WINDOW_DIR"; then
        echo "ERROR: Cannot cd to $WINDOW_DIR"
        FAILED=$((FAILED + 1))
        continue
    fi

    if python run_window.py; then
        echo "[DONE] ${SYSTEM}/window_${NUM}"
        COMPLETED=$((COMPLETED + 1))
    else
        echo "[FAILED] ${SYSTEM}/window_${NUM}"
        FAILED=$((FAILED + 1))
    fi
done

echo ""
echo "=========================================="
echo "GPU $GPU_ID COMPLETE"
echo "Completed: $COMPLETED, Failed: $FAILED"
echo "=========================================="
