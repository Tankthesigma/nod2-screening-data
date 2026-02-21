#!/bin/bash
# Launch FEP windows on a single node with 8 GPUs
# Usage: ./launch_node.sh <start_task> <end_task>
# Tasks 0-19: wt_complex windows 0-19
# Tasks 20-39: mut_complex windows 0-19

set -e

START_TASK=${1:-0}
END_TASK=${2:-19}
BASE_DIR="/root/fep_pmx_natural"

echo "=============================================="
echo "FEP Node Launcher"
echo "Tasks: $START_TASK - $END_TASK"
echo "=============================================="

# Install dependencies if needed
if ! python -c "import openmm" 2>/dev/null; then
    echo "Installing OpenMM..."
    pip install openmm openmmtools numpy
fi

# Function to get system name and window from task ID
get_task_info() {
    local task_id=$1
    if [ $task_id -lt 20 ]; then
        echo "wt_complex $task_id"
    else
        echo "mut_complex $((task_id - 20))"
    fi
}

# Distribute tasks across 8 GPUs
declare -a PIDS=()
declare -a GPU_LOGS=()

for gpu_id in {0..7}; do
    GPU_LOGS[$gpu_id]="$BASE_DIR/gpu_${gpu_id}.log"
    > "${GPU_LOGS[$gpu_id]}"  # Clear log
done

# Assign tasks to GPUs round-robin
task_count=$((END_TASK - START_TASK + 1))
tasks_per_gpu=$(( (task_count + 7) / 8 ))  # Ceiling division

echo "Tasks per GPU: ~$tasks_per_gpu"

# Launch all GPU workers
for gpu_id in {0..7}; do
    (
        for ((i=0; i<tasks_per_gpu; i++)); do
            task_id=$((START_TASK + gpu_id + i * 8))
            if [ $task_id -gt $END_TASK ]; then
                break
            fi

            read sys_name window_idx <<< $(get_task_info $task_id)

            echo "[GPU $gpu_id] Running $sys_name window $window_idx"

            CUDA_VISIBLE_DEVICES=$gpu_id python $BASE_DIR/run_fep_gpu.py \
                $BASE_DIR $sys_name $window_idx $window_idx 0 \
                2>&1 | tee -a "${GPU_LOGS[$gpu_id]}"
        done
    ) &
    PIDS+=($!)
done

echo "Launched ${#PIDS[@]} GPU workers"
echo "Waiting for completion..."

# Wait for all workers
failed=0
for pid in "${PIDS[@]}"; do
    if ! wait $pid; then
        ((failed++))
    fi
done

echo "=============================================="
if [ $failed -eq 0 ]; then
    echo "[SUCCESS] All GPU workers completed!"
else
    echo "[WARNING] $failed GPU workers had failures"
fi
echo "=============================================="

exit $failed
