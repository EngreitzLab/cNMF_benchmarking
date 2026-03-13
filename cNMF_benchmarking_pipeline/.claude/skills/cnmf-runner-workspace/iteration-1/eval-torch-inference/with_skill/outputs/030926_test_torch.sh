#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=030926_test_torch
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030926_test_torch/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030926_test_torch/logs/%j.err
#SBATCH --partition=engreitz
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=256G
#SBATCH --gres=gpu:1
#SBATCH -C GPU_SKU:V100S_PCIE

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu

# Configuration
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result"
RUN_NAME="030926_test_torch"
LOG_DIR="$OUT_DIR/$RUN_NAME/logs"

# Store start time
START_TIME=$(date +%s)

# Print job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "CPUs: $SLURM_CPUS_PER_TASK"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Log directory: $LOG_DIR"

# Create log directory
mkdir -p "$LOG_DIR"

# Start resource monitoring
MONITOR_LOG="$LOG_DIR/resource_monitor_${SLURM_JOB_ID}.log"
monitor_resources() {
    while true; do
        echo "$(date '+%Y-%m-%d %H:%M:%S')" >> "$MONITOR_LOG"
        echo "=== GPU Usage ===" >> "$MONITOR_LOG"
        nvidia-smi --query-gpu=timestamp,name,utilization.gpu,utilization.memory,memory.total,memory.used,memory.free,temperature.gpu --format=csv >> "$MONITOR_LOG" 2>/dev/null
        echo "=== Memory ===" >> "$MONITOR_LOG"
        free -h >> "$MONITOR_LOG"
        echo "---" >> "$MONITOR_LOG"
        sleep 30
    done
}
monitor_resources &
MONITOR_PID=$!

# Activate conda environment
echo "Activating conda environment: torch-cNMF"
eval "$(conda shell.bash hook)"
conda activate torch-cNMF

echo "Active env: $CONDA_DEFAULT_ENV"
echo "Python: $(python --version)"

# Run pipeline
echo "Running: inference-torch"
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Inference/torch-cNMF/Slurm_Version/torch-cNMF_inference_pipeline.py \
        --counts_fn "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Data/raw.h5ad" \
        --output_directory "$OUT_DIR" \
        --run_name "$RUN_NAME" \
        --algo "halsvar" \
        --mode "batch" \
        --init "random" \
        --tol 1e-4 \
        --batch_max_iter 500 \
        --batch_hals_max_iter 200 \
        --batch_hals_tol 0.05 \
        --numiter 100 \
        --species "human" \
        --use_gpu \
        --run_factorize \
        --run_refit \
        --run_complie_annotation \
        --sel_thresh 0.2 2.0 \
        --numhvgenes 17538 \
        --K 50

# Stop resource monitoring
kill $MONITOR_PID 2>/dev/null
echo "Final GPU status:"
nvidia-smi 2>/dev/null || echo "GPU not available"

# Elapsed time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED} seconds)"
