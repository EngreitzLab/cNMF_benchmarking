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
#SBATCH --mem=128G
#SBATCH --gres=gpu:1

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu

START_TIME=$(date +%s)

# Print job and system information for debugging
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Working directory: $(pwd)"

# Configuration
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result"
RUN_NAME="030926_test_torch"
LOG_DIR="$OUT_DIR/$RUN_NAME"

# Create logs directory if it doesn't exist
mkdir -p "$LOG_DIR/logs"

# Environment information
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"

# Activate conda environment
echo "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate torch-cNMF

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"

# Start resource monitoring
MONITOR_LOG="$LOG_DIR/logs/resource_monitor_${SLURM_JOB_ID}.log"

monitor_resources() {
    while true; do
        echo "$(date '+%Y-%m-%d %H:%M:%S')" >> "$MONITOR_LOG"
        echo "=== Memory Usage ===" >> "$MONITOR_LOG"
        free -h >> "$MONITOR_LOG"
        echo "=== GPU Usage ===" >> "$MONITOR_LOG"
        nvidia-smi --query-gpu=timestamp,name,utilization.gpu,utilization.memory,memory.total,memory.used,memory.free,temperature.gpu --format=csv >> "$MONITOR_LOG" 2>/dev/null || echo "GPU monitoring not available" >> "$MONITOR_LOG"
        echo "---" >> "$MONITOR_LOG"
        sleep 30
    done
}

monitor_resources &
MONITOR_PID=$!

echo "Initial GPU status:"
nvidia-smi 2>/dev/null || echo "GPU monitoring not available"

# Run the torch-cNMF inference pipeline
echo "Running torch-cNMF inference pipeline..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Inference/torch-cNMF/Slurm_Version/torch-cNMF_inference_pipeline.py \
    --counts_fn "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Data/raw.h5ad" \
    --output_directory "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --algo "halsvar" \
    --mode "batch" \
    --init "random" \
    --tol 1e-7 \
    --batch_max_iter 1000 \
    --batch_hals_max_iter 1000 \
    --batch_hals_tol 0.005 \
    --online_chunk_size 200000 \
    --online_max_pass 1000 \
    --online_chunk_max_iter 1000 \
    --numiter 100 \
    --online_usage_tol 0.005 \
    --online_spectra_tol 0.005 \
    --species "human" \
    --use_gpu \
    --run_factorize \
    --run_refit \
    --run_complie_annotation \
    --sel_thresh 0.2 2.0 \
    --numhvgenes 17538 \
    --K 50

# Record end time and calculate duration
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))

# Stop resource monitoring
kill $MONITOR_PID 2>/dev/null

echo "========================================="
echo "EXECUTION SUMMARY"
echo "========================================="
echo "Total execution time: ${DURATION} seconds ($(($DURATION / 3600))h $(($DURATION % 3600 / 60))m $(($DURATION % 60))s)"
echo "Final GPU status:"
nvidia-smi 2>/dev/null || echo "GPU monitoring not available"

echo "Job completed at: $(date)"
