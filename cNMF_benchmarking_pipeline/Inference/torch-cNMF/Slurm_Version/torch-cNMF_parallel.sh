#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=091625_100k_cells_10iter_torch_halsvar_batch_e7_par        # Job name
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/torch-cNMF_evaluation/091625_100k_cells_10iter_torch_halsvar_batch_e7_par/logs/%A_%a.out       # Output file (fixed double slash)
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/torch-cNMF_evaluation/091625_100k_cells_10iter_torch_halsvar_batch_e7_par/logs/%A_%a.err     # Error file (fixed double slash)
#SBATCH --partition=gpu                # partition name
#SBATCH --array=1-8                    # Run parallel jobs (array indices 1-#)
#SBATCH --time=15:00:00                # Time limit 
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=1              # CPUs per task (keeping 1 core as requested)
#SBATCH --mem=64G                      # Memory per node
#SBATCH --gres=gpu:1                   # Request 1 GPU (for future use)
##SBATCH -C GPU_SKU:H100_SXM5

# Optional: Request specific GPU type if available
##SBATCH --constraint="GPU_MEM:32GB" # GPU memory constraint

# Email notifications
#SBATCH --mail-type=BEGIN              # Send email when job starts
#SBATCH --mail-type=END                # Send email when job ends
#SBATCH --mail-type=FAIL               # Send email if job fails
#SBATCH --mail-user=ymo@stanford.edu   # the email address sent 

START_TIME=$(date +%s)

# Define K values array
K_VALUES=(30 50 60 80 100 200 250 300)

# Get K value for this array task
K=${K_VALUES[$((SLURM_ARRAY_TASK_ID-1))]}

# Configuration - Set your log directory here (fixed to match SLURM paths)
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/torch-cNMF_evaluation/"
RUN_NAME="091625_100k_cells_10iter_torch_halsvar_batch_e7_par"
LOG_DIR="$OUT_DIR/$RUN_NAME"

# Create logs directory if it doesn't exist
mkdir -p "$LOG_DIR/logs"

# Print job and system information for debugging
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "K value for this task: $K"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "Number of CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Log directory: $LOG_DIR"

# Environment information
echo "PATH: $PATH"
echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"

# Activate conda base environment
echo "Activating conda base environment..."
source activate torch-cNMF

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"

# Start resource monitoring
MONITOR_LOG="$LOG_DIR/logs/resource_monitor_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

# Function to monitor resources
monitor_resources() {
    while true; do
        echo "$(date '+%Y-%m-%d %H:%M:%S')" >> "$MONITOR_LOG"
        echo "=== Memory Usage ===" >> "$MONITOR_LOG"
        free -h >> "$MONITOR_LOG"
        echo "=== CPU Usage ===" >> "$MONITOR_LOG"
        top -bn1 | grep "Cpu(s)" >> "$MONITOR_LOG"
        echo "=== GPU Usage ===" >> "$MONITOR_LOG"
        nvidia-smi --query-gpu=timestamp,name,utilization.gpu,utilization.memory,memory.total,memory.used,memory.free,temperature.gpu --format=csv >> "$MONITOR_LOG" 2>/dev/null || echo "GPU monitoring not available" >> "$MONITOR_LOG"
        echo "=== Process Memory ===" >> "$MONITOR_LOG"
        ps -eo pid,ppid,cmd,%mem,%cpu --sort=-%mem | head -10 >> "$MONITOR_LOG"
        echo "---" >> "$MONITOR_LOG"
        sleep 30  # Monitor every 30 seconds
    done
}

# Start monitoring in background
monitor_resources &
MONITOR_PID=$!

# Record start time and initial memory
SCRIPT_START_TIME=$(date)
echo "Script start time: $SCRIPT_START_TIME"
echo "Initial memory usage:" 
free -h
echo "Initial GPU status:"
nvidia-smi 2>/dev/null || echo "GPU monitoring not available"

# Run the Python script (CPU-only for now)
echo "Running Python script with CPU-only (K=$K)..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Inference/torch-cNMF/Slurm_Version/torch-cNMF_inference_pipeline.py \
        --counts_fn "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Cell_data/100k_250genes.h5ad" \
        --output_directory "$OUT_DIR/$RUN_NAME" \
        --run_name "${RUN_NAME}_${K}" \
        --algo "halsvar" \
        --mode "batch" \
        --tol 1e-7 \
        --batch_max_iter 1000 \
        --batch_hals_max_iter 1000 \
        --batch_hals_tol 0.005 \
        --densify \
        --online_chunk_size 50000 \
        --online_max_pass 1000 \
        --online_chunk_max_iter 1000 \
        --numiter 10 \
        --online_usage_tol 0.005 \
        --online_spectra_tol 0.005 \
        --K $K \
        --use_gpu

# Record end time and calculate duration
END_TIME=$(date +%s)
SCRIPT_END_TIME=$(date)
DURATION=$((END_TIME - START_TIME))

# Stop resource monitoring
kill $MONITOR_PID 2>/dev/null

# Final resource summary
echo "========================================="
echo "EXECUTION SUMMARY"
echo "========================================="
echo "Array Task ID: $SLURM_ARRAY_TASK_ID"
echo "K value: $K"
echo "Script start time: $SCRIPT_START_TIME"
echo "Script end time: $SCRIPT_END_TIME"
echo "Total execution time: ${DURATION} seconds ($(($DURATION / 3600))h $(($DURATION % 3600 / 60))m $(($DURATION % 60))s)"
echo "Final memory usage:"
free -h
echo "Final GPU status:"
nvidia-smi 2>/dev/null || echo "GPU monitoring not available"

# Peak memory usage from SLURM
echo "SLURM reported peak memory usage: $(sacct -j $SLURM_JOB_ID --format=MaxRSS --noheader | head -1 | tr -d ' ')" 2>/dev/null || echo "SLURM memory stats not available"

echo "Resource monitoring log saved to: $MONITOR_LOG"
echo "Job completed at: $(date)"