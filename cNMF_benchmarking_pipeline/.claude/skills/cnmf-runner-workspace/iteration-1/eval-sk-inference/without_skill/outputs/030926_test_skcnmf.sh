#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=030926_test_skcnmf
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030926_test_skcnmf/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030926_test_skcnmf/logs/%j.err
#SBATCH --partition=engreitz,owners,bigmem
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=256G

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu

# Define paths
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result"
RUN_NAME="030926_test_skcnmf"
COUNTS_FN="/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Data/raw.h5ad"
PIPELINE_DIR="/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline"
CONDA_BASE="/oak/stanford/groups/engreitz/Users/ymo/miniforge3"
LOG_DIR="$OUT_DIR/$RUN_NAME/logs"

# Store start time
START_TIME=$(date +%s)

# Print job information
echo "============================================================"
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "Number of CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Run name: $RUN_NAME"
echo "Output directory: $OUT_DIR"
echo "Input file: $COUNTS_FN"
echo "K values: 30 50 80 100"
echo "Iterations: 10"
echo "Species: human"
echo "============================================================"

# Create logs directory if it doesn't exist
mkdir -p "$LOG_DIR"

# Activate conda environment
echo "Activating conda environment..."
source "$CONDA_BASE/etc/profile.d/conda.sh"
eval "$(conda shell.bash hook)"
conda activate sk-cNMF

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"

# Run the sk-cNMF batch inference pipeline
echo "Running sk-cNMF batch inference pipeline..."
python3 "$PIPELINE_DIR/Inference/sk-cNMF/Slurm_Version/sk-cNMF_batch_inference_pipeline.py" \
    --counts_fn "$COUNTS_FN" \
    --output_directory "$OUT_DIR" \
    --run_name "$RUN_NAME" \
    --K 30 50 80 100 \
    --numiter 10 \
    --numhvgenes 2000 \
    --seed 14 \
    --init "random" \
    --algo "mu" \
    --loss "frobenius" \
    --max_NMF_iter 1000 \
    --tol 1e-4 \
    --sel_thresh 0.2 2.0 \
    --species "human" \
    --run_factorize \
    --run_refit \
    --run_complie_annotation

# Calculate and print elapsed time
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "============================================================"
echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
echo "============================================================"
