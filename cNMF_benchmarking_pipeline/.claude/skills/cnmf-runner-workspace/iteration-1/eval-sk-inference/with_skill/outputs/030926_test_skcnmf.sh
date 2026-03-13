#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=030926_test_skcnmf
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030926_test_skcnmf/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030926_test_skcnmf/logs/%j.err
#SBATCH --partition=engreitz,owners
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=256G

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu

# Configuration
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result"
RUN_NAME="030926_test_skcnmf"
LOG_DIR="/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030926_test_skcnmf/logs"

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

# Activate conda environment
echo "Activating conda environment: sk-cNMF"
eval "$(conda shell.bash hook)"
conda activate sk-cNMF

echo "Active env: $CONDA_DEFAULT_ENV"
echo "Python: $(python --version)"

# Run pipeline
echo "Running: inference-sk"
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Inference/sk-cNMF/Slurm_Version/sk-cNMF_batch_inference_pipeline.py \
        --counts_fn /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Data/raw.h5ad \
        --output_directory /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result \
        --run_name 030926_test_skcnmf \
        --species human \
        --K 30 50 80 100 \
        --numiter 10 \
        --run_factorize \
        --run_refit \
        --run_complie_annotation

# Elapsed time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED} seconds)"
