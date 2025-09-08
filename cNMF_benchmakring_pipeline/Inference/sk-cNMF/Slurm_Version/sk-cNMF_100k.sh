#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=090525_100k_10iter_1000batiter_sk_cd_frobenius          # Job name
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation/090525_100k_10iter_1000batiter_sk_cd_frobenius/logs/%A_%a.out      # Output file (%j = job ID)
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation/090525_100k_10iter_1000batiter_sk_cd_frobenius/logs/%A_%a.err       # Error file
#SBATCH --partition=engreitz           # partition name
#SBATCH --array=1-8                    # Run parallel jobs (array indices 1-#)
#SBATCH --time=40:00:00                # Time limit
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=1              # CPUs per task
#SBATCH --mem=32G                       # Memory per node

# Email notifications
#SBATCH --mail-type=BEGIN              # Send email when job starts
#SBATCH --mail-type=END                # Send email when job ends
#SBATCH --mail-type=FAIL               # Send email if job fails
#SBATCH --mail-user=ymo@stanford.edu   # the email address sent 

# Define the cNMF case
LOG_DIR="/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation/090525_100k_10iter_1000batiter_sk_cd_frobenius"

# Store start time
START_TIME=$(date +%s)


# Define K values array
K_VALUES=(30 50 60 80 100 200 250 300)

# Get K value for this array task
K=${K_VALUES[$((SLURM_ARRAY_TASK_ID-1))]}

# Print some job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "Number of CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Log directory: $LOG_DIR"


# Create logs directory if it doesn't exist
mkdir -p "$LOG_DIR/logs"

# Activate conda base environment
echo "Activating conda base environment..."
source activate sk-cNMF

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"


# Run the Python script
echo "Running Python script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Inference/sk-cNMF/Slurm_Version/sk-cNMF_batch_inference_pipeline.py \
        --counts_fn "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Cell_data/100k_250genes.h5ad" \
        --output_directory "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation/090525_100k_10iter_1000batiter_sk_cd_frobenius" \
        --run_name "090525_100k_10iter_1000batiter_sk_cd_frobenius_${K}" \
        --algo "cd" \
        --K $K \
        --max_NMF_iter 1000



# Calculate and print elapsed time at the end
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
