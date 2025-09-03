#!/bin/bash

# Define the cNMF case
LOG_DIR="/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation/082525_100k_10iter_sk_mu_frobenius"


# SLURM job configuration
#SBATCH --job-name=cNMF_evaluation                                # Job name
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation/082525_100k_10iter_sk_mu_frobenius/Eval/logs/cNMF_evaluation_%j.out      # Output file (%j = job ID)
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation/082525_100k_10iter_sk_mu_frobenius/Eval/logs/cNMF_evaluation_%j.err       # Error file
#SBATCH --partition=service              # partition name
#SBATCH --time=25:00:00                 # Time limit (5 minutes)
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=10              # CPUs per task
#SBATCH --mem=32G                       # Memory per node

# Email notifications
#SBATCH --mail-type=BEGIN              # Send email when job starts
#SBATCH --mail-type=END                # Send email when job ends
#SBATCH --mail-type=FAIL               # Send email if job fails
#SBATCH --mail-user=ymo@stanford.edu   # the email address sent 


# Store start time
START_TIME=$(date +%s)

# Print some job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "Number of CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Log directory: $LOG_DIR"


# Create logs directory if it doesn't exist
mkdir -p "$LOG_DIR/Eval/logs"

# Activate conda base environment
echo "Activating conda base environment..."
source activate NMF_Benchmarking

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"


# Run the Python script
echo "Running Python script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Evaluation/Slurm_Version/cNMF_evaluation_pipeline.py\
        --input_folder "$LOG_DIR"\
        --Perform_categorical False\
        --Perform_perturbation True\
        --Perform_geneset False\
        --Perform_trait False\
        --Perform_explained_variance False\
        --Perform_motif True


# Calculate and print elapsed time at the end
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
