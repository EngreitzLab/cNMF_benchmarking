#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=sk-cNMF-100k           # Job name
#SBATCH --output=logs/sk-cNMF-100k_%j.out      # Output file (%j = job ID)
#SBATCH --error=logs/sk-cNMF-100k_%j.err       # Error file
#SBATCH --partition=engreitz           # partition name
#SBATCH --array=1-3                    # Run parallel jobs (array indices 1-#)
#SBATCH --time=35:00:00                # Time limit (5 minutes)
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=1              # CPUs per task
#SBATCH --mem=32G                       # Memory per node

# Email notifications
#SBATCH --mail-type=BEGIN              # Send email when job starts
#SBATCH --mail-type=END                # Send email when job ends
#SBATCH --mail-type=FAIL               # Send email if job fails
#SBATCH --mail-user=ymo@stanford.edu   # the email address sent 

# Define K values array
K_VALUES=(200 250 300)

# Get K value for this array task
K=${K_VALUES[$((SLURM_ARRAY_TASK_ID-1))]}

# Print some job information
echo "Job started at: $(date)"
echo "Main Job ID: $SLURM_JOB_ID"
echo "Array Job ID: $SLURM_ARRAY_JOB_ID"
echo "Running with K=$K"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"


# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda base environment
echo "Activating conda base environment..."
source activate sk-cNMF

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"


# Run the Python script
echo "Running Python script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Script/Python/sk-cNMF_batch_inference_pipeline.py \
        --counts_fn "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Cell_data/100k_250genes.h5ad" \
        --output_directory "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation" \
        --run_name "082525_100k_10iter_sk_cd_frobenius_test_${K}"\
        --algo "cd"\
         --K $K


echo "Job completed at: $(date)"