#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=torch-cNMF-10k           # Job name
#SBATCH --output=logs/torch-cNMF-10k_%j.out      # Output file (%j = job ID)
#SBATCH --error=logs/torch-cNMF-10k_%j.err       # Error file
#SBATCH --partition=engretiz                # partition name
#SBATCH --time=05:00:00                # Time limit (5 minutes)
#SBATCH --nodes=1                      # Number of nodes
#SBATCH --ntasks=1                     # Number of tasks
#SBATCH --cpus-per-task=1              # CPUs per task
#SBATCH --mem=32G                       # Memory per node


# Email notifications
#SBATCH --mail-type=BEGIN              # Send email when job starts
#SBATCH --mail-type=END                # Send email when job ends
#SBATCH --mail-type=FAIL               # Send email if job fails
#SBATCH --mail-user=ymo@stanford.edu   # the email address sent 


# Print some job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"


# Create logs directory if it doesn't exist
mkdir -p logs

# Activate conda base environment
echo "Activating conda base environment..."
source activate torch-cNMF

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"


# Run the Inference script
echo "Running Inference script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Inference/Slurm_Version/sk-cNMF_batch_inference_pipeline.py 

echo "Inference completed at: $(date)"



# Run the Evaluation script
echo "Running Evaluation script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Evaluation/Slurm_Version/cNMF_evaluation_pipeline.py\
        --input_folder "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation/082525_100k_10iter_sk_mu_frobenius"

echo "Evaluation completed at: $(date)"



# Run the Plotting script
echo "Running Plotting script..."

echo "Plotting completed at: $(date)"