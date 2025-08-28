#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=sk-cNMF-10k           # Job name
#SBATCH --output=logs/sk-cNMF-10k_%j.out      # Output file (%j = job ID)
#SBATCH --error=logs/sk-cNMF-10k_%j.err       # Error file
#SBATCH --partition=engreitz           # partition name
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
source activate sk-cNMF

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"


# Run the Python script
echo "Running Python script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Inference/Slurm_Version/sk-cNMF_batch_inference_pipeline.py \
        --counts_fn "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Cell_data/10K_cell.h5ad" \
        --output_directory "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/sk-cNMF_evaluation" \
        --run_name "10k_10iter_sk_cd_frobenius"\
        --algo "cd"\


echo "Job completed at: $(date)"