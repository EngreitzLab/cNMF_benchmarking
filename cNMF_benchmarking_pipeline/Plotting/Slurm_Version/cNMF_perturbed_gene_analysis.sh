#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=cNMF100_plot           # Job name
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/cNMF_100_07102024/consensus_NMF/Plot/Perturb_gene/logs/%j.out      # Output file (%j = job ID)
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/cNMF_100_07102024/consensus_NMF/Plot/Perturb_gene/logs/%j.err       # Error file
#SBATCH --partition=engreitz            # partition name
#SBATCH --time=15:00:00                  # Time limit 
#SBATCH --nodes=1                       # Number of nodes
#SBATCH --ntasks=1                      # Number of tasks
#SBATCH --cpus-per-task=20               # CPUs per task
#SBATCH --mem=240G                       # Memory per node

# Email notifications
#SBATCH --mail-type=BEGIN              # Send email when job starts
#SBATCH --mail-type=END                # Send email when job ends
#SBATCH --mail-type=FAIL               # Send email if job fails
#SBATCH --mail-user=ymo@stanford.edu   # the email address sent 



# Define the cNMF case
LOG_DIR="/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Results/cNMF_100_07102024/consensus_NMF/Plot/Perturb_gene"

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
mkdir -p "$LOG_DIR/logs"

# Activate conda base environment
echo "Activating conda base environment..."
source activate NMF_Benchmarking

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"


# Run the Python script
echo "Running Python script..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Plotting/Slurm_Version/cNMF_perturbed_gene_analysis.py\
        --mdata_path "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Cell_data/Tri_orginal_data.h5mu" \
        --perturb_path "/oak/stanford/groups/engreitz/Users/ymo/NMF_re-inplementing/Script/Revant_code/orginal_code/shared/250110_ipsc_ec_dashboard_setup/cNMF_100/cNMF_100_gene_sample"\
        --top_program 5 \
        --p_value 0.05 \
        --pdf_save_path "$LOG_DIR" \
        --PDF \
        --file_to_dictionary "/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Resources/weissman_guides_with_coordinates.tsv" \
        --n_processes 70 \
        --samples D0 sample_D1 sample_D2 sample_D3 \
        --square_plots \
        --figsize 35 35




# Calculate and print elapsed time at the end
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
