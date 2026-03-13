#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=030526_100k_cells
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells/Eval/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells/Eval/logs/%j.err
#SBATCH --partition=engreitz,owners
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu

# Configuration
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result"
RUN_NAME="030526_100k_cells"
LOG_DIR="/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells/Eval/logs"

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
echo "Activating conda environment: NMF_Benchmarking"
eval "$(conda shell.bash hook)"
conda activate NMF_Benchmarking

echo "Active env: $CONDA_DEFAULT_ENV"
echo "Python: $(python --version)"

# Run pipeline
echo "Running: evaluation"
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Slurm_Version/cNMF_evaluation_pipeline.py \
        --out_dir /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result \
        --run_name 030526_100k_cells \
        --K 50 \
        --sel_thresh 0.2 2.0 \
        --Perform_categorical \
        --Perform_perturbation \
        --Perform_geneset \
        --gwas_data_path /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz \
        --guide_annotation_path /oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Data/guide/guide_metadata_v43.tsv

# Elapsed time
END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
HOURS=$((ELAPSED / 3600))
MINUTES=$(((ELAPSED % 3600) / 60))
SECONDS=$((ELAPSED % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED} seconds)"
