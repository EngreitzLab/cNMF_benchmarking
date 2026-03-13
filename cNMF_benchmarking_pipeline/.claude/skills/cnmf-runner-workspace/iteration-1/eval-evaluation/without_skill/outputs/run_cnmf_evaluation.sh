#!/bin/bash

# SLURM job configuration
#SBATCH --job-name=030526_100k_cells_eval
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells/Eval/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells/Eval/logs/%j.err
#SBATCH --partition=engreitz
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=96G

# Email notifications
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=ymo@stanford.edu

# Define the cNMF case
OUT_DIR="/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result"
RUN_NAME="030526_100k_cells"
RUN_DIR="$OUT_DIR/$RUN_NAME"

# Store start time
START_TIME=$(date +%s)

# Print some job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURMD_NODENAME"
echo "Working directory: $(pwd)"
echo "Number of CPUs allocated: $SLURM_CPUS_PER_TASK"
echo "Partition: $SLURM_JOB_PARTITION"
echo "Run directory: $RUN_DIR"

# Create logs directory if it doesn't exist
mkdir -p "$RUN_DIR/Eval/logs"

# Activate conda environment
echo "Activating conda environment..."
eval "$(conda shell.bash hook)"
conda activate NMF_Benchmarking

echo "Active conda environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"
echo "Python path: $(which python)"

# Run the evaluation Python script
echo "Running cNMF evaluation pipeline..."
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Slurm_Version/cNMF_evaluation_pipeline.py \
        --out_dir "$OUT_DIR" \
        --run_name "$RUN_NAME" \
        --Perform_categorical \
        --Perform_perturbation \
        --Perform_geneset \
        --data_key 'rna' \
        --prog_key 'cNMF' \
        --categorical_key 'batch' \
        --organism 'human' \
        --guide_annotation_path "/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Data/guide/guide_metadata_v43.tsv" \
        --gwas_data_path '/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz' \
        --reference_gtf_path "/oak/stanford/groups/engreitz/Users/opushkar/genome/IGVFFI9573KOZR.gtf.gz" \
        --sel_thresh 0.2 2.0 \
        --K 50 \
        --FDR_method "StoreyQ"

# Calculate and print elapsed time at the end
END_TIME=$(date +%s)
ELAPSED_TIME=$((END_TIME - START_TIME))
HOURS=$((ELAPSED_TIME / 3600))
MINUTES=$(((ELAPSED_TIME % 3600) / 60))
SECONDS=$((ELAPSED_TIME % 60))

echo "Job completed at: $(date)"
echo "Total elapsed time: ${HOURS}h ${MINUTES}m ${SECONDS}s (${ELAPSED_TIME} seconds)"
