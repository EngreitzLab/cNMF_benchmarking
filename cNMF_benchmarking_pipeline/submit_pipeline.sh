#!/bin/bash
# Submit cNMF pipeline steps with SLURM dependency chaining
# Usage: bash submit_pipeline.sh <step> <script_dir>
#
# Steps:
#   inference      - Run cNMF factorization
#   evaluation     - Run evaluation pipeline
#   plotting       - Run program analysis plotting
#   calibration    - Run perturbation calibration (prompts for method)
#   all            - Chain inference -> evaluation -> plotting
#
# Calibration methods (use with "calibration" step):
#   --crt          - Conditional Randomization Test
#   --utest        - U-test perturbation calibration
#   --matched      - Matched cell program DE
#
# Examples:
#   bash submit_pipeline.sh all /path/to/your/Script
#   bash submit_pipeline.sh calibration --crt /path/to/your/Script
#   bash submit_pipeline.sh plotting /path/to/your/Script

PIPELINE_DIR="$(cd "$(dirname "$0")" && pwd)"
CALIBRATION_DIR="$PIPELINE_DIR/Perturbation_association_calibration/Slurm_version"

STEP="${1:?Error: provide step as first argument (inference|evaluation|plotting|calibration|all)}"
shift

# Parse calibration method flag if present
CALIB_METHOD=""
if [[ "$1" == --crt || "$1" == --utest || "$1" == --matched ]]; then
    CALIB_METHOD="$1"
    shift
fi

SCRIPT_DIR="${1:?Error: provide script directory as last argument}"

INFERENCE_SCRIPT="$SCRIPT_DIR/Inference/sk-cNMF_batch.sh"
EVALUATION_SCRIPT="$SCRIPT_DIR/Evaluation/cNMF_evaluation_pipeline.sh"
PLOTTING_SCRIPT="$SCRIPT_DIR/Plotting/cNMF_program_analysis.sh"

submit_calibration() {
    local dep_flag="$1"
    case "$CALIB_METHOD" in
        --crt)
            CALIB_SCRIPT="$CALIBRATION_DIR/CRT/CRT.sh"
            ;;
        --utest)
            CALIB_SCRIPT="$CALIBRATION_DIR/U-test_perturbation_calibration/U-test_perturbation_calibration.sh"
            ;;
        --matched)
            CALIB_SCRIPT="$CALIBRATION_DIR/Matched_cell_programDE/run.sh"
            ;;
        *)
            echo "Error: specify calibration method: --crt, --utest, or --matched"
            echo "  e.g. bash submit_pipeline.sh calibration --crt /path/to/Script"
            exit 1
            ;;
    esac
    if [[ -n "$dep_flag" ]]; then
        JOB=$(sbatch --parsable --dependency=afterok:"$dep_flag" "$CALIB_SCRIPT")
    else
        JOB=$(sbatch --parsable "$CALIB_SCRIPT")
    fi
    echo "Submitted calibration ($CALIB_METHOD): $JOB"
}

case "$STEP" in
    inference)
        JOB1=$(sbatch --parsable "$INFERENCE_SCRIPT")
        echo "Submitted inference: $JOB1"
        ;;
    evaluation)
        JOB2=$(sbatch --parsable "$EVALUATION_SCRIPT")
        echo "Submitted evaluation: $JOB2"
        ;;
    plotting)
        JOB3=$(sbatch --parsable "$PLOTTING_SCRIPT")
        echo "Submitted plotting: $JOB3"
        ;;
    calibration)
        submit_calibration ""
        ;;
    all)
        JOB1=$(sbatch --parsable "$INFERENCE_SCRIPT")
        echo "Submitted inference: $JOB1"

        JOB2=$(sbatch --parsable --dependency=afterok:"$JOB1" "$EVALUATION_SCRIPT")
        echo "Submitted evaluation: $JOB2 (waits for $JOB1)"

        JOB3=$(sbatch --parsable --dependency=afterok:"$JOB2" "$PLOTTING_SCRIPT")
        echo "Submitted plotting: $JOB3 (waits for $JOB2)"

        echo ""
        echo "Pipeline chain: $JOB1 -> $JOB2 -> $JOB3"
        echo "Monitor: squeue -u $USER"
        ;;
    *)
        echo "Usage: bash submit_pipeline.sh <step> [--crt|--utest|--matched] <script_dir>"
        echo "Steps: inference | evaluation | plotting | calibration | all"
        exit 1
        ;;
esac
