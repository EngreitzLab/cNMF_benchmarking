#!/bin/bash
#SBATCH -J torch-cNMF
#SBATCH -p owners                 
#SBATCH -c 1                  
#SBATCH --mem=16G              
#SBATCH -t 02:00:00
#SBATCH --array=0-3%4          
#SBATCH -o logs/%x_%A_%a.out
#SBATCH -e logs/%x_%A_%a.err

#SBATCH --mail-user ymo@stanford.edu
#SBATCH --mail-type BEGIN
#SBATCH --mail-type END
#SBATCH --mail-type FAIL

set -euo pipefail
mkdir -p logs

# --- Conda env ---
conda activate torch-cnmf

# --- actual code par ---
echo "Running combo: beta_loss=${b}, solver=${s}"
srun -N1 -n1 python /path/to/test.py --beta_loss "${b}" --solver "${s}"