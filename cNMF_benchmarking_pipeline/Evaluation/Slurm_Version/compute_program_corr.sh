#!/bin/bash
#SBATCH --job-name=prog_corr
#SBATCH --partition=engreitz,owners,bigmem
#SBATCH --mem=512G
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --output=/oak/stanford/groups/engreitz/Users/ymo/Project/combined_final_merged_hvg10k/Results/combined_final_merged_hvg10k/Eval/50_0_5/logs/%j.out
#SBATCH --error=/oak/stanford/groups/engreitz/Users/ymo/Project/combined_final_merged_hvg10k/Results/combined_final_merged_hvg10k/Eval/50_0_5/logs/%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=ymo@stanford.edu

EVAL_DIR="/oak/stanford/groups/engreitz/Users/ymo/Project/combined_final_merged_hvg10k/Results/combined_final_merged_hvg10k/Eval/50_0_5"
PERTURB_BASE="${EVAL_DIR}/50_perturbation_assocition_results"
CONDITIONS="il1b oxldl resting"

eval "$(conda shell.bash hook)" && conda activate NMF_Benchmarking

python -c "
import pandas as pd
import numpy as np

conditions = '${CONDITIONS}'.split()
perturb_base = '${PERTURB_BASE}'
eval_dir = '${EVAL_DIR}'

for cond in conditions:
    print(f'Computing program correlation matrix for {cond}...')
    perturb_path = f'{perturb_base}_{cond}.txt'
    df = pd.read_csv(perturb_path, sep='\t', index_col=0)
    pivot_df = df.pivot_table(index='program_name', columns='target_name', values='log2FC')
    corr_matrix = pivot_df.T.corr()
    np.fill_diagonal(corr_matrix.values, np.nan)
    out_path = f'{eval_dir}/corr_program_matrix_{cond}.txt'
    corr_matrix.to_csv(out_path, sep='\t')
    print(f'Saved to {out_path} (shape: {corr_matrix.shape})')

print('Done.')
"
