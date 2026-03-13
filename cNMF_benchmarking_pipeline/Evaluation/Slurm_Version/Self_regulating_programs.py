#%%
import pandas as pd
import numpy as np
import muon as mu
import os

# %%

def get_top_genes_per_program(spectra_score_path, top_n=300, gene_name_dict=None):
    """
    Read gene spectra score file and return top N genes per program.

    Args:
        gene_name_dict: optional dict mapping ensembl_id -> gene_name.
                        If provided, column names are converted before ranking.

    Returns dict: {program_number (str): set of gene names}
    """
    df = pd.read_csv(spectra_score_path, sep='\t', index_col=0)

    # Rename columns from Ensembl IDs to gene names if dict provided
    if gene_name_dict is not None:
        df.columns = [gene_name_dict.get(x, x) for x in df.columns]

    # Rows = programs (numbered), columns = genes
    top_genes = {}
    for prog_idx in df.index:
        ranked = df.loc[prog_idx].sort_values(ascending=False)
        top_genes[str(prog_idx - 1)] = set(ranked.head(top_n).index)
    return top_genes


def find_autoregulatory_hits(perturb_path, top_genes, pval_col='adj_pval', pval_thresh=0.05):
    """
    Find regulators that appear in the top genes of programs they significantly regulate.

    Returns a DataFrame with columns:
        target_name, program_name, log2FC, p-value, adj_pval
    filtered to rows where:
        1. adj_pval < threshold
        2. target_name is in the top 300 genes of that program
    """
    crt = pd.read_csv(perturb_path, sep='\t')
    # Filter to significant hits
    sig = crt[crt[pval_col] < pval_thresh].copy()
    # Check if regulator is in the top genes of the program it regulates
    mask = sig.apply(
        lambda row: row['target_name'] in top_genes.get(str(row['program_name']), set()),
        axis=1
    )
    return sig[mask].reset_index(drop=True)


#%%
K_values = [50, 150]
thresholds = [0.5]
conditions = ['il1b', 'oxldl', 'resting']
perturb_prefix = 'perturbation_assocition_results'
base_results = '/oak/stanford/groups/engreitz/Users/ymo/Project/combined_final_merged_hvg10k/Results/combined_final_merged_hvg10k'
base_eval = '/oak/stanford/groups/engreitz/Users/ymo/Project/combined_final_merged_hvg10k/Results/combined_final_merged_hvg10k/Eval'

#%%

# Load mdata and build gene name dict
mdata_path = '/oak/stanford/groups/engreitz/Users/ymo/Project/combined_final_merged_hvg10k/Results/combined_final_merged_hvg10k/adata/cNMF_50_0_5.h5mu'
mdata = mu.read(mdata_path)
gene_name_dict = dict(zip(mdata['rna'].var['var_names'], mdata['rna'].var['gene_name']))

#%%
all_hits = []
for K in K_values:
    for thresh in thresholds:
        thresh_str = str(thresh).replace('.', '_')  # 0.2 -> 0_2

        # Get top genes
        spectra_path = f'{base_results}/combined_final_merged_hvg10k.gene_spectra_score.k_{K}.dt_{thresh_str}.txt'
        top_genes = get_top_genes_per_program(spectra_path, top_n=10, gene_name_dict=gene_name_dict)

        # Find autoregulatory hits per condition
        for Condition in conditions:
            crt_path = f'{base_eval}/{K}_{thresh_str}/{K}_{perturb_prefix}_{Condition}.txt'                                                                                         

            if not os.path.exists(crt_path):
                print(f'Warning: {crt_path} not found, skipping')
                continue

            hits = find_autoregulatory_hits(crt_path, top_genes,
                                            pval_col='adj_pval',
                                            pval_thresh=0.05)
            hits['K'] = K
            hits['threshold'] = thresh
            hits['Condition'] = Condition
            all_hits.append(hits)
            print(f'  K={K}, thresh={thresh}, {Condition}: {len(hits)} autoregulatory hits')

# %%
if all_hits:
    result = pd.concat(all_hits, ignore_index=True)
    result = result[['K', 'threshold', 'Condition', 'target_name', 'program_name', 'log2FC', 'pval', 'adj_pval']]
    result = result.sort_values(['K', 'threshold', 'Condition', 'program_name', 'target_name']).reset_index(drop=True)

    print(f'\nTotal autoregulatory hits: {len(result)}')
    print(f'\nSummary per K/threshold/Condition:')
    print(result.groupby(['K', 'threshold', 'Condition']).size().to_string())
    print(f'\nAll hits:')
    print(result.to_string(index=False))
# %%

'''

Total autoregulatory hits: 65

Summary per K/threshold/Condition:
K    threshold  Condition
50   0.5        il1b          9
                oxldl         6
                resting      11
150  0.5        il1b         15
                oxldl         7
                resting      17

All hits:
  K  threshold Condition target_name  program_name    log2FC         pval     adj_pval
 50        0.5      il1b        RPS8             5  0.054901 3.163649e-05 1.410178e-02
 50        0.5      il1b       EGFL7             8  0.033241 1.392668e-04 4.444559e-02
 50        0.5      il1b       CREB5            16 -0.046644 1.472925e-04 4.625178e-02
 50        0.5      il1b        EXT1            16 -0.028513 6.399019e-05 2.475107e-02
 50        0.5      il1b       MEF2A            25 -0.030555 1.084139e-04 3.690184e-02
 50        0.5      il1b    SERPINE1            29 -0.089760 4.245648e-12 9.370200e-09
 50        0.5      il1b      CDKN1A            39  0.081048 9.260574e-33 1.218815e-28
 50        0.5      il1b      DNAJB9            43  0.038182 1.601615e-06 1.126838e-03
 50        0.5      il1b       HSPA5            43  0.186528 2.468044e-24 1.964683e-20
 50        0.5     oxldl       ISG15             2  0.072131 2.155552e-15 1.084413e-11
 50        0.5     oxldl    SERPINE1            21 -0.048272 8.680384e-06 6.240723e-03
 50        0.5     oxldl    SERPINE1            29 -0.078204 9.729111e-11 2.395445e-07
 50        0.5     oxldl      CDKN1A            39  0.067077 2.648329e-26 4.693186e-22
 50        0.5     oxldl        MDM2            39  0.033716 1.826556e-06 1.635545e-03
 50        0.5     oxldl       HSPA5            43  0.150631 9.814610e-27 1.865814e-22
 50        0.5   resting       ISG15             2  0.079431 5.988214e-17 1.639677e-13
 50        0.5   resting        TGM2            10 -0.034663 2.012324e-04 3.815669e-02
 50        0.5   resting        EXT1            16 -0.021554 7.173986e-05 1.710507e-02
 50        0.5   resting        PSAP            26  0.029037 3.135720e-05 8.709015e-03
 50        0.5   resting    SERPINE1            29 -0.045251 3.187161e-08 2.135236e-05
 50        0.5   resting        H1-0            37 -0.034200 3.228321e-06 1.273890e-03
 50        0.5   resting      HMGCS1            38  0.024568 2.160531e-05 6.380076e-03
 50        0.5   resting      CDKN1A            39  0.052515 7.287678e-23 3.661781e-19
 50        0.5   resting        MDM2            39  0.024721 1.862294e-05 5.647647e-03
 50        0.5   resting      DNAJB9            43  0.028376 9.410071e-06 3.184465e-03
 50        0.5   resting       HSPA5            43  0.252984 6.905350e-39 9.695405e-35
150        0.5      il1b       CCNL1             5  0.024730 1.244501e-04 4.717837e-02
150        0.5      il1b       EGFL7            12  0.040449 9.378668e-05 3.814455e-02
150        0.5      il1b    SERPINE1            26 -0.109929 3.823764e-18 2.209427e-14
150        0.5      il1b      ZNF736            26 -0.029530 1.216613e-04 4.641531e-02
150        0.5      il1b       CREB5            36 -0.047464 5.614779e-05 2.576439e-02
150        0.5      il1b     COLEC10            48  0.041988 8.377393e-05 3.507171e-02
150        0.5      il1b      DNAJB9            53  0.052358 2.116817e-13 6.828464e-10
150        0.5      il1b       HSPA5            53  0.185486 7.667291e-25 7.657723e-21
150        0.5      il1b      CDKN1A            55  0.075867 7.391932e-28 9.082875e-24
150        0.5      il1b     TNFAIP3            60  0.033693 3.334904e-05 1.701097e-02
150        0.5      il1b     TNFAIP3            63  0.035335 4.471170e-06 3.229384e-03
150        0.5      il1b        LDLR            72  0.038517 1.512771e-06 1.281032e-03
150        0.5      il1b       STAT1            73  0.031603 9.935857e-06 6.290332e-03
150        0.5      il1b     TMEM45B           119 -0.034049 6.028198e-06 4.157751e-03
150        0.5      il1b       CXCL1           136 -0.033838 3.946633e-05 1.947871e-02
150        0.5     oxldl       ISG15             9  0.059852 1.780206e-11 6.211138e-08
150        0.5     oxldl    SERPINE1            26 -0.087855 5.269612e-13 2.385927e-09
150        0.5     oxldl    SERPINE1            29 -0.049204 1.722209e-05 1.285869e-02
150        0.5     oxldl       HSPA5            53  0.176328 2.330572e-33 9.734732e-29
150        0.5     oxldl      CDKN1A            55  0.067801 3.663293e-28 1.005930e-23
150        0.5     oxldl        MDM2            55  0.030529 9.573698e-06 7.938963e-03
150        0.5     oxldl       STAT1            73  0.062402 2.240774e-15 1.435281e-11
150        0.5   resting       MARK3             3 -0.024845 4.360546e-05 1.331297e-02
150        0.5   resting      ARGLU1             5 -0.018757 1.286378e-05 4.861243e-03
150        0.5   resting       ISG15             9  0.073876 1.728309e-16 5.527518e-13
150        0.5   resting       CDC5L            14  0.025363 8.326359e-06 3.383584e-03
150        0.5   resting       MEF2A            22 -0.023249 1.782607e-04 4.027498e-02
150        0.5   resting        FTH1            23  0.028105 8.414465e-05 2.254041e-02
150        0.5   resting    SERPINE1            26 -0.065748 9.999832e-15 2.592318e-11
150        0.5   resting        EXT1            36 -0.021348 1.352220e-04 3.261572e-02
150        0.5   resting        CALR            53  0.047152 1.992789e-08 1.652524e-05
150        0.5   resting      DNAJB9            53  0.039053 9.765882e-11 1.302213e-07
150        0.5   resting       HSPA5            53  0.274293 3.526634e-43 8.181654e-39
150        0.5   resting        MANF            53  0.032668 5.518141e-08 4.114083e-05
150        0.5   resting      CDKN1A            55  0.053234 1.680677e-24 1.230790e-20
150        0.5   resting        MDM2            55  0.023303 4.834092e-05 1.446793e-02
150        0.5   resting     TNFAIP3            60  0.022923 5.527887e-05 1.612271e-02
150        0.5   resting      HMGCS1            72  0.027117 2.078194e-06 1.023433e-03
150        0.5   resting       STAT1            73  0.043169 4.159476e-10 4.922986e-07

'RPS8', 'EGFL7', 'CREB5', 'EXT1', 'MEF2A', 'SERPINE1', 'CDKN1A',
       'DNAJB9', 'HSPA5', 'ISG15', 'MDM2', 'TGM2', 'PSAP', 'H1-0',
       'HMGCS1', 'CCNL1', 'ZNF736', 'COLEC10', 'TNFAIP3', 'LDLR', 'STAT1',
       'TMEM45B', 'CXCL1', 'MARK3', 'ARGLU1', 'CDC5L', 'FTH1', 'CALR',
       'MANF']

29 in total 

'''
# %%
