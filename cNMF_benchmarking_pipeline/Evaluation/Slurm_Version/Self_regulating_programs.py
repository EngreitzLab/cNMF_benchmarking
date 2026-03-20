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
K_values = [50]
thresholds = [0.2]
conditions = ['D0', 'D1', 'D2', 'D3']
perturb_prefix = 'CRT'
base_results = '/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50'
base_eval = '/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/Eval'

#%%

# Load mdata and build gene name dict
mdata_path = '/oak/stanford/groups/engreitz/Users/ymo/IGVF_ccperturbseq/Result/030526_100k_cells_100iter_allHVG_torch_halsvar_batch_e7_50/adata/cNMF_50_0_2.h5mu'
mdata = mu.read(mdata_path)
#gene_name_dict = dict(zip(mdata['rna'].var['var_names'], mdata['rna'].var['gene_name']))
gene_name_dict = None

#%%
all_hits = []
for K in K_values:
    for thresh in thresholds:
        thresh_str = str(thresh).replace('.', '_')  # 0.2 -> 0_2

        # Get top genes
        spectra_path = f'{base_results}.gene_spectra_score.k_{K}.dt_{thresh_str}.txt'
        top_genes = get_top_genes_per_program(spectra_path, top_n=300, gene_name_dict=gene_name_dict)

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
    result = result[['K', 'threshold', 'Condition', 'target_name', 'program_name', 'log2FC', 'p-value', 'adj_pval']]
    result = result.sort_values(['K', 'threshold', 'Condition', 'program_name', 'target_name']).reset_index(drop=True)

    print(f'\nTotal autoregulatory hits: {len(result)}')
    print(f'\nSummary per K/threshold/Condition:')
    print(result.groupby(['K', 'threshold', 'Condition']).size().to_string())
    print(f'\nAll hits:')
    print(result.to_string(index=False))
# %%

'''


Total autoregulatory hits: 113

Summary per K/threshold/Condition:
K   threshold  Condition
50  0.2        D0           24
               D1           13
               D2           21
               D3           55

All hits:
 K  threshold Condition target_name  program_name    log2FC      p-value     adj_pval
50        0.2        D0       ARIH1             0 -0.910582 3.753340e-04 1.090420e-02
50        0.2        D0        CTR9             0  1.133175 3.789590e-06 2.938252e-04
50        0.2        D0     SMARCA5             0 -0.631355 1.178786e-03 2.601689e-02
50        0.2        D0      EIF4G2             3 -1.085091 4.618692e-04 1.267790e-02
50        0.2        D0       RACK1             4 -1.134124 5.368456e-07 5.386665e-05
50        0.2        D0       EIF5A             7 -0.977236 1.039127e-04 4.015388e-03
50        0.2        D0      EIF4G1             9 -1.227517 1.975334e-09 3.804227e-07
50        0.2        D0       SRPRB             9 -0.631248 6.895156e-05 2.993842e-03
50        0.2        D0       RESF1            14 -0.784283 1.313929e-03 2.801570e-02
50        0.2        D0       HSPA5            15  1.209983 1.302040e-03 2.781186e-02
50        0.2        D0       PSMD4            16  1.050106 2.739126e-03 4.691619e-02
50        0.2        D0      EIF4G1            18 -1.212165 1.833791e-05 1.052698e-03
50        0.2        D0       PSMD4            20  1.308195 2.622540e-04 8.176000e-03
50        0.2        D0       ERCC2            25 -0.903896 6.990602e-04 1.764700e-02
50        0.2        D0       SRSF6            25 -0.854924 1.192428e-03 2.622105e-02
50        0.2        D0      EIF4G1            28 -1.707446 3.392344e-10 7.501077e-08
50        0.2        D0       HUWE1            28 -0.901900 8.992825e-04 2.113732e-02
50        0.2        D0      EIF4G1            29  0.877413 5.037209e-05 2.396259e-03
50        0.2        D0       ASCC3            36 -0.877257 1.355141e-04 4.963459e-03
50        0.2        D0        CTR9            38 -1.544033 1.133527e-08 1.935871e-06
50        0.2        D0      EIF4G2            38  0.751719 1.237717e-03 2.691947e-02
50        0.2        D0      HNRNPM            38  0.614416 2.196316e-03 4.051306e-02
50        0.2        D0       RACK1            38 -0.763075 1.819437e-03 3.624442e-02
50        0.2        D0      ATXN2L            39  0.783958 5.037451e-04 1.354707e-02
50        0.2        D1       RACK1             4 -0.964513 8.676175e-04 4.752493e-02
50        0.2        D1      EIF4G1             9 -1.524582 2.180857e-09 1.168958e-06
50        0.2        D1      BMPR1A            14 -1.562296 1.192970e-04 1.145269e-02
50        0.2        D1        TBX3            14 -1.193532 9.402270e-04 4.977477e-02
50        0.2        D1       PSMD4            16  1.470465 1.711848e-04 1.508326e-02
50        0.2        D1      BMPR1A            17 -1.483385 9.354045e-04 4.977477e-02
50        0.2        D1      BMPR1A            19 -1.675168 8.153317e-05 9.025545e-03
50        0.2        D1       MIXL1            30 -1.881451 1.091718e-09 6.585587e-07
50        0.2        D1        TBX3            30  1.301231 8.300025e-05 9.025545e-03
50        0.2        D1        TBXT            30 -1.525991 5.463281e-06 1.033540e-03
50        0.2        D1       EOMES            33 -2.462652 2.856751e-14 4.593728e-11
50        0.2        D1        SKA3            36  1.896763 7.688503e-05 8.752784e-03
50        0.2        D1       EOMES            48 -2.529282 8.083377e-13 9.453291e-10
50        0.2        D2      BMPR1A             1  0.995491 1.305929e-03 3.394527e-02
50        0.2        D2      EIF4G2             3 -1.233992 6.534385e-04 2.060452e-02
50        0.2        D2      EIF4G1             9 -1.538917 3.823657e-09 7.555675e-07
50        0.2        D2      BMPR1A            14 -3.248052 1.104922e-19 1.017759e-16
50        0.2        D2        TBX3            14 -1.275488 7.424773e-05 3.947516e-03
50        0.2        D2         KDR            15 -1.370310 6.952058e-04 2.163386e-02
50        0.2        D2      TGFBR2            15 -1.317723 7.820724e-05 4.102835e-03
50        0.2        D2      TRIM24            15 -1.679644 8.582694e-06 6.738950e-04
50        0.2        D2      BMPR1A            17 -1.812582 1.608737e-06 1.563403e-04
50        0.2        D2       HAND1            17  1.753230 5.925838e-08 8.782241e-06
50        0.2        D2       RESF1            17  1.999635 2.678566e-09 5.475452e-07
50        0.2        D2       SMAD3            17  1.024700 4.850250e-05 2.889726e-03
50        0.2        D2        TBX3            17 -3.538271 1.942672e-23 3.290387e-20
50        0.2        D2        WWC3            17  0.771836 2.396807e-03 4.984826e-02
50        0.2        D2      BMPR1A            19  2.003454 1.171343e-10 3.019061e-08
50        0.2        D2        TBX3            19  2.034768 5.374253e-14 2.450701e-11
50        0.2        D2      EIF4G1            26  1.274163 3.186871e-05 2.064709e-03
50        0.2        D2        TBX3            30  2.515145 2.047073e-17 1.650510e-14
50        0.2        D2      AMOTL2            34  1.329497 4.977074e-06 4.371051e-04
50        0.2        D2       MED12            37  1.120488 3.387019e-05 2.158988e-03
50        0.2        D2        CTR9            38 -1.198173 3.993409e-04 1.434747e-02
50        0.2        D3        UPF2             0  0.447383 2.186945e-03 1.535917e-02
50        0.2        D3        CDH5             2 -0.991554 5.641431e-03 3.177133e-02
50        0.2        D3        ETS1             2 -0.868688 5.823474e-03 3.246922e-02
50        0.2        D3         KDR             2 -3.034996 4.007331e-17 4.653610e-15
50        0.2        D3      EIF4G2             3 -1.739055 2.423439e-05 3.467073e-04
50        0.2        D3        ETS1             3 -1.090543 3.273046e-04 3.254474e-03
50        0.2        D3         KDR             3 -3.693359 2.562666e-30 8.716229e-28
50        0.2        D3       SMAD1             3 -2.493210 2.786962e-29 8.348084e-27
50        0.2        D3       RACK1             4 -1.148356 4.929701e-06 8.483422e-05
50        0.2        D3      RPS15A             5 -0.854885 7.560625e-03 3.961676e-02
50        0.2        D3        CDH5             8 -1.043087 1.061860e-04 1.265950e-03
50        0.2        D3        EEF2             8 -0.815766 1.329282e-03 1.035058e-02
50        0.2        D3       LENG8             8 -0.746908 1.192229e-03 9.444522e-03
50        0.2        D3       MED12             8 -1.506203 7.588539e-05 9.434296e-04
50        0.2        D3       MED15             8 -1.416065 1.420404e-05 2.171451e-04
50        0.2        D3       PLCG1             8 -2.159853 3.687253e-13 2.356228e-11
50        0.2        D3      SETD1B             8 -1.523487 1.465223e-04 1.668004e-03
50        0.2        D3       SRRM2             8 -1.115100 1.398548e-04 1.601227e-03
50        0.2        D3     ATP6AP1             9 -0.566999 8.200889e-03 4.186932e-02
50        0.2        D3      EIF4G1             9 -1.254606 4.428385e-06 7.731767e-05
50        0.2        D3        CDH5            11 -0.860348 8.442645e-03 4.268211e-02
50        0.2        D3         KDR            11 -2.800128 5.974995e-18 7.636282e-16
50        0.2        D3       LATS2            11 -1.352074 7.946541e-07 1.670394e-05
50        0.2        D3       SMAD1            11 -2.521477 5.116658e-33 2.043529e-30
50        0.2        D3      BMPR1A            14 -1.237129 5.874100e-06 9.912854e-05
50        0.2        D3       PRKD1            14  0.755341 3.540235e-03 2.242843e-02
50        0.2        D3        ETV2            15 -1.561161 6.361791e-06 1.062364e-04
50        0.2        D3         KDR            15 -1.550760 1.688316e-05 2.532554e-04
50        0.2        D3      TRIM24            15 -1.423330 4.555779e-04 4.265923e-03
50        0.2        D3     SMARCA4            16 -0.607950 2.418789e-04 2.517353e-03
50        0.2        D3       HAND1            17  1.517151 2.513121e-06 4.760676e-05
50        0.2        D3       RESF1            17  1.000488 8.122479e-03 4.181332e-02
50        0.2        D3       SMAD3            17  1.347109 4.292028e-09 1.469299e-07
50        0.2        D3        TBX3            17 -0.661805 5.247373e-04 4.799394e-03
50        0.2        D3        WWC3            17  0.675593 3.855482e-03 2.404421e-02
50        0.2        D3      BMPR1A            19  1.412254 9.423821e-06 1.520710e-04
50        0.2        D3        TBX3            19  0.805331 3.951989e-04 3.795494e-03
50        0.2        D3        ETV2            23  1.188579 1.093179e-05 1.737719e-04
50        0.2        D3       SMAD3            24  1.317374 3.876439e-07 8.846864e-06
50        0.2        D3      EIF4G2            28 -0.942098 2.096336e-03 1.491740e-02
50        0.2        D3       EOMES            28 -0.879685 4.118384e-03 2.515994e-02
50        0.2        D3       MIXL1            28 -0.866745 2.153899e-03 1.519776e-02
50        0.2        D3      EIF4G1            29 -0.571889 2.487887e-03 1.707025e-02
50        0.2        D3        TBX3            30  0.516649 5.665296e-03 3.184960e-02
50        0.2        D3        EEF2            31  0.958317 1.784097e-04 1.961136e-03
50        0.2        D3       EWSR1            31 -0.656129 9.472956e-04 7.847582e-03
50        0.2        D3       RACK1            31 -0.522526 8.391746e-03 4.251629e-02
50        0.2        D3       PTBP1            32 -0.800838 8.934782e-03 4.448964e-02
50        0.2        D3        SKA3            32 -1.320857 3.912541e-03 2.430517e-02
50        0.2        D3      AMOTL2            34  3.335193 7.506015e-43 4.232200e-40
50        0.2        D3       SMAD3            34  0.889978 9.254647e-04 7.707087e-03
50        0.2        D3         KDR            35  0.845744 6.925778e-03 3.694249e-02
50        0.2        D3        SKA3            36  1.623695 2.870871e-03 1.917112e-02
50        0.2        D3      EIF4G2            38  0.826548 6.338624e-03 3.457292e-02
50        0.2        D3       SMAD1            42 -1.190789 1.296282e-11 6.731711e-10
'''
# %%
