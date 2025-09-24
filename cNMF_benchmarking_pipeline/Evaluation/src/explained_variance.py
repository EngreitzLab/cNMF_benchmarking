import numpy as np
import scipy
import pandas as pd
import scanpy as sc
import argparse
import os
import re
from cnmf import cNMF
from sklearn.decomposition import PCA

# for explained variance
def compute_Var(X):
        return np.sum(np.var(X, axis=0, ddof=1))

# for explained variance
def computeVarianceExplained(X, H, Var_X, i):
    if not isinstance(H, (pd.DataFrame)):
        B_k = X @ H[i,:].T / np.sqrt((np.asarray(H[i,:])**2).sum())
        numerator = compute_Var(X - np.outer(B_k, H[i,:]))
    else:
        B_k = X @ H.iloc[i,:] / np.sqrt((H.iloc[i,:]**2).sum())
        numerator = compute_Var(X - np.outer(B_k, H.iloc[i,:]))
    return (1 - numerator / Var_X)


# given cnmf_obj, normalized X array and components k, calculate explained variance
def compute_explained_variance(cnmf_obj, X, k, output_folder, thre = '2.0'):

    # check X is sparce 
    if hasattr(X, 'toarray'):
        X = X.toarray()

    # read median for W and H
    thre_name = (thre).replace('.', '_')
    H_path = cnmf_obj.paths['consensus_spectra__txt'] % (k, thre_name) ## median_spectra_file
    H_df = pd.read_csv(H_path, sep='\t', index_col=0).T
    H = H_df.to_numpy()
    H = (H/H.sum(0))

    W_path = cnmf_obj.paths['consensus_usages__txt'] % (k, thre_name) ## median_spectra_file
    W_df = pd.read_csv(W_path, sep='\t', index_col=0)
    W = W_df.to_numpy()

    WH = W @ H.T
    diff = X - WH
    diff_sumOfSquaresError = (np.asarray(diff)**2).sum()
    Var_diff = compute_Var(diff)
    Var_X = compute_Var(X)
    TotalVarianceExplained = 1 - Var_diff / Var_X
    V_k = np.empty([k])


    # calculate each program's variance 
    for i in range(k):
        V_k[i] = computeVarianceExplained(X, H.T, Var_X, i)

    ProgramID = ['K' + str(k) + '_' + str(i+1) for i in range(k)]

    # save data 
    metrics_df = pd.DataFrame({'VarianceExplained': V_k,
                                'ProgramID': ProgramID })

    metrics_summary = pd.DataFrame({'Sum' : metrics_df['VarianceExplained'].sum(),
                                    'Median' : metrics_df['VarianceExplained'].median(),
                                    'Max' : metrics_df['VarianceExplained'].max(),
                                    'Total' : TotalVarianceExplained},
                                    index = [0])


    metrics_df.to_csv(os.path.join(output_folder, f"{k}_Explained_Variance.txt"), index = None, sep="\t")
    metrics_summary.to_csv(os.path.join(output_folder, f"{k}_Explained_Variance_Summary.txt"), index = None, sep="\t")
 