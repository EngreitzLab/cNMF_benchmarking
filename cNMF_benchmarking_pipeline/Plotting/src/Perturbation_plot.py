import os
import math

import mudata
import scanpy as scp

import numpy as np
import pandas as pd

import seaborn as sns
from matplotlib import pyplot as plt

from statsmodels.stats.multitest import fdrcorrection
from statsmodels.regression.mixed_linear_model import MixedLM

from joblib import Parallel, delayed
from tqdm.auto import tqdm

# load perturbation data 
def load_perturbation_data(folder, pval = 0.000335, components = [30, 50, 60, 80, 100, 200, 250, 300]):
    
    # Compute no. of unique regulators
    test_stats_df = []

    for k in components:  
        # Run perturbation assocation
        for samp in ['D0', 'sample_D1', 'sample_D2', 'sample_D3']:
            test_stats_df_ = pd.read_csv('{}/{}/{}_perturbation_association_results_{}.txt'.format(folder,k,k,samp), sep='\t')
            test_stats_df_['sample'] = samp
            test_stats_df_['K'] = k
            test_stats_df_['fdr'] = fdrcorrection(test_stats_df_['pval'])[1]
            test_stats_df.append(test_stats_df_)
            
    test_stats_df = pd.concat(test_stats_df, ignore_index=True)

    # pring some stats
    plotting_df = test_stats_df.loc[test_stats_df.pval< pval, ['K','target_name']].drop_duplicates().groupby(['K']).count().reset_index()

    print("min regulators is", plotting_df["target_name"].min())
    print("max regulators is", plotting_df["target_name"].max())

    return test_stats_df

# plot perturbation data
def plot_perturbation(test_stats_df, folder_name = None, file_name = None):

    fig, axs = plt.subplots(ncols=1, nrows=2,figsize=(9, 9))

    axs[0].set_title('Unique regulators of progams per sample (pval <=0.000335)', fontsize=10)
    plotting_df = test_stats_df.loc[test_stats_df.pval<=0.000335, ['K', 'sample','target_name']].drop_duplicates().groupby(['K', 'sample']).count().reset_index()
    sns.lineplot(x='K', y='target_name', hue='sample', data=plotting_df, ax=axs[0])
    axs[0].set_ylabel('# Regulators', fontsize=10)
    axs[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left')


    axs[1].set_title('Unique regulators of progams (pval <=0.000335)', fontsize=10)
    plotting_df = test_stats_df.loc[test_stats_df.pval<=0.000335, ['K','target_name']].drop_duplicates().groupby(['K']).count().reset_index()
    sns.lineplot(x='K', y='target_name', data=plotting_df, color='black', ax=axs[1])
    axs[1].set_ylabel('# Regulators', fontsize=10)

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.png")
