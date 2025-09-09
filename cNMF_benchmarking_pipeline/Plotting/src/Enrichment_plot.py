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


# Load data for differeent enrichment test
def load_enrichment_data(folder, components = [30, 50, 60, 80, 100, 200, 250, 300]):

    # loading function 
    def load(k, term, folder):
        # read evaluation results
        df = pd.read_csv(folder, sep='\t')
        df = df.loc[df['Adjusted P-value']<=0.05]
        df['num_programs'] = k
        df['test_term'] = term

        return df

    # collect all the results for each k
    term_df = []

    for k in components:

        term_df.append(load(k, 'go_terms', '{}/{}/{}_GO_term_enrichment.txt'.format(folder,k,k)))
        term_df.append(load(k, 'genesets', '{}/{}/{}_geneset_enrichment.txt'.format(folder,k,k)))
        term_df.append(load(k, 'traits', '{}/{}/{}_trait_enrichment.txt'.format(folder,k,k)))
        #term_df.append(load(k, 'motifs', '{}/{}_motif_enrichment.txt'.format(k),k))   

    term_df = pd.concat(term_df, ignore_index=True)

    # Count unique terms per k
    count_df = pd.DataFrame(index=components, columns=term_df['test_term'].unique()) 

    for k in components:
        for col in count_df.columns:
            count_df.loc[k, col] = term_df.loc[(term_df['num_programs']==k) & (term_df['test_term']==col), 'Term'].unique().shape[0]


    #print out some stats
    print("min go_terms is", count_df['go_terms'].min())
    print("max go_terms is", count_df['go_terms'].max())

    print("min genesets is", count_df['genesets'].min())
    print("max genesets is", count_df['genesets'].max())

    print("min traits is", count_df['traits'].min())
    print("max traits is", count_df['traits'].max())


    return count_df

# plot loaded df
def plot_enrichment_data(count_df, folder_name = None, file_name = None):

    # Create the plot with two subplots
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 8))

    # Top subplot - Stability (using silhouette as stability metric)
    ax1.plot(count_df.index, count_df['go_terms'], 'k-', linewidth=2)
    ax1.set_ylabel('GO Terms', fontsize=12)
    ax1.set_xlim(0, 310)
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both', which='major', labelsize=10)

    # Bottom subplot - Error (using prediction_error)
    ax2.plot(count_df.index, count_df['genesets'], 'k-', linewidth=2)
    ax2.set_xlabel('k', fontsize=12)
    ax2.set_ylabel('Genesets', fontsize=12)
    ax2.set_xlim(0, 310)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both', which='major', labelsize=10)

    # Bottom subplot - Error (using prediction_error)
    ax3.plot(count_df.index, count_df['traits'], 'k-', linewidth=2)
    ax3.set_xlabel('k', fontsize=12)
    ax3.set_ylabel('Traits', fontsize=12)
    ax3.set_xlim(0, 310)
    ax3.grid(True, alpha=0.3)
    ax3.tick_params(axis='both', which='major', labelsize=10)

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.png")