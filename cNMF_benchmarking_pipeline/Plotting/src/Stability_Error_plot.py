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

import cnmf
import scanpy as sc


# collect all NMF runs 
# important: i must use the same package inference (torch-nmf or sk-nmf) to call combine here
def load_stablity_error_data(output_directory,run_name, components = [30, 50, 60, 80, 100, 200, 250, 300]):

    cnmf_obj = cnmf.cNMF(output_dir=output_directory, name=run_name)

    stats = []
    norm_counts = sc.read(cnmf_obj.paths['normalized_counts'])
    for k in components:
        stats.append(cnmf_obj.consensus(k, skip_density_and_return_after_stats=True, show_clustering=False, close_clustergram_fig=True, norm_counts=norm_counts).stats)

    stats = pd.DataFrame(stats)

    print("min stablity is", stats['silhouette'].min())
    print("max stablity is", stats['silhouette'].max())

    print("min error is",stats['prediction_error'].min())
    print("max error is",stats['prediction_error'].max())

    return stats

# plot NMF stability and error
def plot_stablity_error(stats, folder_name = None, file_name = None):
    # Create the plot with two subplots
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8))

    # Top subplot - Stability (using silhouette as stability metric)
    ax1.plot(stats['k'], stats['silhouette'], 'k-', linewidth=2)
    ax1.set_ylabel('Stability', fontsize=12)
    ax1.set_xlim(0, 310)
    #ax1.set_ylim(0.1, 0.4)  # Adjust based on your data range
    ax1.grid(True, alpha=0.3)
    ax1.tick_params(axis='both', which='major', labelsize=10)

    # Bottom subplot - Error (using prediction_error)
    ax2.plot(stats['k'], stats['prediction_error'], 'k-', linewidth=2)
    ax2.set_xlabel('k', fontsize=12)
    ax2.set_ylabel('Error', fontsize=12)
    ax2.set_xlim(0, 310)
    ax2.grid(True, alpha=0.3)
    ax2.tick_params(axis='both', which='major', labelsize=10)

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.png")
