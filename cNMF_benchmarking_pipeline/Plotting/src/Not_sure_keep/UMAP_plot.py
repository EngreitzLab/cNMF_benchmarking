import matplotlib.colors as mcolors
import muon as mu 
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import xarray as xr
import scanpy as sc
import anndata as ad
from pathlib import Path
import mygene

import sys
# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from .utilities import convert_adata_with_mygene


def plot_umap_per_program(mdata_path, usage_path, program_num = 0, color = 'purple', folder_name = None, file_name = None, figsize = (8,6)):

    # read cell and program data
    mdata = mu.read_h5mu(mdata_path)
    cell_usage = pd.read_csv(usage_path, sep='\t', index_col=0)

    adata_plot = mdata['rna'].copy()
    adata_plot.obs['cell_program_loading'] = cell_usage.iloc[:,program_num]

    
    colors = ['lightgrey', color]
    n_bins = 100
    cmap = mcolors.LinearSegmentedColormap.from_list('custom', colors, N=n_bins)

    if file_name is None:
        file_name = f'Program {program_num} Cell Loading'

    fig, ax = plt.subplots(figsize=figsize)

    sc.pl.umap(adata_plot, color = 'cell_program_loading', title = file_name, cmap = cmap, ax=ax)

    if folder_name and file_name:
        plt.savefig(f"{folder_name}/{file_name}.png")

    return fig


def plot_umap_per_gene(mdata_path, gene_name, color = 'purple', save_path = None, save_name = None, figsize = (8,6)):

    # read cell and program data
    mdata = mu.read_h5mu(mdata_path)
    renamed = convert_adata_with_mygene(mdata['rna'])

    # check if gene exist
    gene_name_list = renamed.var_names.tolist()
    if gene_name not in gene_name_list:
        print("gene name is not found in mdata")
        return 
    
    # set color
    colors = ['lightgrey', color]
    n_bins = 100
    cmap = mcolors.LinearSegmentedColormap.from_list('custom', colors, N=n_bins)

    if save_name is None:
        save_name = f'Gene Expression for {gene_name}'

    fig, ax = plt.subplots(figsize=figsize)


    sc.pl.umap(renamed, color = gene_name, title = save_name, cmap = cmap,ax=ax)

    if save_path and save_name:
        fig.savefig(f"{save_path}/{save_name}.png")

    return fig

    