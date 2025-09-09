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

def plot_umap_per_program(mdata_path, usage_path, program_num = 0, color = 'purple', folder_name = None, file_name = None):

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

    sc.pl.umap(adata_plot, color = 'cell_program_loading', title = file_name, cmap = cmap )

    if folder_name and file_name:
        plt.savefig(f"{folder_name}/{file_name}.png")
