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

# plot barplot for top/bottom genes  
def plot_top_bottom_genes(path, gene_col='gene', log2fc_col='log2FC', num_gene = 5, save_path=None):

    # read df
    df = pd.read_csv(path, sep = "\t")

    # Sort by log2FC
    df_sorted = df.sort_values(by=log2fc_col, ascending=False)
    
    # Get top and bottom gene
    top_gene = df_sorted.head(num_gene)
    bottom_gene = df_sorted.tail(num_gene)
    
    # Combine and add category
    top = top_gene.copy()
    top['category'] = f'Top {num_gene} (Upregulated)'
    
    bottom = bottom_gene.copy()
    bottom['category'] = f'Bottom {num_gene} (Downregulated)'
    
    # Combine data
    plot_data = pd.concat([bottom,top])
    
    # Create the plot
    plt.figure(figsize=(10, 8))
    
    # Create horizontal bar plot
    colors = ['red' if x > 0 else 'blue' for x in plot_data[log2fc_col]]
    
    bars = plt.barh(range(len(plot_data)), plot_data[log2fc_col], color=colors, alpha=0.7)
    
    # Customize the plot
    plt.yticks(range(len(plot_data)), plot_data[gene_col])
    plt.xlabel('log2FC', fontsize=12)
    plt.ylabel('Genes', fontsize=12)

    title = f'Top {num_gene} and Bottom {num_gene} Genes by log2FC'
    plt.title(title, fontsize=14, fontweight='bold')
    
    # Add a vertical line at x=0
    plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    
    # Add value labels on bars
    #for i, (bar, value) in enumerate(zip(bars, plot_data[log2fc_col])):
    #    plt.text(value + (0.1 if value > 0 else -0.1), i, f'{value:.2f}', 
    #            va='center', ha='left' if value > 0 else 'right', fontsize=10)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='red', alpha=0.7, label='Upregulated'),
                      Patch(facecolor='blue', alpha=0.7, label='Downregulated')]
    plt.legend(handles=legend_elements, loc='lower right')
    
    # Adjust layout
    plt.tight_layout()
    plt.grid(axis='x', alpha=0.3)
    
    # Save if path provided
    if save_path:
        plt.savefig(f'{title}/{save_path}.png', dpi=300, bbox_inches='tight')
    
    plt.show()
