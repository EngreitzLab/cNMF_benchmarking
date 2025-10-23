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

import sys
# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from .utilities import convert_with_mygene



# read the cNMF gene matrix (in txt), plot top x gene in the program 
def plot_top_gene_per_program(path, program_num = 0, num_gene = 10, species = "human", save_path = None, save_name = None):

    # read cNMF gene program matrix
    df = pd.read_csv(path, sep='\t', index_col=0)
    df_sorted = df.iloc[program_num].nlargest(num_gene) #local the top x gene 

    # rename gene
    df_renamed = convert_with_mygene(df_sorted, species = species)

    # Create the plot
    fig, ax = plt.subplots(figsize=(5, 8))

    # Create horizontal bar plot
    bars = ax.barh(range(len(df_renamed)), df_renamed.values, color='#808080', alpha=0.8)

    # Customize the plot
    ax.set_yticks(range(len(df_renamed)))
    ax.set_yticklabels(df_renamed.index, fontsize=10)
    ax.set_xlabel('Program Specificity (z-score)', fontsize=11)

    # Format x-axis to match your reference
    ax.set_xlim(0, max(df_renamed.values) * 1.1)
    ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

    # Add grid
    ax.grid(axis='x', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # set title
    ax.set_title(f"Top {num_gene} in program {program_num}" , fontsize=14, fontweight='bold')

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.png")

    return fig


# read the cNMF gene matrix (in txt), plot top x gene in the program 
def plot_top_program_per_gene(path, Target_Gene, top_program = 10, species = "human", save_path = None, save_name = None, figsize = (5,8)):

    # read cNMF gene program matrix
    df = pd.read_csv(path, sep='\t', index_col=0)

    # rename gene
    df_renamed = convert_with_mygene(df, species = species, index = False)

    # sort top x program
    df_sorted = df_renamed[Target_Gene].nlargest(top_program) #local the top x gene 


    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)

    # Create horizontal bar plot
    bars = ax.barh(range(len(df_sorted)), df_sorted.values, color='#808080', alpha=0.8)

    # Customize the plot
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted.index, fontsize=10)
    ax.set_xlabel('Gene Loading Score (z-score)', fontsize=11)

    # Format x-axis to match your reference
    ax.set_xlim(0, max(df_sorted.values) * 1.1)
    ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

    # Add grid
    ax.grid(axis='x', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # set title
    ax.set_title(f"Top {top_program} program for {Target_Gene}" , fontsize=14, fontweight='bold')

    if save_path and save_name:
        fig.savefig(f"{save_path}/{save_name}.png")

    return fig 


# plot top x GO terms with path to txt file 
def top_GO_per_program(path, num_term = 10,  program_num = 0, p_value_name = "Adjusted P-value", title = "GO Term", folder_name = None, file_name = None):

    # read txt file
    df = pd.read_csv(path, sep='\t', index_col=0)

    # local to a program and isolate Term
    df_program = df.loc[df['program_name'] == program_num]

    # rename index
    df_program.index = df_program['Term']

    # sort by the smallest p value
    df_sort = df_program[p_value_name].nsmallest(num_term)

    # -log10 tranform
    df_sort_log = -np.log10(df_sort)

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 8))

    # Create horizontal bar plot
    bars = ax.barh(range(len(df_sort_log)), df_sort_log.values, color='#808080', alpha=0.8)

    # Customize the plot
    ax.set_yticks(range(len(df_sort_log)))
    ax.set_yticklabels(df_sort_log.index, fontsize=10)
    ax.set_xlabel('Adjusted P-value(-10)', fontsize=11)

    # Format x-axis to match your reference
    ax.set_xlim(0, max(df_sort_log.values))
    ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

    # Add grid
    ax.grid(axis='x', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add title
    ax.set_title(title, fontsize=14, fontweight='bold')

    # Adjust layout
    plt.tight_layout()

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.png")


# plot top enriched enhancer or promoter
def plot_motif_per_program(path, num_term = 10,  program_num = 0, day = None, title = "Enhancer", p_value_name = "Adjusted P-value", folder_name = None, file_name = None):

    # read txt file
    df = pd.read_csv(path, sep='\t', index_col=0)

    # local to a program and isolate Term
    df_program = df.loc[df['program_name'] == program_num]

    if day is not None:
        df_program = df_program.loc[df_program['sample'] == day]

    # rename index
    df_program.index = df_program['motif']

    # sort by the smallest p value
    df_sort = df_program[p_value_name].nsmallest(num_term)

    # -log10 tranform
    df_sort_log = -np.log10(df_sort)

    # Create the plot
    fig, ax = plt.subplots(figsize=(5, 8))

    # Create horizontal bar plot
    bars = ax.barh(range(len(df_sort_log)), df_sort_log.values, color='#808080', alpha=0.8)

    # Customize the plot
    ax.set_yticks(range(len(df_sort_log)))
    ax.set_yticklabels(df_sort_log.index, fontsize=10)
    ax.set_xlabel('Adjusted P-value(-10)', fontsize=11)

    # Format x-axis to match your reference
    ax.set_xlim(0, max(df_sort_log.values))
    ax.ticklabel_format(style='scientific', axis='x', scilimits=(0,0))

    # Add grid
    ax.grid(axis='x', alpha=0.3, linestyle='-', linewidth=0.5)
    ax.set_axisbelow(True)

    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


    # Add title
    if day is not None:
        ax.set_title(f"{title} at {day} and program {program_num}", fontsize=14, fontweight='bold')
    else:
         ax.set_title(f"{title} at program {program_num}", fontsize=14, fontweight='bold')



    # Adjust layout
    plt.tight_layout()    

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.png")


# plot volcano plots by date
def plot_all_days_motif(path, num_term = 10,  program_num = 0, title = "Enhancer", p_value_name = "Adjusted P-value", folder_name = None, file_name = None):
    for samp in ['D0', 'D1', 'D2', 'D3']:
        plot_motif_per_program(path,num_term=num_term, program_num=program_num,day = samp, p_value_name ="adj_pval",  title = title,folder_name=folder_name,file_name=file_name)
