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
import matplotlib.colors as mcolors
from cnmf import cNMF


import sys
# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from .utilities import convert_with_mygene


# read the cNMF gene matrix (in txt), plot top x gene in the program 
def plot_top_gene_per_program(gene_loading_path, target_program = 0, num_gene = 10, species = "human", save_path = None, save_name = None , figsize =(5, 8) ):

    # read cNMF gene program matrix
    df = pd.read_csv(gene_loading_path, sep='\t', index_col=0)
    df_sorted = df.iloc[target_program].nlargest(num_gene) #local the top x gene 

    # rename gene
    df_renamed = convert_with_mygene(df_sorted, species = species)

    # Create the plot
    fig, ax = plt.subplots(figsize=figsize)

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
    ax.set_title(f"Top {num_gene} in program {target_program}" , fontsize=14, fontweight='bold')

    if save_path and save_name:
        fig.savefig(f"{save_path}/{save_name}.png", dpi=300, bbox_inches='tight')

    return fig


# make UMAP of program score
def plot_umap_per_program(mdata_path, program_loading_path, target_program = 0, color = 'purple', save_path = None, save_name = None, figsize = (8,6)):

    # read cell and program data
    mdata = mu.read_h5mu(mdata_path)
    cell_usage = pd.read_csv(program_loading_path, sep='\t', index_col=0)

    adata_plot = mdata['rna'].copy()
    adata_plot.obs['cell_program_loading'] = cell_usage.iloc[:,target_program]
    
    colors = ['lightgrey', color]
    n_bins = 100
    cmap = mcolors.LinearSegmentedColormap.from_list('custom', colors, N=n_bins)

    if save_name is None:
        save_name = f'Program {target_program} Cell Loading'

    fig, ax = plt.subplots(figsize=figsize)

    sc.pl.umap(adata_plot, color = 'cell_program_loading', title = save_name, cmap = cmap, ax=ax)

    if save_path and save_name:
        fig.savefig(f"{save_path}/{save_name}.png", dpi=300, bbox_inches='tight')

    return fig


# find most and least simliar programs
def analyze_program_correlations(program_loading_path, target_progarm = 0, num_program = 5, save_path=None, save_name = None, figsize = (5, 4)):


    df =  pd.read_csv(program_loading_path, sep='\t', index_col=0)

    # Calculate correlation matrix
    correlation_matrix = df.corr()
    
    # Get correlations with the target program
    target_correlations = correlation_matrix.iloc[target_progarm]
    target_correlations = target_correlations.drop(target_correlations.index[target_progarm])# Remove self-correlation
    
    # Sort correlations
    sorted_correlations = target_correlations.sort_values(ascending=False)
    
    # Get top and bottom gene
    top = sorted_correlations.head(num_program)
    bottom = sorted_correlations.tail(num_program)


    # Get the middle gene
    abs_correlations = sorted_correlations.abs().sort_values()
    middle = abs_correlations.head(num_program)
    
    # Combine for plotting
    combined_correlations = pd.concat([bottom, middle, top])

    # Create the plot
    fig = plt.figure(figsize=figsize)
    
    # Create horizontal bar plot
    colors = ['blue' if i < num_program else 'grey' if i < num_program*2 else 'red' for i in range(len(combined_correlations))]
    
    bars = plt.barh(range(len(combined_correlations)), combined_correlations.values, color=colors, alpha=0.7)
    
    # Customize the plot
    title = f'Program Correlations on Top, Middle, Bottom {num_program} on Program {target_progarm}'
    plt.title( title, fontsize=14, fontweight='bold', pad=20)
    plt.ylabel('Programs', fontsize=12, fontweight='bold')
    plt.xlabel('Correlation Coefficient', fontsize=12, fontweight='bold')
    
    # Set x-axis labels
    plt.yticks(range(len(combined_correlations)), 
               combined_correlations.index, 
               rotation=45, 
               ha='right')

    
    # Add a vertical line at x=0
    plt.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
    
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='red', alpha=0.7, label='Positive correlation'),
                      Patch(facecolor='grey', alpha=0.7, label='Low correlation'),
                      Patch(facecolor='blue', alpha=0.7, label='Negative correlation')]
    plt.legend(handles=legend_elements, loc='lower right')
    
    # Adjust layout
    plt.tight_layout()
    plt.grid(axis='x', alpha=0.3)
    plt.show()
    
    # Save if path provided
    if save_path and save_name:
        fig.savefig(f'{save_path}/{save_name}.png', dpi=300, bbox_inches='tight')
   

# make violin plot on cell loadings 
def plot_violin(mdata_path, program_loading_path, target_program, figsize = (4,3), save_path = None, save_name = None):

    # read data
    mdata = mu.read_h5mu(mdata_path)
    adata = mdata['rna']
    usage_norm =  pd.read_csv(program_loading_path, sep='\t', index_col=0)


    # Build dataframe
    df = pd.DataFrame({
        "expression": usage_norm.iloc[:,0],
        "cell_type": adata.obs["sample"].values
    })

    # Compute summary stats per cell_type
    summary = (
        df.groupby("cell_type")
        .agg(
            mean_expr=("expression", "mean"),
            frac_cells=("expression", lambda x: (x > 0).mean())
        )
        .reset_index()
    )

    # Make violinplot
    fig = plt.figure(figsize=figsize)
    ax = sns.violinplot(
        data=df,
        x="cell_type",
        y="expression",
        inner="quartile",
        density_norm='width',
        cut = 0)

    # Annotate mean & fraction above violins
    for i, row in summary.iterrows():
        ax.text(
            i,                      # x position
            df["expression"].max()*1.05,  # just above max y
            f"Mean={row['mean_expr']:.2f}\nFrac={row['frac_cells']:.2f}",
            ha="center", va="bottom", fontsize=9, color="black"
        )

    plt.tight_layout()
    plt.show()

    #ax.set_title("Expression Distribution per Day", fontsize=14, weight="bold")


    if save_path and save_name:
        fig.savefig(f"{save_path}/{save_name}.png", dpi=300, bbox_inches='tight')

    return fig


# plot one volcone plot given file path
def plot_volcano(path, target_program, down_thred_log = -0.05, up_thred_log = 0.05, p_value_thred = 0.05,  save_path = None, save_name = None):

    df = pd.read_csv(path, sep = "\t")
    df_program = df.loc[df["program_name"] == target_program]

    plt.scatter(x=df_program['log2FC'],y=df_program['adj_pval'].apply(lambda x:-np.log10(x)),s=1,label="Not significant")

    # highlight down- or up- regulated genes
    down = df_program[(df_program['log2FC']<=down_thred_log)&(df_program['adj_pval']<=p_value_thred)]
    up = df_program[(df_program['log2FC']>=up_thred_log)&(df_program['adj_pval']<=p_value_thred)]

    plt.scatter(x=down['log2FC'],y=down['adj_pval'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
    plt.scatter(x=up['log2FC'],y=up['adj_pval'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")

    for i,r in up.iterrows():
        plt.text(x=r['log2FC'],y=-np.log10(r['adj_pval']),s=r['target_name'])

    for i,r in down.iterrows():
        plt.text(x=r['log2FC'],y=-np.log10(r['adj_pval']),s=r['target_name'])
        
    plt.xlabel("log2FC")
    plt.ylabel("adj_pval")
    plt.axvline(down_thred_log,color="grey",linestyle="--")
    plt.axvline(up_thred_log,color="grey",linestyle="--")
    plt.axhline(-np.log10(p_value_thred),color="grey",linestyle="--")
    plt.legend()
    plt.title(f"volcano plot for program {target_program} on {file_name}")

    if save_path and save_name:
        fig.savefig(f"{save_path}/{save_name}.png")

    plt.show()
    plt.close()


# plot volcano plots by date
def plot_all_days_valcano(in_folder_name, in_file_name,  down_thred_log = -0.05, up_thred_log = 0.05, p_value_thred = 0.05, program_num = 0, save_path = None, save_name = None):
    for samp in ['D0', 'sample_D1', 'sample_D2', 'sample_D3']:
        path = f"{in_folder_name}/{in_file_name}_{samp}_perturbation_association.txt"
        plt = plot_volcano(path, down_thred_log=down_thred_log, up_thred_log=up_thred_log, p_value_thred=p_value_thred, program_num=program_num,save_path=save_path, save_name = f"{save_name}_{samp}")










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
