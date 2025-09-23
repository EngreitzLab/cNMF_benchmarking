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
from PIL import Image
import mygene
import os
from PyPDF2 import PdfMerger
import glob
import matplotlib.colors as mcolors
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rcParams["axes.spines.top"] = False
matplotlib.rcParams["axes.spines.right"] = False


import sys
# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from .utilities import convert_adata_with_mygene, convert_with_mygene, merge_pdfs_in_folder

# plot gene expression UMAP given mdata
def plot_umap_per_gene(mdata, Target_Gene, color = 'purple', save_path = None, save_name = None, figsize = (8,6), show = False):

    # read cell and program data
    # mdata = mu.read_h5mu(mdata_path)
    renamed = convert_adata_with_mygene(mdata['rna'])

    # check if gene exist
    gene_name_list = renamed.var_names.tolist()
    if Target_Gene not in gene_name_list:
        print("gene name is not found in mdata")
        return 
    
    # set color
    colors = ['lightgrey', color]
    n_bins = 100
    cmap = mcolors.LinearSegmentedColormap.from_list('custom', colors, N=n_bins)

    if save_name is None:
        save_name = f'Gene Expression for {Target_Gene}'

    fig, ax = plt.subplots(figsize=figsize)

    sc.pl.umap(renamed, color = Target_Gene, title = save_name, cmap = cmap,ax=ax, show=False)

    if save_path and save_name:
        fig.savefig(f"{save_path}/{save_name}.png")

    # Control whether to display the plot
    if show:
        plt.show()
    else:
        plt.close(fig) 

# read the cNMF gene matrix (in txt), plot top x gene in the program 
def plot_top_program_per_gene(program_loading_path, Target_Gene, top_program = 10, species = "human", save_path = None, save_name = None, figsize = (5,8), show = False):

    # read cNMF gene program matrix
    df = pd.read_csv(program_loading_path, sep='\t', index_col=0)

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

    # Control whether to display the plot
    if show:
        plt.show()
    else:
        plt.close(fig)  # Close the figure to free memory if not showing
 
# plot the gene expresion dotplot split by days given mdata
def perturbed_gene_dotplot(mdata, Target_Gene, groupby='sample', save_name = None, save_path = None, figsize = (3,2), show = False):

    # read in adata
    #mdata = mu.read_h5mu(mdata_path)
    renamed = convert_adata_with_mygene(mdata['rna'])
    renamed.var_names_make_unique()

    if save_name is None:
        save_name = f"{Target_Gene} Expression by days"

    # Create the dotplot
    dp = sc.pl.dotplot(renamed, Target_Gene, groupby=groupby, title = save_name, figsize = figsize,
                     swap_axes=True, dendrogram = False,show = False, return_fig = True)

    # make fig object for pdf plotting
    (dp.make_figure())
    fig = dp.fig

    # save fig 
    if save_name and save_path:
        fig.savefig(f"{save_path}/{save_name}.png", bbox_inches='tight')

    # Control whether to display the plot
    if show:
        plt.show()
    else:
        plt.close(fig)

# plot barplot for up/down regulated genes log2FC given results from perturbation analysis, return genes in df
def plot_log2FC(perturb_path, Target, tagert_col_name="target_name", plot_col_name = "program_name", log2fc_col='log2FC', num_item = 5, p_value = 0.05, save_path=None, save_name = None, figsize = (5, 4), show = False):

    # read df
    df = pd.read_csv(perturb_path, sep = "\t")

    # Sort by log2FC
    df_sorted = df.loc[df[tagert_col_name] == Target]
    df_sorted = df_sorted[df_sorted['adj_pval'] < p_value]
    df_sorted = df_sorted.sort_values(by=log2fc_col, ascending=False) 

    
    # Get top and bottom gene
    top_df = df_sorted.head(num_item)
    bottom_df = df_sorted.tail(num_item)
    
    # Combine and add category
    top = top_df.copy()
    #top['category'] = f'Top {num_item} (Upregulated)'
    
    bottom = bottom_df.copy()
    #bottom['category'] = f'Bottom {num_item} (Downregulated)'
    
    # Combine data
    #plot_data = pd.concat([bottom,top])
    #merge_columns = [col for col in bottom.columns if col != 'category']
    plot_data = pd.merge(top, bottom, how='outer')
    plot_data = plot_data.sort_values(by=log2fc_col, ascending=False)

    # Create the plot
    fig = plt.figure(figsize=figsize)
    
    # Create horizontal bar plot
    colors = ['red' if x > 0 else 'blue' for x in plot_data[log2fc_col]]    
    bars = plt.barh(range(len(plot_data)), plot_data[log2fc_col], color=colors, alpha=0.7)
    
    # Customize the plot
    plt.yticks(range(len(plot_data)), plot_data[plot_col_name])
    plt.xlabel('log2FC', fontsize=12)
    plt.ylabel(plot_col_name, fontsize=12)

    #title = f'Top {num_item} and Bottom {num_item} for {tagert}'
    plt.title(save_name, fontweight='bold')
    
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
        fig.savefig(f'{save_path}/{save_name}.png', dpi=300, bbox_inches='tight')
    
    # Control whether to display the plot
    if show:
        plt.show()
    else:
        plt.close(fig)

    return plot_data

# plot one volcone plot given perturbation analysis, return genes in df
def plot_volcano(perturb_path, Target, tagert_col_name="target_name",plot_col_name = "program_name", down_thred_log = -0.05, up_thred_log = 0.05, p_value_thred = 0.05, save_path = None, save_name = None, figsize = (5,4), show = False):

    df = pd.read_csv(perturb_path, sep = "\t")
    df_program = df.loc[df[tagert_col_name] == Target]

    fig, ax = plt.subplots(figsize = figsize)
    
    ax.scatter(x=df_program['log2FC'],y=df_program['adj_pval'].apply(lambda x:-np.log10(x)),s=1,label="Not significant" ,color="grey")

    # highlight down- or up- regulated genes
    down = df_program[(df_program['log2FC']<=down_thred_log)&(df_program['adj_pval']<=p_value_thred)]
    up = df_program[(df_program['log2FC']>=up_thred_log)&(df_program['adj_pval']<=p_value_thred)]

    ax.scatter(x=down['log2FC'],y=down['adj_pval'].apply(lambda x:-np.log10(x)),s=3,label="Down-regulated",color="blue")
    ax.scatter(x=up['log2FC'],y=up['adj_pval'].apply(lambda x:-np.log10(x)),s=3,label="Up-regulated",color="red")

    for i,r in up.iterrows():
        ax.text(x=r['log2FC'],y=-np.log10(r['adj_pval']),s=r[plot_col_name])

    for i,r in down.iterrows():
        ax.text(x=r['log2FC'],y=-np.log10(r['adj_pval']),s=r[plot_col_name])
        
    ax.set_xlabel("log2FC")
    ax.set_ylabel("-log(adj_pval)")
    ax.axvline(down_thred_log,color="grey",linestyle="--")
    ax.axvline(up_thred_log,color="grey",linestyle="--")
    ax.axhline(-np.log10(p_value_thred),color="grey",linestyle="--")
    ax.legend()
    ax.set_title(f"volcano plot for program {Target} on {save_name}")

    if save_path and save_name:
        fig.savefig(f"{save_path}/{save_name}.png")

    # Control whether to display the plot
    if show:
        plt.show()
    else:
        plt.close(fig)

    #return pd.merge(up, down, how='outer')

# given mdata, list of programs to plot, plot dotplot for programs, split by days
def programs_dotplots(mdata, program_loading_path, groupby = "sample", program_list = None, save_path = None, save_name = None, figsize = (5,4), show = False):

    # load data
    df = pd.read_csv(program_loading_path, sep = '\t', index_col = 0 )
    #adata = mu.read_h5mu(mdata_path)

    # make anndata from program loadings
    adata_new = ad.AnnData(X=df.values)
    adata_new.obs[groupby] = mdata['rna'].obs[groupby].values

    if save_name is None:
        save_name = "Program Loadings by Days"
        
    gene_list = adata_new.var_names.tolist()

    if program_list is not None:
        gene_list = list(map(str, program_list))

    if not gene_list:
        blank_img = np.ones((400, 400, 3), dtype=np.uint8) * 255 
        img = Image.fromarray(blank_img)
        full_path = os.path.join(save_path, f"{save_name}.png")
        img.save(full_path)

        return None

    gene_list = gene_list[::-1]

    dp = sc.pl.dotplot(adata_new, gene_list, groupby=groupby,swap_axes=True, dendrogram = False, show = False,
                       return_fig = True) # title = save_name)

    (dp.make_figure())
    fig = dp.fig

    plt.subplots_adjust(top=0.8, bottom=0.15, left=0.15, right=0.85)

    # save fig 
    if save_name and save_path:
        fig.savefig(f"{save_path}/{save_name}.png", bbox_inches='tight')

    # Control whether to display the plot
    if show:
        plt.show()
    else:
        plt.close(fig)     

# plot top x programs for perturbed gene given Target gene and gene loading file path
def analyze_correlations(gene_loading_path, Target, title, top_num = 5, save_path=None, save_name = None, figsize = (10, 8), show = False):

    # read data
    df = pd.read_csv(gene_loading_path, sep = '\t', index_col = 0 )
    df_rename = convert_with_mygene(df, index = False) # reanme

    # Calculate correlation matrix
    correlation_matrix = df_rename.corr()
    
    # Get correlations with the target program
    target_correlations = correlation_matrix[Target]
    target_correlations = target_correlations.drop(Target) # Remove self-correlation
    
    # Sort correlations
    sorted_correlations = target_correlations.sort_values(ascending=False)
    
    # Get top and bottom gene
    top = sorted_correlations.head(top_num)
    bottom = sorted_correlations.tail(top_num)

    # Get the middle gene
    #abs_correlations = sorted_correlations.abs().sort_values()
    #middle = abs_correlations.head(top_num)
    
    # Combine for plotting
    combined_correlations = pd.concat([bottom, top])
    combined_correlations = combined_correlations.sort_values(ascending = False)

    # Create the plot
    fig = plt.figure(figsize=figsize)
    
    # Create horizontal bar plot
    colors = ['blue' if x < 0 else 'red' for x in combined_correlations]
    
    bars = plt.barh(range(len(combined_correlations)), combined_correlations.values, color=colors, alpha=0.7)
    
    # Customize the plot
    plt.title( title, fontsize=14, fontweight='bold', pad=20)
    plt.ylabel('Genes', fontsize=12, fontweight='bold')
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
    legend_elements =[Patch(facecolor='red', alpha=0.7, label='Positive correlation'),
                      Patch(facecolor='blue', alpha=0.7, label='Negative correlation')]
    plt.legend(handles=legend_elements, loc='lower right')
    
    # Adjust layout
    plt.tight_layout()
    plt.grid(axis='x', alpha=0.3)
    
    # Save if path provided
    if save_path:
        fig.savefig(f'{save_path}/{save_name}.png', dpi=300, bbox_inches='tight')

    # Control whether to display the plot
    if show:
        plt.show()
    else:
        plt.close(fig)

    return combined_correlations
   
# graph pdf 
def graph_pdf_gene_QC(mdata, gene_loading_path, program_loading_path, perturb_path, Target_Gene, top_program, save_path, save_name, p_value = 0.05, remove_png = True):

    save_path = f'{save_path}/Perturbed_Gene_QC'

    os.makedirs(f'{save_path}/', exist_ok=True)

    with PdfPages(f"{save_path}/{save_name}.pdf") as pdf:

        plot_umap_per_gene(mdata, Target_Gene, figsize = (4,3), save_path = save_path, save_name = f"{Target_Gene}_UMAP" )
        plot_top_program_per_gene(gene_loading_path, Target_Gene, top_program, figsize = (3,6), save_path = save_path, save_name = f"{Target_Gene}_Top_program")
        perturbed_gene_dotplot(mdata, Target_Gene, figsize = (5,2),  save_path = save_path, save_name = f"{Target_Gene}_DotPlot")
        out = analyze_correlations(gene_loading_path, Target_Gene, f"Gene Loading Correlation for {Target_Gene}",  figsize = (5, 4),  save_path = save_path, save_name = f"{Target_Gene}_Correlated_gene")

        for samp in ['D0', 'sample_D1', 'sample_D2','sample_D3']:
            df = plot_log2FC(f"{perturb_path}_{samp}_perturbation_association.txt", Target = Target_Gene, tagert_col_name = "target_name", plot_col_name = "program_name", p_value = p_value , save_path = save_path, save_name = f"{Target_Gene}_log2FC_{samp}")
            plot_volcano(f"{perturb_path}_{samp}_perturbation_association.txt", Target = Target_Gene, down_thred_log = -0.05, up_thred_log = 0.05, p_value_thred = p_value, save_path = save_path, save_name = f"{Target_Gene}_volcano_{samp}")
            programs_dotplots(mdata, program_loading_path, program_list = df["program_name"].tolist(), save_path = save_path, save_name = f"{Target_Gene}_program_{samp}")
            

        fig, axes = plt.subplots(4, 4, figsize=(30, 30))
        fig.suptitle(save_name, 
                    fontsize=24, fontweight='bold', y=0.95)

        img1 = plt.imread(f"{save_path}/{Target_Gene}_UMAP.png")
        axes[0,0].imshow(img1)
        axes[0,0].axis('off')
        
        img2 = plt.imread(f"{save_path}/{Target_Gene}_Top_program.png") 
        axes[0,1].imshow(img2)
        axes[0,1].axis('off')

        img3 = plt.imread(f"{save_path}/{Target_Gene}_DotPlot.png") 
        axes[0,2].imshow(img3)
        axes[0,2].axis('off')

        img4 = plt.imread(f"{save_path}/{Target_Gene}_Correlated_gene.png") 
        axes[0,3].imshow(img4)
        axes[0,3].axis('off')

        for row, plot_type in enumerate(['log2FC','volcano', 'program'], 1):
            for col, day in enumerate(['D0', 'D1', 'D2', 'D3']):
                suffix = '' if day == 'D0' else '_sample'
                img = plt.imread(f"{save_path}/{Target_Gene}_{plot_type}{suffix}_{day}.png")
                axes[row, col].imshow(img)
                axes[row, col].axis('off')
                
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches='tight')
        plt.close()


        if remove_png:
            os.remove(f"{save_path}/{Target_Gene}_UMAP.png")
            os.remove(f"{save_path}/{Target_Gene}_Top_program.png")
            os.remove(f"{save_path}/{Target_Gene}_DotPlot.png")
            os.remove(f"{save_path}/{Target_Gene}_Correlated_gene.png")

            for plot_type in ['log2FC', 'volcano', 'program']:
                for day in ['D0', 'D1', 'D2', 'D3']:
                    suffix = '' if day == 'D0' else '_sample'
                    os.remove(f"{save_path}/{Target_Gene}_{plot_type}{suffix}_{day}.png")

     