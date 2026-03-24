import os
import math
import mudata
import scanpy as sc
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import matplotlib.ticker as ticker
from statsmodels.stats.multitest import fdrcorrection
from statsmodels.regression.mixed_linear_model import MixedLM
from joblib import Parallel, delayed
from tqdm.auto import tqdm
import cnmf
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples
import scipy.sparse as sp
import anndata as ad


# Publication-quality color palette — one distinct color per plot
_COLORS = {
    'primary': '#2c3e50',        # dark slate for line strokes
    'stability': '#3498db',      # blue
    'error': '#e74c3c',          # red
    'go_terms': '#2ecc71',       # green
    'genesets': '#d35400',       # burnt orange
    'traits': '#9b59b6',         # purple
    'perturbation': '#1abc9c',   # teal
    'explained_var': '#f1c40f',  # yellow
}

# Qualitative palette for multi-line plots
_PALETTE = ['#3498db', '#e74c3c', '#2ecc71', '#9b59b6', '#f39c12', '#1abc9c', '#e67e22', '#34495e']


def _style_ax(ax, xlabel=None, ylabel=None, title=None, hide_top_right=True):
    """Apply publication-quality styling to an axis."""
    if hide_top_right:
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)
    ax.tick_params(axis='both', which='major', labelsize=9, width=0.8, length=4,
                   direction='out', color='#333333')
    ax.tick_params(axis='both', which='minor', width=0.5, length=2,
                   direction='out', color='#333333')
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=11, fontweight='medium', labelpad=6)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=11, fontweight='medium', labelpad=6)
    if title:
        ax.set_title(title, fontsize=12, fontweight='semibold', pad=8)


# collect all NMF runs
# important: i must use the same package inference (torch-nmf or sk-nmf) to call combine here
def load_stablity_error_data(output_directory, run_name, components=[30, 50, 60, 80, 100, 200, 250, 300]):

    cache_path = os.path.join(output_directory, run_name, f'{run_name}_stability_error_cache.tsv')

    if os.path.exists(cache_path):
        print(f"Loading cached stability/error data from {cache_path}")
        stats = pd.read_csv(cache_path, sep='\t')
    else:
        cnmf_obj = cnmf.cNMF(output_dir=output_directory, name=run_name)

        stats = []
        norm_counts = sc.read(cnmf_obj.paths['normalized_counts'])
        for k in components:
            stats.append(cnmf_obj.consensus(k, skip_density_and_return_after_stats=True, show_clustering=False,
            close_clustergram_fig=True, norm_counts=norm_counts, density_threshold = 2.0,local_neighborhood_size = 0.3).stats)

        stats = pd.DataFrame(stats)
        stats.to_csv(cache_path, sep='\t', index=False)
        print(f"Saved stability/error data to {cache_path}")

    print("min stablity is", stats['silhouette'].min())
    print("max stablity is", stats['silhouette'].max())

    print("min error is",stats['prediction_error'].min())
    print("max error is",stats['prediction_error'].max())

    return stats


# plot NMF stability and error — 2 independent square figures
def plot_stablity_error(stats, folder_name=None, file_name=None, selected_k=None):
    # --- Stability ---
    fig1, ax1 = plt.subplots(figsize=(4, 3.5))
    ax1.plot(stats['k'], stats['silhouette'], color=_COLORS['primary'], linewidth=1.5,
             marker='o', markersize=4, markerfacecolor=_COLORS['stability'], markeredgecolor=_COLORS['primary'],
             markeredgewidth=0.8, zorder=3)
    ax1.fill_between(stats['k'], stats['silhouette'], alpha=0.08, color=_COLORS['stability'])
    if selected_k is not None:
        ax1.axvline(x=selected_k, color='red', linestyle='--', linewidth=1, zorder=2)
    _style_ax(ax1, xlabel='Number of components (k)', ylabel='Stability (silhouette)',
              title='Program stability')

    if folder_name and file_name:
        fig1.savefig(f"{folder_name}/{file_name}_stability.svg", bbox_inches="tight")
        fig1.savefig(f"{folder_name}/{file_name}_stability.png", dpi=300, bbox_inches="tight")

    plt.show()
    plt.close(fig1)

    # --- Error ---
    fig2, ax2 = plt.subplots(figsize=(4, 3.5))
    ax2.plot(stats['k'], stats['prediction_error'], color=_COLORS['primary'], linewidth=1.5,
             marker='s', markersize=4, markerfacecolor=_COLORS['error'], markeredgecolor=_COLORS['primary'],
             markeredgewidth=0.8, zorder=3)
    ax2.fill_between(stats['k'], stats['prediction_error'], alpha=0.08, color=_COLORS['error'])
    if selected_k is not None:
        ax2.axvline(x=selected_k, color='red', linestyle='--', linewidth=1, zorder=2)
    _style_ax(ax2, xlabel='Number of components (k)', ylabel='Prediction error',
              title='Reconstruction error')

    if folder_name and file_name:
        fig2.savefig(f"{folder_name}/{file_name}_error.svg", bbox_inches="tight")
        fig2.savefig(f"{folder_name}/{file_name}_error.png", dpi=300, bbox_inches="tight")

    plt.show()
    plt.close(fig2)





# Load data for differeent enrichment test
def load_enrichment_data(folder, components = [30, 50, 60, 80, 100, 200, 250, 300], sel_thresh = 2.0, pval = 0.05):

    # loading function 
    def load(k, term, folder):
        # read evaluation results
        df = pd.read_csv(folder, sep='\t')
        df = df.loc[df['Adjusted P-value']<=pval]
        df['num_programs'] = k
        df['test_term'] = term

        return df

    # collect all the results for each k
    term_df = []

    for k in components:

        term_df.append(load(k, 'go_terms', '{}/{}_{}/{}_GO_term_enrichment.txt'.format(folder,k,str(sel_thresh).replace('.','_'),k)))
        term_df.append(load(k, 'genesets', '{}/{}_{}/{}_geneset_enrichment.txt'.format(folder,k,str(sel_thresh).replace('.','_'),k)))
        term_df.append(load(k, 'traits', '{}/{}_{}/{}_trait_enrichment.txt'.format(folder,k,str(sel_thresh).replace('.','_'),k)))
        #term_df.append(load(k, 'motifs', '{}/{}_{}/{}_motif_enrichment.txt'.format(folder,k,str(sel_thresh).replace('.','_'),k)))

    term_df = pd.concat(term_df, ignore_index=True)

    # Count unique terms per k
    count_df = pd.DataFrame(index=components, columns=term_df['test_term'].unique()) 

    for k in components:
        for col in count_df.columns:
            count_df.loc[k, col] = term_df.loc[(term_df['num_programs']==k) & (term_df['test_term']==col), 'Term'].unique().shape[0]


    #print out some stats
    print(f"min go_terms for {sel_thresh} is", count_df['go_terms'].min())
    print(f"max go_terms for {sel_thresh} is", count_df['go_terms'].max())

    print(f"min genesets for {sel_thresh} is", count_df['genesets'].min())
    print(f"max genesets for {sel_thresh} is", count_df['genesets'].max())

    print(f"min traits for {sel_thresh}  is", count_df['traits'].min())
    print(f"max traits for {sel_thresh}  is", count_df['traits'].max())

    #print("min motif is", count_df['motifs'].min())
    #print("max motif is", count_df['motifs'].max())


    return count_df


# plot loaded df — 3 independent square figures
def plot_enrichment(count_df, folder_name=None, file_name=None, selected_k=None):
    enrichment_config = [
        ('go_terms',  'Unique GO terms',   'GO term enrichment',   _COLORS['go_terms'],  'o'),
        ('genesets',  'Unique gene sets',   'Gene set enrichment',  _COLORS['genesets'],  's'),
        ('traits',    'Unique traits', 'Trait enrichment', _COLORS['traits'],   'D'),
    ]

    for col, ylabel, title, color, marker in enrichment_config:
        fig, ax = plt.subplots(figsize=(4, 3.5))

        vals = count_df[col].astype(float)
        ax.plot(count_df.index, vals, color=_COLORS['primary'], linewidth=1.5,
                marker=marker, markersize=4, markerfacecolor=color,
                markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
        ax.fill_between(count_df.index, vals, alpha=0.08, color=color)
        if selected_k is not None:
            ax.axvline(x=selected_k, color='red', linestyle='--', linewidth=1, zorder=2)
        _style_ax(ax, xlabel='Number of components (k)', ylabel=ylabel, title=title)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

        if folder_name and file_name:
            fig.savefig(f"{folder_name}/{file_name}_{col}.svg", bbox_inches="tight")
            fig.savefig(f"{folder_name}/{file_name}_{col}.png", dpi=300, bbox_inches="tight")

        plt.show()
        plt.close(fig)


# load perturbation data 
def load_perturbation_data(folder, pval = 0.000335, components = [30, 50, 60, 80, 100, 200, 250, 300], sel_thresh = 2.0,
 samples = ['D0', 'sample_D1', 'sample_D2', 'sample_D3']):

    # Compute no. of unique regulators
    test_stats_df = []

    for k in components:
        # Run perturbation assocation
        for samp in samples:
            test_stats_df_ = pd.read_csv('{}/{}_{}/{}_perturbation_association_results_{}.txt'.format(folder,k,str(sel_thresh).replace('.','_'),k,samp), sep='\t')
            test_stats_df_['sample'] = samp
            test_stats_df_['K'] = k
            #test_stats_df_['fdr'] = fdrcorrection(test_stats_df_['pval'])[1]
            test_stats_df.append(test_stats_df_)

    test_stats_df = pd.concat(test_stats_df, ignore_index=True)

    # pring some stats
    plotting_df = test_stats_df.loc[test_stats_df.adj_pval < pval, ['K','target_name']].drop_duplicates().groupby(['K']).count().reset_index()

    print("min regulators is", plotting_df["target_name"].min())
    print("max regulators is", plotting_df["target_name"].max())

    return test_stats_df


# plot perturbation data — 2 independent square figures
def plot_perturbation(test_stats_df, pval=0.000335, folder_name=None, file_name=None, selected_k=None):
    pval_str = f'{pval:.1e}' if pval < 0.01 else str(pval)

    # --- Per-sample plot ---
    fig1, ax1 = plt.subplots(figsize=(4, 3.5))

    plotting_df_sample = (test_stats_df
        .loc[test_stats_df.adj_pval <= pval, ['K', 'sample', 'target_name']]
        .drop_duplicates()
        .groupby(['K', 'sample']).count().reset_index())
    sns.lineplot(x='K', y='target_name', hue='sample', data=plotting_df_sample,
                 palette=_PALETTE, linewidth=1.5, marker='o', markersize=4,
                 ax=ax1, legend='brief')
    if selected_k is not None:
        ax1.axvline(x=selected_k, color='red', linestyle='--', linewidth=1, zorder=2)
    _style_ax(ax1, xlabel='Number of components (k)', ylabel='No. unique regulators',
              title=f'Unique regulators per sample') # (adj. p-value \u2264 {pval_str})')
    ax1.legend(title='Sample', fontsize=8, title_fontsize=9,
               frameon=True, fancybox=False, edgecolor='#cccccc',
               bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

    if folder_name and file_name:
        fig1.savefig(f"{folder_name}/{file_name}_per_sample.svg", bbox_inches="tight")
        fig1.savefig(f"{folder_name}/{file_name}_per_sample.png", dpi=300, bbox_inches="tight")

    plt.show()
    plt.close(fig1)

    # --- Aggregated plot ---
    fig2, ax2 = plt.subplots(figsize=(4, 3.5))

    plotting_df = (test_stats_df
        .loc[test_stats_df.adj_pval <= pval, ['K', 'target_name']]
        .drop_duplicates()
        .groupby(['K']).count().reset_index())
    ax2.plot(plotting_df['K'], plotting_df['target_name'],
             color=_COLORS['primary'], linewidth=1.5,
             marker='o', markersize=4, markerfacecolor=_COLORS['perturbation'],
             markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax2.fill_between(plotting_df['K'], plotting_df['target_name'],
                     alpha=0.08, color=_COLORS['perturbation'])
    if selected_k is not None:
        ax2.axvline(x=selected_k, color='red', linestyle='--', linewidth=1, zorder=2)
    _style_ax(ax2, xlabel='Number of components (k)', ylabel='No. unique regulators',
              title=f'Unique regulators for all samples') #(adj. p-value \u2264 {pval_str})')
    ax2.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

    if folder_name and file_name:
        fig2.savefig(f"{folder_name}/{file_name}_all_samples.svg", bbox_inches="tight")
        fig2.savefig(f"{folder_name}/{file_name}_all_samples.png", dpi=300, bbox_inches="tight")

    plt.show()
    plt.close(fig2)

    return plotting_df



# load total explained variance
def load_explained_variance_data(folder, components = [30, 50, 60, 80, 100, 200, 250, 300], sel_thresh = 2.0):

    stats = {}
    for k in components:
    
        input_path = f"{folder}/{k}_{str(sel_thresh).replace('.','_')}/{k}_Explained_Variance_Summary.txt"
        df = pd.read_csv(input_path, sep = '\t', index_col = 0)

        stats[k] = df['Total'].values[0]

    
    print("min Explained_variance is", min(stats.values()))
    print("max Explained_variance is", max(stats.values()))

    return stats


# plot NMF explained variance
def plot_explained_variance(stats, folder_name=None, file_name=None, selected_k=None):
    ks = list(stats.keys())
    vals = list(stats.values())

    fig, ax = plt.subplots(figsize=(4, 3.5))

    ax.plot(ks, vals, color=_COLORS['primary'], linewidth=1.5,
            marker='o', markersize=4, markerfacecolor=_COLORS['explained_var'],
            markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax.fill_between(ks, vals, alpha=0.08, color=_COLORS['explained_var'])
    if selected_k is not None:
        ax.axvline(x=selected_k, color='red', linestyle='--', linewidth=1, zorder=2)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='Total explained variance',
              title='Explained variance')

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.svg", bbox_inches="tight")
        fig.savefig(f"{folder_name}/{file_name}.png", dpi=300, bbox_inches="tight")

    plt.show()
    plt.close(fig)



# Combined panel figure: 3 rows x 3 columns
# Row 1: Stability, Error, Explained variance
# Row 2: GO terms, Gene sets, Traits
# Row 3: Regulators (all samples), Regulators (per sample)
def plot_k_selection_panel(stability_stats, count_df, test_stats_df, explained_var_stats,
                           pval=0.05, folder_name=None, file_name=None, selected_k=None):

    # Keep text as editable text in SVG (for Adobe Illustrator)
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['pdf.fonttype'] = 42

    fig, axes = plt.subplots(3, 3, figsize=(13, 10.5))
    fig.subplots_adjust(hspace=0.45, wspace=0.35)

    pval_str = f'{pval:.1e}' if pval < 0.01 else str(pval)

    def _add_vline(ax):
        if selected_k is not None:
            ax.axvline(x=selected_k, color='red', linestyle='--', linewidth=1, zorder=2)

    # --- Row 1, Col 0: Stability ---
    ax = axes[0, 0]
    ax.plot(stability_stats['k'], stability_stats['silhouette'], color=_COLORS['primary'], linewidth=1.5,
            marker='o', markersize=4, markerfacecolor=_COLORS['stability'],
            markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax.fill_between(stability_stats['k'], stability_stats['silhouette'], alpha=0.08, color=_COLORS['stability'])
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='Stability (silhouette)', title='Program stability')

    # --- Row 1, Col 1: Error ---
    ax = axes[0, 1]
    ax.plot(stability_stats['k'], stability_stats['prediction_error'], color=_COLORS['primary'], linewidth=1.5,
            marker='s', markersize=4, markerfacecolor=_COLORS['error'],
            markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax.fill_between(stability_stats['k'], stability_stats['prediction_error'], alpha=0.08, color=_COLORS['error'])
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='Prediction error', title='Reconstruction error')

    # --- Row 1, Col 2: Explained variance ---
    ax = axes[0, 2]
    ev_ks = list(explained_var_stats.keys())
    ev_vals = list(explained_var_stats.values())
    ax.plot(ev_ks, ev_vals, color=_COLORS['primary'], linewidth=1.5,
            marker='o', markersize=4, markerfacecolor=_COLORS['explained_var'],
            markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax.fill_between(ev_ks, ev_vals, alpha=0.08, color=_COLORS['explained_var'])
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='Total explained variance', title='Explained variance')

    # --- Row 2: Enrichment (GO, genesets, traits) ---
    enrichment_config = [
        (0, 'go_terms',  'Unique GO terms',   'GO term enrichment',   _COLORS['go_terms'],  'o'),
        (1, 'genesets',  'Unique gene sets',   'Gene set enrichment',  _COLORS['genesets'],  's'),
        (2, 'traits',    'Unique traits',      'Trait enrichment',     _COLORS['traits'],    'D'),
    ]
    for col_idx, col, ylabel, title, color, marker in enrichment_config:
        ax = axes[1, col_idx]
        vals = count_df[col].astype(float)
        ax.plot(count_df.index, vals, color=_COLORS['primary'], linewidth=1.5,
                marker=marker, markersize=4, markerfacecolor=color,
                markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
        ax.fill_between(count_df.index, vals, alpha=0.08, color=color)
        _add_vline(ax)
        _style_ax(ax, xlabel='Number of components (k)', ylabel=ylabel, title=title)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

    # --- Row 3, Col 0: Regulators (all samples) ---
    ax = axes[2, 0]
    plotting_df = (test_stats_df
        .loc[test_stats_df.adj_pval <= pval, ['K', 'target_name']]
        .drop_duplicates()
        .groupby(['K']).count().reset_index())
    ax.plot(plotting_df['K'], plotting_df['target_name'], color=_COLORS['primary'], linewidth=1.5,
            marker='o', markersize=4, markerfacecolor=_COLORS['perturbation'],
            markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax.fill_between(plotting_df['K'], plotting_df['target_name'], alpha=0.08, color=_COLORS['perturbation'])
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='No. unique regulators',
              title='Unique regulators for all samples')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

    # --- Row 3, Col 1: Regulators (per sample) ---
    ax = axes[2, 1]
    plotting_df_sample = (test_stats_df
        .loc[test_stats_df.adj_pval <= pval, ['K', 'sample', 'target_name']]
        .drop_duplicates()
        .groupby(['K', 'sample']).count().reset_index())
    sns.lineplot(x='K', y='target_name', hue='sample', data=plotting_df_sample,
                 palette=_PALETTE, linewidth=1.5, marker='o', markersize=4,
                 ax=ax, legend='brief')
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='No. unique regulators',
              title='Unique regulators per sample')
    ax.legend(title='Sample', fontsize=7, title_fontsize=8,
              frameon=True, fancybox=False, edgecolor='#cccccc',
              bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

    # --- Row 3, Col 2: hide empty panel ---
    axes[2, 2].set_visible(False)

    # Add panel labels (A, B, C, ...)
    panel_labels = 'ABCDEFGH'
    visible_axes = [axes[r, c] for r in range(3) for c in range(3) if axes[r, c].get_visible()]
    for label, ax in zip(panel_labels, visible_axes):
        ax.text(-0.15, 1.12, label, transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top', ha='left')

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.svg", bbox_inches="tight")
        fig.savefig(f"{folder_name}/{file_name}.png", dpi=300, bbox_inches="tight")

    plt.show()
    plt.close(fig)


# Combined panel figure: 3 rows x 3 columns (no trait enrichment)
# Row 1: Stability, Error, Explained variance
# Row 2: GO terms, Gene sets
# Row 3: Regulators (all samples), Regulators (per sample)
def plot_k_selection_panel_no_traits(stability_stats, count_df, test_stats_df, explained_var_stats,
                                     pval=0.05, folder_name=None, file_name=None, selected_k=None):

    # Keep text as editable text in SVG (for Adobe Illustrator)
    plt.rcParams['svg.fonttype'] = 'none'
    plt.rcParams['pdf.fonttype'] = 42

    fig, axes = plt.subplots(3, 3, figsize=(13, 10.5))
    fig.subplots_adjust(hspace=0.45, wspace=0.35)

    pval_str = f'{pval:.1e}' if pval < 0.01 else str(pval)

    def _add_vline(ax):
        if selected_k is not None:
            ax.axvline(x=selected_k, color='red', linestyle='--', linewidth=1, zorder=2)

    # --- Row 1, Col 0: Stability ---
    ax = axes[0, 0]
    ax.plot(stability_stats['k'], stability_stats['silhouette'], color=_COLORS['primary'], linewidth=1.5,
            marker='o', markersize=4, markerfacecolor=_COLORS['stability'],
            markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax.fill_between(stability_stats['k'], stability_stats['silhouette'], alpha=0.08, color=_COLORS['stability'])
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='Stability (silhouette)', title='Program stability')

    # --- Row 1, Col 1: Error ---
    ax = axes[0, 1]
    ax.plot(stability_stats['k'], stability_stats['prediction_error'], color=_COLORS['primary'], linewidth=1.5,
            marker='s', markersize=4, markerfacecolor=_COLORS['error'],
            markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax.fill_between(stability_stats['k'], stability_stats['prediction_error'], alpha=0.08, color=_COLORS['error'])
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='Prediction error', title='Reconstruction error')

    # --- Row 1, Col 2: Explained variance ---
    ax = axes[0, 2]
    ev_ks = list(explained_var_stats.keys())
    ev_vals = list(explained_var_stats.values())
    ax.plot(ev_ks, ev_vals, color=_COLORS['primary'], linewidth=1.5,
            marker='o', markersize=4, markerfacecolor=_COLORS['explained_var'],
            markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax.fill_between(ev_ks, ev_vals, alpha=0.08, color=_COLORS['explained_var'])
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='Total explained variance', title='Explained variance')

    # --- Row 2: Enrichment (GO, genesets) — no traits ---
    enrichment_config = [
        (0, 'go_terms',  'Unique GO terms',   'GO term enrichment',   _COLORS['go_terms'],  'o'),
        (1, 'genesets',  'Unique gene sets',   'Gene set enrichment',  _COLORS['genesets'],  's'),
    ]
    for col_idx, col, ylabel, title, color, marker in enrichment_config:
        ax = axes[1, col_idx]
        vals = count_df[col].astype(float)
        ax.plot(count_df.index, vals, color=_COLORS['primary'], linewidth=1.5,
                marker=marker, markersize=4, markerfacecolor=color,
                markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
        ax.fill_between(count_df.index, vals, alpha=0.08, color=color)
        _add_vline(ax)
        _style_ax(ax, xlabel='Number of components (k)', ylabel=ylabel, title=title)
        ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

    # --- Row 2, Col 2: hide empty panel ---
    axes[1, 2].set_visible(False)

    # --- Row 3, Col 0: Regulators (all samples) ---
    ax = axes[2, 0]
    plotting_df = (test_stats_df
        .loc[test_stats_df.adj_pval <= pval, ['K', 'target_name']]
        .drop_duplicates()
        .groupby(['K']).count().reset_index())
    ax.plot(plotting_df['K'], plotting_df['target_name'], color=_COLORS['primary'], linewidth=1.5,
            marker='o', markersize=4, markerfacecolor=_COLORS['perturbation'],
            markeredgecolor=_COLORS['primary'], markeredgewidth=0.8, zorder=3)
    ax.fill_between(plotting_df['K'], plotting_df['target_name'], alpha=0.08, color=_COLORS['perturbation'])
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='No. unique regulators',
              title='Unique regulators for all samples')
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

    # --- Row 3, Col 1: Regulators (per sample) ---
    ax = axes[2, 1]
    plotting_df_sample = (test_stats_df
        .loc[test_stats_df.adj_pval <= pval, ['K', 'sample', 'target_name']]
        .drop_duplicates()
        .groupby(['K', 'sample']).count().reset_index())
    sns.lineplot(x='K', y='target_name', hue='sample', data=plotting_df_sample,
                 palette=_PALETTE, linewidth=1.5, marker='o', markersize=4,
                 ax=ax, legend='brief')
    _add_vline(ax)
    _style_ax(ax, xlabel='Number of components (k)', ylabel='No. unique regulators',
              title='Unique regulators per sample')
    ax.legend(title='Sample', fontsize=7, title_fontsize=8,
              frameon=True, fancybox=False, edgecolor='#cccccc',
              bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True, nbins=5))

    # --- Row 3, Col 2: hide empty panel ---
    axes[2, 2].set_visible(False)

    # Add panel labels (A, B, C, ...)
    panel_labels = 'ABCDEFG'
    visible_axes = [axes[r, c] for r in range(3) for c in range(3) if axes[r, c].get_visible()]
    for label, ax in zip(panel_labels, visible_axes):
        ax.text(-0.15, 1.12, label, transform=ax.transAxes,
                fontsize=14, fontweight='bold', va='top', ha='left')

    if folder_name and file_name:
        fig.savefig(f"{folder_name}/{file_name}.svg", bbox_inches="tight")
        fig.savefig(f"{folder_name}/{file_name}.png", dpi=300, bbox_inches="tight")

    plt.show()
    plt.close(fig)


''' Faster stability + error calculation & error calculation is a bit off 
because the refit_usages is not the same as the saved one, it needs to refit but it takes time 
def load_stablity_error_data_(k, cnmf_obj, norm_counts, rf_usages_path, density_threshold, local_neighborhood_size):

    # load files 
    merged_spectra = cnmf.load_df_from_npz(cnmf_obj.paths['merged_spectra']%k)
    rf_usages = pd.read_csv(rf_usages_path, sep="\t")
    
    # calculate stability 
    n_neighbors = int(local_neighborhood_size * merged_spectra.shape[0]/k)
    l2_spectra = (merged_spectra.T/np.sqrt((merged_spectra**2).sum(axis=1))).T

    kmeans_model = KMeans(n_clusters=k, n_init=10, random_state=1)
    kmeans_model.fit(l2_spectra)
    kmeans_cluster_labels = pd.Series(kmeans_model.labels_+1, index=l2_spectra.index)
    silhouette = silhouette_score(l2_spectra.values, kmeans_cluster_labels, metric='euclidean')

    # Find median usage for each gene across cluster
    median_spectra = l2_spectra.groupby(kmeans_cluster_labels).median()

    # Normalize median spectra to probability distributions.
    median_spectra = (median_spectra.T/median_spectra.sum(1)).T

    # reset index and col names 
    median_spectra.columns = range(len(median_spectra.columns))
    median_spectra = median_spectra.reset_index(drop=True)

    

    rf_usages = rf_usages.drop("bc_wells", axis = 1)
    rf_usages.columns = range(len(rf_usages.columns))
    rf_usages = rf_usages.reset_index(drop=True)


    # Compute prediction error as a frobenius norm
    rf_pred_norm_counts = rf_usages.dot(median_spectra)        
    if sp.issparse(norm_counts.X):
        prediction_error = ((norm_counts.X.todense() - rf_pred_norm_counts)**2).sum().sum()
    else:
        prediction_error = ((norm_counts.X - rf_pred_norm_counts)**2).sum().sum()    
        
    return pd.DataFrame([k, density_threshold, silhouette,  prediction_error],
            index = ['k', 'local_density_threshold', 'silhouette', 'prediction_error'],
            columns = ['stats'])


def load_stablity_error_data(output_directory, run_name, local_neighborhood_size=0.30, 
   density_threshold=2.0, components = [30, 50, 60, 80, 100, 200, 250, 300]):

    cnmf_obj = cnmf.cNMF(output_dir=output_directory, name=run_name)

    stats = []
    norm_counts = sc.read(cnmf_obj.paths['normalized_counts'])
    
    for k in components:
        rf_usages_path = '{output_directory}/{run_name}/{run_name}.usages.k_{k}.dt_{sel_thresh}.consensus.txt'.format(
                                                                                output_directory=output_directory,
                                                                                run_name = run_name,
                                                                                k=k,
                                                                                sel_thresh = str(sel_thresh).replace('.','_')

        stats.append(load_stablity_error_data_(k, cnmf_obj, norm_counts, rf_usages_path, density_threshold, local_neighborhood_size))

    arr = np.array(stats).squeeze()
    df = pd.DataFrame(arr, columns=['k', 'local_density_threshold', 'silhouette', 'prediction_error'])

    print("min stablity is", df['silhouette'].min())
    print("max stablity is", df['silhouette'].max())

    print("min error is",df['prediction_error'].min())
    print("max error is",df['prediction_error'].max())

    return df
'''