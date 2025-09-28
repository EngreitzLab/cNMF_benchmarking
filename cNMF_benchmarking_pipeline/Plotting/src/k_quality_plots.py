import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.sparse import coo_matrix
import scipy.sparse as sp
import os
import seaborn as sns
import xarray as xr
from matplotlib.backends.backend_pdf import PdfPages



# make a program dotplots by days
def programs_dotplots(mdata_path, program_loading_path, save_path = None, save_name = None):

    # load data
    df = pd.read_csv(program_loading_path, sep = '\t', index_col = 0 )
    adata = mu.read_h5mu(mdata_path)

    # make anndata from program loadings
    adata_new = ad.AnnData(X=df.values)
    adata_new.obs["sample"] = adata['rna'].obs["sample"].values


    if save_name is None:
        save_name = "Program Loadings by Days"
        
    gene_list = adata_new.var_names.tolist()
    dp = sc.pl.dotplot(adata_new, gene_list, groupby="sample",swap_axes=True, dendrogram = True, show = False,
                       return_fig = True, title = save_name)


    # make fig object for pdf plotting
    (dp.make_figure())
    fig = dp.fig

    # save fig 
    if save_name and save_path:
        fig.savefig(f"{save_path}/{save_name}.png", bbox_inches='tight')

    return fig

# compute correlation coeffcient between two program x gene matrices
def program_corr(matrix1,matrix2):
    
    if matrix2.shape != matrix2.shape:
        print("dim is different")
        return 
    
    num_columns = matrix1.shape[0]
    column_correlations = []
    
    for i in range(num_columns):
        for j in range(num_columns):
            
            # Calculate the correlation coefficient between the i-th row of matrix1 and j-th row of matrix2
            correlation = np.corrcoef(matrix1.iloc[i,:], matrix2.iloc[j,:])[0, 1]
            column_correlations.append(correlation)
            
    df = pd.DataFrame(np.array(column_correlations).reshape(num_columns, num_columns), 
                     columns = matrix2.index,
                     index = matrix1.index)

    return df

# compute euclidea distance btween two program x gene matrices
def program_euclidean(matrix1, matrix2):
    
    if matrix1.shape != matrix2.shape:
        print("dim is different")
        return 
    
    num_rows = matrix1.shape[0]
    euclidean_distances = []
    
    for i in range(num_rows):
        for j in range(num_rows):
            
            # Calculate the Euclidean distance between the i-th row of matrix1 and j-th row of matrix2
            distance = np.sqrt(np.sum((matrix1.iloc[i,:] - matrix2.iloc[j,:])**2))
            euclidean_distances.append(distance)
            
    df = pd.DataFrame(np.array(euclidean_distances).reshape(num_rows, num_rows), 
                     columns=matrix2.index,
                     index=matrix1.index)

    return df

# compute top x overlapped genes in percentages btween two program x gene matrices
def top_genes_overlap(matrix1, matrix2, percentage = True, gene_num = 300):
    
    # check dim 
    if matrix1.shape != matrix2.shape:
        print("Different dim")
        return 
    
    # find out top x genes in each k
    top_genes_1 = matrix1.apply(lambda row: row.nlargest(gene_num).index.tolist(), axis=1)
    top_genes_2 = matrix2.apply(lambda row: row.nlargest(gene_num).index.tolist(), axis=1)

    
    n_k = (top_genes_1.shape[0])
    overlap = np.zeros((n_k, n_k), dtype=float)
    gene_shared_list = pd.Series(dtype=object)

    # generate overlap matrix 
    for i in range(n_k):
        for j in range(n_k):
            if percentage:
                s = len(set(top_genes_1.iloc[i]) & set(top_genes_2.iloc[j]))/gene_num
            else: 
                s = len(set(top_genes_1.iloc[i]) & set(top_genes_2.iloc[j])) 

            # compose a shared gene matrix
            #gene_shared = list(set(top_genes_1.iloc[i]) & set(top_genes_2.iloc[j]))
            #name = f"{top_genes_1.index[i]} VS {top_genes_2.index[j]}"
            #gene_shared_list.at[name] = gene_shared 
            
            overlap[i, j] = int(s)

    
    overlap_df = pd.DataFrame(overlap,
                              index=matrix1.index,
                              columns=matrix2.index)
    
    return overlap_df #, top_genes_1, top_genes_2, gene_shared_list

# for each row, sort the values on the columns, max value is on diagnol
def sort_corr_matrix(Matrix):

    n = Matrix.shape[0]
    Matrix_reordered = Matrix.copy()

    for i in range(n):
        # find index of row maximum
        j = Matrix.iloc[i].values.argmax()

        # move that column to diagonal position
        Matrix_reordered.iloc[i] = np.roll(Matrix.iloc[i].values, i - j)
    
    return Matrix_reordered

# given a list of shared genes across two programs, plot bar plot
def max_gene_values_barplot(data_list, title="Maxinum shared genes between sk-cd and torch-halsvar", x_label="torch-halsvar Program", y_label="Shared Genes" , figsize = (3,5)):

    # Create x-axis positions (indices)
    x_positions = range(len(data_list))
    
    # Create the bar plot
    plt.figure(figsize=figsize)
    bars = plt.bar(x_positions, data_list, color='skyblue', alpha=0.7)
    
    # Customize the plot
    plt.title(title, fontweight='bold')
    plt.xlabel(x_label, fontsize=12)
    plt.ylabel(y_label, fontsize=12)
    plt.grid(axis='y', alpha=0.3)
    
    
    # Set x-axis labels to show the actual numbers from the list
    #plt.xticks(x_positions, [f'{i}' for i in range(len(data_list))])
    
    plt.tight_layout()
    plt.show()


# For graphing clustermap given corr/eu/overlap matrix
def graph_cluster(matrix, method = 'single', save = False, save_folder_name = None , save_file_name = None , show = True, figsize = (4,4) ):
 
    # label color for each run 
    #palette = sns.color_palette("tab10", n_colors=matrix.columns.nunique())
    #lut = dict(zip(matrix.columns.unique(), palette))
    #row_colors = matrix.columns.map(lut)
    
    sns.set_theme(style="white")               # optional aesthetics
    g = sns.clustermap(
            matrix,
            cmap="vlag",                       # diverging palette centred at 0
            #linewidths=.5,                     # grid lines
            center=0,                          # keep 0 in the middle of the colour range
            metric="euclidean",                # distance metric for clustering
            method= method,                    # linkage method
            figsize=figsize,                  # size in inches
            #row_colors=row_colors,            # color axis 
            #col_colors = row_colors,          # color axis 
            row_cluster=True,
            col_cluster=True
    )

    # drop names, too many names now
    # Remove axis labels and tick labels
    #g.ax_heatmap.set_xlabel('')
    #g.ax_heatmap.set_ylabel('')
    g.ax_heatmap.set_xticklabels([])
    g.ax_heatmap.set_yticklabels([])

    g.fig.suptitle(save_file_name)


    # Remove dendrograms
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    # Remove tick marks/lines
    g.ax_heatmap.tick_params(left=False, right= False, bottom=False)

    if show:
        plt.show()

    if save:
        g.savefig(f"{save_folder_name}/{save_file_name}.png", bbox_inches='tight', dpi=150)
    
    return g

'''# graph heatmap base on the clustermap -> TODOs: need more thinking in this
def graph_heatmap(g, r, c, folder_name, file_name, num_gene = 300, sorted = False):
    # g = clustermap
    # r,c dimension of calculating averages

    mat = g.data2d.to_numpy() * num_gene

    assert mat.shape[0] % r == 0 and mat.shape[1] % c == 0, \
           "Rows/cols must divide evenly by block size."

    n_row_blocks = mat.shape[0] // r
    n_col_blocks = mat.shape[1] // c

    blocks = (mat.reshape(n_row_blocks, r, n_col_blocks, c)
            .swapaxes(1, 2)              
            .reshape(-1, r, c))

    block_means = (blocks.mean(axis=(1, 2))).astype(int)
    plt.figure(figsize=(12, 8))
    sns.heatmap(block_means.reshape(10,10),annot=True, cmap='inferno_r',fmt='d')        
    plt.title("Heatmap for matching programs " + file_name)

    if folder_name and file_name:
        g.savefig(f"{folder_name}/{file_name}.png")
 
    plt.show()

    # Sorted heatmap
    if sorted: 
        matrix = block_means.reshape(10,10).tolist()

        for i in range(len(matrix)):
            max_index = matrix[i].index(max(matrix[i]))
            # Swap max element with the diagonal element
            matrix[i][i], matrix[i][max_index] = matrix[i][max_index], matrix[i][i]

        plt.figure(figsize=(12, 8))
        sns.heatmap(np.array(matrix).reshape(r,c),annot=True, cmap='inferno_r',fmt='d')        
        plt.title("Sorted Heatmap for matching programs " + file_name)

        if folder_name and file_name:
            g.savefig(f"{folder_name}/{file_name}_sorted.png")
 
        plt.show()
'''
# graph pdf for 3 clustermap
def graph_pdf_clustermap(cor, distance, overlap, save_path, filename):

    with PdfPages(f"{save_path}/{filename}") as pdf:
        

        g1 = graph_cluster(cor, save = True, save_folder_name = save_path , save_file_name = "ClusterMap_for_correlation")
        g2 = graph_cluster(distance, save = True, save_folder_name = save_path , save_file_name = "ClusterMap_for_Euclidean_distance")
        g3 = graph_cluster(overlap, save = True, save_folder_name = save_path , save_file_name = "ClusterMap_for_top_300_genes")

        fig, axes = plt.subplots(1, 3, figsize=(30, 10))
        fig.suptitle(filename, 
                    fontsize=24, fontweight='bold', y=0.95)

        img1 = plt.imread(f"{save_path}/ClusterMap_for_correlation.png")
        axes[0].imshow(img1)
        axes[0].axis('off')
        
        img2 = plt.imread(f"{save_path}/ClusterMap_for_Euclidean_distance.png") 
        axes[1].imshow(img2)
        axes[1].axis('off')

        img3 = plt.imread(f"{save_path}/ClusterMap_for_top_300_genes.png") 
        axes[2].imshow(img3)
        axes[2].axis('off')
        
        plt.tight_layout()
        pdf.savefig(fig,bbox_inches='tight')
        plt.close()
        
        plt.close(g1.fig)
        plt.close(g2.fig)
        plt.close(g3.fig)
        os.remove(f"{save_path}/ClusterMap_for_correlation.png")
        os.remove(f"{save_path}/ClusterMap_for_Euclidean_distance.png")
        os.remove(f"{save_path}/ClusterMap_for_top_300_genes.png")
                              