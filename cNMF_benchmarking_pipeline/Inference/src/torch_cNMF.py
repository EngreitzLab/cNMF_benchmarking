import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image
from matplotlib import gridspec
import scanpy as sc
from pathlib import Path
import scanpy as sc
import anndata
import muon 
import cnmf
import os
import sys
import pandas as pd
import mygene


from tqdm.auto import tqdm
# running multiple k consensus at the same time 
def run_cnmf_consensus(cnmf_obj=None, output_dir=None, name=None, 
                       components=[7,8,9,10], density_thresholds=[0.01, 0.05, 2.0]):

    if cnmf_obj is None:
        cnmf_obj = init_cnmf_obj(output_dir=output_dir, name=name)

    for k in tqdm(components, desc='Running cNMF'):
        for thresh in density_thresholds:
            cnmf_obj.consensus(k=k, density_threshold=thresh, show_clustering=True)

# given a df from cNMF, return top 300 gene for each program in df
def get_top_indices_fast(df, gene_num = 300):
 
    # Get column names
    col_names = df.columns.values
    
    # Use argsort to get indices of top 300 values per row
    # argsort sorts in ascending order, so we use [:, -300:] and reverse
    top_indices = np.argsort(df.values, axis=1)[:, -gene_num:][:, ::-1]
    
    # Map indices to column names
    top_col_names = col_names[top_indices]
    
    # Create result DataFrame
    result_df = pd.DataFrame(
        top_col_names,
        index=df.index,
        columns=[f'top_{i+1}' for i in range(gene_num)]
    )

    result_df.index = [f'Program_{i}' for i in range(1,len(result_df)+1)]
    
    return result_df
    
# annotate genes in excel given a df with rows for each program, cols for genes    
def annotate_genes_to_excel(df, output_file='gene_annotations.xlsx'):
    
    # Initialize MyGene
    mg = mygene.MyGeneInfo()
    
    # Dictionary to store results for each column
    all_annotations = {}
    
    # Process each column
    for row_idx in df.index:
        # Get unique genes from column (remove NaN)
        genes = df.loc[row_idx].dropna().unique().tolist()
        
        if len(genes) == 0:
            print(f"Column '{row_idx}': No genes found")
            continue
        
        print(f"Annotating column '{row_idx}': {len(genes)} genes...")
        
        # Query MyGene for annotations
        results = mg.querymany(
            genes, 
            scopes='symbol,alias,ensembl.gene',  # Multiple search scopes
            fields='symbol,name,entrezgene,summary,type_of_gene',
            species='human',
            returnall=True
        )
        
        # Process results
        annotation_list = []
        for query_gene, gene_info in zip(genes, results['out']):
            # Handle cases where gene is not found or multiple matches
            if 'notfound' in gene_info and gene_info['notfound']:
                annotation_list.append({
                    'Input_Gene': query_gene,
                    'Gene_Symbol': 'NOT FOUND',
                    'Gene_Name': 'NOT FOUND',
                    'Entrez_ID': '',
                    'Type': '',
                    'Summary': ''
                })
            else:
                annotation_list.append({
                    'Input_Gene': query_gene,
                    'Gene_Symbol': gene_info.get('symbol', query_gene),
                    'Gene_Name': gene_info.get('name', ''),
                    'Entrez_ID': gene_info.get('entrezgene', ''),
                    'Type': gene_info.get('type_of_gene', ''),
                    'Summary': gene_info.get('summary', '')
                })
        
        # Create DataFrame for this column
        all_annotations[row_idx] = pd.DataFrame(annotation_list)

    # Export to Excel with multiple sheets
    with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
        for row_idx, annotations_df in all_annotations.items():
            sheet_name = str(row_idx)
            # Truncate sheet name if too long (Excel limit is 31 chars)
            sheet_name = sheet_name[:31] if len(sheet_name) > 31 else sheet_name
            annotations_df.to_excel(writer, sheet_name=sheet_name, index=False)
    
    return all_annotations
