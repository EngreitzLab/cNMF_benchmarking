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
import mygene
import shutil
from pathlib import Path
from tqdm.auto import tqdm

# running multiple k consensus at the same time 
def run_cnmf_consensus(cnmf_obj=None, output_dir=None, name=None, 
                       components=[7,8,9,10], density_thresholds=[0.01, 0.05, 2.0]):

    if cnmf_obj is None:
        cnmf_obj = init_cnmf_obj(output_dir=output_dir, name=name)

    for k in tqdm(components, desc='Running cNMF'):
        for thresh in density_thresholds:
            cnmf_obj.consensus(k=k, density_threshold=thresh, show_clustering=True)


# compile results into correct format for downstream Evaluation pipeline
def compile_results(output_directory, run_name, sel_thresh = [2.0], components = [30, 50, 60, 80, 100, 200, 250, 300] ):
       
 for i in sel_thresh:
    for k in K:

        scores = pd.read_csv('{output_directory}/{run_name}/{run_name}.usages.k_{k}.dt_{sel_thresh}.consensus.txt'.format(
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k,
                                                                                        sel_thresh = str(i).replace('.','_')),
                                                                                        sep='\t', index_col=0)

        loadings = pd.read_csv('{output_directory}/{run_name}/{run_name}.spectra.k_{k}.dt_{sel_thresh}.consensus.txt'.format(
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k,
                                                                                        sel_thresh = str(i).replace('.','_')),
                                                                                        sep='\t', index_col=0)
        

        os.makedirs((f'{output_directory}/{run_name}/loading'), exist_ok=True)


        scores.to_csv('{output_directory}/{run_name}/loading/cNMF_scores_{k}_{sel_thresh}.txt'.format(
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k,
                                                                                        sel_thresh = i), sep='\t')
        loadings.T.to_csv('{output_directory}/{run_name}/loading/cNMF_loadings_{k}_{sel_thresh}.txt'.format(     
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k,
                                                                                        sel_thresh = i), sep='\t')

        adata_ = anndata.read_h5ad('{output_directory}/{run_name}/cnmf_tmp/{run_name}.tpm.h5ad'.format(
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k ))
        adata_.var_names_make_unique()
        adata_.obs_names_make_unique()

        prog_data = anndata.AnnData(X=scores.values, obs=adata_.obs)
        prog_data.varm['loadings'] = loadings.values
        prog_data.uns['var_names'] = loadings.columns.values


        # Make adata
        os.makedirs((f'{output_directory}/{run_name}/prog_data'), exist_ok=True)
        prog_data.write(f'{output_directory}/{run_name}/prog_data/NMF_{k}_{sel_thresh}.h5ad'.format(
                                                                                output_directory=output_directory,
                                                                                run_name = run_name,
                                                                                k=k,
                                                                                sel_thresh = str(i).replace('.','_')))

        # Make mdata
        mdata = muon.MuData({'rna': adata_, 'cNMF': prog_data})

        os.makedirs((f'{output_directory}/{run_name}/adata'), exist_ok=True)
        mdata.write(f'{output_directory}/{run_name}/adata/cNMF_{k}_{sel_thresh}.h5mu'.format(
                                                                                output_directory=output_directory,
                                                                                run_name = run_name,
                                                                                k=k,
                                                                                sel_thresh = str(i).replace('.','_')))
    


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


# combine parallel inferred cNMF component results into one run
def rename_and_move_files_NMF(file_name_input, file_name_output, source_folder, destination_folder, len = 10):

    # Create destination folder if it doesn't exist
    Path(destination_folder).mkdir(parents=True, exist_ok=True)
    
    processed_files = []
    
    # Process files for i 
    for i in range(len):
        source_filename = f"{file_name_input}_{i}.df.npz"
        source_path = os.path.join(source_folder, source_filename)
        
        # Check if source file exists
        if os.path.exists(source_path):
            
            # Create new filename
            new_filename = f"{file_name_output}_{i}.df.npz"
            destination_path = os.path.join(destination_folder, new_filename)
            
            try:
                # Copy and rename file to destination
                shutil.copy2(source_path, destination_path)
                print(f"Successfully processed: {source_filename} -> {new_filename}")
                processed_files.append((source_filename, new_filename))
                
            except Exception as e:
                print(f"Error processing {source_filename}: {e}")
        else:
            print(f"File not found: {source_filename}")
    
    return processed_files


# combine all Ks parallel run
def rename_all_NMF(file_name_input, file_name_output, source_folder, destination_folder, len = 10, components = [30, 50, 60, 80, 100, 200, 250, 300]):

    for k in components:

        file_name_input_new = f"{file_name_input}_{k}.spectra.k_{k}.iter"
        file_name_output_new = f"{file_name_output}.spectra.k_{k}.iter"

        source_folder_new = f"{source_folder}_{k}/cnmf_tmp"

        rename_and_move_files_NMF(file_name_input_new, file_name_output_new, source_folder_new, destination_folder)


# combine 2 cNMF results together for one K
def rename_and_move_files(k, file_name_input, file_name_output, source_dir, dest_dir, len = 10, second = False):
    
    # Create destination folder if it doesn't exist
    Path(dest_dir).mkdir(parents=True, exist_ok=True)
    
    processed_files = []
    
    # Process files for i 
    for i in range(len):

        source_file = f"{file_name_input}.spectra.k_{k}.iter_{i}.df.npz"
        source_path = os.path.join(source_dir, source_file )

        # Check if source file exists
        if os.path.exists(source_dir): 

            if second: 
                i = i + 10
            
            # Create new filename
            new_filename = f"{file_name_output}.spectra.k_{k}.iter_{i}.df.npz"
            destination_path = os.path.join(dest_dir, new_filename)
            
            try:
                # Copy and rename file to destination
                shutil.copy2(source_path, destination_path)
                print(f"Successfully processed: {source_file} -> {new_filename}")
                processed_files.append((source_file , new_filename))
                
            except Exception as e:
                print(f"Error processing {file_name_input}: {e}")
        else:
            print(f"File not found: {file_name_input}")
            print(source_path)
    
    return processed_files


def rename_all(file_name_input, file_name_output, source_dir, dest_dir, components = [30, 50, 60, 80, 100, 200, 250, 300], second = False):

    for k in components:

        data = rename_and_move_files(k, file_name_input, file_name_output, source_dir, dest_dir, second = second)
