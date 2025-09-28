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
import os
from PyPDF2 import PdfMerger
import glob


# replace EnsemblID by gene name given the dataframe with EnsemblID as the index
def convert_with_mygene(dataframe, species='human', index = True):

    mg = mygene.MyGeneInfo()
    
    if index:
        # Query multiple genes at once
        result = mg.querymany(dataframe.index, 
                            scopes='ensembl.gene', 
                            fields='symbol,name', 
                            species='human')
    else:
         # Query multiple genes at once
        result = mg.querymany(dataframe.columns, 
                            scopes='ensembl.gene', 
                            fields='symbol,name', 
                            species='human')

    
    # Create mapping dictionary
    mapping = {}
    for item in result:
        if 'symbol' in item and 'query' in item:
            mapping[item['query']] = item['symbol']

        elif 'query' in item:
            mapping[item['query']] = item['query']  # Keep original if no symbol

    if index: 
        new_dataframe = dataframe.rename(index=mapping)
    else:
        new_dataframe = dataframe.rename(columns=mapping)
    
    return new_dataframe

# same, input is adata
def convert_adata_with_mygene(adata, species='human'):

    mg = mygene.MyGeneInfo()
    gene_list = adata.var_names.tolist()
    annotations = mg.querymany(gene_list, scopes='ensembl.gene', 
                            fields='symbol', species='human')

    # Process the results to create mapping
    gene_dict = {}
    for item in annotations:
        if 'symbol' in item and 'query' in item:
            gene_dict[item['query']] = item['symbol']

        elif 'query' in item:
            gene_dict[item['query']] = item['query']  # Keep original if no symbol

    adata_new = adata.copy()


    adata_new.var['gene_name'] = [gene_dict.get(x, x) for x in adata_new.var_names]
    adata_new.var_names = adata_new.var['gene_name']

    return adata_new

# read cNMF programs
def read_npz(path):

    # Load the NPZ file with pickle enabled
    npz_data = np.load(path, allow_pickle=True)

    # Reconstruct the DataFrame
    df = pd.DataFrame(
        data=npz_data['data'],
        index=npz_data['index'],
        columns=npz_data['columns']
    )
     
    return df


# merge all PDFs into one, save pdf in the same folder_path
def merge_pdfs_in_folder(folder_path, output_filename="merged_perturbed_gene_QC.pdf"):

    
    # Create PdfMerger object
    merger = PdfMerger()
    
    # Get all PDF files in the folder
    pdf_files = glob.glob(os.path.join(folder_path, "*.pdf"))
    print(f"Found {len(pdf_files)} PDF files")
    
    # Merge each PDF
    for pdf_file in pdf_files:
        try:
            merger.append(pdf_file)
        except Exception as e:
            print(f"Error processing {pdf_file}: {str(e)}")
            continue
    
    # Write the merged PDF to output file
    output_path = os.path.join(folder_path, output_filename)
    
    try:
        with open(output_path, 'wb') as output_file:
            merger.write(output_file)
        merger.close()
        
    except Exception as e:
        print(f"Error saving merged PDF: {str(e)}")

def merge_svgs_to_pdf(folder_path, output_filename="merged_perturbed_gene_QC.pdf"):

    svg_files = glob.glob(os.path.join(folder_path, "*.svg"))
    print(f"Found {len(svg_files)} SVG files")
    
    merger = PdfMerger()
    temp_pdfs = []
    
    for svg_file in svg_files:
        try:
            # Convert SVG to PDF
            drawing = svg2rlg(svg_file)
            temp_pdf = svg_file.replace('.svg', '_temp.pdf')
            renderPDF.drawToFile(drawing, temp_pdf)
            temp_pdfs.append(temp_pdf)
            merger.append(temp_pdf)
        except Exception as e:
            print(f"Error processing {svg_file}: {str(e)}")
    
    # Save merged PDF
    output_path = os.path.join(folder_path, output_filename)
    with open(output_path, 'wb') as f:
        merger.write(f)
    merger.close()
    
    # Clean up temp files
    for temp_pdf in temp_pdfs:
        os.remove(temp_pdf)
    
    print(f"PDF created with {len(svg_files)} pages: {output_path}")