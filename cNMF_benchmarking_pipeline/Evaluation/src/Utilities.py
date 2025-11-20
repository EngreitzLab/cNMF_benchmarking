import os
import sys
import yaml
import logging
import mudata as mu
import pandas as pd
import cnmf
import scanpy as sc
import mygene



def rename_gene_mygene(mdata):

    mg = mygene.MyGeneInfo()
    gene_list = mdata['cNMF'].uns['var_names'].tolist()
    annotations = mg.querymany(gene_list, scopes='ensembl.gene', fields='symbol', species='human')

    # Process the results to create mapping
    gene_dict = {}
    
    for item in annotations:
        if 'symbol' in item:
            gene_dict[item['query']] = item['symbol']

    mdata['cNMF'].uns['var_names'] = [gene_dict.get(x, x) for x in mdata['cNMF'].uns['var_names']]

def rename_gene_dictionary(mdata, dictionary_file_path):

    # Convert mapping result to list before assignment
    df = pd.read_csv(dictionary_file_path, sep='\t')
    ensemble_to_gene = dict(zip(df['ensembl_id'], df['gene']))
    new_names = [ensemble_to_gene.get(x, x) for x in mdata['rna'].var_names]
    mdata['rna'].var_names = (new_names)

def check_evaluation_pipeline_format(mdata, prog_key = 'cNMF'):
    """
    Validate that the evaluation pipeline MuData object has the correct format.
    
    Checks for:
    - 'cNMguideF' modality exists
    - mdata[prog_key].uns has 'guide_names' column
    - mdata[prog_key].uns has 'guide_targetes' column
    - mdata[prog_key].obsm has 'guide_assignment' obsm that is dense
    - mdata[data_key]
    
    Parameters
    ----------
    mdata : mudata.MuData
        MuData object to validate
    prog_key : str
        index for the gene program anndata object (mdata[prog_key]) in the mudata object.   

    
    Returns
    -------
    bool
        True if all checks pass, False otherwise
    
    Raises
    ------
    Prints warnings for each validation failure
    """
    is_valid = True
    
    # Check if 'guide' modality exists
    if prog_key not in mdata.mod:
        print("WARNING: 'guide' modality not found in mdata")
        is_valid = False
    else:
        guide_adata = mdata[prog_key]
        
        # Check if 'guide_names' exists in uns
        if 'guide_names' not in guide_adata.uns:
            print("WARNING: 'guide_names' column not found in mdata['guide'].uns")
            is_valid = False
        else:
            print("guide_names' found in mdata['guide'].var")
        
        # Check if 'guide_targetes' exists in uns
        if 'guide_targets' not in guide_adata.uns:
            print("WARNING: 'guide_targets' column not found in mdata['guide'].uns")
            is_valid = False
        else:
            print("'guide_targets' found in mdata['guide'].var")
        
        # Check if 'guide_assignment' layer exists
        if 'guide_assignment' not in guide_adata.obsm:
            print("WARNING: 'guide_assignment' obsm not found in mdata['guide'].obsm")
            is_valid = False
        else:
            guide_assignment = guide_adata.obsm['guide_assignment']
            
            # Check if guide_assignment is sparse
            try:
                import scipy.sparse as sp
                is_sparse = sp.issparse(guide_assignment)
                
                if is_sparse:
                    print("WARNING: 'guide_assignment' is sparse. Converting to dense array...")
                    dense_array = guide_assignment.toarray()
                    print(f"'guide_assignment' converted to dense array (shape: {dense_array.shape})")
                else:
                    print(f"'guide_assignment' is already dense (shape: {guide_assignment.shape})")
            except Exception as e:
                print(f"WARNING: Error checking 'guide_assignment' sparsity: {e}")
                is_valid = False
    
    return is_valid

def check_gene_names(mdata, prog_key = 'cNMF', data_key = 'rna', gene_name_key = 'symbol' , categorical_key= 'batch', reference_gtf_path=None):
    """
    Validate gene names consistency across RNA modality and cNMF guide names,
    with optional reference GTF validation.
    
    Parameters
    ----------
    mdata : mudata.MuData
        MuData object containing 'rna' and 'cNMF' modalities
    prog_key : str
        index for the gene program anndata object (mdata[prog_key]) in the mudata object.   
    data_key : str
        index for the gene expression anndata object (mdata[prog_key]) in the mudata object.   
    categorical_key : str
        index for the gene name anndata object (mdata[data_key]) in the mudata object.
    reference_gtf_path : str, optional
        Path to reference GTF file (e.g., IGVF GENCODE v43 GTF).
        If provided, checks that gene names in mdata['rna'].var match the GTF.
        GTF should have a gene_id or gene_name column.
    
    Returns
    -------
    dict
        Dictionary containing:
        - 'is_valid' (bool): True if all checks pass
        - 'rna_vs_guide_mismatches' (list): Gene names in RNA but not in guide_names
        - 'missing_from_rna' (list): Gene names in guide_names but not in RNA
        - 'gtf_validation' (dict or None): GTF validation results if reference provided
    
    Raises
    ------
    Prints warnings for each validation failure
    """
    results = {
        'is_valid': True,
        'missing_from_rna': [],
        'gtf_validation': None
    }
    
    # check rna expression mod
    if data_key not in mdata.mod:
        print("WARNING: rna modality not found in mdata")
        results['is_valid'] = False
        return results
    else:
        print("mdata[data_key] modality found")

    rna_adata = mdata[data_key]

    # Check if categorical_key exists
    if categorical_key not in rna_adata.obs:
        print("WARNING: categorical_key column not found in mdata[data_key].obs")
        results['is_valid'] = False
        return results
    else:
        print("mdata[data_key].obs[categorical_key] batch key found")
    
    # check gene name in rna mod
    if gene_name_key not in rna_adata.var.columns:
        print("WARNING: gene_name_key column not found in mdata[gene_name_key].var")
        results['is_valid'] = False
        return results
    else:
        print("mdata[data_key].var[gene_name_key] gene names found")
    
    rna_gene_names = set(rna_adata.var[gene_name_key].values)
    print(f"Found {len(rna_gene_names)} gene names in mdata[data_key].var[gene_name_key']")
    
    # check program mod
    if prog_key not in mdata.mod:
        print("WARNING: program modality not found in mdata")
        results['is_valid'] = False
        return results
    else:
        print("mdata[prog_key] modality found")
    
    # check guide targets in program mod
    if 'guide_targets' not in mdata[prog_key].uns:
        print("WARNING: 'guide_targets' key not found in mdata[prog_key].uns")
        results['is_valid'] = False
        return results
    else:
        print("mdata[prog_key].uns['guide_targets'] guide names found")


    guide_targets = mdata[prog_key].uns['guide_targets']
    
    # Handle different data types for guide_names
    if isinstance(guide_targets, (list, tuple)):
        guide_gene_targets = set(guide_targets)
    elif isinstance(guide_targets, pd.Index):
        guide_gene_targets = set(guide_targets.values)
    elif isinstance(guide_targets, pd.Series):
        guide_gene_targets = set(guide_targets.values)
    else:
        guide_gene_targets = set(guide_targets)
    
    print(f"Found {len(guide_gene_targets)} gene names in mdata[prog_key].uns['guide_targets']")
    
    # Compare RNA vs guide names
    guide_not_in_rna = guide_gene_targets - rna_gene_names

    if guide_not_in_rna:
        print(f"\n WARNING: {len(guide_not_in_rna)} gene names in guide_names but NOT in adata.var for RNA. \n \
            This might be caused by mismatch between gene symbol vs Ensembl ID or some perturbed genes are not expressed in the dataset. ")
        print(f" Examples: {list(guide_not_in_rna)[:5]}")
        results['missing_from_rna'] = sorted(list(guide_not_in_rna))
    else:
        print("All guide_names found in RNA gene names")
    
    # compare GTF
    if reference_gtf_path is None:
        print("No reference GTF provided (optional check skipped)")
        print("To validate against a reference, provide reference_gtf_path parameter")
        results['gtf_validation'] = None
    else:
        gtf_validation = _validate_against_reference_gtf(
            rna_gene_names, 
            reference_gtf_path
        )
        results['gtf_validation'] = gtf_validation
        
        if not gtf_validation['is_valid']:
            results['is_valid'] = False
            print("test")

    
    return results

def _validate_against_reference_gtf(rna_gene_names, gtf_path):
    """
    Validate RNA gene names against a reference GTF file. Gene names must be in gene symbols, not in gene ids 
    
    Parameters
    ----------
    rna_gene_names : set
        Set of gene names from RNA data
    gtf_path : str
        Path to GTF file
    
    Returns
    -------
    dict
        Validation results with keys:
        - 'is_valid' (bool)
        - 'gtf_genes' (set): Gene names from GTF
        - 'not_in_gtf' (list): RNA genes not in GTF
    """
    validation = {
        'is_valid': True,
        'gtf_genes': set(),
        'not_in_gtf': [],
     }
    
    try:
        # Read GTF file - handling common gene name columns
        gtf_df = pd.read_csv(
            gtf_path, 
            sep='\t', 
            comment='#', 
            header=None,
            dtype={0: str, 1: str, 2: str, 3: int, 4: int, 5: str, 6: str, 7: str, 8: str}
        )
        
        # Extract gene names from GTF attributes (column 8)
        extracted = gtf_df[8].str.extract(r'gene_name "([^"]+)"')[0]

        if extracted.notna().sum() == 0:
            print("gene_name not found in any rows")
        else:
            print(f"Successfully extracted gene_name from {extracted.notna().sum()} rows")
            extracted = set(extracted)
            validation['gtf_genes'] = extracted
                
        
        # Compare with RNA names
        rna_not_in_gtf = rna_gene_names - extracted
        
        if rna_not_in_gtf:
            print(f"\nWARNING: {len(rna_not_in_gtf)} adata.var RNA NOT found in reference GTF. \
                \n This might be caused by the adata.var RNA being Ensembl ID instead of gene symbols.")
            print(f"   Examples: {list(rna_not_in_gtf)[:5]} \n")
            validation['not_in_gtf'] = sorted(list(rna_not_in_gtf))
        else:
            print("All RNA genes found in reference GTF")
    

    except FileNotFoundError:
        print(f"ERROR: Reference GTF file not found at {gtf_path}")
        validation['is_valid'] = False
    except Exception as e:
        print(f"ERROR: Failed to read reference GTF: {e}")
        validation['is_valid'] = False
    
    return validation
      