import sys
import muon as mu 
import numpy as np
import pandas as pd
import argparse
import yaml
import os

# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from Plotting.src import plot_umap_per_gene, plot_top_program_per_gene, perturbed_gene_dotplot,\
                         plot_log2FC, plot_volcano, programs_dotplot, analyze_correlations, \
                         create_gene_correlation_waterfall, \
                         convert_with_mygene, convert_adata_with_mygene, read_npz, \
                         merge_pdfs_in_folder, merge_svgs_to_pdf, create_comprehensive_plot, rename_adata_gene_dictionary, \
                         rename_list_gene_dictionary, plot_umap_per_gene_guide, process_single_gene, parallel_gene_processing


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--mdata_path', type=str, required=True)  
    parser.add_argument('--perturb_path', type=str, required=True) # partial path for each day
    parser.add_argument('--top_program',  type=int, default=5)
    parser.add_argument('--p_value',  type=float, default=0.05)
    parser.add_argument('--pdf_save_path',  type=str, required=True)
    parser.add_argument('--PDF',  action="store_true")  
    parser.add_argument('--file_to_dictionary',  type=str, default = None)  
    parser.add_argument('--n_processes', type=int, default = -1)  

    args = parser.parse_args()


    # save comfigs used         
    args_dict = vars(args)
    with open(f'{args.pdf_save_path}/config.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)
        

    #read mdata
    mdata = mu.read_h5mu(args.mdata_path)

    # found detected perturbed gene
    perturbed_gene = np.unique(mdata['cNMF'].uns["guide_targets"])
    gene_list = rename_list_gene_dictionary(mdata['rna'].var_names.tolist(), args.file_to_dictionary) # convert gene id to geene name
    perturbed_gene_found = list(set(gene_list) & set(perturbed_gene.tolist()))

    # sort list by alphabetical order 
    perturbed_gene_found = sorted(perturbed_gene_found)
    print(f"there are {len(perturbed_gene_found)} perturbed gene found")
    
    print("Starting parallel gene processing...")                   

    try: 
        # Graph all pdf 
        result = parallel_gene_processing( 
            mdata = mdata,
            perturb_path = args.perturb_path,
            perturbed_gene_list = perturbed_gene_found,
            file_to_dictionary = args.file_to_dictionary,
            save_path = args.pdf_save_path,
            figsize = (30, 30),
            show=False,
            PDF=args.PDF,
            n_processes = args.n_processes
        )
        print(f"Parallel processing completed successfully. Results: {len(result) if result else 'None'}")    

    except Exception as e:                                                                                 
        print(f"ERROR in parallel_gene_processing: {e}")  
    
    # merge pdf 
    if args.PDF:
        merge_pdfs_in_folder(args.pdf_save_path)
    else:
        merge_svgs_to_pdf(args.pdf_save_path)

