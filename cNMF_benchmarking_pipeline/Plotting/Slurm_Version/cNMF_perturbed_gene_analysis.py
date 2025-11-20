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
                         rename_list_gene_dictionary, plot_umap_per_gene_guide, process_single_gene, parallel_gene_processing,\
                         compute_gene_correlation_matrix, compute_gene_waterfall_cor


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--mdata_path', type=str, required=True)  
    parser.add_argument('--perturb_path', type=str, required=True) # partial path for each day
    parser.add_argument('--file_to_dictionary',  type=str, default = None)  
    parser.add_argument('--groupby', type=str, default = "sample")  
    parser.add_argument('--tagert_col_name', type=str, default = "target_name")  
    parser.add_argument('--plot_col_name', type=str, default = "program_name")  
    parser.add_argument('--log2fc_col', type=str, default = "log2FC") 
    parser.add_argument('--top_num',  type=int, default=5)
    parser.add_argument('--p_value',  type=float, default=0.05)
    parser.add_argument('--down_thred_log',  type=float, default=-0.05)
    parser.add_argument('--up_thred_log',  type=float, default=0.05)
    parser.add_argument('--pdf_save_path',  type=str, required=True)
    parser.add_argument('--samples', nargs='+', default=['D0', 'sample_D1', 'sample_D2', 'sample_D3'])
    parser.add_argument('--square_plots',  action="store_true")  
    parser.add_argument('--figsize', type=float, nargs=2, default=(35, 35))
    parser.add_argument('--show',  action="store_true")  
    parser.add_argument('--PDF',  action="store_true")  
    parser.add_argument('--n_processes', type=int, default = -1) 



    args = parser.parse_args()


    # save comfigs used         
    args_dict = vars(args)
    with open(f'{args.pdf_save_path}/config.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)
        

    #read mdata
    mdata = mu.read_h5mu(args.mdata_path)

    # compute corr once
    correlation_matrix = compute_gene_correlation_matrix(mdata, file_to_dictionary = args.file_to_dictionary)


    # found detected perturbed gene
    perturbed_gene = np.unique(mdata['cNMF'].uns["guide_targets"])
    gene_list = rename_list_gene_dictionary(mdata['rna'].var_names.tolist(), args.file_to_dictionary) # convert gene id to geene name
    perturbed_gene_found = list(set(gene_list) & set(perturbed_gene.tolist()))

    # sort list by alphabetical order 
    perturbed_gene_found = sorted(perturbed_gene_found)
    print(f"there are {len(perturbed_gene_found)} perturbed gene found")


   
    # Graph all pdf 
    for gene in perturbed_gene_found:

        create_comprehensive_plot(
            mdata = mdata,
            perturb_path = args.perturb_path,
            Target_Gene = gene,
            correlation_matrix = correlation_matrix,
            file_to_dictionary = args.file_to_dictionary,
            groupby=args.groupby,
            tagert_col_name=args.tagert_col_name,
            plot_col_name=args.plot_col_name,
            log2fc_col=args.log2fc_col,
            top_num=args.top_num,
            p_value=args.p_value,
            down_thred_log=args.down_thred_log,
            up_thred_log=args.up_thred_log,
            save_path = args.pdf_save_path,
            save_name = gene,
            square_plots=args.square_plots,
            figsize = args.figsize,
            show=args.show,
            PDF = True,
            samples= ['D0', 'sample_D1', 'sample_D2', 'sample_D3'],
        )
        
    '''
    print("Starting parallel gene processing...")                   

    try: 
        # Graph all pdf 
        result = parallel_gene_processing( 
            mdata = mdata,
            perturb_path = args.perturb_path,
            perturbed_gene_list = perturbed_gene_found,
            file_to_dictionary = args.file_to_dictionary,
            groupby=args.groupby,
            tagert_col_name=args.tagert_col_name,
            plot_col_name=args.plot_col_name,
            log2fc_col=args.log2fc_col,
            top_num=args.top_num,
            p_value=args.p_value,
            down_thred_log=args.down_thred_log,
            up_thred_log=args.up_thred_log,
            pdf_save_path=args.pdf_save_path,
            samples=args.samples,
            square_plots=args.square_plots,
            save_path = args.pdf_save_path,
            figsize = args.figsize,
            show=args.show,
            PDF=args.PDF,
            n_processes = args.n_processes
        )
        print(f"Parallel processing completed successfully. Results: {len(result) if result else 'None'}")    

    except Exception as e:                                                                                 
        print(f"ERROR in parallel_gene_processing: {e}")  
    
    '''

    # merge pdf 
    if args.PDF:
        merge_pdfs_in_folder(args.pdf_save_path)
    else:
        merge_svgs_to_pdf(args.pdf_save_path)

