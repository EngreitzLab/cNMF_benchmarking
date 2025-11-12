import sys
import muon as mu 
import numpy as np
import pandas as pd
import argparse
import yaml


# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from Plotting.src import merge_pdfs_in_folder, merge_svgs_to_pdf

from Plotting.src import plot_umap_per_program, plot_top_gene_per_program, top_GO_per_program, compute_program_correlation_matrix,\
                              analyze_program_correlations, plot_violin, plot_program_log2FC, plot_program_heatmap, plot_program_volcano, \
                              perturbed_gene_dotplot, compute_program_waterfall_cor, create_program_correlation_waterfall, create_comprehensive_program_plot
                              
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--mdata_path', type=str, required=True)  
    parser.add_argument('--perturb_path_base', type=str, required=True) # partial path for each day
    parser.add_argument('--file_to_dictionary',  type=str, default = None) 
    parser.add_argument('--GO_path', type=str, required=True)  
    parser.add_argument('--groupby', type=str, default = "sample")  
    parser.add_argument('--tagert_col_name', type=str, default = "target_name")  
    parser.add_argument('--plot_col_name', type=str, default = "program_name")  
    parser.add_argument('--top_program',  type=int, default=5)    
    parser.add_argument('--top_enrichned_term',  type=int, default=10)
    parser.add_argument('--log2fc_col', type=str, default = "log2FC") 
    parser.add_argument('--p_value',  type=float, default=0.05)
    parser.add_argument('--down_thred_log',  type=float, default=-0.05)
    parser.add_argument('--up_thred_log',  type=float, default=0.05)
    parser.add_argument('--pdf_save_path',  type=str, required=True)
    parser.add_argument('--samples', nargs='+', default=['D0', 'sample_D1', 'sample_D2', 'sample_D3'])
    parser.add_argument('--square_plots',  action="store_true")  
    parser.add_argument('--figsize', type=float, nargs=2, default=(35, 35))
    parser.add_argument('--show',  action="store_true")  
    parser.add_argument('--PDF',  action="store_true")  


    args = parser.parse_args()


    # save comfigs used         
    args_dict = vars(args)
    with open(f'{args.pdf_save_path}/config.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)
        

    #read mdata
    mdata = mu.read_h5mu(args.mdata_path)
    program_len = len(mdata['cNMF'].var) # find out list of programs to process 
    print(f"there are {program_len} Program found")



    # compute correlations
    Sample = ["D0", "sample_D1", "sample_D2", "sample_D3"]
    waterfall_correlation = {}

    for samp in Sample:
        df = compute_program_waterfall_cor(f"{args.perturb_path_base}_{samp}_perturbation_association.txt")
        waterfall_correlation[samp] = (df)


    program_correlation = compute_program_correlation_matrix(mdata)
        
    
    for program in range(program_len):

        create_comprehensive_program_plot(
            mdata=mdata,
            perturb_path=args.perturb_path_base,
            GO_path=args.GO_path,
            waterfall_correlation=waterfall_correlation,
            program_correlation=program_correlation,
            file_to_dictionary=args.file_to_dictionary,
            Target_Program=program,
            top_program=args.top_program,
            groupby=args.groupby,
            tagert_col_name=args.tagert_col_name,
            plot_col_name=args.plot_col_name,
            log2fc_col=args.log2fc_col,
            top_enrichned_term=args.top_enrichned_term,
            down_thred_log=args.down_thred_log,
            up_thred_log=args.down_thred_log,
            p_value=args.p_value, 
            save_path=args.pdf_save_path,
            save_name= str(program),
            figsize=args.figsize,
            samples=args.samples,
            square_plots=args.square_plots,
            show=args.show,
            PDF=args.PDF
        )


    # merge pdf 
    if args.PDF:
        merge_pdfs_in_folder(args.pdf_save_path)
    else:
        merge_svgs_to_pdf(args.pdf_save_path)

