import sys
import argparse
import cnmf
import yaml
import os
import pandas as pd

# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from Inference.src import (
    run_cnmf_consensus, get_top_indices_fast, annotate_genes_to_excel, \
    rename_and_move_files_NMF, rename_all_NMF, compile_results
)

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--counts_fn', type=str, required=True)  
    parser.add_argument('--output_directory', type=str, required=True)
    parser.add_argument('--run_name', type=str, required=True)

    parser.add_argument('--numiter', type = int, default = 10)
    parser.add_argument('--numhvgenes', type = int, default = 5451)
    parser.add_argument('--seed', type = int, default = 14)
    parser.add_argument('--K', nargs='*', type=int, default=None) # allow zero input 
    parser.add_argument('--init', type = str, default = 'random')
    parser.add_argument('--loss', default = 'frobenius')
    parser.add_argument('--algo', type = str, default = 'mu')
    parser.add_argument('--max_NMF_iter', type = int , default = 500)
    parser.add_argument('--tol', type = float , default = 1e-4)

    parser.add_argument('--sel_thresh', nargs='*', type=float, default=[2.0])  

    parser.add_argument('--parallel_running', action="store_true")
    parser.add_argument('--num_gene', type = int, default = 300)

    args = parser.parse_args()

    # either change the array here or run each component in parallel
    if args.K is None:
        k_value = [30, 50, 60, 80, 100, 200, 250, 300]
    else:
        k_value = args.k
    
    if args.sel_thresh is None:
        sel_thresh_value = [0.2, 2.0]
    else:
        sel_thresh_value = args.sel_thresh

    # save comfigs used      
    os.makedirs((f'{args.output_directory}/{args.run_name}'), exist_ok=True)
    args_dict = vars(args)
    with open(f'{args.output_directory}/{args.run_name}/config.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)


    # running cnmf 
    cnmf_obj = cnmf.cNMF(output_dir=args.output_directory, name=args.run_name)

    cnmf_obj.prepare(counts_fn= args.counts_fn, components= k_value, n_iter= args.numiter,  densify=False, tpm_fn=None, seed= args.seed,
                     beta_loss = args.loss,num_highvar_genes=args.numhvgenes, genes_file=None,
                     alpha_usage=0.0, alpha_spectra=0.0, init=args.init, max_NMF_iter=args.max_NMF_iter, algo = args.algo, tol = args.tol)


    cnmf_obj.factorize(total_workers = 1)

    cnmf_obj.combine()

    cnmf_obj.k_selection_plot()

    # Consensus plots with all k to choose thresh
    run_cnmf_consensus(cnmf_obj, 
                        components=k_value, 
                        density_thresholds=sel_thresh_value)

    # Save all cNMF scores in separate mudata objects
    compile_results(args.output_directory, args.run_name, components = k_value, sel_thresh = sel_thresh_value)


    # annotation for all K
    os.makedirs((f'{args.output_directory}/{args.run_name}/Annotation'), exist_ok=True)

    for i in sel_thresh_value:
        for k in k_value:
            df = pd.read_csv('{output_directory}/{run_name}/{run_name}.gene_spectra_scores.k_{k}.dt_{sel_thresh}.txt'.format(
                                                                                    output_directory=args.output_directory,
                                                                                    run_name = args.run_name,
                                                                                    k=k,
                                                                                    sel_thresh = str(i).replace('.','_')),
                                                                                    sep='\t', index_col=0)   
            overlap = get_top_indices_fast(df, gene_num=args.num_gene)
            annotate_genes_to_excel(overlap, f'{args.output_directory}/{args.run_name}/Annotation/{k}.xlsx')


    # combine the parallel ran K value into "run_name_all" file 
    if args.parallel_running and isinstance(k_value, int):

        file_name_input_new = f"{file_name_input}_{k_value}.spectra.k_{k_value}.iter"
        file_name_output_new = f"{args.run_name}_all.spectra.k_{k_value}.iter"

        source_folder_new = f"{source_folder}_{k_value}/cnmf_tmp"

        rename_all_NMF(source_folder = args.output_directory ,
                                destination_folder = f"{args.utput_directory}_all/cnmf_tmp",
                                file_name_input = args.run_name,
                                file_name_output = f"{args.run_name}_all",
                                len = args.numiter, 
                                components = k_value)

