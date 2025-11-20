#%%
import sys
from statsmodels.stats.multitest import fdrcorrection
import argparse
import yaml
import os
#%%

# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from Plotting.src import (load_stablity_error_data, plot_stablity_error,\
                         load_enrichment_data, plot_enrichment,\
                         load_perturbation_data, plot_perturbation,\
                         load_explained_variance_data,plot_explained_variance, programs_dotplots
                          )


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_directory', type=str, required=True)  
    parser.add_argument('--run_name', type=str, required=True)
    parser.add_argument('--groupby', type=str, default="sample")
    parser.add_argument('--K', nargs='*', type=int, default=None) # allow zero input 
    parser.add_argument('--save_folder_name',  type=str, required=True)
    parser.add_argument('--pval',  type=float, default=0.05)
    parser.add_argument('--eval_folder_name',  type=str, required=True)
    parser.add_argument('--sel_threshs', nargs='*', type=float, default=None) # allow zero input 
    parser.add_argument('--samples', nargs='*', type=str, default = None) 


    args = parser.parse_args()

    # either change the array here or run each component in parallel
    if args.K is None:
        k_value = [30, 50, 60, 80, 100, 200, 250, 300]
    else:
        k_value = args.K

    if args.sel_threshs is None:
        sel_thresh_value = [0.4, 0.8, 2.0]
    else:
        sel_thresh_value = args.sel_threshs

    if args.samples is None:
        samples_value = ['D0', 'sample_D1', 'sample_D2', 'sample_D3']
    else:
        samples_value = args.samples


    # save comfigs used         
    args.sel_threshs= sel_thresh_value
    args.K=k_value
    args.samples=samples_value

    args_dict = vars(args)    
    job_id = os.environ.get('SLURM_JOB_ID')

    os.makedirs((f'{args.output_directory}/{args.run_name}/Plot/k_selection'), exist_ok=True)
    with open(f'{args.save_folder_name}/config_{job_id}.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)


    # Stability & Error
    stats = load_stablity_error_data(output_directory = args.output_directory, run_name = args.run_name, components = k_value)
    plot_stablity_error(stats = stats,folder_name = args.save_folder_name, file_name = "Stability_Error")

    for sel_thresh in sel_thresh_value:

        # Enrichement 
        count_df = load_enrichment_data(folder = args.eval_folder_name, components = k_value, sel_thresh = sel_thresh)
        plot_enrichment(count_df,folder_name = args.save_folder_name, file_name = f"Enrichment_{sel_thresh}")

        # Perturbation
        test_stats_df = load_perturbation_data(folder = args.eval_folder_name, components = k_value, sel_thresh = sel_thresh,
        samples = samples_value, pval = args.pval)
        plot_perturbation(test_stats_df, folder_name = args.save_folder_name, pval=args.pval,file_name = f"Perturbation_{sel_thresh}")

        # Explained Variance
        stats = load_explained_variance_data(folder = args.eval_folder_name, components=k_value, sel_thresh = sel_thresh)
        plot_explained_variance(stats, folder_name = args.save_folder_name, file_name = f"Explained_Variance_{sel_thresh}")

        # Motif (working in progress)

    for sel_thresh in sel_thresh_value:
        for k in k_value:
            fig = programs_dotplots(k, args.output_directory, args.run_name, sel_thresh = sel_thresh, groupby=args.groupby, figsize=(4, 15),
            show = False, save_name=f"Program_dotplot_{k}_{sel_thresh}", save_path = args.save_folder_name, ax = None)
