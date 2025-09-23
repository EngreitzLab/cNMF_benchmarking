import sys
from statsmodels.stats.multitest import fdrcorrection


# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from Plotting.src import (load_stablity_error_data, plot_stablity_error, 
                          load_perturbation_data, plot_perturbation,
                          load_enrichment_data,plot_enrichment_data
                          )


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--output_directory', type=str, required=True)  
    parser.add_argument('--run_name', type=str, required=True)
    parser.add_argument('--K', nargs='*', type=int, default=None) # allow zero input 
    parser.add_argument('--save_folder_name',  type=str, required=True)
    parser.add_argument('--eval_folder_name',  type=str, required=True)


    args = parser.parse_args()

    # either change the array here or run each component in parallel
    if args.K is None:
        k_value = [30, 50, 60, 80, 100, 200, 250, 300]
    else:
        k_value = args.K


    # Stability & Error
    stats = load_stablity_error_data(output_directory = args.output_directory, run_name = args.run_name, components = k_value)
    plot_stablity_error(stats = stats,folder_name = args.save_folder_name, file_name = "Stability_Error")

    # Enrichement 
    count_df = load_enrichment_data(folder = args.eval_folder_name, components = k_value)
    plot_enrichment_data(count_df,folder_name = args.save_folder_name, file_name = "Enrichment")

    # Perturbation
    test_stats_df = load_perturbation_data(folder = args.eval_folder_name, components = k_value)
    plot_perturbation(test_stats_df, folder_name = args.save_folder_name, file_name = "Perturbation")


    # Explained Variance (working in progress)



    # Motif (working in progress for complie the results)