import sys
import argparse
import cnmf
import yaml
import os
import pandas as pd
import muon as mu
import numpy as np

# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from Inference.src import (
    run_cnmf_consensus, get_top_indices_fast, annotate_genes_to_excel, \
    rename_and_move_files_NMF, rename_all_NMF, compile_results
)

from Inference.src import (
    check_data_format, check_guide_names, _validate_against_reference_gtf, check_mdata_format )

    
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    #IO
    parser.add_argument('--counts_fn', type=str, required=True)  
    parser.add_argument('--output_directory', type=str, required=True)
    parser.add_argument('--run_name', type=str, required=True)
    parser.add_argument('--nmf_seeds_path', type=str, required=True)


    # cNMF parameters 
    parser.add_argument('--numiter', type = int, default = 10)
    parser.add_argument('--numhvgenes', type = int, default = 5451)
    parser.add_argument('--seed', type = int, default = 14)
    parser.add_argument('--K', nargs='*', type=int, default=None) # allow zero input 
    parser.add_argument('--init', type = str, default = 'random')
    parser.add_argument('--loss', default = 'frobenius')
    parser.add_argument('--algo', type = str, default = 'mu')
    parser.add_argument('--max_NMF_iter', type = int , default = 500)
    parser.add_argument('--tol', type = float , default = 1e4)
    parser.add_argument('--sel_thresh', nargs='*', type=float, default=[2.0])  

    # annotation parameters 
    parser.add_argument('--species', type=str, required=True)
    parser.add_argument('--check_format', action="store_true")
    parser.add_argument('--parallel_running', action="store_true")
    parser.add_argument('--num_gene', type = int, default = 300)
    parser.add_argument('--run_refit', action="store_true")
    parser.add_argument('--run_complie_annotation', action="store_true")
    parser.add_argument('--run_factorize', action="store_true")

    # resourses 
    parser.add_argument('--guide_annotation_path', type=str,  help='path to file with guide informations',  default=None)
    parser.add_argument('--reference_gtf_path', type=str,  help='path to reference GTF file for checking gene names', default=None)

    # keys
    parser.add_argument('--data_key', help='access gene expression in mdata',type=str, default="rna") 
    parser.add_argument('--prog_key', help='access cNMF program in mdata',type=str,  default="cNMF") 
    parser.add_argument('--categorical_key', help='access cell condition in obs',type=str, default="sample")  
    parser.add_argument('--guide_names_key', help='guide names in uns',type=str, default="guide_names")  
    parser.add_argument('--guide_targets_key', help='guide targets in uns',type=str, default="guide_targets") 
    parser.add_argument('--guide_assignment_key', help='guide assignment in obsm',type=str, default="guide_assignment_key") 



    args = parser.parse_args()

    # either change the array here or run each component in parallel
    if args.K is None:
        args.K = [30, 50, 60, 80, 100, 200, 250, 300]
    
    if args.sel_thresh is None:
        args.sel_thresh = [0.4, 0.8, 2.0]


    # read seeds
    if args.nmf_seeds_path is not None:
        nmf_seeds = np.load(args.nmf_seeds_path)
    else:
        nmf_seeds = None



    # create output directory       
    os.makedirs((f'{args.output_directory}/{args.run_name}'), exist_ok=True)

    args_dict = vars(args)

   # --- Capture Slurm environment info ---
    slurm_info = {
            'job_id': os.environ.get('SLURM_JOB_ID'),
            'job_name': os.environ.get('SLURM_JOB_NAME'),
            'partition': os.environ.get('SLURM_JOB_PARTITION'),
            'node_list': os.environ.get('SLURM_JOB_NODELIST'),
            'num_nodes': os.environ.get('SLURM_JOB_NUM_NODES'),
            'ntasks': os.environ.get('SLURM_NTASKS'),
            'cpus_per_task': os.environ.get('SLURM_CPUS_PER_TASK'),
            'mem_per_node': os.environ.get('SLURM_MEM_PER_NODE'),
            'mem_per_cpu': os.environ.get('SLURM_MEM_PER_CPU'),
            'gres': os.environ.get('SLURM_JOB_GRES'),
            'gpu_device_ids': os.environ.get('SLURM_GPUS_ON_NODE'),
            'gpu_type': os.environ.get('SLURM_JOB_GPUS'),
            'time_limit': os.environ.get('SLURM_JOB_TIMELIMIT'),
            'time_remaining': os.environ.get('SLURM_TIMELIMIT'),
            'submit_dir': os.environ.get('SLURM_SUBMIT_DIR'),
            'submit_host': os.environ.get('SLURM_SUBMIT_HOST'),
            'constraint': os.environ.get('SLURM_JOB_CONSTRAINT'),
            'array_task_id': os.environ.get('SLURM_ARRAY_TASK_ID'),
        }
    job_id = slurm_info['job_id'] or 'no_jobid'

    # Merge them into your config
    config_to_save = {
        'script_args': args_dict,
        'slurm_info': slurm_info
    }

    with open(f'{args.output_directory}/{args.run_name}/config_{job_id}.yml', 'w') as f:
        yaml.dump(config_to_save, f, default_flow_style=False, width=1000)

    # check data format 
    if args.check_format:
        adata = mu.read(args.counts_fn)
        valid = check_guide_names(adata, guide_names_key = args.guide_names_key, guide_targets_key = args.guide_targets_key, 
        categorical_key= args.categorical_key, reference_gtf_path=args.reference_gtf_path, guide_annotation_path = args.guide_annotation_path)
        if not valid['is_valid']:
            raise ValueError("Format is incorrect")

    # running cnmf 
    cnmf_obj = cnmf.cNMF(output_dir=args.output_directory, name=args.run_name)

    cnmf_obj.prepare(counts_fn= args.counts_fn, components= args.K, n_iter= args.numiter,  densify=False, tpm_fn=None, seed= args.seed,
                     beta_loss = args.loss,num_highvar_genes=args.numhvgenes, genes_file=None,
                     alpha_usage=0.0, alpha_spectra=0.0, init=args.init, max_NMF_iter=args.max_NMF_iter, algo = args.algo, tol = args.tol, nmf_seeds = nmf_seeds)

    if args.run_factorize:

        cnmf_obj.factorize(total_workers = 1,skip_completed_runs=True) 

    if args.run_refit:

        cnmf_obj.combine()

        cnmf_obj.k_selection_plot()

        # Consensus plots with all k to choose thresh
        run_cnmf_consensus(cnmf_obj, 
                            components=args.K, 
                            density_thresholds=args.sel_thresh)


    if args.run_complie_annotation:

        # Save all cNMF scores in separate mudata objects
        compile_results(args.output_directory, args.run_name, components= args.K, sel_threshs = args.sel_thresh,
        guide_names_key = args.guide_names_key, guide_targets_key = args.guide_targets_key, categorical_key= args.categorical_key, 
        guide_assignment_key = args.guide_assignment_key )

        # annotation for all K
        os.makedirs((f'{args.output_directory}/{args.run_name}/Annotation'), exist_ok=True)

        # annotation for all K
        for i in args.sel_thresh:
            for k in args.K:
                df = pd.read_csv('{output_directory}/{run_name}/{run_name}.gene_spectra_scores.k_{k}.dt_{sel_thresh}.txt'.format(
                                                                                        output_directory=args.output_directory,
                                                                                        run_name = args.run_name,
                                                                                        k=k,
                                                                                        sel_thresh = str(i).replace('.','_')),
                                                                                        sep='\t', index_col=0)   
                overlap = get_top_indices_fast(df, gene_num=args.num_gene)
                annotate_genes_to_excel(overlap, species = args.species, output_file= f'{args.output_directory}/{args.run_name}/Annotation/{k}.xlsx')



    # combine the parallel ran K value into "run_name_all" file 
    if args.parallel_running and isinstance(args.K, int):

        rename_all_NMF(source_folder = f"{args.utput_directory}/{args.run_name}" ,
                                destination_folder = f"{args.utput_directory}/{args.run_name}_all/cnmf_tmp",
                                file_name_input = args.run_name,
                                file_name_output = f"{args.run_name}_all",
                                len = args.numiter, 
                                components = [args.K])





