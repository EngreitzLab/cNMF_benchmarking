import sys
import argparse
import cnmf
import yaml
import os
import pandas as pd
import muon as mu


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

    # cNMF parameters 
    parser.add_argument('--K', nargs='*', type=int, default=None)
    parser.add_argument('--numiter', type = int, default = 10)
    parser.add_argument('--densify', action='store_true') # Flag: include --densify for True, omit for False
    parser.add_argument('--tpm_fn', type = str, default = None)
    parser.add_argument('--numhvgenes', type = int, default = 5451)
    parser.add_argument('--genes_file', type = str, default = None)
    parser.add_argument('--init', type = str, default = 'random')
    parser.add_argument('--loss', default = 'frobenius')
    parser.add_argument('--algo', type = str, default = 'mu')
    parser.add_argument('--mode', type = str, default = 'batch')
    parser.add_argument('--tol', type = float, default = 1e-4)
    parser.add_argument('--n_jobs', type = int, default = -1)
    parser.add_argument('--seed', type = int, default = 14)
    parser.add_argument('--use_gpu', action="store_true") # Flag: include --use_gpu for True, omit for False
    parser.add_argument('--alpha_usage', type = float, default = 0.0)
    parser.add_argument('--alpha_spectra', type = float, default = 0.0)
    parser.add_argument('--l1_ratio_usage', type = float, default = 0.0)
    parser.add_argument('--l1_ratio_spectra', type = float, default = 0.0)
    parser.add_argument('--online_usage_tol', type = float, default = 0.05)
    parser.add_argument('--online_spectra_tol', type = float, default = 0.05)
    parser.add_argument('--fp_precision', type = str, default = 'float')
    parser.add_argument('--batch_max_iter', type = int, default = 500)
    parser.add_argument('--batch_hals_tol', type = float, default = 0.05)
    parser.add_argument('--batch_hals_max_iter', type = int, default = 200)
    parser.add_argument('--online_max_pass', type = int, default = 20)
    parser.add_argument('--online_chunk_size', type = int, default = 5000)
    parser.add_argument('--online_chunk_max_iter', type = int, default = 200)
    parser.add_argument('--shuffle_cells', action="store_true")
    parser.add_argument('--sel_thresh', nargs='*', type=float, default=None) 

    # annotation parameters 
    parser.add_argument('--sk_cd_refit', action="store_true")
    parser.add_argument('--species', type=str, required=True)
    parser.add_argument('--check_format', action="store_true")
    parser.add_argument('--parallel_running', action="store_true")
    parser.add_argument('--num_gene', type = int, default = 300)
    parser.add_argument('--run_refit_only', action="store_true")

    # resourses 
    parser.add_argument('--guide_annotation_path', type=str,  help='path to file with guide informations', default=None)
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


    # save comfigs used  

    # create output directory       
    os.makedirs((f'{args.output_directory}/{args.run_name}'), exist_ok=True)

    # Get args as dict
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


    # running 
    cnmf_obj = cnmf.cNMF(output_dir=args.output_directory, name=args.run_name)


    cnmf_obj.prepare(counts_fn=args.counts_fn, components=args.K, n_iter=args.numiter, densify=args.densify, tpm_fn=args.tpm_fn, num_highvar_genes=args.numhvgenes, genes_file=args.genes_file,
                init = args.init,  beta_loss = args.loss, 
                algo = args.algo, mode = args.mode, tol=args.tol, n_jobs=args.n_jobs, 
                seed = args.seed,  use_gpu = args.use_gpu, 
                alpha_usage = args.alpha_usage, alpha_spectra = args.alpha_spectra, 
                l1_ratio_usage = args.l1_ratio_usage, l1_ratio_spectra = args.l1_ratio_spectra,
                online_usage_tol = args.online_usage_tol, online_spectra_tol = args.online_spectra_tol,
                fp_precision = args.fp_precision, 
                batch_max_iter = args.batch_max_iter, batch_hals_tol = args.batch_hals_tol, batch_hals_max_iter = args.batch_hals_max_iter,
                online_max_pass = args.online_max_pass, online_chunk_size = args.online_chunk_size, online_chunk_max_iter = args.online_chunk_max_iter,
                shuffle_cells = args.shuffle_cells, sk_cd_refit =args.sk_cd_refit )


    if not args.run_refit_only:

        cnmf_obj.factorize(total_workers=1)

    cnmf_obj.combine()

    cnmf_obj.k_selection_plot()

    # Consensus plots with all k to choose thresh
    run_cnmf_consensus(cnmf_obj, 
                       components=args.K, 
                       density_thresholds=args.sel_thresh)

    # annotation for all K
    os.makedirs((f'{args.output_directory}/{args.run_name}/Annotation'), exist_ok=True)
    for i in args.sel_thresh:
        for k in args.K:
            df = pd.read_csv('{output_directory}/{run_name}/{run_name}.gene_spectra_score.k_{k}.dt_{sel_thresh}.txt'.format(
                                                                                    output_directory=args.output_directory,
                                                                                    run_name = args.run_name,
                                                                                    k=k,
                                                                                    sel_thresh = str(i).replace('.','_')),
                                                                                    sep='\t', index_col=0)   
            overlap = get_top_indices_fast(df, gene_num=args.num_gene)
            annotate_genes_to_excel(overlap, species = args.species, output_file = f'{args.output_directory}/{args.run_name}/Annotation/{k}_{i}.xlsx')



    # combine the parallel ran K value into "run_name_all" file 
    if args.parallel_running and isinstance(args.K, int):

        rename_all_NMF(source_folder = f"{args.utput_directory}/{args.run_name}" ,
                                destination_folder = f"{args.utput_directory}/{args.run_name}_all/cnmf_tmp",
                                file_name_input = args.run_name,
                                file_name_output = f"{args.run_name}_all",
                                len = args.numiter, 
                                components = [args.K])

    # Save all cNMF scores in separate mudata objects
    compile_results(args.output_directory, args.run_name, components= args.K, sel_thresh = args.sel_thresh,
     guide_names_key = args.guide_names_key, guide_targets_key = args.guide_targets_key, categorical_key= args.categorical_key, 
     guide_assignment_key = args.guide_assignment_key )



