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
                         compute_gene_correlation_matrix, compute_gene_waterfall_cor,perturbed_program_dotplot


def _assign_guide(mdata, mdata_guide):
        mdata['cNMF'].obsm['guide_assignment'] = mdata_guide['cNMF'].obsm['guide_assignment'].toarray()


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    #io path
    parser.add_argument('--mdata_path', type=str, required=True, help='path to the MuData (.h5mu) file')
    parser.add_argument('--perturb_path_base', type=str, required=True, help='base path for perturbation result files (sample suffix appended automatically)')
    parser.add_argument('--file_to_dictionary', type=str, default=None, help='path to gene name mapping dictionary file for ID-to-name conversion')
    parser.add_argument('--reference_gtf_path', type=str, default=None, help='path to reference GTF file for checking gene names')

    # plotting variables
    parser.add_argument('--tagert_col_name', type=str, default="target_name", help='column name for target genes in perturbation results')
    parser.add_argument('--plot_col_name', type=str, default="program_name", help='column name for programs in perturbation results')
    parser.add_argument('--log2fc_col', type=str, default="log2FC", help='column name for log2 fold change values')
    parser.add_argument('--top_num', type=int, default=5, help='number of top genes to display per program')
    parser.add_argument('--top_program', type=int, default=10, help='number of top programs to display per gene')
    parser.add_argument('--p_value', type=float, default=0.05, help='p-value threshold for significance')
    parser.add_argument('--down_thred_log', type=float, default=-0.05, help='lower log2FC threshold for volcano plot')
    parser.add_argument('--up_thred_log', type=float, default=0.05, help='upper log2FC threshold for volcano plot')
    parser.add_argument('--pdf_save_path', type=str, required=True, help='directory path to save output plots')
    parser.add_argument('--square_plots', action="store_true", help='use square aspect ratio for plots')
    parser.add_argument('--figsize', type=float, nargs=2, default=(35, 35), help='figure size as width height')
    parser.add_argument('--show', action="store_true", help='display plots interactively')
    parser.add_argument('--PDF', action="store_true", help='save plots as PDF (default is SVG)')
    parser.add_argument('--n_processes', type=int, default=-1, help='number of parallel processes (-1 for all available cores)')
    parser.add_argument('--sample', nargs='*', type=str, default=None, help='list of sample names (default: D0 sample_D1 sample_D2 sample_D3)')
    parser.add_argument('--dot_size', type=int, default=10, help='dot size use to plot UMAP')
    parser.add_argument('--expressed_only', action="store_true", help='only plot perturbed genes found in the gene expression matrix (default: plot all perturbed genes)')
    parser.add_argument('--subsample_frac', type=float, default=None, help='fraction of cells to subsample for UMAP plots (e.g. 0.1 for 10%%). Default: None (plot all cells)')

    # keys
    parser.add_argument('--data_key', type=str, default="rna", help='key to access gene expression data in MuData')
    parser.add_argument('--prog_key', type=str, default="cNMF", help='key to access cNMF programs in MuData')
    parser.add_argument('--gene_name_key', type=str, default="gene_names", help='key to access gene names in var')
    parser.add_argument('--categorical_key', type=str, default="sample", help='key to access sample/condition labels in obs')

    
    args = parser.parse_args()

    if args.sample is None:
        args.sample = ['D0', 'sample_D1', 'sample_D2', 'sample_D3']



    # save comfigs used         
    args_dict = vars(args)
    job_id = os.environ.get('SLURM_JOB_ID')
    os.makedirs(f'{args.pdf_save_path}', exist_ok=True)
    with open(f'{args.pdf_save_path}/config_{job_id}.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)



    #read mdata
    mdata = mu.read_h5mu(args.mdata_path)
    _assign_guide(mdata, mdata)

 
    # check umap exist 
    if 'X_umap' not in mdata['cNMF'].obsm:
        import scanpy as sc
        sc.tl.pca(mdata['rna'], n_comps=50)
        sc.pp.neighbors(mdata['rna'])
        sc.tl.umap(mdata['rna'])
        mdata['cNMF'].obsm['X_pca'] = mdata['rna'].obsm['X_pca'] 
        mdata['cNMF'].obsm['X_umap'] = mdata['rna'].obsm['X_umap'] 




    # found detected perturbed gene
    perturbed_gene = np.unique(mdata['cNMF'].uns["guide_targets"])
    gene_list = mdata['rna'].var_names.tolist()
    # gene_list = rename_list_gene_dictionary(mdata['rna'].var_names.tolist(), args.file_to_dictionary) # convert gene id to gene name
    perturbed_gene_found = sorted(set(gene_list) & set(perturbed_gene.tolist()))
    perturbed_gene_not_found = sorted(set(perturbed_gene.tolist()) - set(gene_list))
    print(f"there are {len(perturbed_gene_found)} perturbed genes found in expression matrix")
    print(f"there are {len(perturbed_gene_not_found)} perturbed genes NOT found in expression matrix: {perturbed_gene_not_found}")

    if args.expressed_only:
        genes_to_plot = perturbed_gene_found
    else:
        genes_to_plot = sorted(perturbed_gene.tolist())


    # compute corr once
    correlation_matrix = compute_gene_correlation_matrix(mdata, file_to_dictionary = args.file_to_dictionary)

    waterfall_correlation = {}
    for samp in args.sample:
        df = compute_gene_waterfall_cor(f"{args.perturb_path_base}_{samp}.txt")
        waterfall_correlation[samp] = (df)



   
    # Graph all pdf
    for gene in genes_to_plot:

        create_comprehensive_plot(
            mdata = mdata,
            perturb_path_base = args.perturb_path_base,
            file_to_dictionary = args.file_to_dictionary,
            Target_Gene = gene,
            gene_correlation = correlation_matrix,
            waterfall_correlation=waterfall_correlation,
            top_program=args.top_program,
            groupby=args.categorical_key,
            tagert_col_name=args.tagert_col_name,
            plot_col_name=args.plot_col_name,
            log2fc_col=args.log2fc_col,
            top_num=args.top_num,
            down_thred_log=args.down_thred_log,
            up_thred_log=args.up_thred_log,
            p_value=args.p_value,
            save_path = args.pdf_save_path,
            save_name = gene,
            figsize = args.figsize,
            sample= args.sample,
            square_plots=args.square_plots,
            show=args.show,
            PDF = True,
            dot_size = args.dot_size,
            subsample_frac = args.subsample_frac
        )
        
    '''
    print("Starting parallel gene processing...")                   

    try: 
        # Graph all pdf 
        result = parallel_gene_processing( 
            mdata = mdata,
            perturb_path_base = args.perturb_path_base,
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
            sample=args.sample,
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
        merge_pdfs_in_folder(args.pdf_save_path, output_filename = "gene.pdf")
    else:
        merge_svgs_to_pdf(args.pdf_save_path)

