import sys
import muon as mu 
import numpy as np
import pandas as pd
import argparse

# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from Plotting.src import plot_umap_per_gene, plot_top_program_per_gene, perturbed_gene_dotplot,\
                                      plot_log2FC, plot_volcano, programs_dotplots, analyze_correlations, \
                                      create_comprehensive_plot, create_gene_correlation_waterfall, convert_with_mygene, \
                                      convert_adata_with_mygene, read_npz, merge_svgs_to_pdf, merge_pdfs_in_folder



if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--mdata_path', type=str, required=True)  
    parser.add_argument('--gene_loading_path', type=str, required=True)
    parser.add_argument('--program_loading_path', type=str, required=True)
    parser.add_argument('--perturb_path', type=str, required=True) # partial path for each day
    parser.add_argument('--top_program',  type=int, default=5)
    parser.add_argument('--p_value',  type=float, default=0.05)
    parser.add_argument('--pdf_save_path',  type=str, required=True)
    parser.add_argument('--PDF',  action="store_true")  



    args = parser.parse_args()


    #read mdata
    mdata = mu.read_h5mu(args.mdata_path)

    # replace name
    replacement_dict = {
        'sample_D3': 'D3',
        'sample_D2': 'D2', 
        'sample_D1': 'D1'
    }  

    mdata['rna'].obs['sample'] = mdata['rna'].obs['sample'].replace(replacement_dict)

    perturbed_gene = np.unique(mdata['cNMF'].uns["guide_targets"])

    # found detected perturbed gene
    gene_list = convert_with_mygene(pd.DataFrame(index=mdata.var_names.tolist())) # convert gene id to geene name
    perturbed_gene_found = list(set(gene_list.index) & set(perturbed_gene))
    print(f"there are {len(perturbed_gene_found)} perturbed gene found")


    # Graph all pdf 
    for gene in perturbed_gene_found:

        create_comprehensive_plot(
            mdata,
            args.perturb_path,
            args.gene_loading_path,
            args.program_loading_path,
            gene,
            save_path=args.pdf_save_path,
            save_name=f"{gene}",
            figsize=(25, 25),
            show=False,
            PDF = args.PDF
        )

    # merge pdf 
    if args.PDF:
        merge_pdfs_in_folder(args.pdf_save_path)
    else:
        merge_svgs_to_pdf(args.pdf_save_path)

    # save comfigs used         
    args_dict = vars(args)
    with open(f'{args.input_folder}/Eval/config.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)
        