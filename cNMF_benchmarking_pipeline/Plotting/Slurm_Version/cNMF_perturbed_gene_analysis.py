import sys
import muon as mu 
import numpy as np
import pandas as pd

# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from Plotting.src import plot_umap_per_gene, plot_top_program_per_gene, perturbed_gene_dotplot,\
                         plot_log2FC, plot_volcano, programs_dotplots, analyze_correlations, \
                         graph_pdf_gene_QC, convert_with_mygene, convert_adata_with_mygene, read_npz, merge_pdfs_in_folder


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--mdata_path', type=str, required=True)  
    parser.add_argument('--gene_loading_path', type=str, required=True)
    parser.add_argument('--program_loading_path', type=str, required=True)
    parser.add_argument('--perturb_path', type=str, required=True) # partial path for each day
    parser.add_argument('--top_program',  type=int, default=5)
    parser.add_argument('--p_value',  type=float, default=0.05)
    parser.add_argument('--pdf_save_path',  type=str, required=True)


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

        save_name = f"{gene} QC"

        graph_pdf_gene_QC(mdata, args.gene_loading_path, args.program_loading_path, args.perturb_path, Target_Gene = gene, top_program= args.top_program,\
             save_path=args.pdf_save_path, save_name=save_name, p_value = args.p_value)

    # merge pdf 
    merge_path = f"{args.pdf_save_path}/Perturbed_Gene_QC"
    merge_pdfs_in_folder(merge_path)
        