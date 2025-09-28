import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image
from matplotlib import gridspec
import scanpy as sc
from pathlib import Path
from tqdm.auto import tqdm
import scanpy as sc
import anndata
import muon 
import cnmf
import os
import argparse



def run_cnmf_consensus(cnmf_obj=None, output_dir=None, name=None, 
                       components=[7,8,9,10], density_thresholds=[0.01, 0.05, 2.0]):

    if cnmf_obj is None:
        cnmf_obj = init_cnmf_obj(output_dir=output_dir, name=name)

    for k in tqdm(components, desc='Running cNMF'):
        for thresh in density_thresholds:
            cnmf_obj.consensus(k=k, density_threshold=thresh, show_clustering=True)

def compile_results(output_directory, run_name, sel_thresh = 2.0, components = [30, 50, 60, 80, 100, 200, 250, 300] ):
       
    for k in components:

        scores = pd.read_csv('{output_directory}/{run_name}/{run_name}.usages.k_{k}.dt_{sel_thresh}.consensus.txt'.format(
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k,
                                                                                        sel_thresh = str(sel_thresh).replace('.','_')),
                                                                                        sep='\t', index_col=0)

        loadings = pd.read_csv('{output_directory}/{run_name}/{run_name}.spectra.k_{k}.dt_{sel_thresh}.consensus.txt'.format(
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k,
                                                                                        sel_thresh = str(sel_thresh).replace('.','_')),
                                                                                        sep='\t', index_col=0)
        

        os.makedirs((f'{output_directory}/{run_name}/loading'), exist_ok=True)


        scores.to_csv('{output_directory}/{run_name}/loading/cNMF_scores_{k}_{sel_thresh}.txt'.format(
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k,
                                                                                        sel_thresh = sel_thresh), sep='\t')
        loadings.T.to_csv('{output_directory}/{run_name}/loading/cNMF_loadings_{k}_{sel_thresh}.txt'.format(     
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k,
                                                                                        sel_thresh = sel_thresh), sep='\t')

        adata_ = anndata.read_h5ad('{output_directory}/{run_name}/cnmf_tmp/{run_name}.tpm.h5ad'.format(
                                                                                        output_directory=output_directory,
                                                                                        run_name = run_name,
                                                                                        k=k ))
        adata_.var_names_make_unique()
        adata_.obs_names_make_unique()

        prog_data = anndata.AnnData(X=scores.values, obs=adata_.obs)
        prog_data.varm['loadings'] = loadings.values
        prog_data.uns['var_names'] = loadings.columns.values


        # Make adata
        os.makedirs((f'{output_directory}/{run_name}/prog_data'), exist_ok=True)
        prog_data.write(f'{output_directory}/{run_name}/prog_data/NMF_{k}_{sel_thresh}.h5ad'.format(
                                                                                output_directory=output_directory,
                                                                                run_name = run_name,
                                                                                k=k,
                                                                                sel_thresh = str(sel_thresh).replace('.','_')))

        # Make mdata
        mdata = muon.MuData({'rna': adata_, 'cNMF': prog_data})

        os.makedirs((f'{output_directory}/{run_name}/adata'), exist_ok=True)
        mdata.write(f'{output_directory}/{run_name}/adata/cNMF_{k}_{sel_thresh}.h5mu'.format(
                                                                                output_directory=output_directory,
                                                                                run_name = run_name,
                                                                                k=k,
                                                                                sel_thresh = str(sel_thresh).replace('.','_')))

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
    parser.add_argument('--n_iter', type = int, default = 10)
    parser.add_argument('--loss', default = 'frobenius')
    parser.add_argument('--algo', type = str, default = 'mu')
    parser.add_argument('--max_NMF_iter', type = int , default = 500)
    parser.add_argument('--tol', type = float , default = 1e-4)


    parser.add_argument('--sel_thresh', nargs='*', type=float, default=[2.0])  



    args = parser.parse_args()

    # either change the array here or run each component in parallel
    if args.K is None:
        k_value = [30, 50, 60, 80, 100, 200, 250, 300]
    else:
        k_value = args.K


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
                        density_thresholds=args.sel_thresh)

    # Save all cNMF scores in separate mudata objects
    compile_results(args.output_directory, args.run_name)

    # save comfigs used         
    args_dict = vars(args)
    with open(f'{args.output_directory}/{run_name}/config.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)