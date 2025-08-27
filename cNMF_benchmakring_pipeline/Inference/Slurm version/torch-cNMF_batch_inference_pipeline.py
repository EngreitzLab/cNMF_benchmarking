import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import Image
from matplotlib import gridspec
import scanpy as sc
from pathlib import Path
import anndata
import muon 
from tqdm.auto import tqdm
import cnmf
import os
import argparse
import gin



def run_cnmf_consensus(cnmf_obj=None, output_dir=None, name=None, 
                       components=[7,8,9,10], density_thresholds=[0.01, 0.05, 2.0]):

    if cnmf_obj is None:
        cnmf_obj = init_cnmf_obj(output_dir=output_dir, name=name)

    for k in tqdm(components, desc='Running cNMF'):
        for thresh in density_thresholds:
            cnmf_obj.consensus(k=k, density_threshold=thresh, show_clustering=True)
    



if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('count_fn')
    parser.add_argument('--output_directory', type = str)
    parser.add_argument('--run_name', type = str)

    parser.add_argument('--numiter', type = int, default = 10)
    parser.add_argument('--numhvgenes', type = int, default = 5451)
    parser.add_argument('--seed', type = int, default = 14)
    parser.add_argument('--K', type = list, default = [30, 50, 60, 80, 100, 200, 250, 300])
    parser.add_argument('--inti', type = str, default = 'random')
    parser.add_argument('--n_iter', type = int, default = 10)
    parser.add_argument('--loss', type = str | int, default = 'frobenius')
    parser.add_argument('--algo', type = str, defult = 'mu')
    parser.add_argument('--sel_thresh', type = list , default = [2.0])


    # running 
    cnmf_obj = cnmf.cNMF(output_dir=output_directory, name=run_name)


    cnmf_obj.prepare(counts_fn=counts_fn, components=K, n_iter=numiter, seed=seed, init = init, total_workers=1, num_highvar_genes=numhvgenes,
                beta_loss = loss, use_gpu = True, algo = algo, mode = mode, tol=tol, batch_max_iter=batch_max_iter,online_max_pass=online_max_pass,
                online_chunk_size=online_chunk_size, online_chunk_max_iter=online_chunk_max_iter, online_usage_tol=online_usage_tol, online_spectra_tol=online_spectra_tol)


    cnmf_obj.factorize(total_workers = 1)

    cnmf_obj.combine()

    cnmf_obj.k_selection_plot()

    # Consensus plots with all k to choose thresh
    run_cnmf_consensus(cnmf_obj, 
                    components=K, 
                    density_thresholds=sel_thresh)

    # Save all cNMF scores in separate mudata objects
    for k in K:
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

        os.makedirs((f'{output_directory}/{run_name}/prog_data'), exist_ok=True)
        prog_data.write(f'{output_directory}/{run_name}/cNMF_{k}_{sel_thresh}.h5ad')

        # Make mdata
        mdata = muon.MuData({'rna': adata_, 'cNMF': prog_data})

        os.makedirs((f'{output_directory}/{run_name}/adata'), exist_ok=True)
        mdata.write('{output_directory}/{run_name}/adata/cNMF_{k}_{sel_thresh}.h5mu'.format(
                                                                                output_directory=output_directory,
                                                                                run_name = run_name,
                                                                                k=k,
                                                                                sel_thresh = str(sel_thresh).replace('.','_')))

