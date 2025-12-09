import os
import sys

import mudata
import pandas as pd
import numpy as np

import seaborn as sns
from matplotlib import pyplot as plt

from qmplot import qqplot

import yaml
import logging
import mudata as mu
import pandas as pd
import cnmf
import scanpy as sc
import logging
import muon as mu
import pandas as pd
import numpy as np
import argparse
import mygene


from Evaluation.src import (
    compute_perturbation_association
)


from Evaluation.src import (
    rename_gene_mygene,
    rename_gene_dictionary,
    check_evaluation_pipeline_format,
    check_gene_names
)


# helper method 
def _assign_guide(mdata, mdata_guide):
        mdata['rna'].var_names = mdata['rna'].var['symbol']
        mdata['cNMF'].uns['guide_names'] = mdata_guide['guide'].var['guide_id']
        mdata['cNMF'].uns['guide_targets'] = mdata_guide['guide'].var['intended_target_name']
        mdata['cNMF'].obsm['guide_assignment'] = mdata_guide['guide'].layers['guide_assignment'].toarray()

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    #IO info
    parser.add_argument('--out_dir', help='path to make cnmf object', type=str, required=True)  
    parser.add_argument('--run_name', help='name of cnmf obj',type=str, required=True)  
    parser.add_argument('--save_path', help='fig saved dirctory',type=str, required=True)  

    # resourses 
    parser.add_argument('--X_normalized_path', type=str,  help='path to normalized input cell x gene matrix from cNMF pipeline', required=True)
    parser.add_argument('--mdata_guide_path', type=str,  help='path to mdata with correct guide assignments', required=True)
    parser.add_argument('--guide_annotation_path', type=str,  help='path to file with guide informations')
    parser.add_argument('--guide_annotation_key', type=str,  help='key to file with guide informations', default = "non-targeting")
    parser.add_argument('--gwas_data_path', type=str,  help='path to file with gwas information for trait enrichment test', required=True)
    parser.add_argument('--reference_gtf_path', type=str,  help='path to reference GTF file for checking gene names')

    # keys
    parser.add_argument('--data_key', help='access gene expression in mdata',type=str, default="rna") 
    parser.add_argument('--prog_key', help='access cNMF program in mdata',type=str,  default="cNMF") 
    parser.add_argument('--categorical_key', help='access cell condition in obs',type=str, default="sample")  
    parser.add_argument('--organism', help='data species',type=str, default="human") 
    parser.add_argument('--FDR_method',type=str, default="StoreyQ")  


    # variables
    parser.add_argument('--num_runs', nargs='*', type=int, default=10)  
    parser.add_argument('--num_guides', nargs='*', type=int, default=5) 
    parser.add_argument('--components', nargs='*', type=int, default=None) # allow zero input 
    parser.add_argument('--sel_threshs', nargs='*', type=float, default=None) # allow zero input 

    # compute real / fake p-value
    parser.add_argument('--compute_real_pertubration_test', action="store_true")
    parser.add_argument('--compute_fake_pertubration_test', action="store_true")


    args = parser.parse_args()

    # read guide
    mdata_guide = mu.read(args.mdata_guide_path)

    # read reference targets
    df_target = pd.read_csv(args.guide_annotation_path, sep = "\t", index_col = 0)
    df_target_non = df_target[df_target["targeting"] == False]
    reference_targets = df_target_non.index.values.tolist()

    # read guide target
    guide_target = pd.read_csv(args.guide_annotation_path, sep='\t', index_col = 0)
    non_targeting_idx = guide_target.index[guide_target['targeting']==False]
    guide_target = guide_target.loc[non_targeting_idx]


    # Compute real p-value, store perturbation test results 
    test_stats_real_df = []
    if args.compute_real_pertubration_test:
        for sel_thresh in args.sel_threshs:
            for k in args.K:  

                output_folder = f"{args.out_dir}/{args.run_name}/Eval/{k}_{str(sel_thresh).replace('.','_')}"

                os.makedirs(output_folder, exist_ok=True)

                # Load mdata
                mdata = mu.read('{out_dir}/{run_name}/adata/cNMF_{k}_{sel_thresh}.h5mu'.format(out_dir = args.out_dir,
                                                                                        run_name = args.run_name,
                                                                                        k=k,
                                                                                        sel_thresh = str(sel_thresh).replace('.','_')))                                                                   
                # assign guide
                _assign_guide(mdata, mdata_guide)   

                # checks for correct mdata format
                if not check_evaluation_pipeline_format(mdata,prog_key=args.prog_key):
                    raise ValueError("mdata format is incorrect")


                valid = check_gene_names(mdata,prog_key=args.prog_key, data_key=args.data_key,categorical_key=args.categorical_key,reference_gtf_path=args.reference_gtf_path)
                if not valid["is_valid"]:
                    raise ValueError("mdata gene naming is incorrect")

                # Run perturbation assocation
                for samp in mdata['rna'].obs[args.categorical_key].unique():
                    mdata_ = mdata[mdata['rna'].obs[args.categorical_key]==samp]
                    test_stats_df = compute_perturbation_association(mdata_, prog_key=args.prog_key, 
                                                                    collapse_targets=True,
                                                                    pseudobulk=False,
                                                                    reference_targets=reference_targets,
                                                                    FDR_method = args.FDR_method,
                                                                    n_jobs=-1, inplace=False)

                    test_stats_df.to_csv('{}/{}_perturbation_association_results_{}.txt'.format(output_folder,k,samp), sep='\t', index=False)

    # read perturbation test results 
    for sel_thresh in args.sel_threshs:
        for k in args.components:  
            for samp in args.sample:
                thresh_str = str(sel_thresh).replace('.', '_')
                test_stats_df_ = pd.read_csv(f'{args.out_dir}/{args.run_name}/Eval/{k}_{thresh_str}/{k}_perturbation_association_results_{samp}.txt', sep='\t')
                test_stats_df_['sample'] = samp
                test_stats_df_['K'] = k
                test_stats_df_['sel_thresh'] = sel_thresh
                #test_stats_df_['fdr'] = fdrcorrection(test_stats_df_['pval'])[1]
                test_stats_real_df.append(test_stats_df_)

    # Concatenate AFTER all loops are done
    test_stats_real_df = pd.concat(test_stats_real_df, ignore_index=True)
    test_stats_real_df['real'] = True



    # Compute fake p-value, store perturbation test results 
    test_stats_fake_dfs = []
    if args.compute_fake_pertubration_test:
        for sel_thresh in args.sel_threshs:
            for k in args.components:  

                # Load mdata
                output_folder = f"{args.out_dir}/{args.run_name}/Eval/{k}_{str(sel_thresh).replace('.','_')}"
                os.makedirs(output_folder, exist_ok=True)
                mdata = mu.read('{out_dir}/{run_name}/adata/cNMF_{k}_{sel_thresh}.h5mu'.format(out_dir = args.out_dir,
                                                                                        run_name =args.run_name,
                                                                                        k=k,
                                                                                        sel_thresh = str(sel_thresh).replace('.','_')))  


                # assign guide
                _assign_guide(mdata, mdata_guide) 
 
                for i in range(args.number_run):
                    print(f"Running iteration {i+1}/{args.number_run} for {k} and {sel_thresh}")       

                    # randomly make number_guide non-targeting guides targeting 
                    guide_target_ = guide_target.copy()
                    selected_guides = np.random.choice(guide_target_.guide_name, args.number_guide, replace=False)
                    guide_target_.loc[guide_target_.guide_name.isin(selected_guides), 'transcript'] = 'targeting'

                    # Filter to only non-targeting guides that exist in both datasets
                    valid_guide_mask = np.isin(mdata['cNMF'].uns['guide_names'], guide_target_.guide_name.values)
                    valid_indices = np.where(valid_guide_mask)[0]
                    
                    print(f"Found {len(valid_indices)} valid guides out of {len(mdata['cNMF'].uns['guide_names'])} total")
                    
                    if len(valid_indices) == 0:
                        raise ValueError("No valid guides found")

                    _mdata = mdata.copy()
                    _mdata['cNMF'].obsm['guide_assignment'] = mdata['cNMF'].obsm['guide_assignment'][:, valid_indices]
                    _mdata['cNMF'].uns['guide_names'] = mdata['cNMF'].uns['guide_names'][valid_indices]
                    _mdata['cNMF'].uns['guide_targets'] = guide_target_.loc[guide_target_.guide_name.isin(mdata['cNMF'].uns['guide_names']), 'transcript'].values


                     # checks for correct mdata format
                    if not check_evaluation_pipeline_format(mdata,prog_key=args.prog_key):
                        raise ValueError("mdata format is incorrect")


                    valid = check_gene_names(mdata,prog_key=args.prog_key, data_key=args.data_key,categorical_key=args.categorical_key,reference_gtf_path=args.reference_gtf_path)
                    if not valid["is_valid"]:
                        raise ValueError("mdata gene naming is incorrect")
 
                    # Run perturbation assocation
                    for samp in _mdata['rna'].obs[args.categorical_key].unique():
                        mdata_samp = _mdata[_mdata['rna'].obs[args.categorical_key]==samp]
                        test_stats_df = compute_perturbation_association(mdata_samp, prog_key=args.prog_key, 
                                                                        collapse_targets=True,
                                                                        pseudobulk=False,
                                                                        reference_targets=reference_targets,
                                                                        FDR_method = args.FDR_method,
                                                                        n_jobs=-1, inplace=False)

                        test_stats_df[args.categorical_key] = samp
                        test_stats_df['K'] = k
                        test_stats_df['run'] = i
                        test_stats_df['sel_thresh'] = sel_thresh
                        test_stats_fake_dfs.append(test_stats_df)

        # save fake test results
        test_stats_fake_dfs = pd.concat(test_stats_fake_dfs)
        test_stats_fake_dfs.to_csv('perturbation_association_calibration.txt', sep='\t', index=False)
    
    # if already exist, read fake pvalues 
    else:
        test_stats_fake_dfs = pd.read_csv('perturbation_association_calibration.txt', sep='\t')
        test_stats_fake_dfs['real'] = False


    # compare 
    test_stats_dfs = pd.concat([test_stats_real_df, test_stats_fake_dfs]) # concate real and fake

    # QQ plot real
    fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(20,8))
    for i,k in enumerate(args.components):

        qqplot(data=test_stats_dfs.loc[(test_stats_dfs.K==k) & (test_stats_dfs.real==True), 'pval'], ax=axs.flat[i])
        axs.flat[i].set_title(str(k))

    axs.flat[-1].axis('off')
    plt.tight_layout()

    plt.savefig(f'{args.save_path}/perturbation_association_qqplot_real.png', bbox_inches='tight', dpi=300)


    # QQ plot fake
    fig, axs = plt.subplots(ncols=4, nrows=2, figsize=(20,8))
    for i,k in enumerate(args.components):

        qqplot(data=test_stats_dfs.loc[(test_stats_dfs.K==k) & (test_stats_dfs.real==False)& (test_stats_dfs.target_name=='targeting'), 'pval'], ax=axs.flat[i])
        axs.flat[i].set_title(str(k))

    axs.flat[-1].axis('off')
    plt.tight_layout()

    plt.savefig(f'{args.save_path}/perturbation_association_qqplot_null.png', bbox_inches='tight', dpi=300)