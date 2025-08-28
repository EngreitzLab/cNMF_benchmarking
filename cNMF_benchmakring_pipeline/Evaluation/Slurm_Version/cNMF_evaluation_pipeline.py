import os
import sys
import yaml
import logging
import mudata
import pandas as pd

# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline')

from Evaluation.src import (
    compute_categorical_association,
    compute_geneset_enrichment,
    compute_trait_enrichment,
    compute_perturbation_association,
    compute_explained_variance_ratio,
    compute_motif_enrichment
)
from Evaluation.src.enrichment_trait import process_enrichment_data

def assign_guide(mdatam,file):

    # read mdata with guide
    mdata_guide = mudata.read(file)

    mdata['cNMF'].uns["guide_names"] = mdata_guide["cNMF_100"].uns["guide_names"]
    mdata['cNMF'].uns["guide_targets"] = mdata_guide["cNMF_100"].uns["guide_targets"]
    mdata['cNMF'].obsm["guide_assignment"] = mdata_guide["cNMF_100"].obsm["guide_assignment"]


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--input_folder', type=str, required=True)  
    parser.add_argument('--mdata_guide', type=str, default="/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Evaluation/Resources/mdata_guide.h5mu")
    parser.add_argument('--K', default = [30, 60, 80, 100, 200, 250, 300])
    parser.add_argument('--thre', default = '2_0')


    for k in args.K:  

        os.makedirs(f"{args.input_folder}/Eval/{k}", exist_ok=True)

        output_folder = f"{args.input_folder}/Eval/{k}"

        # Load mdata
        mdata = mudata.read('{}/adata/cNMF_{}_{}.h5mu'.format(args.input_folder,args.thre,k))
        assign_guide(mdata, args.mdata_guide ) 


        # Run categorical assocation
        results_df, posthoc_df = compute_categorical_association(mdata, prog_key='cNMF', categorical_key='sample', 
                                                                pseudobulk_key=None, test='dunn', n_jobs=-1, inplace=False)

        results_df.to_csv('{}/{}_categorical_association_results.txt'.format(output_folder,k), sep='\t', index=False) # This was made wide form to insert into .var of the program anndata.
        posthoc_df.to_csv('{}/{}_categorical_association_posthoc.txt'.format(output_folder,k), sep='\t', index=False)

    
        # Run perturbation assocation
        for samp in mdata['rna'].obs['sample'].unique():
            mdata_ = mdata[mdata['rna'].obs['sample']==samp]
            test_stats_df = compute_perturbation_association(mdata_, prog_key='cNMF', 
                                                            collapse_targets=True,
                                                            pseudobulk=False,
                                                            reference_targets=('non-targeting'),
                                                            n_jobs=-1, inplace=False)

            test_stats_df.to_csv('{}/{}_perturbation_association_results_{}.txt'.format(folder,k,samp), sep='\t', index=False)
    

        # Gene-set enrichment
        pre_res = compute_geneset_enrichment(mdata, prog_key='cNMF', data_key='rna', prog_nam=None,
                                            organism='human', library='Reactome_2022', method="fisher",
                                            database='enrichr', loading_rank_thresh=300, n_jobs=-1, 
                                            inplace=False, user_geneset=None)
        pre_res.to_csv('{}/{}_geneset_enrichment.txt'.format(output_folder,k), sep='\t', index=False)

        # GO Term enrichment
        pre_res = compute_geneset_enrichment(mdata, prog_key='cNMF', data_key='rna', prog_nam=None,
                                            organism='human', library='GO_Biological_Process_2023', method="fisher",
                                            database='enrichr', loading_rank_thresh=300, n_jobs=-1, 
                                            inplace=False, user_geneset=None)
        pre_res.to_csv('{}/{}_GO_term_enrichment.txt'.format(output_folder,k), sep='\t', index=False)

        # Run trait enrichment
        pre_res_trait = compute_trait_enrichment(mdata, gwas_data='/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz', 
                                                prog_key='cNMF', prog_nam=None, data_key='rna', 
                                                library='OT_GWAS', n_jobs=-1, inplace=False, 
                                                key_column='trait_efos', gene_column='gene_name', 
                                                method='fisher', loading_rank_thresh=300)
        pre_res_trait.to_csv('{}/{}_trait_enrichment.txt'.format(output_folder,k), sep='\t', index=False)

    