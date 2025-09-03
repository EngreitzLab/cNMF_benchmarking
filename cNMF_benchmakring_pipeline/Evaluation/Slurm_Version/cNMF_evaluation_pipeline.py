import os
import sys
import yaml
import logging
import mudata
import pandas as pd
import argparse

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
    parser.add_argument('--thre', default = '2.0')

    # running different tests
    parser.add_argument('--Perform_categorical', type=bool, default = True)
    parser.add_argument('--Perform_perturbation', type=bool, default = True)
    parser.add_argument('--Perform_geneset', type=bool, default = True)
    parser.add_argument('--Perform_trait', type=bool, default = True)
    parser.add_argument('--Perform_explained_variance', type=bool, default = True)
    parser.add_argument('--Perform_motif', type=bool, default = True)


    args = parser.parse_args()


    for k in args.K:  

        os.makedirs(f"{args.input_folder}/Eval/{k}", exist_ok=True)
        output_folder = f"{args.input_folder}/Eval/{k}"

        # Load mdata
        mdata = mudata.read('{}/adata/cNMF_{}_{}.h5mu'.format(args.input_folder,k,args.thre))
        

        # Run categorical assocation
        if args.Perform_categorical: 
            results_df, posthoc_df = compute_categorical_association(mdata, prog_key='cNMF', categorical_key='sample', 
                                                                    pseudobulk_key=None, test='dunn', n_jobs=-1, inplace=False)

            results_df.to_csv('{}/{}_categorical_association_results.txt'.format(output_folder,k), sep='\t', index=False) # This was made wide form to insert into .var of the program anndata.
            posthoc_df.to_csv('{}/{}_categorical_association_posthoc.txt'.format(output_folder,k), sep='\t', index=False)


        # Run perturbation assocation
        if args.Perform_perturbation: 

            assign_guide(mdata, args.mdata_guide) 

            for samp in mdata['rna'].obs['sample'].unique():
                mdata_ = mdata[mdata['rna'].obs['sample']==samp]
                test_stats_df = compute_perturbation_association(mdata_, prog_key='cNMF', 
                                                                collapse_targets=True,
                                                                pseudobulk=False,
                                                                reference_targets=('non-targeting'),
                                                                n_jobs=-1, inplace=False)

                test_stats_df.to_csv('{}/{}_perturbation_association_results_{}.txt'.format(output_folder,k,samp), sep='\t', index=False)
    

        # Gene-set enrichment
        if args.Perform_geneset:
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
        if args.Perform_trait:
            pre_res_trait = compute_trait_enrichment(mdata, gwas_data='/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz', 
                                                    prog_key='cNMF', prog_nam=None, data_key='rna', 
                                                    library='OT_GWAS', n_jobs=-1, inplace=False, 
                                                    key_column='trait_efos', gene_column='gene_name', 
                                                    method='fisher', loading_rank_thresh=300)
            pre_res_trait.to_csv('{}/{}_trait_enrichment.txt'.format(output_folder,k), sep='\t', index=False)

        
        # Run explained variance
        if args.Perform_explained_variance:
            pass


        # Run motif analysis 
        if args.Perform_motif:
            fimo_thresh_enhancer = 1e-6
            fimo_thresh_promoter = 1e-4

            for i in range(4):
                for class_, thresh in [('enhancer', fimo_thresh_enhancer), 
                                    ('promoter', fimo_thresh_promoter)]:

                    loci_file = '/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Evaluation/Resources/scE2G_links/EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp_D{}_{}.tsv'.format(i, class_)
                    motif_match_df, motif_count_df, motif_enrichment_df = compute_motif_enrichment(
                        mdata, 
                        prog_key='cNMF',
                        data_key='rna',
                        motif_file='/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Evaluation/Resources/hocomoco_meme.meme',
                        seq_file='/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmakring_pipeline/Evaluation/Resources/hg38.fa',
                        loci_file=loci_file,
                        window=1000,
                        sig=thresh,
                        eps=1e-4,
                        n_top=2000,
                        n_jobs=-1,
                        inplace=False
                    )

                    motif_match_df.to_csv(os.path.join(output_folder, f'cNMF_{class_}_pearson_topn2000_sample_D{i}_motif_match.txt'), sep='\t', index=False)
                    motif_count_df.to_csv(os.path.join(output_folder, f'cNMF_{class_}_pearson_topn2000_sample_D{i}_motif_count.txt'), sep='\t', index=False)
                    motif_enrichment_df.to_csv(os.path.join(output_folder, f'cNMF_{class_}_pearson_topn2000_sample_D{i}_motif_enrichment.txt'), sep='\t', index=False)


    