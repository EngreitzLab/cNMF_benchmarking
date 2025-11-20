import os
import sys
import yaml
import logging
import mudata
import muon as mu
import pandas as pd
import numpy as np
import argparse
import mygene
import cnmf
import scanpy as sc


# Change path to wherever you have repo locally
sys.path.append('/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline')

from Evaluation.src import (
    compute_categorical_association,
    compute_geneset_enrichment,
    compute_trait_enrichment,
    compute_perturbation_association,
    compute_explained_variance,
    compute_motif_enrichment
)


from Evaluation.src import (
    rename_gene_mygene,
    rename_gene_dictionary,
    check_evaluation_pipeline_format,
    check_gene_names
)

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
    parser.add_argument('--K', nargs='*', type=int, default=None) # allow zero input 
    parser.add_argument('--sel_thresh', nargs='*', type=int, default=None) # allow zero input 

    # running different tests
    parser.add_argument('--Perform_categorical', action="store_true")
    parser.add_argument('--Perform_perturbation', action="store_true")
    parser.add_argument('--Perform_geneset', action="store_true")
    parser.add_argument('--Perform_trait', action="store_true")
    parser.add_argument('--Perform_explained_variance', action="store_true")
    parser.add_argument('--Perform_motif', action="store_true")

    # resourses 
    parser.add_argument('--X_normalized_path', type=str,  help='path to normalized input cell x gene matrix from cNMF pipeline', required=True)
    parser.add_argument('--mdata_guide_path', type=str,  help='path to mdata with correct guide assignments', required=True)
    parser.add_argument('--guide_annotation_path', type=str,  help='path to file with guide informations', required=True)
    parser.add_argument('--gwas_data_path', type=str,  help='path to file with gwas information for trait enrichment test', required=True)
    parser.add_argument('--reference_gtf_path', type=str,  help='path to reference GTF file for checking gene names')

    # keys
    parser.add_argument('--data_key', help='access gene expression in mdata',type=str, default="rna") 
    parser.add_argument('--prog_key', help='access cNMF program in mdata',type=str,  default="cNMF") 
    parser.add_argument('--categorical_key', help='access cell condition in obs',type=str, default="sample")  
    parser.add_argument('--organism', help='data species',type=str, default="human")  


    args = parser.parse_args()


    # either change the array here or run each component in parallel
    if args.K is None:
        args.K = [30, 50, 60, 80, 100, 200, 250, 300]

    if args.sel_thresh is None:
        args.sel_thresh = [0.4, 0.8, 2.0]


    # create output directory       
    os.makedirs((f'{args.out_dir}/{args.run_name}/Eval'), exist_ok=True)

    # save comfigs used         
    args_dict = vars(args)
    job_id = os.environ.get('SLURM_JOB_ID')
    with open(f'{args.out_dir}/{args.run_name}/Eval/config_{job_id}.yml', 'w') as f:
        yaml.dump(args_dict, f, default_flow_style=False, width=1000)


    # defind objects used in explained variance
    cnmf_obj = cnmf.cNMF(output_dir=args.out_dir, name=args.run_name)
    X_norm = sc.read_h5ad(args.X_normalized_path)
    X = X_norm.X

    # list of non-targeting guides
    df_target = pd.read_csv(args.guide_annotation_path, sep = "\t", index_col = 0)
    df_target_non = df_target[df_target["targeting"] == False]
    reference_targets = df_target_non.index.values.tolist()

    # read guide
    mdata_guide = mu.read(args.mdata_guide_path)

    for sel_thresh in args.sel_thresh:
        for k in args.K:  
            
            output_folder = f"{args.out_dir}/{args.run_name}/Eval/{k}_{str(sel_thresh).replace('.','_')}"

            os.makedirs(output_folder, exist_ok=True)

            # Load mdata
            mdata = mu.read('{out_dir}/{run_name}/adata/cNMF_{k}_{sel_thresh}.h5mu'.format(out_dir = args.out_dir,
                                                                                    run_name =args.run_name,
                                                                                    k=k,
                                                                                    sel_thresh = str(sel_thresh).replace('.','_'))) 
            # assign guide
            _assign_guide(mdata, mdata_guide)   

            # checks for correct mdata format
            if not check_evaluation_pipeline_format(mdata,prog_key=args.prog_key):
                raise ValueError("mdata format is incorrect")

            valid = check_gene_names(mdata,prog_key=args.prog_key, data_key=args.data_key,categorical_key=args.categorical_key,reference_gtf_path=args.reference_gtf_path)
            if not  valid["is_valid"]:
                raise ValueError("mdata gene naming is incorrect")


            # Run categorical assocation
            if args.Perform_categorical: 
                results_df, posthoc_df = compute_categorical_association(mdata, prog_key='cNMF', categorical_key=args.categorical_key, 
                                                                        pseudobulk_key=None, test='dunn', n_jobs=-1, inplace=False)

                results_df.to_csv('{}/{}_categorical_association_results.txt'.format(output_folder,k), sep='\t', index=False) # This was made wide form to insert into .var of the program anndata.
                posthoc_df.to_csv('{}/{}_categorical_association_posthoc.txt'.format(output_folder,k), sep='\t', index=False)


            # Run perturbation assocation
            if args.Perform_perturbation: 
                for samp in mdata['rna'].obs[args.categorical_key].unique():
                    mdata_ = mdata[mdata['rna'].obs[args.categorical_key]==samp]
                    test_stats_df = compute_perturbation_association(mdata_, prog_key=args.prog_key, 
                                                                    collapse_targets=True,
                                                                    pseudobulk=False,
                                                                    reference_targets=reference_targets,
                                                                    n_jobs=-1, inplace=False)

                    test_stats_df.to_csv('{}/{}_perturbation_association_results_{}.txt'.format(output_folder,k,samp), sep='\t', index=False)
        

            # Gene-set enrichment
            if args.Perform_geneset:
                pre_res = compute_geneset_enrichment(mdata, prog_key=args.prog_key, data_key=args.data_key, prog_nam=None,
                                                    organism=args.organism, library='Reactome_2022', method="fisher",
                                                    database='enrichr', loading_rank_thresh=300, n_jobs=-1, 
                                                    inplace=False, user_geneset=None)
                pre_res.to_csv('{}/{}_geneset_enrichment.txt'.format(output_folder,k), sep='\t', index=False)


                # GO Term enrichment
                pre_res = compute_geneset_enrichment(mdata, prog_key=args.prog_key, data_key=args.data_key, prog_nam=None,
                                                    organism=args.organism, library='GO_Biological_Process_2023', method="fisher",
                                                    database='enrichr', loading_rank_thresh=300, n_jobs=-1, 
                                                    inplace=False, user_geneset=None)
                pre_res.to_csv('{}/{}_GO_term_enrichment.txt'.format(output_folder,k), sep='\t', index=False)


            # Run trait enrichment
            if args.Perform_trait:
                pre_res_trait = compute_trait_enrichment(mdata, gwas_data=args.gwas_data_path, 
                                                        prog_key=args.prog_key, prog_nam=None, data_key=args.data_key, 
                                                        library='OT_GWAS', n_jobs=-1, inplace=False, 
                                                        key_column='trait_efos', gene_column='gene_name', 
                                                        method='fisher', loading_rank_thresh=300)
                pre_res_trait.to_csv('{}/{}_trait_enrichment.txt'.format(output_folder,k), sep='\t', index=False)


            # Run motif analysis 
            if args.Perform_motif:

                fimo_thresh_enhancer = 1e-6
                fimo_thresh_promoter = 1e-4

                for samp in mdata['rna'].obs[args.categorical_key].unique():
                    for class_, thresh in [('enhancer', fimo_thresh_enhancer), 
                                        ('promoter', fimo_thresh_promoter)]:

                    
                    
                        loci_file = '/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Resources/scE2G_links/EnhancerPredictionsAllPutative.ForVariantOverlap.shrunk150bp_D{}_{}.tsv'.format(i, class_)
                        motif_match_df, motif_count_df, motif_enrichment_df = compute_motif_enrichment(
                            mdata, 
                            prog_key='cNMF',
                            data_key='rna',
                            motif_file='/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Resources/hocomoco_meme.meme',
                            seq_file='/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Resources/hg38.fa',
                            loci_file=loci_file,
                            window=1000,
                            sig=thresh,
                            eps=1e-4,
                            n_top=2000,
                            n_jobs=-1,
                            inplace=False
                        )

                        motif_match_df.to_csv(os.path.join(args.out_dir, f'cNMF_{class_}_pearson_topn2000_{samp}_motif_match.txt'), sep='\t', index=False)
                        motif_count_df.to_csv(os.path.join(args.out_dir, f'cNMF_{class_}_pearson_topn2000_{samp}_motif_count.txt'), sep='\t', index=False)
                        motif_enrichment_df.to_csv(os.path.join(args.out_dir, f'cNMF_{class_}_pearson_topn2000_{samp}_motif_enrichment.txt'), sep='\t', index=False)

            # Run explained variance
            if args.Perform_explained_variance:
                compute_explained_variance(cnmf_obj, X, k, output_folder = output_folder)

