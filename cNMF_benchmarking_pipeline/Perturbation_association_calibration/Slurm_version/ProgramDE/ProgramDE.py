from src.sceptre import (
    add_burden_covariate,
    build_ntc_group_inputs,
    compute_ntc_group_null_pvals_parallel,
    crt_pvals_for_ntc_groups_ensemble,
    crt_pvals_for_ntc_groups_ensemble_skew,
    make_ntc_groups_ensemble,
    prepare_crt_inputs,
    run_all_genes_union_crt,
)
from src.visualization import qq_plot_ntc_pvals
import scanpy as sc
import argparse
import numpy as np


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()

    parser.add_argument('--mdata_path', help='adata object', type=str, required=True)  
    parser.add_argument('--n_jobs', help='jobs for parallelism', type=int, default=-1)  
    parser.add_argument('--n_permutations', help='number of permutations', type=int, default=1023)
    parser.add_argument('--group_size', help='number of guides to group', type=int, default=6)

    parser.add_argument('--covariates', help='lists of covariates to include', nargs = '+', type = str, required = True)  
    parser.add_argument('--NTC_guides', help='names of non-targeting guides', nargs = '+', type = str, default = ["SAFE", "non-targeting", "NTC"])



    args = parser.parse_args()


    # take adata from mdata 
    mdata = sc.read(args.mdata_path)
    adata = mdata["cNMF"].copy()
    adata.obsm["cnmf_usage"] = np.asarray(adata.X)  # ensure dense float

    # Program names
    adata.uns["program_names"] = list(adata.var_names)

    # Guide names must match columns of guide_assignment
    guide_names = list(mdata["cNMF"].uns["guide_names"])
    adata.uns["guide_names"] = guide_names

    # guide2gene must map guide_name -> gene name (keys must match guide_names)
    guide2gene = (
        mdata["cNMF"].uns["guide_id"]
        .loc[guide_names]
        .to_dict()
    )
    adata.uns["guide2gene"] = guide2gene


    # Make covariates a DataFrame 
    covar_df = mdata['rna'].obs.copy()
    missing_covariates = set(args.covariates) - set(covar_df.columns)

    if missing_covariates:
        raise ValueError(f"The following covariates are not found in the DataFrame columns: {missing_covariates}")

    adata.obsm["covar"] = covar_df



    # 1) when cNMF have too many zeros, this can fix log transformation error 
    U = adata.obsm["cnmf_usage"].copy()
    U = np.maximum(U, 1e-8)
    U /= U.sum(axis=1, keepdims=True)
    adata.obsm["cnmf_usage"] = U


    # 2) Prepare CRT inputs
    inputs = prepare_crt_inputs(
        adata=adata, 
        usage_key="cnmf_usage",
        covar_key="covar",
        guide_assignment_key="guide_assignment",
        guide2gene_key="guide2gene")

    # 3) Run gene-level CRT with stratified permutation
    out = run_all_genes_union_crt(
        inputs=inputs,
        B=args.n_permutations,
        n_jobs=args.n_jobs,
        calibrate_skew_normal=True,
        return_raw_pvals=True,
        return_skew_normal=True,
    )

    # 4) Build NTC guide groups (6-guide units) + compute NTC p-values
    ntc_guides, guide_freq, guide_to_bin, real_sigs = build_ntc_group_inputs(
        inputs=inputs,
        ntc_label=args.NTC_guides,
        group_size=args.group_size,
        n_bins=10,
    )

    ntc_groups_ens = make_ntc_groups_ensemble(
        ntc_guides=ntc_guides,
        ntc_freq=guide_freq,
        real_gene_bin_sigs=real_sigs,
        guide_to_bin=guide_to_bin,
        n_ensemble=10,
        seed0=7,
        group_size=args.group_size,
        max_groups=None,
    )

    ntc_group_pvals_ens = crt_pvals_for_ntc_groups_ensemble(
        inputs=inputs,
        ntc_groups_ens=ntc_groups_ens,
        B=args.n_permutations,
        seed0=23,
        # If you ran the main CRT with test_stat="utest", pass the same here.
        # test_stat="utest",
        # test_stat_kwargs={"use": "clr", "rank_method": "average"},
    )

    ntc_group_pvals_skew_ens = crt_pvals_for_ntc_groups_ensemble_skew(
        inputs=inputs,
        ntc_groups_ens=ntc_groups_ens,
        B=args.n_permutations,
        seed0=23,
    )

    # 5) CRT-null p-values matched to the same NTC units
    null_pvals = compute_ntc_group_null_pvals_parallel(
        inputs=inputs,
        ntc_groups_ens=ntc_groups_ens,
        B=args.n_permutations,
        n_jobs=args.n_jobs,
        backend="threading",
        # Keep test_stat/test_stat_kwargs in sync with your main run if needed.
        # test_stat="utest",
        # test_stat_kwargs={"use": "clr", "rank_method": "average"},
    )

    # 6) QQ plot
    ax = qq_plot_ntc_pvals(
        pvals_raw_df=out["pvals_raw_df"],
        guide2gene=args.adata.uns["guide2gene"],
        ntc_genes=args.NTC_guides,
        pvals_skew_df=out["pvals_df"],
        null_pvals=null_pvals,
        ntc_group_pvals_ens=ntc_group_pvals_ens,
        ntc_group_pvals_skew_ens=ntc_group_pvals_skew_ens,
        show_ntc_ensemble_band=True,
        show_all_pvals=True,
        title="QQ plot: S-CRT (grouped NTC controls)",
    )