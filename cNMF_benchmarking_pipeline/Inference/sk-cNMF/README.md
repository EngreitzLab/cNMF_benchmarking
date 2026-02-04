# sk-cNMF

* Individual NMF inference using: [sklearn.decomposition.non_negative_factorization](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.non_negative_factorization.html)
* consensus NMF using [sk-cNMF](https://github.com/EngreitzLab/sk_cNMF) which is a slightly modified version from the [Orginal cNMF](https://github.com/dylkot/cNMF/tree/main) with more flexiblity to choose solver and loss function. 
* To run sk-cNMF, create a new conda environment with `conda env create -f environment.yml --name sk-cNMF` with the provided yml file, then run `pip install git+https://github.com/EngreitzLab/sk_cNMF.git` in the terminal

## cNMF Parameters 


| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| counts_fn | str | - | Path to input counts matrix. If extension is .h5ad, .mtx, mtx.gz, or .npz, data is loaded accordingly. Otherwise assumed to be tab-delimited text. If .mtx or .mtx.gz, assumed to be in 10x-Genomics-formatted mtx directory. |
| run_name | str | - | Name of the current analyisis. |
| components | list or numpy array | - | Values of K to run NMF for |
| n_iter | int | 100 | Number of iterations for factorization. If several k are specified, this many iterations will be run for each value of k. |
| densify | bool | False | Convert sparse data to dense |
| tpm_fn | str or None | None | If provided, load TPM data from file. Otherwise will compute it from the counts file |
| seed | int or None | None | Seed for sklearn random state |
| beta_loss | str or float | "frobenius" | Beta loss metric for approximation. Options: "frobenius" (L2), "kullback-leibler" (KL), "itakura-saito" (IS), or float value |
| num_highvar_genes | int or None | 2000 | Number of highly variable genes to use for factorization if genes_file is None |
| genes_file | str or None | None | Load high-variance genes from a list file |
| alpha_usage | float | 0.0 | Regularization parameter for NMF corresponding to alpha_W in scikit-learn |
| alpha_spectra | float | 0.0 | Regularization parameter for NMF corresponding to alpha_H in scikit-learn |
| max_NMF_iter | int | 1000 | Maximum number of iterations per individual NMF run |
| algo | str | "halsvar" | Algorithm choice: "mu", "halsvar" |
| sel_thresh | int | - | Threshold for filtering NMF runs during consensus step|
| tol | float | 1e-4 | Tolerance for convergence check |

## Additional Parameters for Annotation and Compilation

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| species | str | - | Species for gene annotation (required) |
| check_format | bool | False | Validate input data format against reference files |
| run_refit_only | bool | False | Skip factorization and only run consensus/refit steps |
| parallel_running | bool | False | Compile files after parallel mode for multiple K values |
| num_gene | int | 300 | Number of top genes to include in annotation |
| guide_annotation_path | str | None | Path to file containing guide information |
| reference_gtf_path | str | None | Path to reference GTF file for gene name validation |
| data_key | str | "rna" | Key to access gene expression data in MuData object |
| prog_key | str | "cNMF" | Key to access cNMF programs in MuData object |
| categorical_key | str | "sample" | Key to access cell condition information in obs |
| guide_names_key | str | "guide_names" | Key to access guide names in uns |
| guide_targets_key | str | "guide_targets" | Key to access guide targets in uns |
| guide_assignment_key | str | "guide_assignment_key" | Key to access guide assignments in obsm |
| nmf_seeds_path | str | None |path to .npy file containing custom NMF seeds for reproducibility |


