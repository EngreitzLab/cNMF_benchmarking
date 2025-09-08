# sk-cNMF

* Individual NMF inference using: [sklearn.decomposition.non_negative_factorization](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.non_negative_factorization.html)
* consensus NMF using [sk-cNMF](https://github.com/EngreitzLab/sk_cNMF) which is a slightly modified version from the [Orginal cNMF](https://github.com/dylkot/cNMF/tree/main) with more flexiblity to choose solver and loss function. 
* To run sk-cNMF, create a new conda enviornment, then run `pip install git+https://github.com/EngreitzLab/sk_cNMF.git` in the terminal


    | Parameter | Type | Default | Description |
    |-----------|------|---------|-------------|
    | counts_fn | str | - | Path to input counts matrix. If extension is .h5ad, .mtx, mtx.gz, or .npz, data is loaded accordingly. Otherwise assumed to be tab-delimited text. If .mtx or .mtx.gz, assumed to be in 10x-Genomics-formatted mtx directory. |
    | run_name | str | - | Name of the current analyisis. |
    | components | list or numpy array | - | Values of K to run NMF for |
    | n_iter | int | 100 | Number of iterations for factorization. If several k are specified, this many iterations will be run for each value of k. |
    | densify | bool | False | Convert sparse data to dense |
    | tpm_fn | str or None | None | If provided, load TPM data from file. Otherwise will compute it from the counts file |
    | seed | int or None | None | Seed for sklearn random state |
    | beta_loss | str or None | 'frobenius' | Beta loss function for NMF |
    | num_highvar_genes | int or None | 2000 | Number of highly variable genes to use for factorization if genes_file is None |
    | genes_file | str or None | None | Load high-variance genes from a list file |
    | alpha_usage | float | 0.0 | Regularization parameter for NMF corresponding to alpha_W in scikit-learn |
    | alpha_spectra | float | 0.0 | Regularization parameter for NMF corresponding to alpha_H in scikit-learn |
    | max_NMF_iter | int | 1000 | Maximum number of iterations per individual NMF run |
    | sel_thresh | Threshold for filtering NMF runs during consensus step |

