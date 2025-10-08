# torch-cNMF

* GPU powered cNMF
* Individual NMF inference using: [NMF-Torch](https://github.com/lilab-bcb/nmf-torch)
* consensus NMF using: [torch-cNMF](https://github.com/ymo6/torch_based_cNMF) 
* To run torch-cNMF, create a new conda environment with `conda env create -f environment.yml --name torch-cNMF` with the provided yml file, then run `pip install git+https://github.com/ymo6/torch_based_cNMF.git` in the terminal


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
    | alpha_usage | float | 0.0 | Regularization parameter for NMF corresponding to alpha_W |
    | alpha_spectra | float | 0.0 | Regularization parameter for NMF corresponding to alpha_H |
    | use_gpu | bool | False | Whether to use GPU |
    | mode | str | "batch" | Learning mode: "batch" or "online". Online only works when beta=2.0 |
    | batch_size | int | 5000 | Batch size for online NMF learning |
    | algo | str | "halsvar" | Algorithm choice: "mu", "halsvar" |
    | init | str | "nndsvdar" | Initialization method: "random", "nndsvd", "nndsvda", "nndsvdar" |
    | tol | float | 1e-4 | Tolerance used for convergence check |
    | n_jobs | int | -1 | Number of CPU threads to use. If -1, use PyTorch's default |
    | l1_ratio_usage | float | 0.0 | Ratio of L1 penalty on W (0-1). L2 penalty ratio is (1 - l1_ratio_usage) |
    | l1_ratio_spectra | float | 0.0 | Ratio of L1 penalty on H (0-1). L2 penalty ratio is (1 - l1_ratio_spectra) |
    | fp_precision | str | "float" | Numeric precision: "float" (torch.float) or "double" (torch.double) |
    | batch_max_iter | int | 500 | Maximum iterations for batch learning |
    | batch_hals_tol | float | 0.05 | Tolerance for HALS to mimic BPP - maximal relative change threshold |
    | batch_hals_max_iter | int | 200 | Maximum iterations for updating H & W to mimic BPP. Set to 1 for standard HALS |
    | online_max_pass | int | 20 | Maximum number of online passes through all data |
    | online_chunk_size | int | 5000 | Chunk/mini-batch size for online learning |
    | online_chunk_max_iter | int | 200 | Maximum iterations for updating H or W in online learning |
    | online_usage_tol | float | 0.05 | Tolerance for updating W in each chunk during online learning |
    | shuffle_cells | bool | False | Shuffle cells in obs and guide assignment in obsm to do online learning is recommanded |
    | sel_thresh | int | - | Threshold for filtering NMF runs during consensus step|

