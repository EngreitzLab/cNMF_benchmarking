# torch-cNMF

* GPU powered cNMF
* Individual NMF inference using: [NMF-Torch](https://github.com/lilab-bcb/nmf-torch)
* consensus NMF using: [torch-cNMF](https://github.com/ymo6/torch_based_cNMF) 
* To run torch-cNMF, create a new conda environment with `conda env create -f environment.yml --name torch-cNMF` with the provided yml file, then run `pip install git+https://github.com/ymo6/torch_based_cNMF.git` and `pip install git+https://github.com/ymo6/nmf-torch.git` in the terminal
* Singularity Container can be found in https://hub.docker.com/r/igvf/torch-cnmf/tags

## Required I/O Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| counts_fn | str | Path to input counts matrix (.h5ad, .mtx, .mtx.gz, .npz, or tab-delimited text) |
| output_directory | str | Path to output directory where results will be saved |
| run_name | str | Name for this cNMF run (used for output file naming) |

## cNMF Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| K | list of int | [30, 50, 60, 80, 100, 200, 250, 300] | Values of K (number of components) to run NMF for |
| numiter | int | 10 | Number of NMF iterations per K value |
| densify | flag | False | Densify sparse matrix before factorization |
| tpm_fn | str | None | Path to TPM normalized data (optional) |
| seed | int | 14 | Random seed for reproducibility |
| loss | str | "frobenius" | Loss function: "frobenius" (L2), "kullback-leibler" (KL), "itakura-saito" (IS), or float |
| numhvgenes | int | 5451 | Number of highly variable genes to use for factorization |
| genes_file | str | None | Path to file containing list of genes to use instead of HVG selection |
| alpha_usage | float | 0.0 | Regularization parameter for usage matrix (alpha_W) |
| alpha_spectra | float | 0.0 | Regularization parameter for spectra matrix (alpha_H) |
| use_gpu | flag | False | Use GPU acceleration if available |
| mode | str | "batch" | Learning mode: "batch", "online", or "dataloader". Online/dataloader only works when beta=2.0 |
| algo | str | "mu" | Algorithm choice: "mu" (multiplicative update), "halsvar" |
| init | str | "random" | Initialization method: "random", "nndsvd", "nndsvda", "nndsvdar" |
| tol | float | 1e-4 | Tolerance for convergence check |
| n_jobs | int | -1 | Number of CPU threads (-1 uses all available cores) |
| l1_ratio_usage | float | 0.0 | L1 ratio for usage regularization (0=L2 only, 1=L1 only) |
| l1_ratio_spectra | float | 0.0 | L1 ratio for spectra regularization (0=L2 only, 1=L1 only) |
| fp_precision | str | "float" | Numeric precision: "float" (32-bit) or "double" (64-bit) |
| batch_max_iter | int | 500 | Maximum iterations for batch learning |
| batch_hals_tol | float | 0.05 | Tolerance for HALS - maximal relative change threshold |
| batch_hals_max_iter | int | 200 | Maximum HALS iterations. Set to 1 for standard HALS |
| online_max_pass | int | 20 | Maximum number of online passes through all data |
| online_chunk_size | int | 5000 | Chunk/mini-batch size for online learning |
| online_chunk_max_iter | int | 200 | Maximum iterations per chunk in online mode |
| online_usage_tol | float | 0.05 | Tolerance for usage updates in online learning |
| online_spectra_tol | float | 0.05 | Tolerance for spectra updates in online learning |
| shuffle_cells | flag | False | Shuffle cells before factorization (recommended for online learning) |
| sk_cd_refit | flag | False | Use scikit-learn coordinate descent for refitting |
| sel_thresh | list of float | [0.4, 0.8, 2.0] | Density threshold(s) for consensus selection |
| nmf_seeds_path | str | None | Path to .npy file containing custom NMF seeds for reproducibility |

## Preprocessing Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| remove_noncoding | flag | False | Remove non-coding genes whose symbol starts with an Ensembl ID prefix before factorization |
| ensembl_prefix | str | "ENSG" | Ensembl ID prefix used to identify non-coding genes |
| gene_symbol_key | str | "symbol" | Column in adata.var containing gene symbols |

## Annotation and Compilation Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| species | str | - | Species for gene annotation (required) |
| check_format | flag | False | Validate input data format before running |
| parallel_running | flag | False | Compile files after parallel mode for multiple K values |
| num_gene | int | 300 | Number of top genes to include in annotation |
| run_refit | flag | False | Run the refit step (combine, k_selection_plot, and consensus) |
| run_complie_annotation | flag | False | Run the compilation and annotation step |
| run_factorize | flag | False | Run the NMF factorization step |
| guide_annotation_path | str | None | Path to file containing guide information |
| reference_gtf_path | str | None | Path to reference GTF file for gene name validation |

## Data Access Keys

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| data_key | str | "rna" | Key to access gene expression data in MuData object |
| prog_key | str | "cNMF" | Key to access cNMF programs in MuData object |
| categorical_key | str | "sample" | Key to access cell condition information in obs |
| guide_names_key | str | "guide_names" | Key to access guide names in uns |
| guide_targets_key | str | "guide_targets" | Key to access guide targets in uns |
| guide_assignment_key | str | "guide_assignment" | Key to access guide assignments in obsm |
