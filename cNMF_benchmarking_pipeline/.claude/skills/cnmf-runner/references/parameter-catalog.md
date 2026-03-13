# Parameter Catalog

Complete parameter reference for all pipeline stages, extracted from argparse definitions.

**Note**: Parameter names with typos (e.g., `--tagert_col_name`, `--run_complie_annotation`, `--top_enrichned_term`, `--down_thred_log`) are intentional -- they match the actual argparse definitions in the pipeline scripts and must be used as-is.

---

## 1. sk-cNMF Inference

**Script**: `Inference/sk-cNMF/Slurm_Version/sk-cNMF_batch_inference_pipeline.py`
**Conda**: `sk-cNMF`

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| `--counts_fn` | str | Path to input counts file (.h5ad or .h5mu) |
| `--output_directory` | str | Directory for all outputs |
| `--run_name` | str | Name for this cNMF run |
| `--species` | str | Species for gene annotation (`human`, `mouse`) |

### Core cNMF Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--K` | int (nargs=\*) | `[30,50,60,80,100,200,250,300]` | K values (number of components) |
| `--numiter` | int | `10` | Number of NMF replicates |
| `--numhvgenes` | int | `5451` | Number of highly variable genes |
| `--sel_thresh` | float (nargs=\*) | `[2.0]` | Density thresholds for consensus filtering |
| `--seed` | int | `14` | Random seed |
| `--init` | str | `random` | NMF initialization (`random`, `nndsvd`, `nndsvda`, `nndsvdar`) |
| `--loss` | str | `frobenius` | NMF loss function |
| `--algo` | str | `mu` | NMF algorithm (`mu` or `cd`) |
| `--max_NMF_iter` | int | `500` | Maximum NMF iterations |
| `--tol` | float | `1e4` | Convergence tolerance (**WARNING**: likely a bug, should be `1e-4` -- see torch-cNMF) |

### Workflow Flags

| Flag | Description |
|------|-------------|
| `--run_factorize` | Run the NMF factorization step |
| `--run_refit` | Run combine, k_selection_plot, and consensus steps |
| `--run_complie_annotation` | Compile results and generate gene annotations |
| `--check_format` | Validate input data format before running |
| `--parallel_running` | Enable parallel processing for multiple K values |

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--nmf_seeds_path` | str | None | Path to .npy file with custom NMF seeds |
| `--num_gene` | int | `300` | Top genes for annotation |
| `--guide_annotation_path` | str | None | Guide annotation TSV path |
| `--reference_gtf_path` | str | None | Reference GTF path |

### Keys

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--data_key` | `rna` | RNA modality key in MuData |
| `--prog_key` | `cNMF` | cNMF modality key in MuData |
| `--categorical_key` | `sample` | Categorical variable key in obs |
| `--guide_names_key` | `guide_names` | Guide names key in uns |
| `--guide_targets_key` | `guide_targets` | Guide targets key in uns |
| `--guide_assignment_key` | `guide_assignment_key` | Guide assignment key in obsm (**NOTE**: default is the literal string `"guide_assignment_key"`, not `"guide_assignment"` -- this differs from torch-cNMF) |

---

## 2. torch-cNMF Inference

**Script**: `Inference/torch-cNMF/Slurm_Version/torch-cNMF_inference_pipeline.py`
**Conda**: `torch-cNMF`

torch-cNMF shares most parameters with sk-cNMF but has these important differences:
- Does NOT have `--max_NMF_iter` (use `--batch_max_iter` instead)
- `--sel_thresh` default is `None` (code fallback: `[0.4, 0.8, 2.0]`)
- `--tol` default is `1e-4` (correct, unlike sk-cNMF's `1e4`)
- `--guide_assignment_key` default is `guide_assignment` (not `guide_assignment_key`)

### Required Parameters

Same as sk-cNMF: `--counts_fn`, `--output_directory`, `--run_name`, `--species`

### Core cNMF Parameters (shared with sk-cNMF, same defaults unless noted)

| Parameter | Type | Default | Notes |
|-----------|------|---------|-------|
| `--K` | int (nargs=\*) | `[30,50,60,80,100,200,250,300]` | |
| `--numiter` | int | `10` | |
| `--numhvgenes` | int | `5451` | |
| `--sel_thresh` | float (nargs=\*) | `[0.4,0.8,2.0]` | **Different from sk-cNMF** |
| `--seed` | int | `14` | |
| `--init` | str | `random` | |
| `--loss` | str | `frobenius` | |
| `--algo` | str | `mu` | |
| `--tol` | float | `1e-4` | **Different from sk-cNMF** |

### Workflow Flags

Same as sk-cNMF: `--run_factorize`, `--run_refit`, `--run_complie_annotation`, `--check_format`, `--parallel_running`

### Optional Parameters (shared with sk-cNMF)

Same as sk-cNMF: `--nmf_seeds_path`, `--num_gene`, `--guide_annotation_path`, `--reference_gtf_path`

### Keys

Same as sk-cNMF except:
| Parameter | Default | Notes |
|-----------|---------|-------|
| `--guide_assignment_key` | `guide_assignment` | **Different from sk-cNMF** (which defaults to `guide_assignment_key`) |

### Additional torch-cNMF Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--mode` | str | `batch` | NMF mode: `batch`, `online`, or `dataloader` |
| `--n_jobs` | int | `-1` | Parallel jobs (-1 = all cores) |
| `--use_gpu` | flag | False | Enable GPU acceleration |
| `--densify` | flag | False | Densify sparse matrix before factorization |
| `--tpm_fn` | str | None | Path to TPM normalized data |
| `--genes_file` | str | None | Custom gene list file |

### Regularization

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--alpha_usage` | float | `0.0` | Regularization for usage matrix |
| `--alpha_spectra` | float | `0.0` | Regularization for spectra matrix |
| `--l1_ratio_usage` | float | `0.0` | L1 ratio for usage (0=L2, 1=L1) |
| `--l1_ratio_spectra` | float | `0.0` | L1 ratio for spectra |

### Batch Mode Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--batch_max_iter` | int | `500` | Max iterations for batch NMF |
| `--batch_hals_tol` | float | `0.05` | HALS tolerance in batch mode |
| `--batch_hals_max_iter` | int | `200` | Max HALS iterations |

### Online Mode Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--online_max_pass` | int | `20` | Max passes through data |
| `--online_chunk_size` | int | `5000` | Chunk size for online learning |
| `--online_chunk_max_iter` | int | `200` | Max iterations per chunk |
| `--online_usage_tol` | float | `0.05` | Usage update tolerance |
| `--online_spectra_tol` | float | `0.05` | Spectra update tolerance |

### Other torch-cNMF Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--fp_precision` | str | `float` | Floating point: `float` (32-bit) or `double` (64-bit) |
| `--shuffle_cells` | flag | False | Shuffle cells before factorization |
| `--sk_cd_refit` | flag | False | Use sklearn coordinate descent for refitting |
| `--remove_noncoding` | flag | False | Remove non-coding genes |
| `--ensembl_prefix` | str | `ENSG` | Prefix for identifying non-coding genes |
| `--gene_symbol_key` | str | `symbol` | Column in var with gene symbols |

---

## 3. Evaluation

**Script**: `Evaluation/Slurm_Version/cNMF_evaluation_pipeline.py`
**Conda**: `NMF_Benchmarking`

### Required

| Parameter | Type | Description |
|-----------|------|-------------|
| `--out_dir` | str | Directory containing cNMF output |
| `--run_name` | str | cNMF run name |
| `--gwas_data_path` | str | Path to GWAS data (OpenTargets L2G) |

### Test Flags

| Flag | Description |
|------|-------------|
| `--Perform_categorical` | Categorical association (Kruskal-Wallis + Dunn's test) |
| `--Perform_perturbation` | Perturbation association |
| `--Perform_geneset` | Gene set enrichment (Reactome + GO terms) |
| `--Perform_trait` | GWAS trait enrichment |
| `--Perform_explained_variance` | Explained variance per K |
| `--Perform_motif` | TF motif enrichment (WIP) |

### Optional

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--K` | int (nargs=\*) | `[30,50,60,80,100,200,250,300]` | K values to evaluate |
| `--sel_thresh` | float (nargs=\*) | `[0.4,0.8,2.0]` | Density thresholds |
| `--X_normalized_path` | str | None | Normalized counts h5ad (needed for explained variance) |
| `--guide_annotation_path` | str | None | Guide annotation TSV |
| `--reference_gtf_path` | str | None | Reference GTF path |
| `--data_guide_path` | str | None | MuData with additional guide info |
| `--organism` | str | `human` | Species for enrichment |
| `--FDR_method` | str | `StoreyQ` | FDR method: `StoreyQ` or `BH` |
| `--n_top` | int | `300` | Top genes for enrichment tests |
| `--check_format` | flag | | Validate format before running |
| `--guide_annotation_key` | str (nargs=\*) | `["non-targeting"]` | Non-targeting guide identifiers (accepts multiple values) |

### Keys

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--data_key` | `rna` | RNA modality key |
| `--prog_key` | `cNMF` | cNMF modality key |
| `--categorical_key` | `sample` | Categorical variable key |
| `--guide_names_key` | `guide_names` | Guide names key |
| `--guide_targets_key` | `guide_targets` | Guide targets key |
| `--guide_assignment_key` | `guide_assignment` | Guide assignment key |

---

## 4. K-Selection Plot

**Script**: `Plotting/Slurm_Version/cNMF_k_selection.py`
**Conda**: `torch-cNMF`

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--output_directory` | str | — | Yes | Directory with cNMF output |
| `--run_name` | str | — | Yes | cNMF run name |
| `--save_folder_name` | str | — | Yes | Where to save plots |
| `--eval_folder_name` | str | — | Yes | Path to Eval results |
| `--K` | int (nargs=\*) | `[30,50,60,80,100,200,250,300]` | No | K values |
| `--sel_threshs` | float (nargs=\*) | `[0.4,0.8,2.0]` | No | Density thresholds |
| `--groupby` | str | `sample` | No | Grouping variable |
| `--pval` | float | `0.05` | No | P-value threshold |
| `--samples` | str (nargs=\*) | None | No | Sample names |

---

## 5. Program Analysis Plot

**Script**: `Plotting/Slurm_Version/cNMF_program_analysis.py`
**Conda**: `NMF_Benchmarking`

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--mdata_path` | str | — | Yes | Path to .h5mu file |
| `--perturb_path_base` | str | — | Yes | Base path for perturbation results |
| `--GO_path` | str | — | Yes | Path to GO enrichment results |
| `--pdf_save_path` | str | — | Yes | Output directory for plots |
| `--file_to_dictionary` | str | None | No | Gene name mapping file |
| `--reference_gtf_path` | str | None | No | Reference GTF |
| `--tagert_col_name` | str | `program_name` | No | Target column in perturbation results |
| `--plot_col_name` | str | `target_name` | No | Gene column |
| `--log2fc_col` | str | `log2FC` | No | Log2FC column |
| `--top_program` | int | `10` | No | Top programs to display |
| `--top_enrichned_term` | int | `10` | No | Top GO terms per program |
| `--p_value` | float | `0.05` | No | Significance threshold |
| `--down_thred_log` | float | `-0.05` | No | Lower volcano threshold |
| `--up_thred_log` | float | `0.05` | No | Upper volcano threshold |
| `--figsize` | float (nargs=2) | `35 35` | No | Figure size |
| `--sample` | str (nargs=\*) | None | No | Sample names |
| `--square_plots` | flag | | No | Square aspect ratio |
| `--show` | flag | | No | Display interactively |
| `--PDF` | flag | | No | Save as PDF (else SVG) |

### Keys

`--data_key` (rna), `--prog_key` (cNMF), `--gene_name_key` (gene_names), `--categorical_key` (sample)

---

## 6. Perturbed Gene Analysis Plot

**Script**: `Plotting/Slurm_Version/cNMF_perturbed_gene_analysis.py`
**Conda**: `NMF_Benchmarking`

**Note**: This script does NOT have `--GO_path` or `--top_enrichned_term` (unlike Program Analysis).

| Parameter | Type | Default | Required | Description |
|-----------|------|---------|----------|-------------|
| `--mdata_path` | str | — | Yes | Path to .h5mu file |
| `--perturb_path_base` | str | — | Yes | Base path for perturbation results |
| `--pdf_save_path` | str | — | Yes | Output directory for plots |
| `--file_to_dictionary` | str | None | No | Gene name mapping file |
| `--reference_gtf_path` | str | None | No | Reference GTF |
| `--tagert_col_name` | str | `target_name` | No | Target column (reversed from program analysis) |
| `--plot_col_name` | str | `program_name` | No | Plot column (reversed from program analysis) |
| `--log2fc_col` | str | `log2FC` | No | Log2FC column |
| `--top_num` | int | `5` | No | Top genes per program |
| `--top_program` | int | `10` | No | Top programs to display |
| `--p_value` | float | `0.05` | No | Significance threshold |
| `--down_thred_log` | float | `-0.05` | No | Lower volcano threshold |
| `--up_thred_log` | float | `0.05` | No | Upper volcano threshold |
| `--figsize` | float (nargs=2) | `35 35` | No | Figure size |
| `--sample` | str (nargs=\*) | None | No | Sample names |
| `--n_processes` | int | `-1` | No | Parallel processes |
| `--dot_size` | int | `10` | No | UMAP dot size |
| `--square_plots` | flag | | No | Square aspect ratio |
| `--show` | flag | | No | Display interactively |
| `--PDF` | flag | | No | Save as PDF (else SVG) |
| `--expressed_only` | flag | | No | Only plot expressed perturbed genes |

### Keys

`--data_key` (rna), `--prog_key` (cNMF), `--gene_name_key` (gene_names), `--categorical_key` (sample)

---

## 7. U-test Calibration

**Script**: `Perturbation_association_calibration/Slurm_version/U-test_perturbation_calibration/U-test_perturbation_calibration.py`
**Conda**: `NMF_Benchmarking`

### Required

| Parameter | Type | Description |
|-----------|------|-------------|
| `--out_dir` | str | Directory containing cNMF output |
| `--run_name` | str | cNMF run name |
| `--mdata_guide_path` | str | Path to MuData with guide assignments |

### Optional

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--components` | int (nargs=\*) | `[30,50,60,80,100,200,250,300]` | K values |
| `--sel_thresh` | float (nargs=\*) | `[0.4,0.8,2.0]` | Density thresholds |
| `--guide_annotation_path` | str | None | Guide annotation TSV |
| `--guide_annotation_key` | str | `['non-targeting']` | Non-targeting guide label (default is a list) |
| `--reference_gtf_path` | str | None | Reference GTF |
| `--number_run` | int | `300` | Calibration iterations |
| `--number_guide` | int | `6` | Fake targeting guides per iteration |
| `--organism` | str | `human` | Species |
| `--FDR_method` | str | `StoreyQ` | FDR correction method |

### Workflow Flags

| Flag | Description |
|------|-------------|
| `--compute_real_perturbation_tests` | Run real perturbation tests |
| `--compute_fake_perturbation_tests` | Run calibration null distribution |
| `--visualizations` | Generate QQ and violin plots |
| `--check_format` | Validate format first |

### Keys

| Parameter | Default |
|-----------|---------|
| `--data_key` | `rna` |
| `--prog_key` | `cNMF` |
| `--categorical_key` | `sample` |
| `--guide_names_key` | `guide_names` |
| `--guide_targets_key` | `guide_targets` |
| `--guide_assignment_key` | `guide_assignment` |

---

## 8. CRT Calibration

**Script**: `Perturbation_association_calibration/Slurm_version/CRT/CRT.py`
**Conda**: `programDE`

### Required

| Parameter | Type | Description |
|-----------|------|-------------|
| `--out_dir` | str | Directory containing cNMF output |
| `--run_name` | str | cNMF run name |
| `--mdata_guide_path` | str | Path to MuData with guide assignments |

### Optional

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `--components` | int (nargs=\*) | `[30,50,60,80,100,200,250,300]` | K values |
| `--sel_thresh` | float (nargs=\*) | `[0.4,0.8,2.0]` | Density thresholds |
| `--categorical_key` | str | `sample` | Sample/condition key |
| `--covariates` | str (nargs=\*) | None | Covariate keys in obs (used as-is) |
| `--log_covariates` | str (nargs=\*) | None | Covariate keys to log1p-transform |
| `--number_guide` | int | `6` | Fake targeting guides per iteration |
| `--number_permutations` | int | `1024` | CRT permutations |
| `--guide_annotation_key` | str (nargs=\*) | `non-targeting` | Non-targeting label (accepts multiple values) |
| `--FDR_method` | str | `BH` | `BH` or `StoreyQ` (choices enforced) |
