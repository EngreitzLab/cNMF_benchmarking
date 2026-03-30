---
name: cnmf-runner
description: Interactive cNMF pipeline runner. Guides users through configuring and submitting cNMF inference, evaluation, plotting, and calibration jobs on SLURM. Triggers on keywords like cNMF, NMF, inference, evaluation, calibration, SLURM, submit job, run pipeline, K selection, perturbation, gene programs.
user_invocable: true
---

# cNMF Pipeline Runner

You are an interactive assistant for running the cNMF benchmarking pipeline on Stanford Sherlock HPC. Guide the user step-by-step through data validation, parameter selection, resource estimation, SLURM script generation, and job submission.

## Constants

```
PIPELINE_ROOT=/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline
SKILL_DIR=/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/.claude/skills/cnmf-runner
CONDA_BASE=/oak/stanford/groups/engreitz/Users/ymo/miniforge3
EMAIL=ymo@stanford.edu
GWAS_DATA=/oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/Evaluation/Resources/OpenTargets_L2G_Filtered.csv.gz
REFERENCE_GTF=/oak/stanford/groups/engreitz/Users/opushkar/genome/IGVFFI9573KOZR.gtf.gz
```

## Workflow

When the user wants to run any pipeline stage, follow this conversational flow:

### Step 1: Identify the Stage

Ask the user which stage they want to run (or infer from context):

| Stage | generate_slurm.py `--stage` value | Conda Env | GPU |
|-------|-----------------------------------|-----------|-----|
| sk-cNMF Inference | `inference-sk` | `sk-cNMF` | No |
| torch-cNMF Inference | `inference-torch` | `torch-cNMF` | Yes (V100S) |
| Evaluation | `evaluation` | `NMF_Benchmarking` | No |
| K-selection Plot | `k-selection` | `torch-cNMF` | No |
| Program Analysis Plot | `program-analysis` | `NMF_Benchmarking` | No |
| Perturbed Gene Plot | `perturbed-gene` | `NMF_Benchmarking` | No |
| U-test Calibration | `u-test-calibration` | `NMF_Benchmarking` | No |
| CRT Calibration | `crt-calibration` | `programDE` | No |

### Step 2: Get Data Path & Validate

Ask for the input data file path. Then run validation:

```bash
eval "$(conda shell.bash hook)" && conda activate sk-cNMF && python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/.claude/skills/cnmf-runner/scripts/validate_data.py --counts_fn "<path>"
```

Add `--categorical_key`, `--guide_names_key`, etc. if user specifies non-default keys.

Report the validation results to the user. If issues are found, help resolve them before proceeding.

### Step 3: Collect Parameters

#### Inference (sk-cNMF or torch-cNMF)

Present parameters in two tiers. First collect the common ones, then explicitly offer the advanced ones.

**Tier 1 -- Common parameters** (always ask about these, show defaults):

| Parameter | Default | Description |
|-----------|---------|-------------|
| **Output directory** | (required) | Where results are saved |
| **Run name** | (required) | Convention: `MMDDYY_<description>` |
| **Species** | (required) | `human` or `mouse` |
| `--K` | `30 50 60 80 100 200 250 300` | K values (number of programs) |
| `--numiter` | `10` | Number of NMF replicates per K |
| `--numhvgenes` | `5451` | Highly variable genes to use |
| `--sel_thresh` | sk-cNMF: `2.0` / torch: `0.4 0.8 2.0` | Density thresholds |
| `--seed` | `14` | Random seed |

For **torch-cNMF**, also ask about these common torch-specific parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--algo` | `mu` | Algorithm: `halsvar` (recommended), `mu`, or `cd` |
| `--mode` | `batch` | NMF mode: `batch`, `online`, or `dataloader` |
| `--tol` | `1e-4` | Convergence tolerance |

For **sk-cNMF**, also ask about:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--algo` | `mu` | Algorithm: `mu` or `cd` |
| `--init` | `random` | Initialization: `random`, `nndsvd`, `nndsvda`, `nndsvdar` |
| `--tol` | `1e4` | Convergence tolerance (**WARNING**: likely a bug, recommend overriding to `1e-4`) |
| `--max_NMF_iter` | `500` | Maximum NMF iterations |

**Workflow flags** -- ask which steps to run:
- `--run_factorize` — run the NMF factorization (required for first run)
- `--run_refit` — run combine, k_selection_plot, and consensus (required for first run)
- `--run_complie_annotation` — compile results and generate gene annotations (required for first run)
- `--check_format` — validate input data format before running (recommended)
- `--parallel_running` — enable parallel processing for multiple K values

For a **first/complete run**, include all three: `--run_factorize --run_refit --run_complie_annotation`
For a **rerun after factorization** (e.g., changing sel_thresh), use only: `--run_refit --run_complie_annotation`

**Tier 2 -- Advanced parameters** (after Tier 1 is collected, present these in a table and ask if the user wants to change any):

For **both sk-cNMF and torch-cNMF**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--loss` | `frobenius` | NMF loss function |
| `--nmf_seeds_path` | None | Custom NMF seeds (.npy file) |
| `--num_gene` | `300` | Top genes for annotation |
| `--guide_annotation_path` | None | Guide annotation TSV |
| `--reference_gtf_path` | None | Reference GTF for gene validation |

For **torch-cNMF** only -- additional advanced parameters:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--n_jobs` | `-1` | Parallel jobs (-1 = all cores) |
| `--densify` | off | Densify sparse matrix before factorization |
| `--fp_precision` | `float` | `float` (32-bit) or `double` (64-bit) |
| `--alpha_usage` | `0.0` | Regularization for usage matrix (W) |
| `--alpha_spectra` | `0.0` | Regularization for spectra matrix (H) |
| `--l1_ratio_usage` | `0.0` | L1 vs L2 ratio for usage (0=L2, 1=L1) |
| `--l1_ratio_spectra` | `0.0` | L1 vs L2 ratio for spectra |
| `--shuffle_cells` | off | Shuffle cells before factorization |
| `--remove_noncoding` | off | Remove non-coding genes |
| `--sk_cd_refit` | off | Use sklearn coordinate descent for refitting |

If **torch-cNMF batch mode** is selected:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--batch_max_iter` | `500` | Max iterations for batch NMF |
| `--batch_hals_tol` | `0.05` | HALS tolerance |
| `--batch_hals_max_iter` | `200` | Max HALS inner iterations |

If **torch-cNMF online mode** is selected:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--online_max_pass` | `20` | Max passes through data |
| `--online_chunk_size` | `5000` | Chunk size for online learning |
| `--online_chunk_max_iter` | `200` | Max iterations per chunk |
| `--online_usage_tol` | `0.05` | Usage update tolerance |
| `--online_spectra_tol` | `0.05` | Spectra update tolerance |

Present the Tier 2 table and ask: "Do you want to change any of these advanced parameters, or use the defaults?"

#### Evaluation

Required parameters:
- `--out_dir` — directory containing cNMF output
- `--run_name` — cNMF run name
- `--gwas_data_path` — path to GWAS data (use GWAS_DATA constant above)

Test flags (ask which to enable):
- `--Perform_categorical` — Kruskal-Wallis + Dunn's test
- `--Perform_perturbation` — perturbation association
- `--Perform_geneset` — gene set enrichment (Reactome + GO)
- `--Perform_trait` — GWAS trait enrichment
- `--Perform_explained_variance` — explained variance per K (needs `--X_normalized_path`)
- `--Perform_motif` — TF motif enrichment (WIP)

#### Plotting and Calibration

For these stages, collect the appropriate parameters (refer to parameter catalog). Key required parameters:
- **K-selection**: `--output_directory`, `--run_name`, `--save_folder_name`, `--eval_folder_name`
- **Program analysis**: `--mdata_path`, `--perturb_path_base`, `--GO_path`, `--pdf_save_path`
- **Perturbed gene analysis**: `--mdata_path`, `--perturb_path_base`, `--pdf_save_path`
- **U-test calibration**: `--out_dir`, `--run_name`, `--mdata_guide_path`
- **CRT calibration**: `--out_dir`, `--run_name`, `--mdata_guide_path`

### Step 4: Estimate Resources

Use the data statistics from validation (cell count, gene count, file size) plus chosen parameters to estimate SLURM resources.

#### Memory Estimation

**Inference (sk-cNMF / CPU):**
- Base: `cells x numhvgenes x 8 bytes x 3` (matrix + W + H), plus 50% overhead
- Tiers:
  - <50k cells: 64G
  - 50k-200k cells: 128-256G
  - 200k-500k cells: 256-512G
  - >500k cells: 512-786G

**Inference (torch-cNMF / GPU):**
- GPU bottleneck: V100S = 32GB VRAM
- System RAM: `max(64G, cells x numhvgenes x 8 bytes x 2)`
- Tiers:
  - <100k cells: 64G
  - 100k-300k cells: 128G
  - >300k cells: 256G (consider sk-cNMF instead)

**Evaluation / Plotting / Calibration:**
- `max(64G, 2 x file_size_on_disk)`
- Calibration with many K x sel_thresh combinations: multiply by num_combinations

#### Time Estimation

**Inference (sk-cNMF):** scales with `num_K x numiter x cells x numhvgenes`
- <50k cells: 4-6h | 50k-200k: 8-14h | >200k: 14-24h
- (for numiter=10, 8 K values; double for numiter=20)

**Inference (torch-cNMF):** ~5-10x faster than CPU
- <100k cells: 2-4h | 100k-300k: 6-10h | >300k: 10-15h

**Evaluation:** 1-5h | **Plotting:** 1-3h | **Calibration:** 1-2h per K x sel_thresh

#### CPUs
- sk-cNMF inference: 1-10
- torch-cNMF inference: 1 (GPU-bound)
- Evaluation/Plotting/Calibration: 10-20

#### Partitions
- Default: `engreitz,owners`
- If memory >256G: `engreitz,owners,bigmem`
- For torch-cNMF (GPU): `engreitz` only

Present the recommendation and ask user to confirm or adjust.

### Step 5: Generate SLURM Script

Run the generator. Use `--` to separate generator args from passthrough pipeline args:

```bash
python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/.claude/skills/cnmf-runner/scripts/generate_slurm.py \
  --stage <stage_value> \
  --job_name <run_name> \
  --output_dir <out_dir> \
  --run_name <run_name> \
  --cpus <N> --mem <MG> --time <HH:MM:SS> --partition <parts> \
  [--gpu] \
  --script_output_path <out_dir>/<run_name>/<run_name>.sh \
  -- \
  [all stage-specific args for the pipeline script...]
```

**Stage value mapping**: Use the `--stage` value from the table in Step 1 (e.g., `inference-sk`, `inference-torch`, `evaluation`, etc.)

Show the generated script to the user for review. Read the generated .sh file and display it.

### Step 6: Submit

After user confirms, submit the job:

```bash
sbatch <script_path>
```

Report the job ID and explain monitoring:
- `squeue -u ymo` to check job status
- `sacct -j <job_id> --format=JobID,JobName,State,Elapsed,MaxRSS,MaxVMSize` for completed job details
- Log file locations:
  - Inference/Plotting: `<out_dir>/<run_name>/logs/<job_id>.out` and `.err`
  - Evaluation/Calibration: `<out_dir>/<run_name>/Eval/logs/<job_id>.out` and `.err`

### Step 7: Generate h5mu Structure Files

After inference completes (or when the user has existing .h5mu output), generate structure summary files for all output MuData files. These provide a quick human-readable overview of the h5mu contents (modalities, X matrices, obs/var columns, uns, obsm, layers, varm, etc.) in a tree format.

**For all .h5mu files in the adata directory (recommended):**

```bash
eval "$(conda shell.bash hook)" && conda activate sk-cNMF && python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/.claude/skills/cnmf-runner/scripts/generate_h5mu_structure.py dummy --adata_dir <out_dir>/<run_name>/adata
```

**For a single .h5mu file:**

```bash
eval "$(conda shell.bash hook)" && conda activate sk-cNMF && python3 /oak/stanford/groups/engreitz/Users/ymo/Tools/cNMF_benchmarking/cNMF_benchmarking_pipeline/.claude/skills/cnmf-runner/scripts/generate_h5mu_structure.py <out_dir>/<run_name>/adata/cNMF_<K>_<sel_thresh>.h5mu
```

This produces `_structure.txt` files alongside each `.h5mu` file (e.g., `cNMF_50_0_5_structure.txt` next to `cNMF_50_0_5.h5mu`).

**When to run this step:**
- After inference completes and .h5mu files are generated in `<out_dir>/<run_name>/adata/`
- When the user asks to inspect or summarize a .h5mu file
- Before evaluation or plotting, so the user can verify the output structure

Always offer to run this step after inference completes or when the user has .h5mu output files.

### Step 8: Post-Submission Guidance

After inference completes, guide the user through the next stages:
1. **Generate h5mu structure files** — summarize output MuData contents (Step 7 above)
2. **Evaluation** — run statistical tests on programs
3. **K-selection plot** — visualize stability and error across K values
4. **Program analysis** — per-program detailed plots
5. **Perturbed gene analysis** — per-gene detailed plots
6. **Calibration** — U-test or CRT calibration

## Important Notes

- Always use `eval "$(conda shell.bash hook)"` before conda activate in any Bash command
- For torch-cNMF, the generator auto-adds `--gres=gpu:1` and `-C GPU_SKU:V100S_PCIE`
- The pipeline saves config YAMLs automatically; no need to create them manually
- Run name convention: `MMDDYY_<short_description>`
- Output MuData files are at `<out_dir>/<run_name>/adata/cNMF_<K>_<sel_thresh>.h5mu` with corresponding `_structure.txt` summaries
- Generated SLURM scripts for inference stages automatically run h5mu structure generation after the pipeline completes. Pass `--no_structure` to `generate_slurm.py` to skip this
- The generated script creates log directories automatically via `mkdir -p`
- Known typos in pipeline argparse (use as-is): `--run_complie_annotation` (not compile), `--tagert_col_name` (not target), `--top_enrichned_term` (not enriched)
- sk-cNMF `--tol` default is `1e4` (likely a bug; the torch-cNMF default `1e-4` is correct). Consider passing `--tol 1e-4` for sk-cNMF runs.
- sk-cNMF `--guide_assignment_key` default is `guide_assignment_key` (the string), while torch-cNMF default is `guide_assignment`. If switching between the two, explicitly pass this parameter.
