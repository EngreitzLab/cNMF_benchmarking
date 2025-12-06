# cNMF Plotting Pipeline

Comprehensive visualization tools for cNMF (consensus Non-negative Matrix Factorization) analysis results.

## Overview

This module provides three main plotting pipelines for different aspects of cNMF analysis:

1. **K-Selection Analysis** - Optimal K value selection and model comparison
2. **Perturbed Gene Analysis** - Individual gene perturbation effects visualization  
3. **Program Analysis** - cNMF program characterization and annotation

## Components

### 1. K-Selection Plotting (`cNMF_k_selection.py`)

Generates plots to help select optimal K values and compare model performance.

#### Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| output_directory | str | - | Directory containing cNMF output files |
| run_name | str | - | Name of the cNMF run |
| save_folder_name | str | - | Output directory for plots |
| eval_folder_name | str | - | Directory containing evaluation results |
| groupby | str | "sample" | Column for grouping in analysis |
| K | list of int | [30, 50, 60, 80, 100, 200, 250, 300] | K values to analyze |
| sel_threshs | list of float | [0.4, 0.8, 2.0] | Selection thresholds |
| samples | list of str | ['D0', 'sample_D1', 'sample_D2', 'sample_D3'] | Sample names |
| pval | float | 0.05 | P-value threshold |

#### Outputs
- Stability & error plots across K values
- Program enrichment plots for each threshold
- Perturbation effect analysis
- Explained variance plots
- Program dotplots for detailed visualization

### 2. Perturbed Gene Analysis (`cNMF_perturbed_gene_analysis.py`)

Creates comprehensive analysis plots for each perturbed gene.

#### Parameters

**I/O Paths:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| mdata_path | str | - | Path to multimodal data file (.h5mu) |
| perturb_path_base | str | - | Base path for perturbation data files |
| pdf_save_path | str | - | Output path for PDF files |
| file_to_dictionary | str | None | Gene ID to name conversion file |
| reference_gtf_path | str | None | Reference GTF for validation |

**Visualization Options:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| tagert_col_name | str | "target_name" | Column name for targets |
| plot_col_name | str | "program_name" | Column name for plotting |
| log2fc_col | str | "log2FC" | Log2 fold change column |
| top_num | int | 5 | Number of top items to show |
| top_program | int | 10 | Number of top programs |
| p_value | float | 0.05 | P-value threshold |
| down_thred_log | float | -0.05 | Downregulation threshold |
| up_thred_log | float | 0.05 | Upregulation threshold |
| figsize | list | [35, 35] | Figure size |
| square_plots | flag | False | Create square plots |
| show | flag | False | Display plots |
| PDF | flag | False | Save as PDF |
| n_processes | int | -1 | Number of parallel processes |

#### Outputs
- Individual gene analysis PDFs with UMAP, correlation plots, volcano plots
- Gene correlation matrices
- Waterfall correlation plots
- Merged comprehensive PDF reports

### 3. Program Analysis (`cNMF_program_analysis.py`)

Generates detailed characterization plots for each cNMF program.

#### Parameters

**I/O Paths:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| mdata_path | str | - | Path to multimodal data file (.h5mu) |
| perturb_path_base | str | - | Base path for perturbation data files |
| GO_path | str | - | Path to Gene Ontology data |
| pdf_save_path | str | - | Output path for PDF files |

**Analysis Options:**
| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| tagert_col_name | str | "program_name" | Target column name |
| plot_col_name | str | "target_name" | Plot column name |
| top_program | int | 10 | Number of top programs |
| top_enrichned_term | int | 10 | Number of top enriched terms |
| figsize | list | [35, 35] | Figure size |

#### Outputs
- Program-specific UMAP visualizations
- Top genes per program analysis
- Gene Ontology enrichment plots
- Program correlation matrices
- Violin plots and heatmaps
- Comprehensive PDF reports

## Common Data Keys

All scripts use these configurable keys for data access:

| Key | Default | Description |
|-----|---------|-------------|
| data_key | "rna" | Access gene expression in mdata |
| prog_key | "cNMF" | Access cNMF programs in mdata |
| gene_name_key | "gene_names" | Gene names key |
| categorical_key | "sample" | Cell condition key |

 

## Output Organization

```
save_folder/
├── K_selection_plots/
│   ├── stability_plots.png
│   ├── enrichment_plots/
│   └── variance_plots/
├── gene_analysis/
│   ├── individual_genes/
│   └── merged_reports.pdf
└── program_analysis/
    ├── individual_programs/
    └── merged_reports.pdf
```