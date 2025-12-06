# Perturbation Association Calibration

Statistical validation tool for perturbation association tests in cNMF analysis.

## Purpose

This module validates the statistical properties of perturbation association tests by generating null distributions and comparing them against real perturbation results. It ensures that p-values are properly calibrated and false positive rates are controlled.

## Methodology

The calibration performs two main analyses:

1. **Real Perturbation Test**: Analyzes actual perturbation effects on cNMF programs
2. **Null/Fake Perturbation Test**: Generates null distribution by randomly reassigning non-targeting guides as "targeting"

The results are compared using QQ plots to assess statistical calibration.

## Parameters

### Required Parameters

| Parameter | Type | Description |
|-----------|------|-------------|
| out_dir | str | Path to cNMF object directory |
| run_name | str | Name of cNMF object |
| save_path | str | Directory to save output figures |
| X_normalized_path | str | Path to normalized input cell x gene matrix |
| mdata_guide_path | str | Path to mdata with guide assignments |
| gwas_data_path | str | Path to GWAS information file |

### Optional Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| guide_annotation_path | str | None | Path to guide annotation file |
| guide_annotation_key | str | "non-targeting" | Key for guide annotation |
| reference_gtf_path | str | None | Path to reference GTF file |
| data_key | str | "rna" | Key to access gene expression in mdata |
| prog_key | str | "cNMF" | Key to access cNMF programs in mdata |
| categorical_key | str | "sample" | Key to access cell conditions in obs |
| organism | str | "human" | Data species |
| FDR_method | str | "StoreyQ" | FDR correction method |
| num_runs | list of int | [10] | Number of calibration iterations |
| num_guides | list of int | [5] | Number of guides for calibration |
| components | list of int | None | K values to test |
| sel_threshs | list of float | None | Selection thresholds |

### Flags

| Flag | Description |
|------|-------------|
| --compute_real_pertubration_test | Compute real perturbation associations |
| --compute_fake_pertubration_test | Compute null/fake perturbation associations |

## Usage

### Basic Usage
```bash
python Simple_perturbation_calibration.py \
  --out_dir /path/to/cnmf/results \
  --run_name my_analysis \
  --save_path /path/to/figures \
  --X_normalized_path /path/to/normalized_matrix.h5ad \
  --mdata_guide_path /path/to/guide_data.h5mu \
  --gwas_data_path /path/to/gwas_data.csv \
  --compute_real_pertubration_test \
  --compute_fake_pertubration_test
```

### With Custom Parameters
```bash
python Simple_perturbation_calibration.py \
  --out_dir /path/to/cnmf/results \
  --run_name my_analysis \
  --save_path /path/to/figures \
  --X_normalized_path /path/to/normalized_matrix.h5ad \
  --mdata_guide_path /path/to/guide_data.h5mu \
  --gwas_data_path /path/to/gwas_data.csv \
  --components 50 100 200 \
  --sel_threshs 0.4 0.8 2.0 \
  --num_runs 20 \
  --compute_real_pertubration_test \
  --compute_fake_pertubration_test
```

## Output Files

1. **Perturbation Results**: `perturbation_association_*_k{K}_sample{sample}.txt` - Statistical test results
2. **QQ Plots**: 
   - `perturbation_association_qqplot_real.png` - QQ plot of real p-values
   - `perturbation_association_qqplot_null.png` - QQ plot of null p-values
3. **Calibration Summary**: `perturbation_association_calibration.txt` - Combined calibration results

## Interpretation

- **Well-calibrated test**: Real p-values should deviate from diagonal in QQ plot, while null p-values should align with diagonal
- **Inflation**: If null p-values deviate below diagonal, test statistics may be inflated
- **Under-powered**: If real p-values align with diagonal, perturbation effects may be weak or undetectable

## Dependencies

- Standard: os, sys, yaml, logging, argparse
- Scientific: numpy, pandas, scanpy, muon
- Custom: Pipeline-specific modules for cNMF analysis

## Version Requirements

Compatible with the cNMF benchmarking pipeline v1.0+