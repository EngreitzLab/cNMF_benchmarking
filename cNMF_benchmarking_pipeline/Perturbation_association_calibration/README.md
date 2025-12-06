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

