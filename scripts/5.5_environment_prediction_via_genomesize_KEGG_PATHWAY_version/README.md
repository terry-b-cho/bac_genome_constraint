# KEGG Environment Prediction Pipeline

## Overview

This pipeline adapts the GO-based environment prediction pipeline to work with KEGG functional data (reactions, KOs, pathways) for predicting 8-category environment labels using genome size and KEGG term counts.

## Directory Structure

```
scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/
├── 5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version.py  # Main script
├── run_prediction.sh                                                    # Bash wrapper
├── 00_Figure_descriptions_n_methods.md                                  # Documentation
└── README.md                                                            # This file

results/5.5_environment_prediction/
├── prev99_env_prediction_metrics_reactions.tsv
├── prev99_env_prediction_metrics_ko.tsv
├── prev99_env_prediction_metrics_pathway.tsv
├── prev99_fig01_model_performance_comparison_reactions.png
├── prev99_fig01_model_performance_comparison_ko.png
├── prev99_fig01_model_performance_comparison_pathway.png
└── ... (all outputs with appropriate suffixes)
```

## Quick Start

### Single Data Type

```bash
python3 scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version.py \
    --data-type reactions \
    --prevalence-threshold 99 \
    --model all \
    --plot
```

### All Data Types (Bash Wrapper)

```bash
bash scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/run_prediction.sh 99 "" all
```

### SLURM Job

```bash
sbatch scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/run_prediction.sh
```

## Arguments

### Required
- `--data-type`: KEGG data type (`reactions`, `ko`, `pathway`)

### Optional
- `--prevalence-threshold`: Prevalence threshold (0-100, e.g., 95 or 99). If not specified, uses all KEGG terms.
- `--test-mode`: Run on small subset (50 genomes per environment) for testing
- `--normalize`: Feature normalization (`none`, `per_gene`, `log`, `both`) - default: `none`
- `--model`: Which model(s) to train (`all`, `rf`, `gb`, `lr`, `xgb`, `baseline`) - default: `all`
- `--use-gpu`: Enable GPU acceleration for XGBoost (requires CUDA)
- `--plot`: Generate visualization plots (default: True)
- `--output-dir`: Custom output directory (default: `results/5.5_environment_prediction/`)

## Input Files

The pipeline expects the following input files (automatically located based on `--data-type`):

1. **Master Table:**
   - `results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts/master_table_env_filtered{suffix}.tsv`
   - Contains: Genome metadata, environment labels, KEGG term counts

2. **Valid Environments:**
   - `results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts/valid_environments_min20{suffix}.tsv`
   - Contains: List of 8 valid environment categories

3. **KEGG Labels:**
   - `results/4.5_statistical_analyses_KEGG_PATHWAY_version/05_kegg_labels/kegg_term_labels{suffix}.tsv`
   - Contains: Human-readable names for KEGG terms

Where `{suffix}` = `_reactions`, `_ko`, or `_pathway` based on `--data-type`.

## Output Files

All outputs are saved to `results/5.5_environment_prediction/` with format:
`{prefix}{description}{suffix}.{ext}`

### Metrics Files
- `prev99_env_prediction_metrics{suffix}.tsv` - Overall model performance metrics
- `prev99_env_prediction_per_class_metrics{suffix}.tsv` - Per-environment performance
- `prev99_env_prediction_confusion_matrix_{model}{suffix}.tsv` - Confusion matrices
- `prev99_env_prediction_feature_importances{suffix}.tsv` - Feature importance rankings

### Figures (PNG + PDF)
- `prev99_fig01_model_performance_comparison{suffix}.png/pdf`
- `prev99_fig02_confusion_matrices{suffix}.png/pdf`
- `prev99_fig03_roc_curves{suffix}.png/pdf`
- `prev99_fig04_per_class_metrics{suffix}.png/pdf`
- `prev99_fig05_feature_importance{suffix}.png/pdf`

### Models
- `prev99_env_prediction_model_{model}{suffix}.pkl` - Trained model files
- `prev99_env_prediction_scaler{suffix}.pkl` - Feature scaler (for LR)

### Logs
- `prev99_qc_05_env_prediction{suffix}.log` - QC log file

## Data Types

### Reactions
- **Format:** R-numbers (e.g., R00130)
- **Count:** ~485 reactions
- **Granularity:** Biochemical reaction level
- **Example:** "ATP:dephospho-CoA 3'-phosphotransferase"

### KOs
- **Format:** K-numbers (e.g., K00005)
- **Count:** ~739 KOs
- **Granularity:** Functional ortholog level
- **Example:** "glycerol dehydrogenase [EC:1.1.1.6]"

### Pathways
- **Format:** map-numbers (e.g., map00010)
- **Count:** ~101 pathways (only "map", no "rn" duplicates)
- **Granularity:** Pathway level
- **Example:** "Glycolysis / Gluconeogenesis"

## Testing Protocol

### Step 1: Test Mode (Single Data Type)

```bash
python3 scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version.py \
    --data-type reactions \
    --prevalence-threshold 99 \
    --test-mode \
    --model baseline \
    --plot
```

**Verify:**
- Script runs without errors
- Output files have `_reactions` suffix
- Log file is created
- Test mode uses ~50 genomes per environment

### Step 2: Test Mode (All Models)

```bash
python3 scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version.py \
    --data-type reactions \
    --prevalence-threshold 99 \
    --test-mode \
    --model all \
    --plot
```

**Verify:**
- All models train successfully
- All figures are generated
- Feature importance plots show KEGG names (not IDs)
- Confusion matrices are generated
- ROC curves are generated

### Step 3: Full Run (Single Data Type)

```bash
python3 scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version.py \
    --data-type reactions \
    --prevalence-threshold 99 \
    --model all \
    --plot
```

**Verify:**
- Full dataset processes successfully
- All outputs are generated
- Performance metrics are reasonable

### Step 4: Repeat for Other Data Types

Repeat Steps 1-3 for `ko` and `pathway` data types.

### Step 5: All Data Types (Bash Wrapper)

```bash
bash scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/run_prediction.sh 99 "" all
```

**Verify:**
- All three data types process successfully
- All outputs have correct suffixes
- No conflicts between data types

## Dependencies

### Python Packages
- pandas
- numpy
- scikit-learn
- xgboost (optional, for XGBoost model)
- matplotlib
- seaborn
- joblib or pickle (for model saving)

### O2 Modules
```bash
module load conda/miniforge3/24.11.3-0
module load gcc/14.2.0
module load cuda/12.8  # Optional, for GPU support
```

### Conda Environment
```bash
conda activate genome_constraint_envs_O2_py
```

## Key Differences from GO Pipeline

1. **Data Type Argument:** Required `--data-type` argument (reactions/ko/pathway)
2. **KEGG IDs:** Uses R-numbers, K-numbers, map-numbers (not 7-digit GO format)
3. **File Paths:** Uses KEGG-specific paths from `4.5_statistical_analyses_KEGG_PATHWAY_version`
4. **Feature Names:** `{data_type}_total` instead of `go_total`
5. **Label Files:** Uses `kegg_term_labels{suffix}.tsv` instead of GO labels
6. **Output Suffixes:** All outputs include `_reactions`, `_ko`, or `_pathway` suffix

## Troubleshooting

### Import Errors
- Ensure `prevalence_utils.py` is in `scripts/4.5_statistical_analyses_KEGG_PATHWAY_version/`
- Check that Python path includes the correct directory

### File Not Found Errors
- Verify input files exist in `results/4.5_statistical_analyses_KEGG_PATHWAY_version/`
- Check that suffixes match data type (e.g., `_reactions` for reactions)

### Memory Issues
- Use `--test-mode` for initial testing
- Reduce number of models with `--model rf` (single model)
- Consider reducing prevalence threshold

### GPU Issues
- XGBoost will fall back to CPU if GPU is not available
- Check CUDA installation if `--use-gpu` is specified

## Performance Notes

- **Training Time:** Varies by data type and model
  - Baseline: < 1 minute
  - Random Forest: ~5-10 minutes
  - Gradient Boosting: ~10-15 minutes
  - XGBoost: ~5-10 minutes (CPU) or ~2-5 minutes (GPU)
  - Logistic Regression: ~2-5 minutes

- **Memory Usage:** ~8-16 GB depending on data type and model

- **Expected Accuracy:** Varies by data type, typically 50-65% test accuracy

## Citation

If using this pipeline, please cite:
- The original GO-based pipeline (Script 05)
- KEGG database
- scikit-learn, XGBoost, and other ML libraries used

---

**Version:** 1.0  
**Last Updated:** December 2025  
**Script:** 5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version



