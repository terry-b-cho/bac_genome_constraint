# Code Workflow and Figure Production Summary

## Overview

This document summarizes the complete workflow of the statistical analyses pipeline and details how each plot/figure is produced.

## Execution Order

The scripts must be run in this order:

1. **01_build_master_table_env.py** - Build master table with QC filters
2. **02_define_env_cohorts.py** - Select environments with ≥20 genomes
3. **03_fit_global_scaling.py** - Fit global scaling laws
4. **04_fit_env_scaling_and_Z.py** - Fit per-environment scaling and compute Z-scores
5. **05_map_go_labels.py** - Map GO IDs to names (can run in parallel with 03-04)
6. **04.5_plot_intermediate_figures.py** - Generate intermediate figures (requires 01-04)
7. **06_make_scaling_figures.py** - Generate publication figures (requires 03, 04, 05)
8. **07_scaling_extensions_tf_mobile_nutrient.py** - Extensions for TF and mobile elements (optional)

## Script-by-Script Breakdown

### Script 01: Build Master Table (`01_build_master_table_env.py`)

**Purpose**: Join GO counts with genome metadata and apply QC filters

**Inputs**:
- `results/3_GO_analyses/ubiquitous_counts_table.txt` (GO term counts)
- `results/3_GO_analyses/2_JGIgold_KEGG_anayses_out/05_genome_feature_matrix.tsv` (genome metadata)

**Outputs**:
- `results/4_statistical_analyses/01_master_table/master_table_raw.tsv` (before QC)
- `results/4_statistical_analyses/01_master_table/master_table_high_quality.tsv` (after QC)
- `results/4_statistical_analyses/01_master_table/environment_counts_all.tsv` (environment summary)
- `results/4_statistical_analyses/01_master_table/qc_01_master_table.log` (QC log)

**QC Filters Applied**:
1. CheckM completeness > 90% (NaN treated as 0)
2. CheckM contamination < 5% (NaN treated as 100)
3. genes_total > 0
4. Environment not null/empty/"Unclassified"

**Key Operations**:
- Inner join on `accession` = `Genome`
- Validates data integrity (no duplicates, numeric checks)
- Computes environment-level statistics

**No figures produced** - Data preparation only

---

### Script 02: Define Environment Cohorts (`02_define_env_cohorts.py`)

**Purpose**: Select environments with sufficient sample size (≥20 genomes)

**Inputs**:
- `results/4_statistical_analyses/01_master_table/master_table_high_quality.tsv`

**Outputs**:
- `results/4_statistical_analyses/02_env_cohorts/valid_environments_min20.tsv` (list of valid environments)
- `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv` (filtered master table)
- `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet` (Parquet format for speed)
- `results/4_statistical_analyses/02_env_cohorts/qc_02_env_cohorts.log` (QC log)

**Key Operations**:
- Filters environments with < 20 genomes
- Optional prevalence filtering for GO terms (via `--prevalence-threshold` argument)
- Creates analysis-ready master table

**No figures produced** - Data filtering only

---

### Script 03: Fit Global Scaling (`03_fit_global_scaling.py`)

**Purpose**: Fit power-law scaling across all environments (global fit)

**Inputs**:
- `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet` (or .tsv)

**Model**:
- Power law: `nc(g) = β × n(g)^α`
- Log-log form: `log(nc(g)) = log(β) + α × log(n(g))`
- OLS regression on log-transformed data

**Outputs**:
- `results/4_statistical_analyses/03_global_scaling/global_scaling_params.tsv`
  - Columns: `category`, `alpha_global`, `alpha_global_se`, `alpha_global_ci99_low`, `alpha_global_ci99_high`, `beta_global_log`, `beta_global_log_se`, `beta_global_log_ci99_low`, `beta_global_log_ci99_high`, `n_genomes_used`, `r_squared`, `p_value`
- `results/4_statistical_analyses/03_global_scaling/global_scaling_params.parquet` (optional)
- `results/4_statistical_analyses/03_global_scaling/qc_03_global_scaling.log` (QC log)

**Key Operations**:
- For each GO category:
  - Filter genomes with `genes_total > 0` and `category_count > 0`
  - Compute `x = log(genes_total)` and `y = log(category_count)`
  - Fit OLS regression: `y = β + α × x`
  - Calculate standard errors and 99% confidence intervals
  - Skip categories with < 10 genomes or insufficient variance (std < 0.05)

**No figures produced** - Parameter estimation only

---

### Script 04: Fit Environment-Specific Scaling and Z-Scores (`04_fit_env_scaling_and_Z.py`)

**Purpose**: Fit per-environment scaling laws and compute Z-scores comparing to global

**Inputs**:
- `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet`
- `results/4_statistical_analyses/03_global_scaling/global_scaling_params.tsv`
- `results/4_statistical_analyses/02_env_cohorts/valid_environments_min20.tsv`

**Model**:
- Same power-law model as Script 03, but fitted separately for each environment
- Z-scores: `Z = (α_env - α_global) / sqrt(SE_env² + SE_global²)`

**Outputs**:
- `results/4_statistical_analyses/04_env_scaling/env_scaling_params.tsv`
  - Columns: `environment`, `category`, `alpha_env`, `alpha_env_se`, `alpha_env_ci99_low`, `alpha_env_ci99_high`, `beta_env_log`, `beta_env_log_se`, `n_genomes_used`, `r_squared`, `p_value`
- `results/4_statistical_analyses/04_env_scaling/env_vs_global_Z_scores.tsv`
  - Columns: `environment`, `category`, `Z_alpha`, `Z_beta`, `n_genomes_used`
- `results/4_statistical_analyses/04_env_scaling/category_Z_summary.tsv`
  - Columns: `category`, `Z_alpha_category`, `Z_beta_category`, `n_envs_used`
  - Aggregated Z-scores: `Z_category = sqrt(mean(Z²))` across environments
- `results/4_statistical_analyses/04_env_scaling/qc_04_env_scaling.log` (QC log)

**Key Operations**:
- For each environment × category combination:
  - Filter genomes with `genes_total > 0` and `category_count > 0`
  - Require ≥ 10 genomes (skip if fewer)
  - Fit OLS regression
  - Compute Z-scores comparing to global parameters
- Aggregate Z-scores across environments for each category

**No figures produced** - Parameter estimation only

---

### Script 05: Map GO Labels (`05_map_go_labels.py`)

**Purpose**: Map GO term IDs to human-readable names

**Inputs**:
- `results/3_GO_analyses/ubiquitous_terms.txt` (list of GO term IDs)
- `data/go/go-basic.obo` (GO ontology file)

**Outputs**:
- `results/4_statistical_analyses/05_go_labels/go_term_labels.tsv`
  - Columns: `category`, `go_id`, `name`, `namespace`, `definition`
- `results/4_statistical_analyses/05_go_labels/go_term_labels_for_plots.tsv`
  - Columns: `category`, `short_label`, `full_label`, `super_category`
- `results/4_statistical_analyses/05_go_labels/qc_05_go_labels.log` (QC log)

**Key Operations**:
- Parse OBO file to extract term names and definitions
- Map 7-digit GO IDs (e.g., "0000015") to full GO IDs (e.g., "GO:0000015")
- Create plot-friendly labels

**No figures produced** - Label mapping only

---

### Script 04.5: Plot Intermediate Figures (`04.5_plot_intermediate_figures.py`)

**Purpose**: Generate intermediate figures for QC and validation

**Inputs**: Outputs from scripts 01-04

**Figures Produced**:

#### Script 01 Figures (`script_01/` directory):

1. **Fig 01A: QC Filtering Flowchart** (`fig_01A_QC_filtering_flowchart.png/pdf`)
   - Flowchart showing genome counts at each QC step
   - Created using matplotlib rectangles and arrows
   - Shows: Start → Completeness filter → Contamination filter → Genes filter → Environment filter → Final

2. **Fig 01B: Completeness Distribution** (`fig_01B_completeness_distribution.png/pdf`)
   - Histogram of CheckM completeness scores
   - Shows raw vs. high-quality genomes
   - Vertical line at 90% cutoff

3. **Fig 01C: Contamination Distribution** (`fig_01C_contamination_distribution.png/pdf`)
   - Histogram of CheckM contamination scores
   - Shows raw vs. high-quality genomes
   - Vertical line at 5% cutoff

4. **Fig 01D: Genome Size Distribution** (`fig_01D_genome_size_distribution.png/pdf`)
   - Histogram of `genes_total` (genome size)
   - Shows before vs. after QC

5. **Fig 02A: Environment Counts (QC)** (`fig_02A_environment_counts_QC.png/pdf`)
   - Bar plot comparing environment counts before and after QC
   - Side-by-side bars for each environment

6. **Fig 02D: Genome Size by Environment** (`fig_02D_genome_size_by_environment.png/pdf`)
   - Box plot of `genes_total` stratified by environment
   - Shows environment-specific genome size distributions

#### Script 02 Figures (`script_02/` directory):

1. **Fig 02B: Final Environments** (`fig_02B_final_environments.png/pdf`)
   - Bar plot of environments with ≥20 genomes
   - Shows final retained environments

2. **Fig 02C: Environment Contribution** (`fig_02C_environment_contribution.png/pdf`)
   - Two-panel figure: bar plot and pie chart
   - Shows proportion of genomes from each environment

3. **Fig 02D: Genome Size (Retained Environments)** (`fig_02D_genome_size_retained_envs.png/pdf`)
   - Box plot of genome size for retained environments only

#### Script 03 Figures (`script_03/` directory):

1. **Fig 04A: Global Exponent Histogram** (`fig_04A_global_exponents_histogram.png/pdf`)
   - Histogram of `alpha_global` values across all categories
   - Vertical lines: median and α=1 (linear scaling)

2. **Fig 04B: Exponent vs R²** (`fig_04B_alpha_vs_r2.png/pdf`)
   - Scatter plot: `alpha_global` vs `r_squared`
   - Identifies categories with poor fit quality

3. **Fig 04C: Exponent vs Mean Count** (`fig_04C_alpha_vs_mean_count.png/pdf`)
   - Scatter plot: `alpha_global` vs mean GO category count
   - Log-log scale

4. **Fig 04D: Representative Scaling Plots** (`fig_04D_representative_scaling.png/pdf` and `fig_04D_representative_scaling_linear.png/pdf`)
   - Grid of 20 subplots (4×5) showing scaling for top 20 categories
   - Selected by highest `Z_alpha_category` (or `alpha_global` if Z-scores unavailable)
   - Each panel: scatter plot with fitted line
   - Log scale and linear scale versions
   - Metabolic version available if metabolic terms file exists

#### Script 04 Figures (`script_04/` directory):

1. **Fig 05A: Z-Score Heatmap (Exponents)** (`fig_05A_Z_alpha_heatmap.png/pdf`)
   - Heatmap of `Z_alpha` values (categories × environments)
   - Hierarchical clustering (Euclidean distance, Ward linkage)
   - Top 50 categories by `Z_alpha_category`
   - Color scale: RdBu_r (blue = negative, red = positive)

2. **Fig 05B: Z-Score Heatmap (Offsets)** (`fig_05B_Z_beta_heatmap.png/pdf`)
   - Same as 05A but for `Z_beta` values

3. **Fig 05C: Absolute Z-Score by Environment** (`fig_05C_abs_Z_alpha_by_env.png/pdf`)
   - Box plot of `|Z_alpha|` values stratified by environment
   - Horizontal line at |Z| = 2

4. **Fig 05D: Significant Categories by Environment** (`fig_05D_significant_categories_by_env.png/pdf`)
   - Bar plot: number of categories with `|Z_alpha| > 2` per environment

5. **Fig 06A: Z-Score Distribution (Exponents)** (`fig_06A_Z_alpha_category_histogram.png/pdf`)
   - Histogram of `Z_alpha_category` values
   - Vertical line at Z = 2

6. **Fig 06B: Z-Score Distribution (Offsets)** (`fig_06B_Z_beta_category_histogram.png/pdf`)
   - Histogram of `Z_beta_category` values

7. **Fig 06C: Exception Categories** (`fig_06C_exception_categories.png/pdf`)
   - Horizontal bar plot of top 20 categories with `Z_alpha_category > 2`
   - Labels show GO term names

8. **Fig 07: Environment-Stratified Scaling** (`fig_07_env_stratified_scaling.png/pdf` and `fig_07_env_stratified_scaling_linear.png/pdf`)
   - Grid of 20 subplots (4×5) for top 20 variable categories
   - Each panel: scatter plot with points colored by environment
   - Environment-specific fit lines (colored) and global fit (gray)
   - Log scale and linear scale versions
   - Metabolic version available

**Command-line Arguments**:
- `--script`: Which script figures to generate ('01', '02', '03', '04', 'all')
- `--test-mode`: Generate subset of figures for testing
- `--prevalence-threshold`: Filter by prevalence threshold (e.g., 95)

---

### Script 06: Make Scaling Figures (`06_make_scaling_figures.py`)

**Purpose**: Generate publication-quality figures (main results)

**Inputs**:
- `results/4_statistical_analyses/04_env_scaling/category_Z_summary.tsv`
- `results/4_statistical_analyses/04_env_scaling/env_scaling_params.tsv`
- `results/4_statistical_analyses/03_global_scaling/global_scaling_params.tsv`
- `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv`
- `results/4_statistical_analyses/05_go_labels/go_term_labels_for_plots.tsv`

**Figures Produced** (`results/4_statistical_analyses/06_figures/`):

1. **Panel 1a: Z-Statistics for Exponents by Category** (`fig1a_Z_exponents_by_category_env.png/pdf`)
   - Bar plot of `Z_alpha_category` for all categories, sorted by magnitude
   - Horizontal line at Z = 2
   - Top categories labeled with GO term names
   - Metabolic version available

2. **Panel 1b: Z-Statistics for Offsets by Category** (`fig1b_Z_offsets_by_category_env.png/pdf`)
   - Same as 1a but for `Z_beta_category`

3. **Panels 1c-1e: Exponent Comparisons** (`fig1cde_env_exponents_selected_categories.png/pdf`)
   - Grid of subplots (up to 10 categories)
   - Each panel: environment-specific exponents with 99% CI error bars
   - Horizontal dashed line: global exponent
   - Categories selected: low Z, medium Z, high Z, quartiles

4. **Panels 1f-1k: Scatter Plots with Fits** (`fig1f_to_k_env_scatter_scaling.png/pdf` and `fig1f_to_k_env_scatter_scaling_linear.png/pdf`)
   - Grid of subplots (5 categories × 3 environments = 15 panels)
   - Each panel:
     - Gray dots: all genomes (all environments)
     - Colored dots: genomes from specific environment
     - Colored line: environment-specific fit
     - Gray line: global fit
   - Log scale and linear scale versions
   - Panel metadata saved to `fig1_panel_metadata.tsv`

5. **Panels: GO Category vs Total Annotated Domains** (`fig1_domains_vs_total_domains.png/pdf` and `fig1_domains_vs_total_domains_linear.png/pdf`)
   - Similar to panels 1f-1k but X-axis is total annotated domains (sum of all GO counts)
   - Tests whether scaling is driven by genome size or functional diversity

**Key Operations**:
- Selects categories based on Z-score variance
- Selects top 3 environments by genome count for each category
- Creates both log-scale and linear-scale versions of scatter plots
- Generates metabolic-focused versions if metabolic terms file exists

---

### Script 07: Scaling Extensions (`07_scaling_extensions_tf_mobile_nutrient.py`)

**Purpose**: Extend scaling analysis to transcription factors and mobile elements

**Inputs**:
- `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet`

**Model**:
- Same power-law model applied to `tf_count` and `mobile_element_count`

**Outputs**:
- `results/4_statistical_analyses/07_extensions/tf_mobile_scaling_params.tsv`
  - Columns: `category` ("tf_count" or "mobile_element_count"), `environment` ("global" or env name), `alpha`, `alpha_se`, `beta_log`, `beta_log_se`, `n_genomes_used`, `r_squared`, `p_value`
- `results/4_statistical_analyses/07_extensions/tf_mobile_env_Z_scores.tsv`
  - Columns: `category`, `environment`, `Z_alpha`, `Z_beta`, `n_genomes_used`
- `results/4_statistical_analyses/07_extensions/qc_07_extensions.log` (QC log)

**Key Operations**:
- Global fits for `tf_count` and `mobile_element_count`
- Per-environment fits
- Z-scores comparing environment-specific to global
- Outlier detection using Cook's distance

**Expected Results**:
- TF scaling: α ≈ 2 (quadratic scaling, regulatory complexity hypothesis)
- Mobile elements: Variable scaling

**No figures produced** - Parameter estimation only (figures could be added)

---

## Utility Functions

### `prevalence_utils.py`

**Functions**:
- `get_prevalence_prefix(prevalence_threshold)`: Returns filename prefix (e.g., "prev95_")
- `filter_go_columns_by_prevalence(df, prevalence_threshold)`: Filters GO columns by prevalence threshold
- `load_prevalence_filtered_terms(prevalence_threshold, base_dir)`: Loads pre-computed filtered terms

**Usage**: Used by scripts 02-07 to support optional prevalence filtering

### `create_prevalence_filtered_terms.py`

**Purpose**: Create prevalence-filtered GO terms file

**Command**: `python create_prevalence_filtered_terms.py --prevalence-threshold 95`

**Output**: `results/3_GO_analyses/prev95_ubiquitous_terms.txt`

### `extract_metabolic_go_terms.py`

**Purpose**: Extract metabolism-related GO terms from ontology

**Method**:
- Parses GO ontology (go-basic.obo)
- Traverses DAG starting from metabolism root terms (GO:0008152, etc.)
- Collects all descendant terms
- Filters to terms present in dataset

**Output**: `results/4_statistical_analyses/05_go_labels/metabolic_go_terms.txt`

**Usage**: Used by figure scripts to generate metabolism-focused versions

---

## Command-Line Arguments

All scripts support:
- `--test-mode`: Run on small subset for testing
- `--prevalence-threshold`: Filter GO terms by prevalence (0-100, e.g., 95 for 95%)

Script 04.5 also supports:
- `--script`: Which script figures to generate ('01', '02', '03', '04', 'all')

---

## File Naming Conventions

**Base files**: `{description}.{ext}`

**With prevalence filter**: `prev{P}_{description}.{ext}` (e.g., `prev95_global_scaling_params.tsv`)

**Metabolic focus**: `metabolic_{description}.{ext}` (e.g., `metabolic_fig1a_Z_exponents_by_category_env.png`)

**Linear scale**: `{description}_linear.{ext}` (e.g., `fig1f_to_k_env_scatter_scaling_linear.png`)

**Combined**: `prev{P}_metabolic_{description}_linear.{ext}`

---

## Key Data Structures

### Master Table Structure
- **Primary key**: `accession` (GCF_* format)
- **Environment**: `environment` (from GOLD metadata)
- **Genome size**: `genes_total` (n(g) for scaling laws)
- **GO counts**: 334 columns with 7-digit GO IDs (nc(g) for each category)
- **Quality metrics**: `checkm_completeness`, `checkm_contamination`
- **Additional**: `tf_count`, `mobile_element_count`, etc.

### Scaling Parameters Structure
- **Global**: One row per GO category
- **Environment-specific**: One row per environment × category combination
- **Z-scores**: One row per environment × category combination
- **Category summary**: One row per category (aggregated Z-scores)

---

## Mathematical Framework

### Power-Law Model
```
nc(g) = β × n(g)^α
```

### Log-Log Form (for regression)
```
log(nc(g)) = log(β) + α × log(n(g))
```

### Z-Score Formula
```
Z = (α_env - α_global) / sqrt(SE_env² + SE_global²)
```

### Category-Level Z-Score (aggregated)
```
Z_category = sqrt(mean(Z² across environments))
```

---

## Quality Control Checks

### Minimum Sample Sizes
- **Per-environment fit**: ≥ 10 genomes
- **Environment selection**: ≥ 20 genomes
- **Variance check**: std(log(n(g))) ≥ 0.05

### Data Validation
- No negative counts
- No duplicate accessions
- Valid environment labels (not null/empty/"Unclassified")
- CheckM completeness > 90%
- CheckM contamination < 5%

### Fit Quality
- R² reported (no hard threshold, but low R² < 0.3 flagged)
- Standard errors calculated
- 99% confidence intervals computed
- P-values from OLS regression

---

## Output Directory Structure

```
results/4_statistical_analyses/
├── 01_master_table/
│   ├── master_table_raw.tsv
│   ├── master_table_high_quality.tsv
│   ├── environment_counts_all.tsv
│   └── qc_01_master_table.log
├── 02_env_cohorts/
│   ├── valid_environments_min20.tsv
│   ├── master_table_env_filtered.tsv
│   ├── master_table_env_filtered.parquet
│   └── qc_02_env_cohorts.log
├── 03_global_scaling/
│   ├── global_scaling_params.tsv
│   ├── global_scaling_params.parquet
│   └── qc_03_global_scaling.log
├── 04_env_scaling/
│   ├── env_scaling_params.tsv
│   ├── env_vs_global_Z_scores.tsv
│   ├── category_Z_summary.tsv
│   └── qc_04_env_scaling.log
├── 04.5_intermediate_figures/
│   ├── script_01/
│   ├── script_02/
│   ├── script_03/
│   ├── script_04/
│   └── qc_04.5_figures.log
├── 05_go_labels/
│   ├── go_term_labels.tsv
│   ├── go_term_labels_for_plots.tsv
│   ├── metabolic_go_terms.txt
│   └── qc_05_go_labels.log
├── 06_figures/
│   ├── fig1a_Z_exponents_by_category_env.png/pdf
│   ├── fig1b_Z_offsets_by_category_env.png/pdf
│   ├── fig1cde_env_exponents_selected_categories.png/pdf
│   ├── fig1f_to_k_env_scatter_scaling.png/pdf
│   ├── fig1f_to_k_env_scatter_scaling_linear.png/pdf
│   ├── fig1_domains_vs_total_domains.png/pdf
│   ├── fig1_domains_vs_total_domains_linear.png/pdf
│   ├── fig1_panel_metadata.tsv
│   └── qc_06_figures.log
└── 07_extensions/
    ├── tf_mobile_scaling_params.tsv
    ├── tf_mobile_env_Z_scores.tsv
    └── qc_07_extensions.log
```

---

## Dependencies

### Python Packages
- `pandas`: Data manipulation
- `numpy`: Numerical operations
- `scipy`: Statistical functions (optional, has fallback)
- `matplotlib`: Plotting
- `seaborn`: Enhanced plotting
- `pyarrow`: Parquet file support (optional)

### Data Dependencies
- Script 01 requires: GO counts and feature matrix from previous analyses
- Script 02 requires: Output from Script 01
- Script 03 requires: Output from Script 02
- Script 04 requires: Outputs from Scripts 02 and 03
- Script 05 requires: GO ontology file and ubiquitous terms list
- Script 04.5 requires: Outputs from Scripts 01-04
- Script 06 requires: Outputs from Scripts 03, 04, and 05
- Script 07 requires: Output from Script 02

---

## Notes

1. **Prevalence Filtering**: Optional feature to filter GO terms by prevalence (e.g., 95% threshold). When used, all output files are prefixed with `prev{P}_`.

2. **Metabolic Focus**: Scripts 04.5 and 06 can generate metabolism-focused versions of figures if `metabolic_go_terms.txt` exists (created by `extract_metabolic_go_terms.py`).

3. **Test Mode**: All scripts support `--test-mode` for quick testing on small subsets.

4. **File Formats**: TSV for human-readable outputs, Parquet for fast I/O in downstream scripts.

5. **QC Logs**: Every script generates a QC log file with validation checks, summary statistics, and warnings.

6. **Figure Formats**: All figures saved as both PNG (300 DPI) and PDF for publication.

---

## Summary

The pipeline follows a clear workflow:
1. **Data Preparation** (Scripts 01-02): Join data, apply QC, filter environments
2. **Parameter Estimation** (Scripts 03-04): Fit scaling laws globally and per-environment, compute Z-scores
3. **Label Mapping** (Script 05): Map GO IDs to names for plotting
4. **Visualization** (Scripts 04.5, 06): Generate QC figures and publication figures
5. **Extensions** (Script 07): Analyze TF and mobile elements

Each script is self-contained with clear inputs/outputs and comprehensive QC logging. The figure generation scripts (04.5 and 06) are the primary visualization tools, creating both intermediate QC figures and final publication-ready figures.

