# KEGG Pathway Statistical Analyses Pipeline - Execution Guide

## Script Classification

### ✅ MAIN SCRIPTS (Required for Core Analysis)

These scripts are **essential** and must be run in sequence to generate statistical results:

1. **Script 01**: `01_build_master_table_env.py`
2. **Script 02**: `02_define_env_cohorts.py`
3. **Script 03**: `03_fit_global_scaling.py`
4. **Script 04**: `04_fit_env_scaling_and_Z_single.py` (via `04_fit_env_scaling_and_Z.sh`)
5. **Script 05**: `05_map_kegg_labels_single.py` (via `05_map_kegg_labels.sh`) - **Recommended but optional**

### 📊 OPTIONAL SCRIPTS (Visualization & Extensions)

These scripts generate figures and extended analyses but are not required for core statistics:

1. **Script 04.5**: `04.5_plot_intermediate_figures_single.py` (via `04.5_plot_intermediate_figures.sh`)
   - Generates QC and diagnostic plots
   - Alternative: `04.5_plot_minimal_single.py` (minimal version)

2. **Script 06**: `06_make_scaling_figures_single.py` (via `06_make_scaling_figures.sh`)
   - Generates publication-quality figures

3. **Script 07**: `07_scaling_extensions_single.py` (via `07_scaling_extensions.sh`)
   - Analyzes transcription factors and mobile elements

---

## Execution Order & Sequence

### Core Pipeline (Required)

```bash
cd scripts/4.5_statistical_analyses_KEGG_PATHWAY_version

# Step 1: Build master tables (processes all 3 data types automatically)
python 01_build_master_table_env.py

# Step 2: Filter environments (≥20 genomes)
python 02_define_env_cohorts.py

# Step 3: Fit global scaling laws
python 03_fit_global_scaling.py

# Step 4: Fit per-environment scaling and compute Z-scores
bash 04_fit_env_scaling_and_Z.sh

# Step 5: Map KEGG IDs to human-readable names (recommended for plots)
bash 05_map_kegg_labels.sh
```

### With Test Mode (Fast Testing)

```bash
python 01_build_master_table_env.py --test-mode
python 02_define_env_cohorts.py --test-mode
python 03_fit_global_scaling.py --test-mode
bash 04_fit_env_scaling_and_Z.sh --test-mode
bash 05_map_kegg_labels.sh --test-mode
```

### With Prevalence Filtering

```bash
python 01_build_master_table_env.py --prevalence-threshold 99
python 02_define_env_cohorts.py --prevalence-threshold 99
python 03_fit_global_scaling.py --prevalence-threshold 99
bash 04_fit_env_scaling_and_Z.sh --prevalence-threshold 99
bash 05_map_kegg_labels.sh --prevalence-threshold 99
```

### Optional Visualization (After Core Pipeline)

```bash
# Generate intermediate QC figures
bash 04.5_plot_intermediate_figures.sh

# Generate publication figures
bash 06_make_scaling_figures.sh

# Analyze extensions (TF, mobile elements)
bash 07_scaling_extensions.sh
```

---

## Key Methods by Script

### Script 01: Build Master Table (`01_build_master_table_env.py`)

**Purpose**: Integrate KEGG counts with genome metadata and apply QC filters

**Key Methods/Steps**:
1. **Input Validation**
   - Loads KEGG count tables (reactions, KOs, pathways) from `results/3.5_KEGG_n_reaction_analyses/`
   - Loads genome feature matrix with metadata
   - Validates column structure and data types

2. **Join Integrity Checks**
   - Performs inner join on `accession` = `Genome`
   - Verifies join size matches expected (~3088 genomes)
   - Checks for duplicate accessions

3. **QC Filtering**
   - **Completeness filter**: `checkm_completeness > 90` (NaN → 0)
   - **Contamination filter**: `checkm_contamination < 5` (NaN → 100)
   - **Genes filter**: `genes_total > 0`
   - **Environment filter**: Not null, not empty, ≠ "Unclassified"

4. **Environment-Level Sanity Checks**
   - Computes genome counts per environment
   - Flags environments with < 5 or > 1000 genomes

5. **Output Generation**
   - `master_table_raw_{suffix}.tsv` - Before QC (after join)
   - `master_table_high_quality_{suffix}.tsv` - After QC filters
   - `environment_counts_all_{suffix}.tsv` - Environment summary
   - `qc_01_master_table_{suffix}.log` - QC log

**Data Types Processed**: All three (reactions, ko, pathway) automatically

---

### Script 02: Define Environment Cohorts (`02_define_env_cohorts.py`)

**Purpose**: Select environments with sufficient sample size (≥20 genomes)

**Key Methods/Steps**:
1. **Consistency Verification**
   - Loads high-quality master table from Script 01
   - Recomputes environment counts as sanity check
   - Verifies `env_n_genomes` column matches

2. **Environment Threshold Filter**
   - Keeps environments with ≥ 20 genomes
   - Creates `valid_environments_min20_{suffix}.tsv`
   - Optionally creates `valid_environments_min10_{suffix}.tsv` if < 5 environments

3. **KEGG Column Filtering**
   - For pathways: Explicitly drops "rn" columns (keeps only "map" pathways)
   - Applies prevalence filtering if `--prevalence-threshold` specified
   - Uses `filter_kegg_columns_by_prevalence()` from `prevalence_utils.py`

4. **Output Generation**
   - `valid_environments_min20_{suffix}.tsv` - Valid environments list
   - `master_table_env_filtered_{suffix}.tsv` - Filtered master table
   - `master_table_env_filtered_{suffix}.parquet` - Parquet format (faster loading)
   - `qc_02_env_cohorts_{suffix}.log` - QC log

**Data Types Processed**: All three (reactions, ko, pathway) automatically

---

### Script 03: Fit Global Scaling (`03_fit_global_scaling.py`)

**Purpose**: Fit power-law scaling relationships: `nc(g) = β × n(g)^α` across all environments

**Key Methods/Steps**:
1. **Data Loading**
   - Loads filtered master table from Script 02
   - Identifies KEGG columns using `get_kegg_columns()` or `filter_kegg_columns_by_prevalence()`
   - Optionally loads KEGG labels for better category names

2. **Pre-Regression Validation**
   - For each category, drops rows with `genes_total <= 0` or `nc(g) <= 0`
   - Tracks statistics: zero drop counts, genomes per category

3. **OLS Regression Fitting**
   - **Model**: `log(nc(g)) = log(β) + α × log(n(g))`
   - **Method**: Ordinary Least Squares (OLS) in log-log space
   - **Implementation**: Uses `scipy.stats.linregress()` if available, else manual OLS
   - **Validation**:
     - Minimum sample size: ≥ 10 genomes
     - Minimum variance: `std(log(n(g))) >= 0.05`
   - **Outputs per category**:
     - `alpha_global` (scaling exponent)
     - `alpha_global_se` (standard error)
     - `alpha_global_ci99_low/high` (99% confidence intervals)
     - `beta_global_log` (log intercept)
     - `beta_global_log_se` (standard error)
     - `beta_global_ci99_low/high` (99% confidence intervals)
     - `r_squared` (coefficient of determination)
     - `p_value` (t-test p-value)
     - `n_genomes_used` (sample size)

4. **Global Summary QC**
   - Distribution statistics for `alpha_global`
   - QC flags: negative alpha, high alpha (>4), high SE (>1), low R² (<0.2)
   - Top 10 categories by alpha

5. **Output Generation**
   - `global_scaling_params_{suffix}.tsv` - All fitted parameters
   - `global_scaling_params_{suffix}.parquet` - Parquet format
   - `qc_03_global_scaling_{suffix}.log` - QC log

**Data Types Processed**: All three (reactions, ko, pathway) automatically

**Mathematical Model**:
- **Linear space**: `nc(g) = β × n(g)^α`
- **Log-log space**: `log(nc(g)) = log(β) + α × log(n(g))`
- **Interpretation**:
  - `α = 1`: Linear scaling
  - `α < 1`: Sub-linear scaling (core genes)
  - `α > 1`: Super-linear scaling (regulatory genes)

---

### Script 04: Environment-Specific Fits & Z-Scores (`04_fit_env_scaling_and_Z_single.py`)

**Purpose**: Fit per-environment scaling laws and compute Z-scores comparing environment-specific to global parameters

**Key Methods/Steps**:
1. **Data Loading**
   - Loads filtered master table from Script 02
   - Loads global scaling parameters from Script 03
   - Loads valid environments list from Script 02

2. **Per-Environment Scaling Fits**
   - For each environment × category combination:
     - Filters: `genes_total > 0` and `nc(g) > 0`
     - Minimum sample size: ≥ 10 genomes
     - Minimum variance: `std(log(n(g))) >= 0.05`
     - Fits OLS regression: `log(nc(g)) = log(β_env) + α_env × log(n(g))`
     - Computes 99% confidence intervals using t-distribution
   - **Outputs per fit**:
     - `alpha_env`, `alpha_env_se`, `alpha_env_ci99_low/high`
     - `beta_env_log`, `beta_env_log_se`
     - `n_genomes_used`, `r_squared`, `p_value`

3. **Z-Score Computation**
   - **Z-Score for Exponents**:
     ```
     Z_alpha = (alpha_env - alpha_global) / sqrt(SE(alpha_env)² + SE(alpha_global)²)
     ```
   - **Z-Score for Offsets**:
     ```
     Z_beta = (beta_env_log - beta_global_log) / sqrt(SE(beta_env_log)² + SE(beta_global_log)²)
     ```
   - **Interpretation**:
     - `|Z| > 2`: Significant deviation (p < 0.05)
     - `|Z| > 3`: Highly significant deviation (p < 0.01)

4. **Category-Level Aggregation**
   - For each category, computes summary Z-statistics:
     - `Z_alpha_category = sqrt(mean(Z_alpha²))` (root mean square across environments)
     - `Z_beta_category = sqrt(mean(Z_beta²))`
     - `n_envs_used` (number of environments with valid fits)

5. **Output Generation**
   - `env_scaling_params_{suffix}.tsv` - All environment×category fits
   - `env_vs_global_Z_scores_{suffix}.tsv` - Z-scores for each env×category
   - `category_Z_summary_{suffix}.tsv` - Category-level aggregated Z-scores
   - `qc_04_env_scaling_{suffix}.log` - QC log

**Data Types Processed**: One at a time (use wrapper `04_fit_env_scaling_and_Z.sh` for all)

**Wrapper Script**: `04_fit_env_scaling_and_Z.sh` loops over all three data types

---

### Script 05: Map KEGG Labels (`05_map_kegg_labels_single.py`)

**Purpose**: Fetch human-readable names and definitions for KEGG IDs from KEGG API

**Key Methods/Steps**:
1. **KEGG Column Identification**
   - Loads master table or global scaling params to get KEGG category IDs
   - Uses `get_kegg_columns()` to identify relevant columns

2. **KEGG API Fetching**
   - **Method**: `fetch_kegg_info(kegg_id, data_type, use_cache=True)`
   - **Caching**: Checks `results/3.5_KEGG_n_reaction_analyses/kegg_api_cache/` first
   - **API Endpoints**:
     - Reactions: `http://rest.kegg.jp/get/rn:{kegg_id}`
     - KOs: `http://rest.kegg.jp/get/ko:{kegg_id}`
     - Pathways: `http://rest.kegg.jp/get/{kegg_id}` (already has map/rn prefix)
   - **Rate Limiting**: 0.5 second delay between API calls
   - **Error Handling**: Retries on timeout, skips on 404

3. **Label Parsing**
   - Extracts `NAME` field from KEGG response
   - Extracts `DEFINITION` field (if available)
   - Handles multi-line definitions

4. **Output Generation**
   - `kegg_term_labels_{suffix}.tsv` - Two columns: `category`, `name`
   - Optionally includes `definition` column if available
   - Cached API responses saved for future runs

**Data Types Processed**: One at a time (use wrapper `05_map_kegg_labels.sh` for all)

**Wrapper Script**: `05_map_kegg_labels.sh` loops over all three data types

**Note**: This script is **recommended** but not strictly required. Without it, plots will show KEGG IDs instead of names.

---

## Optional Scripts Summary

### Script 04.5: Plot Intermediate Figures

**Purpose**: Generate QC and diagnostic plots for each analysis step

**Key Outputs**:
- Script 01 figures: QC filtering flowchart, completeness/contamination distributions, genome size distributions
- Script 02 figures: Final environments, environment contribution, genome size by environment
- Script 03 figures: Alpha distribution, alpha vs R², representative scaling plots
- Script 04 figures: Z-score heatmaps, Z-score distributions, exception categories

**Alternative**: `04.5_plot_minimal_single.py` - Minimal version with essential plots only

### Script 06: Make Scaling Figures

**Purpose**: Generate publication-quality figures following van Nimwegen (2003) structure

**Key Outputs**:
- Z-statistics by category (panels 1a-1b)
- Exponent comparisons for selected categories (panels 1c-1e)
- Scatter plots with fits (panels 1f-1k)
- Domain vs total domains plots

### Script 07: Scaling Extensions

**Purpose**: Analyze transcription factors and mobile elements scaling

**Key Methods**:
- Fits scaling laws for `tf_count` and `mobile_element_count`
- Tests hypothesis: `α_TF ≈ 2` (quadratic scaling for regulatory complexity)
- Computes environment-specific parameters and Z-scores

---

## Utility Functions (`prevalence_utils.py`)

Key helper functions used across scripts:

1. **`get_prevalence_prefix(prevalence_threshold)`**: Returns prefix string (e.g., "prev99_")
2. **`detect_kegg_data_type(df)`**: Auto-detects data type from column patterns
3. **`get_kegg_columns(df, data_type)`**: Returns list of KEGG columns for given type
4. **`filter_kegg_columns_by_prevalence(df, threshold, data_type)`**: Filters columns by prevalence
5. **`get_kegg_input_files(data_type, threshold, base_dir)`**: Returns correct input file path
6. **`get_data_type_suffix(data_type)`**: Returns suffix string (e.g., "_reactions")

---

## Data Flow Summary

```
Raw KEGG Count Tables (results/3.5_KEGG_n_reaction_analyses/)
    ↓
Script 01: Build Master Table
    ↓
Script 02: Filter Environments (≥20 genomes)
    ↓
Script 03: Fit Global Scaling
    ↓
Script 04: Fit Per-Environment Scaling + Z-Scores
    ↓
Script 05: Map KEGG Labels (optional, for better plots)
    ↓
Script 04.5/06: Generate Figures (optional)
    ↓
Script 07: Extensions (optional)
```

---

## Expected Outputs

After running the core pipeline (Scripts 01-04), you should have:

- **2,164 high-quality genomes** across **8 environments**
- **~448 reactions**, **~684 KOs**, **~198 pathways** with successful global fits
- **~2,668 reaction fits**, **~4,617 KO fits**, **~1,420 pathway fits** for environment-specific scaling
- **Z-scores** quantifying environment-specific variation for each category

---

## Notes

- All scripts process **all three data types** (reactions, KOs, pathways) automatically, except Scripts 04-07 which use wrapper scripts
- Use `--test-mode` for fast testing on small subsets
- Use `--prevalence-threshold 99` to focus on ubiquitous terms (present in ≥99% of genomes)
- Script 05 can be run at any time after Script 01 (it only needs category IDs)
- Scripts 04.5, 06, and 07 require Script 05 to be run first for human-readable labels in plots


