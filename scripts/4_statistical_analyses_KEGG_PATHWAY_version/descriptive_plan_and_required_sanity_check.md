# 4_statistical_analyses: Plan and Required Sanity Checks

This folder contains **only new analysis code** for the environment-stratified scaling law project.  

All existing code/results from earlier stages are treated as **read-only inputs**.

## Global Rules

- All scripts in this folder:

  - MUST live under `scripts/4_statistical_analyses/`.

  - MUST only write outputs under `results/4_statistical_analyses/`.

  - MUST NOT modify or overwrite anything under:

    - `scripts/3_*`

    - `results/3_*`

    - any other legacy directory.

- Each script:

  - Has its own subdirectory under `results/4_statistical_analyses/NN_task_name/`.

  - Writes a small QC log (`qc_*.log`) summarizing key counts, filters, and sanity checks.

  - Should not silently overwrite important files. If something may be rerun, add a suffix such as `_v1`, `_v2`, or a date stamp.

## Files and Scripts

### 01_build_master_table_env.py

**Location:**  

`scripts/4_statistical_analyses/01_build_master_table_env.py`

**Purpose:**  

Join GO counts with genome + environment metadata and apply basic QC filters.

**Reads (read-only):**

- `results/3_GO_analyses/ubiquitous_counts_table.txt`
  - Format: TSV with `Genome` column (index) + 334 GO term columns (7-digit IDs as strings)
  - Values: Integer counts ≥ 0
  - Rows: 3,089 (header + 3,088 genomes)
  - Key: `Genome` column contains accessions (GCF_* format)

- `results/3_GO_analyses/2_JGIgold_KEGG_anayses_out/05_genome_feature_matrix.tsv`
  - Format: TSV with 21 columns
  - Rows: 3,102 (header + 3,101 genomes)
  - Required columns (exact names):
    - `accession` (string) - Primary key, GCF_* format
    - `environment` (string) - Environment category from GOLD
    - `genes_total` (int) - Total genes (**n(g)** for scaling laws)
    - `genes_proteinCoding` (int) - Protein-coding genes
    - `checkm_completeness` (float) - Quality metric (0-100)
    - `checkm_contamination` (float) - Quality metric
    - `organism_name` (string) - Organism scientific name
    - `organism_taxId` (int) - NCBI taxonomy ID
    - `genome_size_mb` (float) - Genome size in megabases
    - `amino_n_burden` (float) - Amino acid nitrogen burden
    - `tf_count` (int) - Transcription factor count
    - `mobile_element_count` (int) - Mobile element count
  - Note: KEGG columns (`ko_count`, `module_count`, etc.) are 100% missing - ignore them

**Writes:**

Directory: `results/4_statistical_analyses/01_master_table/`

- `master_table_raw.tsv`  

  All genomes with both GO counts and feature matrix joined on `accession` = `Genome`.
  
  Columns:
  - All columns from `05_genome_feature_matrix.tsv` (21 columns)
  - All 334 GO term columns from `ubiquitous_counts_table.txt`
  - Key: `accession` (unique, no duplicates)

- `master_table_high_quality.tsv`  

  Subset of `master_table_raw.tsv` with filters:
  - `checkm_completeness > 90` (filter out low-quality genomes)
  - `checkm_contamination < 5` (filter out contaminated genomes)
  - `genes_total > 0` (ensure valid gene counts)
  - `environment` is not null/empty AND `environment != "Unclassified"` (require valid environment)
  
  Additional column:
  - `env_n_genomes` (int) - Per row, the total number of genomes in that environment after filtering

- `environment_counts_all.tsv`  

  Summary table with columns:
  - `environment` (string) - Environment name
  - `n_genomes` (int) - Number of genomes in this environment after QC filters

- `qc_01_master_table.log`  

  Human-readable text log with:
  - Number of genomes in GO counts file (should be 3,088)
  - Number of genomes in feature matrix (should be 3,101)
  - Number after join (`master_table_raw.tsv` row count)
  - Number after QC filters (`master_table_high_quality.tsv` row count)
  - Top 10 environments with their `n_genomes` counts
  - Any warnings about missing data or join issues

**Required sanity checks:**

- **Join verification**: Count of rows in `master_table_raw.tsv` must equal number of rows in GO table (3,088) or be within a small difference (explain any difference in log)
- **No duplicates**: Confirm no duplicated `accession` in `master_table_raw.tsv`
- **Data validity**:
  - Minimum `genes_total > 0` (verify all values are positive)
  - All GO columns are integers ≥ 0 (verify no negative counts)
  - All `checkm_completeness` values are in valid range (0-100)
- **Environment reporting**: Report any environments with < 5 genomes after QC (may want to exclude these later)
- **Missing data check**: Report how many genomes were lost due to each filter (completeness, contamination, missing environment, etc.)

---

### 02_define_env_cohorts.py

**Location:**  

`scripts/4_statistical_analyses/02_define_env_cohorts.py`

**Purpose:**  

Select environments with sufficient sample size (≥20 genomes) and create an analysis-ready master table.

**Reads:**

- `results/4_statistical_analyses/01_master_table/master_table_high_quality.tsv`
  - Expected columns: All from previous step, including `env_n_genomes`

**Writes:**

Directory: `results/4_statistical_analyses/02_env_cohorts/`

- `valid_environments_min20.tsv`  

  Columns:
  - `environment` (string) - Environment name
  - `n_genomes` (int) - Number of genomes in this environment
  
  Only rows where `n_genomes ≥ 20` (threshold is part of the spec, can be adjusted if needed)

- `master_table_env_filtered.tsv`  

  Subset of `master_table_high_quality.tsv` where `environment` is in the list from `valid_environments_min20.tsv`.
  
  Same columns as input, but only genomes from valid environments.

- `master_table_env_filtered.parquet`  

  Same data as TSV but in Parquet format for fast downstream I/O.
  
  Use pandas: `df.to_parquet(path, engine='pyarrow')`

- `qc_02_env_cohorts.log`  

  Text file summarizing:
  - Total genomes in input (`master_table_high_quality.tsv`)
  - Total genomes after environment filter (`master_table_env_filtered.tsv`)
  - Number of environments kept (from `valid_environments_min20.tsv`)
  - Full list of valid environments and `n_genomes` each (sorted by count, descending)
  - Note that all selected environments have ≥ 20 genomes

**Required sanity checks:**

- **Environment verification**: Every environment in `valid_environments_min20.tsv` must appear in `master_table_env_filtered.tsv` with matching `n_genomes`
- **Threshold enforcement**: Confirm there is **no environment** in the filtered table with `n_genomes < 20`
- **Data consistency**: Verify that `env_n_genomes` column values match the counts in `valid_environments_min20.tsv` for each environment
- **Row count check**: Sum of `n_genomes` in `valid_environments_min20.tsv` should equal row count of `master_table_env_filtered.tsv`

---

### 03_fit_global_scaling.py

**Location:**  

`scripts/4_statistical_analyses/03_fit_global_scaling.py`

**Purpose:**  

Fit the power-law scaling for each GO category using all genomes in all valid environments (global fit, not environment-specific).

**Reads:**

- `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet`
  - Required columns: `genes_total` (n(g)), all 334 GO term columns (nc(g) for each category c)

**Model:**

For each GO category column `c` (one of the 334 GO term IDs):

- **x_g** = log(n(g)) where n(g) = `genes_total` (integer)
- **y_g** = log(nc(g)) where nc(g) = count for category c in genome g
- **Zero handling**: If nc(g) = 0, either:
  - Drop that genome for this category (recommended), OR
  - Add a small pseudocount (e.g., 0.5) before log transformation
  - **Choice must be documented in QC log**

- **Fit linear regression** (OLS or Bayesian):
  ```
  y_g = β_c + α_c * x_g + ε_g
  ```
  
  Where:
  - α_c = scaling exponent (slope)
  - β_c = log-offset (intercept)
  - ε_g = error term

- **Store results**:
  - Exponent α_c with standard error (SE) or posterior standard deviation (SD)
  - 99% confidence/posterior intervals for α_c
  - Log-offset β_c with SE/SD
  - Number of genomes used (n_genomes_used)
  - R-squared (r_squared)
  - P-value (if using OLS)

**Writes:**

Directory: `results/4_statistical_analyses/03_global_scaling/`

- `global_scaling_params.tsv`  

  Columns:
  - `category` (string) - GO ID as in column name (e.g., "0000015")
  - `alpha_global` (float) - Scaling exponent α_c
  - `alpha_global_se` (float) - Standard error of exponent
  - `alpha_global_ci99_low` (float) - Lower bound of 99% confidence interval
  - `alpha_global_ci99_high` (float) - Upper bound of 99% confidence interval
  - `beta_global_log` (float) - Log-offset β_c (intercept in log space)
  - `beta_global_log_se` (float) - Standard error of log-offset
  - `n_genomes_used` (int) - Number of genomes used in fit (after filtering zeros if applicable)
  - `r_squared` (float) - R² of the linear fit
  - `p_value` (float) - P-value from OLS regression (or NaN if Bayesian)

- `global_scaling_params.parquet` (optional, for speed)

- `qc_03_global_scaling.log`  

  Summaries:
  - Number of genomes per category used (should be similar across categories, ~2,000-3,000)
  - Statistics of `alpha_global`:
    - Min, max, median, IQR (interquartile range)
  - Count of categories with:
    - `alpha_global < 0` (negative scaling - unusual)
    - `alpha_global > 4` (very steep scaling - unusual)
    - `alpha_global_se > 1` (high uncertainty - may indicate poor fit)
  - Simple text note about how zeros in nc(g) were handled (dropped vs pseudocount)
  - Sample of 5-10 categories with their exponents and R² values

**Required sanity checks:**

- **Minimum sample size**: At least 10 genomes used for each category; if fewer, log which categories were skipped
- **Exponent bounds**: No category should have absurd exponents without a comment (e.g., `alpha_global = 20` → must be flagged in log)
- **Data quality**: Verify that `n_genomes_used` is reasonable (should be close to total genomes minus those with zero counts)
- **R-squared check**: Most categories should have R² > 0.5 (good fit); log categories with R² < 0.3 (poor fit)

---

### 04_fit_env_scaling_and_Z.py

**Location:**  

`scripts/4_statistical_analyses/04_fit_env_scaling_and_Z.py`

**Purpose:**  

Fit per-environment scaling laws, then compute Z-scores comparing environment-specific parameters to global parameters.

**Reads:**

1. `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet`
   - Contains: `environment`, `genes_total`, all 334 GO term columns

2. `results/4_statistical_analyses/03_global_scaling/global_scaling_params.tsv`
   - Contains: `category`, `alpha_global`, `alpha_global_se`, `beta_global_log`, `beta_global_log_se`

3. `results/4_statistical_analyses/02_env_cohorts/valid_environments_min20.tsv`
   - Contains: `environment`, `n_genomes` (list of valid environments)

**Per-environment fit:**

For each environment `i` and GO category `c`:

- Subset rows where `environment == i`
- Filter: Use only genomes with `n(g) > 0` AND `nc(g) > 0` (drop zeros)
- **Require `n_genomes_used ≥ 10`**. If fewer, skip this env×category combination and log it
- Fit: `y = β_{i,c} + α_{i,c} * x + ε` (same model as global fit)
- Store: `alpha_env`, `alpha_env_se`, `beta_env_log`, `beta_env_log_se`, `n_genomes_used`, `r_squared`, `p_value`

**Z-statistics:**

For each env×category pair with a successful fit:

- From global params, get: `alpha_global`, `alpha_global_se`, `beta_global_log`, `beta_global_log_se`
- From env fit, get: `alpha_env`, `alpha_env_se`, `beta_env_log`, `beta_env_log_se`

- **Compute Z-scores**:
  ```
  Z_{i,c}^{(α)} = (α_{i,c} - α_c) / sqrt(σ²_{α_{i,c}} + σ²_{α_c})
  
  Z_{i,c}^{(β)} = (β_{i,c}^{(log)} - β_c^{(log)}) / sqrt(σ²_{β_{i,c}} + σ²_{β_c})
  ```
  
  Where σ² = SE² (variance = squared standard error)

- **Aggregate across environments** (for each category c):
  ```
  Z_c^{(α)} = sqrt((1/N_c) * Σ_i (Z_{i,c}^{(α)})²)
  
  Z_c^{(β)} = sqrt((1/N_c) * Σ_i (Z_{i,c}^{(β)})²)
  ```
  
  Where N_c = number of environments where category c had a successful fit

**Writes:**

Directory: `results/4_statistical_analyses/04_env_scaling/`

- `env_scaling_params.tsv`  

  Columns:
  - `environment` (string) - Environment name
  - `category` (string) - GO term ID (e.g., "0000015")
  - `alpha_env` (float) - Environment-specific exponent
  - `alpha_env_se` (float) - Standard error of exponent
  - `alpha_env_ci99_low` (float) - Lower bound of 99% CI
  - `alpha_env_ci99_high` (float) - Upper bound of 99% CI
  - `beta_env_log` (float) - Environment-specific log-offset
  - `beta_env_log_se` (float) - Standard error of log-offset
  - `n_genomes_used` (int) - Number of genomes used in fit
  - `r_squared` (float) - R² of fit
  - `p_value` (float) - P-value from regression

- `env_vs_global_Z_scores.tsv`  

  Columns:
  - `environment` (string)
  - `category` (string)
  - `Z_alpha` (float) - Z-score for exponent deviation
  - `Z_beta` (float) - Z-score for offset deviation
  - `n_genomes_used` (int) - Number of genomes used

- `category_Z_summary.tsv`  

  Columns:
  - `category` (string) - GO term ID
  - `Z_alpha_category` (float) - Aggregated Z-score for exponent variation
  - `Z_beta_category` (float) - Aggregated Z-score for offset variation
  - `n_envs_used` (int) - Number of environments where this category had a fit

- `qc_04_env_scaling.log`  

  Text summary:
  - Number of env×category fits successfully completed
  - Number of env×category combinations skipped because `n_genomes_used < 10` (list a few examples)
  - Summary statistics (min, max, median) for `Z_alpha` and `Z_beta` across all env×category pairs
  - Top 10 categories by `Z_alpha_category` (most variable exponents)
  - Top 10 categories by `Z_beta_category` (most variable offsets)
  - Note: Most categories should have `Z_alpha_category < 2` (consistent across environments)

**Required sanity checks:**

- **No invalid values**: Confirm no `Z_alpha` or `Z_beta` is NaN or infinite (check for division by zero, missing data)
- **Environment count**: Confirm `n_envs_used ≤ number of environments` and `≥ 1` for categories present in summary
- **Z-score distribution**: At least some categories should have `Z_alpha_category > 2` (variable) and most should be < 2 (consistent). If everything > 2 or everything ≈ 0, something is wrong
- **Fit quality**: Check that `n_genomes_used` values are reasonable (≥ 10, typically 20-100+ per environment)
- **Coverage**: Report how many env×category combinations were successfully fitted vs skipped

---

### 05_map_go_labels.py

**Location:**  

`scripts/4_statistical_analyses/05_map_go_labels.py`

**Purpose:**  

Map GO IDs (7-digit format) to human-readable names and definitions for plotting and reporting.

**Reads:**

1. `results/3_GO_analyses/ubiquitous_terms.txt`
   - Format: One GO term ID per line (7-digit format, e.g., "0000015")
   - Rows: 334 (one per ubiquitous term)

2. `data/go/go-basic.obo`
   - Format: OBO (Open Biomedical Ontology) format
   - Size: ~550,491 lines
   - Contains: Complete GO ontology with term definitions, relationships, names

**Parsing OBO file:**

- Look for `[Term]` sections
- Extract:
  - `id: GO:0000015` (convert to 7-digit format by removing "GO:" prefix)
  - `name: mitotic cell cycle` (term name)
  - `namespace: biological_process` (or `molecular_function`, `cellular_component`)
  - `def: "..."` (definition text)

**Writes:**

Directory: `results/4_statistical_analyses/05_go_labels/`

- `go_term_labels.tsv`  

  Columns:
  - `category` (string) - 7-digit GO ID (e.g., "0000015")
  - `go_id` (string) - Full GO ID with prefix (e.g., "GO:0000015")
  - `name` (string) - GO term name (e.g., "mitotic cell cycle")
  - `namespace` (string) - One of: "biological_process", "molecular_function", "cellular_component"
  - `definition` (string) - Full definition text from OBO file

- `go_term_labels_for_plots.tsv`  

  Columns:
  - `category` (string) - 7-digit GO ID
  - `short_label` (string) - Trimmed or cleaned name, ≤ ~40 characters (for plot labels)
  - `super_category` (string, optional) - Manual grouping (e.g., "Metabolism", "Regulation", "Transport"). Can be empty if not used

- `qc_05_go_labels.log`  

  Text file with:
  - Number of GO IDs in input (`ubiquitous_terms.txt` - should be 334)
  - Number of GO IDs successfully mapped (should be 334)
  - Listing of any IDs not found in OBO file (should be zero - all ubiquitous terms should be in ontology)
  - 5-10 sample mappings printed in human-readable form (e.g., "0000015 → GO:0000015: mitotic cell cycle")

**Required sanity checks:**

- **Count match**: Number of rows in `go_term_labels.tsv` must equal number of lines in `ubiquitous_terms.txt` (334)
- **No duplicates**: No duplicate `category` values in `go_term_labels.tsv`
- **Complete mapping**: All 334 GO IDs should be found in OBO file; if any are missing, log them and investigate
- **Format validation**: Verify that `go_id` column has format "GO:XXXXXXX" (7 digits after "GO:")

---

### 06_make_scaling_figures.py

**Location:**  

`scripts/4_statistical_analyses/06_make_scaling_figures.py`

**Purpose:**  

Generate environment-stratified version of Fig. 1 (panels a-k) from van Nimwegen (2003), showing Z-statistics, exponent comparisons, and scatter plots.

**Reads:**

1. `results/4_statistical_analyses/04_env_scaling/category_Z_summary.tsv`
   - Contains: `category`, `Z_alpha_category`, `Z_beta_category`, `n_envs_used`

2. `results/4_statistical_analyses/04_env_scaling/env_scaling_params.tsv`
   - Contains: Per-environment exponents and offsets for each category

3. `results/4_statistical_analyses/03_global_scaling/global_scaling_params.tsv`
   - Contains: Global exponents and offsets for each category

4. `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet`
   - Contains: Raw data for scatter plots (`genes_total`, GO term counts, `environment`)

5. `results/4_statistical_analyses/05_go_labels/go_term_labels_for_plots.tsv`
   - Contains: `category`, `short_label` for plot labels

**Figure panels:**

- **Panel 1a**: Z-statistics for exponents by category
  - X-axis: GO categories (sorted by Z_alpha_category)
  - Y-axis: Z_alpha_category values
  - Horizontal line at Z=2 (threshold for significance)
  - Similar to van Nimwegen Fig. 1a

- **Panel 1b**: Z-statistics for offsets by category
  - X-axis: GO categories (sorted by Z_beta_category)
  - Y-axis: Z_beta_category values
  - Horizontal line at Z=2
  - Similar to van Nimwegen Fig. 1b

- **Panels 1c-1e**: Exponent comparisons for selected categories
  - Show fitted exponents (dots) with 99% CI (vertical bars) across all environments
  - Select 3 interesting categories (e.g., one with low Z, one with medium Z, one with high Z)
  - Horizontal dotted line shows global exponent
  - Similar to van Nimwegen Fig. 1c-e

- **Panels 1f-1k**: Scatter plots with fits
  - For each selected category (from panels c-e), show scatter plots for 2-3 representative environments
  - X-axis: log(n(g)) = log(genes_total)
  - Y-axis: log(nc(g)) = log(GO term count)
  - Gray dots: All genomes (all environments)
  - Colored dots: Genomes from specific environment
  - Colored line: Environment-specific fit
  - Gray line: Global fit
  - Similar to van Nimwegen Fig. 1f-k

**Writes:**

Directory: `results/4_statistical_analyses/06_figures/`

- `fig1a_Z_exponents_by_category_env.png` / `.pdf`
  - Panel 1a: Z-statistics for exponents

- `fig1b_Z_offsets_by_category_env.png` / `.pdf`
  - Panel 1b: Z-statistics for offsets

- `fig1cde_env_exponents_selected_categories.png` / `.pdf`
  - Panels 1c, 1d, 1e: Exponent comparisons

- `fig1f_to_k_env_scatter_scaling.png` / `.pdf`
  - Panels 1f through 1k: Scatter plots with fits

- `fig1_panel_metadata.tsv`
  
  Columns:
  - `panel` (string) - Panel identifier (e.g., "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k")
  - `category` (string) - GO term ID used in this panel
  - `short_label` (string) - Human-readable label for category
  - `environment` (string, optional) - Environment name (for scatter panels f-k)
  - `notes` (string, optional) - Any extra notes about this panel

- `qc_06_figures.log`
  
  Text summary:
  - Categories used for panels c-e (list the 3 selected categories and why)
  - Environment & category combinations used for panels f-k (list which environments shown for each category)
  - Confirmation that all figure files were written and are non-empty (check file sizes)
  - Any warnings about missing data or plotting issues

**Required sanity checks:**

- **Data consistency**: For one category, verify that points in panels c-e match values from `env_scaling_params.tsv` (spot check 2-3 values)
- **Fit verification**: For one panel (e.g., panel f), visually confirm that the global vs environment-specific regression lines match numeric exponents and offsets from parameter tables
- **File creation**: All figure files should be created and have non-zero file size
- **Panel metadata**: All panels should be documented in `fig1_panel_metadata.tsv`

---

### 07_scaling_extensions_tf_mobile_nutrient.py

**Location:**  

`scripts/4_statistical_analyses/07_scaling_extensions_tf_mobile_nutrient.py`

**Purpose:**  

Extend scaling analyses to transcription factors (`tf_count`), mobile elements (`mobile_element_count`), and optionally nutrient-related GO categories (amino acid metabolism, etc.).

**Reads:**

- `results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet`
  - Required columns: `environment`, `genes_total`, `tf_count`, `mobile_element_count`, all GO term columns

**Analysis:**

Treat `tf_count` and `mobile_element_count` as additional "functional categories" and analyze their scaling:

- **Global fits**: Fit power laws for `tf_count` and `mobile_element_count` vs `genes_total` (same model as GO categories)
- **Environment-specific fits**: Fit per-environment scaling laws
- **Z-scores**: Compute Z-scores comparing environment-specific to global parameters
- **Optional**: Identify nutrient-related GO categories (e.g., amino acid metabolism) and analyze their scaling separately

**Writes:**

Directory: `results/4_statistical_analyses/07_extensions/`

- `tf_mobile_scaling_params.tsv`

  Columns:
  - `category` (string) - Either "tf_count" or "mobile_element_count"
  - `environment` (string) - Environment name, or "global" for global fit
  - `alpha` (float) - Scaling exponent
  - `alpha_se` (float) - Standard error
  - `beta_log` (float) - Log-offset
  - `beta_log_se` (float) - Standard error of offset
  - `n_genomes_used` (int)
  - `r_squared` (float)
  - `p_value` (float)

- `tf_mobile_env_Z_scores.tsv`

  Columns:
  - `category` (string) - "tf_count" or "mobile_element_count"
  - `environment` (string)
  - `Z_alpha` (float) - Z-score for exponent
  - `Z_beta` (float) - Z-score for offset
  - `n_genomes_used` (int)

- `tf_scaling_by_environment.png` / `.pdf` (optional)
  - Scatter plot showing TF count scaling across environments

- `mobile_scaling_by_environment.png` / `.pdf` (optional)
  - Scatter plot showing mobile element count scaling across environments

- `qc_07_extensions.log`

  Text summary:
  - Global exponent for `tf_count` (should be ~>1, often ~2 in literature)
  - Global exponent for `mobile_element_count`
  - Summary of environment-specific variation (Z-scores)
  - Any notes about outliers or unusual patterns

**Required sanity checks:**

- **Data validity**: Ensure `tf_count` and `mobile_element_count` are non-negative integers
- **Biological reasonableness**: Global exponent for TF counts should be biologically reasonable (>1, typically ~2 for quadratic scaling). Log observed values and comment in QC log
- **Outlier check**: Check that exponents aren't driven solely by a few outlier genomes (plot residuals or check for high-leverage points)
- **Coverage**: Verify that all environments have sufficient data (≥10 genomes) for fits

---

## Minimal Execution Order

1. `01_build_master_table_env.py` - Join data and apply QC filters
2. `02_define_env_cohorts.py` - Select environments with ≥20 genomes
3. `03_fit_global_scaling.py` - Fit global scaling laws for all GO categories
4. `04_fit_env_scaling_and_Z.py` - Fit per-environment scaling and compute Z-scores
5. `05_map_go_labels.py` - Map GO IDs to names (can run in parallel with 03-04)
6. `06_make_scaling_figures.py` - Generate figures
7. `07_scaling_extensions_tf_mobile_nutrient.py` - Extensions (optional)

**Dependencies:**
- Scripts 01-04 must run in order
- Script 05 can run anytime after GO terms are identified (can run in parallel)
- Script 06 requires outputs from 03, 04, and 05
- Script 07 requires output from 02 (can run in parallel with 03-06)

Running these scripts in order and reviewing each `qc_*.log` should reproduce the environment-stratified scaling analyses end-to-end.

## Notes

- Thresholds (≥20 genomes per environment, ≥10 genomes per env×category fit) can be adjusted if needed, but document any changes in the QC logs
- The "≥10 genomes" threshold for env×category fits is a minimum; most should have 20-100+ genomes
- Z-score threshold of 2 is conventional but can be adjusted for multiple testing correction if needed
- All file paths are relative to the repository root (`/n/scratch/users/b/byc014/github/bac_genome_constraint/`)

