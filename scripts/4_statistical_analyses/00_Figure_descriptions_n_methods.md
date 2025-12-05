# Figure Descriptions and Methods

This document provides a comprehensive description of all figures, panels, mathematical models, and methodological justifications for the statistical analyses pipeline.

---

## Table of Contents

1. [Mathematical Framework](#mathematical-framework)
2. [Data Processing Pipeline](#data-processing-pipeline)
3. [Script 01: Master Table & QC Figures](#script-01-master-table--qc-figures)
4. [Script 02: Environment Cohort Figures](#script-02-environment-cohort-figures)
5. [Script 03: Global Scaling Figures](#script-03-global-scaling-figures)
6. [Script 04: Environment-Specific Scaling & Z-Score Figures](#script-04-environment-specific-scaling--z-score-figures)
7. [Script 06: Publication Figures (Main Results)](#script-06-publication-figures-main-results)
8. [Script 07: Extensions (TF & Mobile Elements)](#script-07-extensions-tf--mobile-elements)

---

## Mathematical Framework

### Power-Law Scaling Model

The core analysis models the relationship between genome size and gene content using a power-law scaling relationship:

**Linear Space:**
\[
n_c(g) = \beta \times n(g)^{\alpha}
\]

**Log-Log Space (for regression):**
\[
\log(n_c(g)) = \log(\beta) + \alpha \times \log(n(g))
\]

Where:
- \(n(g)\) = total number of genes in genome \(g\) (`genes_total`)
- \(n_c(g)\) = number of genes in functional category \(c\) in genome \(g\) (GO term count)
- \(\alpha\) = scaling exponent (slope in log-log space)
- \(\beta\) = scaling prefactor/offset (intercept in log-log space, stored as \(\log(\beta)\))

### Interpretation of Scaling Exponents

- **\(\alpha = 1\)**: Linear scaling — category grows proportionally with genome size
- **\(\alpha < 1\)**: Sub-linear scaling — category grows slower than genome size (e.g., core genes)
- **\(\alpha > 1\)**: Super-linear scaling — category grows faster than genome size (e.g., regulatory genes)
- **\(\alpha \approx 2\)**: Quadratic scaling — category scales with genome size squared (e.g., transcription factors)

### Ordinary Least Squares (OLS) Regression

For each GO category, we fit the model in log-log space using OLS:

\[
\hat{\alpha} = \frac{\sum_{i=1}^{n} (x_i - \bar{x})(y_i - \bar{y})}{\sum_{i=1}^{n} (x_i - \bar{x})^2}
\]

\[
\hat{\beta}_{\text{log}} = \bar{y} - \hat{\alpha} \times \bar{x}
\]

Where:
- \(x_i = \log(n(g_i))\), \(y_i = \log(n_c(g_i))\)
- \(\bar{x}, \bar{y}\) are sample means
- \(n\) = number of genomes

**Standard Errors:**
\[
\text{SE}(\hat{\alpha}) = \sqrt{\frac{\text{MSE}}{S_{xx}}}, \quad \text{SE}(\hat{\beta}_{\text{log}}) = \sqrt{\text{MSE} \times \left(\frac{1}{n} + \frac{\bar{x}^2}{S_{xx}}\right)}
\]

Where:
- \(\text{MSE} = \frac{\text{SS}_{\text{res}}}{n-2}\) (mean squared error)
- \(S_{xx} = \sum_{i=1}^{n} (x_i - \bar{x})^2\)
- \(\text{SS}_{\text{res}} = \sum_{i=1}^{n} (y_i - \hat{y}_i)^2\) (residual sum of squares)

**99% Confidence Intervals:**
\[
\alpha \in [\hat{\alpha} - t_{0.995} \times \text{SE}(\hat{\alpha}), \hat{\alpha} + t_{0.995} \times \text{SE}(\hat{\alpha})]
\]

Where \(t_{0.995}\) is the 99th percentile of the t-distribution (approximately 2.576 for large \(n\)).

### Z-Statistics for Environment Comparison

To quantify whether environment-specific scaling parameters differ significantly from global parameters, we compute Z-scores:

**Z-Score for Exponents:**
\[
Z_{\alpha} = \frac{\alpha_{\text{env}} - \alpha_{\text{global}}}{\sqrt{\text{SE}(\alpha_{\text{env}})^2 + \text{SE}(\alpha_{\text{global}})^2}}
\]

**Z-Score for Offsets:**
\[
Z_{\beta} = \frac{\beta_{\text{env,log}} - \beta_{\text{global,log}}}{\sqrt{\text{SE}(\beta_{\text{env,log}})^2 + \text{SE}(\beta_{\text{global,log}})^2}}
\]

**Interpretation:**
- \(|Z| > 2\): Significant deviation from global parameter (approximately \(p < 0.05\))
- \(|Z| > 3\): Highly significant deviation (approximately \(p < 0.01\))

**Category-Level Z-Statistics:**
For each GO category, we compute summary statistics across environments:
- \(Z_{\alpha,\text{category}} = \sqrt{\frac{1}{n_{\text{envs}}} \sum_{j=1}^{n_{\text{envs}}} Z_{\alpha,j}^2}\) (root mean square)
- \(Z_{\beta,\text{category}} = \sqrt{\frac{1}{n_{\text{envs}}} \sum_{j=1}^{n_{\text{envs}}} Z_{\beta,j}^2}\)

These quantify the overall variability of environment-specific parameters relative to the global fit.

---

## Data Processing Pipeline

### Quality Control Filters

1. **CheckM Completeness:** ≥ 90% (high-quality genomes)
2. **CheckM Contamination:** ≤ 5%
3. **Gene Count:** `genes_total > 0` and GO category count > 0 (for regression)
4. **Environment Filtering:** Only environments with ≥20 genomes (for reliable per-environment fits)
5. **GO Term Prevalence:** Optional filter for terms present in ≥X% of genomes (e.g., 95%, 99%)

### Data Flow

```
Raw Data (NCBI genomes)
    ↓
Script 01: Build master table + QC
    ↓
Script 02: Filter environments (≥20 genomes)
    ↓
Script 03: Fit global scaling (all environments combined)
    ↓
Script 04: Fit per-environment scaling + compute Z-scores
    ↓
Script 05: Map GO term IDs to labels
    ↓
Script 06: Generate publication figures
```

---

## Script 01: Master Table & QC Figures

**Purpose:** Visualize data quality and filtering steps.

### Figure 01A: QC Filtering Flowchart

**File:** `fig_01A_QC_filtering_flowchart.png/pdf`

**Description:** Flowchart showing the number of genomes at each QC step.

**Justification:** Provides transparency about data filtering and sample size reduction at each step. Essential for reproducibility and understanding the final dataset composition.

**Data Flow:**
1. Raw genomes (from Script 01 input)
2. After CheckM completeness filter (≥90%)
3. After CheckM contamination filter (≤5%)
4. After gene count filter (`genes_total > 0`)
5. Final high-quality genomes

### Figure 01B: CheckM Completeness Distribution

**File:** `fig_01B_completeness_distribution.png/pdf`

**Description:** Histogram showing CheckM completeness scores for raw vs. high-quality genomes.

**X-axis:** CheckM completeness (0-100%)  
**Y-axis:** Number of genomes (frequency)

**Justification:** Visualizes the distribution of completeness scores and the effect of the ≥90% threshold. Helps identify potential biases in the filtered dataset.

### Figure 01C: CheckM Contamination Distribution

**File:** `fig_01C_contamination_distribution.png/pdf`

**Description:** Histogram showing CheckM contamination scores for raw vs. high-quality genomes.

**X-axis:** CheckM contamination (0-100%)  
**Y-axis:** Number of genomes (frequency)

**Justification:** Shows contamination levels and validates the ≤5% threshold. High contamination can indicate mixed cultures or assembly errors.

### Figure 01D: Genome Size Distribution

**File:** `fig_01D_genome_size_distribution.png/pdf`

**Description:** Histogram showing total gene count (`genes_total`) for raw vs. high-quality genomes.

**X-axis:** `genes_total` (number of genes)  
**Y-axis:** Number of genomes (frequency)

**Justification:** Displays the distribution of genome sizes (as gene counts) and ensures QC filtering doesn't introduce systematic biases in genome size.

### Figure 02A: Environment Counts (QC)

**File:** `fig_02A_environment_counts_QC.png/pdf`

**Description:** Bar plot showing number of genomes per environment before and after QC filtering.

**X-axis:** Environment name (GOLD metadata)  
**Y-axis:** Number of genomes

**Justification:** Identifies which environments are most represented and how QC filtering affects each environment. Helps understand the final environment distribution.

### Figure 02D: Genome Size by Environment

**File:** `fig_02D_genome_size_by_environment.png/pdf`

**Description:** Box plot or violin plot showing `genes_total` distribution stratified by environment.

**X-axis:** Environment name  
**Y-axis:** `genes_total` (number of genes)

**Justification:** Reveals environment-specific genome size patterns. If environments show distinct genome size distributions, this supports the hypothesis that environment constrains genome size.

---

## Script 02: Environment Cohort Figures

**Purpose:** Visualize environment filtering and final cohort selection.

### Figure 02B: Final Environments

**File:** `fig_02B_final_environments.png/pdf`

**Description:** Bar plot showing environments retained after filtering (≥20 genomes).

**X-axis:** Environment name  
**Y-axis:** Number of genomes

**Justification:** Shows which environments meet the minimum sample size requirement for reliable per-environment scaling fits. Environments with too few genomes would have unreliable parameter estimates.

### Figure 02C: Environment Contribution

**File:** `fig_02C_environment_contribution.png/pdf`

**Description:** Pie chart or bar plot showing the proportion of genomes contributed by each environment.

**Justification:** Highlights potential imbalances in environment representation. Highly imbalanced datasets may bias global scaling fits toward dominant environments.

### Figure 02D: Genome Size (Retained Environments)

**File:** `fig_02D_genome_size_retained_envs.png/pdf`

**Description:** Box plot showing `genes_total` distribution for retained environments only.

**X-axis:** Environment name  
**Y-axis:** `genes_total`

**Justification:** Final view of genome size distribution after all filtering. Used to identify environments with distinct genome size characteristics.

---

## Script 03: Global Scaling Figures

**Purpose:** Visualize global scaling laws fitted across all environments.

### Figure 04A: Global Scaling Exponent Distribution

**File:** `fig_04A_alpha_distribution.png/pdf`

**Description:** Histogram of global scaling exponents (\(\alpha_{\text{global}}\)) across all GO categories.

**X-axis:** \(\alpha_{\text{global}}\) (scaling exponent)  
**Y-axis:** Number of GO categories (frequency)

**Vertical lines:**
- Red dashed: Median \(\alpha\)
- Green dashed: \(\alpha = 1\) (linear scaling reference)

**Justification:** Reveals the distribution of scaling exponents across functional categories. Most categories should show \(\alpha \approx 1\) (linear scaling), with some showing sub-linear (\(\alpha < 1\)) or super-linear (\(\alpha > 1\)) patterns. This provides a baseline for understanding how different functional categories scale with genome size.

### Figure 04B: Exponent vs. R²

**File:** `fig_04B_alpha_vs_r2.png/pdf`

**Description:** Scatter plot of global scaling exponent vs. R² (coefficient of determination).

**X-axis:** \(\alpha_{\text{global}}\)  
**Y-axis:** R² (from OLS regression)

**Justification:** Identifies categories with poor fit quality (low R²) or unusual scaling exponents. Low R² indicates high variance around the power-law fit, suggesting the model may not capture the relationship well. Categories with high R² and extreme \(\alpha\) values are particularly interesting.

### Figure 04C: Exponent vs. Mean Count

**File:** `fig_04C_alpha_vs_mean_count.png/pdf`

**Description:** Scatter plot of global scaling exponent vs. mean GO category count per genome.

**X-axis:** \(\alpha_{\text{global}}\)  
**Y-axis:** Mean count (average \(n_c(g)\) across all genomes)

**Justification:** Tests whether scaling exponent correlates with category abundance. Rare categories (low mean count) may have less reliable exponent estimates due to sparse data. This helps identify potential artifacts from low-count categories.

### Figure 04D: Representative Global Scaling Plots

**File:** `fig_04D_representative_scaling.png/pdf` (log scale)  
**File:** `fig_04D_representative_scaling_linear.png/pdf` (linear scale)

**Description:** Faceted scatter plots showing log-log and linear-scale fits for top 20 GO categories (selected by Z-score variance or exponent magnitude).

**Each panel:**
- **X-axis (log):** \(\log(n(g))\) = \(\log(\text{genes\_total})\)  
- **Y-axis (log):** \(\log(n_c(g))\) = \(\log(\text{GO category count})\)  
- **X-axis (linear):** \(n(g)\) = `genes_total`  
- **Y-axis (linear):** \(n_c(g)\) = GO category count

**Data points:** Individual genomes (gray dots)  
**Red line:** Fitted power-law: \(y = \beta + \alpha \times x\) (log) or \(y = \exp(\beta) \times x^{\alpha}\) (linear)

**Justification:** Provides visual validation of the power-law model for representative categories. Log-log plots should show linear relationships if the power-law model holds. Linear-scale plots show the actual scaling in natural units, which is more interpretable for non-specialists. Selecting top categories by Z-score variance highlights categories with environment-specific variation.

**Metabolic Version:** `metabolic_fig_04D_representative_scaling.png/pdf`  
Filters to only metabolism-related GO terms (descendants of "metabolic process" GO:0008152).

---

## Script 04: Environment-Specific Scaling & Z-Score Figures

**Purpose:** Visualize environment-specific scaling and deviations from global patterns.

### Figure 05A: Z-Score Heatmap (Exponents)

**File:** `fig_05A_Z_alpha_heatmap.png/pdf`

**Description:** Heatmap showing \(Z_{\alpha}\) (Z-scores for exponents) for each environment × GO category combination.

**Rows:** GO categories (hierarchically clustered by Euclidean distance)  
**Columns:** Environments (hierarchically clustered)  
**Color scale:** \(Z_{\alpha}\) values (blue = negative, red = positive, white = near zero)

**Justification:** Identifies environment-category combinations with significant deviations from global scaling. Clustering reveals groups of categories or environments with similar scaling patterns. This is the primary visualization for detecting environment-specific scaling.

### Figure 05B: Z-Score Heatmap (Offsets)

**File:** `fig_05B_Z_beta_heatmap.png/pdf`

**Description:** Heatmap showing \(Z_{\beta}\) (Z-scores for offsets) for each environment × GO category combination.

**Rows:** GO categories (hierarchically clustered)  
**Columns:** Environments (hierarchically clustered)  
**Color scale:** \(Z_{\beta}\) values

**Justification:** Reveals environment-specific differences in baseline gene content (offset) independent of scaling exponent. Categories with high \(|Z_{\beta}|\) have environment-specific absolute abundances, even if scaling exponents are similar.

### Figure 05C: Absolute Z-Score by Environment

**File:** `fig_05C_abs_Z_alpha_by_env.png/pdf`

**Description:** Box plot showing distribution of \(|Z_{\alpha}|\) (absolute Z-scores) across GO categories, stratified by environment.

**X-axis:** Environment name  
**Y-axis:** \(|Z_{\alpha}|\) (absolute Z-score for exponents)

**Justification:** Identifies which environments show the most deviation from global scaling patterns. Environments with high median \(|Z_{\alpha}|\) have distinct scaling relationships compared to the global average.

### Figure 05D: Significant Categories by Environment

**File:** `fig_05D_significant_categories_by_env.png/pdf`

**Description:** Bar plot or heatmap showing number of GO categories with \(|Z_{\alpha}| > 2\) (significant deviation) per environment.

**X-axis:** Environment name  
**Y-axis:** Number of categories with \(|Z_{\alpha}| > 2\)

**Justification:** Quantifies how many functional categories show environment-specific scaling in each environment. Environments with many significant categories are candidates for detailed investigation.

### Figure 06A: Z-Score Distribution (Exponents)

**File:** `fig_06A_Z_alpha_category_histogram.png/pdf`

**Description:** Histogram of \(Z_{\alpha,\text{category}}\) (category-level Z-statistics for exponents) across all GO categories.

**X-axis:** \(Z_{\alpha,\text{category}}\) (root mean square Z-score across environments)  
**Y-axis:** Number of categories (frequency)

**Vertical line:** \(|Z| = 2\) threshold

**Justification:** Shows the distribution of environment-specific variation across categories. Categories with high \(Z_{\alpha,\text{category}}\) show consistent environment-specific scaling across multiple environments.

### Figure 06B: Z-Score Distribution (Offsets)

**File:** `fig_06B_Z_beta_category_histogram.png/pdf`

**Description:** Histogram of \(Z_{\beta,\text{category}}\) (category-level Z-statistics for offsets) across all GO categories.

**X-axis:** \(Z_{\beta,\text{category}}\)  
**Y-axis:** Number of categories (frequency)

**Justification:** Similar to Figure 06A but for offset parameters. Identifies categories with environment-specific baseline abundances.

### Figure 06C: Exception Categories

**File:** `fig_06C_exception_categories.png/pdf`

**Description:** Scatter plot or bar plot highlighting GO categories with extreme Z-scores or unusual scaling patterns.

**Justification:** Highlights categories that deviate most strongly from global patterns. These "exception" categories may represent functional groups under strong environment-specific selection.

### Figure 07: Environment-Stratified Scaling

**File:** `fig_07_env_stratified_scaling.png/pdf` (log scale)  
**File:** `fig_07_env_stratified_scaling_linear.png/pdf` (linear scale)

**Description:** Faceted scatter plots showing scaling relationships for top 20 GO categories (by Z-score variance), with points colored by environment and separate fit lines per environment.

**Each panel:**
- **X-axis:** \(\log(n(g))\) or \(n(g)\)  
- **Y-axis:** \(\log(n_c(g))\) or \(n_c(g)\)  
- **Points:** Colored by environment  
- **Lines:** Environment-specific fits (colored) and global fit (gray)

**Justification:** Visualizes environment-specific scaling for categories with high Z-score variance. If environment-specific lines diverge, this confirms environment-dependent scaling. This is the most direct visualization of the core hypothesis: that scaling relationships differ between environments.

**Metabolic Version:** `metabolic_fig_07_env_stratified_scaling.png/pdf`  
Filters to metabolism-related GO terms only.

---

## Script 06: Publication Figures (Main Results)

**Purpose:** Generate publication-quality figures following the structure of van Nimwegen (2003).

### Panel 1a: Z-Statistics for Exponents by Category

**File:** `fig1a_Z_exponents_by_category_env.png/pdf`

**Description:** Bar plot or scatter plot showing \(Z_{\alpha,\text{category}}\) for each GO category, sorted by Z-score magnitude.

**X-axis:** GO categories (sorted by \(Z_{\alpha,\text{category}}\))  
**Y-axis:** \(Z_{\alpha,\text{category}}\) (category-level Z-statistic)

**Horizontal line:** \(|Z| = 2\) (significance threshold)  
**Text labels:** Top categories labeled with GO term names

**Justification:** Provides a ranked view of categories by environment-specific variation in scaling exponents. Categories with high \(|Z_{\alpha,\text{category}}|\) are prime candidates for environment-specific scaling. This panel directly addresses the question: "Which functional categories show the strongest environment-specific scaling?"

**Metabolic Version:** `metabolic_fig1a_Z_exponents_by_category_env.png/pdf`  
Shows only metabolism-related categories.

### Panel 1b: Z-Statistics for Offsets by Category

**File:** `fig1b_Z_offsets_by_category_env.png/pdf`

**Description:** Bar plot showing \(Z_{\beta,\text{category}}\) for each GO category, sorted by Z-score magnitude.

**X-axis:** GO categories (sorted by \(Z_{\beta,\text{category}}\))  
**Y-axis:** \(Z_{\beta,\text{category}}\) (category-level Z-statistic for offsets)

**Horizontal line:** \(|Z| = 2\)  
**Text labels:** Top categories labeled

**Justification:** Similar to Panel 1a but for offset parameters. Identifies categories with environment-specific baseline abundances, independent of scaling exponent.

### Panels 1c-1e: Exponent Comparisons for Selected Categories

**File:** `fig1cde_env_exponents_selected_categories.png/pdf`

**Description:** Faceted plots showing fitted exponents (\(\alpha_{\text{env}}\)) with 99% confidence intervals for selected GO categories across all environments.

**Each panel:**
- **X-axis:** Environment name  
- **Y-axis:** \(\alpha_{\text{env}}\) (environment-specific exponent)  
- **Dots:** Point estimates  
- **Vertical bars:** 99% confidence intervals  
- **Horizontal dashed line:** \(\alpha_{\text{global}}\) (global exponent)

**Selected categories:** Low Z, medium Z, high Z, and additional interesting categories (up to 8-10 total)

**Justification:** Shows how scaling exponents vary across environments for representative categories. If confidence intervals don't overlap with the global exponent, this confirms environment-specific scaling. Categories are selected to represent the full range of Z-scores, providing both negative and positive examples.

### Panels 1f-1k: Scatter Plots with Fits

**File:** `fig1f_to_k_env_scatter_scaling.png/pdf` (log scale)  
**File:** `fig1f_to_k_env_scatter_scaling_linear.png/pdf` (linear scale)

**Description:** Faceted scatter plots showing scaling relationships for selected categories and environments, with both global and environment-specific fits.

**Each panel:**
- **X-axis:** \(\log(n(g))\) or \(n(g)\)  
- **Y-axis:** \(\log(n_c(g))\) or \(n_c(g)\)  
- **Gray points:** All genomes (all environments)  
- **Colored points:** Genomes from specific environment  
- **Colored line:** Environment-specific fit  
- **Gray line:** Global fit

**Selected combinations:** Top 5 categories × 3 representative environments = 15 panels

**Justification:** Provides direct visual evidence of environment-specific scaling. If the colored line (environment-specific fit) deviates from the gray line (global fit), this confirms environment-dependent scaling. The linear-scale version is more interpretable for non-specialists, while the log-scale version better shows the power-law relationship.

### Panels: GO Category Count vs. Total Annotated Domains

**File:** `fig1_domains_vs_total_domains.png/pdf` (log scale)  
**File:** `fig1_domains_vs_total_domains_linear.png/pdf` (linear scale)

**Description:** Faceted scatter plots showing relationship between GO category count and total annotated domains (sum of all GO counts per genome).

**Each panel:**
- **X-axis:** Total annotated domains = \(\sum_c n_c(g)\) (sum of all GO category counts)  
- **Y-axis:** \(n_c(g)\) (count for specific GO category)  
- **Points:** Colored by environment  
- **Lines:** Environment-specific fits

**Justification:** Tests whether category abundance scales with total functional diversity (total annotated domains) rather than just genome size. This addresses whether scaling is driven by genome size per se or by overall functional complexity. If scaling is similar to the `genes_total` plots, this suggests genome size is the primary driver. If scaling differs, this suggests functional diversity may be more relevant.

---

## Script 07: Extensions (TF & Mobile Elements)

**Purpose:** Extend scaling analyses to transcription factors and mobile elements.

### Transcription Factor Scaling

**Model:** Same power-law model applied to `tf_count`:
\[
n_{\text{TF}}(g) = \beta_{\text{TF}} \times n(g)^{\alpha_{\text{TF}}}
\]

**Expected:** \(\alpha_{\text{TF}} \approx 2\) (quadratic scaling), as transcription factors scale with the square of genome size (regulatory complexity hypothesis).

### Mobile Element Scaling

**Model:** Applied to `mobile_element_count`:
\[
n_{\text{mobile}}(g) = \beta_{\text{mobile}} \times n(g)^{\alpha_{\text{mobile}}}
\]

**Expected:** Variable scaling, as mobile elements may be more common in certain environments or genome size ranges.

**Output Files:**
- `tf_mobile_scaling_params.tsv`: Global and environment-specific scaling parameters
- `tf_mobile_env_Z_scores.tsv`: Z-scores comparing environment-specific to global parameters

**Justification:** Transcription factors and mobile elements represent distinct functional categories that may scale differently than GO-annotated protein domains. TF scaling tests the regulatory complexity hypothesis, while mobile element scaling may reveal environment-specific patterns of horizontal gene transfer.

---

## Prevalence Filtering

**Purpose:** Filter GO terms by prevalence (percentage of genomes with non-zero count).

**Method:** For a prevalence threshold \(P\) (e.g., 95%), keep only GO terms where:
\[
\frac{\text{Number of genomes with } n_c(g) > 0}{\text{Total number of genomes}} \geq \frac{P}{100}
\]

**Justification:** 
- **High prevalence (e.g., 95-99%):** Focuses on "ubiquitous" terms present in most genomes. This reduces noise from rare, environment-specific terms and provides a more conservative analysis.
- **Lower prevalence:** Includes more terms but may introduce sparsity issues in per-environment fits.

**Output Prefix:** Files generated with prevalence filtering are prefixed with `prev{P}_` (e.g., `prev95_global_scaling_params.tsv`).

---

## Metabolic Focus Analysis

**Purpose:** Generate metabolism-specific plots by filtering to GO terms related to metabolism.

**Method:** 
1. Parse GO ontology (go-basic.obo)
2. Identify metabolism-related root terms:
   - GO:0008152 (metabolic process)
   - GO:0044281 (small-molecule metabolic process)
   - GO:0006520 (amino acid metabolic process)
   - GO:0009056 (catabolic process)
   - GO:0009058 (biosynthetic process)
3. Traverse GO DAG to collect all descendant terms
4. Filter GO categories to only those in the metabolic sub-ontology

**Justification:** Metabolism is central to the hypothesis that environment constrains genome size through nutrient availability. Focusing on metabolic terms allows for a more targeted analysis of metabolic complexity vs. nutrient limitation hypotheses.

**Output Prefix:** Metabolic plots are prefixed with `metabolic_` (e.g., `metabolic_fig1a_Z_exponents_by_category_env.png`).

---

## Statistical Validation

### Regression Assumptions

1. **Linearity:** Log-log transformation assumes power-law relationship
2. **Independence:** Genomes are assumed independent (may be violated by phylogenetic structure)
3. **Homoscedasticity:** Variance of residuals should be constant (checked via R² and residual plots)
4. **Normality:** Residuals should be approximately normal (large sample sizes mitigate violations)

### Quality Control Checks

- **Minimum sample size:** ≥10 genomes per fit (per category or per environment×category)
- **Variance threshold:** Standard deviation of \(\log(n(g))\) must be ≥ 0.05 (ensures sufficient genome size variation)
- **R² threshold:** No hard threshold, but low R² (< 0.3) indicates poor fit quality
- **Z-score significance:** \(|Z| > 2\) used as approximate significance threshold

### Limitations

1. **Phylogenetic non-independence:** Genomes are not independent due to shared evolutionary history. This may inflate significance of Z-scores.
2. **Environment classification:** GOLD metadata may not capture all relevant environmental variation.
3. **GO annotation quality:** Incomplete or incorrect GO annotations may introduce noise.
4. **Power-law assumption:** Not all categories may follow power-law scaling; alternative models (e.g., exponential, logarithmic) are not tested.

---

## References

- van Nimwegen, E. (2003). Scaling laws in the functional content of genomes. *Trends in Genetics*, 19(9), 479-484.
- Gene Ontology Consortium. Gene Ontology annotations and ontology structure.
- GOLD (Genomes OnLine Database) for environment metadata.

---

## Figure File Naming Convention

- **Base files:** `fig{number}{letter}_{description}.{ext}`
- **With prevalence filter:** `prev{P}_fig{number}{letter}_{description}.{ext}`
- **Metabolic focus:** `metabolic_fig{number}{letter}_{description}.{ext}`
- **Combined:** `prev{P}_metabolic_fig{number}{letter}_{description}.{ext}`
- **Linear scale:** `{filename}_linear.{ext}`

**Example:** `prev95_metabolic_fig1a_Z_exponents_by_category_env_linear.png`

