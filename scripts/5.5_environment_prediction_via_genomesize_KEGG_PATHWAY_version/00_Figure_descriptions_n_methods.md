# Figure Descriptions and Methods: KEGG Environment Prediction Pipeline

This document provides a comprehensive description of all figures, panels, mathematical models, methodological justifications, and results interpretation for the KEGG-based environment prediction pipeline (Script 5.5).

---

## Table of Contents

1. [Mathematical Framework](#mathematical-framework)
2. [Data Processing Pipeline](#data-processing-pipeline)
3. [Feature Engineering](#feature-engineering)
4. [Model Selection and Training](#model-selection-and-training)
5. [Evaluation Metrics](#evaluation-metrics)
6. [Figure Descriptions](#figure-descriptions)
7. [Results Summary](#results-summary)
8. [Robustness Checks and Limitations](#robustness-checks-and-limitations)

---

## Mathematical Framework

### Supervised Multi-Class Classification

The environment prediction task is formulated as a **multi-class classification problem**, where we predict one of 8 environment categories based on genomic features:

**Problem Formulation:**
\[
\hat{y} = f(\mathbf{x})
\]

Where:
- \(\mathbf{x} \in \mathbb{R}^{d}\) = feature vector (KEGG term counts + genome size)
- \(y \in \{1, 2, ..., 8\}\) = environment label (one of 8 categories)
- \(f: \mathbb{R}^{d} \to \{1, 2, ..., 8\}\) = classification function (learned model)
- \(\hat{y}\) = predicted environment

**Environment Categories:**
1. Aquatic
2. Terrestrial
3. Mammals: Human
4. Plants
5. Mammals
6. Food production
7. Wastewater
8. Birds

**KEGG Data Types:**
- **Reactions** (R-numbers): Biochemical reactions (e.g., R00130 = "ATP:dephospho-CoA 3'-phosphotransferase")
- **KOs** (K-numbers): Functional orthologs (e.g., K00005 = "glycerol dehydrogenase [EC:1.1.1.6]")
- **Pathways** (map-numbers): Metabolic pathways (e.g., map00010 = "Glycolysis / Gluconeogenesis")

### Model Architectures

#### 1. Random Forest (RF)

**Algorithm:** Ensemble of decision trees trained on bootstrapped samples with random feature subsets.

**Prediction:**
\[
\hat{y}_{\text{RF}} = \text{mode}\left(\{T_1(\mathbf{x}), T_2(\mathbf{x}), ..., T_B(\mathbf{x})\}\right)
\]

Where:
- \(T_i(\mathbf{x})\) = prediction from tree \(i\)
- \(B\) = number of trees (100 in our implementation)
- \(\text{mode}(\cdot)\) = majority vote

**Hyperparameters:**
- `n_estimators=100`: Number of trees
- `max_depth=20`: Maximum tree depth
- `min_samples_split=5`: Minimum samples to split a node
- `min_samples_leaf=2`: Minimum samples in a leaf
- `class_weight='balanced'`: Adjusts class weights to handle imbalance

**Feature Importance:**
\[
\text{Importance}(f_j) = \frac{1}{B} \sum_{i=1}^{B} \sum_{t \in T_i} \frac{n_t}{N} \Delta_t(j)
\]

Where:
- \(n_t\) = number of samples in node \(t\)
- \(N\) = total number of samples
- \(\Delta_t(j)\) = impurity decrease when feature \(j\) is used to split node \(t\)

#### 2. Gradient Boosting (GB)

**Algorithm:** Sequential ensemble where each tree corrects errors of previous trees.

**Prediction:**
\[
\hat{y}_{\text{GB}} = \text{argmax}_k \sum_{m=1}^{M} \alpha_m h_m(\mathbf{x}, k)
\]

Where:
- \(h_m(\mathbf{x}, k)\) = prediction of tree \(m\) for class \(k\)
- \(\alpha_m\) = learning rate (0.1)
- \(M\) = number of trees (100)

**Hyperparameters:**
- `n_estimators=100`: Number of boosting stages
- `learning_rate=0.1`: Shrinkage parameter
- `max_depth=5`: Maximum tree depth

#### 3. XGBoost (XGB)

**Algorithm:** Optimized gradient boosting with regularization and parallel processing.

**Objective Function:**
\[
\mathcal{L} = \sum_{i=1}^{n} l(y_i, \hat{y}_i) + \sum_{m=1}^{M} \Omega(f_m)
\]

Where:
- \(l(y_i, \hat{y}_i)\) = loss function (multi-class log-loss)
- \(\Omega(f_m)\) = regularization term (L1 + L2)

**Hyperparameters:**
- `n_estimators=100`
- `learning_rate=0.1`
- `max_depth=5`
- `tree_method='hist'` or `'gpu_hist'` (GPU acceleration if available)
- `objective='multi:softprob'`: Multi-class classification with probability outputs

#### 4. Logistic Regression (LR)

**Algorithm:** Linear model with multinomial logistic loss.

**Prediction Probability:**
\[
P(y=k|\mathbf{x}) = \frac{\exp(\mathbf{w}_k^T \mathbf{x} + b_k)}{\sum_{j=1}^{8} \exp(\mathbf{w}_j^T \mathbf{x} + b_j)}
\]

Where:
- \(\mathbf{w}_k\) = weight vector for class \(k\)
- \(b_k\) = bias term for class \(k\)

**Prediction:**
\[
\hat{y}_{\text{LR}} = \text{argmax}_k P(y=k|\mathbf{x})
\]

**Hyperparameters:**
- `max_iter=1000`: Maximum iterations
- `C=1.0`: Inverse regularization strength
- `class_weight='balanced'`: Handles class imbalance

**Feature Preprocessing:** StandardScaler (zero mean, unit variance):
\[
\mathbf{x}_{\text{scaled}} = \frac{\mathbf{x} - \mu}{\sigma}
\]

#### 5. Baseline (Dummy Classifier)

**Algorithm:** Predicts the majority class (most frequent environment).

**Prediction:**
\[
\hat{y}_{\text{baseline}} = \text{argmax}_k \sum_{i=1}^{n} \mathbb{1}[y_i = k]
\]

**Purpose:** Provides a baseline for comparison. Any model should outperform this.

---

## Data Processing Pipeline

### Input Data

**Source:** Master table from Script 02 (`master_table_env_filtered{suffix}.tsv`)

**Filters Applied:**
1. **Environment Filtering:** Only genomes with environment labels in the 8 valid categories
2. **KEGG Term Prevalence:** Optional filtering by prevalence threshold (e.g., 95% or 99%)
3. **Quality Control:** High-quality genomes (CheckM completeness ≥90%, contamination ≤5%)

**Final Dataset (varies by data type):**
- **Reactions:** ~485 KEGG reaction terms
- **KOs:** ~739 KEGG KO terms
- **Pathways:** ~101 KEGG pathway terms (map only, no rn duplicates)
- **Total genomes:** 2,164
- **Environments:** 8 categories
- **Features:** KEGG terms + `genes_total` + `{data_type}_total`

### Class Distribution

| Environment | Count | Percentage |
|------------|-------|------------|
| Aquatic | 598 | 27.6% |
| Terrestrial | 479 | 22.1% |
| Mammals: Human | 406 | 18.8% |
| Plants | 281 | 13.0% |
| Mammals | 198 | 9.1% |
| Food production | 103 | 4.8% |
| Wastewater | 63 | 2.9% |
| Birds | 36 | 1.7% |

**Note:** Class imbalance is addressed using `class_weight='balanced'` in tree-based models and Logistic Regression.

---

## Feature Engineering

### Base Features

1. **KEGG Term Counts:**
   - **Reactions:** Raw counts of each reaction per genome (R-numbers, e.g., R00130)
   - **KOs:** Raw counts of each KO per genome (K-numbers, e.g., K00005)
   - **Pathways:** Raw counts of each pathway per genome (map-numbers, e.g., map00010)
   - Optional prevalence filtering (e.g., 95% or 99% threshold)
   - Format: KEGG IDs (not 7-digit GO format)

2. **Genome Size (`genes_total`):**
   - Total number of genes in the genome
   - Direct measure of genome size

3. **Total Annotated KEGG Terms (`{data_type}_total`):**
   - Sum of all KEGG term counts per genome
   - \(\text{reactions\_total} = \sum_{c=1}^{n} n_c(g)\) for reactions
   - \(\text{ko\_total} = \sum_{c=1}^{n} n_c(g)\) for KOs
   - \(\text{pathway\_total} = \sum_{c=1}^{n} n_c(g)\) for pathways
   - Measures overall functional annotation density

### Feature Preprocessing

1. **Missing Value Handling:**
   - KEGG term counts: Filled with 0 (assumes unannotated = absent)
   - `genes_total`: Dropped if missing (critical feature)

2. **Zero Variance Filtering:**
   - Removed features with zero variance (constant across all genomes)

3. **Scaling (for Logistic Regression only):**
   - StandardScaler: \(\mathbf{x}_{\text{scaled}} = \frac{\mathbf{x} - \mu}{\sigma}\)
   - Applied only to training set; validation/test sets transformed using training parameters

### Feature Statistics

Feature statistics vary by data type and prevalence threshold. Typical ranges:
- **Mean:** Varies by data type
- **Standard Deviation:** Varies by data type
- **Range:** [0, varies] depending on KEGG term abundance

---

## Model Selection and Training

### Data Splitting

**Stratified Split (80/10/10):**
- **Training Set:** ~1,731 genomes (80%)
- **Validation Set:** ~216 genomes (10%)
- **Test Set:** ~217 genomes (10%)

**Stratification:** Ensures proportional representation of each environment class in all splits.

**Random Seed:** 42 (for reproducibility)

### Training Procedure

1. **Baseline:** Trained on full training set (majority class prediction)
2. **Tree-Based Models (RF, GB, XGB):**
   - Trained on raw (unscaled) features
   - No preprocessing required
3. **Logistic Regression:**
   - Features standardized using `StandardScaler`
   - Scaler fit on training set, applied to validation/test sets

### Model Hyperparameters

All hyperparameters were set to default/reasonable values (not optimized via cross-validation):

| Model | Key Hyperparameters |
|-------|---------------------|
| Random Forest | `n_estimators=100`, `max_depth=20`, `class_weight='balanced'` |
| Gradient Boosting | `n_estimators=100`, `learning_rate=0.1`, `max_depth=5` |
| XGBoost | `n_estimators=100`, `learning_rate=0.1`, `max_depth=5`, `tree_method='hist'` |
| Logistic Regression | `max_iter=1000`, `C=1.0`, `class_weight='balanced'` |

**Note:** Hyperparameter optimization (e.g., via grid search or Bayesian optimization) was not performed but could improve performance.

---

## Evaluation Metrics

### Accuracy Metrics

1. **Accuracy:**
   \[
   \text{Accuracy} = \frac{\text{Correct Predictions}}{\text{Total Predictions}} = \frac{1}{n} \sum_{i=1}^{n} \mathbb{1}[\hat{y}_i = y_i]
   \]

2. **Balanced Accuracy:**
   \[
   \text{Balanced Accuracy} = \frac{1}{K} \sum_{k=1}^{K} \frac{\text{TP}_k}{\text{TP}_k + \text{FN}_k}
   \]
   Where:
   - \(\text{TP}_k\) = True Positives for class \(k\)
   - \(\text{FN}_k\) = False Negatives for class \(k\)
   - \(K\) = number of classes (8)

   **Justification:** Accounts for class imbalance. Equal weight to each class, regardless of sample size.

### Per-Class Metrics

For each environment class \(k\):

1. **Precision:**
   \[
   \text{Precision}_k = \frac{\text{TP}_k}{\text{TP}_k + \text{FP}_k}
   \]

2. **Recall (Sensitivity):**
   \[
   \text{Recall}_k = \frac{\text{TP}_k}{\text{TP}_k + \text{FN}_k}
   \]

3. **F1-Score:**
   \[
   \text{F1}_k = \frac{2 \times \text{Precision}_k \times \text{Recall}_k}{\text{Precision}_k + \text{Recall}_k}
   \]

4. **Support:** Number of true instances of class \(k\) in test set

### Aggregate Metrics

1. **Macro F1-Score:**
   \[
   \text{Macro F1} = \frac{1}{K} \sum_{k=1}^{K} \text{F1}_k
   \]
   **Interpretation:** Average F1 across all classes (equal weight to each class)

2. **Weighted F1-Score:**
   \[
   \text{Weighted F1} = \sum_{k=1}^{K} w_k \times \text{F1}_k, \quad w_k = \frac{n_k}{n}
   \]
   **Interpretation:** F1 weighted by class frequency (favors majority classes)

### ROC Curves and AUC

**One-vs-Rest Approach:**
For each class \(k\), compute:
- **True Positive Rate (TPR):**
  \[
  \text{TPR}_k = \frac{\text{TP}_k}{\text{TP}_k + \text{FN}_k}
  \]
- **False Positive Rate (FPR):**
  \[
  \text{FPR}_k = \frac{\text{FP}_k}{\text{FP}_k + \text{TN}_k}
  \]

**AUC (Area Under ROC Curve):**
\[
\text{AUC}_k = \int_0^1 \text{TPR}_k(\text{FPR}_k^{-1}(t)) \, dt
\]

**Mean AUC:**
\[
\text{Mean AUC} = \frac{1}{K} \sum_{k=1}^{K} \text{AUC}_k
\]

### Overfitting Check

**Train-Validation Gap:**
\[
\text{Overfit} = \text{Accuracy}_{\text{train}} - \text{Accuracy}_{\text{val}}
\]

**Warning Threshold:** If overfit > 0.1, model may be overfitting.

---

## Figure Descriptions

### Figure 1: Model Performance Comparison

**File:** `prev99_fig01_model_performance_comparison{suffix}.png/pdf`

**Description:** Four-panel comparison of model performance metrics.

**Panels:**

1. **Test Accuracy by Model:**
   - Horizontal bar chart showing test accuracy for each model
   - Models sorted by accuracy (highest to lowest)
   - Values annotated on bars

2. **Test Balanced Accuracy by Model:**
   - Horizontal bar chart showing balanced accuracy
   - Accounts for class imbalance
   - More informative for imbalanced datasets

3. **Macro F1-Score by Model:**
   - Horizontal bar chart showing macro-averaged F1-score
   - Equal weight to each class

4. **Train vs Test Accuracy (Overfitting Check):**
   - Grouped bar chart comparing training and test accuracy
   - Large gap indicates overfitting

**Justification:** Provides a comprehensive overview of model performance across multiple metrics. Balanced accuracy and macro F1 are particularly important given class imbalance. The overfitting check helps identify models that may not generalize well.

---

### Figure 2: Confusion Matrices

**File:** `prev99_fig02_confusion_matrices{suffix}.png/pdf`

**Description:** Normalized confusion matrices for all models (baseline, RF, GB, XGB, LR).

**Format:**
- **Rows:** True environment labels
- **Columns:** Predicted environment labels
- **Values:** Normalized proportions (percentages)
- **Color scale:** Blue gradient (darker = higher proportion)

**Justification:** Reveals which environments are confused with each other. Diagonal elements show correct predictions; off-diagonal elements show common misclassifications.

**Common Confusions:**
- "Terrestrial" ↔ "Aquatic" (ecologically similar)
- "Mammals" ↔ "Mammals: Human" (taxonomically related)
- "Plants" ↔ "Terrestrial" (environmental overlap)

---

### Figure 3: ROC Curves (AUC)

**File:** `prev99_fig03_roc_curves{suffix}.png/pdf`

**Description:** ROC curves for each model, showing one-vs-rest classification performance for each environment class.

**Format:**
- **X-axis:** False Positive Rate (FPR)
- **Y-axis:** True Positive Rate (TPR)
- **Lines:** One curve per environment class (8 classes)
- **Diagonal dashed line:** Random classifier (AUC = 0.5)
- **Legend:** Class name with AUC score

**Justification:** ROC curves show the trade-off between sensitivity (TPR) and specificity (1-FPR) across different classification thresholds. AUC summarizes overall discriminative ability.

**Interpretation:**
- **AUC = 1.0:** Perfect classifier
- **AUC = 0.5:** Random classifier (no discriminative power)
- **AUC > 0.7:** Good discriminative ability
- **AUC > 0.9:** Excellent discriminative ability

---

### Figure 4: Per-Class Performance Metrics

**File:** `prev99_fig04_per_class_metrics{suffix}.png/pdf`

**Description:** Four-panel heatmap showing per-class metrics across all models.

**Panels:**

1. **Precision by Environment and Model:**
   - Heatmap: Rows = environments, Columns = models
   - Color scale: Yellow-Orange-Red (higher = better)
   - Values annotated in cells

2. **Recall by Environment and Model:**
   - Heatmap: Rows = environments, Columns = models
   - Color scale: Yellow-Green (higher = better)

3. **F1-Score by Environment and Model:**
   - Heatmap: Rows = environments, Columns = models
   - Color scale: Blue gradient (higher = better)

4. **Test Set Sample Size by Environment:**
   - Horizontal bar chart showing number of test samples per environment
   - Indicates class imbalance in test set

**Justification:** Provides detailed per-class performance breakdown. Helps identify which models perform best for specific environments and which environments are most challenging to predict.

---

### Figure 5: Feature Importance

**File:** `prev99_fig05_feature_importance{suffix}.png/pdf`

**Description:** Top 20 most important features for each model (RF, GB, XGB, LR), with annotated KEGG term labels.

**Format:**
- **Subplots:** One per model (4 panels)
- **Y-axis:** Top 20 features (sorted by importance, highest at top)
- **X-axis:** Feature importance score
- **Labels:** Human-readable KEGG term names (e.g., "glycerol dehydrogenase [EC:1.1.1.6]" for K00005) instead of numeric IDs
- **Special labels:**
  - "Total genes (genome size)" for `genes_total`
  - "Total annotated {data_type}" for `{data_type}_total` (e.g., "Total annotated reactions")

**Justification:** Identifies which genomic features are most predictive of environment. Provides biological insight into environment-specific functional requirements.

**Feature Importance Methods:**
- **Tree-based models (RF, GB, XGB):** Gini impurity decrease
- **Logistic Regression:** Mean absolute coefficient magnitude across all classes

**Biological Interpretation:**
- **Reactions:** Environment-specific biochemical reactions (e.g., anaerobic vs. aerobic pathways)
- **KOs:** Functional orthologs that vary by environment (e.g., transporters, metabolic enzymes)
- **Pathways:** Complete metabolic pathways that differ across environments
- **Genome size:** Larger genomes may indicate more complex environments or host-associated lifestyles

---

## Results Summary

### Overall Performance

Performance metrics vary by KEGG data type (reactions, KOs, pathways). Results will be generated separately for each data type.

**Key Observations:**
1. **All models should outperform baseline**, confirming predictive signal
2. **Tree-based models may show overfitting** (train accuracy >> test accuracy)
3. **Logistic Regression may have better balanced accuracy**, indicating better performance on minority classes
4. **Macro F1 reflects challenges with rare classes**

### Per-Environment Performance

**Best-Performing Environments (typically):**
- **Aquatic:** Usually highest accuracy (largest class, distinct features)
- **Mammals: Human:** Moderate accuracy (moderate class size)
- **Terrestrial:** Moderate accuracy (large class, but overlaps with Aquatic)

**Challenging Environments (typically):**
- **Wastewater:** Low recall (too rare)
- **Birds:** Low recall (too rare)
- **Food production:** Low recall (small class, overlaps with other environments)

---

## Robustness Checks and Limitations

### Robustness Checks

1. **Taxonomic Distribution:**
   - Each environment has diverse taxonomic representation
   - Total unique taxa: 2,164 (one per genome)
   - **Finding:** No single taxon dominates any environment (good for generalization)

2. **KEGG Annotation Completeness:**
   - Mean total KEGG counts per genome varies by environment and data type
   - **Potential Issue:** Annotation bias may confound predictions (some environments may have better annotations)

3. **Feature Statistics:**
   - Varies by data type and prevalence threshold
   - **Finding:** Wide range suggests good feature diversity

### Limitations

1. **Class Imbalance:**
   - Largest class (Aquatic): 27.6%
   - Smallest class (Birds): 1.7%
   - **Impact:** Models struggle with rare classes (Wastewater, Birds)
   - **Mitigation:** Used `class_weight='balanced'`, but limited training data for rare classes

2. **Overfitting:**
   - Tree-based models may show significant overfitting
   - **Causes:**
     - High model complexity (deep trees, many features)
     - Limited regularization
   - **Potential Solutions:**
     - Hyperparameter tuning (reduce `max_depth`, increase `min_samples_split`)
     - Cross-validation for model selection
     - Feature selection/dimensionality reduction

3. **Phylogenetic Confounding:**
   - Genomes are not independent (shared evolutionary history)
   - Certain taxa may be overrepresented in specific environments
   - **Impact:** Model may learn taxonomic signals rather than pure environmental signals
   - **Mitigation:** Could use phylogenetic-aware splitting or include taxonomic features explicitly

4. **Annotation Bias:**
   - KEGG annotation completeness varies by environment
   - Some environments may have better-curated annotations
   - **Impact:** Predictions may reflect annotation quality rather than true biological differences

5. **Environment Classification:**
   - GOLD metadata may not capture all relevant environmental variation
   - Some environments are broad categories (e.g., "Terrestrial" includes diverse habitats)
   - **Impact:** Reduces prediction accuracy and interpretability

6. **Feature Selection:**
   - Prevalence threshold (if used) ensures ubiquitous terms but may exclude environment-specific terms
   - **Trade-off:** Ubiquitous terms are more reliable but may be less informative for environment prediction

7. **No Hyperparameter Optimization:**
   - Hyperparameters set to default/reasonable values
   - **Impact:** Performance could potentially be improved with tuning

### Future Improvements

1. **Hyperparameter Tuning:**
   - Grid search or Bayesian optimization
   - Cross-validation for model selection

2. **Feature Engineering:**
   - Include lower-prevalence KEGG terms (e.g., 95% threshold)
   - Add taxonomic features (phylum, class, order)
   - Include additional genomic features (GC content, mobile elements, etc.)

3. **Class Balancing:**
   - Oversampling minority classes (SMOTE)
   - Undersampling majority classes
   - Cost-sensitive learning

4. **Ensemble Methods:**
   - Combine predictions from multiple models
   - Stacking or voting ensembles

5. **Phylogenetic Awareness:**
   - Phylogenetic cross-validation (split by clade, not randomly)
   - Include phylogenetic features or constraints

6. **Interpretability:**
   - SHAP values for feature attribution
   - Partial dependence plots
   - Model-agnostic explainability methods

---

## References

- Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.
- Friedman, J. H. (2001). Greedy function approximation: A gradient boosting machine. *Annals of Statistics*, 29(5), 1189-1232.
- Chen, T., & Guestrin, C. (2016). XGBoost: A scalable tree boosting system. *Proceedings of the 22nd ACM SIGKDD International Conference on Knowledge Discovery and Data Mining*, 785-794.
- Pedregosa, F., et al. (2011). Scikit-learn: Machine learning in Python. *Journal of Machine Learning Research*, 12, 2825-2830.
- KEGG (Kyoto Encyclopedia of Genes and Genomes) database.
- GOLD (Genomes OnLine Database) for environment metadata.

---

## File Naming Convention

All outputs use format: `{prefix}{description}{suffix}.{ext}`

**Examples:**
- `prev99_env_prediction_metrics_reactions.tsv`
- `prev99_fig01_model_performance_comparison_ko.png`
- `prev99_env_prediction_feature_importances_pathway.tsv`
- `prev99_qc_05_env_prediction_reactions.log`

Where:
- `prefix` = `prev99_` or `prev95_` or empty (based on prevalence threshold)
- `suffix` = `_reactions`, `_ko`, or `_pathway` (based on data type)
- `description` = descriptive name (e.g., `env_prediction_metrics`)

---

## Appendix: Command-Line Usage

### Basic Usage (Single Data Type)

```bash
python scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version.py \
    --data-type reactions \
    --prevalence-threshold 99 \
    --model all \
    --normalize none \
    --plot
```

### Process All Data Types (Bash Wrapper)

```bash
bash scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/run_prediction.sh 99 "" all
```

### SLURM Job Script

```bash
sbatch scripts/5.5_environment_prediction_via_genomesize_KEGG_PATHWAY_version/run_prediction.sh
```

### Key Arguments

- `--data-type`: KEGG data type (`reactions`, `ko`, `pathway`) - **REQUIRED**
- `--prevalence-threshold`: KEGG term prevalence threshold (0-100, e.g., 95 for 95%)
- `--test-mode`: Run on small subset for testing
- `--normalize`: Feature normalization (`none`, `per_gene`, `log`, `both`)
- `--model`: Which model(s) to train (`all`, `rf`, `gb`, `lr`, `xgb`, `baseline`)
- `--use-gpu`: Enable GPU acceleration for XGBoost (requires CUDA)
- `--output-dir`: Custom output directory (default: `results/5.5_environment_prediction/`)

### Data Type-Specific Notes

**Reactions:**
- R-number format (e.g., R00130)
- Typically ~485 reactions
- Biochemical reaction-level granularity

**KOs:**
- K-number format (e.g., K00005)
- Typically ~739 KOs
- Functional ortholog-level granularity

**Pathways:**
- map-number format (e.g., map00010)
- Only "map" pathways (no "rn" duplicates)
- Typically ~101 pathways
- Pathway-level granularity

---

**Document Version:** 1.0  
**Last Updated:** December 2025  
**Script Version:** Script 5.5 (KEGG Environment Prediction Pipeline)
