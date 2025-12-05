#!/usr/bin/env python3
"""
Script 03: Global Scaling Fits

Fit global log-log OLS scaling laws across all environments for each GO category.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import argparse
from prevalence_utils import get_prevalence_prefix, filter_go_columns_by_prevalence

# Try to import scipy.stats, fall back to manual implementation
try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not available, using manual OLS implementation")

# Parse arguments
parser = argparse.ArgumentParser(description='Fit global scaling laws for all GO categories')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset (first 10 GO categories)')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for GO term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
INPUT_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet"
INPUT_FILE_TSV = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv"
GO_LABELS_FILE = BASE_DIR / "results/4_statistical_analyses/05_go_labels/go_term_labels.tsv"
OUTPUT_DIR = BASE_DIR / "results/4_statistical_analyses/03_global_scaling"
QC_LOG_FILE = OUTPUT_DIR / "qc_03_global_scaling.log"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Initialize QC log
qc_log = []

def log_message(msg):
    """Add message to QC log and print to stdout."""
    qc_log.append(msg)
    print(msg)

log_message("=" * 80)
log_message("Script 03: Global Scaling Fits")
if args.test_mode:
    log_message("  TEST MODE: Processing small subset only")
if args.prevalence_threshold is not None:
    log_message(f"  PREVALENCE FILTER: {args.prevalence_threshold}%% threshold")
log_message("=" * 80)
log_message("")

# ============================================================================
# Load Data
# ============================================================================

log_message("Loading data...")
try:
    prefix = get_prevalence_prefix(args.prevalence_threshold)
    base_name_parquet = "master_table_env_filtered.parquet"
    INPUT_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts" / (f"{prefix}{base_name_parquet}" if prefix else base_name_parquet)
    base_name_tsv = "master_table_env_filtered.tsv"
    INPUT_FILE_TSV = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts" / (f"{prefix}{base_name_tsv}" if prefix else base_name_tsv)
    
    df = pd.read_parquet(INPUT_FILE, engine='pyarrow')
    log_message(f"  ✓ Loaded {len(df)} genomes from {INPUT_FILE.name}")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load data: {e}")
    log_message(f"    Trying TSV format instead...")
    try:
        df = pd.read_csv(INPUT_FILE_TSV, sep='\t')
        log_message(f"  ✓ Loaded {len(df)} genomes from TSV ({INPUT_FILE_TSV.name})")
    except Exception as e2:
        log_message(f"  ✗ ERROR: Failed to load TSV: {e2}")
        sys.exit(1)

# Identify GO category columns (7-digit IDs) and filter by prevalence if needed
go_cols = filter_go_columns_by_prevalence(df, args.prevalence_threshold)
if args.test_mode:
    # Use first 10 GO categories for testing
    go_cols = go_cols[:10]
    log_message(f"  Found {len(go_cols)} GO category columns (TEST MODE: using first 10)")
else:
    if args.prevalence_threshold is not None:
        log_message(f"  Found {len(go_cols)} GO category columns (filtered by {args.prevalence_threshold}%% prevalence)")
    else:
        log_message(f"  Found {len(go_cols)} GO category columns (no prevalence filter)")
log_message("")

# Load GO labels if available
go_labels = None
if GO_LABELS_FILE.exists():
    try:
        go_labels = pd.read_csv(GO_LABELS_FILE, sep='\t')
        go_labels['category'] = go_labels['category'].astype(str).str.zfill(7)
        log_message(f"  ✓ Loaded GO labels for {len(go_labels)} categories")
    except Exception as e:
        log_message(f"  ⚠ Could not load GO labels: {e}")
        log_message("    Will use category IDs only")
else:
    log_message(f"  ⚠ GO labels file not found: {GO_LABELS_FILE}")
    log_message("    Run Script 05 first for better labels. Will use category IDs only")

def get_category_label(category_id):
    """Get human-readable label for a category ID."""
    if go_labels is not None:
        cat_str = str(category_id).zfill(7)
        label_row = go_labels[go_labels['category'] == cat_str]
        if len(label_row) > 0:
            name = label_row.iloc[0]['name']
            go_id = label_row.iloc[0]['go_id']
            return f"{name} ({go_id})"
    return f"Category {str(category_id).zfill(7)}"

log_message("")

# ============================================================================
# A. Pre-Regression Validation
# ============================================================================

log_message("A. Pre-Regression Validation")
log_message("-" * 80)

# Track statistics across categories
zero_drop_counts = []
n_genomes_per_category = []

for category in go_cols:
    # Drop rows with genes_total <= 0 or nc(g) <= 0
    subset = df[(df['genes_total'] > 0) & (df[category] > 0)].copy()
    
    n_dropped = len(df) - len(subset)
    zero_drop_counts.append(n_dropped)
    n_genomes_per_category.append(len(subset))

log_message(f"  Genomes removed due to zero counts:")
log_message(f"    Min: {min(zero_drop_counts)}")
log_message(f"    Median: {np.median(zero_drop_counts):.1f}")
log_message(f"    Max: {max(zero_drop_counts)}")

log_message(f"  Genomes remaining per category after filtering:")
log_message(f"    Min: {min(n_genomes_per_category)}")
log_message(f"    Median: {np.median(n_genomes_per_category):.1f}")
log_message(f"    Max: {max(n_genomes_per_category)}")

log_message("")

# ============================================================================
# B. Regression Validity Checks & Fitting
# ============================================================================

log_message("B. Fitting Global Scaling Laws")
log_message("-" * 80)

results = []
skipped_low_variance = []
skipped_low_n = []

        # 99% confidence interval: t_critical for large df (use df=1000 as approximation)
if HAS_SCIPY:
    t_critical_99 = stats.t.ppf(0.995, df=1000)  # 99% CI = 0.5% on each tail
else:
    # Approximate t-critical for 99% CI with large df (≈2.576 for normal distribution)
    t_critical_99 = 2.576

for category in go_cols:
    # Filter: genes_total > 0 and nc(g) > 0
    subset = df[(df['genes_total'] > 0) & (df[category] > 0)].copy()
    
    if len(subset) < 10:
        skipped_low_n.append(category)
        continue
    
    # Compute x = log(n(g)) and y = log(nc(g))
    x = np.log(subset['genes_total'].values)
    y = np.log(subset[category].values)
    
    # Check variance of x
    x_std = np.std(x)
    if x_std < 0.05:
        skipped_low_variance.append(category)
        log_message(f"  SKIPPING CATEGORY {category}: insufficient variation in genome size (std={x_std:.4f})")
        continue
    
    # Fit OLS regression
    try:
        if HAS_SCIPY:
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            r_squared = r_value ** 2
            alpha_se = std_err
        else:
            # Manual OLS implementation
            n = len(x)
            x_mean = np.mean(x)
            y_mean = np.mean(y)
            
            # Calculate slope and intercept
            numerator = np.sum((x - x_mean) * (y - y_mean))
            denominator = np.sum((x - x_mean) ** 2)
            
            if denominator == 0:
                log_message(f"  ✗ WARNING: Category {category} has zero variance in x, skipping")
                continue
            
            slope = numerator / denominator
            intercept = y_mean - slope * x_mean
            
            # Calculate residuals and R-squared
            y_pred = slope * x + intercept
            ss_res = np.sum((y - y_pred) ** 2)
            ss_tot = np.sum((y - y_mean) ** 2)
            r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
            
            # Calculate standard errors
            mse = ss_res / (n - 2) if n > 2 else 0
            sxx = denominator
            alpha_se = np.sqrt(mse / sxx) if sxx > 0 else 0
            
            # Calculate p-value (approximate using t-test)
            if alpha_se > 0:
                t_stat = slope / alpha_se
                # Use large-sample approximation for p-value
                p_value = 2 * (1 - 0.5 * (1 + np.sign(t_stat) * (1 - np.exp(-2 * t_stat**2 / np.pi))))
            else:
                p_value = np.nan
        
        # Check for NaN or inf
        if np.isnan(slope) or np.isinf(slope) or np.isnan(intercept) or np.isinf(intercept):
            log_message(f"  ✗ WARNING: Category {category} produced NaN/inf values, skipping")
            continue
        
        # Calculate standard error of intercept (beta)
        n = len(x)
        x_mean = np.mean(x)
        sxx = np.sum((x - x_mean) ** 2)
        mse = np.sum((y - (slope * x + intercept)) ** 2) / (n - 2) if n > 2 else 0
        beta_se = np.sqrt(mse * (1/n + x_mean**2 / sxx)) if sxx > 0 and n > 0 else 0
        
        # Calculate 99% confidence intervals
        alpha_ci99_low = slope - t_critical_99 * alpha_se
        alpha_ci99_high = slope + t_critical_99 * alpha_se
        
        beta_ci99_low = intercept - t_critical_99 * beta_se
        beta_ci99_high = intercept + t_critical_99 * beta_se
        
        # Degrees of freedom
        df_regression = n - 2
        
        # Store results
        results.append({
            'category': category,
            'alpha_global': slope,
            'alpha_global_se': alpha_se,
            'alpha_global_ci99_low': alpha_ci99_low,
            'alpha_global_ci99_high': alpha_ci99_high,
            'beta_global_log': intercept,
            'beta_global_log_se': beta_se,
            'beta_global_ci99_low': beta_ci99_low,
            'beta_global_ci99_high': beta_ci99_high,
            'n_genomes_used': n,
            'r_squared': r_squared,
            'p_value': p_value
        })
        
    except Exception as e:
        log_message(f"  ✗ ERROR: Failed to fit category {category}: {e}")
        continue

log_message(f"  Successfully fitted {len(results)} categories")
log_message(f"  Skipped {len(skipped_low_n)} categories due to low sample size (<10)")
log_message(f"  Skipped {len(skipped_low_variance)} categories due to low variance")

log_message("")

# ============================================================================
# C. Global Summary QC
# ============================================================================

log_message("C. Global Summary QC")
log_message("-" * 80)

if len(results) == 0:
    log_message("  ✗ ERROR: No successful fits!")
    sys.exit(1)

results_df = pd.DataFrame(results)

# Distribution summary for alpha
alpha_values = results_df['alpha_global'].values
log_message(f"  Distribution of alpha_global:")
log_message(f"    Mean: {np.mean(alpha_values):.4f}")
log_message(f"    Median: {np.median(alpha_values):.4f}")
log_message(f"    SD: {np.std(alpha_values):.4f}")
log_message(f"    Min: {np.min(alpha_values):.4f}")
log_message(f"    Max: {np.max(alpha_values):.4f}")

# QC Flags
alpha_negative = (alpha_values < 0).sum()
alpha_high = (alpha_values > 4).sum()
se_high = (results_df['alpha_global_se'] > 1).sum()
r2_low = (results_df['r_squared'] < 0.2).sum()

log_message("")
log_message(f"  QC Flags:")
log_message(f"    Categories with alpha < 0: {alpha_negative}")
log_message(f"    Categories with alpha > 4: {alpha_high}")
log_message(f"    Categories with SE > 1: {se_high}")
log_message(f"    Categories with R² < 0.2: {r2_low}")

# Top 10 categories with highest alpha
log_message("")
log_message(f"  Top 10 categories with highest alpha_global:")
top_10 = results_df.nlargest(10, 'alpha_global')[['category', 'alpha_global', 'r_squared']]
for _, row in top_10.iterrows():
    cat_label = get_category_label(row['category'])
    log_message(f"    {cat_label}: alpha={row['alpha_global']:.4f}, R²={row['r_squared']:.4f}")

log_message("")
log_message(f"  Categories skipped:")
log_message(f"    Low variation: {len(skipped_low_variance)}")
log_message(f"    Low sample size (<10): {len(skipped_low_n)}")

log_message("")

# ============================================================================
# D. Outputs
# ============================================================================

log_message("D. Writing Outputs")
log_message("-" * 80)

# Write global_scaling_params.tsv
prefix = get_prevalence_prefix(args.prevalence_threshold)
base_name = "global_scaling_params_test.tsv" if args.test_mode else "global_scaling_params.tsv"
output_file_tsv = f"{prefix}{base_name}" if prefix else base_name
results_df.to_csv(OUTPUT_DIR / output_file_tsv, sep='\t', index=False)
log_message(f"  ✓ {output_file_tsv}: {len(results_df)} categories")

# Write global_scaling_params.parquet (optional)
base_name_parquet = "global_scaling_params_test.parquet" if args.test_mode else "global_scaling_params.parquet"
output_file_parquet = f"{prefix}{base_name_parquet}" if prefix else base_name_parquet
try:
    results_df.to_parquet(OUTPUT_DIR / output_file_parquet, engine='pyarrow')
    log_message(f"  ✓ {output_file_parquet}: {len(results_df)} categories")
except Exception as e:
    log_message(f"  ⚠ NOTE: Could not write Parquet file: {e}")

# Write QC log
base_name_log = "qc_03_global_scaling_test.log" if args.test_mode else "qc_03_global_scaling.log"
qc_log_file = f"{prefix}{base_name_log}" if prefix else base_name_log
with open(OUTPUT_DIR / qc_log_file, 'w') as f:
    f.write('\n'.join(qc_log))
log_message(f"  ✓ {qc_log_file} written")

log_message("")
log_message("=" * 80)
log_message("Script 03 completed successfully!")
log_message("=" * 80)

