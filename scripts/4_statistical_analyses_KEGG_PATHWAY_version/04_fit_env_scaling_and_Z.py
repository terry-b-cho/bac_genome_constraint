#!/usr/bin/env python3
"""
Script 04: Environment-Specific Fits & Z-Scores

Fit per-environment scaling laws, then compute Z-scores comparing environment-specific 
parameters to global parameters.
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
parser = argparse.ArgumentParser(description='Fit per-environment scaling laws and compute Z-scores')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset (first 5 environments, first 10 GO categories)')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for GO term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
INPUT_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.parquet"
INPUT_FILE_TSV = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv"
# GLOBAL_PARAMS_FILE will be set dynamically based on prevalence threshold
VALID_ENV_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/valid_environments_min20.tsv"
GO_LABELS_FILE = BASE_DIR / "results/4_statistical_analyses/05_go_labels/go_term_labels.tsv"
OUTPUT_DIR = BASE_DIR / "results/4_statistical_analyses/04_env_scaling"
QC_LOG_FILE = OUTPUT_DIR / "qc_04_env_scaling.log"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Initialize QC log
qc_log = []

def log_message(msg):
    """Add message to QC log and print to stdout."""
    qc_log.append(msg)
    print(msg)

def fit_ols(x, y):
    """Fit OLS regression and return slope, intercept, R², SEs, and p-value."""
    n = len(x)
    if n < 2:
        return None
    
    x_mean = np.mean(x)
    y_mean = np.mean(y)
    
    # Calculate slope and intercept
    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sum((x - x_mean) ** 2)
    
    if denominator == 0:
        return None
    
    slope = numerator / denominator
    intercept = y_mean - slope * x_mean
    
    # Calculate residuals and R²
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y_mean) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    # Calculate standard errors
    mse = ss_res / (n - 2) if n > 2 else 0
    sxx = denominator
    alpha_se = np.sqrt(mse / sxx) if sxx > 0 else 0
    
    # Calculate standard error of intercept
    beta_se = np.sqrt(mse * (1/n + x_mean**2 / sxx)) if sxx > 0 and n > 0 else 0
    
    # Calculate p-value
    if HAS_SCIPY and alpha_se > 0:
        t_stat = slope / alpha_se
        p_value = stats.t.sf(np.abs(t_stat), n - 2) * 2
    elif alpha_se > 0:
        t_stat = slope / alpha_se
        # Approximate p-value
        p_value = 2 * (1 - 0.5 * (1 + np.sign(t_stat) * (1 - np.exp(-2 * t_stat**2 / np.pi))))
    else:
        p_value = np.nan
    
    return {
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_squared,
        'alpha_se': alpha_se,
        'beta_se': beta_se,
        'p_value': p_value,
        'n': n
    }

log_message("=" * 80)
log_message("Script 04: Environment-Specific Fits & Z-Scores")
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

# Load master table (with prevalence prefix if applicable)
prefix = get_prevalence_prefix(args.prevalence_threshold)
base_name_parquet = "master_table_env_filtered.parquet"
INPUT_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts" / (f"{prefix}{base_name_parquet}" if prefix else base_name_parquet)
base_name_tsv = "master_table_env_filtered.tsv"
INPUT_FILE_TSV = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts" / (f"{prefix}{base_name_tsv}" if prefix else base_name_tsv)

try:
    df = pd.read_parquet(INPUT_FILE, engine='pyarrow')
    log_message(f"  ✓ Loaded {len(df)} genomes from Parquet ({INPUT_FILE.name})")
except Exception as e:
    log_message(f"  ⚠ Could not load Parquet: {e}")
    log_message(f"    Trying TSV format...")
    try:
        df = pd.read_csv(INPUT_FILE_TSV, sep='\t')
        log_message(f"  ✓ Loaded {len(df)} genomes from TSV ({INPUT_FILE_TSV.name})")
    except Exception as e2:
        log_message(f"  ✗ ERROR: Failed to load data: {e2}")
        sys.exit(1)

# Load global parameters (with prevalence prefix if applicable)
prefix = get_prevalence_prefix(args.prevalence_threshold)
base_name_global = "global_scaling_params.tsv"
global_params_file = BASE_DIR / "results/4_statistical_analyses/03_global_scaling" / (f"{prefix}{base_name_global}" if prefix else base_name_global)
try:
    global_params = pd.read_csv(global_params_file, sep='\t', dtype={'category': str})
    # Convert category to string format with leading zeros if needed
    global_params['category'] = global_params['category'].astype(str).str.zfill(7)
    log_message(f"  ✓ Loaded {len(global_params)} global scaling parameters from {global_params_file.name}")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load global parameters from {global_params_file}: {e}")
    log_message(f"    Make sure Script 03 has been run with --prevalence-threshold {args.prevalence_threshold}")
    sys.exit(1)

# Load valid environments (with prevalence prefix if applicable)
base_name_valid = "valid_environments_min20.tsv"
VALID_ENV_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts" / (f"{prefix}{base_name_valid}" if prefix else base_name_valid)
try:
    valid_envs = pd.read_csv(VALID_ENV_FILE, sep='\t')
    log_message(f"  ✓ Loaded {len(valid_envs)} valid environments from {VALID_ENV_FILE.name}")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load valid environments from {VALID_ENV_FILE}: {e}")
    sys.exit(1)

# Load GO labels if available (with prevalence prefix if applicable)
go_labels = None
base_name_go_labels = "go_term_labels.tsv"
go_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}{base_name_go_labels}" if prefix else base_name_go_labels)
if go_labels_file.exists():
    try:
        go_labels = pd.read_csv(go_labels_file, sep='\t')
        go_labels['category'] = go_labels['category'].astype(str).str.zfill(7)
        log_message(f"  ✓ Loaded GO labels for {len(go_labels)} categories from {go_labels_file.name}")
    except Exception as e:
        log_message(f"  ⚠ Could not load GO labels: {e}")
else:
    log_message(f"  ⚠ GO labels file not found: {go_labels_file}")
    log_message("    Will use category IDs only")

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

# Identify GO category columns and filter by prevalence if needed
go_cols = filter_go_columns_by_prevalence(df, args.prevalence_threshold)

if args.test_mode:
    # Use first 5 environments and first 10 GO categories for testing
    test_envs = valid_envs['environment'].head(5).tolist()
    df = df[df['environment'].isin(test_envs)].copy()
    go_cols = go_cols[:10]
    log_message(f"  TEST MODE: Using {len(test_envs)} environments and {len(go_cols)} GO categories")
else:
    test_envs = valid_envs['environment'].tolist()
    df = df[df['environment'].isin(test_envs)].copy()
    if args.prevalence_threshold is not None:
        log_message(f"  Using {len(test_envs)} environments and {len(go_cols)} GO categories (filtered by {args.prevalence_threshold}%% prevalence)")
    else:
        log_message(f"  Using {len(test_envs)} environments and {len(go_cols)} GO categories")

log_message("")

# ============================================================================
# A. Fit QC per env×category
# ============================================================================

log_message("A. Fitting Per-Environment Scaling Laws")
log_message("-" * 80)

env_scaling_results = []
skipped_low_n = []
skipped_low_variance = []

# 99% confidence interval: t_critical
if HAS_SCIPY:
    t_critical_99 = stats.t.ppf(0.995, df=1000)
else:
    t_critical_99 = 2.576

total_combos = len(test_envs) * len(go_cols)
fitted_combos = 0

for env in test_envs:
    env_df = df[df['environment'] == env].copy()
    env_fits = 0
    
    for category in go_cols:
        # Drop genomes with genes_total <= 0 or nc(g) <= 0
        subset = env_df[(env_df['genes_total'] > 0) & (env_df[category] > 0)].copy()
        
        if len(subset) < 10:
            skipped_low_n.append((env, category))
            continue
        
        # Compute x = log(n(g)) and y = log(nc(g))
        x = np.log(subset['genes_total'].values)
        y = np.log(subset[category].values)
        
        # Check variance of x
        x_std = np.std(x)
        if x_std < 0.05:
            skipped_low_variance.append((env, category))
            continue
        
        # Fit OLS regression
        fit_result = fit_ols(x, y)
        if fit_result is None:
            continue
        
        # Check for NaN or inf
        if (np.isnan(fit_result['slope']) or np.isinf(fit_result['slope']) or
            np.isnan(fit_result['intercept']) or np.isinf(fit_result['intercept'])):
            continue
        
        # Calculate 99% confidence intervals
        alpha_ci99_low = fit_result['slope'] - t_critical_99 * fit_result['alpha_se']
        alpha_ci99_high = fit_result['slope'] + t_critical_99 * fit_result['alpha_se']
        
        env_scaling_results.append({
            'environment': env,
            'category': category,
            'alpha_env': fit_result['slope'],
            'alpha_env_se': fit_result['alpha_se'],
            'alpha_env_ci99_low': alpha_ci99_low,
            'alpha_env_ci99_high': alpha_ci99_high,
            'beta_env_log': fit_result['intercept'],
            'beta_env_log_se': fit_result['beta_se'],
            'n_genomes_used': fit_result['n'],
            'r_squared': fit_result['r_squared'],
            'p_value': fit_result['p_value']
        })
        
        env_fits += 1
        fitted_combos += 1
    
    log_message(f"  {env}: {env_fits} / {len(go_cols)} categories fitted")

log_message("")
log_message(f"  Total env×category fits: {fitted_combos} / {total_combos}")
log_message(f"  Skipped due to low sample size (<10): {len(skipped_low_n)}")
log_message(f"  Skipped due to low variance: {len(skipped_low_variance)}")

log_message("")

# ============================================================================
# B. Z-Score Formula Correctness
# ============================================================================

log_message("B. Computing Z-Scores")
log_message("-" * 80)

env_scaling_df = pd.DataFrame(env_scaling_results)
z_scores = []
z_score_debug_count = 0

for _, row in env_scaling_df.iterrows():
    env = row['environment']
    category = row['category']
    
    # Get global parameters (ensure category is string for comparison)
    category_str = str(category).zfill(7) if isinstance(category, (int, np.integer)) else str(category)
    global_row = global_params[global_params['category'] == category_str]
    if len(global_row) == 0:
        z_score_debug_count += 1
        if z_score_debug_count <= 3:
            log_message(f"  ⚠ No global params for category {category} (searched as {category_str})")
        continue
    
    global_row = global_row.iloc[0]
    
    # Get values
    alpha_env = row['alpha_env']
    alpha_env_se = row['alpha_env_se']
    beta_env_log = row['beta_env_log']
    beta_env_log_se = row['beta_env_log_se']
    
    alpha_global = global_row['alpha_global']
    alpha_global_se = global_row['alpha_global_se']
    beta_global_log = global_row['beta_global_log']
    beta_global_log_se = global_row['beta_global_log_se']
    
    # Check SEs > 0 before computing Z
    z_alpha = np.nan
    z_beta = np.nan
    
    if alpha_env_se > 0 and alpha_global_se > 0:
        denominator_alpha = np.sqrt(alpha_env_se**2 + alpha_global_se**2)
        if denominator_alpha > 0:
            z_alpha = (alpha_env - alpha_global) / denominator_alpha
    else:
        log_message(f"  ⚠ Z-score undefined for {env}×{category}: SE=0 (alpha)")
    
    if beta_env_log_se > 0 and beta_global_log_se > 0:
        denominator_beta = np.sqrt(beta_env_log_se**2 + beta_global_log_se**2)
        if denominator_beta > 0:
            z_beta = (beta_env_log - beta_global_log) / denominator_beta
    else:
        log_message(f"  ⚠ Z-score undefined for {env}×{category}: SE=0 (beta)")
    
    # Reject NaN or infinite Z-scores
    if np.isnan(z_alpha) or np.isinf(z_alpha) or np.isnan(z_beta) or np.isinf(z_beta):
        z_score_debug_count += 1
        # Only log first few to avoid spam
        if z_score_debug_count <= 5:
            log_message(f"  ⚠ Skipping {env}×{category}: invalid Z (alpha={z_alpha}, beta={z_beta}, alpha_se={alpha_env_se:.6f}, beta_se={beta_env_log_se:.6f})")
        continue
    
    z_scores.append({
        'environment': env,
        'category': category,
        'Z_alpha': z_alpha,
        'Z_beta': z_beta,
        'n_genomes_used': row['n_genomes_used']
    })

z_scores_df = pd.DataFrame(z_scores)
log_message(f"  Computed {len(z_scores_df)} valid Z-scores")

if len(z_scores_df) == 0:
    log_message("  ✗ WARNING: No valid Z-scores computed!")
    log_message("    This may indicate issues with SE values or data quality")

log_message("")

# ============================================================================
# C. Category-Level Aggregation QC
# ============================================================================

log_message("C. Category-Level Aggregation")
log_message("-" * 80)

category_z_summary = []

if len(z_scores_df) > 0:
    for category in go_cols:
        # Get all valid Z-scores for this category
        cat_z_scores = z_scores_df[z_scores_df['category'] == category].copy()
    
        if len(cat_z_scores) < 2:
            continue
        
        # Aggregate Z-scores
        z_alpha_squared = cat_z_scores['Z_alpha'] ** 2
        z_beta_squared = cat_z_scores['Z_beta'] ** 2
        
        z_alpha_category = np.sqrt(np.mean(z_alpha_squared))
        z_beta_category = np.sqrt(np.mean(z_beta_squared))
        
        n_envs_used = len(cat_z_scores)
        
        category_z_summary.append({
            'category': category,
            'Z_alpha_category': z_alpha_category,
            'Z_beta_category': z_beta_category,
            'n_envs_used': n_envs_used
        })

category_z_summary_df = pd.DataFrame(category_z_summary) if len(category_z_summary) > 0 else pd.DataFrame()

# Summary statistics
if len(category_z_summary_df) > 0:
    log_message(f"  Aggregated Z-scores for {len(category_z_summary_df)} categories")
    log_message(f"  Z_alpha_category: min={category_z_summary_df['Z_alpha_category'].min():.4f}, "
                f"median={category_z_summary_df['Z_alpha_category'].median():.4f}, "
                f"max={category_z_summary_df['Z_alpha_category'].max():.4f}")
    log_message(f"  Z_beta_category: min={category_z_summary_df['Z_beta_category'].min():.4f}, "
                f"median={category_z_summary_df['Z_beta_category'].median():.4f}, "
                f"max={category_z_summary_df['Z_beta_category'].max():.4f}")
    
    # Warnings
    low_n_envs = (category_z_summary_df['n_envs_used'] < 3).sum()
    high_z_alpha = (category_z_summary_df['Z_alpha_category'] > 2).sum()
    high_z_beta = (category_z_summary_df['Z_beta_category'] > 2).sum()
    
    if low_n_envs > 0:
        log_message(f"  ⚠ {low_n_envs} categories with n_envs_used < 3")
    if high_z_alpha > 0:
        log_message(f"  ⚠ {high_z_alpha} categories with Z_alpha_category > 2")
    if high_z_beta > 0:
        log_message(f"  ⚠ {high_z_beta} categories with Z_beta_category > 2")
    
    # Top categories
    log_message("")
    log_message(f"  Top 10 categories by Z_alpha_category:")
    top_z_alpha = category_z_summary_df.nlargest(10, 'Z_alpha_category')
    for _, row in top_z_alpha.iterrows():
        cat_label = get_category_label(row['category'])
        log_message(f"    {cat_label}: Z_alpha={row['Z_alpha_category']:.4f}, n_envs={row['n_envs_used']}")
    
    log_message("")
    log_message(f"  Top 10 categories by Z_beta_category:")
    top_z_beta = category_z_summary_df.nlargest(10, 'Z_beta_category')
    for _, row in top_z_beta.iterrows():
        cat_label = get_category_label(row['category'])
        log_message(f"    {cat_label}: Z_beta={row['Z_beta_category']:.4f}, n_envs={row['n_envs_used']}")

log_message("")

# ============================================================================
# D. Outputs
# ============================================================================

log_message("D. Writing Outputs")
log_message("-" * 80)

# Write env_scaling_params.tsv
prefix = get_prevalence_prefix(args.prevalence_threshold)
base_name_env = "env_scaling_params_test.tsv" if args.test_mode else "env_scaling_params.tsv"
output_file_env = f"{prefix}{base_name_env}" if prefix else base_name_env
env_scaling_df.to_csv(OUTPUT_DIR / output_file_env, sep='\t', index=False)
log_message(f"  ✓ {output_file_env}: {len(env_scaling_df)} env×category fits")

# Write env_vs_global_Z_scores.tsv
base_name_z = "env_vs_global_Z_scores_test.tsv" if args.test_mode else "env_vs_global_Z_scores.tsv"
output_file_z = f"{prefix}{base_name_z}" if prefix else base_name_z
z_scores_df.to_csv(OUTPUT_DIR / output_file_z, sep='\t', index=False)
log_message(f"  ✓ {output_file_z}: {len(z_scores_df)} env×category Z-scores")

# Write category_Z_summary.tsv
base_name_summary = "category_Z_summary_test.tsv" if args.test_mode else "category_Z_summary.tsv"
output_file_summary = f"{prefix}{base_name_summary}" if prefix else base_name_summary
category_z_summary_df.to_csv(OUTPUT_DIR / output_file_summary, sep='\t', index=False)
log_message(f"  ✓ {output_file_summary}: {len(category_z_summary_df)} categories")

# Write QC log
base_name_log = "qc_04_env_scaling_test.log" if args.test_mode else "qc_04_env_scaling.log"
qc_log_file = f"{prefix}{base_name_log}" if prefix else base_name_log
with open(OUTPUT_DIR / qc_log_file, 'w') as f:
    f.write('\n'.join(qc_log))
log_message(f"  ✓ {qc_log_file} written")

log_message("")
log_message("=" * 80)
log_message("Script 04 completed successfully!")
log_message("=" * 80)

