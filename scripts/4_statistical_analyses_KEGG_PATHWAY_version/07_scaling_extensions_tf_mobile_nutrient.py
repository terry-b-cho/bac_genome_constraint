#!/usr/bin/env python3
"""
Script 07: Scaling Extensions (TF, Mobile, Nutrient)

Extend scaling analyses to transcription factors (tf_count), mobile elements 
(mobile_element_count), and optionally nutrient-related GO categories.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import argparse
from prevalence_utils import get_prevalence_prefix

# Try to import scipy.stats, fall back to manual implementation
try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not available, using manual OLS implementation")

# Parse arguments
parser = argparse.ArgumentParser(description='Analyze scaling of TF and mobile elements')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for GO term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths (will be set dynamically based on prevalence threshold)
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
OUTPUT_DIR = BASE_DIR / "results/4_statistical_analyses/07_extensions"
QC_LOG_FILE = OUTPUT_DIR / "qc_07_extensions.log"

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
    
    numerator = np.sum((x - x_mean) * (y - y_mean))
    denominator = np.sum((x - x_mean) ** 2)
    
    if denominator == 0:
        return None
    
    slope = numerator / denominator
    intercept = y_mean - slope * x_mean
    
    y_pred = slope * x + intercept
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - y_mean) ** 2)
    r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0
    
    mse = ss_res / (n - 2) if n > 2 else 0
    sxx = denominator
    alpha_se = np.sqrt(mse / sxx) if sxx > 0 else 0
    beta_se = np.sqrt(mse * (1/n + x_mean**2 / sxx)) if sxx > 0 and n > 0 else 0
    
    if HAS_SCIPY and alpha_se > 0:
        t_stat = slope / alpha_se
        p_value = stats.t.sf(np.abs(t_stat), n - 2) * 2
    elif alpha_se > 0:
        t_stat = slope / alpha_se
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
log_message("Script 07: Scaling Extensions (TF, Mobile, Nutrient)")
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

# Verify required columns
required_cols = ['tf_count', 'mobile_element_count', 'genes_total', 'environment']
missing_cols = [col for col in required_cols if col not in df.columns]
if missing_cols:
    log_message(f"  ✗ ERROR: Missing required columns: {missing_cols}")
    sys.exit(1)

# Data validity checks
log_message("")
log_message("Data validation:")
tf_negative = (df['tf_count'] < 0).sum()
mobile_negative = (df['mobile_element_count'] < 0).sum()
log_message(f"  tf_count: min={df['tf_count'].min()}, max={df['tf_count'].max()}, negative={tf_negative}")
log_message(f"  mobile_element_count: min={df['mobile_element_count'].min()}, max={df['mobile_element_count'].max()}, negative={mobile_negative}")

if tf_negative > 0 or mobile_negative > 0:
    log_message(f"  ✗ WARNING: Found negative values in TF or mobile counts")

log_message("")

# ============================================================================
# Global Fits for TF and Mobile
# ============================================================================

log_message("=" * 80)
log_message("Global Scaling Fits")
log_message("=" * 80)
log_message("")

global_results = []
categories_to_fit = ['tf_count', 'mobile_element_count']

# 99% confidence interval
if HAS_SCIPY:
    t_critical_99 = stats.t.ppf(0.995, df=1000)
else:
    t_critical_99 = 2.576

for category in categories_to_fit:
    log_message(f"Fitting global scaling for {category}...")
    
    # Drop zeros
    subset = df[(df['genes_total'] > 0) & (df[category] > 0)].copy()
    
    if len(subset) < 10:
        log_message(f"  ✗ SKIPPING: Only {len(subset)} genomes with {category} > 0")
        continue
    
    # Compute x = log(n(g)) and y = log(count)
    x = np.log(subset['genes_total'].values)
    y = np.log(subset[category].values)
    
    # Check variance
    x_std = np.std(x)
    if x_std < 0.05:
        log_message(f"  ✗ SKIPPING: Insufficient variation (std={x_std:.4f})")
        continue
    
    # Fit OLS
    fit_result = fit_ols(x, y)
    if fit_result is None:
        log_message(f"  ✗ SKIPPING: Fit failed")
        continue
    
    # Calculate 99% CI
    alpha_ci99_low = fit_result['slope'] - t_critical_99 * fit_result['alpha_se']
    alpha_ci99_high = fit_result['slope'] + t_critical_99 * fit_result['alpha_se']
    
    global_results.append({
        'category': category,
        'environment': 'global',
        'alpha': fit_result['slope'],
        'alpha_se': fit_result['alpha_se'],
        'alpha_ci99_low': alpha_ci99_low,
        'alpha_ci99_high': alpha_ci99_high,
        'beta_log': fit_result['intercept'],
        'beta_log_se': fit_result['beta_se'],
        'n_genomes_used': fit_result['n'],
        'r_squared': fit_result['r_squared'],
        'p_value': fit_result['p_value']
    })
    
    log_message(f"  ✓ Fitted: α={fit_result['slope']:.4f} ± {fit_result['alpha_se']:.4f}, "
               f"R²={fit_result['r_squared']:.4f}, n={fit_result['n']}")

log_message("")

# ============================================================================
# Environment-Specific Fits
# ============================================================================

log_message("=" * 80)
log_message("Environment-Specific Scaling Fits")
log_message("=" * 80)
log_message("")

env_results = []
valid_envs = sorted(df['environment'].unique())

for category in categories_to_fit:
    log_message(f"Fitting {category} per environment...")
    
    for env in valid_envs:
        env_df = df[df['environment'] == env].copy()
        
        # Drop zeros
        subset = env_df[(env_df['genes_total'] > 0) & (env_df[category] > 0)].copy()
        
        if len(subset) < 10:
            continue
        
        x = np.log(subset['genes_total'].values)
        y = np.log(subset[category].values)
        
        x_std = np.std(x)
        if x_std < 0.05:
            continue
        
        fit_result = fit_ols(x, y)
        if fit_result is None:
            continue
        
        if (np.isnan(fit_result['slope']) or np.isinf(fit_result['slope']) or
            np.isnan(fit_result['intercept']) or np.isinf(fit_result['intercept'])):
            continue
        
        alpha_ci99_low = fit_result['slope'] - t_critical_99 * fit_result['alpha_se']
        alpha_ci99_high = fit_result['slope'] + t_critical_99 * fit_result['alpha_se']
        
        env_results.append({
            'category': category,
            'environment': env,
            'alpha': fit_result['slope'],
            'alpha_se': fit_result['alpha_se'],
            'alpha_ci99_low': alpha_ci99_low,
            'alpha_ci99_high': alpha_ci99_high,
            'beta_log': fit_result['intercept'],
            'beta_log_se': fit_result['beta_se'],
            'n_genomes_used': fit_result['n'],
            'r_squared': fit_result['r_squared'],
            'p_value': fit_result['p_value']
        })
    
    log_message(f"  ✓ {category}: {len([r for r in env_results if r['category'] == category])} environment fits")

log_message("")

# ============================================================================
# Z-Scores (Environment vs Global)
# ============================================================================

log_message("=" * 80)
log_message("Computing Z-Scores")
log_message("=" * 80)
log_message("")

z_scores = []
all_results = pd.DataFrame(global_results + env_results)

for _, row in all_results[all_results['environment'] != 'global'].iterrows():
    category = row['category']
    env = row['environment']
    
    # Get global parameters
    global_row = all_results[(all_results['category'] == category) & 
                             (all_results['environment'] == 'global')]
    if len(global_row) == 0:
        continue
    
    global_row = global_row.iloc[0]
    
    # Compute Z-scores
    z_alpha = np.nan
    z_beta = np.nan
    
    if row['alpha_se'] > 0 and global_row['alpha_se'] > 0:
        denominator = np.sqrt(row['alpha_se']**2 + global_row['alpha_se']**2)
        if denominator > 0:
            z_alpha = (row['alpha'] - global_row['alpha']) / denominator
    
    if row['beta_log_se'] > 0 and global_row['beta_log_se'] > 0:
        denominator = np.sqrt(row['beta_log_se']**2 + global_row['beta_log_se']**2)
        if denominator > 0:
            z_beta = (row['beta_log'] - global_row['beta_log']) / denominator
    
    if not (np.isnan(z_alpha) or np.isinf(z_alpha) or np.isnan(z_beta) or np.isinf(z_beta)):
        z_scores.append({
            'category': category,
            'environment': env,
            'Z_alpha': z_alpha,
            'Z_beta': z_beta,
            'n_genomes_used': row['n_genomes_used']
        })

z_scores_df = pd.DataFrame(z_scores)
log_message(f"  ✓ Computed {len(z_scores_df)} valid Z-scores")

log_message("")

# ============================================================================
# Outlier Detection (Cook's Distance)
# ============================================================================

log_message("=" * 80)
log_message("Outlier Detection (Cook's Distance)")
log_message("=" * 80)
log_message("")

# For global fits, compute Cook's distance
for category in categories_to_fit:
    subset = df[(df['genes_total'] > 0) & (df[category] > 0)].copy()
    
    if len(subset) < 10:
        continue
    
    x = np.log(subset['genes_total'].values)
    y = np.log(subset[category].values)
    
    # Get fit
    fit_result = fit_ols(x, y)
    if fit_result is None:
        continue
    
    # Compute residuals
    y_pred = fit_result['slope'] * x + fit_result['intercept']
    residuals = y - y_pred
    mse = np.sum(residuals**2) / (len(x) - 2)
    
    # Compute leverage (hat matrix diagonal)
    x_mean = np.mean(x)
    leverage = 1/len(x) + (x - x_mean)**2 / np.sum((x - x_mean)**2)
    
    # Cook's distance
    cooks_d = (residuals**2 / (2 * mse)) * (leverage / (1 - leverage)**2)
    
    # Flag high Cook's distance (> 4/n threshold)
    threshold = 4 / len(x)
    high_cooks = (cooks_d > threshold).sum()
    
    log_message(f"  {category}: {high_cooks} genomes with Cook's D > {threshold:.4f} "
               f"(max Cook's D = {cooks_d.max():.4f})")

log_message("")

# ============================================================================
# QC Summary
# ============================================================================

log_message("=" * 80)
log_message("QC Summary")
log_message("=" * 80)

global_df = pd.DataFrame(global_results)
if len(global_df) > 0:
    log_message("")
    log_message("Global scaling exponents:")
    for _, row in global_df.iterrows():
        log_message(f"  {row['category']}: α={row['alpha']:.4f} ± {row['alpha_se']:.4f}, "
                   f"R²={row['r_squared']:.4f}, n={row['n_genomes_used']}")
        
        # Check if TF exponent is ~2 (expected)
        if row['category'] == 'tf_count':
            if 1.5 < row['alpha'] < 2.5:
                log_message(f"    ✓ TF exponent ({row['alpha']:.4f}) is in expected range (~2)")
            else:
                log_message(f"    ⚠ TF exponent ({row['alpha']:.4f}) deviates from expected ~2")

log_message("")

# ============================================================================
# Outputs
# ============================================================================

log_message("=" * 80)
log_message("Writing Outputs")
log_message("=" * 80)

# Write tf_mobile_scaling_params.tsv
all_results_df = pd.DataFrame(global_results + env_results)
prefix = get_prevalence_prefix(args.prevalence_threshold)
base_name_params = "tf_mobile_scaling_params_test.tsv" if args.test_mode else "tf_mobile_scaling_params.tsv"
output_file_params = f"{prefix}{base_name_params}" if prefix else base_name_params
all_results_df.to_csv(OUTPUT_DIR / output_file_params, sep='\t', index=False)
log_message(f"  ✓ {output_file_params}: {len(all_results_df)} fits")

# Write tf_mobile_env_Z_scores.tsv
base_name_z = "tf_mobile_env_Z_scores_test.tsv" if args.test_mode else "tf_mobile_env_Z_scores.tsv"
output_file_z = f"{prefix}{base_name_z}" if prefix else base_name_z
z_scores_df.to_csv(OUTPUT_DIR / output_file_z, sep='\t', index=False)
log_message(f"  ✓ {output_file_z}: {len(z_scores_df)} Z-scores")

# Write QC log
base_name_log = "qc_07_extensions_test.log" if args.test_mode else "qc_07_extensions.log"
qc_log_file = f"{prefix}{base_name_log}" if prefix else base_name_log
with open(OUTPUT_DIR / qc_log_file, 'w') as f:
    f.write('\n'.join(qc_log))
log_message(f"  ✓ {qc_log_file} written")

log_message("")
log_message("=" * 80)
log_message("Script 07 completed successfully!")
log_message("=" * 80)

