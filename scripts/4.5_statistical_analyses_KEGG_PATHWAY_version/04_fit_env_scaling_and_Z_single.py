#!/usr/bin/env python3
"""
Script 04: Environment-Specific Fits & Z-Scores (Single Data Type)

This script processes ONE data type at a time.
Use 04_fit_env_scaling_and_Z.py wrapper to process all three.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import argparse
from prevalence_utils import get_prevalence_prefix, filter_kegg_columns_by_prevalence, get_data_type_suffix

# Try to import scipy.stats, fall back to manual implementation
try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("WARNING: scipy not available, using manual OLS implementation")

# Parse arguments
parser = argparse.ArgumentParser(description='Fit per-environment scaling laws and compute Z-scores')
parser.add_argument('--data-type', type=str, required=True, choices=['reactions', 'ko', 'pathway'],
                    help='Data type to process: reactions, ko, or pathway')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset (first 5 environments, first 10 categories)')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for KEGG term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
INPUT_BASE_DIR = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts"
GLOBAL_PARAMS_BASE_DIR = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/03_global_scaling"
OUTPUT_BASE_DIR = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling"

# Create output directory
OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

# Initialize QC log
qc_log = []

def log_message(msg):
    """Add message to QC log and print to stdout."""
    qc_log.append(msg)
    print(msg, flush=True)

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

data_type = args.data_type
suffix = get_data_type_suffix(data_type)

log_message("=" * 80)
log_message(f"Script 04: Environment-Specific Fits & Z-Scores - {data_type.upper()}")
if args.test_mode:
    log_message("  TEST MODE: Processing small subset only")
if args.prevalence_threshold is not None:
    log_message(f"  PREVALENCE FILTER: {args.prevalence_threshold}%% threshold")
log_message("=" * 80)
log_message("")

log_message("Loading data...")
prefix = get_prevalence_prefix(args.prevalence_threshold)

# Load master table
base_name_parquet = "master_table_env_filtered.parquet"
parquet_name = base_name_parquet.replace('.parquet', f'{suffix}.parquet')
INPUT_FILE = INPUT_BASE_DIR / (f"{prefix}{parquet_name}" if prefix else parquet_name)

base_name_tsv = "master_table_env_filtered.tsv"
tsv_name = base_name_tsv.replace('.tsv', f'{suffix}.tsv')
INPUT_FILE_TSV = INPUT_BASE_DIR / (f"{prefix}{tsv_name}" if prefix else tsv_name)

try:
    df = pd.read_parquet(INPUT_FILE, engine='pyarrow')
    log_message(f"  ✓ Loaded {len(df)} genomes from Parquet")
except Exception as e:
    log_message(f"  ⚠ Could not load Parquet: {e}")
    log_message(f"    Trying TSV format...")
    try:
        df = pd.read_csv(INPUT_FILE_TSV, sep='\t')
        log_message(f"  ✓ Loaded {len(df)} genomes from TSV")
    except Exception as e2:
        log_message(f"  ✗ ERROR: Failed to load data: {e2}")
        sys.exit(1)

# Load global parameters
base_name_global = "global_scaling_params.tsv"
global_params_name = base_name_global.replace('.tsv', f'{suffix}.tsv')
global_params_file = GLOBAL_PARAMS_BASE_DIR / (f"{prefix}{global_params_name}" if prefix else global_params_name)
try:
    global_params = pd.read_csv(global_params_file, sep='\t', dtype={'category': str})
    log_message(f"  ✓ Loaded {len(global_params)} global scaling parameters")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load global parameters: {e}")
    sys.exit(1)

# Load valid environments
base_name_valid = "valid_environments_min20.tsv"
valid_env_name = base_name_valid.replace('.tsv', f'{suffix}.tsv')
VALID_ENV_FILE = INPUT_BASE_DIR / (f"{prefix}{valid_env_name}" if prefix else valid_env_name)
try:
    valid_envs = pd.read_csv(VALID_ENV_FILE, sep='\t')
    log_message(f"  ✓ Loaded {len(valid_envs)} valid environments")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load valid environments: {e}")
    sys.exit(1)

kegg_cols = filter_kegg_columns_by_prevalence(df, args.prevalence_threshold, data_type)

if args.test_mode:
    test_envs = valid_envs['environment'].head(5).tolist()
    df = df[df['environment'].isin(test_envs)].copy()
    kegg_cols = kegg_cols[:10]
    log_message(f"  TEST MODE: Using {len(test_envs)} environments and {len(kegg_cols)} {data_type} categories")
else:
    test_envs = valid_envs['environment'].tolist()
    df = df[df['environment'].isin(test_envs)].copy()
    log_message(f"  Using {len(test_envs)} environments and {len(kegg_cols)} {data_type} categories")

log_message("")

log_message("A. Fitting Per-Environment Scaling Laws")
log_message("-" * 80)

env_scaling_results = []
skipped_low_n = []
skipped_low_variance = []

t_critical_99 = stats.t.ppf(0.995, df=1000) if HAS_SCIPY else 2.576

total_combos = len(test_envs) * len(kegg_cols)
fitted_combos = 0

for env in test_envs:
    env_df = df[df['environment'] == env].copy()
    env_fits = 0
    
    for category in kegg_cols:
        subset = env_df[(env_df['genes_total'] > 0) & (env_df[category] > 0)].copy()
        
        if len(subset) < 10:
            skipped_low_n.append((env, category))
            continue
        
        x = np.log(subset['genes_total'].values)
        y = np.log(subset[category].values)
        
        x_std = np.std(x)
        if x_std < 0.05:
            skipped_low_variance.append((env, category))
            continue
        
        fit_result = fit_ols(x, y)
        if fit_result is None:
            continue
        
        if (np.isnan(fit_result['slope']) or np.isinf(fit_result['slope']) or
            np.isnan(fit_result['intercept']) or np.isinf(fit_result['intercept'])):
            continue
        
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
    
    log_message(f"  {env}: {env_fits} / {len(kegg_cols)} categories fitted")

log_message("")
log_message(f"  Total env×category fits: {fitted_combos} / {total_combos}")
log_message(f"  Skipped due to low sample size (<10): {len(skipped_low_n)}")
log_message(f"  Skipped due to low variance: {len(skipped_low_variance)}")
log_message("")

log_message("B. Computing Z-Scores")
log_message("-" * 80)

env_scaling_df = pd.DataFrame(env_scaling_results)
z_scores = []
z_score_debug_count = 0

for _, row in env_scaling_df.iterrows():
    env = row['environment']
    category = row['category']
    category_str = str(category)
    
    global_row = global_params[global_params['category'] == category_str]
    if len(global_row) == 0:
        z_score_debug_count += 1
        if z_score_debug_count <= 3:
            log_message(f"  ⚠ No global params for category {category}")
        continue
    
    global_row = global_row.iloc[0]
    
    alpha_env = row['alpha_env']
    alpha_env_se = row['alpha_env_se']
    beta_env_log = row['beta_env_log']
    beta_env_log_se = row['beta_env_log_se']
    
    alpha_global = global_row['alpha_global']
    alpha_global_se = global_row['alpha_global_se']
    beta_global_log = global_row['beta_global_log']
    beta_global_log_se = global_row['beta_global_log_se']
    
    z_alpha = np.nan
    z_beta = np.nan
    
    if alpha_env_se > 0 and alpha_global_se > 0:
        denominator_alpha = np.sqrt(alpha_env_se**2 + alpha_global_se**2)
        if denominator_alpha > 0:
            z_alpha = (alpha_env - alpha_global) / denominator_alpha
    
    if beta_env_log_se > 0 and beta_global_log_se > 0:
        denominator_beta = np.sqrt(beta_env_log_se**2 + beta_global_log_se**2)
        if denominator_beta > 0:
            z_beta = (beta_env_log - beta_global_log) / denominator_beta
    
    if np.isnan(z_alpha) or np.isinf(z_alpha) or np.isnan(z_beta) or np.isinf(z_beta):
        z_score_debug_count += 1
        if z_score_debug_count <= 8:
            log_message(f"  ⚠ Skipping {env}×{category}: invalid Z")
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

log_message("")

log_message("C. Category-Level Aggregation")
log_message("-" * 80)

category_z_summary = []

if len(z_scores_df) > 0:
    for category in kegg_cols:
        cat_z_scores = z_scores_df[z_scores_df['category'] == category].copy()
    
        if len(cat_z_scores) < 2:
            continue
        
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

if len(category_z_summary_df) > 0:
    log_message(f"  Aggregated Z-scores for {len(category_z_summary_df)} categories")
    log_message(f"  Z_alpha: min={category_z_summary_df['Z_alpha_category'].min():.4f}, "
                f"median={category_z_summary_df['Z_alpha_category'].median():.4f}, "
                f"max={category_z_summary_df['Z_alpha_category'].max():.4f}")

log_message("")

log_message("D. Writing Outputs")
log_message("-" * 80)

base_name_env = "env_scaling_params_test.tsv" if args.test_mode else "env_scaling_params.tsv"
output_file_env_name = base_name_env.replace('.tsv', f'{suffix}.tsv')
output_file_env = f"{prefix}{output_file_env_name}" if prefix else output_file_env_name
env_scaling_df.to_csv(OUTPUT_BASE_DIR / output_file_env, sep='\t', index=False)
log_message(f"  ✓ {output_file_env}: {len(env_scaling_df)} env×category fits")

base_name_z = "env_vs_global_Z_scores_test.tsv" if args.test_mode else "env_vs_global_Z_scores.tsv"
output_file_z_name = base_name_z.replace('.tsv', f'{suffix}.tsv')
output_file_z = f"{prefix}{output_file_z_name}" if prefix else output_file_z_name
z_scores_df.to_csv(OUTPUT_BASE_DIR / output_file_z, sep='\t', index=False)
log_message(f"  ✓ {output_file_z}: {len(z_scores_df)} env×category Z-scores")

base_name_summary = "category_Z_summary_test.tsv" if args.test_mode else "category_Z_summary.tsv"
output_file_summary_name = base_name_summary.replace('.tsv', f'{suffix}.tsv')
output_file_summary = f"{prefix}{output_file_summary_name}" if prefix else output_file_summary_name
category_z_summary_df.to_csv(OUTPUT_BASE_DIR / output_file_summary, sep='\t', index=False)
log_message(f"  ✓ {output_file_summary}: {len(category_z_summary_df)} categories")

base_name_log = "qc_04_env_scaling_test.log" if args.test_mode else "qc_04_env_scaling.log"
qc_log_file_name = base_name_log.replace('.log', f'{suffix}.log')
qc_log_file = f"{prefix}{qc_log_file_name}" if prefix else qc_log_file_name
with open(OUTPUT_BASE_DIR / qc_log_file, 'w') as f:
    f.write('\n'.join(qc_log))
log_message(f"  ✓ {qc_log_file} written")

log_message("")
log_message("=" * 80)
log_message(f"Script 04 completed successfully for {data_type}!")
log_message("=" * 80)


