#!/usr/bin/env python3
"""
Script 02: Define Environment Cohorts

Select environments with sufficient sample size (≥20 genomes) and create an analysis-ready master table.
"""

import pandas as pd
import sys
import argparse
from pathlib import Path
from prevalence_utils import get_prevalence_prefix, filter_go_columns_by_prevalence

# Parse arguments
parser = argparse.ArgumentParser(description='Define environment cohorts with sufficient sample size')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for GO term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
INPUT_FILE = BASE_DIR / "results/4_statistical_analyses/01_master_table/master_table_high_quality.tsv"
OUTPUT_DIR = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts"
QC_LOG_FILE = OUTPUT_DIR / "qc_02_env_cohorts.log"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Initialize QC log
qc_log = []

def log_message(msg):
    """Add message to QC log and print to stdout."""
    qc_log.append(msg)
    print(msg)

log_message("=" * 80)
log_message("Script 02: Define Environment Cohorts")
if args.test_mode:
    log_message("  TEST MODE: Processing small subset only")
if args.prevalence_threshold is not None:
    log_message(f"  PREVALENCE FILTER: {args.prevalence_threshold}%% threshold")
log_message("=" * 80)
log_message("")

# ============================================================================
# A. Verify Consistency with Script 01
# ============================================================================

log_message("A. Verify Consistency with Script 01")
log_message("-" * 80)

# Load input file
log_message("Loading master_table_high_quality.tsv...")
try:
    master_table = pd.read_csv(INPUT_FILE, sep='\t')
    if args.test_mode:
        # For test mode, take a small subset but keep environments together
        master_table = (master_table.groupby('environment', group_keys=False)
                       .apply(lambda x: x.head(10) if len(x) >= 10 else x)
                       .reset_index(drop=True))
        log_message(f"  ✓ File loaded successfully (TEST MODE: using up to 10 genomes per environment)")
    else:
        log_message(f"  ✓ File loaded successfully")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load input file: {e}")
    sys.exit(1)

log_message(f"  Total genomes in input: {len(master_table)}")

# Check for missing environments
if master_table['environment'].isna().any():
    log_message(f"  ✗ WARNING: Found missing environments in input")
else:
    log_message(f"  ✓ No missing environments")

# Recompute environment counts as sanity check
log_message("  Recomputing environment counts...")
recomputed_counts = master_table.groupby('environment').size().reset_index(name='n_genomes')
recomputed_counts = recomputed_counts.sort_values('n_genomes', ascending=False)

# Check against env_n_genomes column
if 'env_n_genomes' in master_table.columns:
    env_n_genomes_from_col = master_table.groupby('environment')['env_n_genomes'].first().reset_index()
    env_n_genomes_from_col = env_n_genomes_from_col.rename(columns={'env_n_genomes': 'n_genomes_from_col'})
    env_n_genomes_from_col = env_n_genomes_from_col.sort_values('n_genomes_from_col', ascending=False)
    
    # Merge and compare
    comparison = pd.merge(
        recomputed_counts.rename(columns={'n_genomes': 'n_genomes_recomputed'}),
        env_n_genomes_from_col,
        on='environment'
    )
    
    mismatches = comparison[abs(comparison['n_genomes_recomputed'] - comparison['n_genomes_from_col']) > 1]
    
    if len(mismatches) > 0:
        log_message(f"  ✗ WARNING: Found {len(mismatches)} environments with count mismatches > 1:")
        for _, row in mismatches.iterrows():
            log_message(f"    {row['environment']}: recomputed={row['n_genomes_recomputed']}, from_col={row['n_genomes_from_col']}")
    else:
        log_message(f"  ✓ Environment counts match env_n_genomes column (within ±1 tolerance)")
else:
    log_message(f"  ⚠ NOTE: env_n_genomes column not found, skipping consistency check")

log_message("")

# ============================================================================
# B. Environment Threshold Filter
# ============================================================================

log_message("B. Environment Threshold Filter")
log_message("-" * 80)

# Keep environments with ≥ 20 genomes
threshold = 20
valid_environments = recomputed_counts[recomputed_counts['n_genomes'] >= threshold].copy()
dropped_environments = recomputed_counts[recomputed_counts['n_genomes'] < threshold].copy()

log_message(f"  Threshold: ≥ {threshold} genomes per environment")
log_message(f"  Environments kept: {len(valid_environments)}")
log_message(f"  Environments dropped: {len(dropped_environments)}")

log_message("")
log_message("  All environments with counts:")
for _, row in recomputed_counts.iterrows():
    env_name = row['environment']
    n_genomes = row['n_genomes']
    status = "✓ KEPT" if n_genomes >= threshold else "✗ DROPPED"
    log_message(f"    {status}: {env_name}: {n_genomes} genomes")

log_message("")
log_message("  Final retained environments:")
for _, row in valid_environments.iterrows():
    log_message(f"    {row['environment']}: {row['n_genomes']} genomes")

log_message("")

# ============================================================================
# C. Optional Flexibility
# ============================================================================

log_message("C. Optional Flexibility Check")
log_message("-" * 80)

if len(valid_environments) < 5:
    log_message(f"  ⚠ WARNING: Too few environments ({len(valid_environments)}) for robust scaling analysis.")
    log_message(f"    Consider lowering threshold to 15 or 10.")
    
    # Produce min10 variant
    valid_environments_min10 = recomputed_counts[recomputed_counts['n_genomes'] >= 10].copy()
    valid_environments_min10.to_csv(OUTPUT_DIR / "valid_environments_min10.tsv", sep='\t', index=False)
    log_message(f"    Created valid_environments_min10.tsv with {len(valid_environments_min10)} environments")
    log_message(f"    NOTE: This file exists but is not used by default")
else:
    log_message(f"  ✓ Sufficient environments ({len(valid_environments)}) for analysis")

log_message("")

# ============================================================================
# D. Outputs
# ============================================================================

log_message("D. Writing Outputs")
log_message("-" * 80)

# Write valid_environments_min20.tsv
prefix = get_prevalence_prefix(args.prevalence_threshold)
base_name_valid = "valid_environments_min20_test.tsv" if args.test_mode else "valid_environments_min20.tsv"
output_file_valid = f"{prefix}{base_name_valid}" if prefix else base_name_valid
valid_environments.to_csv(OUTPUT_DIR / output_file_valid, sep='\t', index=False)
log_message(f"  ✓ {output_file_valid}: {len(valid_environments)} environments")

# Filter master table to only include valid environments
master_filtered = master_table[master_table['environment'].isin(valid_environments['environment'])].copy()

# Filter GO columns by prevalence threshold if specified
if args.prevalence_threshold is not None:
    log_message("")
    log_message("Filtering GO columns by prevalence threshold...")
    go_cols_all = [col for col in master_filtered.columns 
                   if col.startswith('0') and len(col) == 7 and col.isdigit()]
    go_cols_filtered = filter_go_columns_by_prevalence(master_filtered, args.prevalence_threshold)
    cols_to_drop = set(go_cols_all) - set(go_cols_filtered)
    if len(cols_to_drop) > 0:
        master_filtered = master_filtered.drop(columns=list(cols_to_drop))
        log_message(f"  ✓ Filtered GO columns: {len(go_cols_all)} → {len(go_cols_filtered)} "
                   f"(removed {len(cols_to_drop)} columns below {args.prevalence_threshold}%% threshold)")
    else:
        log_message(f"  ✓ All {len(go_cols_all)} GO columns meet {args.prevalence_threshold}%% threshold")

# Write master_table_env_filtered.tsv
base_name_tsv = "master_table_env_filtered_test.tsv" if args.test_mode else "master_table_env_filtered.tsv"
output_file_tsv = f"{prefix}{base_name_tsv}" if prefix else base_name_tsv
master_filtered.to_csv(OUTPUT_DIR / output_file_tsv, sep='\t', index=False)
log_message(f"  ✓ {output_file_tsv}: {len(master_filtered)} genomes, {len(go_cols_filtered) if args.prevalence_threshold is not None else len(go_cols_all)} GO columns")

# Write master_table_env_filtered.parquet
base_name_parquet = "master_table_env_filtered_test.parquet" if args.test_mode else "master_table_env_filtered.parquet"
output_file_parquet = f"{prefix}{base_name_parquet}" if prefix else base_name_parquet
try:
    master_filtered.to_parquet(OUTPUT_DIR / output_file_parquet, engine='pyarrow')
    log_message(f"  ✓ {output_file_parquet}: {len(master_filtered)} genomes")
except Exception as e:
    log_message(f"  ✗ WARNING: Failed to write Parquet file: {e}")
    log_message(f"    Install pyarrow: conda install pyarrow")

# Sanity check: verify row count matches sum of n_genomes
sum_n_genomes = valid_environments['n_genomes'].sum()
if len(master_filtered) == sum_n_genomes:
    log_message(f"  ✓ Row count verification: {len(master_filtered)} genomes = sum of n_genomes")
else:
    log_message(f"  ✗ WARNING: Row count mismatch: {len(master_filtered)} != {sum_n_genomes}")

# Write QC log
base_name_log = "qc_02_env_cohorts_test.log" if args.test_mode else "qc_02_env_cohorts.log"
qc_log_file = f"{prefix}{base_name_log}" if prefix else base_name_log
with open(OUTPUT_DIR / qc_log_file, 'w') as f:
    f.write('\n'.join(qc_log))
log_message(f"  ✓ {qc_log_file} written")

log_message("")
log_message("=" * 80)
log_message("Script 02 completed successfully!")
log_message("=" * 80)

