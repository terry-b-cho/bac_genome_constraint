#!/usr/bin/env python3
"""
Script 02: Define Environment Cohorts

Select environments with sufficient sample size (≥20 genomes) and create an analysis-ready master table.
Processes reactions, KOs, and pathways.
"""

import pandas as pd
import sys
import argparse
from pathlib import Path
from prevalence_utils import get_prevalence_prefix, filter_kegg_columns_by_prevalence, get_kegg_columns, get_data_type_suffix

# Parse arguments
parser = argparse.ArgumentParser(description='Define environment cohorts with sufficient sample size')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for KEGG term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
INPUT_BASE_DIR = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/01_master_table"
OUTPUT_BASE_DIR = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts"

# Create output directory
OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

# Data types to process
DATA_TYPES = ['reactions', 'ko', 'pathway']

# ============================================================================
# Process each data type
# ============================================================================

for data_type in DATA_TYPES:
    # Initialize QC log for this data type
    qc_log = []
    
    def log_message(msg):
        """Add message to QC log and print to stdout."""
        qc_log.append(msg)
        print(msg)
    
    log_message("=" * 80)
    log_message(f"Script 02: Define Environment Cohorts - {data_type.upper()}")
    if args.test_mode:
        log_message("  TEST MODE: Processing small subset only")
    if args.prevalence_threshold is not None:
        log_message(f"  PREVALENCE FILTER: {args.prevalence_threshold}%% threshold")
    log_message("=" * 80)
    log_message("")
    
    # Get suffix for this data type
    suffix = get_data_type_suffix(data_type)
    
    # ============================================================================
    # A. Verify Consistency with Script 01
    # ============================================================================
    
    log_message("A. Verify Consistency with Script 01")
    log_message("-" * 80)
    
    # Construct input filename
    base_name_input = "master_table_high_quality_test.tsv" if args.test_mode else "master_table_high_quality.tsv"
    input_file_name = base_name_input.replace('.tsv', f'{suffix}.tsv')
    INPUT_FILE = INPUT_BASE_DIR / input_file_name
    
    # Load input file
    log_message(f"Loading {input_file_name}...")
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
        log_message(f"    Skipping {data_type}")
        continue
    
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
        base_name_min10 = "valid_environments_min10.tsv"
        output_file_min10 = base_name_min10.replace('.tsv', f'{suffix}.tsv')
        valid_environments_min10.to_csv(OUTPUT_BASE_DIR / output_file_min10, sep='\t', index=False)
        log_message(f"    Created {output_file_min10} with {len(valid_environments_min10)} environments")
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
    output_file_valid_name = base_name_valid.replace('.tsv', f'{suffix}.tsv')
    output_file_valid = f"{prefix}{output_file_valid_name}" if prefix else output_file_valid_name
    valid_environments.to_csv(OUTPUT_BASE_DIR / output_file_valid, sep='\t', index=False)
    log_message(f"  ✓ {output_file_valid}: {len(valid_environments)} environments")
    
    # Filter master table to only include valid environments
    master_filtered = master_table[master_table['environment'].isin(valid_environments['environment'])].copy()
    
    # Get KEGG columns for this data type
    kegg_cols_all = get_kegg_columns(master_filtered, data_type)
    
    # For pathways, explicitly drop "rn" columns if they exist (keep only "map")
    if data_type == 'pathway':
        rn_cols = [col for col in master_filtered.columns if col.startswith('rn') and col[2:].isdigit()]
        if len(rn_cols) > 0:
            master_filtered = master_filtered.drop(columns=rn_cols)
            log_message(f"  ✓ Removed {len(rn_cols)} 'rn' duplicate pathway columns (keeping only 'map' pathways)")
    
    # Refresh KEGG columns list after dropping "rn" columns
    kegg_cols_all = get_kegg_columns(master_filtered, data_type)
    
    # Filter KEGG columns by prevalence threshold if specified
    if args.prevalence_threshold is not None:
        log_message("")
        log_message(f"Filtering {data_type} columns by prevalence threshold...")
        kegg_cols_filtered = filter_kegg_columns_by_prevalence(master_filtered, args.prevalence_threshold, data_type)
        cols_to_drop = set(kegg_cols_all) - set(kegg_cols_filtered)
        if len(cols_to_drop) > 0:
            master_filtered = master_filtered.drop(columns=list(cols_to_drop))
            log_message(f"  ✓ Filtered {data_type} columns: {len(kegg_cols_all)} → {len(kegg_cols_filtered)} "
                       f"(removed {len(cols_to_drop)} columns below {args.prevalence_threshold}%% threshold)")
        else:
            log_message(f"  ✓ All {len(kegg_cols_all)} {data_type} columns meet {args.prevalence_threshold}%% threshold")
        kegg_cols_count = len(kegg_cols_filtered)
    else:
        kegg_cols_count = len(kegg_cols_all)
    
    # Write master_table_env_filtered.tsv
    base_name_tsv = "master_table_env_filtered_test.tsv" if args.test_mode else "master_table_env_filtered.tsv"
    output_file_tsv_name = base_name_tsv.replace('.tsv', f'{suffix}.tsv')
    output_file_tsv = f"{prefix}{output_file_tsv_name}" if prefix else output_file_tsv_name
    master_filtered.to_csv(OUTPUT_BASE_DIR / output_file_tsv, sep='\t', index=False)
    log_message(f"  ✓ {output_file_tsv}: {len(master_filtered)} genomes, {kegg_cols_count} {data_type} columns")
    
    # Write master_table_env_filtered.parquet
    base_name_parquet = "master_table_env_filtered_test.parquet" if args.test_mode else "master_table_env_filtered.parquet"
    output_file_parquet_name = base_name_parquet.replace('.parquet', f'{suffix}.parquet')
    output_file_parquet = f"{prefix}{output_file_parquet_name}" if prefix else output_file_parquet_name
    try:
        master_filtered.to_parquet(OUTPUT_BASE_DIR / output_file_parquet, engine='pyarrow')
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
    qc_log_file_name = base_name_log.replace('.log', f'{suffix}.log')
    qc_log_file = f"{prefix}{qc_log_file_name}" if prefix else qc_log_file_name
    with open(OUTPUT_BASE_DIR / qc_log_file, 'w') as f:
        f.write('\n'.join(qc_log))
    log_message(f"  ✓ {qc_log_file} written")
    
    log_message("")
    log_message("=" * 80)
    log_message(f"Script 02 completed successfully for {data_type}!")
    log_message("=" * 80)
    log_message("")

print("")
print("=" * 80)
print("Script 02: All data types processed!")
print("=" * 80)
