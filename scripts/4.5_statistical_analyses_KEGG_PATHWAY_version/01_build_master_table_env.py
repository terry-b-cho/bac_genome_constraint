#!/usr/bin/env python3
"""
Script 01: Build Master Table

Integrate KEGG term counts (reactions, KOs, pathways) with genome + environment metadata 
and apply basic QC filters.
"""

import pandas as pd
import numpy as np
import sys
import argparse
from pathlib import Path
from prevalence_utils import get_kegg_input_files, detect_kegg_data_type, get_kegg_columns, get_data_type_suffix

# Parse arguments
parser = argparse.ArgumentParser(description='Build master table with KEGG counts and genome metadata')
parser.add_argument('--test-mode', action='store_true', 
                    help='Run on small test subset (first 50 genomes)')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for KEGG term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
FEATURE_MATRIX_FILE = BASE_DIR / "results/2_JGIgold_KEGG_anayses_out/05_genome_feature_matrix.tsv"
OUTPUT_BASE_DIR = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/01_master_table"

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
    log_message(f"Script 01: Build Master Table - {data_type.upper()}")
    if args.test_mode:
        log_message("  TEST MODE: Processing small subset only")
    if args.prevalence_threshold is not None:
        log_message(f"  PREVALENCE FILTER: {args.prevalence_threshold}%% threshold")
    log_message("=" * 80)
    log_message("")
    
    # ============================================================================
    # A. Input Validation
    # ============================================================================
    
    log_message("A. Input Validation")
    log_message("-" * 80)
    
    # Get appropriate KEGG counts file
    try:
        KEGG_COUNTS_FILE = get_kegg_input_files(data_type, args.prevalence_threshold, BASE_DIR)
        log_message(f"Using input file: {KEGG_COUNTS_FILE}")
    except Exception as e:
        log_message(f"  ✗ ERROR: Failed to determine input file: {e}")
        sys.exit(1)
    
    # Load KEGG counts table
    log_message(f"Loading {data_type} counts table...")
    try:
        kegg_counts = pd.read_csv(KEGG_COUNTS_FILE, sep='\t', index_col=0)
        if args.test_mode:
            # Take first 50 genomes for testing
            kegg_counts = kegg_counts.head(50)
            log_message(f"  ✓ File loaded successfully (TEST MODE: using first 50 genomes)")
        else:
            log_message(f"  ✓ File loaded successfully")
    except Exception as e:
        log_message(f"  ✗ ERROR: Failed to load {data_type} counts file: {e}")
        continue  # Skip this data type
    
    # Verify KEGG counts structure
    log_message(f"  Number of genomes with {data_type} data: {len(kegg_counts)}")
    log_message(f"  Number of {data_type} categories: {len(kegg_counts.columns)}")
    
    # Get KEGG columns for this data type
    kegg_cols = get_kegg_columns(kegg_counts, data_type)
    log_message(f"  Identified {len(kegg_cols)} {data_type} columns")
    
    # Check that all KEGG columns are numeric
    non_numeric = []
    for col in kegg_cols:
        if not pd.api.types.is_numeric_dtype(kegg_counts[col]):
            non_numeric.append(col)
    
    if non_numeric:
        log_message(f"  ✗ WARNING: Non-numeric {data_type} columns found: {non_numeric[:10]}")
    else:
        log_message(f"  ✓ All {data_type} columns are numeric")
    
    # Verify all values are >= 0
    if (kegg_counts[kegg_cols] >= 0).all().all():
        log_message(f"  ✓ All {data_type} values are >= 0")
    else:
        negative_count = (kegg_counts[kegg_cols] < 0).sum().sum()
        log_message(f"  ✗ WARNING: Found {negative_count} negative values in {data_type} counts")
    
    # Check for unique genome identifiers
    if kegg_counts.index.is_unique:
        log_message(f"  ✓ Index contains unique genome identifiers")
    else:
        duplicates = kegg_counts.index.duplicated().sum()
        log_message(f"  ✗ WARNING: Found {duplicates} duplicate genome identifiers")
    
    log_message("")
    
    # Load feature matrix (only once, reuse for all data types)
    if data_type == DATA_TYPES[0]:  # Only log this once
        log_message("Loading feature matrix...")
        try:
            feature_matrix = pd.read_csv(
                FEATURE_MATRIX_FILE, 
                sep='\t',
                na_values=['']  # Treat empty strings as NaN
            )
            if args.test_mode:
                # Take first 50 accessions for testing
                feature_matrix = feature_matrix.head(50)
                log_message(f"  ✓ File loaded successfully (TEST MODE: using first 50 genomes)")
            else:
                log_message(f"  ✓ File loaded successfully")
        except Exception as e:
            log_message(f"  ✗ ERROR: Failed to load feature matrix file: {e}")
            sys.exit(1)
        
        log_message(f"  Number of genomes in feature matrix: {len(feature_matrix)}")
        
        # Verify required columns exist
        required_cols = ['accession', 'environment', 'genes_total', 'checkm_completeness', 
                         'checkm_contamination']
        missing_cols = [col for col in required_cols if col not in feature_matrix.columns]
        if missing_cols:
            log_message(f"  ✗ ERROR: Missing required columns: {missing_cols}")
            sys.exit(1)
        else:
            log_message(f"  ✓ All required columns present")
        
        # Check genes_total
        if pd.api.types.is_integer_dtype(feature_matrix['genes_total']):
            log_message(f"  ✓ genes_total is integer type")
        else:
            log_message(f"  ✗ WARNING: genes_total is not integer type")
        
        genes_positive = (feature_matrix['genes_total'] > 0).sum()
        log_message(f"  Number of genomes with genes_total > 0: {genes_positive} / {len(feature_matrix)}")
        
        # Log presence of empty strings and NaN in each column
        log_message("  Missing data summary:")
        for col in required_cols:
            nan_count = feature_matrix[col].isna().sum()
            empty_count = (feature_matrix[col] == '').sum() if feature_matrix[col].dtype == 'object' else 0
            if nan_count > 0 or empty_count > 0:
                log_message(f"    {col}: {nan_count} NaN, {empty_count} empty strings")
        
        log_message("")
    
    # ============================================================================
    # B. Join Integrity Checks
    # ============================================================================
    
    log_message("B. Join Integrity Checks")
    log_message("-" * 80)
    
    # Join on Genome (index) = accession
    log_message(f"Joining {data_type} counts with feature matrix...")
    log_message(f"  {data_type} counts index name: {kegg_counts.index.name}")
    log_message(f"  Feature matrix accession column: accession")
    
    # Reset index to make Genome a column for joining
    kegg_counts_reset = kegg_counts.reset_index()
    kegg_counts_reset.rename(columns={kegg_counts.index.name: 'Genome'}, inplace=True)
    
    # Perform inner join
    master_raw = pd.merge(
        feature_matrix,
        kegg_counts_reset,
        left_on='accession',
        right_on='Genome',
        how='inner'
    )
    
    expected_rows = 3088
    actual_rows = len(master_raw)
    difference = abs(actual_rows - expected_rows)
    percent_diff = (difference / expected_rows) * 100
    
    log_message(f"  Expected ≈{expected_rows} rows (all KEGG genomes)")
    log_message(f"  Actual: {actual_rows} rows")
    log_message(f"  Difference: {difference} genomes")
    
    if percent_diff > 2:
        log_message(f"  ✗ WARNING: Unexpected join size; investigate missing or duplicated keys.")
        log_message(f"    Difference ({percent_diff:.2f}%) exceeds 2% threshold")
    else:
        log_message(f"  ✓ Join size within expected range")
        if difference > 0:
            log_message(f"    Difference likely due to {difference} genomes with missing {data_type} annotations")
    
    # Check for duplicate accessions after join
    duplicate_accessions = master_raw['accession'].duplicated().sum()
    if duplicate_accessions > 0:
        log_message(f"  ✗ ERROR: Found {duplicate_accessions} duplicate accessions after join")
        continue  # Skip this data type
    else:
        log_message(f"  ✓ No duplicate accessions after join")
    
    # Check for missing environments
    missing_env = master_raw['environment'].isna().sum()
    empty_env = (master_raw['environment'] == '').sum()
    unclassified_env = (master_raw['environment'] == 'Unclassified').sum()
    log_message(f"  Missing environments: {missing_env} NaN, {empty_env} empty, {unclassified_env} 'Unclassified'")
    
    log_message("")
    
    # ============================================================================
    # C. QC Filtering
    # ============================================================================
    
    log_message("C. QC Filtering")
    log_message("-" * 80)
    
    # Start with all rows
    n_before_filter = len(master_raw)
    log_message(f"  Starting with {n_before_filter} genomes")
    
    # Apply filters
    # 1. Completeness: > 90, treat NaN as 0
    completeness_nan = master_raw['checkm_completeness'].isna().sum()
    master_raw['checkm_completeness'] = master_raw['checkm_completeness'].fillna(0)
    n_lost_completeness = (master_raw['checkm_completeness'] <= 90).sum()
    master_raw = master_raw[master_raw['checkm_completeness'] > 90].copy()
    log_message(f"  Completeness filter (>90):")
    log_message(f"    NaN values treated as 0: {completeness_nan}")
    log_message(f"    Lost: {n_lost_completeness} genomes")
    
    # 2. Contamination: < 5, treat NaN as 100
    contamination_nan = master_raw['checkm_contamination'].isna().sum()
    master_raw['checkm_contamination'] = master_raw['checkm_contamination'].fillna(100)
    n_lost_contamination = (master_raw['checkm_contamination'] >= 5).sum()
    master_raw = master_raw[master_raw['checkm_contamination'] < 5].copy()
    log_message(f"  Contamination filter (<5):")
    log_message(f"    NaN values treated as 100: {contamination_nan}")
    log_message(f"    Lost: {n_lost_contamination} genomes")
    
    # 3. Genes: genes_total > 0
    n_lost_genes = (master_raw['genes_total'] <= 0).sum()
    master_raw = master_raw[master_raw['genes_total'] > 0].copy()
    log_message(f"  Genes filter (>0): Lost {n_lost_genes} genomes")
    
    # 4. Environment: not null, not empty, ≠ "Unclassified"
    n_lost_env_null = master_raw['environment'].isna().sum()
    n_lost_env_empty = (master_raw['environment'] == '').sum()
    n_lost_env_unclassified = (master_raw['environment'] == 'Unclassified').sum()
    master_raw = master_raw[
        master_raw['environment'].notna() & 
        (master_raw['environment'] != '') & 
        (master_raw['environment'] != 'Unclassified')
    ].copy()
    n_lost_env = n_lost_env_null + n_lost_env_empty + n_lost_env_unclassified
    log_message(f"  Environment filter (not null/empty/'Unclassified'):")
    log_message(f"    Lost: {n_lost_env} genomes ({n_lost_env_null} null, {n_lost_env_empty} empty, {n_lost_env_unclassified} 'Unclassified')")
    
    n_after_filter = len(master_raw)
    log_message(f"  Total retained: {n_after_filter} genomes")
    log_message(f"  Total lost: {n_before_filter - n_after_filter} genomes")
    
    log_message("")
    
    # ============================================================================
    # D. Environment-Level Sanity Checks
    # ============================================================================
    
    log_message("D. Environment-Level Sanity Checks")
    log_message("-" * 80)
    
    # Compute count of genomes per environment
    env_counts = master_raw.groupby('environment').size().reset_index(name='n_genomes')
    env_counts = env_counts.sort_values('n_genomes', ascending=False)
    
    # Add env_n_genomes column to master table
    master_raw['env_n_genomes'] = master_raw.groupby('environment')['accession'].transform('count')
    
    log_message("  Environment counts (sorted by genome count):")
    for _, row in env_counts.iterrows():
        env_name = row['environment']
        n_genomes = row['n_genomes']
        log_message(f"    {env_name}: {n_genomes} genomes")
        
        # Check for warnings
        if n_genomes < 5:
            log_message(f"      ⚠ WARNING: Environment has < 5 genomes")
        if n_genomes > 1000:
            log_message(f"      ⚠ WARNING: Environment has > 1000 genomes (possible metadata mapping issue)")
    
    log_message("")
    
    # ============================================================================
    # E. Outputs
    # ============================================================================
    
    log_message("E. Writing Outputs")
    log_message("-" * 80)
    
    # Get suffix for this data type
    suffix = get_data_type_suffix(data_type)
    
    # Write master_table_raw.tsv (before QC filters, but after join)
    master_raw_before_qc = pd.merge(
        feature_matrix,
        kegg_counts_reset,
        left_on='accession',
        right_on='Genome',
        how='inner'
    )
    base_name_raw = "master_table_raw_test.tsv" if args.test_mode else "master_table_raw.tsv"
    output_file_raw = base_name_raw.replace('.tsv', f'{suffix}.tsv')
    master_raw_before_qc.to_csv(OUTPUT_BASE_DIR / output_file_raw, sep='\t', index=False)
    log_message(f"  ✓ {output_file_raw}: {len(master_raw_before_qc)} rows")
    
    # Write master_table_high_quality.tsv (after QC filters)
    base_name_hq = "master_table_high_quality_test.tsv" if args.test_mode else "master_table_high_quality.tsv"
    output_file_hq = base_name_hq.replace('.tsv', f'{suffix}.tsv')
    master_raw.to_csv(OUTPUT_BASE_DIR / output_file_hq, sep='\t', index=False)
    log_message(f"  ✓ {output_file_hq}: {len(master_raw)} rows")
    
    # Write environment_counts_all.tsv
    base_name_env = "environment_counts_all_test.tsv" if args.test_mode else "environment_counts_all.tsv"
    output_file_env = base_name_env.replace('.tsv', f'{suffix}.tsv')
    env_counts.to_csv(OUTPUT_BASE_DIR / output_file_env, sep='\t', index=False)
    log_message(f"  ✓ {output_file_env}: {len(env_counts)} environments")
    
    # Write QC log
    base_name_log = "qc_01_master_table_test.log" if args.test_mode else "qc_01_master_table.log"
    qc_log_file = base_name_log.replace('.log', f'{suffix}.log')
    with open(OUTPUT_BASE_DIR / qc_log_file, 'w') as f:
        f.write('\n'.join(qc_log))
    log_message(f"  ✓ {qc_log_file} written")
    
    log_message("")
    log_message("=" * 80)
    log_message(f"Script 01 completed successfully for {data_type}!")
    log_message("=" * 80)
    log_message("")

print("")
print("=" * 80)
print("Script 01: All data types processed!")
print("=" * 80)
