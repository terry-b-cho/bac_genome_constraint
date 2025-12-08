#!/usr/bin/env python3
"""
Script 05: Map KEGG Labels (Single Data Type)

Map KEGG IDs (reactions, KOs, pathways) to human-readable names for plotting and reporting.
"""

import pandas as pd
from pathlib import Path
import sys
import argparse
import time
import requests
from prevalence_utils import get_prevalence_prefix, get_kegg_columns, get_data_type_suffix

# Parse arguments
parser = argparse.ArgumentParser(description='Map KEGG IDs to names and definitions')
parser.add_argument('--data-type', type=str, required=True, choices=['reactions', 'ko', 'pathway'],
                    help='Data type to process: reactions, ko, or pathway')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset (first 10 terms)')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for KEGG term filtering (e.g., 95 for 95%%)')
parser.add_argument('--use-cache', action='store_true', default=True,
                    help='Use cached KEGG data if available')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
KEGG_CACHE_DIR = BASE_DIR / "results/3.5_KEGG_n_reaction_analyses/kegg_api_cache"
OUTPUT_BASE_DIR = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/05_kegg_labels"

# Create output directory
OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

# Initialize QC log
qc_log = []

def log_message(msg):
    """Add message to QC log and print to stdout."""
    qc_log.append(msg)
    print(msg, flush=True)

def fetch_kegg_info(kegg_id, data_type, use_cache=True):
    """
    Fetch KEGG information for a given ID.
    
    Args:
        kegg_id: KEGG identifier (e.g., R00130, K00001, map00010)
        data_type: 'reactions', 'ko', or 'pathway'
        use_cache: Whether to use cached data
    
    Returns:
        Dict with 'name' and 'definition' keys, or None if not found
    """
    # Check cache first
    if use_cache and KEGG_CACHE_DIR.exists():
        cache_file = KEGG_CACHE_DIR / f"{kegg_id}.txt"
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    content = f.read()
                    # Parse cached content
                    name = None
                    definition = None
                    for line in content.split('\n'):
                        if line.startswith('NAME'):
                            name = line.split(maxsplit=1)[1] if len(line.split(maxsplit=1)) > 1 else kegg_id
                        elif line.startswith('DEFINITION'):
                            definition = line.split(maxsplit=1)[1] if len(line.split(maxsplit=1)) > 1 else ""
                    return {'name': name or kegg_id, 'definition': definition or ""}
            except:
                pass
    
    # If not in cache, try KEGG API (with rate limiting)
    try:
        # Determine KEGG database
        if data_type == 'reactions':
            db = 'rn'
        elif data_type == 'ko':
            db = 'ko'
        elif data_type == 'pathway':
            # Pathways already have prefix (map/rn)
            db = None
        
        # Construct URL
        if db:
            url = f"http://rest.kegg.jp/get/{db}:{kegg_id}"
        else:
            url = f"http://rest.kegg.jp/get/{kegg_id}"
        
        response = requests.get(url, timeout=10)
        if response.status_code == 200:
            content = response.text
            
            # Save to cache
            if KEGG_CACHE_DIR.exists():
                cache_file = KEGG_CACHE_DIR / f"{kegg_id}.txt"
                with open(cache_file, 'w') as f:
                    f.write(content)
            
            # Parse response
            name = None
            definition = None
            for line in content.split('\n'):
                if line.startswith('NAME'):
                    name = line.split(maxsplit=1)[1] if len(line.split(maxsplit=1)) > 1 else kegg_id
                elif line.startswith('DEFINITION'):
                    definition = line.split(maxsplit=1)[1] if len(line.split(maxsplit=1)) > 1 else ""
            
            # Rate limiting
            time.sleep(0.1)
            
            return {'name': name or kegg_id, 'definition': definition or ""}
    except Exception as e:
        log_message(f"  ⚠ Could not fetch {kegg_id}: {e}")
    
    # Return default if all else fails
    return {'name': kegg_id, 'definition': ""}

data_type = args.data_type
suffix = get_data_type_suffix(data_type)

log_message("=" * 80)
log_message(f"Script 05: Map KEGG Labels - {data_type.upper()}")
if args.test_mode:
    log_message("  TEST MODE: Processing small subset only")
if args.prevalence_threshold is not None:
    log_message(f"  PREVALENCE FILTER: {args.prevalence_threshold}%% threshold")
log_message("=" * 80)
log_message("")

# ============================================================================
# Load KEGG Term IDs from master table
# ============================================================================

log_message("Loading KEGG term IDs from master table...")
try:
    prefix = get_prevalence_prefix(args.prevalence_threshold)
    base_name_tsv = "master_table_env_filtered.tsv"
    tsv_name = base_name_tsv.replace('.tsv', f'{suffix}.tsv')
    INPUT_FILE = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts" / (f"{prefix}{tsv_name}" if prefix else tsv_name)
    
    df = pd.read_csv(INPUT_FILE, sep='\t', nrows=10)  # Just need column names
    kegg_term_ids = get_kegg_columns(df, data_type)
    
    if args.test_mode:
        kegg_term_ids = kegg_term_ids[:10]
        log_message(f"  ✓ Loaded {len(kegg_term_ids)} {data_type} term IDs (TEST MODE: using first 10)")
    else:
        log_message(f"  ✓ Loaded {len(kegg_term_ids)} {data_type} term IDs")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load KEGG term IDs: {e}")
    sys.exit(1)

log_message("")

# ============================================================================
# Fetch KEGG Information
# ============================================================================

log_message("Fetching KEGG information...")
log_message(f"  Using cache: {args.use_cache}")
log_message("")

term_labels = []
not_found = []
fetch_count = 0

for term_id in kegg_term_ids:
    fetch_count += 1
    if fetch_count % 50 == 0:
        log_message(f"  Progress: {fetch_count}/{len(kegg_term_ids)} terms processed")
    
    info = fetch_kegg_info(term_id, data_type, use_cache=args.use_cache)
    
    if info:
        term_labels.append({
            'category': term_id,
            'name': info['name'],
            'definition': info['definition']
        })
    else:
        not_found.append(term_id)
        term_labels.append({
            'category': term_id,
            'name': term_id,
            'definition': ""
        })

log_message(f"  ✓ Processed {len(term_labels)} terms")
log_message(f"  Terms not found: {len(not_found)}")
if len(not_found) > 0 and len(not_found) <= 10:
    log_message(f"    {', '.join(not_found)}")

log_message("")

# ============================================================================
# Outputs
# ============================================================================

log_message("Writing Outputs")
log_message("-" * 80)

term_labels_df = pd.DataFrame(term_labels)

# Write kegg_term_labels.tsv
base_name = "kegg_term_labels_test.tsv" if args.test_mode else "kegg_term_labels.tsv"
output_file_name = base_name.replace('.tsv', f'{suffix}.tsv')
output_file = f"{prefix}{output_file_name}" if prefix else output_file_name
term_labels_df.to_csv(OUTPUT_BASE_DIR / output_file, sep='\t', index=False)
log_message(f"  ✓ {output_file}: {len(term_labels_df)} terms")

# Write QC log
base_name_log = "qc_05_kegg_labels_test.log" if args.test_mode else "qc_05_kegg_labels.log"
qc_log_file_name = base_name_log.replace('.log', f'{suffix}.log')
qc_log_file = f"{prefix}{qc_log_file_name}" if prefix else qc_log_file_name
with open(OUTPUT_BASE_DIR / qc_log_file, 'w') as f:
    f.write('\n'.join(qc_log))
log_message(f"  ✓ {qc_log_file} written")

log_message("")
log_message("=" * 80)
log_message(f"Script 05 completed successfully for {data_type}!")
log_message("=" * 80)

