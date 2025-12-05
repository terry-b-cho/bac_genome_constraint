#!/usr/bin/env python3
"""
Script 05: Map GO Labels

Map GO IDs (7-digit format) to human-readable names and definitions for plotting and reporting.
"""

import pandas as pd
from pathlib import Path
import sys
import argparse
import re
from prevalence_utils import get_prevalence_prefix

# Parse arguments
parser = argparse.ArgumentParser(description='Map GO IDs to names and definitions')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset (first 10 GO terms)')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for GO term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
# GO_TERMS_FILE will be set dynamically based on prevalence threshold
OBO_FILE = BASE_DIR / "data/go/go-basic.obo"
OUTPUT_DIR = BASE_DIR / "results/4_statistical_analyses/05_go_labels"
QC_LOG_FILE = OUTPUT_DIR / "qc_05_go_labels.log"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Initialize QC log
qc_log = []

def log_message(msg):
    """Add message to QC log and print to stdout."""
    qc_log.append(msg)
    print(msg)

def parse_obo_file(obo_file, term_ids):
    """
    Parse OBO file and extract term information.
    
    Returns dict mapping 7-digit GO ID to {go_id, name, namespace, definition}
    """
    term_data = {}
    current_term = None
    current_id = None
    
    with open(obo_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Start of a new term
            if line == "[Term]":
                current_term = {}
                current_id = None
                continue
            
            # End of term section
            if line == "" and current_term is not None and current_id is not None:
                # Convert GO:XXXXXXX to 7-digit format
                if current_id.startswith("GO:"):
                    go_id_7digit = current_id[3:]  # Remove "GO:" prefix
                    if go_id_7digit in term_ids:
                        term_data[go_id_7digit] = current_term
                current_term = None
                current_id = None
                continue
            
            # Parse term fields
            if current_term is not None:
                if line.startswith("id: "):
                    current_id = line[4:].strip()
                    current_term['go_id'] = current_id
                elif line.startswith("name: "):
                    current_term['name'] = line[6:].strip()
                elif line.startswith("namespace: "):
                    current_term['namespace'] = line[11:].strip()
                elif line.startswith("def: "):
                    # Definition is in quotes, extract text
                    def_text = line[5:].strip()
                    # Remove quotes and [source] references
                    def_text = re.sub(r'^"', '', def_text)
                    def_text = re.sub(r'"\s*\[.*$', '', def_text)
                    current_term['definition'] = def_text
    
    # Handle last term if file doesn't end with blank line
    if current_term is not None and current_id is not None:
        if current_id.startswith("GO:"):
            go_id_7digit = current_id[3:]
            if go_id_7digit in term_ids:
                term_data[go_id_7digit] = current_term
    
    return term_data

log_message("=" * 80)
log_message("Script 05: Map GO Labels")
if args.test_mode:
    log_message("  TEST MODE: Processing small subset only")
log_message("=" * 80)
log_message("")

# ============================================================================
# Load GO Term IDs
# ============================================================================

log_message("Loading GO term IDs...")
try:
    # Use prevalence-filtered terms if threshold is set, otherwise use all ubiquitous terms
    prefix = get_prevalence_prefix(args.prevalence_threshold)
    if prefix:
        # Try to load prevalence-filtered file
        filtered_file = BASE_DIR / f"results/3_GO_analyses/{prefix}ubiquitous_terms.txt"
        if filtered_file.exists():
            GO_TERMS_FILE = filtered_file
        else:
            # Fallback: compute on-the-fly from master table
            log_message(f"  ⚠ {filtered_file} not found, will compute prevalence filter from master table")
            GO_TERMS_FILE = BASE_DIR / "results/3_GO_analyses/ubiquitous_terms.txt"
    else:
        GO_TERMS_FILE = BASE_DIR / "results/3_GO_analyses/ubiquitous_terms.txt"
    
    with open(GO_TERMS_FILE, 'r') as f:
        go_term_ids = [line.strip() for line in f if line.strip()]
    
    if args.test_mode:
        go_term_ids = go_term_ids[:10]
        log_message(f"  ✓ Loaded {len(go_term_ids)} GO term IDs (TEST MODE: using first 10)")
    else:
        log_message(f"  ✓ Loaded {len(go_term_ids)} GO term IDs")
    
    expected_count = 334
    if len(go_term_ids) != expected_count and not args.test_mode:
        log_message(f"  ⚠ WARNING: Expected {expected_count} terms, found {len(go_term_ids)}")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load GO term IDs: {e}")
    sys.exit(1)

log_message("")

# ============================================================================
# Parse OBO File
# ============================================================================

log_message("Parsing OBO file...")
try:
    term_data = parse_obo_file(OBO_FILE, set(go_term_ids))
    log_message(f"  ✓ Parsed OBO file")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to parse OBO file: {e}")
    sys.exit(1)

log_message("")

# ============================================================================
# Map Terms and Validate
# ============================================================================

log_message("Mapping GO terms...")

go_labels = []
missing_terms = []

for term_id in go_term_ids:
    if term_id in term_data:
        data = term_data[term_id]
        go_labels.append({
            'category': term_id,
            'go_id': data.get('go_id', f'GO:{term_id}'),
            'name': data.get('name', 'Unknown'),
            'namespace': data.get('namespace', 'unknown'),
            'definition': data.get('definition', 'No definition available')
        })
    else:
        missing_terms.append(term_id)
        # Add fallback entry
        go_labels.append({
            'category': term_id,
            'go_id': f'GO:{term_id}',
            'name': f'unknown_{term_id}',
            'namespace': 'unknown',
            'definition': 'No definition available'
        })

go_labels_df = pd.DataFrame(go_labels)

log_message(f"  ✓ Mapped {len(go_labels_df)} GO terms")
if missing_terms:
    log_message(f"  ✗ WARNING: {len(missing_terms)} terms not found in OBO file:")
    for term in missing_terms[:10]:  # Show first 10
        log_message(f"    {term}")
    if len(missing_terms) > 10:
        log_message(f"    ... and {len(missing_terms) - 10} more")
else:
    log_message(f"  ✓ All {len(go_term_ids)} terms successfully mapped")

log_message("")

# ============================================================================
# Create Plot Labels
# ============================================================================

log_message("Creating plot labels...")

plot_labels = []
for _, row in go_labels_df.iterrows():
    # Create short label (trim to ~40 characters)
    name = row['name']
    if len(name) > 40:
        short_label = name[:37] + "..."
    else:
        short_label = name
    
    # Create label with GO ID: "name (GO:XXXXXXX)"
    go_id = row['go_id']
    full_label = f"{name} ({go_id})"
    if len(full_label) > 60:  # Truncate if too long
        full_label = name[:50] + f" ({go_id})"
    
    plot_labels.append({
        'category': row['category'],
        'short_label': short_label,
        'full_label': full_label,  # "name (GO:XXXXXXX)" format
        'super_category': ''  # Can be filled manually later
    })

plot_labels_df = pd.DataFrame(plot_labels)
log_message(f"  ✓ Created {len(plot_labels_df)} plot labels")

log_message("")

# ============================================================================
# Sample Mappings
# ============================================================================

log_message("Sample mappings (first 5-10):")
for i, (_, row) in enumerate(go_labels_df.head(10).iterrows()):
    log_message(f"  {row['category']} → {row['go_id']}: {row['name']}")

log_message("")

# ============================================================================
# Sanity Checks
# ============================================================================

log_message("Sanity Checks")
log_message("-" * 80)

# Count match
if len(go_labels_df) == len(go_term_ids):
    log_message(f"  ✓ Row count matches input: {len(go_labels_df)} rows")
else:
    log_message(f"  ✗ WARNING: Row count mismatch: {len(go_labels_df)} != {len(go_term_ids)}")

# No duplicates
duplicates = go_labels_df['category'].duplicated().sum()
if duplicates == 0:
    log_message(f"  ✓ No duplicate categories")
else:
    log_message(f"  ✗ ERROR: Found {duplicates} duplicate categories")

# All have definitions
no_def = (go_labels_df['definition'] == 'No definition available').sum()
if no_def == 0:
    log_message(f"  ✓ All terms have definitions")
else:
    log_message(f"  ⚠ {no_def} terms have fallback definitions")

log_message("")

# ============================================================================
# Outputs
# ============================================================================

log_message("Writing Outputs")
log_message("-" * 80)

# Write go_term_labels.tsv
prefix = get_prevalence_prefix(args.prevalence_threshold)
base_name_labels = "go_term_labels_test.tsv" if args.test_mode else "go_term_labels.tsv"
output_file_labels = f"{prefix}{base_name_labels}" if prefix else base_name_labels
go_labels_df.to_csv(OUTPUT_DIR / output_file_labels, sep='\t', index=False)
log_message(f"  ✓ {output_file_labels}: {len(go_labels_df)} terms")

# Write go_term_labels_for_plots.tsv
base_name_plots = "go_term_labels_for_plots_test.tsv" if args.test_mode else "go_term_labels_for_plots.tsv"
output_file_plots = f"{prefix}{base_name_plots}" if prefix else base_name_plots
plot_labels_df.to_csv(OUTPUT_DIR / output_file_plots, sep='\t', index=False)
log_message(f"  ✓ {output_file_plots}: {len(plot_labels_df)} terms")

# Write QC log
base_name_log = "qc_05_go_labels_test.log" if args.test_mode else "qc_05_go_labels.log"
qc_log_file = f"{prefix}{base_name_log}" if prefix else base_name_log
with open(OUTPUT_DIR / qc_log_file, 'w') as f:
    f.write('\n'.join(qc_log))
log_message(f"  ✓ {qc_log_file} written")

log_message("")
log_message("=" * 80)
log_message("Script 05 completed successfully!")
log_message("=" * 80)

