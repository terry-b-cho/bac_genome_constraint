#!/usr/bin/env python3
"""
Create prevalence-filtered GO terms file.

This script filters ubiquitous_terms.txt by prevalence threshold and saves
the result to prev{threshold}_ubiquitous_terms.txt.
"""

import pandas as pd
import argparse
from pathlib import Path
from prevalence_utils import get_prevalence_prefix, filter_go_columns_by_prevalence

# Parse arguments
parser = argparse.ArgumentParser(description='Create prevalence-filtered GO terms file')
parser.add_argument('--prevalence-threshold', type=float, required=True,
                    help='Prevalence threshold (0-100) for filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
MASTER_TABLE_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv"
UBIQUITOUS_TERMS_FILE = BASE_DIR / "results/3_GO_analyses/ubiquitous_terms.txt"
OUTPUT_DIR = BASE_DIR / "results/3_GO_analyses"

print(f"Creating prevalence-filtered GO terms file (threshold: {args.prevalence_threshold}%%)...")

# Load master table
print(f"Loading master table from {MASTER_TABLE_FILE}...")
df = pd.read_csv(MASTER_TABLE_FILE, sep='\t')

# Filter GO columns by prevalence
print(f"Filtering GO columns by {args.prevalence_threshold}%% prevalence...")
go_cols_filtered = filter_go_columns_by_prevalence(df, args.prevalence_threshold)

print(f"  Found {len(go_cols_filtered)} GO columns meeting {args.prevalence_threshold}%% threshold")

# Save filtered terms
prefix = get_prevalence_prefix(args.prevalence_threshold)
output_file = OUTPUT_DIR / f"{prefix}ubiquitous_terms.txt"
with open(output_file, 'w') as f:
    for term in sorted(go_cols_filtered):
        f.write(f"{term}\n")

print(f"✓ Saved {len(go_cols_filtered)} terms to {output_file}")

