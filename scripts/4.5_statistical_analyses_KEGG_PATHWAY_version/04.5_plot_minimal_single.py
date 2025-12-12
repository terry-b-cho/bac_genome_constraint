#!/usr/bin/env python3
"""
Script 04.5: Minimal Plotting for KEGG Data (Single Type)

Generates essential QC and diagnostic plots.
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import argparse
from prevalence_utils import get_prevalence_prefix, get_data_type_suffix

# Parse arguments
parser = argparse.ArgumentParser(description='Generate minimal QC plots for KEGG data')
parser.add_argument('--data-type', type=str, required=True, choices=['reactions', 'ko', 'pathway'],
                    help='Data type to process')
parser.add_argument('--test-mode', action='store_true',
                    help='Run on small test subset')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
OUTPUT_BASE_DIR = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/04.5_intermediate_figures"
OUTPUT_BASE_DIR.mkdir(parents=True, exist_ok=True)

data_type = args.data_type
suffix = get_data_type_suffix(data_type)
prefix = get_prevalence_prefix(args.prevalence_threshold)

print(f"Generating minimal plots for {data_type.upper()}...")

# ============================================================================
# Load Data
# ============================================================================

try:
    # Global scaling parameters
    global_params = pd.read_csv(
        BASE_DIR / f"results/4.5_statistical_analyses_KEGG_PATHWAY_version/03_global_scaling/{prefix}global_scaling_params{suffix}.tsv",
        sep='\t'
    )
    print(f"  ✓ Loaded {len(global_params)} global parameters")
    
    # Z-scores
    z_scores = pd.read_csv(
        BASE_DIR / f"results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/{prefix}category_Z_summary{suffix}.tsv",
        sep='\t'
    )
    print(f"  ✓ Loaded {len(z_scores)} Z-score summaries")
    
    # Master table
    master = pd.read_csv(
        BASE_DIR / f"results/4.5_statistical_analyses_KEGG_PATHWAY_version/02_env_cohorts/{prefix}master_table_env_filtered{suffix}.tsv",
        sep='\t',
        nrows=1000 if args.test_mode else None
    )
    print(f"  ✓ Loaded master table ({len(master)} rows)")
    
except Exception as e:
    print(f"  ✗ ERROR loading data: {e}")
    import sys
    sys.exit(1)

# ============================================================================
# Plot 1: Distribution of Alpha Values
# ============================================================================

print("  Generating Plot 1: Alpha distribution...")
fig, ax = plt.subplots(figsize=(10, 6))
ax.hist(global_params['alpha_global'], bins=50, edgecolor='black', alpha=0.7)
ax.axvline(global_params['alpha_global'].median(), color='red', linestyle='--', label=f"Median: {global_params['alpha_global'].median():.3f}")
ax.set_xlabel('Alpha (Scaling Exponent)', fontsize=12)
ax.set_ylabel('Count', fontsize=12)
ax.set_title(f'Distribution of Global Scaling Exponents - {data_type.upper()}', fontsize=14)
ax.legend()
ax.grid(alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_BASE_DIR / f"{prefix}alpha_distribution{suffix}.pdf", dpi=300)
plt.close()

# ============================================================================
# Plot 2: Alpha vs R²
# ============================================================================

print("  Generating Plot 2: Alpha vs R²...")
fig, ax = plt.subplots(figsize=(10, 6))
scatter = ax.scatter(global_params['alpha_global'], global_params['r_squared'], 
                     alpha=0.5, s=50, c=global_params['r_squared'], cmap='viridis')
ax.set_xlabel('Alpha (Scaling Exponent)', fontsize=12)
ax.set_ylabel('R² (Fit Quality)', fontsize=12)
ax.set_title(f'Scaling Exponent vs Fit Quality - {data_type.upper()}', fontsize=14)
ax.grid(alpha=0.3)
plt.colorbar(scatter, label='R²')
plt.tight_layout()
plt.savefig(OUTPUT_BASE_DIR / f"{prefix}alpha_vs_r2{suffix}.pdf", dpi=300)
plt.close()

# ============================================================================
# Plot 3: Distribution of Z-scores
# ============================================================================

print("  Generating Plot 3: Z-score distributions...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

ax1.hist(z_scores['Z_alpha_category'], bins=30, edgecolor='black', alpha=0.7, color='steelblue')
ax1.axvline(z_scores['Z_alpha_category'].median(), color='red', linestyle='--', label=f"Median: {z_scores['Z_alpha_category'].median():.3f}")
ax1.set_xlabel('Z_alpha (Environmental Variation)', fontsize=12)
ax1.set_ylabel('Count', fontsize=12)
ax1.set_title(f'Z-Score Distribution (Alpha) - {data_type.upper()}', fontsize=13)
ax1.legend()
ax1.grid(alpha=0.3)

ax2.hist(z_scores['Z_beta_category'], bins=30, edgecolor='black', alpha=0.7, color='coral')
ax2.axvline(z_scores['Z_beta_category'].median(), color='red', linestyle='--', label=f"Median: {z_scores['Z_beta_category'].median():.3f}")
ax2.set_xlabel('Z_beta (Environmental Variation)', fontsize=12)
ax2.set_ylabel('Count', fontsize=12)
ax2.set_title(f'Z-Score Distribution (Beta) - {data_type.upper()}', fontsize=13)
ax2.legend()
ax2.grid(alpha=0.3)

plt.tight_layout()
plt.savefig(OUTPUT_BASE_DIR / f"{prefix}z_score_distributions{suffix}.pdf", dpi=300)
plt.close()

# ============================================================================
# Plot 4: Top Variable Categories
# ============================================================================

print("  Generating Plot 4: Top variable categories...")
top_n = 20
top_variable = z_scores.nlargest(top_n, 'Z_alpha_category')

fig, ax = plt.subplots(figsize=(12, 8))
y_pos = np.arange(len(top_variable))
ax.barh(y_pos, top_variable['Z_alpha_category'], color='steelblue', alpha=0.7)
ax.set_yticks(y_pos)
ax.set_yticklabels(top_variable['category'], fontsize=9)
ax.set_xlabel('Z_alpha (Environmental Variation)', fontsize=12)
ax.set_title(f'Top {top_n} Most Environmentally Variable Categories - {data_type.upper()}', fontsize=14)
ax.invert_yaxis()
ax.grid(axis='x', alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_BASE_DIR / f"{prefix}top_variable_categories{suffix}.pdf", dpi=300)
plt.close()

# ============================================================================
# Plot 5: Environment Counts
# ============================================================================

print("  Generating Plot 5: Environment genome counts...")
env_counts = master.groupby('environment').size().reset_index(name='n_genomes')
env_counts = env_counts.sort_values('n_genomes', ascending=True)

fig, ax = plt.subplots(figsize=(10, 6))
ax.barh(env_counts['environment'], env_counts['n_genomes'], color='forestgreen', alpha=0.7)
ax.set_xlabel('Number of Genomes', fontsize=12)
ax.set_ylabel('Environment', fontsize=12)
ax.set_title(f'Genome Counts per Environment - {data_type.upper()}', fontsize=14)
ax.grid(axis='x', alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_BASE_DIR / f"{prefix}environment_counts{suffix}.pdf", dpi=300)
plt.close()

# ============================================================================
# Plot 6: Alpha vs Number of Environments Used
# ============================================================================

print("  Generating Plot 6: Alpha vs environments used...")
fig, ax = plt.subplots(figsize=(10, 6))
scatter = ax.scatter(global_params['alpha_global'], global_params['n_genomes_used'], 
                     alpha=0.5, s=50, c=global_params['r_squared'], cmap='plasma')
ax.set_xlabel('Alpha (Scaling Exponent)', fontsize=12)
ax.set_ylabel('Number of Genomes Used', fontsize=12)
ax.set_title(f'Scaling Exponent vs Sample Size - {data_type.upper()}', fontsize=14)
ax.grid(alpha=0.3)
plt.colorbar(scatter, label='R²')
plt.tight_layout()
plt.savefig(OUTPUT_BASE_DIR / f"{prefix}alpha_vs_sample_size{suffix}.pdf", dpi=300)
plt.close()

print("")
print("=" * 80)
print(f"Minimal plots generated successfully for {data_type}!")
print(f"Output directory: {OUTPUT_BASE_DIR}")
print("=" * 80)



