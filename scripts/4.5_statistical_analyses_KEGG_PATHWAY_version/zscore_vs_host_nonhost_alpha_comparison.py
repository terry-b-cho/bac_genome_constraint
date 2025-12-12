#!/usr/bin/env python3
"""
Compare Z-scores with host vs non-host alpha differences.

This script analyzes the relationship between category-level Z-scores
and the difference in average scaling exponents (alpha) between
host environments (Mammals, Mammals: Human, Plants) and non-host environments
(Aquatic, Terrestrial).

Author: Generated from terminal exploration code
Date: December 2025
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats
from pathlib import Path
import sys

# ============================================================================
# Configuration
# ============================================================================

BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
DATA_TYPE = "pathway"  # Can be changed to "ko" or "reactions"

# Input files
ENV_SCALING_FILE = BASE_DIR / f"results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/env_scaling_params_{DATA_TYPE}.tsv"
RSQFILTERED_Z_SUMMARY_FILE = BASE_DIR / f"results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/rsqfiltered_category_Z_summary_{DATA_TYPE}.tsv"
KEGG_LABELS_FILE = BASE_DIR / f"results/4.5_statistical_analyses_KEGG_PATHWAY_version/05_kegg_labels/kegg_term_labels_{DATA_TYPE}.tsv"

# Output directory
OUTPUT_DIR = BASE_DIR / f"results/4.5_statistical_analyses_KEGG_PATHWAY_version/zscore_host_comparison"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Environment definitions
HOST_ENVIRONMENTS = ["Mammals: Human", "Mammals", "Plants"]
NONHOST_ENVIRONMENTS = ["Aquatic", "Terrestrial"]

# R-squared threshold
R_SQUARED_THRESHOLD = 0.05

# ============================================================================
# Load Data
# ============================================================================

print("=" * 80)
print("Z-Score vs Host/Non-Host Alpha Comparison Analysis")
print("=" * 80)
print(f"Data Type: {DATA_TYPE}")
print(f"R² Threshold: {R_SQUARED_THRESHOLD}")
print()

# Load environment scaling parameters (filtered by R²)
print("Loading environment scaling parameters...")
df = pd.read_csv(ENV_SCALING_FILE, sep="\t")
df = df[df["r_squared"] >= R_SQUARED_THRESHOLD].copy()
print(f"  ✓ Loaded {len(df)} environment×category combinations with R² >= {R_SQUARED_THRESHOLD}")

# Load rsqfiltered Z-scores
print("Loading rsqfiltered category Z-scores...")
zdf = pd.read_csv(RSQFILTERED_Z_SUMMARY_FILE, sep="\t", index_col="category")
zdf.sort_values("Z_alpha_category", ascending=False, inplace=True)
print(f"  ✓ Loaded {len(zdf)} categories with rsqfiltered Z-scores")

# Load KEGG labels
print("Loading KEGG labels...")
try:
    namedf = pd.read_csv(KEGG_LABELS_FILE, sep="\t", index_col="category", usecols=["category", "name"])
    print(f"  ✓ Loaded {len(namedf)} KEGG labels")
except Exception as e:
    print(f"  ⚠ Could not load KEGG labels: {e}")
    namedf = pd.DataFrame()

# ============================================================================
# Compute Statistics
# ============================================================================

print()
print("Computing host vs non-host alpha averages...")

z_scores = []
alpha_avg_host = []
alpha_avg_nonhost = []
differences = []
categories = []

for cat, z_row in zdf.iterrows():
    # Get environment-specific data for this category
    cat_data = df[df["category"] == cat].copy()
    
    if len(cat_data) == 0:
        continue
    
    # Separate host and non-host environments
    host_rows = cat_data[cat_data["environment"].isin(HOST_ENVIRONMENTS)]
    nonhost_rows = cat_data[cat_data["environment"].isin(NONHOST_ENVIRONMENTS)]
    
    # Compute averages
    z_score = z_row["Z_alpha_category"]
    host_alpha = host_rows["alpha_env"].mean() if len(host_rows) > 0 else np.nan
    nonhost_alpha = nonhost_rows["alpha_env"].mean() if len(nonhost_rows) > 0 else np.nan
    difference = host_alpha - nonhost_alpha if (not np.isnan(host_alpha) and not np.isnan(nonhost_alpha)) else np.nan
    
    z_scores.append(z_score)
    alpha_avg_host.append(host_alpha)
    alpha_avg_nonhost.append(nonhost_alpha)
    differences.append(difference)
    categories.append(cat)

print(f"  ✓ Computed statistics for {len(categories)} categories")

# Convert to numpy arrays for easier handling
z_scores = np.array(z_scores)
alpha_avg_host = np.array(alpha_avg_host)
alpha_avg_nonhost = np.array(alpha_avg_nonhost)
differences = np.array(differences)

# Filter out NaN values for regression
valid_mask = ~np.isnan(differences)
z_scores_valid = z_scores[valid_mask]
differences_valid = differences[valid_mask]
categories_valid = [cat for i, cat in enumerate(categories) if valid_mask[i]]

print(f"  ✓ {len(z_scores_valid)} categories have valid differences (non-NaN)")

# ============================================================================
# Linear Regression
# ============================================================================

print()
print("Computing linear regression...")
if len(z_scores_valid) > 1:
    regress_result = scipy.stats.linregress(z_scores_valid, differences_valid)
    print(f"  Slope: {regress_result.slope:.6f}")
    print(f"  Intercept: {regress_result.intercept:.6f}")
    print(f"  R²: {regress_result.rvalue**2:.6f}")
    print(f"  p-value: {regress_result.pvalue:.6f}")
else:
    regress_result = None
    print("  ⚠ Insufficient data for regression")

# ============================================================================
# Create Plots
# ============================================================================

print()
print("Generating plots...")

# Set style
plt.style.use('default')
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['font.size'] = 10

# ----------------------------------------------------------------------------
# Plot 1: Z-score vs Alpha Difference (Scatter)
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 6))

ax.scatter(z_scores_valid, differences_valid, alpha=0.6, s=50, color='steelblue', edgecolors='black', linewidth=0.5)
ax.axhline(y=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Zero difference line')

# Add regression line if available
if regress_result is not None and not np.isnan(regress_result.slope):
    x_line = np.linspace(z_scores_valid.min(), z_scores_valid.max(), 100)
    y_line = regress_result.slope * x_line + regress_result.intercept
    ax.plot(x_line, y_line, 'r-', linewidth=2, alpha=0.8, 
            label=f'Linear fit (R²={regress_result.rvalue**2:.3f}, p={regress_result.pvalue:.3f})')
    
    # Calculate perpendicular distances from points to regression line
    # Distance = |mx - y + b| / sqrt(m^2 + 1)
    # But we need signed distance (positive = above line, negative = below line)
    m = regress_result.slope
    b = regress_result.intercept
    distances = []
    for x, y in zip(z_scores_valid, differences_valid):
        # Perpendicular distance with sign
        # Points above the line have positive distance, below have negative
        distance = (m * x - y + b) / np.sqrt(m**2 + 1)
        distances.append(distance)
    
    distances = np.array(distances)
    
    # Find top 3 and bottom 3 points (furthest orthogonally from line)
    top3_idx = np.argsort(distances)[-3:][::-1]  # Top 3 (largest positive distances)
    bottom3_idx = np.argsort(distances)[:3]      # Bottom 3 (most negative distances)
    
    # Annotate top 3 points (furthest above the regression line)
    offsets_top = [(15, 20), (15, 40), (15, 60)]  # Vary vertical offset to avoid overlap
    for i, idx in enumerate(top3_idx):
        x_val = z_scores_valid[idx]
        y_val = differences_valid[idx]
        cat = categories_valid[idx]
        
        # Get pathway name
        if len(namedf) > 0 and cat in namedf.index:
            label = namedf.at[cat, "name"]
            # Truncate if too long
            if len(label) > 35:
                label = label[:32] + "..."
        else:
            label = cat
        
        # Add annotation with arrow
        ax.annotate(label, 
                   xy=(x_val, y_val),
                   xytext=offsets_top[i],
                   textcoords='offset points',
                   fontsize=9,
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.8, edgecolor='black', linewidth=0.8),
                   arrowprops=dict(arrowstyle='->', color='black', lw=1.2, alpha=0.7, connectionstyle='arc3,rad=0.1'),
                   ha='left',
                   va='bottom')
    
    # Annotate bottom 3 points (furthest below the regression line)
    offsets_bottom = [(15, -20), (15, -40), (15, -60)]  # Vary vertical offset to avoid overlap
    for i, idx in enumerate(bottom3_idx):
        x_val = z_scores_valid[idx]
        y_val = differences_valid[idx]
        cat = categories_valid[idx]
        
        # Get pathway name
        if len(namedf) > 0 and cat in namedf.index:
            label = namedf.at[cat, "name"]
            # Truncate if too long
            if len(label) > 35:
                label = label[:32] + "..."
        else:
            label = cat
        
        # Add annotation with arrow
        ax.annotate(label, 
                   xy=(x_val, y_val),
                   xytext=offsets_bottom[i],
                   textcoords='offset points',
                   fontsize=9,
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.8, edgecolor='black', linewidth=0.8),
                   arrowprops=dict(arrowstyle='->', color='black', lw=1.2, alpha=0.7, connectionstyle='arc3,rad=0.1'),
                   ha='left',
                   va='top')

ax.set_xlabel('Category Z-Score (Z_alpha_category)', fontsize=12, weight='bold')
ax.set_ylabel('Alpha Difference (Host - Non-Host)', fontsize=12, weight='bold')
ax.set_title(f'Z-Score vs Host/Non-Host Alpha Difference\n({DATA_TYPE.title()} Data, R² ≥ {R_SQUARED_THRESHOLD})', 
             fontsize=14, weight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)

plt.tight_layout()
output_file = OUTPUT_DIR / f"zscore_vs_alpha_difference_{DATA_TYPE}.png"
fig.savefig(output_file, bbox_inches='tight')
print(f"  ✓ Saved: {output_file.name}")

output_file_pdf = OUTPUT_DIR / f"zscore_vs_alpha_difference_{DATA_TYPE}.pdf"
fig.savefig(output_file_pdf, bbox_inches='tight')
print(f"  ✓ Saved: {output_file_pdf.name}")

plt.close(fig)

# ----------------------------------------------------------------------------
# Plot 2: Z-score vs Average Alphas (Host and Non-Host)
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 6))

# Filter valid data for both host and non-host
valid_both_mask = ~(np.isnan(alpha_avg_host) | np.isnan(alpha_avg_nonhost))
z_scores_both = z_scores[valid_both_mask]
alpha_host_both = alpha_avg_host[valid_both_mask]
alpha_nonhost_both = alpha_avg_nonhost[valid_both_mask]

ax.scatter(z_scores_both, alpha_host_both, alpha=0.6, s=50, color='red', 
          edgecolors='darkred', linewidth=0.5, label='Host environments (avg)')
ax.scatter(z_scores_both, alpha_nonhost_both, alpha=0.6, s=50, color='blue', 
          edgecolors='darkblue', linewidth=0.5, label='Non-host environments (avg)')

ax.set_xlabel('Category Z-Score (Z_alpha_category)', fontsize=12, weight='bold')
ax.set_ylabel('Average Scaling Exponent (α)', fontsize=12, weight='bold')
ax.set_title(f'Z-Score vs Average Alpha by Environment Type\n({DATA_TYPE.title()} Data, R² ≥ {R_SQUARED_THRESHOLD})', 
             fontsize=14, weight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)

plt.tight_layout()
output_file = OUTPUT_DIR / f"zscore_vs_avg_alpha_by_env_{DATA_TYPE}.png"
fig.savefig(output_file, bbox_inches='tight')
print(f"  ✓ Saved: {output_file.name}")

output_file_pdf = OUTPUT_DIR / f"zscore_vs_avg_alpha_by_env_{DATA_TYPE}.pdf"
fig.savefig(output_file_pdf, bbox_inches='tight')
print(f"  ✓ Saved: {output_file_pdf.name}")

plt.close(fig)

# ----------------------------------------------------------------------------
# Plot 3: Z-score vs Alpha Difference (Line Plot)
# ----------------------------------------------------------------------------

fig, ax = plt.subplots(figsize=(10, 6))

# Sort by Z-score for line plot
sort_idx = np.argsort(z_scores_valid)
z_sorted = z_scores_valid[sort_idx]
diff_sorted = differences_valid[sort_idx]

ax.plot(z_sorted, diff_sorted, 'o-', markersize=4, linewidth=1.5, alpha=0.7, color='steelblue')
ax.axhline(y=0, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Zero difference line')

ax.set_xlabel('Category Z-Score (Z_alpha_category)', fontsize=12, weight='bold')
ax.set_ylabel('Alpha Difference (Host - Non-Host)', fontsize=12, weight='bold')
ax.set_title(f'Z-Score vs Host/Non-Host Alpha Difference (Line Plot)\n({DATA_TYPE.title()} Data, R² ≥ {R_SQUARED_THRESHOLD})', 
             fontsize=14, weight='bold')
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)

plt.tight_layout()
output_file = OUTPUT_DIR / f"zscore_vs_alpha_difference_line_{DATA_TYPE}.png"
fig.savefig(output_file, bbox_inches='tight')
print(f"  ✓ Saved: {output_file.name}")

output_file_pdf = OUTPUT_DIR / f"zscore_vs_alpha_difference_line_{DATA_TYPE}.pdf"
fig.savefig(output_file_pdf, bbox_inches='tight')
print(f"  ✓ Saved: {output_file_pdf.name}")

plt.close(fig)

# ============================================================================
# Summary Statistics
# ============================================================================

print()
print("=" * 80)
print("Summary Statistics")
print("=" * 80)
print(f"Total categories analyzed: {len(categories)}")
print(f"Categories with valid differences: {len(z_scores_valid)}")
print(f"Categories with negative differences (host < non-host): {np.sum(differences_valid < 0)}")
print(f"Categories with positive differences (host > non-host): {np.sum(differences_valid > 0)}")
print()

if len(differences_valid) > 0:
    print(f"Mean difference: {np.mean(differences_valid):.6f}")
    print(f"Median difference: {np.median(differences_valid):.6f}")
    print(f"Std difference: {np.std(differences_valid):.6f}")
    print(f"Min difference: {np.min(differences_valid):.6f}")
    print(f"Max difference: {np.max(differences_valid):.6f}")
    print()

# List categories with negative differences (if labels available)
if len(namedf) > 0:
    print("Categories with negative differences (host alpha < non-host alpha):")
    print("-" * 80)
    for i, (cat, diff, z) in enumerate(zip(categories_valid, differences_valid, z_scores_valid)):
        if diff < 0:
            cat_name = namedf.at[cat, "name"] if cat in namedf.index else cat
            print(f"  {cat:12s}  Z={z:6.3f}  diff={diff:7.4f}  {cat_name}")
    print()

# ============================================================================
# Save Summary Table
# ============================================================================

summary_df = pd.DataFrame({
    'category': categories,
    'Z_alpha_category': z_scores,
    'alpha_avg_host': alpha_avg_host,
    'alpha_avg_nonhost': alpha_avg_nonhost,
    'alpha_difference': differences
})

# Add pathway names if available
if len(namedf) > 0:
    summary_df['pathway_name'] = summary_df['category'].apply(
        lambda x: namedf.at[x, "name"] if x in namedf.index else x
    )

summary_file = OUTPUT_DIR / f"zscore_host_comparison_summary_{DATA_TYPE}.tsv"
summary_df.to_csv(summary_file, sep='\t', index=False)
print(f"✓ Saved summary table: {summary_file.name}")

print()
print("=" * 80)
print("Analysis complete!")
print("=" * 80)

