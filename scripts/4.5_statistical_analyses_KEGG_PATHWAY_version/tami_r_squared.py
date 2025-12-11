#!/usr/bin/env python3
"""
Script: tami_r_squared

Generate histograms of R² values from global and environment-specific scaling fits.
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import numpy as np

# Try to import scipy for density plots
try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
GLOBAL_SCALING_FILE = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/03_global_scaling/global_scaling_params_pathway.tsv"
ENV_SCALING_FILE = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/env_scaling_params_pathway.tsv"
OUTPUT_DIR = BASE_DIR / "results/tami_request"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

print("=" * 80)
print("Script: tami_r_squared - R² Histogram Analysis")
print("=" * 80)
print("")

# Load data
print("Loading data...")
try:
    global_df = pd.read_csv(GLOBAL_SCALING_FILE, sep='\t')
    print(f"  ✓ Loaded global scaling: {len(global_df)} categories")
except Exception as e:
    print(f"  ✗ ERROR loading global scaling file: {e}")
    exit(1)

try:
    env_df = pd.read_csv(ENV_SCALING_FILE, sep='\t')
    print(f"  ✓ Loaded environment scaling: {len(env_df)} env×category fits")
except Exception as e:
    print(f"  ✗ ERROR loading environment scaling file: {e}")
    exit(1)

# Check for r_squared columns
if 'r_squared' not in global_df.columns:
    print("  ✗ ERROR: 'r_squared' column not found in global scaling file")
    exit(1)

if 'r_squared' not in env_df.columns:
    print("  ✗ ERROR: 'r_squared' column not found in environment scaling file")
    exit(1)

if 'environment' not in env_df.columns:
    print("  ✗ ERROR: 'environment' column not found in environment scaling file")
    exit(1)

print("")

# ============================================================================
# Figure 1: Global R² Histogram
# ============================================================================

print("Generating Figure 1: Global R² histogram...")
try:
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Filter out any NaN or invalid R² values
    r2_global = global_df['r_squared'].dropna()
    r2_global = r2_global[r2_global >= 0]  # R² should be >= 0
    
    ax.hist(r2_global, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax.axvline(r2_global.median(), color='red', linestyle='--', linewidth=2,
               label=f'Median = {r2_global.median():.3f}')
    ax.axvline(r2_global.mean(), color='green', linestyle='--', linewidth=2,
               label=f'Mean = {r2_global.mean():.3f}')
    
    ax.set_xlabel('R² (Coefficient of Determination)', fontsize=12)
    ax.set_ylabel('Number of Categories', fontsize=12)
    ax.set_title('Distribution of R² Values - Global Scaling Fits', fontsize=14, weight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add text with summary stats
    stats_text = f'n = {len(r2_global)}\nMin = {r2_global.min():.3f}\nMax = {r2_global.max():.3f}\nSD = {r2_global.std():.3f}'
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    output_file = OUTPUT_DIR / "global_r_squared_histogram.png"
    fig.savefig(output_file, bbox_inches='tight', dpi=300)
    output_file_pdf = OUTPUT_DIR / "global_r_squared_histogram.pdf"
    fig.savefig(output_file_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ Saved: {output_file}")
    print(f"  ✓ Saved: {output_file_pdf}")
    
except Exception as e:
    print(f"  ✗ ERROR generating global R² histogram: {e}")
    import traceback
    traceback.print_exc()

print("")

# ============================================================================
# Figure 2: Environment-Specific R² Histograms (Faceted)
# ============================================================================

print("Generating Figure 2: Environment-specific R² histograms (faceted)...")
try:
    # Filter out NaN and invalid R² values
    env_df_clean = env_df[env_df['r_squared'].notna() & (env_df['r_squared'] >= 0)].copy()
    
    # Get unique environments and sort them
    environments = sorted(env_df_clean['environment'].unique())
    n_envs = len(environments)
    
    print(f"  Found {n_envs} environments: {', '.join(environments)}")
    
    # Determine grid layout
    # Try to make it roughly square
    n_cols = int(np.ceil(np.sqrt(n_envs)))
    n_rows = int(np.ceil(n_envs / n_cols))
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(4*n_cols, 3*n_rows))
    
    # Flatten axes array for easier indexing
    if n_rows == 1:
        axes = axes.reshape(1, -1) if n_cols > 1 else [axes]
    elif n_cols == 1:
        axes = axes.reshape(-1, 1)
    axes_flat = axes.flatten()
    
    for idx, env in enumerate(environments):
        ax = axes_flat[idx]
        
        # Get R² values for this environment
        env_r2 = env_df_clean[env_df_clean['environment'] == env]['r_squared']
        
        # Create histogram
        ax.hist(env_r2, bins=30, color='coral', edgecolor='black', alpha=0.7)
        
        # Add median and mean lines
        median_val = env_r2.median()
        mean_val = env_r2.mean()
        ax.axvline(median_val, color='red', linestyle='--', linewidth=1.5,
                  label=f'Median = {median_val:.3f}')
        ax.axvline(mean_val, color='green', linestyle='--', linewidth=1.5,
                  label=f'Mean = {mean_val:.3f}')
        
        ax.set_xlabel('R²', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title(f'{env}\n(n={len(env_r2)})', fontsize=11, weight='bold')
        ax.legend(fontsize=8, loc='upper right')
        ax.grid(True, alpha=0.3, axis='y')
        
        # Set x-axis limits to be consistent (0 to 1, or max if > 1)
        max_r2 = env_r2.max()
        ax.set_xlim(0, min(1.0, max_r2 * 1.1) if max_r2 <= 1.0 else max_r2 * 1.1)
    
    # Hide unused subplots
    for idx in range(len(environments), len(axes_flat)):
        axes_flat[idx].axis('off')
    
    plt.suptitle('Distribution of R² Values - Environment-Specific Scaling Fits', 
                 fontsize=14, weight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.99])  # Leave space for suptitle
    
    output_file = OUTPUT_DIR / "env_r_squared_histograms_faceted.png"
    fig.savefig(output_file, bbox_inches='tight', dpi=300)
    output_file_pdf = OUTPUT_DIR / "env_r_squared_histograms_faceted.pdf"
    fig.savefig(output_file_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ Saved: {output_file}")
    print(f"  ✓ Saved: {output_file_pdf}")
    
except Exception as e:
    print(f"  ✗ ERROR generating environment R² histograms: {e}")
    import traceback
    traceback.print_exc()

print("")

# ============================================================================
# Figure 3: Side-by-Side Comparison (Global vs Environment-Specific)
# ============================================================================

print("Generating Figure 3: Side-by-side comparison (Global vs Environment-Specific)...")
try:
    # Filter data
    r2_global = global_df['r_squared'].dropna()
    r2_global = r2_global[r2_global >= 0]
    env_df_clean = env_df[env_df['r_squared'].notna() & (env_df['r_squared'] >= 0)].copy()
    
    # Create figure with two subplots side by side (30% wider than 8)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 3))
    
    # ========================================================================
    # Left plot: Global R² Histogram
    # ========================================================================
    ax1.hist(r2_global, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
    ax1.axvline(r2_global.median(), color='red', linestyle='--', linewidth=2,
               label=f'Median = {r2_global.median():.3f}')
    ax1.axvline(r2_global.mean(), color='green', linestyle='--', linewidth=2,
               label=f'Mean = {r2_global.mean():.3f}')
    
    ax1.set_xlabel('R² (Coefficient of Determination)', fontsize=12)
    ax1.set_ylabel('Number of Categories', fontsize=12)
    ax1.set_title('Global Scaling Fits', fontsize=14, weight='bold')
    ax1.legend(fontsize=10, loc='center right')
    ax1.grid(True, alpha=0.3, axis='y')
    
    # Add text with summary stats
    stats_text = f'n = {len(r2_global)}\nMin = {r2_global.min():.3f}\nMax = {r2_global.max():.3f}\nSD = {r2_global.std():.3f}'
    ax1.text(0.05, 0.95, stats_text, transform=ax1.transAxes, fontsize=10,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # ========================================================================
    # Right plot: Environment-Specific R² Histogram (colored by environment) + Density
    # ========================================================================
    
    # Get unique environments and assign colors
    environments = sorted(env_df_clean['environment'].unique())
    n_envs = len(environments)
    
    # Use a colormap to assign colors to each environment
    colors = plt.cm.tab10(np.linspace(0, 1, n_envs))
    env_color_map = {env: colors[i] for i, env in enumerate(environments)}
    
    # Determine bin edges (use same bins for all environments for consistency)
    all_r2 = env_df_clean['r_squared'].values
    bins = np.linspace(0, min(1.0, all_r2.max() * 1.1), 50)
    
    # Plot histogram bars colored by environment (overlapping, not stacked)
    # Each environment gets its own histogram with transparency
    for env in environments:
        env_r2 = env_df_clean[env_df_clean['environment'] == env]['r_squared'].values
        ax2.hist(env_r2, bins=bins, alpha=0.25, color=env_color_map[env], 
                edgecolor='none', label=env, density=False, histtype='stepfilled')
    
    # Overlay density plots for each environment
    if HAS_SCIPY:
        try:
            bin_width = bins[1] - bins[0]
            x_density = np.linspace(0, min(1.0, all_r2.max() * 1.1), 200)
            
            # Plot density curve for each environment
            for env in environments:
                env_r2 = env_df_clean[env_df_clean['environment'] == env]['r_squared'].values
                if len(env_r2) > 1:  # Need at least 2 points for KDE
                    kde = stats.gaussian_kde(env_r2)
                    y_density = kde(x_density)
                    
                    # Scale density to match histogram scale (multiply by count and bin width)
                    env_count = len(env_r2)
                    y_density_scaled = y_density * env_count * bin_width
                    
                    # Plot density curve overlayed with same color as histogram
                    ax2.plot(x_density, y_density_scaled, color=env_color_map[env], 
                           linewidth=1.5, alpha=1.0, zorder=10)
        except Exception as e:
            print(f"  ⚠ Could not generate density plots: {e}")
    else:
        print("  ⚠ scipy not available, skipping density plot overlay")
    
    ax2.set_xlabel('R² (Coefficient of Determination)', fontsize=12)
    ax2.set_ylabel('Number of Fits', fontsize=12)
    ax2.set_title('Environment-Specific Scaling Fits\n(Colored by Environment)', fontsize=14, weight='bold')
    ax2.legend(fontsize=8, loc='center right', ncol=1)
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add summary stats text
    stats_text_env = f'n = {len(env_df_clean)}\nMean = {env_df_clean["r_squared"].mean():.3f}\nMedian = {env_df_clean["r_squared"].median():.3f}'
    ax2.text(0.05, 0.95, stats_text_env, transform=ax2.transAxes, fontsize=10,
           verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    output_file = OUTPUT_DIR / "r_squared_comparison_side_by_side.png"
    fig.savefig(output_file, bbox_inches='tight', dpi=300)
    output_file_pdf = OUTPUT_DIR / "r_squared_comparison_side_by_side.pdf"
    fig.savefig(output_file_pdf, bbox_inches='tight')
    plt.close(fig)
    print(f"  ✓ Saved: {output_file}")
    print(f"  ✓ Saved: {output_file_pdf}")
    
except Exception as e:
    print(f"  ✗ ERROR generating side-by-side comparison: {e}")
    import traceback
    traceback.print_exc()

print("")

# ============================================================================
# Summary Statistics
# ============================================================================

print("Summary Statistics:")
print("-" * 80)

# Global R² summary
r2_global = global_df['r_squared'].dropna()
r2_global = r2_global[r2_global >= 0]
print(f"\nGlobal Scaling R²:")
print(f"  n = {len(r2_global)}")
print(f"  Mean = {r2_global.mean():.4f}")
print(f"  Median = {r2_global.median():.4f}")
print(f"  SD = {r2_global.std():.4f}")
print(f"  Min = {r2_global.min():.4f}")
print(f"  Max = {r2_global.max():.4f}")
print(f"  Q25 = {r2_global.quantile(0.25):.4f}")
print(f"  Q75 = {r2_global.quantile(0.75):.4f}")

# Environment-specific R² summary
env_df_clean = env_df[env_df['r_squared'].notna() & (env_df['r_squared'] >= 0)].copy()
r2_env = env_df_clean['r_squared']
print(f"\nEnvironment-Specific Scaling R² (all environments combined):")
print(f"  n = {len(r2_env)}")
print(f"  Mean = {r2_env.mean():.4f}")
print(f"  Median = {r2_env.median():.4f}")
print(f"  SD = {r2_env.std():.4f}")
print(f"  Min = {r2_env.min():.4f}")
print(f"  Max = {r2_env.max():.4f}")
print(f"  Q25 = {r2_env.quantile(0.25):.4f}")
print(f"  Q75 = {r2_env.quantile(0.75):.4f}")

# Per-environment summary
print(f"\nEnvironment-Specific Scaling R² (by environment):")
for env in sorted(env_df_clean['environment'].unique()):
    env_r2 = env_df_clean[env_df_clean['environment'] == env]['r_squared']
    print(f"  {env}:")
    print(f"    n = {len(env_r2)}, Mean = {env_r2.mean():.4f}, Median = {env_r2.median():.4f}")

# Save summary statistics to file
summary_file = OUTPUT_DIR / "r_squared_summary_statistics.tsv"
summary_data = []

# Global summary
summary_data.append({
    'dataset': 'global',
    'environment': 'all',
    'n_fits': len(r2_global),
    'mean_r2': r2_global.mean(),
    'median_r2': r2_global.median(),
    'sd_r2': r2_global.std(),
    'min_r2': r2_global.min(),
    'max_r2': r2_global.max(),
    'q25_r2': r2_global.quantile(0.25),
    'q75_r2': r2_global.quantile(0.75)
})

# Per-environment summaries
for env in sorted(env_df_clean['environment'].unique()):
    env_r2 = env_df_clean[env_df_clean['environment'] == env]['r_squared']
    summary_data.append({
        'dataset': 'environment_specific',
        'environment': env,
        'n_fits': len(env_r2),
        'mean_r2': env_r2.mean(),
        'median_r2': env_r2.median(),
        'sd_r2': env_r2.std(),
        'min_r2': env_r2.min(),
        'max_r2': env_r2.max(),
        'q25_r2': env_r2.quantile(0.25),
        'q75_r2': env_r2.quantile(0.75)
    })

summary_df = pd.DataFrame(summary_data)
summary_df.to_csv(summary_file, sep='\t', index=False)
print(f"\n  ✓ Saved summary statistics: {summary_file}")

print("")
print("=" * 80)
print("Script completed successfully!")
print(f"All outputs saved to: {OUTPUT_DIR}")
print("=" * 80)

