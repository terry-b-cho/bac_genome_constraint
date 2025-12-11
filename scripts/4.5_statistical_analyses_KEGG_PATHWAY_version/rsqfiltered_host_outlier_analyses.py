#!/usr/bin/env python3
"""
R²-Filtered Host and Non-Host Outlier Analyses

Identifies outlier KEGG categories (pathway/ko/reactions) in host environments
(Mammals, Mammals: Human, Plants) and non-host environments (Aquatic, Terrestrial)
based on alpha and beta distributions.

Only includes fits with R² >= 0.05.
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# Try to import adjustText for better label positioning
try:
    from adjustText import adjust_text
    HAS_ADJUST_TEXT = True
except ImportError:
    HAS_ADJUST_TEXT = False
    print("  ⚠ adjustText not available. Using manual positioning algorithm.")


def position_labels_no_overlap(ax, texts, points, expand_factor=1.5, max_iterations=200):
    """
    Manually position text labels to avoid overlaps using spiral positioning.
    More aggressive spacing to prevent any overlaps.
    
    Args:
        ax: matplotlib axes
        texts: list of Text/Annotation objects
        points: list of (x, y) tuples for each text
        expand_factor: how much to expand around points
        max_iterations: maximum iterations for positioning
    """
    if len(texts) == 0:
        return
    
    # Get axis limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    x_range = xlim[1] - xlim[0]
    y_range = ylim[1] - ylim[0]
    
    # Estimate text size (more generous to ensure no overlap)
    # Assume labels are wider than tall, with padding
    avg_text_width = x_range * 0.12  # Increased from 0.08
    avg_text_height = y_range * 0.05  # Increased from 0.03
    
    # Minimum distance between labels (much larger to prevent overlap)
    min_dist = max(avg_text_width, avg_text_height) * 1.8  # Increased from 1.2
    
    # Initial positions using spiral pattern with larger radius
    positions = []
    for i, (text, point) in enumerate(zip(texts, points)):
        x, y = point
        
        # Spiral out from point with larger initial radius
        angle = (i * 137.5) % 360  # Golden angle for even distribution
        # Start with larger radius and increase more aggressively
        radius = min(x_range, y_range) * 0.10 * (1 + i * 0.15)  # Increased from 0.06 and 0.1
        
        offset_x = np.cos(np.radians(angle)) * radius
        offset_y = np.sin(np.radians(angle)) * radius
        
        new_x = x + offset_x
        new_y = y + offset_y
        
        # Keep within bounds (more margin)
        new_x = np.clip(new_x, xlim[0] + x_range * 0.05, xlim[1] - x_range * 0.05)
        new_y = np.clip(new_y, ylim[0] + y_range * 0.05, ylim[1] - y_range * 0.05)
        
        positions.append((new_x, new_y))
        text.set_position((new_x, new_y))
    
    # Aggressive overlap detection and adjustment
    fig = ax.figure
    fig.canvas.draw()
    
    for iteration in range(max_iterations):
        overlaps = False
        max_overlap = 0
        
        for i in range(len(texts)):
            for j in range(i + 1, len(texts)):
                pos_i = positions[i]
                pos_j = positions[j]
                
                # Check if positions are too close
                dx = pos_i[0] - pos_j[0]
                dy = pos_i[1] - pos_j[1]
                dist = np.sqrt(dx**2 + dy**2)
                
                if dist < min_dist:
                    overlaps = True
                    overlap_amount = min_dist - dist
                    max_overlap = max(max_overlap, overlap_amount)
                    
                    # Push labels apart more aggressively
                    if dist > 0:
                        # Stronger repulsion
                        push_factor = 0.8  # Increased from 0.5
                        push_x = (dx / dist) * (min_dist - dist) * push_factor
                        push_y = (dy / dist) * (min_dist - dist) * push_factor
                        
                        new_pos_i = (pos_i[0] + push_x, pos_i[1] + push_y)
                        new_pos_j = (pos_j[0] - push_x, pos_j[1] - push_y)
                        
                        # Keep within bounds (more margin)
                        new_pos_i = (
                            np.clip(new_pos_i[0], xlim[0] + x_range * 0.05, xlim[1] - x_range * 0.05),
                            np.clip(new_pos_i[1], ylim[0] + y_range * 0.05, ylim[1] - y_range * 0.05)
                        )
                        new_pos_j = (
                            np.clip(new_pos_j[0], xlim[0] + x_range * 0.05, xlim[1] - x_range * 0.05),
                            np.clip(new_pos_j[1], ylim[0] + y_range * 0.05, ylim[1] - y_range * 0.05)
                        )
                        
                        positions[i] = new_pos_i
                        positions[j] = new_pos_j
                        texts[i].set_position(new_pos_i)
                        texts[j].set_position(new_pos_j)
                    else:
                        # If points are exactly on top of each other, push in opposite directions
                        angle_i = (i * 137.5) % 360
                        angle_j = (j * 137.5) % 360
                        push_dist = min_dist * 0.5
                        
                        new_pos_i = (
                            pos_i[0] + np.cos(np.radians(angle_i)) * push_dist,
                            pos_i[1] + np.sin(np.radians(angle_i)) * push_dist
                        )
                        new_pos_j = (
                            pos_j[0] + np.cos(np.radians(angle_j)) * push_dist,
                            pos_j[1] + np.sin(np.radians(angle_j)) * push_dist
                        )
                        
                        new_pos_i = (
                            np.clip(new_pos_i[0], xlim[0] + x_range * 0.05, xlim[1] - x_range * 0.05),
                            np.clip(new_pos_i[1], ylim[0] + y_range * 0.05, ylim[1] - y_range * 0.05)
                        )
                        new_pos_j = (
                            np.clip(new_pos_j[0], xlim[0] + x_range * 0.05, xlim[1] - x_range * 0.05),
                            np.clip(new_pos_j[1], ylim[0] + y_range * 0.05, ylim[1] - y_range * 0.05)
                        )
                        
                        positions[i] = new_pos_i
                        positions[j] = new_pos_j
                        texts[i].set_position(new_pos_i)
                        texts[j].set_position(new_pos_j)
        
        if not overlaps:
            break
        
        # Redraw to update positions
        if iteration % 10 == 0:
            fig.canvas.draw()
    
    # Final check: ensure labels stay reasonably close to points (but allow more distance)
    max_dist = min(x_range, y_range) * 0.35  # Increased from 0.2 to allow more spacing
    
    for text, point, pos in zip(texts, points, positions):
        x, y = point
        px, py = pos
        
        dx = px - x
        dy = py - y
        dist = np.sqrt(dx**2 + dy**2)
        
        if dist > max_dist:
            # Scale back towards point, but not too aggressively
            scale = max_dist / dist
            new_px = x + dx * scale
            new_py = y + dy * scale
            text.set_position((new_px, new_py))

# Add parent directory to path for imports
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")

# Set style
try:
    plt.style.use('seaborn-v0_8')
except OSError:
    try:
        plt.style.use('seaborn')
    except OSError:
        plt.style.use('default')
sns.set_palette("husl")

# Environment groups
HOST_ENVIRONMENTS = ["Mammals", "Mammals: Human", "Plants"]
NONHOST_ENVIRONMENTS = ["Aquatic", "Terrestrial"]

# R² threshold
R_SQUARED_THRESHOLD = 0.05

# Percentile thresholds (top 5% and bottom 5%)
TOP_PERCENTILE = 95
BOTTOM_PERCENTILE = 5


def load_kegg_labels(data_type):
    """Load KEGG labels for human-readable names."""
    labels_file = BASE_DIR / f"results/4.5_statistical_analyses_KEGG_PATHWAY_version/05_kegg_labels/kegg_term_labels_{data_type}.tsv"
    
    if not labels_file.exists():
        print(f"  ⚠ KEGG labels file not found: {labels_file}")
        return {}
    
    try:
        labels_df = pd.read_csv(labels_file, sep='\t')
        labels_map = {}
        if 'category' in labels_df.columns and 'name' in labels_df.columns:
            for _, row in labels_df.iterrows():
                cat = str(row['category'])
                name = str(row['name'])
                labels_map[cat] = name
        return labels_map
    except Exception as e:
        print(f"  ⚠ Error loading KEGG labels: {e}")
        return {}


def load_env_scaling_data(data_type):
    """Load environment-specific scaling parameters."""
    env_scaling_file = BASE_DIR / f"results/4.5_statistical_analyses_KEGG_PATHWAY_version/04_env_scaling/env_scaling_params_{data_type}.tsv"
    
    if not env_scaling_file.exists():
        print(f"  ✗ Environment scaling file not found: {env_scaling_file}")
        return None
    
    try:
        df = pd.read_csv(env_scaling_file, sep='\t')
        print(f"  ✓ Loaded {len(df)} environment×category combinations")
        return df
    except Exception as e:
        print(f"  ✗ Error loading environment scaling data: {e}")
        return None


def filter_environments(df, environment_list, label):
    """Filter for specific environments."""
    filtered_df = df[df['environment'].isin(environment_list)].copy()
    print(f"  ✓ Filtered to {len(filtered_df)} {label} environment×category combinations")
    return filtered_df


def filter_r_squared(df):
    """Filter for R² >= threshold."""
    filtered_df = df[df['r_squared'] >= R_SQUARED_THRESHOLD].copy()
    print(f"  ✓ Filtered to {len(filtered_df)} combinations with R² >= {R_SQUARED_THRESHOLD}")
    return filtered_df


def calculate_percentiles(df, column):
    """Calculate top and bottom percentiles for a column."""
    values = df[column].dropna()
    if len(values) == 0:
        return None, None
    
    top_percentile_value = np.percentile(values, TOP_PERCENTILE)
    bottom_percentile_value = np.percentile(values, BOTTOM_PERCENTILE)
    
    return top_percentile_value, bottom_percentile_value


def identify_outliers(df, column, top_value, bottom_value, labels_map):
    """Identify outlier categories based on percentile thresholds."""
    top_outliers = df[df[column] >= top_value].copy()
    bottom_outliers = df[df[column] <= bottom_value].copy()
    
    # Add human-readable labels
    top_outliers['category_label'] = top_outliers['category'].map(
        lambda x: labels_map.get(str(x), str(x))
    )
    bottom_outliers['category_label'] = bottom_outliers['category'].map(
        lambda x: labels_map.get(str(x), str(x))
    )
    
    # Sort by the column value
    top_outliers = top_outliers.sort_values(column, ascending=False)
    bottom_outliers = bottom_outliers.sort_values(column, ascending=True)
    
    return top_outliers, bottom_outliers


def create_outlier_table(top_outliers, bottom_outliers, column_name, data_type, output_dir, prefix):
    """Create a summary table of outliers."""
    results = []
    
    # Top outliers
    for _, row in top_outliers.iterrows():
        results.append({
            'outlier_type': f'Top 5% ({column_name})',
            'environment': row['environment'],
            'category': row['category'],
            'category_label': row['category_label'],
            'value': row[column_name],
            'r_squared': row['r_squared'],
            'n_genomes_used': row['n_genomes_used'],
            'alpha_env': row.get('alpha_env', np.nan),
            'beta_env_log': row.get('beta_env_log', np.nan)
        })
    
    # Bottom outliers
    for _, row in bottom_outliers.iterrows():
        results.append({
            'outlier_type': f'Bottom 5% ({column_name})',
            'environment': row['environment'],
            'category': row['category'],
            'category_label': row['category_label'],
            'value': row[column_name],
            'r_squared': row['r_squared'],
            'n_genomes_used': row['n_genomes_used'],
            'alpha_env': row.get('alpha_env', np.nan),
            'beta_env_log': row.get('beta_env_log', np.nan)
        })
    
    results_df = pd.DataFrame(results)
    
    # Save table
    output_file = output_dir / f"{prefix}_outliers_{column_name}_{data_type}.tsv"
    results_df.to_csv(output_file, sep='\t', index=False)
    print(f"  ✓ Saved outlier table: {output_file.name}")
    
    return results_df


def create_outlier_visualization(df, column_name, top_outliers, bottom_outliers, 
                                 top_value, bottom_value, data_type, output_dir, labels_map, prefix, env_label):
    """Create visualization of outliers with improved labeling."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    # Left panel: Distribution with outliers highlighted
    ax1 = axes[0]
    
    # Plot all data
    values = df[column_name].dropna()
    ax1.hist(values, bins=50, alpha=0.5, color='gray', edgecolor='black', label='All data')
    
    # Highlight outliers
    top_values = top_outliers[column_name].values
    bottom_values = bottom_outliers[column_name].values
    
    if len(top_values) > 0:
        ax1.hist(top_values, bins=50, alpha=0.8, color='red', edgecolor='darkred', 
                label=f'Top 5% (≥{top_value:.4f})')
    
    if len(bottom_values) > 0:
        ax1.hist(bottom_values, bins=50, alpha=0.8, color='blue', edgecolor='darkblue', 
                label=f'Bottom 5% (≤{bottom_value:.4f})')
    
    # Add percentile lines
    ax1.axvline(top_value, color='red', linestyle='--', linewidth=2, alpha=0.7, label=f'{TOP_PERCENTILE}th percentile')
    ax1.axvline(bottom_value, color='blue', linestyle='--', linewidth=2, alpha=0.7, label=f'{BOTTOM_PERCENTILE}th percentile')
    
    ax1.set_xlabel(f'{column_name.replace("_", " ").title()}', fontsize=12, weight='bold')
    ax1.set_ylabel('Frequency', fontsize=12, weight='bold')
    ax1.set_title(f'Distribution of {column_name.replace("_", " ").title()}\n({env_label} Environments, R² ≥ {R_SQUARED_THRESHOLD})', 
                  fontsize=14, weight='bold')
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # Right panel: Scatter plot with labeled outliers
    ax2 = axes[1]
    
    # Plot all data
    ax2.scatter(df['r_squared'], df[column_name], alpha=0.3, color='gray', s=20, label='All data')
    
    # Add vertical line at R² = 0.05
    ax2.axvline(R_SQUARED_THRESHOLD, color='green', linestyle='--', linewidth=2, alpha=0.7, 
                label=f'R² = {R_SQUARED_THRESHOLD}')
    
    # Highlight and label top outliers
    texts_top = []
    if len(top_outliers) > 0:
        ax2.scatter(top_outliers['r_squared'], top_outliers[column_name], 
                   alpha=0.8, color='red', s=100, edgecolor='darkred', linewidth=1.5,
                   label=f'Top 5% (n={len(top_outliers)})', zorder=5)
        
        # Label top outliers (show top 10 to avoid overcrowding)
        top_to_label = top_outliers.head(10)
        for _, row in top_to_label.iterrows():
            label = row['category_label'][:40] + '...' if len(row['category_label']) > 40 else row['category_label']
            # Start with offset position, will be adjusted
            texts_top.append(ax2.annotate(label,
                                         xy=(row['r_squared'], row[column_name]),
                                         xytext=(5, 5), textcoords='offset points',
                                         fontsize=8, alpha=0.8,
                                         bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7),
                                         arrowprops=dict(arrowstyle='->', color='black', alpha=0.4, lw=0.5)))
    
    # Highlight and label bottom outliers
    texts_bottom = []
    if len(bottom_outliers) > 0:
        ax2.scatter(bottom_outliers['r_squared'], bottom_outliers[column_name], 
                   alpha=0.8, color='blue', s=100, edgecolor='darkblue', linewidth=1.5,
                   label=f'Bottom 5% (n={len(bottom_outliers)})', zorder=5)
        
        # Label bottom outliers (show bottom 10 to avoid overcrowding)
        bottom_to_label = bottom_outliers.head(10)
        for _, row in bottom_to_label.iterrows():
            label = row['category_label'][:40] + '...' if len(row['category_label']) > 40 else row['category_label']
            # Start with offset position, will be adjusted
            texts_bottom.append(ax2.annotate(label,
                                            xy=(row['r_squared'], row[column_name]),
                                            xytext=(5, -15), textcoords='offset points',
                                            fontsize=8, alpha=0.8,
                                            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.7),
                                            arrowprops=dict(arrowstyle='->', color='black', alpha=0.4, lw=0.5)))
    
    # Use adjust_text or manual positioning to repel overlapping labels
    all_texts = texts_top + texts_bottom
    all_points = []
    # Extract point positions from annotation objects
    for text in all_texts:
        # For annotate objects, get the xy position (the point being annotated)
        if hasattr(text, 'xy'):
            all_points.append(text.xy)
        else:
            # Fallback: use text position
            pos = text.get_position()
            all_points.append(pos)
    
    if len(all_texts) > 0:
        if HAS_ADJUST_TEXT:
            try:
                adjust_text(all_texts, ax=ax2, 
                           arrowprops=dict(arrowstyle='->', color='black', alpha=0.5, lw=0.5),
                           expand_points=(2.0, 2.0),  # Increased from 1.2
                           force_text=(1.0, 1.0),  # Increased from 0.5
                           force_points=(0.3, 0.3),  # Reduced to allow more movement
                           lim=500,  # Maximum iterations
                           precision=0.01)  # Higher precision
            except Exception as e:
                print(f"  ⚠ Warning: Could not adjust text positions: {e}")
                # Fall back to manual positioning
                position_labels_no_overlap(ax2, all_texts, all_points)
        else:
            # Use improved manual positioning
            position_labels_no_overlap(ax2, all_texts, all_points)
    
    # Add percentile lines
    ax2.axhline(top_value, color='red', linestyle='--', linewidth=2, alpha=0.7)
    ax2.axhline(bottom_value, color='blue', linestyle='--', linewidth=2, alpha=0.7)
    
    ax2.set_xlabel('R² (Coefficient of Determination)', fontsize=12, weight='bold')
    ax2.set_ylabel(f'{column_name.replace("_", " ").title()}', fontsize=12, weight='bold')
    ax2.set_title(f'Outlier Categories by {column_name.replace("_", " ").title()}\n({env_label} Environments, R² ≥ {R_SQUARED_THRESHOLD})', 
                  fontsize=14, weight='bold')
    ax2.legend(fontsize=10, loc='best')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = output_dir / f"{prefix}_outliers_{column_name}_{data_type}.png"
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved figure: {output_file.name}")
    
    output_file_pdf = output_dir / f"{prefix}_outliers_{column_name}_{data_type}.pdf"
    fig.savefig(output_file_pdf, bbox_inches='tight')
    print(f"  ✓ Saved figure: {output_file_pdf.name}")
    
    plt.close(fig)


def create_combined_summary_table(alpha_top, alpha_bottom, beta_top, beta_bottom, data_type, output_dir, prefix):
    """Create a combined summary table of all outliers."""
    summary = []
    
    # Collect all unique outlier categories
    all_outliers = set()
    
    for df in [alpha_top, alpha_bottom, beta_top, beta_bottom]:
        if len(df) > 0:
            all_outliers.update(df['category'].unique())
    
    # For each outlier category, collect information
    for category in sorted(all_outliers):
        category_info = {
            'category': category,
            'category_label': None,
            'alpha_outlier_type': None,
            'alpha_value': np.nan,
            'alpha_environment': None,
            'beta_outlier_type': None,
            'beta_value': np.nan,
            'beta_environment': None,
            'r_squared': np.nan,
            'n_genomes_used': np.nan
        }
        
        # Check alpha outliers
        alpha_top_match = alpha_top[alpha_top['category'] == category]
        alpha_bottom_match = alpha_bottom[alpha_bottom['category'] == category]
        
        if len(alpha_top_match) > 0:
            row = alpha_top_match.iloc[0]
            category_info['alpha_outlier_type'] = 'Top 5%'
            category_info['alpha_value'] = row['alpha_env']
            category_info['alpha_environment'] = row['environment']
            if category_info['category_label'] is None:
                category_info['category_label'] = row.get('category_label', category)
            if np.isnan(category_info['r_squared']):
                category_info['r_squared'] = row['r_squared']
                category_info['n_genomes_used'] = row['n_genomes_used']
        
        if len(alpha_bottom_match) > 0:
            row = alpha_bottom_match.iloc[0]
            if category_info['alpha_outlier_type'] is None:
                category_info['alpha_outlier_type'] = 'Bottom 5%'
                category_info['alpha_value'] = row['alpha_env']
                category_info['alpha_environment'] = row['environment']
            if category_info['category_label'] is None:
                category_info['category_label'] = row.get('category_label', category)
            if np.isnan(category_info['r_squared']):
                category_info['r_squared'] = row['r_squared']
                category_info['n_genomes_used'] = row['n_genomes_used']
        
        # Check beta outliers
        beta_top_match = beta_top[beta_top['category'] == category]
        beta_bottom_match = beta_bottom[beta_bottom['category'] == category]
        
        if len(beta_top_match) > 0:
            row = beta_top_match.iloc[0]
            category_info['beta_outlier_type'] = 'Top 5%'
            category_info['beta_value'] = row['beta_env_log']
            category_info['beta_environment'] = row['environment']
            if category_info['category_label'] is None:
                category_info['category_label'] = row.get('category_label', category)
            if np.isnan(category_info['r_squared']):
                category_info['r_squared'] = row['r_squared']
                category_info['n_genomes_used'] = row['n_genomes_used']
        
        if len(beta_bottom_match) > 0:
            row = beta_bottom_match.iloc[0]
            if category_info['beta_outlier_type'] is None:
                category_info['beta_outlier_type'] = 'Bottom 5%'
                category_info['beta_value'] = row['beta_env_log']
                category_info['beta_environment'] = row['environment']
            if category_info['category_label'] is None:
                category_info['category_label'] = row.get('category_label', category)
            if np.isnan(category_info['r_squared']):
                category_info['r_squared'] = row['r_squared']
                category_info['n_genomes_used'] = row['n_genomes_used']
        
        summary.append(category_info)
    
    summary_df = pd.DataFrame(summary)
    
    # Save summary table
    output_file = output_dir / f"{prefix}_outliers_summary_{data_type}.tsv"
    summary_df.to_csv(output_file, sep='\t', index=False)
    print(f"  ✓ Saved summary table: {output_file.name}")
    
    return summary_df


def create_comparison_plot(host_df, nonhost_df, host_alpha_top, host_alpha_bottom, host_beta_top, host_beta_bottom,
                          nonhost_alpha_top, nonhost_alpha_bottom, nonhost_beta_top, nonhost_beta_bottom,
                          data_type, output_dir, labels_map):
    """Create faceted comparison plot: beta vs alpha for host vs non-host."""
    fig, axes = plt.subplots(1, 2, figsize=(18, 8))
    
    # Left panel: Host environments
    ax1 = axes[0]
    
    # Plot all host data
    ax1.scatter(host_df['alpha_env'], host_df['beta_env_log'], 
               alpha=0.3, color='gray', s=20, label='All data', zorder=1)
    
    # Highlight outliers
    if len(host_alpha_top) > 0:
        ax1.scatter(host_alpha_top['alpha_env'], host_alpha_top['beta_env_log'],
                   alpha=0.8, color='red', s=150, edgecolor='darkred', linewidth=2,
                   marker='^', label='Top 5% alpha', zorder=5)
    
    if len(host_alpha_bottom) > 0:
        ax1.scatter(host_alpha_bottom['alpha_env'], host_alpha_bottom['beta_env_log'],
                   alpha=0.8, color='blue', s=150, edgecolor='darkblue', linewidth=2,
                   marker='v', label='Bottom 5% alpha', zorder=5)
    
    if len(host_beta_top) > 0:
        ax1.scatter(host_beta_top['alpha_env'], host_beta_top['beta_env_log'],
                   alpha=0.8, color='orange', s=150, edgecolor='darkorange', linewidth=2,
                   marker='>', label='Top 5% beta', zorder=5)
    
    if len(host_beta_bottom) > 0:
        ax1.scatter(host_beta_bottom['alpha_env'], host_beta_bottom['beta_env_log'],
                   alpha=0.8, color='purple', s=150, edgecolor='darkviolet', linewidth=2,
                   marker='<', label='Bottom 5% beta', zorder=5)
    
    # Label top outliers (combine alpha and beta top outliers)
    texts_host = []
    all_host_outliers = pd.concat([
        host_alpha_top.head(5),
        host_alpha_bottom.head(5),
        host_beta_top.head(5),
        host_beta_bottom.head(5)
    ]).drop_duplicates(subset=['category', 'environment'])
    
    for _, row in all_host_outliers.iterrows():
        label = row['category_label'][:30] + '...' if len(row['category_label']) > 30 else row['category_label']
        # Use annotate with arrow
        ann = ax1.annotate(label,
                          xy=(row['alpha_env'], row['beta_env_log']),
                          xytext=(10, 10), textcoords='offset points',
                          fontsize=7, alpha=0.8,
                          bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.6),
                          arrowprops=dict(arrowstyle='->', color='black', alpha=0.4, lw=0.5))
        texts_host.append(ann)
    
    # Get point positions for host labels
    host_points = [(row['alpha_env'], row['beta_env_log']) for _, row in all_host_outliers.iterrows()]
    
    if len(texts_host) > 0:
        if HAS_ADJUST_TEXT:
            try:
                adjust_text(texts_host, ax=ax1, 
                           expand_points=(2.0, 2.0),
                           force_text=(1.0, 1.0),
                           force_points=(0.3, 0.3),
                           lim=500,
                           precision=0.01)
            except Exception as e:
                print(f"  ⚠ Warning: Could not adjust text positions for host: {e}")
                position_labels_no_overlap(ax1, texts_host, host_points)
        else:
            position_labels_no_overlap(ax1, texts_host, host_points)
    
    ax1.set_xlabel('Alpha (Scaling Exponent)', fontsize=12, weight='bold')
    ax1.set_ylabel('Beta (Log Offset)', fontsize=12, weight='bold')
    ax1.set_title('Host Environments\n(Mammals, Mammals: Human, Plants)', fontsize=14, weight='bold')
    ax1.legend(fontsize=9, loc='best')
    ax1.grid(True, alpha=0.3)
    
    # Right panel: Non-host environments
    ax2 = axes[1]
    
    # Plot all non-host data
    ax2.scatter(nonhost_df['alpha_env'], nonhost_df['beta_env_log'], 
               alpha=0.3, color='gray', s=20, label='All data', zorder=1)
    
    # Highlight outliers
    if len(nonhost_alpha_top) > 0:
        ax2.scatter(nonhost_alpha_top['alpha_env'], nonhost_alpha_top['beta_env_log'],
                   alpha=0.8, color='red', s=150, edgecolor='darkred', linewidth=2,
                   marker='^', label='Top 5% alpha', zorder=5)
    
    if len(nonhost_alpha_bottom) > 0:
        ax2.scatter(nonhost_alpha_bottom['alpha_env'], nonhost_alpha_bottom['beta_env_log'],
                   alpha=0.8, color='blue', s=150, edgecolor='darkblue', linewidth=2,
                   marker='v', label='Bottom 5% alpha', zorder=5)
    
    if len(nonhost_beta_top) > 0:
        ax2.scatter(nonhost_beta_top['alpha_env'], nonhost_beta_top['beta_env_log'],
                   alpha=0.8, color='orange', s=150, edgecolor='darkorange', linewidth=2,
                   marker='>', label='Top 5% beta', zorder=5)
    
    if len(nonhost_beta_bottom) > 0:
        ax2.scatter(nonhost_beta_bottom['alpha_env'], nonhost_beta_bottom['beta_env_log'],
                   alpha=0.8, color='purple', s=150, edgecolor='darkviolet', linewidth=2,
                   marker='<', label='Bottom 5% beta', zorder=5)
    
    # Label top outliers
    texts_nonhost = []
    all_nonhost_outliers = pd.concat([
        nonhost_alpha_top.head(5),
        nonhost_alpha_bottom.head(5),
        nonhost_beta_top.head(5),
        nonhost_beta_bottom.head(5)
    ]).drop_duplicates(subset=['category', 'environment'])
    
    for _, row in all_nonhost_outliers.iterrows():
        label = row['category_label'][:30] + '...' if len(row['category_label']) > 30 else row['category_label']
        # Use annotate with arrow
        ann = ax2.annotate(label,
                          xy=(row['alpha_env'], row['beta_env_log']),
                          xytext=(10, 10), textcoords='offset points',
                          fontsize=7, alpha=0.8,
                          bbox=dict(boxstyle='round,pad=0.2', facecolor='lightblue', alpha=0.6),
                          arrowprops=dict(arrowstyle='->', color='black', alpha=0.4, lw=0.5))
        texts_nonhost.append(ann)
    
    # Get point positions for non-host labels
    nonhost_points = [(row['alpha_env'], row['beta_env_log']) for _, row in all_nonhost_outliers.iterrows()]
    
    if len(texts_nonhost) > 0:
        if HAS_ADJUST_TEXT:
            try:
                adjust_text(texts_nonhost, ax=ax2,
                           expand_points=(2.0, 2.0),
                           force_text=(1.0, 1.0),
                           force_points=(0.3, 0.3),
                           lim=500,
                           precision=0.01)
            except Exception as e:
                print(f"  ⚠ Warning: Could not adjust text positions for non-host: {e}")
                position_labels_no_overlap(ax2, texts_nonhost, nonhost_points)
        else:
            position_labels_no_overlap(ax2, texts_nonhost, nonhost_points)
    
    ax2.set_xlabel('Alpha (Scaling Exponent)', fontsize=12, weight='bold')
    ax2.set_ylabel('Beta (Log Offset)', fontsize=12, weight='bold')
    ax2.set_title('Non-Host Environments\n(Aquatic, Terrestrial)', fontsize=14, weight='bold')
    ax2.legend(fontsize=9, loc='best')
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = output_dir / f"host_vs_nonhost_beta_vs_alpha_{data_type}.png"
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved comparison figure: {output_file.name}")
    
    output_file_pdf = output_dir / f"host_vs_nonhost_beta_vs_alpha_{data_type}.pdf"
    fig.savefig(output_file_pdf, bbox_inches='tight')
    print(f"  ✓ Saved comparison figure: {output_file_pdf.name}")
    
    plt.close(fig)


def analyze_environment_group(df, environment_list, label, prefix, data_type, output_dir, labels_map):
    """Analyze a group of environments (host or non-host)."""
    print(f"\n{'='*80}")
    print(f"Analyzing {label} Environments")
    print(f"{'='*80}")
    
    # Filter for environments
    filtered_df = filter_environments(df, environment_list, label)
    
    if len(filtered_df) == 0:
        print(f"  ✗ No data found for {label} environments.")
        return None, None, None, None, None
    
    # Filter for R² >= threshold
    filtered_df = filter_r_squared(filtered_df)
    
    if len(filtered_df) == 0:
        print(f"  ✗ No data passes R² threshold for {label} environments.")
        return None, None, None, None, None
    
    # Calculate percentiles for alpha
    print(f"\nCalculating percentiles for alpha...")
    alpha_top_value, alpha_bottom_value = calculate_percentiles(filtered_df, 'alpha_env')
    print(f"  Top {TOP_PERCENTILE}% threshold: {alpha_top_value:.4f}")
    print(f"  Bottom {BOTTOM_PERCENTILE}% threshold: {alpha_bottom_value:.4f}")
    
    # Identify alpha outliers
    print(f"\nIdentifying alpha outliers...")
    alpha_top_outliers, alpha_bottom_outliers = identify_outliers(
        filtered_df, 'alpha_env', alpha_top_value, alpha_bottom_value, labels_map
    )
    print(f"  Top 5% outliers: {len(alpha_top_outliers)}")
    print(f"  Bottom 5% outliers: {len(alpha_bottom_outliers)}")
    
    # Create alpha outlier table
    print(f"\nCreating alpha outlier table...")
    create_outlier_table(
        alpha_top_outliers, alpha_bottom_outliers, 'alpha_env', data_type, output_dir, prefix
    )
    
    # Create alpha visualization
    print(f"\nCreating alpha visualization...")
    create_outlier_visualization(
        filtered_df, 'alpha_env', alpha_top_outliers, alpha_bottom_outliers,
        alpha_top_value, alpha_bottom_value, data_type, output_dir, labels_map, prefix, label
    )
    
    # Calculate percentiles for beta
    print(f"\nCalculating percentiles for beta...")
    beta_top_value, beta_bottom_value = calculate_percentiles(filtered_df, 'beta_env_log')
    print(f"  Top {TOP_PERCENTILE}% threshold: {beta_top_value:.4f}")
    print(f"  Bottom {BOTTOM_PERCENTILE}% threshold: {beta_bottom_value:.4f}")
    
    # Identify beta outliers
    print(f"\nIdentifying beta outliers...")
    beta_top_outliers, beta_bottom_outliers = identify_outliers(
        filtered_df, 'beta_env_log', beta_top_value, beta_bottom_value, labels_map
    )
    print(f"  Top 5% outliers: {len(beta_top_outliers)}")
    print(f"  Bottom 5% outliers: {len(beta_bottom_outliers)}")
    
    # Create beta outlier table
    print(f"\nCreating beta outlier table...")
    create_outlier_table(
        beta_top_outliers, beta_bottom_outliers, 'beta_env_log', data_type, output_dir, prefix
    )
    
    # Create beta visualization
    print(f"\nCreating beta visualization...")
    create_outlier_visualization(
        filtered_df, 'beta_env_log', beta_top_outliers, beta_bottom_outliers,
        beta_top_value, beta_bottom_value, data_type, output_dir, labels_map, prefix, label
    )
    
    # Create combined summary table
    print(f"\nCreating combined summary table...")
    create_combined_summary_table(
        alpha_top_outliers, alpha_bottom_outliers,
        beta_top_outliers, beta_bottom_outliers,
        data_type, output_dir, prefix
    )
    
    return filtered_df, alpha_top_outliers, alpha_bottom_outliers, beta_top_outliers, beta_bottom_outliers


def main():
    parser = argparse.ArgumentParser(
        description="R²-Filtered Host and Non-Host Outlier Analyses",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        '--data-type',
        type=str,
        required=True,
        choices=['pathway', 'ko', 'reactions'],
        help='KEGG data type to analyze'
    )
    parser.add_argument(
        '--output-dir',
        type=str,
        default=None,
        help='Output directory (default: results/4.5_statistical_analyses_KEGG_PATHWAY_version/rsqfiltered_host_outlier_analyses)'
    )
    
    args = parser.parse_args()
    
    # Set up output directory
    if args.output_dir is None:
        output_dir = BASE_DIR / "results/4.5_statistical_analyses_KEGG_PATHWAY_version/rsqfiltered_host_outlier_analyses"
    else:
        output_dir = Path(args.output_dir)
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("=" * 80)
    print("R²-Filtered Host and Non-Host Outlier Analyses")
    print("=" * 80)
    print(f"Data Type: {args.data_type}")
    print(f"Host Environments: {', '.join(HOST_ENVIRONMENTS)}")
    print(f"Non-Host Environments: {', '.join(NONHOST_ENVIRONMENTS)}")
    print(f"R² Threshold: >= {R_SQUARED_THRESHOLD}")
    print(f"Percentile Thresholds: Top {TOP_PERCENTILE}%, Bottom {BOTTOM_PERCENTILE}%")
    print("=" * 80)
    
    # Load data
    print("\nLoading data...")
    df = load_env_scaling_data(args.data_type)
    if df is None:
        print("  ✗ Failed to load data. Exiting.")
        return 1
    
    # Load KEGG labels
    print("Loading KEGG labels...")
    labels_map = load_kegg_labels(args.data_type)
    print(f"  ✓ Loaded {len(labels_map)} KEGG labels")
    
    # Analyze host environments
    host_results = analyze_environment_group(
        df, HOST_ENVIRONMENTS, "Host", "host", args.data_type, output_dir, labels_map
    )
    
    # Analyze non-host environments
    nonhost_results = analyze_environment_group(
        df, NONHOST_ENVIRONMENTS, "Non-Host", "nonhost", args.data_type, output_dir, labels_map
    )
    
    # Create comparison plot if both analyses succeeded
    if (host_results[0] is not None and nonhost_results[0] is not None):
        print(f"\n{'='*80}")
        print("Creating Host vs Non-Host Comparison Plot")
        print(f"{'='*80}")
        create_comparison_plot(
            host_results[0], nonhost_results[0],
            host_results[1], host_results[2], host_results[3], host_results[4],
            nonhost_results[1], nonhost_results[2], nonhost_results[3], nonhost_results[4],
            args.data_type, output_dir, labels_map
        )
    
    print("\n" + "=" * 80)
    print("Analysis Complete!")
    print("=" * 80)
    print(f"Output directory: {output_dir}")
    print("=" * 80)
    
    return 0


if __name__ == '__main__':
    sys.exit(main())
