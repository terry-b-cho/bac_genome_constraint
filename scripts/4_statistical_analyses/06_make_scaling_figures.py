#!/usr/bin/env python3
"""
Script 06: Make Scaling Figures

Generate environment-stratified version of Fig. 1 (panels a-k) from van Nimwegen (2003),
showing Z-statistics, exponent comparisons, and scatter plots.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import sys
import argparse
from prevalence_utils import get_prevalence_prefix

# Try to import plotting libraries
try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib/seaborn not available")

# Parse arguments
parser = argparse.ArgumentParser(description='Generate scaling figures (Fig. 1 panels a-k)')
parser.add_argument('--test-mode', action='store_true',
                    help='Generate only a subset of figures for testing')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for GO term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths (will be set dynamically based on prevalence threshold)
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
OUTPUT_DIR = BASE_DIR / "results/4_statistical_analyses/06_figures"
QC_LOG_FILE = OUTPUT_DIR / "qc_06_figures.log"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Set style
if HAS_MATPLOTLIB:
    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300

# Initialize QC log
qc_log = []

def log_message(msg):
    """Add message to QC log and print to stdout."""
    qc_log.append(msg)
    print(msg)

def save_figure(fig, filename):
    """Save figure in both PNG and PDF formats."""
    if not HAS_MATPLOTLIB:
        log_message(f"  ✗ Cannot save {filename}: matplotlib not available")
        return
    prefix = get_prevalence_prefix(args.prevalence_threshold)
    filename_with_prefix = f"{prefix}{filename}" if prefix else filename
    png_path = OUTPUT_DIR / f"{filename_with_prefix}.png"
    pdf_path = OUTPUT_DIR / f"{filename_with_prefix}.pdf"
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    fig.savefig(pdf_path, bbox_inches='tight')
    log_message(f"  ✓ Saved {filename_with_prefix}.png and .pdf")
    plt.close(fig)

log_message("=" * 80)
log_message("Script 06: Make Scaling Figures")
if args.test_mode:
    log_message("  TEST MODE: Generating subset of figures")
if args.prevalence_threshold is not None:
    log_message(f"  PREVALENCE FILTER: {args.prevalence_threshold}%% threshold")
if not HAS_MATPLOTLIB:
    log_message("  ✗ ERROR: matplotlib/seaborn required")
    sys.exit(1)
log_message("=" * 80)
log_message("")

# ============================================================================
# Load Data
# ============================================================================

log_message("Loading data...")
try:
    prefix = get_prevalence_prefix(args.prevalence_threshold)
    CATEGORY_Z_FILE = BASE_DIR / "results/4_statistical_analyses/04_env_scaling" / (f"{prefix}category_Z_summary.tsv" if prefix else "category_Z_summary.tsv")
    ENV_SCALING_FILE = BASE_DIR / "results/4_statistical_analyses/04_env_scaling" / (f"{prefix}env_scaling_params.tsv" if prefix else "env_scaling_params.tsv")
    GLOBAL_PARAMS_FILE = BASE_DIR / "results/4_statistical_analyses/03_global_scaling" / (f"{prefix}global_scaling_params.tsv" if prefix else "global_scaling_params.tsv")
    MASTER_TABLE_FILE = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv"
    GO_LABELS_FILE = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels_for_plots.tsv" if prefix else "go_term_labels_for_plots.tsv")
    
    category_z = pd.read_csv(CATEGORY_Z_FILE, sep='\t')
    env_scaling = pd.read_csv(ENV_SCALING_FILE, sep='\t')
    global_params = pd.read_csv(GLOBAL_PARAMS_FILE, sep='\t')
    master_table = pd.read_csv(MASTER_TABLE_FILE, sep='\t')
    go_labels = pd.read_csv(GO_LABELS_FILE, sep='\t')
    
    # Ensure category columns are strings
    category_z['category'] = category_z['category'].astype(str).str.zfill(7)
    env_scaling['category'] = env_scaling['category'].astype(str).str.zfill(7)
    global_params['category'] = global_params['category'].astype(str).str.zfill(7)
    go_labels['category'] = go_labels['category'].astype(str).str.zfill(7)
    
    # Load metabolic GO terms (with prevalence prefix if applicable)
    base_name_metabolic = "metabolic_go_terms.txt"
    metabolic_terms_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}{base_name_metabolic}" if prefix else base_name_metabolic)
    metabolic_terms = None
    if metabolic_terms_file.exists():
        with open(metabolic_terms_file, 'r') as f:
            metabolic_terms = set(line.strip().zfill(7) for line in f if line.strip())
        log_message(f"  ✓ Loaded {len(metabolic_terms)} metabolism-related GO terms from {metabolic_terms_file.name}")
    
    log_message(f"  ✓ Loaded {len(category_z)} categories, {len(env_scaling)} env×category fits")
    log_message(f"  ✓ Loaded {len(master_table)} genomes, {len(go_labels)} GO labels")
except Exception as e:
    log_message(f"  ✗ ERROR: Failed to load data: {e}")
    sys.exit(1)

log_message("")

# ============================================================================
# Panel 1a: Z-statistics for exponents by category
# ============================================================================

log_message("Generating Panel 1a: Z-statistics for exponents by category...")
try:
    category_z_sorted = category_z.sort_values('Z_alpha_category', ascending=False)
    
    # Load full GO labels for text annotations
    try:
        full_labels = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/05_go_labels/go_term_labels.tsv",
            sep='\t'
        )
        full_labels['category'] = full_labels['category'].astype(str).str.zfill(7)
        category_z_sorted['category'] = category_z_sorted['category'].astype(str).str.zfill(7)
        category_z_sorted = category_z_sorted.merge(
            full_labels[['category', 'name', 'go_id']], 
            on='category', 
            how='left'
        )
        category_z_sorted['label'] = category_z_sorted.apply(
            lambda row: f"{row['name']} ({row['go_id']})" if pd.notna(row['name']) 
            else f"Category {row['category']}", axis=1
        )
    except:
        category_z_sorted['label'] = category_z_sorted['category'].apply(
            lambda x: f"Category {x}"
        )
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    x_pos = np.arange(len(category_z_sorted))
    bars = ax.bar(x_pos, category_z_sorted['Z_alpha_category'], 
                  color='steelblue', edgecolor='black', alpha=0.7)
    
    ax.axhline(y=2, color='red', linestyle='--', linewidth=2, label='Z = 2 (threshold)')
    
    # Add text labels for top categories (selective: top 15 by Z, or Z > 2.0)
    # Use a smarter approach: label top 15, plus any with Z > 2.0 that aren't in top 15
    top_15 = category_z_sorted.head(15)
    high_z = category_z_sorted[category_z_sorted['Z_alpha_category'] > 2.0]
    cats_to_label = pd.concat([top_15, high_z]).drop_duplicates().head(20)
    
    for idx, (_, row) in enumerate(category_z_sorted.iterrows()):
        if row.name in cats_to_label.index:
            pos = list(category_z_sorted.index).index(row.name)
            # Truncate label if too long
            label_text = row['label']
            if len(label_text) > 50:
                # Try to keep GO ID visible
                if '(' in label_text:
                    name_part = label_text.split('(')[0].strip()
                    go_part = label_text.split('(')[1].rstrip(')')
                    label_text = name_part[:40] + '...\n(' + go_part + ')'
                else:
                    label_text = label_text[:50] + '...'
            
            # Position label above bar with slight offset
            y_pos = row['Z_alpha_category'] + 0.15
            ax.text(pos, y_pos, label_text, 
                   rotation=45, ha='left', va='bottom', fontsize=6, 
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5))
    
    ax.set_xlabel('GO Category (sorted by Z_alpha_category)', fontsize=12)
    ax.set_ylabel('Z_alpha_category', fontsize=12)
    ax.set_title('Z-Statistics for Exponents by Category', fontsize=14, weight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_xticks([])  # Too many categories to label on x-axis
    
    save_figure(fig, "fig1a_Z_exponents_by_category_env")
    log_message("  ✓ Panel 1a completed")
    
    # Generate metabolic version
    if metabolic_terms:
        log_message("  Generating metabolic version of Panel 1a...")
        metabolic_category_z = category_z[category_z['category'].isin(metabolic_terms)].copy()
        if len(metabolic_category_z) > 0:
            metabolic_category_z_sorted = metabolic_category_z.sort_values('Z_alpha_category', ascending=False)
            
            # Merge labels
            try:
                full_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels.tsv" if prefix else "go_term_labels.tsv")
                full_labels = pd.read_csv(full_labels_file, sep='\t')
                full_labels['category'] = full_labels['category'].astype(str).str.zfill(7)
                metabolic_category_z_sorted['category'] = metabolic_category_z_sorted['category'].astype(str).str.zfill(7)
                metabolic_category_z_sorted = metabolic_category_z_sorted.merge(
                    full_labels[['category', 'name', 'go_id']], 
                    on='category', 
                    how='left'
                )
                metabolic_category_z_sorted['label'] = metabolic_category_z_sorted.apply(
                    lambda row: f"{row['name']} ({row['go_id']})" if pd.notna(row['name']) 
                    else f"Category {row['category']}", axis=1
                )
            except:
                metabolic_category_z_sorted['label'] = metabolic_category_z_sorted['category'].apply(
                    lambda x: f"Category {x}"
                )
            
            fig_met, ax_met = plt.subplots(figsize=(16, 8))
            x_pos_met = np.arange(len(metabolic_category_z_sorted))
            bars_met = ax_met.bar(x_pos_met, metabolic_category_z_sorted['Z_alpha_category'], 
                                  color='steelblue', edgecolor='black', alpha=0.7)
            ax_met.axhline(y=2, color='red', linestyle='--', linewidth=2, label='Z = 2 (threshold)')
            
            # Add labels for top categories
            top_15_met = metabolic_category_z_sorted.head(15)
            high_z_met = metabolic_category_z_sorted[metabolic_category_z_sorted['Z_alpha_category'] > 2.0]
            cats_to_label_met = pd.concat([top_15_met, high_z_met]).drop_duplicates().head(20)
            
            for idx, (_, row) in enumerate(metabolic_category_z_sorted.iterrows()):
                if row.name in cats_to_label_met.index:
                    pos = list(metabolic_category_z_sorted.index).index(row.name)
                    label_text = row['label']
                    if len(label_text) > 50:
                        if '(' in label_text:
                            name_part = label_text.split('(')[0].strip()
                            go_part = label_text.split('(')[1].rstrip(')')
                            label_text = name_part[:40] + '...\n(' + go_part + ')'
                        else:
                            label_text = label_text[:50] + '...'
                    y_pos = row['Z_alpha_category'] + 0.15
                    ax_met.text(pos, y_pos, label_text, 
                               rotation=45, ha='left', va='bottom', fontsize=6, 
                               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5))
            
            ax_met.set_xlabel('GO Category - Metabolism Only (sorted by Z_alpha_category)', fontsize=12)
            ax_met.set_ylabel('Z_alpha_category', fontsize=12)
            ax_met.set_title('Z-Statistics for Exponents by Category (Metabolism Focus)', fontsize=14, weight='bold')
            ax_met.legend(fontsize=10)
            ax_met.grid(True, alpha=0.3, axis='y')
            ax_met.set_xticks([])
            save_figure(fig_met, "metabolic_fig1a_Z_exponents_by_category_env")
            log_message(f"  ✓ Metabolic Panel 1a completed ({len(metabolic_category_z_sorted)} categories)")
    
except Exception as e:
    log_message(f"  ✗ ERROR generating Panel 1a: {e}")
    import traceback
    log_message(traceback.format_exc())

# ============================================================================
# Panel 1b: Z-statistics for offsets by category
# ============================================================================

log_message("")
log_message("Generating Panel 1b: Z-statistics for offsets by category...")
try:
    category_z_sorted_beta = category_z.sort_values('Z_beta_category', ascending=False)
    
    # Load full GO labels for text annotations
    try:
        full_labels = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/05_go_labels/go_term_labels.tsv",
            sep='\t'
        )
        full_labels['category'] = full_labels['category'].astype(str).str.zfill(7)
        category_z_sorted_beta['category'] = category_z_sorted_beta['category'].astype(str).str.zfill(7)
        category_z_sorted_beta = category_z_sorted_beta.merge(
            full_labels[['category', 'name', 'go_id']], 
            on='category', 
            how='left'
        )
        category_z_sorted_beta['label'] = category_z_sorted_beta.apply(
            lambda row: f"{row['name']} ({row['go_id']})" if pd.notna(row['name']) 
            else f"Category {row['category']}", axis=1
        )
    except:
        category_z_sorted_beta['label'] = category_z_sorted_beta['category'].apply(
            lambda x: f"Category {x}"
        )
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    x_pos = np.arange(len(category_z_sorted_beta))
    bars = ax.bar(x_pos, category_z_sorted_beta['Z_beta_category'],
                  color='steelblue', edgecolor='black', alpha=0.7)
    
    ax.axhline(y=2, color='red', linestyle='--', linewidth=2, label='Z = 2 (threshold)')
    
    # Add text labels for top categories (selective: top 15 by Z, or Z > 2.0)
    top_15_beta = category_z_sorted_beta.head(15)
    high_z_beta = category_z_sorted_beta[category_z_sorted_beta['Z_beta_category'] > 2.0]
    cats_to_label_beta = pd.concat([top_15_beta, high_z_beta]).drop_duplicates().head(20)
    
    for idx, (_, row) in enumerate(category_z_sorted_beta.iterrows()):
        if row.name in cats_to_label_beta.index:
            pos = list(category_z_sorted_beta.index).index(row.name)
            # Truncate label if too long
            label_text = row['label']
            if len(label_text) > 50:
                # Try to keep GO ID visible
                if '(' in label_text:
                    name_part = label_text.split('(')[0].strip()
                    go_part = label_text.split('(')[1].rstrip(')')
                    label_text = name_part[:40] + '...\n(' + go_part + ')'
                else:
                    label_text = label_text[:50] + '...'
            
            # Position label above bar with slight offset
            y_pos = row['Z_beta_category'] + 0.15
            ax.text(pos, y_pos, label_text, 
                   rotation=45, ha='left', va='bottom', fontsize=6, 
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='gray', linewidth=0.5))
    
    ax.set_xlabel('GO Category (sorted by Z_beta_category)', fontsize=12)
    ax.set_ylabel('Z_beta_category', fontsize=12)
    ax.set_title('Z-Statistics for Offsets by Category', fontsize=14, weight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_xticks([])
    
    save_figure(fig, "fig1b_Z_offsets_by_category_env")
    log_message("  ✓ Panel 1b completed")
    
except Exception as e:
    log_message(f"  ✗ ERROR generating Panel 1b: {e}")
    import traceback
    log_message(traceback.format_exc())

# ============================================================================
# Panels 1c-1e: Exponent comparisons for selected categories
# ============================================================================

log_message("")
log_message("Generating Panels 1c-1e: Exponent comparisons for selected categories...")
try:
    # Select more categories: low Z, medium Z, high Z, and additional interesting ones
    low_z_cat = category_z.nsmallest(1, 'Z_alpha_category')['category'].iloc[0]
    med_z_cat = category_z.iloc[len(category_z)//2]['category']
    high_z_cat = category_z.nlargest(1, 'Z_alpha_category')['category'].iloc[0]
    
    # Select additional categories: top 5 by Z, bottom 2, and some in the middle
    top_5_z = category_z.nlargest(5, 'Z_alpha_category')['category'].tolist()
    bottom_2_z = category_z.nsmallest(2, 'Z_alpha_category')['category'].tolist()
    quartile_25 = category_z.iloc[len(category_z)//4]['category']
    quartile_75 = category_z.iloc[3*len(category_z)//4]['category']
    
    # Combine and deduplicate
    selected_cats = list(dict.fromkeys([low_z_cat, med_z_cat, high_z_cat] + 
                                       top_5_z[:3] + bottom_2_z + [quartile_25, quartile_75]))[:10]
    panel_labels = [f'1c-{chr(ord("c")+i)}' for i in range(len(selected_cats))]
    
    # Get environment order (consistent across panels)
    envs = sorted(env_scaling['environment'].unique())
    
    # Create subplot grid: 2 rows, 5 columns for up to 10 categories
    n_cats = len(selected_cats)
    n_cols = min(5, n_cats)
    n_rows = (n_cats + n_cols - 1) // n_cols
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    if n_cats == 1:
        axes = [axes]
    elif n_rows == 1:
        axes = axes if isinstance(axes, list) else [axes]
    else:
        axes = axes.flatten()
    
    for idx, (cat, label) in enumerate(zip(selected_cats, panel_labels)):
        if idx >= len(axes):
            break
        ax = axes[idx]
        
        # Get env-specific exponents for this category
        cat_env_data = env_scaling[env_scaling['category'] == cat].copy()
        cat_env_data = cat_env_data.sort_values('environment')
        
        # Get global exponent
        global_row = global_params[global_params['category'] == cat]
        if len(global_row) > 0:
            global_alpha = global_row.iloc[0]['alpha_global']
        else:
            global_alpha = np.nan
        
        # Plot env-specific exponents with CI
        x_pos = np.arange(len(cat_env_data))
        alphas = cat_env_data['alpha_env'].values
        ci_low = cat_env_data['alpha_env_ci99_low'].values
        ci_high = cat_env_data['alpha_env_ci99_high'].values
        
        ax.errorbar(x_pos, alphas, yerr=[alphas - ci_low, ci_high - alphas],
                   fmt='o', capsize=5, capthick=2, markersize=8,
                   color='steelblue', label='Environment-specific')
        
        # Global exponent line
        if not np.isnan(global_alpha):
            ax.axhline(y=global_alpha, color='red', linestyle='--', linewidth=2,
                      label=f'Global (α={global_alpha:.3f})')
        
        # Get category label - try full_label first, then fallback
        cat_label = f"Category {cat}"  # Default
        try:
            # Try to get full_label from go_term_labels_for_plots.tsv
            plot_labels = pd.read_csv(
                BASE_DIR / "results/4_statistical_analyses/05_go_labels/go_term_labels_for_plots.tsv",
                sep='\t'
            )
            plot_labels['category'] = plot_labels['category'].astype(str).str.zfill(7)
            plot_row = plot_labels[plot_labels['category'] == cat]
            if len(plot_row) > 0 and 'full_label' in plot_row.columns:
                cat_label = plot_row.iloc[0]['full_label']
            else:
                # Fallback to go_term_labels.tsv
                full_labels = pd.read_csv(
                    BASE_DIR / "results/4_statistical_analyses/05_go_labels/go_term_labels.tsv",
                    sep='\t'
                )
                full_labels['category'] = full_labels['category'].astype(str).str.zfill(7)
                full_row = full_labels[full_labels['category'] == cat]
                if len(full_row) > 0:
                    cat_label = f"{full_row.iloc[0]['name']} ({full_row.iloc[0]['go_id']})"
        except:
            pass
        
        ax.set_xlabel('Environment', fontsize=11)
        ax.set_ylabel('Scaling Exponent (α)', fontsize=11)
        ax.set_title(f'{label}\n{cat_label}', fontsize=12, weight='bold')
        ax.set_xticks(x_pos)
        ax.set_xticklabels(cat_env_data['environment'], rotation=45, ha='right', fontsize=9)
        ax.legend(fontsize=9)
        ax.grid(True, alpha=0.3, axis='y')
    
    # Hide unused subplots
    for idx in range(len(selected_cats), len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    save_figure(fig, "fig1cde_env_exponents_selected_categories")
    log_message(f"  ✓ Panels 1c-1e completed ({len(selected_cats)} categories)")
    log_message(f"    Selected categories: {len(selected_cats)} total (low Z, medium Z, high Z, and additional)")
    
except Exception as e:
    log_message(f"  ✗ ERROR generating Panels 1c-1e: {e}")
    import traceback
    log_message(traceback.format_exc())

# ============================================================================
# Panels 1f-1k: Scatter plots with fits
# ============================================================================

log_message("")
log_message("Generating Panels 1f-1k: Scatter plots with fits...")
try:
    # Use more categories from panels c-e (up to 10)
    # For each category, show 2-3 representative environments
    # Select environments with most genomes for each category
    env_counts = master_table.groupby('environment').size().sort_values(ascending=False)
    top_envs = env_counts.head(3).index.tolist()
    
    # Create more panels: up to 10 categories × 3 environments = 30 panels max
    # But we'll limit to 15 panels for readability
    n_cats_to_plot = min(5, len(selected_cats))  # Plot top 5 categories
    n_envs_per_cat = 3
    n_panels = n_cats_to_plot * n_envs_per_cat
    
    n_rows = n_cats_to_plot
    n_cols = n_envs_per_cat
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    if n_rows == 1:
        axes = axes.reshape(1, -1) if n_cols > 1 else [axes]
    axes = axes.flatten()
    
    panel_metadata = []
    panel_idx = 0
    
    for cat_idx, cat in enumerate(selected_cats[:n_cats_to_plot]):
        cat_str = str(cat).zfill(7)
        
        # Get environments that have data for this category
        cat_envs = env_scaling[env_scaling['category'] == cat_str]['environment'].unique()
        # Use top environments that have data
        envs_to_plot = [e for e in top_envs if e in cat_envs][:n_envs_per_cat]
        
        for env_idx, env in enumerate(envs_to_plot):
            if panel_idx >= n_panels:
                break
            
            ax = axes[panel_idx]
            
            # Get data
            env_data = master_table[master_table['environment'] == env].copy()
            subset = env_data[(env_data['genes_total'] > 0) & (env_data[cat_str] > 0)].copy()
            
            if len(subset) == 0:
                continue
            
            # Plot in log-log space (log10)
            # Get all genomes with both genes_total > 0 and category > 0
            all_subset = master_table[(master_table['genes_total'] > 0) & 
                                     (master_table[cat_str] > 0)].copy()
            x_all = np.log10(all_subset['genes_total'].values)
            y_all = np.log10(all_subset[cat_str].values)
            x_env = np.log10(subset['genes_total'].values)
            y_env = np.log10(subset[cat_str].values)
            
            # Gray dots: all genomes
            ax.scatter(x_all, y_all, alpha=0.2, s=10, color='lightgray', label='All genomes')
            
            # Colored dots: specific environment
            colors = plt.cm.tab10(cat_idx)
            ax.scatter(x_env, y_env, alpha=0.6, s=30, color=colors, 
                      edgecolors='black', linewidth=0.5, label=env)
            
            # Get fits
            env_fit = env_scaling[(env_scaling['category'] == cat_str) & 
                                 (env_scaling['environment'] == env)]
            global_row = global_params[global_params['category'] == cat_str]
            
            # Environment-specific fit line
            if len(env_fit) > 0:
                alpha_env = env_fit.iloc[0]['alpha_env']
                beta_env = env_fit.iloc[0]['beta_env_log']
                x_line = np.linspace(x_env.min(), x_env.max(), 100)
                y_line = beta_env / np.log(10) + alpha_env * x_line
                ax.plot(x_line, y_line, color=colors, linewidth=2.5, 
                       linestyle='-', label=f'Env fit (α={alpha_env:.3f})')
            
            # Global fit line
            if len(global_row) > 0:
                alpha_global = global_row.iloc[0]['alpha_global']
                beta_global = global_row.iloc[0]['beta_global_log']
                x_line = np.linspace(x_env.min(), x_env.max(), 100)
                y_line = beta_global / np.log(10) + alpha_global * x_line
                ax.plot(x_line, y_line, color='gray', linewidth=2, 
                       linestyle='--', alpha=0.8, label=f'Global (α={alpha_global:.3f})')
            
            panel_letter = chr(ord('f') + panel_idx)
            
            # Get GO label - try full_label first, then fallback to constructing it
            cat_label = f"Category {cat_str}"  # Default
            try:
                # Try to get full_label from go_term_labels_for_plots.tsv
                plot_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels_for_plots.tsv" if prefix else "go_term_labels_for_plots.tsv")
                plot_labels = pd.read_csv(plot_labels_file, sep='\t')
                plot_labels['category'] = plot_labels['category'].astype(str).str.zfill(7)
                plot_row = plot_labels[plot_labels['category'] == cat_str]
                if len(plot_row) > 0 and 'full_label' in plot_row.columns:
                    cat_label = plot_row.iloc[0]['full_label']
                else:
                    # Fallback to go_term_labels.tsv
                    full_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels.tsv" if prefix else "go_term_labels.tsv")
                    full_labels = pd.read_csv(full_labels_file, sep='\t')
                    full_labels['category'] = full_labels['category'].astype(str).str.zfill(7)
                    full_row = full_labels[full_labels['category'] == cat_str]
                    if len(full_row) > 0:
                        cat_label = f"{full_row.iloc[0]['name']} ({full_row.iloc[0]['go_id']})"
            except Exception as e:
                # If all else fails, use default
                pass
            
            ax.set_xlabel('log10(genes_total)', fontsize=10)
            ax.set_ylabel(f'log10({cat_str} count)', fontsize=10)
            ax.set_title(f'Panel {panel_letter}: {env}\n{cat_label}', 
                        fontsize=11, weight='bold')
            ax.legend(fontsize=8, loc='best')
            ax.grid(True, alpha=0.3)
            
            short_label = cat_label  # Use full label for metadata
            
            panel_metadata.append({
                'panel': panel_letter,
                'category': cat_str,
                'short_label': short_label,
                'environment': env,
                'notes': f'Category {cat_str} in {env}',
                'axes_log_transformed': 'log10'
            })
            
            panel_idx += 1
    
    # Hide unused subplots
    for idx in range(panel_idx, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    save_figure(fig, "fig1f_to_k_env_scatter_scaling")
    log_message(f"  ✓ Panels 1f-1k (log scale) completed ({panel_idx} panels)")
    
    # Create linear scale version
    fig_linear, axes_linear = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    if n_rows == 1:
        axes_linear = axes_linear.reshape(1, -1) if n_cols > 1 else [axes_linear]
    axes_linear = axes_linear.flatten()
    
    panel_idx_linear = 0
    
    for cat_idx, cat in enumerate(selected_cats[:n_cats_to_plot]):
        cat_str = str(cat).zfill(7)
        
        # Get environments that have data for this category
        cat_envs = env_scaling[env_scaling['category'] == cat_str]['environment'].unique()
        # Use top environments that have data
        envs_to_plot = [e for e in top_envs if e in cat_envs][:n_envs_per_cat]
        
        for env_idx, env in enumerate(envs_to_plot):
            if panel_idx_linear >= n_panels:
                break
            
            ax = axes_linear[panel_idx_linear]
            
            # Get data
            env_data = master_table[master_table['environment'] == env].copy()
            subset = env_data[(env_data['genes_total'] > 0) & (env_data[cat_str] > 0)].copy()
            
            if len(subset) == 0:
                continue
            
            # Plot in linear space (raw values)
            all_subset = master_table[(master_table['genes_total'] > 0) & 
                                     (master_table[cat_str] > 0)].copy()
            x_all = all_subset['genes_total'].values
            y_all = all_subset[cat_str].values
            x_env = subset['genes_total'].values
            y_env = subset[cat_str].values
            
            # Gray dots: all genomes
            ax.scatter(x_all, y_all, alpha=0.2, s=10, color='lightgray', label='All genomes')
            
            # Colored dots: specific environment
            colors = plt.cm.tab10(cat_idx)
            ax.scatter(x_env, y_env, alpha=0.6, s=30, color=colors, 
                      edgecolors='black', linewidth=0.5, label=env)
            
            # Get fits
            env_fit = env_scaling[(env_scaling['category'] == cat_str) & 
                                 (env_scaling['environment'] == env)]
            global_row = global_params[global_params['category'] == cat_str]
            
            # Environment-specific fit line (convert from log space)
            if len(env_fit) > 0:
                alpha_env = env_fit.iloc[0]['alpha_env']
                beta_env = env_fit.iloc[0]['beta_env_log']
                x_line = np.linspace(x_env.min(), x_env.max(), 100)
                # Convert: y = exp(beta) * x^alpha
                y_line = np.exp(beta_env) * (x_line ** alpha_env)
                ax.plot(x_line, y_line, color=colors, linewidth=2.5, 
                       linestyle='-', label=f'Env fit (α={alpha_env:.3f})')
            
            # Global fit line (convert from log space)
            if len(global_row) > 0:
                alpha_global = global_row.iloc[0]['alpha_global']
                beta_global = global_row.iloc[0]['beta_global_log']
                x_line = np.linspace(x_env.min(), x_env.max(), 100)
                # Convert: y = exp(beta) * x^alpha
                y_line = np.exp(beta_global) * (x_line ** alpha_global)
                ax.plot(x_line, y_line, color='gray', linewidth=2, 
                       linestyle='--', alpha=0.8, label=f'Global (α={alpha_global:.3f})')
            
            panel_letter = chr(ord('f') + panel_idx_linear)
            
            # Get GO label
            try:
                plot_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels_for_plots.tsv" if prefix else "go_term_labels_for_plots.tsv")
                plot_labels = pd.read_csv(plot_labels_file, sep='\t')
                plot_labels['category'] = plot_labels['category'].astype(str).str.zfill(7)
                plot_row = plot_labels[plot_labels['category'] == cat_str]
                if len(plot_row) > 0 and 'full_label' in plot_row.columns:
                    cat_label = plot_row.iloc[0]['full_label']
                else:
                    full_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels.tsv" if prefix else "go_term_labels.tsv")
                    full_labels = pd.read_csv(full_labels_file, sep='\t')
                    full_labels['category'] = full_labels['category'].astype(str).str.zfill(7)
                    full_row = full_labels[full_labels['category'] == cat_str]
                    if len(full_row) > 0:
                        cat_label = f"{full_row.iloc[0]['name']} ({full_row.iloc[0]['go_id']})"
                    else:
                        cat_label = f"Category {cat_str}"
            except:
                cat_label = f"Category {cat_str}"
            
            ax.set_xlabel('genes_total', fontsize=10)
            ax.set_ylabel(f'{cat_str} count', fontsize=10)
            ax.set_title(f'Panel {panel_letter}: {env}\n{cat_label} (linear scale)', 
                        fontsize=11, weight='bold')
            ax.legend(fontsize=8, loc='best')
            ax.grid(True, alpha=0.3)
            
            panel_idx_linear += 1
    
    # Hide unused subplots
    for idx in range(panel_idx_linear, len(axes_linear)):
        axes_linear[idx].axis('off')
    
    plt.tight_layout()
    save_figure(fig_linear, "fig1f_to_k_env_scatter_scaling_linear")
    log_message(f"  ✓ Panels 1f-1k (linear scale) completed ({panel_idx_linear} panels)")
    
    # Save panel metadata
    metadata_df = pd.DataFrame(panel_metadata)
    metadata_filename = f"{prefix}fig1_panel_metadata.tsv" if prefix else "fig1_panel_metadata.tsv"
    metadata_df.to_csv(OUTPUT_DIR / metadata_filename, sep='\t', index=False)
    log_message(f"  ✓ {metadata_filename} written")
    
except Exception as e:
    log_message(f"  ✗ ERROR generating Panels 1f-1k: {e}")
    import traceback
    log_message(traceback.format_exc())

# ============================================================================
# Panels: GO category count vs total annotated domains
# ============================================================================

log_message("")
log_message("Generating Panels: GO category count vs total annotated domains...")
try:
    # Calculate total annotated domains per genome (sum of all GO category counts)
    go_cols = [col for col in master_table.columns 
              if col.startswith('0') and len(col) == 7 and col.isdigit()]
    master_table['total_annotated_domains'] = master_table[go_cols].sum(axis=1)
    
    # Use same categories as before
    n_cats_to_plot = min(5, len(selected_cats))
    n_envs_per_cat = 3
    n_panels = n_cats_to_plot * n_envs_per_cat
    
    n_rows = n_cats_to_plot
    n_cols = n_envs_per_cat
    fig_domains, axes_domains = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    if n_rows == 1:
        axes_domains = axes_domains.reshape(1, -1) if n_cols > 1 else [axes_domains]
    axes_domains = axes_domains.flatten()
    
    panel_metadata_domains = []
    panel_idx_domains = 0
    
    for cat_idx, cat in enumerate(selected_cats[:n_cats_to_plot]):
        cat_str = str(cat).zfill(7)
        
        # Get environments that have data for this category
        cat_envs = env_scaling[env_scaling['category'] == cat_str]['environment'].unique()
        envs_to_plot = [e for e in top_envs if e in cat_envs][:n_envs_per_cat]
        
        for env_idx, env in enumerate(envs_to_plot):
            if panel_idx_domains >= n_panels:
                break
            
            ax = axes_domains[panel_idx_domains]
            
            # Get data
            env_data = master_table[master_table['environment'] == env].copy()
            subset = env_data[(env_data['total_annotated_domains'] > 0) & (env_data[cat_str] > 0)].copy()
            
            if len(subset) == 0:
                continue
            
            # Plot in log-log space (log10)
            all_subset = master_table[(master_table['total_annotated_domains'] > 0) & 
                                     (master_table[cat_str] > 0)].copy()
            x_all = np.log10(all_subset['total_annotated_domains'].values)
            y_all = np.log10(all_subset[cat_str].values)
            x_env = np.log10(subset['total_annotated_domains'].values)
            y_env = np.log10(subset[cat_str].values)
            
            # Gray dots: all genomes
            ax.scatter(x_all, y_all, alpha=0.2, s=10, color='lightgray', label='All genomes')
            
            # Colored dots: specific environment
            colors = plt.cm.tab10(cat_idx)
            ax.scatter(x_env, y_env, alpha=0.6, s=30, color=colors, 
                      edgecolors='black', linewidth=0.5, label=env)
            
            # Fit scaling: log(nc) = beta + alpha * log(total_domains)
            # Fit for this env×category
            if len(subset) >= 10:
                x_fit = np.log(subset['total_annotated_domains'].values)
                y_fit = np.log(subset[cat_str].values)
                
                # Simple OLS fit
                x_mean = np.mean(x_fit)
                y_mean = np.mean(y_fit)
                numerator = np.sum((x_fit - x_mean) * (y_fit - y_mean))
                denominator = np.sum((x_fit - x_mean) ** 2)
                if denominator > 0:
                    alpha_fit = numerator / denominator
                    beta_fit = y_mean - alpha_fit * x_mean
                    
                    x_line = np.linspace(x_env.min(), x_env.max(), 100)
                    y_line = beta_fit / np.log(10) + alpha_fit * x_line
                    ax.plot(x_line, y_line, color=colors, linewidth=2.5, 
                           linestyle='-', label=f'Fit (α={alpha_fit:.3f})')
            
            panel_letter = chr(ord('f') + panel_idx_domains)
            
            # Get GO label
            try:
                plot_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels_for_plots.tsv" if prefix else "go_term_labels_for_plots.tsv")
                plot_labels = pd.read_csv(plot_labels_file, sep='\t')
                plot_labels['category'] = plot_labels['category'].astype(str).str.zfill(7)
                plot_row = plot_labels[plot_labels['category'] == cat_str]
                if len(plot_row) > 0 and 'full_label' in plot_row.columns:
                    cat_label = plot_row.iloc[0]['full_label']
                else:
                    full_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels.tsv" if prefix else "go_term_labels.tsv")
                    full_labels = pd.read_csv(full_labels_file, sep='\t')
                    full_labels['category'] = full_labels['category'].astype(str).str.zfill(7)
                    full_row = full_labels[full_labels['category'] == cat_str]
                    if len(full_row) > 0:
                        cat_label = f"{full_row.iloc[0]['name']} ({full_row.iloc[0]['go_id']})"
                    else:
                        cat_label = f"Category {cat_str}"
            except:
                cat_label = f"Category {cat_str}"
            
            ax.set_xlabel('log10(total annotated domains)', fontsize=10)
            ax.set_ylabel(f'log10({cat_str} count)', fontsize=10)
            ax.set_title(f'Panel {panel_letter}: {env}\n{cat_label}\n(vs total domains)', 
                        fontsize=11, weight='bold')
            ax.legend(fontsize=8, loc='best')
            ax.grid(True, alpha=0.3)
            
            panel_metadata_domains.append({
                'panel': panel_letter,
                'category': cat_str,
                'short_label': cat_label,
                'environment': env,
                'notes': f'Category {cat_str} vs total annotated domains in {env}',
                'axes_log_transformed': 'log10',
                'x_axis': 'total_annotated_domains'
            })
            
            panel_idx_domains += 1
    
    # Hide unused subplots
    for idx in range(panel_idx_domains, len(axes_domains)):
        axes_domains[idx].axis('off')
    
    plt.tight_layout()
    save_figure(fig_domains, "fig1_domains_vs_total_domains")
    log_message(f"  ✓ Panels (domains vs total domains, log scale) completed ({panel_idx_domains} panels)")
    
    # Linear scale version
    fig_domains_linear, axes_domains_linear = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 5*n_rows))
    if n_rows == 1:
        axes_domains_linear = axes_domains_linear.reshape(1, -1) if n_cols > 1 else [axes_domains_linear]
    axes_domains_linear = axes_domains_linear.flatten()
    
    panel_idx_domains_linear = 0
    
    for cat_idx, cat in enumerate(selected_cats[:n_cats_to_plot]):
        cat_str = str(cat).zfill(7)
        cat_envs = env_scaling[env_scaling['category'] == cat_str]['environment'].unique()
        envs_to_plot = [e for e in top_envs if e in cat_envs][:n_envs_per_cat]
        
        for env_idx, env in enumerate(envs_to_plot):
            if panel_idx_domains_linear >= n_panels:
                break
            
            ax = axes_domains_linear[panel_idx_domains_linear]
            
            env_data = master_table[master_table['environment'] == env].copy()
            subset = env_data[(env_data['total_annotated_domains'] > 0) & (env_data[cat_str] > 0)].copy()
            
            if len(subset) == 0:
                continue
            
            # Linear scale
            all_subset = master_table[(master_table['total_annotated_domains'] > 0) & 
                                     (master_table[cat_str] > 0)].copy()
            x_all = all_subset['total_annotated_domains'].values
            y_all = all_subset[cat_str].values
            x_env = subset['total_annotated_domains'].values
            y_env = subset[cat_str].values
            
            ax.scatter(x_all, y_all, alpha=0.2, s=10, color='lightgray', label='All genomes')
            colors = plt.cm.tab10(cat_idx)
            ax.scatter(x_env, y_env, alpha=0.6, s=30, color=colors, 
                      edgecolors='black', linewidth=0.5, label=env)
            
            # Fit in linear space
            if len(subset) >= 10:
                x_fit = np.log(subset['total_annotated_domains'].values)
                y_fit = np.log(subset[cat_str].values)
                x_mean = np.mean(x_fit)
                y_mean = np.mean(y_fit)
                numerator = np.sum((x_fit - x_mean) * (y_fit - y_mean))
                denominator = np.sum((x_fit - x_mean) ** 2)
                if denominator > 0:
                    alpha_fit = numerator / denominator
                    beta_fit = y_mean - alpha_fit * x_mean
                    
                    x_line = np.linspace(x_env.min(), x_env.max(), 100)
                    y_line = np.exp(beta_fit) * (x_line ** alpha_fit)
                    ax.plot(x_line, y_line, color=colors, linewidth=2.5, 
                           linestyle='-', label=f'Fit (α={alpha_fit:.3f})')
            
            try:
                plot_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels_for_plots.tsv" if prefix else "go_term_labels_for_plots.tsv")
                plot_labels = pd.read_csv(plot_labels_file, sep='\t')
                plot_labels['category'] = plot_labels['category'].astype(str).str.zfill(7)
                plot_row = plot_labels[plot_labels['category'] == cat_str]
                if len(plot_row) > 0 and 'full_label' in plot_row.columns:
                    cat_label = plot_row.iloc[0]['full_label']
                else:
                    full_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}go_term_labels.tsv" if prefix else "go_term_labels.tsv")
                    full_labels = pd.read_csv(full_labels_file, sep='\t')
                    full_labels['category'] = full_labels['category'].astype(str).str.zfill(7)
                    full_row = full_labels[full_labels['category'] == cat_str]
                    if len(full_row) > 0:
                        cat_label = f"{full_row.iloc[0]['name']} ({full_row.iloc[0]['go_id']})"
                    else:
                        cat_label = f"Category {cat_str}"
            except:
                cat_label = f"Category {cat_str}"
            
            ax.set_xlabel('total annotated domains', fontsize=10)
            ax.set_ylabel(f'{cat_str} count', fontsize=10)
            ax.set_title(f'Panel {panel_letter}: {env}\n{cat_label}\n(vs total domains, linear)', 
                        fontsize=11, weight='bold')
            ax.legend(fontsize=8, loc='best')
            ax.grid(True, alpha=0.3)
            
            panel_idx_domains_linear += 1
    
    # Hide unused subplots
    for idx in range(panel_idx_domains_linear, len(axes_domains_linear)):
        axes_domains_linear[idx].axis('off')
    
    plt.tight_layout()
    save_figure(fig_domains_linear, "fig1_domains_vs_total_domains_linear")
    log_message(f"  ✓ Panels (domains vs total domains, linear scale) completed ({panel_idx_domains_linear} panels)")
    
except Exception as e:
    log_message(f"  ✗ ERROR generating domains vs total domains plots: {e}")
    import traceback
    log_message(traceback.format_exc())

# ============================================================================
# QC and Final Output
# ============================================================================

log_message("")
log_message("=" * 80)
log_message("QC Summary")
log_message("=" * 80)

# Verify figure files
figure_files = [
    "fig1a_Z_exponents_by_category_env",
    "fig1b_Z_offsets_by_category_env",
    "fig1cde_env_exponents_selected_categories",
    "fig1f_to_k_env_scatter_scaling",
    "fig1f_to_k_env_scatter_scaling_linear"
]

for fig_file in figure_files:
    png_path = OUTPUT_DIR / f"{fig_file}.png"
    pdf_path = OUTPUT_DIR / f"{fig_file}.pdf"
    if png_path.exists() and pdf_path.exists():
        png_size = png_path.stat().st_size
        pdf_size = pdf_path.stat().st_size
        if png_size > 10000 and pdf_size > 10000:
            log_message(f"  ✓ {fig_file}: PNG ({png_size/1024:.1f} KB), PDF ({pdf_size/1024:.1f} KB)")
        else:
            log_message(f"  ⚠ {fig_file}: File size suspiciously small")
    else:
        log_message(f"  ✗ {fig_file}: Missing files")

log_message("")
log_message("Writing QC log...")
with open(QC_LOG_FILE, 'w') as f:
    f.write('\n'.join(qc_log))
log_message(f"  ✓ qc_06_figures.log written")

log_message("")
log_message("=" * 80)
log_message("Script 06 completed successfully!")
log_message("=" * 80)

