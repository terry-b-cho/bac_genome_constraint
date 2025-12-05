#!/usr/bin/env python3
"""
Script 04.5: Plot Intermediate Figures

Produces and saves all needed figures/panels/tables for each step of the analyses.
Organized by script output (01-05).
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
    matplotlib.use('Agg')  # Non-interactive backend
    import matplotlib.pyplot as plt
    import seaborn as sns
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("WARNING: matplotlib/seaborn not available. Install with: pip install matplotlib seaborn")

# Parse arguments
parser = argparse.ArgumentParser(description='Generate intermediate figures for all analysis steps')
parser.add_argument('--test-mode', action='store_true',
                    help='Generate only a subset of figures for testing')
parser.add_argument('--script', type=str, choices=['01', '02', '03', '04', 'all'],
                    default='all', help='Which script figures to generate')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for GO term filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
OUTPUT_DIR = BASE_DIR / "results/4_statistical_analyses/04.5_intermediate_figures"
QC_LOG_FILE = OUTPUT_DIR / "qc_04.5_figures.log"

# Create output directories
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
(OUTPUT_DIR / "script_01").mkdir(exist_ok=True)
(OUTPUT_DIR / "script_02").mkdir(exist_ok=True)
(OUTPUT_DIR / "script_03").mkdir(exist_ok=True)
(OUTPUT_DIR / "script_04").mkdir(exist_ok=True)

# Set style (if matplotlib available)
if HAS_MATPLOTLIB:
    sns.set_style("whitegrid")
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['savefig.dpi'] = 300

# Helper function to load metabolic GO terms
def load_metabolic_go_terms():
    """Load list of metabolism-related GO term IDs (7-digit format)."""
    prefix = get_prevalence_prefix(args.prevalence_threshold)
    base_name_metabolic = "metabolic_go_terms.txt"
    metabolic_terms_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels" / (f"{prefix}{base_name_metabolic}" if prefix else base_name_metabolic)
    if not metabolic_terms_file.exists():
        return None
    with open(metabolic_terms_file, 'r') as f:
        metabolic_terms = set(line.strip().zfill(7) for line in f if line.strip())
    return metabolic_terms

# Initialize QC log
qc_log = []

def log_message(msg):
    """Add message to QC log and print to stdout."""
    qc_log.append(msg)
    print(msg)

def save_figure(fig, filename, script_dir):
    """Save figure in both PNG and PDF formats."""
    if not HAS_MATPLOTLIB:
        log_message(f"  ✗ Cannot save {filename}: matplotlib not available")
        return
    prefix = get_prevalence_prefix(args.prevalence_threshold)
    filename_with_prefix = f"{prefix}{filename}" if prefix else filename
    png_path = OUTPUT_DIR / script_dir / f"{filename_with_prefix}.png"
    pdf_path = OUTPUT_DIR / script_dir / f"{filename_with_prefix}.pdf"
    fig.savefig(png_path, bbox_inches='tight', dpi=300)
    fig.savefig(pdf_path, bbox_inches='tight')
    log_message(f"  ✓ Saved {filename_with_prefix}.png and .pdf")
    plt.close(fig)

log_message("=" * 80)
log_message("Script 04.5: Plot Intermediate Figures")
if args.test_mode:
    log_message("  TEST MODE: Generating subset of figures")
if not HAS_MATPLOTLIB:
    log_message("  ✗ ERROR: matplotlib/seaborn required for figure generation")
    log_message("    Install with: pip install matplotlib seaborn")
    sys.exit(1)
log_message("=" * 80)
log_message("")

# ============================================================================
# Script 01 Figures
# ============================================================================

if args.script in ['01', 'all']:
    log_message("=" * 80)
    log_message("Script 01: QC & Basic Dataset Figures")
    log_message("=" * 80)
    log_message("")
    
    # Load data
    log_message("Loading Script 01 data...")
    try:
        master_raw = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/01_master_table/master_table_raw.tsv",
            sep='\t'
        )
        master_hq = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/01_master_table/master_table_high_quality.tsv",
            sep='\t'
        )
        env_counts = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/01_master_table/environment_counts_all.tsv",
            sep='\t'
        )
        log_message(f"  ✓ Loaded data: {len(master_raw)} raw, {len(master_hq)} HQ genomes")
    except Exception as e:
        log_message(f"  ✗ ERROR: Failed to load Script 01 data: {e}")
        if args.script == '01':
            sys.exit(1)
        else:
            log_message("  ⚠ Skipping Script 01 figures")
    
    # Fig 1A: QC filtering flowchart
    log_message("")
    log_message("Generating Fig 1A: QC filtering flowchart...")
    try:
        n_start = 3088  # From QC log: "Number of genomes with GO data: 3088"
        n_hq = len(master_hq)
        
        # Read QC log to get exact counts
        qc_log_file = BASE_DIR / "results/4_statistical_analyses/01_master_table/qc_01_master_table.log"
        lost_completeness = 259
        lost_contamination = 543
        lost_genes = 0
        lost_env = 77
        nan_completeness = 1
        nan_contamination = 353
        
        if qc_log_file.exists():
            with open(qc_log_file, 'r') as f:
                content = f.read()
                # Extract numbers from log with more specific patterns
                import re
                comp_match = re.search(r'Completeness filter.*?Lost: (\d+)', content, re.DOTALL)
                cont_match = re.search(r'Contamination filter.*?Lost: (\d+)', content, re.DOTALL)
                genes_match = re.search(r'Genes filter.*?Lost (\d+)', content)
                env_match = re.search(r'Environment filter.*?Lost: (\d+) genomes.*?\((\d+) null, (\d+) empty, (\d+)', content, re.DOTALL)
                nan_comp_match = re.search(r'NaN values treated as 0: (\d+)', content)
                nan_cont_match = re.search(r'NaN values treated as 100: (\d+)', content)
                
                if comp_match:
                    lost_completeness = int(comp_match.group(1))
                if cont_match:
                    lost_contamination = int(cont_match.group(1))
                if genes_match:
                    lost_genes = int(genes_match.group(1))
                if env_match:
                    lost_env = int(env_match.group(1))
                if nan_comp_match:
                    nan_completeness = int(nan_comp_match.group(1))
                if nan_cont_match:
                    nan_contamination = int(nan_cont_match.group(1))
        
        # Calculate intermediate counts
        n_after_completeness = n_start - lost_completeness
        n_after_contamination = n_after_completeness - lost_contamination
        n_after_genes = n_after_contamination - lost_genes
        n_after_env = n_after_genes - lost_env
        
        fig, ax = plt.subplots(figsize=(10, 8))
        ax.axis('off')
        
        # Create flowchart using text boxes
        y_pos = 0.9
        box_height = 0.1
        box_width = 0.3
        
        # Start
        rect_start = plt.Rectangle((0.35, y_pos - box_height), box_width, box_height,
                                   facecolor='lightblue', edgecolor='black', linewidth=2)
        ax.add_patch(rect_start)
        ax.text(0.5, y_pos - box_height/2, f'Start\n{n_start} genomes\nwith GO data',
                ha='center', va='center', fontsize=12, weight='bold')
        
        # Arrow down
        ax.arrow(0.5, y_pos - box_height, 0, -0.05, head_width=0.02, head_length=0.02, fc='black')
        
        y_pos -= 0.15
        
        # Completeness filter
        rect_comp = plt.Rectangle((0.35, y_pos - box_height), box_width, box_height,
                                   facecolor='lightyellow', edgecolor='black', linewidth=1.5)
        ax.add_patch(rect_comp)
        ax.text(0.5, y_pos - box_height/2, 
                f'Completeness >90\nNaN→0: {nan_completeness}\nLost: {lost_completeness}\nRemaining: {n_after_completeness}',
                ha='center', va='center', fontsize=9)
        
        ax.arrow(0.5, y_pos - box_height, 0, -0.05, head_width=0.02, head_length=0.02, fc='black')
        y_pos -= 0.2
        
        # Contamination filter
        rect_cont = plt.Rectangle((0.35, y_pos - box_height), box_width, box_height,
                                   facecolor='lightyellow', edgecolor='black', linewidth=1.5)
        ax.add_patch(rect_cont)
        ax.text(0.5, y_pos - box_height/2,
                f'Contamination <5\nNaN→100: {nan_contamination}\nLost: {lost_contamination}\nRemaining: {n_after_contamination}',
                ha='center', va='center', fontsize=9)
        
        ax.arrow(0.5, y_pos - box_height, 0, -0.05, head_width=0.02, head_length=0.02, fc='black')
        y_pos -= 0.2
        
        # Genes filter
        rect_genes = plt.Rectangle((0.35, y_pos - box_height), box_width, box_height,
                                    facecolor='lightyellow', edgecolor='black', linewidth=1.5)
        ax.add_patch(rect_genes)
        ax.text(0.5, y_pos - box_height/2,
                f'Genes >0\nLost: {lost_genes}\nRemaining: {n_after_genes}',
                ha='center', va='center', fontsize=9)
        
        ax.arrow(0.5, y_pos - box_height, 0, -0.05, head_width=0.02, head_length=0.02, fc='black')
        y_pos -= 0.2
        
        # Environment filter
        rect_env = plt.Rectangle((0.35, y_pos - box_height), box_width, box_height,
                                 facecolor='lightyellow', edgecolor='black', linewidth=1.5)
        ax.add_patch(rect_env)
        ax.text(0.5, y_pos - box_height/2,
                f'Valid Environment\n(not null/empty/"Unclassified")\nLost: {lost_env}\nRemaining: {n_after_env}',
                ha='center', va='center', fontsize=9)
        
        ax.arrow(0.5, y_pos - box_height, 0, -0.05, head_width=0.02, head_length=0.02, fc='black')
        y_pos -= 0.2
        
        # Final
        rect_final = plt.Rectangle((0.35, y_pos - box_height), box_width, box_height,
                                   facecolor='lightgreen', edgecolor='black', linewidth=2)
        ax.add_patch(rect_final)
        ax.text(0.5, y_pos - box_height/2, 
                f'High Quality Dataset\n{n_after_env} genomes\n11 environments',
                ha='center', va='center', fontsize=12, weight='bold')
        
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_title('QC Filtering Flowchart', fontsize=14, weight='bold', pad=20)
        
        save_figure(fig, "fig_01A_QC_filtering_flowchart", "script_01")
        log_message("  ✓ Fig 1A completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 1A: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 1B: Completeness distribution + cutoff
    log_message("")
    log_message("Generating Fig 1B: Completeness distribution + cutoff...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Histogram of completeness
        ax.hist(master_raw['checkm_completeness'].fillna(0), bins=50, alpha=0.7, 
                label='All genomes', color='lightblue', edgecolor='black')
        ax.hist(master_hq['checkm_completeness'], bins=50, alpha=0.7,
                label='High quality (after QC)', color='lightgreen', edgecolor='black')
        
        # Vertical line at 90%
        ax.axvline(x=90, color='red', linestyle='--', linewidth=2, label='Cutoff (90%)')
        
        ax.set_xlabel('CheckM Completeness (%)', fontsize=12)
        ax.set_ylabel('Number of Genomes', fontsize=12)
        ax.set_title('Distribution of CheckM Completeness', fontsize=14, weight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        
        save_figure(fig, "fig_01B_completeness_distribution", "script_01")
        log_message("  ✓ Fig 1B completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 1B: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 1C: Contamination distribution + cutoff
    log_message("")
    log_message("Generating Fig 1C: Contamination distribution + cutoff...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Histogram of contamination
        ax.hist(master_raw['checkm_contamination'].fillna(100), bins=50, alpha=0.7,
                label='All genomes', color='lightcoral', edgecolor='black')
        ax.hist(master_hq['checkm_contamination'], bins=50, alpha=0.7,
                label='High quality (after QC)', color='lightgreen', edgecolor='black')
        
        # Vertical line at 5%
        ax.axvline(x=5, color='red', linestyle='--', linewidth=2, label='Cutoff (5%)')
        
        ax.set_xlabel('CheckM Contamination (%)', fontsize=12)
        ax.set_ylabel('Number of Genomes', fontsize=12)
        ax.set_title('Distribution of CheckM Contamination', fontsize=14, weight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        
        save_figure(fig, "fig_01C_contamination_distribution", "script_01")
        log_message("  ✓ Fig 1C completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 1C: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 1D: Genome size distribution (global)
    log_message("")
    log_message("Generating Fig 1D: Genome size distribution (global)...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Histogram of genes_total
        ax.hist(master_raw['genes_total'], bins=50, alpha=0.7,
                label='Before QC', color='lightblue', edgecolor='black')
        ax.hist(master_hq['genes_total'], bins=50, alpha=0.7,
                label='After QC', color='lightgreen', edgecolor='black')
        
        ax.set_xlabel('Total Genes (genes_total)', fontsize=12)
        ax.set_ylabel('Number of Genomes', fontsize=12)
        ax.set_title('Distribution of Genome Size (Total Genes)', fontsize=14, weight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3)
        
        save_figure(fig, "fig_01D_genome_size_distribution", "script_01")
        log_message("  ✓ Fig 1D completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 1D: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 2A: Environment counts before and after QC
    log_message("")
    log_message("Generating Fig 2A: Environment counts before and after QC...")
    try:
        # Compute pre-QC environment counts
        env_counts_raw = master_raw.groupby('environment').size().reset_index(name='n_genomes_raw')
        env_counts_merged = pd.merge(env_counts, env_counts_raw, on='environment', how='outer').fillna(0)
        env_counts_merged = env_counts_merged.sort_values('n_genomes', ascending=False)
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        x = np.arange(len(env_counts_merged))
        width = 0.35
        
        bars1 = ax.bar(x - width/2, env_counts_merged['n_genomes_raw'], width,
                       label='Before QC', color='lightblue', edgecolor='black')
        bars2 = ax.bar(x + width/2, env_counts_merged['n_genomes'], width,
                       label='After QC', color='lightgreen', edgecolor='black')
        
        ax.set_xlabel('Environment', fontsize=12)
        ax.set_ylabel('Number of Genomes', fontsize=12)
        ax.set_title('Environment Counts Before and After QC', fontsize=14, weight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(env_counts_merged['environment'], rotation=45, ha='right')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        save_figure(fig, "fig_02A_environment_counts_QC", "script_01")
        log_message("  ✓ Fig 2A completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 2A: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 2D: Genome size distributions by environment (pre Script 02 filter)
    log_message("")
    log_message("Generating Fig 2D: Genome size distributions by environment...")
    try:
        fig, ax = plt.subplots(figsize=(14, 6))
        
        # Box plot
        env_order = env_counts.sort_values('n_genomes', ascending=False)['environment'].tolist()
        data_to_plot = [master_hq[master_hq['environment'] == env]['genes_total'].values 
                        for env in env_order]
        
        bp = ax.boxplot(data_to_plot, tick_labels=env_order, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightblue')
            patch.set_edgecolor('black')
        
        ax.set_xlabel('Environment', fontsize=12)
        ax.set_ylabel('Total Genes (genes_total)', fontsize=12)
        ax.set_title('Genome Size Distribution by Environment (After QC)', fontsize=14, weight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        save_figure(fig, "fig_02D_genome_size_by_environment", "script_01")
        log_message("  ✓ Fig 2D completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 2D: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    log_message("")
    log_message("Script 01 figures completed!")
    log_message("")

# ============================================================================
# Script 02 Figures
# ============================================================================

if args.script in ['02', 'all']:
    log_message("=" * 80)
    log_message("Script 02: Final Cohort Environment Figures")
    log_message("=" * 80)
    log_message("")
    
    # Load data
    log_message("Loading Script 02 data...")
    try:
        valid_envs = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/valid_environments_min20.tsv",
            sep='\t'
        )
        master_filtered = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv",
            sep='\t'
        )
        log_message(f"  ✓ Loaded data: {len(valid_envs)} environments, {len(master_filtered)} genomes")
    except Exception as e:
        log_message(f"  ✗ ERROR: Failed to load Script 02 data: {e}")
        if args.script == '02':
            sys.exit(1)
        else:
            log_message("  ⚠ Skipping Script 02 figures")
    
    # Fig 2B: Final 8 environments (≥20 genomes) barplot
    log_message("")
    log_message("Generating Fig 2B: Final 8 environments barplot...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        valid_envs_sorted = valid_envs.sort_values('n_genomes', ascending=False)
        
        bars = ax.bar(range(len(valid_envs_sorted)), valid_envs_sorted['n_genomes'],
                      color='steelblue', edgecolor='black')
        
        # Add value labels on bars
        for i, (idx, row) in enumerate(valid_envs_sorted.iterrows()):
            ax.text(i, row['n_genomes'] + 10, str(int(row['n_genomes'])),
                   ha='center', va='bottom', fontsize=10, weight='bold')
        
        ax.set_xlabel('Environment', fontsize=12)
        ax.set_ylabel('Number of Genomes', fontsize=12)
        ax.set_title('Final Environments (≥20 genomes)', fontsize=14, weight='bold')
        ax.set_xticks(range(len(valid_envs_sorted)))
        ax.set_xticklabels(valid_envs_sorted['environment'], rotation=45, ha='right')
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        save_figure(fig, "fig_02B_final_environments", "script_02")
        log_message("  ✓ Fig 2B completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 2B: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 2C: Environment contribution to final dataset
    log_message("")
    log_message("Generating Fig 2C: Environment contribution to final dataset...")
    try:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        
        env_counts_final = master_filtered.groupby('environment').size().reset_index(name='n_genomes')
        env_counts_final = env_counts_final.sort_values('n_genomes', ascending=False)
        
        # Bar plot
        bars = ax1.bar(range(len(env_counts_final)), env_counts_final['n_genomes'],
                      color='steelblue', edgecolor='black')
        ax1.set_xlabel('Environment', fontsize=12)
        ax1.set_ylabel('Number of Genomes', fontsize=12)
        ax1.set_title('Genome Counts by Environment', fontsize=12, weight='bold')
        ax1.set_xticks(range(len(env_counts_final)))
        ax1.set_xticklabels(env_counts_final['environment'], rotation=45, ha='right')
        ax1.grid(True, alpha=0.3, axis='y')
        
        # Pie chart
        colors = plt.cm.Set3(range(len(env_counts_final)))
        wedges, texts, autotexts = ax2.pie(env_counts_final['n_genomes'], 
                                           labels=env_counts_final['environment'],
                                           autopct='%1.1f%%',
                                           colors=colors,
                                           startangle=90)
        ax2.set_title('Environment Contribution (%)', fontsize=12, weight='bold')
        
        plt.tight_layout()
        save_figure(fig, "fig_02C_environment_contribution", "script_02")
        log_message("  ✓ Fig 2C completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 2C: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 2D: Genome size distributions by retained environment
    log_message("")
    log_message("Generating Fig 2D: Genome size distributions by retained environment...")
    try:
        fig, ax = plt.subplots(figsize=(14, 6))
        
        env_order = env_counts_final.sort_values('n_genomes', ascending=False)['environment'].tolist()
        data_to_plot = [master_filtered[master_filtered['environment'] == env]['genes_total'].values 
                        for env in env_order]
        
        bp = ax.boxplot(data_to_plot, tick_labels=env_order, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightgreen')
            patch.set_edgecolor('black')
        
        ax.set_xlabel('Environment', fontsize=12)
        ax.set_ylabel('Total Genes (genes_total)', fontsize=12)
        ax.set_title('Genome Size Distribution by Retained Environment', fontsize=14, weight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        save_figure(fig, "fig_02D_genome_size_retained_envs", "script_02")
        log_message("  ✓ Fig 2D completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 2D: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    log_message("")
    log_message("Script 02 figures completed!")
    log_message("")

# ============================================================================
# Script 03 Figures
# ============================================================================

if args.script in ['03', 'all']:
    log_message("=" * 80)
    log_message("Script 03: Global Scaling Figures")
    log_message("=" * 80)
    log_message("")
    
    # Load data
    log_message("Loading Script 03 data...")
    try:
        prefix = get_prevalence_prefix(args.prevalence_threshold)
        base_name_global = "global_scaling_params.tsv"
        global_params_file = BASE_DIR / "results/4_statistical_analyses/03_global_scaling" / (f"{prefix}{base_name_global}" if prefix else base_name_global)
        global_params = pd.read_csv(global_params_file, sep='\t')
        master_filtered = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv",
            sep='\t'
        )
        # Load GO labels
        go_labels = None
        go_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels/go_term_labels.tsv"
        if go_labels_file.exists():
            try:
                go_labels = pd.read_csv(go_labels_file, sep='\t')
                go_labels['category'] = go_labels['category'].astype(str).str.zfill(7)
                log_message(f"  ✓ Loaded GO labels for {len(go_labels)} categories")
            except Exception as e:
                log_message(f"  ⚠ Could not load GO labels: {e}")
        else:
            log_message(f"  ⚠ GO labels file not found, will use category IDs only")
        
        def get_category_label(category_id):
            """Get human-readable label for a category ID."""
            if go_labels is not None:
                cat_str = str(category_id).zfill(7)
                label_row = go_labels[go_labels['category'] == cat_str]
                if len(label_row) > 0:
                    name = label_row.iloc[0]['name']
                    go_id = label_row.iloc[0]['go_id']
                    return f"{name} ({go_id})"
            return f"Category {str(category_id).zfill(7)}"
        
        log_message(f"  ✓ Loaded data: {len(global_params)} categories, {len(master_filtered)} genomes")
    except Exception as e:
        log_message(f"  ✗ ERROR: Failed to load Script 03 data: {e}")
        if args.script == '03':
            sys.exit(1)
        else:
            log_message("  ⚠ Skipping Script 03 figures")
    
    # Fig 4A: Histogram of global exponents (α_global)
    log_message("")
    log_message("Generating Fig 4A: Histogram of global exponents...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        alpha_values = global_params['alpha_global'].values
        
        ax.hist(alpha_values, bins=50, color='steelblue', edgecolor='black', alpha=0.7)
        ax.axvline(x=np.median(alpha_values), color='red', linestyle='--', linewidth=2,
                   label=f'Median = {np.median(alpha_values):.3f}')
        ax.axvline(x=1, color='green', linestyle='--', linewidth=2, label='Linear scaling (α=1)')
        
        ax.set_xlabel('Global Scaling Exponent (α_global)', fontsize=12)
        ax.set_ylabel('Number of Categories', fontsize=12)
        ax.set_title('Distribution of Global Scaling Exponents', fontsize=14, weight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')
        
        save_figure(fig, "fig_04A_global_exponents_histogram", "script_03")
        log_message("  ✓ Fig 4A completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 4A: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 4B: Scatter: α_global vs R²
    log_message("")
    log_message("Generating Fig 4B: α_global vs R²...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.scatter(global_params['alpha_global'], global_params['r_squared'],
                  alpha=0.6, s=30, color='steelblue', edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel('Global Scaling Exponent (α_global)', fontsize=12)
        ax.set_ylabel('R²', fontsize=12)
        ax.set_title('Global Exponent vs Fit Quality (R²)', fontsize=14, weight='bold')
        ax.grid(True, alpha=0.3)
        
        # Add text with summary stats
        high_r2 = (global_params['r_squared'] > 0.5).sum()
        ax.text(0.05, 0.95, f'Categories with R² > 0.5: {high_r2} / {len(global_params)}',
               transform=ax.transAxes, fontsize=10, verticalalignment='top',
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
        
        save_figure(fig, "fig_04B_alpha_vs_r2", "script_03")
        log_message("  ✓ Fig 4B completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 4B: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 4C: Scatter: α_global vs mean GO count
    log_message("")
    log_message("Generating Fig 4C: α_global vs mean GO count...")
    try:
        # Compute mean counts per category
        go_cols = [col for col in master_filtered.columns 
                  if col.startswith('0') and len(col) == 7 and col.isdigit()]
        
        mean_counts = []
        categories = []
        for col in go_cols:
            if col in global_params['category'].values:
                mean_counts.append(master_filtered[col].mean())
                categories.append(col)
        
        mean_counts_df = pd.DataFrame({'category': categories, 'mean_count': mean_counts})
        merged = pd.merge(global_params, mean_counts_df, on='category', how='inner')
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        # Filter out zeros/negatives for log scale
        merged_positive = merged[(merged['alpha_global'] > 0) & (merged['mean_count'] > 0)].copy()
        
        if len(merged_positive) > 0:
            ax.scatter(merged_positive['alpha_global'], merged_positive['mean_count'],
                      alpha=0.6, s=30, color='steelblue', edgecolors='black', linewidth=0.5)
            ax.set_xscale('log')
            ax.set_yscale('log')
        else:
            ax.scatter(merged['alpha_global'], merged['mean_count'],
                      alpha=0.6, s=30, color='steelblue', edgecolors='black', linewidth=0.5)
        
        ax.set_xlabel('Global Scaling Exponent (α_global)', fontsize=12)
        ax.set_ylabel('Mean GO Count (log scale)', fontsize=12)
        ax.set_title('Global Exponent vs Mean GO Count', fontsize=14, weight='bold')
        ax.grid(True, alpha=0.3)
        
        save_figure(fig, "fig_04C_alpha_vs_mean_count", "script_03")
        log_message("  ✓ Fig 4C completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 4C: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 4D: Representative global scaling plots (TOP 20 by Z variance)
    log_message("")
    log_message("Generating Fig 4D: Representative global scaling plots (top 20 by Z variance)...")
    try:
        # Load category Z scores to select top 20 by variance
        prefix = get_prevalence_prefix(args.prevalence_threshold)
        base_name_cat_z = "category_Z_summary.tsv"
        category_z_file = BASE_DIR / "results/4_statistical_analyses/04_env_scaling" / (f"{prefix}{base_name_cat_z}" if prefix else base_name_cat_z)
        try:
            category_z_for_selection = pd.read_csv(category_z_file, sep='\t')
            category_z_for_selection['category'] = category_z_for_selection['category'].astype(str).str.zfill(7)
            # Select top 20 by Z_alpha_category
            top_20_cats = category_z_for_selection.nlargest(20, 'Z_alpha_category')['category'].tolist()
        except:
            # Fallback: select by alpha_global
            global_params['category'] = global_params['category'].astype(str).str.zfill(7)
            top_20_cats = global_params.nlargest(20, 'alpha_global')['category'].tolist()
        
        # Ensure categories are strings
        global_params['category'] = global_params['category'].astype(str).str.zfill(7)
        
        # Create grid: 4 rows × 5 columns for 20 categories
        n_rows = 4
        n_cols = 5
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(25, 20))
        axes = axes.flatten()
        
        for idx, cat in enumerate(top_20_cats):
            if idx >= len(axes):
                break
            ax = axes[idx]
            
            # Get data for this category (ensure cat is string)
            cat_str = str(cat).zfill(7)
            if cat_str not in master_filtered.columns:
                log_message(f"  ⚠ Category {cat_str} not found in master_filtered, skipping")
                continue
            
            subset = master_filtered[(master_filtered['genes_total'] > 0) & 
                                    (master_filtered[cat_str] > 0)].copy()
            
            if len(subset) == 0:
                continue
            
            x = np.log(subset['genes_total'].values)
            y = np.log(subset[cat_str].values)
            
            # Scatter plot
            ax.scatter(x, y, alpha=0.3, s=8, color='steelblue')
            
            # Get global fit parameters
            global_row = global_params[global_params['category'] == cat_str]
            if len(global_row) == 0:
                continue
            global_row = global_row.iloc[0]
            alpha = global_row['alpha_global']
            beta = global_row['beta_global_log']
            
            # Plot fitted line
            x_line = np.linspace(x.min(), x.max(), 100)
            y_line = beta + alpha * x_line
            ax.plot(x_line, y_line, 'r-', linewidth=1.5, label=f'α={alpha:.3f}')
            
            cat_label = get_category_label(cat_str)
            # Truncate if too long
            if len(cat_label) > 50:
                cat_label = cat_label[:47] + '...'
            ax.set_xlabel('log(genes_total)', fontsize=8)
            ax.set_ylabel(f'log({cat_str})', fontsize=8)
            ax.set_title(cat_label, fontsize=9, weight='bold')
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for idx in range(len(top_20_cats), len(axes)):
            axes[idx].axis('off')
        
        plt.tight_layout()
        save_figure(fig, "fig_04D_representative_scaling", "script_03")
        log_message(f"  ✓ Fig 4D (log scale) completed - showing top {len(top_20_cats)} categories")
        
        # Create linear scale version (TOP 20)
        fig_linear, axes_linear = plt.subplots(n_rows, n_cols, figsize=(25, 20))
        axes_linear = axes_linear.flatten()
        
        for idx, cat in enumerate(top_20_cats):
            if idx >= len(axes_linear):
                break
            ax = axes_linear[idx]
            
            cat_str = str(cat).zfill(7)
            if cat_str not in master_filtered.columns:
                continue
            
            subset = master_filtered[(master_filtered['genes_total'] > 0) & 
                                    (master_filtered[cat_str] > 0)].copy()
            
            if len(subset) == 0:
                continue
            
            # Use linear scale (raw values)
            x = subset['genes_total'].values
            y = subset[cat_str].values
            
            # Scatter plot
            ax.scatter(x, y, alpha=0.3, s=8, color='steelblue')
            
            # Get global fit parameters
            global_row = global_params[global_params['category'] == cat_str]
            if len(global_row) == 0:
                continue
            global_row = global_row.iloc[0]
            alpha = global_row['alpha_global']
            beta = global_row['beta_global_log']
            
            # Convert fit from log space to linear: y = exp(beta) * x^alpha
            x_line = np.linspace(x.min(), x.max(), 100)
            y_line = np.exp(beta) * (x_line ** alpha)
            ax.plot(x_line, y_line, 'r-', linewidth=1.5, label=f'α={alpha:.3f}')
            
            cat_label = get_category_label(cat_str)
            # Truncate if too long
            if len(cat_label) > 50:
                cat_label = cat_label[:47] + '...'
            ax.set_xlabel('genes_total', fontsize=8)
            ax.set_ylabel(f'{cat_str} count', fontsize=8)
            ax.set_title(f'{cat_label} (linear)', fontsize=9, weight='bold')
            ax.legend(fontsize=7)
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for idx in range(len(top_20_cats), len(axes_linear)):
            axes_linear[idx].axis('off')
        
        plt.tight_layout()
        save_figure(fig_linear, "fig_04D_representative_scaling_linear", "script_03")
        log_message(f"  ✓ Fig 4D (linear scale) completed - showing top {len(top_20_cats)} categories")
        
        # Generate metabolic version
        metabolic_terms = load_metabolic_go_terms()
        if metabolic_terms:
            log_message("  Generating metabolic version of Fig 4D...")
            try:
                category_z_file = BASE_DIR / "results/4_statistical_analyses/04_env_scaling" / (f"{prefix}{base_name_cat_z}" if prefix else base_name_cat_z)
                category_z_for_selection = pd.read_csv(category_z_file, sep='\t')
                category_z_for_selection['category'] = category_z_for_selection['category'].astype(str).str.zfill(7)
                metabolic_category_z = category_z_for_selection[
                    category_z_for_selection['category'].isin(metabolic_terms)
                ]
                if len(metabolic_category_z) > 0:
                    top_20_metabolic = metabolic_category_z.nlargest(20, 'Z_alpha_category')['category'].tolist()
                    
                    # Log scale metabolic version
                    n_rows_met = 4
                    n_cols_met = 5
                    fig_met, axes_met = plt.subplots(n_rows_met, n_cols_met, figsize=(25, 20))
                    axes_met = axes_met.flatten()
                    
                    for idx, cat in enumerate(top_20_metabolic):
                        if idx >= len(axes_met):
                            break
                        ax = axes_met[idx]
                        cat_str = str(cat).zfill(7)
                        if cat_str not in master_filtered.columns:
                            continue
                        subset = master_filtered[(master_filtered['genes_total'] > 0) & 
                                                (master_filtered[cat_str] > 0)].copy()
                        if len(subset) == 0:
                            continue
                        x = np.log(subset['genes_total'].values)
                        y = np.log(subset[cat_str].values)
                        ax.scatter(x, y, alpha=0.3, s=8, color='steelblue')
                        global_row = global_params[global_params['category'] == cat_str]
                        if len(global_row) == 0:
                            continue
                        global_row = global_row.iloc[0]
                        alpha = global_row['alpha_global']
                        beta = global_row['beta_global_log']
                        x_line = np.linspace(x.min(), x.max(), 100)
                        y_line = beta + alpha * x_line
                        ax.plot(x_line, y_line, 'r-', linewidth=1.5, label=f'α={alpha:.3f}')
                        cat_label = get_category_label(cat_str)
                        if len(cat_label) > 50:
                            cat_label = cat_label[:47] + '...'
                        ax.set_xlabel('log(genes_total)', fontsize=8)
                        ax.set_ylabel(f'log({cat_str})', fontsize=8)
                        ax.set_title(cat_label, fontsize=9, weight='bold')
                        ax.legend(fontsize=7)
                        ax.grid(True, alpha=0.3)
                    
                    for idx in range(len(top_20_metabolic), len(axes_met)):
                        axes_met[idx].axis('off')
                    plt.tight_layout()
                    save_figure(fig_met, "metabolic_fig_04D_representative_scaling", "script_03")
                    log_message(f"  ✓ Metabolic Fig 4D (log scale) completed - {len(top_20_metabolic)} categories")
                    
                    # Linear scale metabolic version
                    fig_met_linear, axes_met_linear = plt.subplots(n_rows_met, n_cols_met, figsize=(25, 20))
                    axes_met_linear = axes_met_linear.flatten()
                    
                    for idx, cat in enumerate(top_20_metabolic):
                        if idx >= len(axes_met_linear):
                            break
                        ax = axes_met_linear[idx]
                        cat_str = str(cat).zfill(7)
                        if cat_str not in master_filtered.columns:
                            continue
                        subset = master_filtered[(master_filtered['genes_total'] > 0) & 
                                                (master_filtered[cat_str] > 0)].copy()
                        if len(subset) == 0:
                            continue
                        x = subset['genes_total'].values
                        y = subset[cat_str].values
                        ax.scatter(x, y, alpha=0.3, s=8, color='steelblue')
                        global_row = global_params[global_params['category'] == cat_str]
                        if len(global_row) == 0:
                            continue
                        global_row = global_row.iloc[0]
                        alpha = global_row['alpha_global']
                        beta = global_row['beta_global_log']
                        x_line = np.linspace(x.min(), x.max(), 100)
                        y_line = np.exp(beta) * (x_line ** alpha)
                        ax.plot(x_line, y_line, 'r-', linewidth=1.5, label=f'α={alpha:.3f}')
                        cat_label = get_category_label(cat_str)
                        if len(cat_label) > 50:
                            cat_label = cat_label[:47] + '...'
                        ax.set_xlabel('genes_total', fontsize=8)
                        ax.set_ylabel(f'{cat_str} count', fontsize=8)
                        ax.set_title(f'{cat_label} (linear)', fontsize=9, weight='bold')
                        ax.legend(fontsize=7)
                        ax.grid(True, alpha=0.3)
                    
                    for idx in range(len(top_20_metabolic), len(axes_met_linear)):
                        axes_met_linear[idx].axis('off')
                    plt.tight_layout()
                    save_figure(fig_met_linear, "metabolic_fig_04D_representative_scaling_linear", "script_03")
                    log_message(f"  ✓ Metabolic Fig 4D (linear scale) completed - {len(top_20_metabolic)} categories")
            except Exception as met_e:
                log_message(f"  ⚠ Error generating metabolic Fig 4D: {met_e}")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 4D: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    log_message("")
    log_message("Script 03 figures completed!")
    log_message("")

# ============================================================================
# Script 04 Figures
# ============================================================================

if args.script in ['04', 'all']:
    log_message("=" * 80)
    log_message("Script 04: Environment-Specific Scaling & Z-Score Figures")
    log_message("=" * 80)
    log_message("")
    
    # Load data
    log_message("Loading Script 04 data...")
    try:
        prefix = get_prevalence_prefix(args.prevalence_threshold)
        base_name_env = "env_scaling_params.tsv"
        base_name_z = "env_vs_global_Z_scores.tsv"
        base_name_cat_z = "category_Z_summary.tsv"
        env_scaling = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/04_env_scaling" / (f"{prefix}{base_name_env}" if prefix else base_name_env),
            sep='\t'
        )
        z_scores = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/04_env_scaling" / (f"{prefix}{base_name_z}" if prefix else base_name_z),
            sep='\t'
        )
        category_z = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/04_env_scaling" / (f"{prefix}{base_name_cat_z}" if prefix else base_name_cat_z),
            sep='\t'
        )
        master_filtered = pd.read_csv(
            BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv",
            sep='\t'
        )
        # Load GO labels
        go_labels = None
        go_labels_file = BASE_DIR / "results/4_statistical_analyses/05_go_labels/go_term_labels.tsv"
        if go_labels_file.exists():
            try:
                go_labels = pd.read_csv(go_labels_file, sep='\t')
                go_labels['category'] = go_labels['category'].astype(str).str.zfill(7)
                log_message(f"  ✓ Loaded GO labels for {len(go_labels)} categories")
            except Exception as e:
                log_message(f"  ⚠ Could not load GO labels: {e}")
        else:
            log_message(f"  ⚠ GO labels file not found, will use category IDs only")
        
        def get_category_label(category_id):
            """Get human-readable label for a category ID."""
            if go_labels is not None:
                cat_str = str(category_id).zfill(7)
                label_row = go_labels[go_labels['category'] == cat_str]
                if len(label_row) > 0:
                    name = label_row.iloc[0]['name']
                    go_id = label_row.iloc[0]['go_id']
                    return f"{name} ({go_id})"
            return f"Category {str(category_id).zfill(7)}"
        
        log_message(f"  ✓ Loaded data: {len(env_scaling)} env×category fits, "
                   f"{len(z_scores)} Z-scores, {len(category_z)} categories")
    except Exception as e:
        log_message(f"  ✗ ERROR: Failed to load Script 04 data: {e}")
        if args.script == '04':
            sys.exit(1)
        else:
            log_message("  ⚠ Skipping Script 04 figures")
    
    # Ensure category columns are strings
    env_scaling['category'] = env_scaling['category'].astype(str).str.zfill(7)
    z_scores['category'] = z_scores['category'].astype(str).str.zfill(7)
    category_z['category'] = category_z['category'].astype(str).str.zfill(7)
    
    # Fig 5A: Heatmap of Z_alpha (categories × environments) with clustering
    log_message("")
    log_message("Generating Fig 5A: Heatmap of Z_alpha with hierarchical clustering...")
    try:
        # Create pivot table
        z_alpha_pivot = z_scores.pivot_table(
            values='Z_alpha', 
            index='category', 
            columns='environment',
            aggfunc='mean'
        )
        
        # Fill NaN with 0 for clustering
        z_alpha_pivot_filled = z_alpha_pivot.fillna(0)
        
        # Select top variable categories (by Z_alpha_category) for visualization
        top_cats = category_z.nlargest(50, 'Z_alpha_category')['category'].tolist()
        z_alpha_pivot_subset = z_alpha_pivot_filled.loc[z_alpha_pivot_filled.index.isin(top_cats)]
        
        # Try clustermap for hierarchical clustering (requires scipy)
        # If scipy not available, use simple distance-based sorting
        try:
            from scipy.cluster.hierarchy import linkage, dendrogram
            from scipy.spatial.distance import pdist, squareform
            
            # Compute distances and cluster
            # Cluster rows (categories)
            row_distances = pdist(z_alpha_pivot_subset.values, metric='euclidean')
            row_linkage = linkage(row_distances, method='ward')
            row_order = dendrogram(row_linkage, no_plot=True)['leaves']
            
            # Cluster columns (environments)
            col_distances = pdist(z_alpha_pivot_subset.values.T, metric='euclidean')
            col_linkage = linkage(col_distances, method='ward')
            col_order = dendrogram(col_linkage, no_plot=True)['leaves']
            
            # Reorder data
            z_alpha_clustered = z_alpha_pivot_subset.iloc[row_order, col_order]
            
            fig, ax = plt.subplots(figsize=(12, 14))
            sns.heatmap(z_alpha_clustered, cmap='RdBu_r', center=0, 
                       cbar_kws={'label': 'Z_alpha'}, ax=ax, 
                       xticklabels=True, yticklabels=False)
            ax.set_xlabel('Environment (clustered)', fontsize=12)
            ax.set_ylabel('Category (top 50 by Z_alpha_category, clustered)', fontsize=12)
            ax.set_title('Z_alpha Heatmap (Categories × Environments)\nHierarchical Clustering (Euclidean, Ward)', 
                        fontsize=14, weight='bold')
            plt.tight_layout()
            save_figure(fig, "fig_05A_Z_alpha_heatmap", "script_04")
            log_message("  ✓ Fig 5A (with clustering) completed")
        except ImportError:
            log_message("  ⚠ scipy not available for clustering, using distance-based sorting")
            # Fallback: simple distance-based sorting
            # Sort rows by mean Z_alpha (descending)
            row_means = z_alpha_pivot_subset.mean(axis=1).sort_values(ascending=False)
            z_alpha_sorted = z_alpha_pivot_subset.loc[row_means.index]
            
            # Sort columns by mean Z_alpha (descending)
            col_means = z_alpha_sorted.mean(axis=0).sort_values(ascending=False)
            z_alpha_sorted = z_alpha_sorted[col_means.index]
            
            fig, ax = plt.subplots(figsize=(10, 12))
            sns.heatmap(z_alpha_sorted, cmap='RdBu_r', center=0, 
                       cbar_kws={'label': 'Z_alpha'}, ax=ax, 
                       xticklabels=True, yticklabels=False)
            ax.set_xlabel('Environment (sorted by mean Z)', fontsize=12)
            ax.set_ylabel('Category (top 50 by Z_alpha_category, sorted)', fontsize=12)
            ax.set_title('Z_alpha Heatmap (Categories × Environments)\nSorted by Mean Z Values', 
                        fontsize=14, weight='bold')
            plt.tight_layout()
            save_figure(fig, "fig_05A_Z_alpha_heatmap", "script_04")
            log_message("  ✓ Fig 5A (sorted fallback) completed")
        except Exception as cluster_error:
            log_message(f"  ⚠ Clustering failed ({cluster_error}), using unsorted heatmap")
            fig, ax = plt.subplots(figsize=(10, 12))
            sns.heatmap(z_alpha_pivot_subset, cmap='RdBu_r', center=0, 
                       cbar_kws={'label': 'Z_alpha'}, ax=ax, 
                       xticklabels=True, yticklabels=False)
            ax.set_xlabel('Environment', fontsize=12)
            ax.set_ylabel('Category (top 50 by Z_alpha_category)', fontsize=12)
            ax.set_title('Z_alpha Heatmap (Categories × Environments)', fontsize=14, weight='bold')
            plt.tight_layout()
            save_figure(fig, "fig_05A_Z_alpha_heatmap", "script_04")
            log_message("  ✓ Fig 5A (unsorted fallback) completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 5A: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 5B: Heatmap of Z_beta with clustering
    log_message("")
    log_message("Generating Fig 5B: Heatmap of Z_beta with hierarchical clustering...")
    try:
        z_beta_pivot = z_scores.pivot_table(
            values='Z_beta',
            index='category',
            columns='environment',
            aggfunc='mean'
        )
        
        # Fill NaN with 0 for clustering
        z_beta_pivot_filled = z_beta_pivot.fillna(0)
        
        top_cats_beta = category_z.nlargest(50, 'Z_beta_category')['category'].tolist()
        z_beta_pivot_subset = z_beta_pivot_filled.loc[z_beta_pivot_filled.index.isin(top_cats_beta)]
        
        # Try hierarchical clustering (requires scipy)
        try:
            from scipy.cluster.hierarchy import linkage, dendrogram
            from scipy.spatial.distance import pdist
            
            # Cluster rows (categories)
            row_distances = pdist(z_beta_pivot_subset.values, metric='euclidean')
            row_linkage = linkage(row_distances, method='ward')
            row_order = dendrogram(row_linkage, no_plot=True)['leaves']
            
            # Cluster columns (environments)
            col_distances = pdist(z_beta_pivot_subset.values.T, metric='euclidean')
            col_linkage = linkage(col_distances, method='ward')
            col_order = dendrogram(col_linkage, no_plot=True)['leaves']
            
            # Reorder data
            z_beta_clustered = z_beta_pivot_subset.iloc[row_order, col_order]
            
            fig, ax = plt.subplots(figsize=(12, 14))
            sns.heatmap(z_beta_clustered, cmap='RdBu_r', center=0,
                       cbar_kws={'label': 'Z_beta'}, ax=ax,
                       xticklabels=True, yticklabels=False)
            ax.set_xlabel('Environment (clustered)', fontsize=12)
            ax.set_ylabel('Category (top 50 by Z_beta_category, clustered)', fontsize=12)
            ax.set_title('Z_beta Heatmap (Categories × Environments)\nHierarchical Clustering (Euclidean, Ward)', 
                        fontsize=14, weight='bold')
            plt.tight_layout()
            save_figure(fig, "fig_05B_Z_beta_heatmap", "script_04")
            log_message("  ✓ Fig 5B (with clustering) completed")
        except ImportError:
            log_message("  ⚠ scipy not available for clustering, using distance-based sorting")
            # Fallback: simple distance-based sorting
            row_means = z_beta_pivot_subset.mean(axis=1).sort_values(ascending=False)
            z_beta_sorted = z_beta_pivot_subset.loc[row_means.index]
            col_means = z_beta_sorted.mean(axis=0).sort_values(ascending=False)
            z_beta_sorted = z_beta_sorted[col_means.index]
            
            fig, ax = plt.subplots(figsize=(10, 12))
            sns.heatmap(z_beta_sorted, cmap='RdBu_r', center=0,
                       cbar_kws={'label': 'Z_beta'}, ax=ax,
                       xticklabels=True, yticklabels=False)
            ax.set_xlabel('Environment (sorted by mean Z)', fontsize=12)
            ax.set_ylabel('Category (top 50 by Z_beta_category, sorted)', fontsize=12)
            ax.set_title('Z_beta Heatmap (Categories × Environments)\nSorted by Mean Z Values', 
                        fontsize=14, weight='bold')
            plt.tight_layout()
            save_figure(fig, "fig_05B_Z_beta_heatmap", "script_04")
            log_message("  ✓ Fig 5B (sorted fallback) completed")
        except Exception as cluster_error:
            log_message(f"  ⚠ Clustering failed ({cluster_error}), using unsorted heatmap")
            fig, ax = plt.subplots(figsize=(10, 12))
            sns.heatmap(z_beta_pivot_subset, cmap='RdBu_r', center=0,
                       cbar_kws={'label': 'Z_beta'}, ax=ax,
                       xticklabels=True, yticklabels=False)
            ax.set_xlabel('Environment', fontsize=12)
            ax.set_ylabel('Category (top 50 by Z_beta_category)', fontsize=12)
            ax.set_title('Z_beta Heatmap (Categories × Environments)', fontsize=14, weight='bold')
            plt.tight_layout()
            save_figure(fig, "fig_05B_Z_beta_heatmap", "script_04")
            log_message("  ✓ Fig 5B (unsorted fallback) completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 5B: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 5C: Distribution of |Z_alpha| per environment
    log_message("")
    log_message("Generating Fig 5C: Distribution of |Z_alpha| per environment...")
    try:
        z_scores['abs_Z_alpha'] = z_scores['Z_alpha'].abs()
        
        fig, ax = plt.subplots(figsize=(12, 6))
        
        envs = sorted(z_scores['environment'].unique())
        data_to_plot = [z_scores[z_scores['environment'] == env]['abs_Z_alpha'].values 
                        for env in envs]
        
        bp = ax.boxplot(data_to_plot, tick_labels=envs, patch_artist=True)
        for patch in bp['boxes']:
            patch.set_facecolor('lightcoral')
            patch.set_edgecolor('black')
        
        ax.axhline(y=2, color='red', linestyle='--', linewidth=2, label='|Z| = 2')
        ax.set_xlabel('Environment', fontsize=12)
        ax.set_ylabel('|Z_alpha|', fontsize=12)
        ax.set_title('Distribution of |Z_alpha| per Environment', fontsize=14, weight='bold')
        ax.tick_params(axis='x', rotation=45)
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        save_figure(fig, "fig_05C_abs_Z_alpha_by_env", "script_04")
        log_message("  ✓ Fig 5C completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 5C: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 5D: Number of categories with |Z_alpha| > 2 per environment
    log_message("")
    log_message("Generating Fig 5D: Categories with |Z_alpha| > 2 per environment...")
    try:
        z_scores['abs_Z_alpha'] = z_scores['Z_alpha'].abs()
        significant = z_scores[z_scores['abs_Z_alpha'] > 2].groupby('environment').size().reset_index(name='n_categories')
        significant = significant.sort_values('n_categories', ascending=False)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        bars = ax.bar(range(len(significant)), significant['n_categories'],
                     color='coral', edgecolor='black')
        
        for i, (idx, row) in enumerate(significant.iterrows()):
            ax.text(i, row['n_categories'] + 2, str(int(row['n_categories'])),
                   ha='center', va='bottom', fontsize=10, weight='bold')
        
        ax.set_xlabel('Environment', fontsize=12)
        ax.set_ylabel('Number of Categories with |Z_alpha| > 2', fontsize=12)
        ax.set_title('Categories with Significant Exponent Deviation', fontsize=14, weight='bold')
        ax.set_xticks(range(len(significant)))
        ax.set_xticklabels(significant['environment'], rotation=45, ha='right')
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        save_figure(fig, "fig_05D_significant_categories_by_env", "script_04")
        log_message("  ✓ Fig 5D completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 5D: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 6A: Histogram of Z_alpha_category
    log_message("")
    log_message("Generating Fig 6A: Histogram of Z_alpha_category...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.hist(category_z['Z_alpha_category'], bins=50, color='steelblue', 
               edgecolor='black', alpha=0.7)
        ax.axvline(x=2, color='red', linestyle='--', linewidth=2, label='Z = 2')
        ax.axvline(x=np.median(category_z['Z_alpha_category']), color='green', 
                  linestyle='--', linewidth=2, 
                  label=f'Median = {np.median(category_z["Z_alpha_category"]):.3f}')
        
        ax.set_xlabel('Z_alpha_category', fontsize=12)
        ax.set_ylabel('Number of Categories', fontsize=12)
        ax.set_title('Distribution of Category-Level Z_alpha (Environment Variation)', 
                    fontsize=14, weight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')
        
        save_figure(fig, "fig_06A_Z_alpha_category_histogram", "script_04")
        log_message("  ✓ Fig 6A completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 6A: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 6B: Histogram of Z_beta_category
    log_message("")
    log_message("Generating Fig 6B: Histogram of Z_beta_category...")
    try:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        ax.hist(category_z['Z_beta_category'], bins=50, color='steelblue',
               edgecolor='black', alpha=0.7)
        ax.axvline(x=2, color='red', linestyle='--', linewidth=2, label='Z = 2')
        ax.axvline(x=np.median(category_z['Z_beta_category']), color='green',
                  linestyle='--', linewidth=2,
                  label=f'Median = {np.median(category_z["Z_beta_category"]):.3f}')
        
        ax.set_xlabel('Z_beta_category', fontsize=12)
        ax.set_ylabel('Number of Categories', fontsize=12)
        ax.set_title('Distribution of Category-Level Z_beta (Environment Variation)',
                    fontsize=14, weight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.3, axis='y')
        
        save_figure(fig, "fig_06B_Z_beta_category_histogram", "script_04")
        log_message("  ✓ Fig 6B completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 6B: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 6C: "Exception categories" barplot (Z_alpha_category > 2)
    log_message("")
    log_message("Generating Fig 6C: Exception categories barplot...")
    try:
        exception_cats = category_z[category_z['Z_alpha_category'] > 2].sort_values(
            'Z_alpha_category', ascending=False).head(20)
        
        # Get labels for exception categories
        exception_labels = []
        for _, row in exception_cats.iterrows():
            cat_str = str(row['category']).zfill(7)
            label = get_category_label(cat_str)
            exception_labels.append(label)
        
        fig, ax = plt.subplots(figsize=(14, 8))
        
        bars = ax.barh(range(len(exception_cats)), exception_cats['Z_alpha_category'],
                      color='coral', edgecolor='black')
        
        ax.set_yticks(range(len(exception_cats)))
        ax.set_yticklabels(exception_labels, fontsize=9)
        ax.set_xlabel('Z_alpha_category', fontsize=12)
        ax.set_ylabel('Category', fontsize=12)
        ax.set_title('Top 20 Exception Categories (Z_alpha_category > 2)', 
                    fontsize=14, weight='bold')
        ax.axvline(x=2, color='red', linestyle='--', linewidth=2)
        ax.grid(True, alpha=0.3, axis='x')
        
        plt.tight_layout()
        save_figure(fig, "fig_06C_exception_categories", "script_04")
        log_message("  ✓ Fig 6C completed")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 6C: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    # Fig 7: Environment-stratified scaling scatterplots (TOP 20)
    log_message("")
    log_message("Generating Fig 7: Environment-stratified scaling scatterplots (top 20 categories)...")
    try:
        # Select top 20 variable categories by Z_alpha_category
        top_var_cats = category_z.nlargest(20, 'Z_alpha_category')['category'].tolist()
        
        # Also generate metabolic version
        metabolic_terms = load_metabolic_go_terms()
        if metabolic_terms:
            metabolic_category_z = category_z[category_z['category'].isin(metabolic_terms)]
            if len(metabolic_category_z) > 0:
                top_20_metabolic = metabolic_category_z.nlargest(20, 'Z_alpha_category')['category'].tolist()
                log_message(f"  Generating metabolic version with {len(top_20_metabolic)} categories...")
        
        # Create grid: 4 rows × 5 columns for 20 categories
        n_rows = 4
        n_cols = 5
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 16))
        axes = axes.flatten()
        
        for idx, cat in enumerate(top_var_cats):
            if idx >= len(axes):
                break
            ax = axes[idx]
            cat_str = str(cat).zfill(7)
            
            # Get environments for this category
            envs_for_cat = sorted(z_scores[z_scores['category'] == cat_str]['environment'].unique())
            
            colors = plt.cm.tab10(np.linspace(0, 1, len(envs_for_cat)))
            
            for env_idx, env in enumerate(envs_for_cat):
                # Get data
                env_data = master_filtered[master_filtered['environment'] == env].copy()
                subset = env_data[(env_data['genes_total'] > 0) & (env_data[cat_str] > 0)].copy()
                
                if len(subset) == 0:
                    continue
                
                x = np.log10(subset['genes_total'].values)
                y = np.log10(subset[cat_str].values)
                
                # Scatter
                ax.scatter(x, y, alpha=0.4, s=15, color=colors[env_idx], 
                          label=env, edgecolors='black', linewidth=0.2)
                
                # Get env-specific fit
                env_fit = env_scaling[(env_scaling['category'] == cat_str) & 
                                     (env_scaling['environment'] == env)]
                if len(env_fit) > 0:
                    alpha_env = env_fit.iloc[0]['alpha_env']
                    beta_env = env_fit.iloc[0]['beta_env_log']
                    x_line = np.linspace(x.min(), x.max(), 100)
                    y_line = beta_env / np.log(10) + alpha_env * x_line  # Convert from ln to log10
                    ax.plot(x_line, y_line, color=colors[env_idx], linewidth=1.5, 
                           linestyle='--', alpha=0.7)
            
            cat_label = get_category_label(cat_str)
            # Truncate label if too long for subplot
            if len(cat_label) > 60:
                cat_label = cat_label[:57] + '...'
            ax.set_xlabel('log10(genes_total)', fontsize=8)
            ax.set_ylabel(f'log10({cat_str})', fontsize=8)
            ax.set_title(cat_label, fontsize=9, weight='bold')
            ax.legend(fontsize=6, ncol=2, loc='best')
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots if any
        for idx in range(len(top_var_cats), len(axes)):
            axes[idx].axis('off')
        
        plt.tight_layout()
        save_figure(fig, "fig_07_env_stratified_scaling", "script_04")
        log_message(f"  ✓ Fig 7 (log scale) completed - showing top {len(top_var_cats)} categories")
        
        # Create linear scale version (TOP 20)
        fig_linear, axes_linear = plt.subplots(n_rows, n_cols, figsize=(20, 16))
        axes_linear = axes_linear.flatten()
        
        for idx, cat in enumerate(top_var_cats):
            if idx >= len(axes_linear):
                break
            ax = axes_linear[idx]
            cat_str = str(cat).zfill(7)
            
            # Get environments for this category
            envs_for_cat = sorted(z_scores[z_scores['category'] == cat_str]['environment'].unique())
            
            colors = plt.cm.tab10(np.linspace(0, 1, len(envs_for_cat)))
            
            for env_idx, env in enumerate(envs_for_cat):
                # Get data
                env_data = master_filtered[master_filtered['environment'] == env].copy()
                subset = env_data[(env_data['genes_total'] > 0) & (env_data[cat_str] > 0)].copy()
                
                if len(subset) == 0:
                    continue
                
                # Use linear scale (raw values)
                x = subset['genes_total'].values
                y = subset[cat_str].values
                
                # Scatter
                ax.scatter(x, y, alpha=0.4, s=20, color=colors[env_idx], 
                          label=env, edgecolors='black', linewidth=0.3)
                
                # Get env-specific fit
                env_fit = env_scaling[(env_scaling['category'] == cat_str) & 
                                     (env_scaling['environment'] == env)]
                if len(env_fit) > 0:
                    alpha_env = env_fit.iloc[0]['alpha_env']
                    beta_env = env_fit.iloc[0]['beta_env_log']
                    x_line = np.linspace(x.min(), x.max(), 100)
                    # Convert from log space: y = exp(beta) * x^alpha
                    y_line = np.exp(beta_env) * (x_line ** alpha_env)
                    ax.plot(x_line, y_line, color=colors[env_idx], linewidth=2, 
                           linestyle='--', alpha=0.8)
            
            cat_label = get_category_label(cat_str)
            # Truncate label if too long
            if len(cat_label) > 60:
                cat_label = cat_label[:57] + '...'
            ax.set_xlabel('genes_total', fontsize=8)
            ax.set_ylabel(f'{cat_str} count', fontsize=8)
            ax.set_title(f'{cat_label} (linear)', fontsize=9, weight='bold')
            ax.legend(fontsize=6, ncol=2, loc='best')
            ax.grid(True, alpha=0.3)
        
        # Hide unused subplots
        for idx in range(len(top_var_cats), len(axes_linear)):
            axes_linear[idx].axis('off')
        
        plt.tight_layout()
        save_figure(fig_linear, "fig_07_env_stratified_scaling_linear", "script_04")
        log_message(f"  ✓ Fig 7 (linear scale) completed - showing top {len(top_var_cats)} categories")
        
        # Generate metabolic version if available
        if metabolic_terms and len(metabolic_category_z) > 0:
            # Log scale metabolic version
            fig_met, axes_met = plt.subplots(n_rows, n_cols, figsize=(20, 16))
            axes_met = axes_met.flatten()
            
            for idx, cat in enumerate(top_20_metabolic):
                if idx >= len(axes_met):
                    break
                ax = axes_met[idx]
                cat_str = str(cat).zfill(7)
                envs_for_cat = sorted(z_scores[z_scores['category'] == cat_str]['environment'].unique())
                colors = plt.cm.tab10(np.linspace(0, 1, len(envs_for_cat)))
                
                for env_idx, env in enumerate(envs_for_cat):
                    env_data = master_filtered[master_filtered['environment'] == env].copy()
                    subset = env_data[(env_data['genes_total'] > 0) & (env_data[cat_str] > 0)].copy()
                    if len(subset) == 0:
                        continue
                    x = np.log10(subset['genes_total'].values)
                    y = np.log10(subset[cat_str].values)
                    ax.scatter(x, y, alpha=0.4, s=15, color=colors[env_idx], 
                              label=env, edgecolors='black', linewidth=0.2)
                    env_fit = env_scaling[(env_scaling['category'] == cat_str) & 
                                         (env_scaling['environment'] == env)]
                    if len(env_fit) > 0:
                        alpha_env = env_fit.iloc[0]['alpha_env']
                        beta_env = env_fit.iloc[0]['beta_env_log']
                        x_line = np.linspace(x.min(), x.max(), 100)
                        y_line = beta_env / np.log(10) + alpha_env * x_line
                        ax.plot(x_line, y_line, color=colors[env_idx], linewidth=1.5, 
                               linestyle='--', alpha=0.7)
                
                cat_label = get_category_label(cat_str)
                if len(cat_label) > 60:
                    cat_label = cat_label[:57] + '...'
                ax.set_xlabel('log10(genes_total)', fontsize=8)
                ax.set_ylabel(f'log10({cat_str})', fontsize=8)
                ax.set_title(cat_label, fontsize=9, weight='bold')
                ax.legend(fontsize=6, ncol=2, loc='best')
                ax.grid(True, alpha=0.3)
            
            for idx in range(len(top_20_metabolic), len(axes_met)):
                axes_met[idx].axis('off')
            plt.tight_layout()
            save_figure(fig_met, "metabolic_fig_07_env_stratified_scaling", "script_04")
            log_message(f"  ✓ Metabolic Fig 7 (log scale) completed - {len(top_20_metabolic)} categories")
            
            # Linear scale metabolic version
            fig_met_linear, axes_met_linear = plt.subplots(n_rows, n_cols, figsize=(20, 16))
            axes_met_linear = axes_met_linear.flatten()
            
            for idx, cat in enumerate(top_20_metabolic):
                if idx >= len(axes_met_linear):
                    break
                ax = axes_met_linear[idx]
                cat_str = str(cat).zfill(7)
                envs_for_cat = sorted(z_scores[z_scores['category'] == cat_str]['environment'].unique())
                colors = plt.cm.tab10(np.linspace(0, 1, len(envs_for_cat)))
                
                for env_idx, env in enumerate(envs_for_cat):
                    env_data = master_filtered[master_filtered['environment'] == env].copy()
                    subset = env_data[(env_data['genes_total'] > 0) & (env_data[cat_str] > 0)].copy()
                    if len(subset) == 0:
                        continue
                    x = subset['genes_total'].values
                    y = subset[cat_str].values
                    ax.scatter(x, y, alpha=0.4, s=15, color=colors[env_idx], 
                              label=env, edgecolors='black', linewidth=0.2)
                    env_fit = env_scaling[(env_scaling['category'] == cat_str) & 
                                         (env_scaling['environment'] == env)]
                    if len(env_fit) > 0:
                        alpha_env = env_fit.iloc[0]['alpha_env']
                        beta_env = env_fit.iloc[0]['beta_env_log']
                        x_line = np.linspace(x.min(), x.max(), 100)
                        y_line = np.exp(beta_env) * (x_line ** alpha_env)
                        ax.plot(x_line, y_line, color=colors[env_idx], linewidth=1.5, 
                               linestyle='--', alpha=0.7)
                
                cat_label = get_category_label(cat_str)
                if len(cat_label) > 60:
                    cat_label = cat_label[:57] + '...'
                ax.set_xlabel('genes_total', fontsize=8)
                ax.set_ylabel(f'{cat_str} count', fontsize=8)
                ax.set_title(f'{cat_label} (linear)', fontsize=9, weight='bold')
                ax.legend(fontsize=6, ncol=2, loc='best')
                ax.grid(True, alpha=0.3)
            
            for idx in range(len(top_20_metabolic), len(axes_met_linear)):
                axes_met_linear[idx].axis('off')
            plt.tight_layout()
            save_figure(fig_met_linear, "metabolic_fig_07_env_stratified_scaling_linear", "script_04")
            log_message(f"  ✓ Metabolic Fig 7 (linear scale) completed - {len(top_20_metabolic)} categories")
        
    except Exception as e:
        log_message(f"  ✗ ERROR generating Fig 7: {e}")
        import traceback
        log_message(traceback.format_exc())
    
    log_message("")
    log_message("Script 04 figures completed!")
    log_message("")

# ============================================================================
# Final QC Log
# ============================================================================

log_message("")
log_message("=" * 80)
log_message("Writing QC log...")
prefix = get_prevalence_prefix(args.prevalence_threshold)
base_name_qc = "qc_04.5_figures.log"
qc_log_filename = f"{prefix}{base_name_qc}" if prefix else base_name_qc
qc_log_path = OUTPUT_DIR / qc_log_filename
with open(qc_log_path, 'w') as f:
    f.write('\n'.join(qc_log))
log_message(f"  ✓ {qc_log_filename} written")
log_message("")
log_message("=" * 80)
log_message("Script 04.5 completed successfully!")
log_message("=" * 80)

