#!/usr/bin/env python3
"""
Utility functions for prevalence-based KEGG term filtering.
Adapted for reactions, KOs, and pathways.
"""

import pandas as pd
import numpy as np
from pathlib import Path

def get_prevalence_prefix(prevalence_threshold):
    """Get filename prefix for prevalence threshold (e.g., 'prev95_')."""
    if prevalence_threshold is None:
        return ""
    return f"prev{int(prevalence_threshold)}_"

def detect_kegg_data_type(df):
    """
    Detect KEGG data type based on column names.
    
    Args:
        df: DataFrame with KEGG category columns
    
    Returns:
        String: 'reactions', 'ko', 'pathway', or None if cannot determine
    """
    # Sample some column names (excluding common metadata columns)
    metadata_cols = {'Genome', 'accession', 'environment', 'genes_total', 'checkm_completeness', 
                     'checkm_contamination', 'env_n_genomes'}
    sample_cols = [col for col in df.columns if col not in metadata_cols][:100]
    
    if not sample_cols:
        return None
    
    # Check for reactions (R + 5 digits)
    reaction_pattern = sum(1 for col in sample_cols if col.startswith('R') and len(col) == 6 and col[1:].isdigit())
    
    # Check for KOs (K + 5 digits)
    ko_pattern = sum(1 for col in sample_cols if col.startswith('K') and len(col) == 6 and col[1:].isdigit())
    
    # Check for pathways (map or rn followed by 5 digits)
    pathway_pattern = sum(1 for col in sample_cols if (col.startswith('map') or col.startswith('rn')) and len(col) >= 7)
    
    # Determine which pattern is most prevalent
    if reaction_pattern > max(ko_pattern, pathway_pattern):
        return 'reactions'
    elif ko_pattern > max(reaction_pattern, pathway_pattern):
        return 'ko'
    elif pathway_pattern > max(reaction_pattern, ko_pattern):
        return 'pathway'
    
    return None

def get_kegg_columns(df, data_type=None):
    """
    Get KEGG category columns based on data type.
    
    Args:
        df: DataFrame with KEGG category columns
        data_type: 'reactions', 'ko', 'pathway', or None (auto-detect)
    
    Returns:
        List of KEGG column names
    """
    if data_type is None:
        data_type = detect_kegg_data_type(df)
    
    if data_type == 'reactions':
        # Reactions: R + 5 digits (e.g., R00130)
        return [col for col in df.columns if col.startswith('R') and len(col) == 6 and col[1:].isdigit()]
    elif data_type == 'ko':
        # KOs: K + 5 digits (e.g., K00005)
        return [col for col in df.columns if col.startswith('K') and len(col) == 6 and col[1:].isdigit()]
    elif data_type == 'pathway':
        # Pathways: ONLY 'map' prefix (reference pathways)
        # Exclude 'rn' prefix (those are reaction networks, duplicates of map pathways)
        return [col for col in df.columns if col.startswith('map') and len(col) >= 7]
    else:
        return []

def filter_kegg_columns_by_prevalence(df, prevalence_threshold, data_type=None, min_genomes=10):
    """
    Filter KEGG category columns by prevalence threshold.
    
    Args:
        df: DataFrame with KEGG category columns
        prevalence_threshold: Float (0-100) or None. If None, return all KEGG columns.
        data_type: 'reactions', 'ko', 'pathway', or None (auto-detect)
        min_genomes: Minimum number of genomes required for prevalence calculation
    
    Returns:
        List of KEGG column names that meet the threshold
    """
    if prevalence_threshold is None:
        # Return all KEGG columns
        return get_kegg_columns(df, data_type)
    
    # Identify KEGG columns
    kegg_cols = get_kegg_columns(df, data_type)
    
    if len(kegg_cols) == 0:
        return []
    
    # Calculate prevalence for each KEGG term
    n_genomes = len(df)
    if n_genomes < min_genomes:
        # Not enough genomes, return all
        return kegg_cols
    
    # Prevalence = percentage of genomes with count > 0
    prevalence_threshold_decimal = prevalence_threshold / 100.0
    min_count = int(np.ceil(n_genomes * prevalence_threshold_decimal))
    
    filtered_cols = []
    for col in kegg_cols:
        # Count genomes with this term present (count > 0)
        n_present = (df[col] > 0).sum()
        if n_present >= min_count:
            filtered_cols.append(col)
    
    return filtered_cols

def get_kegg_input_files(data_type, prevalence_threshold, base_dir):
    """
    Get KEGG input file paths based on data type and prevalence threshold.
    
    Args:
        data_type: 'reactions', 'ko', or 'pathway'
        prevalence_threshold: Float (0-100) or None
        base_dir: Base directory path
    
    Returns:
        Path to the appropriate KEGG counts file
    """
    kegg_dir = base_dir / "results/3.5_KEGG_n_reaction_analyses"
    
    # Base filenames for each data type
    base_files = {
        'reactions': 'all_reactions_counts_table_kegg_reaction',
        'ko': 'all_ko_counts_table_kegg_reaction',
        'pathway': 'all_pathway_counts_table_kegg_reaction'
    }
    
    if data_type not in base_files:
        raise ValueError(f"Invalid data_type: {data_type}. Must be 'reactions', 'ko', or 'pathway'")
    
    base_name = base_files[data_type]
    
    # Check for pre-filtered file if threshold specified
    if prevalence_threshold is not None:
        # Try both prev_XX formats (95 and 99 are available)
        threshold_int = int(prevalence_threshold)
        filtered_file = kegg_dir / f"{base_name}_prev_{threshold_int}.txt"
        if filtered_file.exists():
            return filtered_file
    
    # Return base file (no prevalence filter)
    base_file = kegg_dir / f"{base_name}.txt"
    return base_file

def get_data_type_suffix(data_type):
    """
    Get suffix for output files based on data type.
    
    Args:
        data_type: 'reactions', 'ko', or 'pathway'
    
    Returns:
        String suffix (e.g., '_reactions', '_ko', '_pathway')
    """
    return f"_{data_type}"

# Legacy functions for backward compatibility (now redirect to KEGG versions)
def filter_go_columns_by_prevalence(df, prevalence_threshold, min_genomes=10):
    """Legacy function for GO term filtering (redirects to KEGG version)."""
    return filter_kegg_columns_by_prevalence(df, prevalence_threshold, data_type=None, min_genomes=min_genomes)

def load_prevalence_filtered_terms(prevalence_threshold, base_dir):
    """
    Load GO terms filtered by prevalence threshold.
    (Legacy function - kept for compatibility but not used in KEGG pipeline)
    
    Args:
        prevalence_threshold: Float (0-100) or None
        base_dir: Base directory path
    
    Returns:
        Set of GO term IDs (7-digit format) or None if threshold is None
    """
    if prevalence_threshold is None:
        return None
    
    # Try to load from pre-computed file
    prefix = get_prevalence_prefix(prevalence_threshold)
    filtered_file = base_dir / f"results/3_GO_analyses/{prefix}ubiquitous_terms.txt"
    
    if filtered_file.exists():
        with open(filtered_file, 'r') as f:
            terms = set(line.strip().zfill(7) for line in f if line.strip())
        return terms
    
    # If file doesn't exist, return None (will be computed on-the-fly)
    return None
