#!/usr/bin/env python3
"""
Utility functions for prevalence-based GO term filtering.
"""

import pandas as pd
import numpy as np
from pathlib import Path

def get_prevalence_prefix(prevalence_threshold):
    """Get filename prefix for prevalence threshold (e.g., 'prev95_')."""
    if prevalence_threshold is None:
        return ""
    return f"prev{int(prevalence_threshold)}_"

def filter_go_columns_by_prevalence(df, prevalence_threshold, min_genomes=10):
    """
    Filter GO category columns by prevalence threshold.
    
    Args:
        df: DataFrame with GO category columns (7-digit format)
        prevalence_threshold: Float (0-100) or None. If None, return all GO columns.
        min_genomes: Minimum number of genomes required for prevalence calculation
    
    Returns:
        List of GO column names that meet the threshold
    """
    if prevalence_threshold is None:
        # Return all GO columns
        go_cols = [col for col in df.columns 
                  if col.startswith('0') and len(col) == 7 and col.isdigit()]
        return go_cols
    
    # Identify GO columns
    go_cols = [col for col in df.columns 
              if col.startswith('0') and len(col) == 7 and col.isdigit()]
    
    if len(go_cols) == 0:
        return []
    
    # Calculate prevalence for each GO term
    n_genomes = len(df)
    if n_genomes < min_genomes:
        # Not enough genomes, return all
        return go_cols
    
    # Prevalence = percentage of genomes with count > 0
    prevalence_threshold_decimal = prevalence_threshold / 100.0
    min_count = int(np.ceil(n_genomes * prevalence_threshold_decimal))
    
    filtered_cols = []
    for col in go_cols:
        # Count genomes with this term present (count > 0)
        n_present = (df[col] > 0).sum()
        if n_present >= min_count:
            filtered_cols.append(col)
    
    return filtered_cols

def load_prevalence_filtered_terms(prevalence_threshold, base_dir):
    """
    Load GO terms filtered by prevalence threshold.
    
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

