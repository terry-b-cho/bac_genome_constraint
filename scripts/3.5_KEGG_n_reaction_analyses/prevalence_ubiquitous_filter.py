#!/usr/bin/env python3
"""
Prevalence Ubiquitous Filter

This script creates filtered versions of count tables based on prevalence thresholds.
For each input file, it creates:
- prev_95: Terms present in ≥95% of genomes
- prev_99: Terms present in ≥99% of genomes

Input files:
- all_ko_counts_table_kegg_reaction.txt
- all_pathway_counts_table_kegg_reaction.txt
- all_reactions_counts_table_kegg_reaction.txt

Output files are saved to the same directory with _prev_95 and _prev_99 suffixes.
"""

import os
import sys
import logging
import argparse
from collections import defaultdict
from datetime import datetime

# ============================================================================
# Configuration
# ============================================================================

BASE_DIR = "/n/scratch/users/b/byc014/github/bac_genome_constraint"
RESULTS_DIR = f"{BASE_DIR}/results/3.5_KEGG_n_reaction_analyses"

# Input files to process
INPUT_FILES = [
    "all_ko_counts_table_kegg_reaction.txt",
    "all_pathway_counts_table_kegg_reaction.txt",
    "all_reactions_counts_table_kegg_reaction.txt"
]

# ============================================================================
# Safety Functions
# ============================================================================

def validate_output_path(file_path):
    """
    Validate that a file path is within RESULTS_DIR to prevent corruption of upstream data.
    Raises AssertionError if path is outside RESULTS_DIR.
    """
    abs_results_dir = os.path.abspath(RESULTS_DIR)
    abs_file_path = os.path.abspath(file_path)
    if not abs_file_path.startswith(abs_results_dir):
        raise AssertionError(
            f"ERROR: Attempted to write outside RESULTS_DIR!\n"
            f"  File: {file_path}\n"
            f"  RESULTS_DIR: {RESULTS_DIR}\n"
            f"  This would corrupt upstream data. All outputs must be in {RESULTS_DIR}/"
        )
    return abs_file_path

def safe_open_write(file_path, mode='w'):
    """
    Safely open a file for writing, ensuring it's within RESULTS_DIR.
    Usage: with safe_open_write(path) as wf: ...
    """
    validated_path = validate_output_path(file_path)
    # Ensure parent directory exists
    os.makedirs(os.path.dirname(validated_path), exist_ok=True)
    return open(validated_path, mode)

# ============================================================================
# Logging Setup
# ============================================================================

def setup_logging(log_file=None):
    """Setup detailed logging to both file and console."""
    if log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        log_file = f"{RESULTS_DIR}/prevalence_filter_{timestamp}.log"
    
    # Validate log file path
    validate_output_path(log_file)
    
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return log_file

# ============================================================================
# Main Processing Functions
# ============================================================================

def calculate_prevalence(count_table_file):
    """
    Calculate prevalence for each term (column) in the count table.
    
    Args:
        count_table_file: Path to the count table file
        
    Returns:
        tuple: (header, data_rows, term_prevalence_dict, total_genomes)
            - header: List of column names (first is "Genome")
            - data_rows: List of lists, each containing genome_id and counts
            - term_prevalence_dict: Dict mapping term_id -> (count_with_term, total_genomes)
            - total_genomes: Number of genomes (data rows)
    """
    logging.info(f"  Reading {os.path.basename(count_table_file)}...")
    
    header = None
    data_rows = []
    term_prevalence = defaultdict(int)  # term_id -> count of genomes with non-zero value
    
    with open(count_table_file, "r") as rf:
        for line_num, line in enumerate(rf, 1):
            fields = line.strip().split("\t")
            
            if line_num == 1:
                # Header row
                header = fields
                if header[0] != "Genome":
                    raise ValueError(f"Expected first column to be 'Genome', got '{header[0]}'")
                continue
            
            # Data row
            if len(fields) != len(header):
                logging.warning(f"  Line {line_num}: Expected {len(header)} columns, got {len(fields)}")
                continue
            
            genome_id = fields[0]
            data_rows.append(fields)
            
            # Count non-zero values for each term
            for i, count_str in enumerate(fields[1:], 1):  # Skip "Genome" column
                term_id = header[i]
                try:
                    count = int(count_str)
                    if count > 0:
                        term_prevalence[term_id] += 1
                except ValueError:
                    logging.warning(f"  Line {line_num}, term {term_id}: Invalid count '{count_str}'")
    
    total_genomes = len(data_rows)
    logging.info(f"  Found {total_genomes} genomes and {len(header) - 1} terms")
    
    return header, data_rows, term_prevalence, total_genomes

def filter_by_prevalence(count_table_file, output_file, prevalence_threshold):
    """
    Filter count table to only include terms with prevalence ≥ threshold.
    
    Args:
        count_table_file: Path to input count table
        output_file: Path to output filtered table
        prevalence_threshold: Minimum prevalence (0.95 for 95%, 0.99 for 99%)
    """
    logging.info(f"  Processing {os.path.basename(count_table_file)} (threshold: {prevalence_threshold*100}%)...")
    
    # Calculate prevalence
    header, data_rows, term_prevalence, total_genomes = calculate_prevalence(count_table_file)
    
    if total_genomes == 0:
        logging.error(f"  ERROR: No data rows found in {count_table_file}")
        return
    
    # Calculate minimum count for threshold
    min_count = int(total_genomes * prevalence_threshold)
    
    # Filter terms that meet threshold
    filtered_terms = []
    for term_id in header[1:]:  # Skip "Genome" column
        count_with_term = term_prevalence.get(term_id, 0)
        if count_with_term >= min_count:
            filtered_terms.append(term_id)
    
    logging.info(f"  Found {len(filtered_terms)} terms with ≥{prevalence_threshold*100}% prevalence (out of {len(header)-1} total)")
    
    # Get indices of filtered terms (including Genome column at index 0)
    term_indices = [0] + [header.index(term_id) for term_id in filtered_terms]
    
    # Write filtered table
    with safe_open_write(output_file, "w") as wf:
        # Write header
        filtered_header = [header[i] for i in term_indices]
        wf.write("\t".join(filtered_header) + "\n")
        
        # Write data rows
        for row in data_rows:
            filtered_row = [row[i] for i in term_indices]
            wf.write("\t".join(filtered_row) + "\n")
    
    logging.info(f"  ✓ Written: {os.path.basename(output_file)}")

def process_all_files():
    """Process all input files and create filtered versions."""
    logging.info("=" * 80)
    logging.info("Prevalence Ubiquitous Filter")
    logging.info("=" * 80)
    
    for input_file in INPUT_FILES:
        input_path = f"{RESULTS_DIR}/{input_file}"
        
        if not os.path.exists(input_path):
            logging.warning(f"  Skipping {input_file}: File not found")
            continue
        
        # Create output file names
        base_name = input_file.replace(".txt", "")
        
        # 95% prevalence filter
        output_95 = f"{RESULTS_DIR}/{base_name}_prev_95.txt"
        filter_by_prevalence(input_path, output_95, 0.95)
        
        # 99% prevalence filter
        output_99 = f"{RESULTS_DIR}/{base_name}_prev_99.txt"
        filter_by_prevalence(input_path, output_99, 0.99)
        
        logging.info("")
    
    logging.info("=" * 80)
    logging.info("SUCCESS: All filtered tables created!")
    logging.info("=" * 80)

# ============================================================================
# Main Function
# ============================================================================

def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(description="Prevalence Ubiquitous Filter")
    parser.add_argument("--log", type=str, default=None,
                       help="Path to log file (default: auto-generated)")
    args = parser.parse_args()
    
    log_file = setup_logging(log_file=args.log)
    
    try:
        process_all_files()
        logging.info(f"Log file: {log_file}")
    except Exception as e:
        logging.error(f"ERROR: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    main()

