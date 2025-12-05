#!/usr/bin/env python3
"""
Extract Metabolism-Related GO Terms

This script extracts all GO terms related to metabolism by traversing the GO ontology DAG
starting from metabolism-related root terms (e.g., GO:0008152 "metabolic process").

Can optionally filter by prevalence threshold.
"""

import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict, deque
from prevalence_utils import get_prevalence_prefix, filter_go_columns_by_prevalence

# Parse arguments
parser = argparse.ArgumentParser(description='Extract metabolism-related GO terms')
parser.add_argument('--prevalence-threshold', type=float, default=None,
                    help='Prevalence threshold (0-100) for filtering (e.g., 95 for 95%%)')
args = parser.parse_args()

# Define paths
BASE_DIR = Path("/n/scratch/users/b/byc014/github/bac_genome_constraint")
GO_OBO_FILE = BASE_DIR / "data/go/go-basic.obo"
OUTPUT_DIR = BASE_DIR / "results/4_statistical_analyses/05_go_labels"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Metabolism-related root GO terms
METABOLISM_ROOT_TERMS = [
    "GO:0008152",  # metabolic process
    "GO:0044238",  # primary metabolic process
    "GO:0044281",  # small-molecule metabolic process
    "GO:0006520",  # amino acid metabolic process
    "GO:0009056",  # catabolic process
    "GO:0009058",  # biosynthetic process
    "GO:0044249",  # cellular biosynthetic process
    "GO:0044255",  # cellular catabolic process
]

# Also include molecular function terms related to metabolism
METABOLISM_MF_ROOTS = [
    "GO:0003824",  # catalytic activity
    "GO:0016491",  # oxidoreductase activity
    "GO:0003674",  # molecular_function (too broad, but we'll filter)
]

def parse_go_ontology(obo_file):
    """Parse GO ontology OBO file and build DAG structure."""
    print(f"Parsing GO ontology from {obo_file}...")
    
    terms = {}
    is_a_relations = defaultdict(list)  # child -> [parents]
    part_of_relations = defaultdict(list)
    namespace_map = {}  # term_id -> namespace
    
    current_term = None
    current_namespace = None
    
    with open(obo_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line.startswith("id: GO:"):
                current_term = line.split(": ")[1]
                terms[current_term] = {
                    'id': current_term,
                    'name': None,
                    'namespace': None,
                    'def': None,
                    'is_obsolete': False
                }
            
            elif line.startswith("name: ") and current_term:
                terms[current_term]['name'] = line.split("name: ", 1)[1]
            
            elif line.startswith("namespace: ") and current_term:
                namespace = line.split("namespace: ", 1)[1]
                terms[current_term]['namespace'] = namespace
                namespace_map[current_term] = namespace
                current_namespace = namespace
            
            elif line.startswith("def: ") and current_term:
                terms[current_term]['def'] = line.split("def: ", 1)[1]
            
            elif line.startswith("is_obsolete: true") and current_term:
                terms[current_term]['is_obsolete'] = True
            
            elif line.startswith("is_a: GO:") and current_term:
                # Format: is_a: GO:XXXXXXX ! name
                parent = line.split("is_a: ")[1].split(" !")[0].strip()
                is_a_relations[current_term].append(parent)
            
            elif line.startswith("relationship: part_of GO:") and current_term:
                # Format: relationship: part_of GO:XXXXXXX ! name
                parent = line.split("part_of GO:")[1].split(" !")[0].strip()
                part_of_relations[current_term].append(parent)
    
    print(f"  Parsed {len(terms)} GO terms")
    print(f"  Found {len(is_a_relations)} terms with 'is_a' relations")
    print(f"  Found {len(part_of_relations)} terms with 'part_of' relations")
    
    return terms, is_a_relations, part_of_relations, namespace_map

def get_all_descendants(root_terms, is_a_relations, part_of_relations, namespace_map, 
                        include_namespaces=None, exclude_obsolete=True):
    """
    Get all descendant terms from root terms by traversing the DAG.
    
    Args:
        root_terms: List of root GO term IDs
        is_a_relations: Dict mapping child -> [parents] for is_a relations
        part_of_relations: Dict mapping child -> [parents] for part_of relations
        namespace_map: Dict mapping term_id -> namespace
        include_namespaces: List of namespaces to include (e.g., ['biological_process', 'molecular_function'])
        exclude_obsolete: Whether to exclude obsolete terms
    
    Returns:
        Set of descendant GO term IDs
    """
    descendants = set()
    queue = deque(root_terms)
    visited = set()
    
    # Reverse relations: parent -> [children]
    parent_to_children = defaultdict(set)
    for child, parents in is_a_relations.items():
        for parent in parents:
            parent_to_children[parent].add(child)
    for child, parents in part_of_relations.items():
        for parent in parents:
            parent_to_children[parent].add(child)
    
    while queue:
        current = queue.popleft()
        if current in visited:
            continue
        visited.add(current)
        
        # Add current term if it's in the desired namespace
        if include_namespaces is None or namespace_map.get(current) in include_namespaces:
            descendants.add(current)
        
        # Add children to queue
        for child in parent_to_children.get(current, []):
            if child not in visited:
                queue.append(child)
    
    print(f"  Found {len(descendants)} descendant terms from {len(root_terms)} root(s)")
    return descendants

def main():
    print("=" * 80)
    print("Extracting Metabolism-Related GO Terms")
    print("=" * 80)
    
    # Parse ontology
    terms, is_a_relations, part_of_relations, namespace_map = parse_go_ontology(GO_OBO_FILE)
    
    # Get all descendants from metabolism root terms
    # Focus on Biological Process and Molecular Function (exclude Cellular Component)
    print("\nExtracting metabolism-related terms from Biological Process...")
    bp_metabolism_terms = get_all_descendants(
        METABOLISM_ROOT_TERMS,
        is_a_relations,
        part_of_relations,
        namespace_map,
        include_namespaces=['biological_process']
    )
    
    print("\nExtracting metabolism-related terms from Molecular Function...")
    mf_metabolism_terms = get_all_descendants(
        METABOLISM_MF_ROOTS,
        is_a_relations,
        part_of_relations,
        namespace_map,
        include_namespaces=['molecular_function']
    )
    
    # Combine and remove duplicates
    all_metabolism_terms = bp_metabolism_terms | mf_metabolism_terms
    
    print(f"\nTotal unique metabolism-related GO terms: {len(all_metabolism_terms)}")
    print(f"  - Biological Process: {len(bp_metabolism_terms)}")
    print(f"  - Molecular Function: {len(mf_metabolism_terms)}")
    print(f"  - Overlap: {len(bp_metabolism_terms & mf_metabolism_terms)}")
    
    # Filter to only terms that exist in our dataset
    print("\nFiltering to terms present in our dataset...")
    try:
        # Load terms based on prevalence threshold
        prefix = get_prevalence_prefix(args.prevalence_threshold)
        if prefix:
            # Try to load prevalence-filtered file
            filtered_file = BASE_DIR / f"results/3_GO_analyses/{prefix}ubiquitous_terms.txt"
            if filtered_file.exists():
                ubiquitous_terms_file = filtered_file
            else:
                # Fallback: compute from master table
                print(f"  ⚠ {filtered_file} not found, computing from master table...")
                master_table_file = BASE_DIR / "results/4_statistical_analyses/02_env_cohorts/master_table_env_filtered.tsv"
                if master_table_file.exists():
                    master_df = pd.read_csv(master_table_file, sep='\t', nrows=1000)  # Sample for speed
                    filtered_cols = filter_go_columns_by_prevalence(master_df, args.prevalence_threshold)
                    dataset_terms = set(filtered_cols)
                else:
                    # Ultimate fallback: use all ubiquitous terms
                    ubiquitous_terms_file = BASE_DIR / "results/3_GO_analyses/ubiquitous_terms.txt"
                    with open(ubiquitous_terms_file, 'r') as f:
                        dataset_terms = set(line.strip().zfill(7) for line in f if line.strip())
            if 'dataset_terms' not in locals():
                with open(ubiquitous_terms_file, 'r') as f:
                    dataset_terms = set(line.strip().zfill(7) for line in f if line.strip())
        else:
            ubiquitous_terms_file = BASE_DIR / "results/3_GO_analyses/ubiquitous_terms.txt"
            with open(ubiquitous_terms_file, 'r') as f:
                dataset_terms = set(line.strip().zfill(7) for line in f if line.strip())
        
        # Convert GO:XXXXXXX format to 7-digit format
        metabolism_terms_7digit = set()
        for term in all_metabolism_terms:
            if term.startswith("GO:"):
                term_7digit = term.replace("GO:", "").zfill(7)
                metabolism_terms_7digit.add(term_7digit)
        
        # Find intersection
        metabolism_in_dataset = metabolism_terms_7digit & dataset_terms
        
        print(f"  Terms in our dataset: {len(dataset_terms)}")
        print(f"  Metabolism terms in dataset: {len(metabolism_in_dataset)}")
        print(f"  Percentage: {100*len(metabolism_in_dataset)/len(dataset_terms):.1f}%")
        
    except FileNotFoundError:
        print("  ⚠ ubiquitous_terms.txt not found, skipping dataset filtering")
        metabolism_in_dataset = metabolism_terms_7digit
    
    # Save results (with prevalence prefix if applicable)
    prefix = get_prevalence_prefix(args.prevalence_threshold)
    base_name = "metabolic_go_terms.txt"
    output_file = OUTPUT_DIR / (f"{prefix}{base_name}" if prefix else base_name)
    with open(output_file, 'w') as f:
        for term in sorted(metabolism_in_dataset):
            f.write(f"{term}\n")
    
    print(f"\n✓ Saved {len(metabolism_in_dataset)} metabolism-related GO terms to:")
    print(f"  {output_file}")
    
    # Also save with GO: prefix for reference
    base_name_go = "metabolic_go_terms_GO_format.txt"
    output_file_go_format = OUTPUT_DIR / (f"{prefix}{base_name_go}" if prefix else base_name_go)
    with open(output_file_go_format, 'w') as f:
        for term in sorted(metabolism_in_dataset):
            go_id = f"GO:{term}"
            f.write(f"{go_id}\n")
    
    print(f"✓ Also saved in GO: format to:")
    print(f"  {output_file_go_format}")
    
    # Create summary
    base_name_summary = "metabolic_go_terms_summary.txt"
    summary_file = OUTPUT_DIR / (f"{prefix}{base_name_summary}" if prefix else base_name_summary)
    with open(summary_file, 'w') as f:
        f.write("Metabolism-Related GO Terms Summary\n")
        f.write("=" * 80 + "\n\n")
        if args.prevalence_threshold is not None:
            f.write(f"Prevalence threshold: {args.prevalence_threshold}%%\n\n")
        f.write(f"Root terms used:\n")
        for root in METABOLISM_ROOT_TERMS + METABOLISM_MF_ROOTS:
            f.write(f"  - {root}\n")
        f.write(f"\nTotal metabolism-related terms found: {len(all_metabolism_terms)}\n")
        f.write(f"Terms present in our dataset: {len(metabolism_in_dataset)}\n")
        f.write(f"Percentage of dataset: {100*len(metabolism_in_dataset)/len(dataset_terms):.1f}%\n")
    
    print(f"\n✓ Summary saved to: {summary_file}")
    print("\n" + "=" * 80)
    print("Extraction complete!")
    print("=" * 80)

if __name__ == "__main__":
    main()

