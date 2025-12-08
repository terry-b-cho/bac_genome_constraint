#!/usr/bin/env python3
"""
KEGG Reaction and KO Analysis from GO Terms

This script builds count matrices for:
- Reaction terms (KEGG reaction IDs)
- KO terms (KEGG ortholog IDs)
- Pathway terms (KEGG pathway IDs)

based on GO term counts from individual genome files, using KEGG REST API to map
GO → reaction → KO/pathway.

ALL outputs are saved to: results/3.5_KEGG_n_reaction_analyses/
"""

import os
import sys
import glob
import time
import logging
import argparse
import urllib.request
import urllib.parse
from collections import defaultdict
from datetime import datetime
from pathlib import Path

# ============================================================================
# Configuration
# ============================================================================

BASE_DIR = "/n/scratch/users/b/byc014/github/bac_genome_constraint"
RESULTS_DIR = f"{BASE_DIR}/results/3.5_KEGG_n_reaction_analyses"
KEGG_MAPPING_FILE = f"{BASE_DIR}/data/kegg/kegg_reaction2go.txt"
GO_COUNTS_DIR = f"{BASE_DIR}/results/3_GO_analyses/GO_lists_with_counts"

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

def setup_logging(log_file=None, test_mode=False):
    """Setup detailed logging to both file and console."""
    if log_file is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        mode_suffix = "_test" if test_mode else ""
        log_file = f"{RESULTS_DIR}/kegg_go_analyses{mode_suffix}_{timestamp}.log"
    
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
# Phase 1: Parse GO → Reaction Mapping
# ============================================================================

def parse_kegg_reaction2go(mapping_file):
    """
    Parse kegg_reaction2go.txt to build GO → Reaction mapping.
    
    File format: KEGG_REACTION:R00004 > GO:description ; GO:0004427
    
    Returns:
        go_to_reactions: dict mapping GO term (7-digit format) to set of reaction IDs
        reaction_to_gos: dict mapping reaction ID to set of GO terms
    """
    logging.info(f"Phase 1: Parsing GO → Reaction mapping from {mapping_file}")
    go_to_reactions = defaultdict(set)
    reaction_to_gos = defaultdict(set)
    
    line_count = 0
    parsed_count = 0
    
    with open(mapping_file, "r") as rf:
        for line in rf:
            line = line.strip()
            line_count += 1
            
            # Skip comments and empty lines
            if not line or line.startswith("!"):
                continue
            
            # Parse format: KEGG_REACTION:R00004 > GO:description ; GO:0004427
            if ">" in line:
                parts = line.split(">")
                if len(parts) == 2:
                    # Extract reaction ID (e.g., R00004 from KEGG_REACTION:R00004)
                    reaction_part = parts[0].strip()
                    if "KEGG_REACTION:" in reaction_part:
                        reaction_id = reaction_part.split(":")[1].strip()
                        
                        # Extract GO terms from right side
                        go_part = parts[1].strip()
                        # Find all GO:XXXXXXX patterns
                        go_terms = []
                        for item in go_part.split(";"):
                            item = item.strip()
                            if item.startswith("GO:"):
                                go_term = item[3:]  # Remove "GO:" prefix, keep 7-digit ID
                                go_terms.append(go_term)
                        
                        # Build bidirectional mapping
                        for go_term in go_terms:
                            go_to_reactions[go_term].add(reaction_id)
                            reaction_to_gos[reaction_id].add(go_term)
                        parsed_count += 1
    
    logging.info(f"  Parsed {parsed_count} mappings from {line_count} lines")
    logging.info(f"  Found {len(go_to_reactions)} GO terms mapping to {len(reaction_to_gos)} reactions")
    
    # Save intermediate mappings
    logging.info("  Saving intermediate mapping files...")
    with safe_open_write(f"{RESULTS_DIR}/intermediate_mappings/go_to_reactions.tsv", "w") as wf:
        wf.write("GO_term\tReaction_IDs\n")
        for go_term, reactions in sorted(go_to_reactions.items()):
            wf.write(f"{go_term}\t{','.join(sorted(reactions))}\n")
    
    with safe_open_write(f"{RESULTS_DIR}/intermediate_mappings/reaction_to_gos.tsv", "w") as wf:
        wf.write("Reaction_ID\tGO_terms\n")
        for reaction_id, gos in sorted(reaction_to_gos.items()):
            wf.write(f"{reaction_id}\t{','.join(sorted(gos))}\n")
    
    logging.info("  ✓ Phase 1 complete: GO → Reaction mapping built and saved")
    return dict(go_to_reactions), dict(reaction_to_gos)

# ============================================================================
# Phase 2: Load GO Counts from Individual Genome Files
# ============================================================================

def load_go_counts_from_genome_files(go_counts_dir, test_mode=False, max_genomes=None):
    """
    Load GO term counts from individual genome files.
    
    Args:
        go_counts_dir: Directory containing GCF_*.txt files
        test_mode: If True, only process subset
        max_genomes: Maximum number of genomes to process (for testing)
    
    Returns:
        dict mapping genome_id -> dict mapping GO_term -> count
    """
    logging.info(f"Phase 2: Loading GO counts from {go_counts_dir}")
    
    genome_files = sorted(glob.glob(f"{go_counts_dir}/*.txt"))
    total_files = len(genome_files)
    
    if test_mode and max_genomes:
        genome_files = genome_files[:max_genomes]
        logging.info(f"  TEST MODE: Processing first {len(genome_files)} of {total_files} genomes")
    else:
        logging.info(f"  Found {total_files} genome files")
    
    go_counts_per_genome = {}
    missing_files = []
    empty_files = []
    
    for i, genome_file in enumerate(genome_files, 1):
        genome_id = os.path.basename(genome_file).replace(".txt", "")
        
        if not os.path.exists(genome_file):
            missing_files.append(genome_id)
            continue
        
        go_counts = {}
        with open(genome_file, "r") as rf:
            for line in rf:
                fields = line.strip().split("\t")
                if len(fields) == 2:
                    go_term = fields[0]
                    try:
                        count = int(fields[1])
                        go_counts[go_term] = count
                    except ValueError:
                        logging.warning(f"  Invalid count in {genome_id}: {line.strip()}")
        
        if len(go_counts) == 0:
            empty_files.append(genome_id)
        else:
            go_counts_per_genome[genome_id] = go_counts
        
        if i % 100 == 0:
            logging.info(f"  Processed {i}/{len(genome_files)} genomes...")
    
    logging.info(f"  ✓ Loaded GO counts for {len(go_counts_per_genome)} genomes")
    if missing_files:
        logging.warning(f"  Missing files: {len(missing_files)}")
        with safe_open_write(f"{RESULTS_DIR}/missing_genomes_go.txt", "w") as wf:
            for gid in missing_files:
                wf.write(f"{gid}\n")
    if empty_files:
        logging.warning(f"  Empty files: {len(empty_files)}")
    
    # Get all unique GO terms
    all_go_terms = set()
    for go_counts in go_counts_per_genome.values():
        all_go_terms.update(go_counts.keys())
    logging.info(f"  Found {len(all_go_terms)} unique GO terms across all genomes")
    
    logging.info("  ✓ Phase 2 complete: GO counts loaded")
    return go_counts_per_genome

# ============================================================================
# Phase 3: KEGG API Queries
# ============================================================================

def batch_kegg_query(operation, db, ids, batch_size=100, cache_dir=None):
    """
    Make batch KEGG REST API queries with caching.
    
    Args:
        operation: KEGG operation (e.g., 'link', 'get')
        db: Target database (e.g., 'ko', 'pathway')
        ids: List of IDs to query
        batch_size: Number of IDs per batch (KEGG has limits)
        cache_dir: Directory to cache API responses
    
    Returns:
        List of response lines
    """
    all_results = []
    base_url = "https://rest.kegg.jp"
    total_batches = (len(ids) + batch_size - 1) // batch_size
    
    logging.info(f"  Querying KEGG API: {operation}/{db} for {len(ids)} IDs ({total_batches} batches)")
    
    # Process in batches
    for i in range(0, len(ids), batch_size):
        batch = ids[i:i+batch_size]
        batch_num = i // batch_size + 1
        batch_str = "+".join(batch)
        
        # Check cache
        cache_file = None
        if cache_dir:
            cache_key = f"{operation}_{db}_{hash(batch_str)}.txt"
            cache_file = os.path.join(cache_dir, cache_key)
            if os.path.exists(cache_file):
                logging.debug(f"    Batch {batch_num}/{total_batches}: Using cache")
                with open(cache_file, "r") as cf:
                    cached_results = cf.read().strip().split('\n')
                    if cached_results and cached_results[0]:  # Non-empty
                        all_results.extend(cached_results)
                        continue
        
        # Build URL
        if operation == "link":
            url = f"{base_url}/link/{db}/{batch_str}"
        elif operation == "get":
            url = f"{base_url}/get/{batch_str}"
        else:
            raise ValueError(f"Unknown operation: {operation}")
        
        try:
            logging.debug(f"    Batch {batch_num}/{total_batches}: Querying {len(batch)} IDs...")
            with urllib.request.urlopen(url, timeout=30) as response:
                result = response.read().decode('utf-8')
                if result.strip():
                    results = result.strip().split('\n')
                    all_results.extend(results)
                    
                    # Cache result
                    if cache_file:
                        with safe_open_write(cache_file, "w") as cf:
                            cf.write(result)
            
            # Be polite to KEGG API
            time.sleep(0.1)
            
            if batch_num % 10 == 0:
                logging.info(f"    Progress: {batch_num}/{total_batches} batches completed")
                
        except Exception as e:
            logging.error(f"    Error querying batch {batch_num}/{total_batches}: {e}")
            continue
    
    logging.info(f"  ✓ Retrieved {len(all_results)} results from KEGG API")
    return all_results

def get_reaction_to_ko_mappings(all_reaction_ids, cache_dir=None):
    """Get Reaction → KO mappings from KEGG API."""
    logging.info("Phase 3A: Getting Reaction → KO mappings")
    
    # Use "rn:" prefix for API calls
    reaction_ids_with_prefix = [f"rn:{r}" for r in all_reaction_ids]
    results = batch_kegg_query("link", "ko", reaction_ids_with_prefix, batch_size=100, cache_dir=cache_dir)
    
    # Parse results: format is "rn:R00004\tko:K00001"
    reaction_to_kos = defaultdict(set)
    for line in results:
        if "\t" in line:
            parts = line.split("\t")
            if len(parts) == 2:
                reaction_part = parts[0].strip()
                ko_part = parts[1].strip()
                # Extract IDs (remove prefixes like "rn:" and "ko:")
                if ":" in reaction_part and ":" in ko_part:
                    reaction_id = reaction_part.split(":")[1]
                    ko_id = ko_part.split(":")[1]
                    reaction_to_kos[reaction_id].add(ko_id)
    
    logging.info(f"  Mapped {len(reaction_to_kos)} reactions to KO terms")
    
    # Save mapping
    with safe_open_write(f"{RESULTS_DIR}/intermediate_mappings/reaction_to_kos.tsv", "w") as wf:
        wf.write("Reaction_ID\tKO_IDs\n")
        for reaction_id, kos in sorted(reaction_to_kos.items()):
            wf.write(f"{reaction_id}\t{','.join(sorted(kos))}\n")
    
    logging.info("  ✓ Phase 3A complete: Reaction → KO mappings saved")
    return dict(reaction_to_kos)

def get_reaction_to_pathway_mappings(all_reaction_ids, cache_dir=None):
    """Get Reaction → Pathway mappings from KEGG API."""
    logging.info("Phase 3B: Getting Reaction → Pathway mappings")
    
    # Use "rn:" prefix for API calls
    reaction_ids_with_prefix = [f"rn:{r}" for r in all_reaction_ids]
    results = batch_kegg_query("link", "pathway", reaction_ids_with_prefix, batch_size=100, cache_dir=cache_dir)
    
    # Parse results: format is "rn:R00004\tpath:map00010"
    reaction_to_pathways = defaultdict(set)
    for line in results:
        if "\t" in line:
            parts = line.split("\t")
            if len(parts) == 2:
                reaction_part = parts[0].strip()
                pathway_part = parts[1].strip()
                # Extract IDs (pathway format can be "path:map00010" or just "map00010")
                if ":" in reaction_part:
                    reaction_id = reaction_part.split(":")[1]
                    # Pathway can have "path:" prefix or just be "map00010"
                    if ":" in pathway_part:
                        pathway_id = pathway_part.split(":")[1]
                    else:
                        pathway_id = pathway_part
                    reaction_to_pathways[reaction_id].add(pathway_id)
    
    logging.info(f"  Mapped {len(reaction_to_pathways)} reactions to pathways")
    
    # Save mapping
    with safe_open_write(f"{RESULTS_DIR}/intermediate_mappings/reaction_to_pathways.tsv", "w") as wf:
        wf.write("Reaction_ID\tPathway_IDs\n")
        for reaction_id, pathways in sorted(reaction_to_pathways.items()):
            wf.write(f"{reaction_id}\t{','.join(sorted(pathways))}\n")
    
    logging.info("  ✓ Phase 3B complete: Reaction → Pathway mappings saved")
    return dict(reaction_to_pathways)

# ============================================================================
# Phase 4: Build Count Matrices
# ============================================================================

def build_reaction_counts(go_counts_per_genome, go_to_reactions):
    """
    Build reaction-level counts by propagating GO counts to reactions.
    
    Returns:
        dict mapping genome_id -> dict mapping reaction_id -> count
    """
    logging.info("Phase 4A: Building Reaction-level counts")
    reaction_counts_per_genome = defaultdict(dict)
    
    total_genomes = len(go_counts_per_genome)
    genomes_processed = 0
    
    for genome_id, go_counts in go_counts_per_genome.items():
        # For each GO term with counts, propagate to reactions
        for go_term, go_count in go_counts.items():
            # Get all reactions mapped to this GO term
            reactions = go_to_reactions.get(go_term, set())
            for reaction_id in reactions:
                reaction_counts_per_genome[genome_id][reaction_id] = \
                    reaction_counts_per_genome[genome_id].get(reaction_id, 0) + go_count
        
        genomes_processed += 1
        if genomes_processed % 100 == 0:
            logging.info(f"  Processed {genomes_processed}/{total_genomes} genomes...")
    
    # Get all unique reactions that have counts
    all_reactions_with_counts = set()
    for genome_counts in reaction_counts_per_genome.values():
        all_reactions_with_counts.update(genome_counts.keys())
    all_reactions_with_counts = sorted(all_reactions_with_counts)
    
    logging.info(f"  ✓ Processed {len(reaction_counts_per_genome)} genomes")
    logging.info(f"  ✓ Found {len(all_reactions_with_counts)} reactions with counts across genomes")
    logging.info("  ✓ Phase 4A complete: Reaction counts built")
    
    return dict(reaction_counts_per_genome), all_reactions_with_counts

def build_ko_counts(reaction_counts_per_genome, reaction_to_kos):
    """
    Build KO-level counts by aggregating reaction counts.
    
    Returns:
        dict mapping genome_id -> dict mapping ko_id -> count
    """
    logging.info("Phase 4B: Building KO-level counts")
    ko_counts_per_genome = defaultdict(dict)
    
    total_genomes = len(reaction_counts_per_genome)
    genomes_processed = 0
    
    for genome_id, reaction_counts in reaction_counts_per_genome.items():
        for reaction_id, reaction_count in reaction_counts.items():
            # Get all KOs mapped to this reaction
            kos = reaction_to_kos.get(reaction_id, set())
            for ko_id in kos:
                ko_counts_per_genome[genome_id][ko_id] = \
                    ko_counts_per_genome[genome_id].get(ko_id, 0) + reaction_count
        
        genomes_processed += 1
        if genomes_processed % 100 == 0:
            logging.info(f"  Processed {genomes_processed}/{total_genomes} genomes...")
    
    # Get all unique KOs
    all_kos_with_counts = set()
    for genome_counts in ko_counts_per_genome.values():
        all_kos_with_counts.update(genome_counts.keys())
    all_kos_with_counts = sorted(all_kos_with_counts)
    
    logging.info(f"  ✓ Found {len(all_kos_with_counts)} KOs with counts across genomes")
    logging.info("  ✓ Phase 4B complete: KO counts built")
    
    return dict(ko_counts_per_genome), all_kos_with_counts

def build_pathway_counts(reaction_counts_per_genome, reaction_to_pathways):
    """
    Build Pathway-level counts by aggregating reaction counts.
    Multiple reactions can map to the same pathway - counts are cumulative (summed).
    
    Returns:
        dict mapping genome_id -> dict mapping pathway_id -> count
    """
    logging.info("Phase 4C: Building Pathway-level counts (cumulative)")
    pathway_counts_per_genome = defaultdict(dict)
    
    total_genomes = len(reaction_counts_per_genome)
    genomes_processed = 0
    
    for genome_id, reaction_counts in reaction_counts_per_genome.items():
        for reaction_id, reaction_count in reaction_counts.items():
            # Get all pathways mapped to this reaction
            pathways = reaction_to_pathways.get(reaction_id, set())
            for pathway_id in pathways:
                # Cumulative: sum all reaction counts that map to this pathway
                pathway_counts_per_genome[genome_id][pathway_id] = \
                    pathway_counts_per_genome[genome_id].get(pathway_id, 0) + reaction_count
        
        genomes_processed += 1
        if genomes_processed % 100 == 0:
            logging.info(f"  Processed {genomes_processed}/{total_genomes} genomes...")
    
    # Get all unique pathways
    all_pathways_with_counts = set()
    for genome_counts in pathway_counts_per_genome.values():
        all_pathways_with_counts.update(genome_counts.keys())
    all_pathways_with_counts = sorted(all_pathways_with_counts)
    
    logging.info(f"  ✓ Found {len(all_pathways_with_counts)} pathways with counts across genomes")
    logging.info("  ✓ Phase 4C complete: Pathway counts built (cumulative)")
    
    return dict(pathway_counts_per_genome), all_pathways_with_counts

# ============================================================================
# Phase 5: Write Output Files
# ============================================================================

def write_output_files(reaction_counts_per_genome, ko_counts_per_genome, 
                       pathway_counts_per_genome, all_reactions_with_counts,
                       all_kos_with_counts, all_pathways_with_counts):
    """Write all output files to RESULTS_DIR."""
    logging.info("Phase 5: Writing output files")
    
    # 5A: Identify ubiquitous reactions (≥95% prevalence)
    logging.info("  5A: Identifying ubiquitous reactions...")
    reaction_prevalence = defaultdict(int)
    for genome_counts in reaction_counts_per_genome.values():
        for reaction_id in genome_counts.keys():
            reaction_prevalence[reaction_id] += 1
    
    abundance_cutoff = len(reaction_counts_per_genome) * 0.95
    ubiquitous_reactions = [r for r in all_reactions_with_counts 
                            if reaction_prevalence[r] >= abundance_cutoff]
    ubiquitous_reactions = sorted(ubiquitous_reactions)
    
    logging.info(f"  Found {len(ubiquitous_reactions)} ubiquitous reactions (≥95% prevalence)")
    
    # Write ubiquitous reaction counts table
    logging.info("  Writing ubiquitous_counts_table_kegg_reaction.txt...")
    with safe_open_write(f"{RESULTS_DIR}/ubiquitous_counts_table_kegg_reaction.txt", "w") as wf:
        wf.write("Genome\t" + "\t".join(ubiquitous_reactions) + "\n")
        for genome_id in sorted(reaction_counts_per_genome.keys()):
            counts_list = [str(reaction_counts_per_genome[genome_id].get(r, 0)) 
                          for r in ubiquitous_reactions]
            wf.write(f"{genome_id}\t" + "\t".join(counts_list) + "\n")
    logging.info(f"  ✓ Written: ubiquitous_counts_table_kegg_reaction.txt")
    
    # 5B: Write individual reaction files
    logging.info(f"  5B: Writing {len(all_reactions_with_counts)} reaction files...")
    files_written = 0
    for reaction_id in all_reactions_with_counts:
        reaction_file = f"{RESULTS_DIR}/REACTION_lists_with_counts_/{reaction_id}.txt"
        with safe_open_write(reaction_file, "w") as wf:
            for genome_id in sorted(reaction_counts_per_genome.keys()):
                count = reaction_counts_per_genome[genome_id].get(reaction_id, 0)
                if count > 0:
                    wf.write(f"{genome_id}\t{count}\n")
        files_written += 1
        if files_written % 100 == 0:
            logging.info(f"    Written {files_written}/{len(all_reactions_with_counts)} reaction files...")
    logging.info(f"  ✓ Written {len(all_reactions_with_counts)} reaction files")
    
    # 5C: Write individual KO files
    logging.info(f"  5C: Writing {len(all_kos_with_counts)} KO files...")
    files_written = 0
    for ko_id in all_kos_with_counts:
        ko_file = f"{RESULTS_DIR}/KO_lists_with_counts_/{ko_id}.txt"
        with safe_open_write(ko_file, "w") as wf:
            for genome_id in sorted(ko_counts_per_genome.keys()):
                count = ko_counts_per_genome[genome_id].get(ko_id, 0)
                if count > 0:
                    wf.write(f"{genome_id}\t{count}\n")
        files_written += 1
        if files_written % 100 == 0:
            logging.info(f"    Written {files_written}/{len(all_kos_with_counts)} KO files...")
    logging.info(f"  ✓ Written {len(all_kos_with_counts)} KO files")
    
    # 5D: Write individual pathway files
    logging.info(f"  5D: Writing {len(all_pathways_with_counts)} pathway files...")
    files_written = 0
    for pathway_id in all_pathways_with_counts:
        pathway_file = f"{RESULTS_DIR}/PATHWAY_lists_with_counts_/{pathway_id}.txt"
        with safe_open_write(pathway_file, "w") as wf:
            for genome_id in sorted(pathway_counts_per_genome.keys()):
                count = pathway_counts_per_genome[genome_id].get(pathway_id, 0)
                if count > 0:
                    wf.write(f"{genome_id}\t{count}\n")
        files_written += 1
        if files_written % 100 == 0:
            logging.info(f"    Written {files_written}/{len(all_pathways_with_counts)} pathway files...")
    logging.info(f"  ✓ Written {len(all_pathways_with_counts)} pathway files")
    
    # 5E: Write all KO counts table
    logging.info("  5E: Writing all_ko_counts_table_kegg_reaction.txt...")
    with safe_open_write(f"{RESULTS_DIR}/all_ko_counts_table_kegg_reaction.txt", "w") as wf:
        wf.write("Genome\t" + "\t".join(all_kos_with_counts) + "\n")
        for genome_id in sorted(ko_counts_per_genome.keys()):
            counts_list = [str(ko_counts_per_genome[genome_id].get(ko_id, 0)) 
                          for ko_id in all_kos_with_counts]
            wf.write(f"{genome_id}\t" + "\t".join(counts_list) + "\n")
    logging.info(f"  ✓ Written: all_ko_counts_table_kegg_reaction.txt")
    
    # 5F: Write all Pathway counts table
    logging.info("  5F: Writing all_pathway_counts_table_kegg_reaction.txt...")
    with safe_open_write(f"{RESULTS_DIR}/all_pathway_counts_table_kegg_reaction.txt", "w") as wf:
        wf.write("Genome\t" + "\t".join(all_pathways_with_counts) + "\n")
        for genome_id in sorted(pathway_counts_per_genome.keys()):
            counts_list = [str(pathway_counts_per_genome[genome_id].get(pathway_id, 0)) 
                          for pathway_id in all_pathways_with_counts]
            wf.write(f"{genome_id}\t" + "\t".join(counts_list) + "\n")
    logging.info(f"  ✓ Written: all_pathway_counts_table_kegg_reaction.txt")
    
    logging.info("  ✓ Phase 5 complete: All output files written")

# ============================================================================
# Phase 6: Assemble Tables from Existing Files (Skip CPU-heavy parts)
# ============================================================================

def assemble_tables_from_files():
    """
    Assemble all_ko_counts_table and all_pathway_counts_table from existing
    individual files. This allows skipping CPU-heavy parts and just assembling
    the final tables.
    """
    logging.info("Phase 6: Assembling tables from existing files")
    
    # 6A: Assemble KO counts table
    logging.info("  6A: Assembling all_ko_counts_table_kegg_reaction.txt from KO_lists_with_counts_/...")
    ko_dir = f"{RESULTS_DIR}/KO_lists_with_counts_"
    
    if not os.path.exists(ko_dir):
        logging.error(f"  ERROR: Directory not found: {ko_dir}")
        return
    
    # Get all KO files
    ko_files = sorted(glob.glob(f"{ko_dir}/*.txt"))
    if not ko_files:
        logging.error(f"  ERROR: No KO files found in {ko_dir}")
        return
    
    # Extract KO IDs from filenames
    all_kos = sorted([os.path.basename(f).replace(".txt", "") for f in ko_files])
    logging.info(f"  Found {len(all_kos)} KO files")
    
    # Collect all genomes and counts
    all_genomes = set()
    ko_counts_per_genome = defaultdict(dict)
    
    for ko_file in ko_files:
        ko_id = os.path.basename(ko_file).replace(".txt", "")
        with open(ko_file, "r") as rf:
            for line in rf:
                fields = line.strip().split("\t")
                if len(fields) == 2:
                    genome_id = fields[0]
                    try:
                        count = int(fields[1])
                        all_genomes.add(genome_id)
                        ko_counts_per_genome[genome_id][ko_id] = count
                    except ValueError:
                        logging.warning(f"  Invalid count in {ko_file}: {line.strip()}")
    
    all_genomes = sorted(all_genomes)
    logging.info(f"  Found {len(all_genomes)} genomes")
    
    # Write KO counts table
    with safe_open_write(f"{RESULTS_DIR}/all_ko_counts_table_kegg_reaction.txt", "w") as wf:
        wf.write("Genome\t" + "\t".join(all_kos) + "\n")
        for genome_id in all_genomes:
            counts_list = [str(ko_counts_per_genome[genome_id].get(ko_id, 0)) 
                          for ko_id in all_kos]
            wf.write(f"{genome_id}\t" + "\t".join(counts_list) + "\n")
    logging.info(f"  ✓ Written: all_ko_counts_table_kegg_reaction.txt")
    
    # 6B: Assemble Pathway counts table
    logging.info("  6B: Assembling all_pathway_counts_table_kegg_reaction.txt from PATHWAY_lists_with_counts_/...")
    pathway_dir = f"{RESULTS_DIR}/PATHWAY_lists_with_counts_"
    
    if not os.path.exists(pathway_dir):
        logging.error(f"  ERROR: Directory not found: {pathway_dir}")
        return
    
    # Get all pathway files
    pathway_files = sorted(glob.glob(f"{pathway_dir}/*.txt"))
    if not pathway_files:
        logging.error(f"  ERROR: No pathway files found in {pathway_dir}")
        return
    
    # Extract pathway IDs from filenames
    all_pathways = sorted([os.path.basename(f).replace(".txt", "") for f in pathway_files])
    logging.info(f"  Found {len(all_pathways)} pathway files")
    
    # Collect all genomes and counts
    all_genomes_pathway = set()
    pathway_counts_per_genome = defaultdict(dict)
    
    for pathway_file in pathway_files:
        pathway_id = os.path.basename(pathway_file).replace(".txt", "")
        with open(pathway_file, "r") as rf:
            for line in rf:
                fields = line.strip().split("\t")
                if len(fields) == 2:
                    genome_id = fields[0]
                    try:
                        count = int(fields[1])
                        all_genomes_pathway.add(genome_id)
                        pathway_counts_per_genome[genome_id][pathway_id] = count
                    except ValueError:
                        logging.warning(f"  Invalid count in {pathway_file}: {line.strip()}")
    
    all_genomes_pathway = sorted(all_genomes_pathway)
    logging.info(f"  Found {len(all_genomes_pathway)} genomes")
    
    # Write Pathway counts table
    with safe_open_write(f"{RESULTS_DIR}/all_pathway_counts_table_kegg_reaction.txt", "w") as wf:
        wf.write("Genome\t" + "\t".join(all_pathways) + "\n")
        for genome_id in all_genomes_pathway:
            counts_list = [str(pathway_counts_per_genome[genome_id].get(pathway_id, 0)) 
                          for pathway_id in all_pathways]
            wf.write(f"{genome_id}\t" + "\t".join(counts_list) + "\n")
    logging.info(f"  ✓ Written: all_pathway_counts_table_kegg_reaction.txt")
    
    logging.info("  ✓ Phase 6 complete: Tables assembled from existing files")

# ============================================================================
# Main Function
# ============================================================================

def main(test_mode=False, max_genomes=10):
    """Main execution function."""
    start_time = time.time()
    
    # Setup
    log_file = setup_logging(test_mode=test_mode)
    logging.info("=" * 80)
    logging.info("KEGG GO Analyses Script")
    logging.info(f"Mode: {'TEST' if test_mode else 'FULL'}")
    if test_mode:
        logging.info(f"Max genomes: {max_genomes}")
    logging.info("=" * 80)
    
    # Create output directories
    logging.info("Creating output directories...")
    os.makedirs(f"{RESULTS_DIR}/REACTION_lists_with_counts_", exist_ok=True)
    os.makedirs(f"{RESULTS_DIR}/KO_lists_with_counts_", exist_ok=True)
    os.makedirs(f"{RESULTS_DIR}/PATHWAY_lists_with_counts_", exist_ok=True)
    os.makedirs(f"{RESULTS_DIR}/intermediate_mappings", exist_ok=True)
    os.makedirs(f"{RESULTS_DIR}/kegg_api_cache", exist_ok=True)
    
    # Verify input files
    logging.info("Verifying input files...")
    assert os.path.exists(KEGG_MAPPING_FILE), f"Missing: {KEGG_MAPPING_FILE}"
    assert os.path.exists(GO_COUNTS_DIR), f"Missing: {GO_COUNTS_DIR}"
    logging.info("  ✓ All required input files found")
    
    try:
        # Phase 1: Parse GO → Reaction mapping
        go_to_reactions, reaction_to_gos = parse_kegg_reaction2go(KEGG_MAPPING_FILE)
        
        # Phase 2: Load GO counts
        go_counts_per_genome = load_go_counts_from_genome_files(
            GO_COUNTS_DIR, test_mode=test_mode, max_genomes=max_genomes
        )
        
        # Phase 3: Get KEGG API mappings
        all_reaction_ids = sorted(set(reaction_to_gos.keys()))
        logging.info(f"Total unique reactions to query: {len(all_reaction_ids)}")
        
        reaction_to_kos = get_reaction_to_ko_mappings(
            all_reaction_ids, 
            cache_dir=f"{RESULTS_DIR}/kegg_api_cache"
        )
        
        reaction_to_pathways = get_reaction_to_pathway_mappings(
            all_reaction_ids,
            cache_dir=f"{RESULTS_DIR}/kegg_api_cache"
        )
        
        # Phase 4: Build count matrices
        reaction_counts_per_genome, all_reactions_with_counts = build_reaction_counts(
            go_counts_per_genome, go_to_reactions
        )
        
        ko_counts_per_genome, all_kos_with_counts = build_ko_counts(
            reaction_counts_per_genome, reaction_to_kos
        )
        
        pathway_counts_per_genome, all_pathways_with_counts = build_pathway_counts(
            reaction_counts_per_genome, reaction_to_pathways
        )
        
        # Phase 5: Write outputs
        write_output_files(
            reaction_counts_per_genome, ko_counts_per_genome, pathway_counts_per_genome,
            all_reactions_with_counts, all_kos_with_counts, all_pathways_with_counts
        )
        
        # Summary
        elapsed_time = time.time() - start_time
        logging.info("=" * 80)
        logging.info("SUCCESS: All phases completed!")
        logging.info(f"  Genomes processed: {len(go_counts_per_genome)}")
        logging.info(f"  Reactions found: {len(all_reactions_with_counts)}")
        logging.info(f"  KOs found: {len(all_kos_with_counts)}")
        logging.info(f"  Pathways found: {len(all_pathways_with_counts)}")
        logging.info(f"  Total time: {elapsed_time:.2f} seconds ({elapsed_time/60:.2f} minutes)")
        logging.info(f"  Log file: {log_file}")
        logging.info("=" * 80)
        
    except Exception as e:
        logging.error(f"ERROR: {e}", exc_info=True)
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="KEGG GO Analyses Script")
    parser.add_argument("--test", action="store_true", 
                       help="Run in test mode on subset of genomes")
    parser.add_argument("--max-genomes", type=int, default=10,
                       help="Maximum number of genomes for test mode (default: 10)")
    parser.add_argument("--assemble-only", action="store_true",
                       help="Skip CPU-heavy parts and only assemble tables from existing files")
    args = parser.parse_args()
    
    if args.assemble_only:
        log_file = setup_logging()
        logging.info("=" * 80)
        logging.info("KEGG GO Analyses Script - Assembly Only Mode")
        logging.info("=" * 80)
        assemble_tables_from_files()
    else:
        main(test_mode=args.test, max_genomes=args.max_genomes)

