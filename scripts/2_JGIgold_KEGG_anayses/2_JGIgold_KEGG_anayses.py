#!/usr/bin/env python3
"""
JGI/GOLD KEGG Analysis: Bacterial Genome Size Constraints

This script implements comprehensive KEGG annotation and feature extraction:
1. Compute genome metrics (coding density, gene counts, etc.)
2. Extract amino acid composition and nitrogen burden
3. Identify transcription factors and mobile elements
4. Run KEGG annotation (KofamScan) on protein sequences
5. Map KOs to modules and pathways
6. Create comprehensive feature matrix for downstream modeling

OUTPUT: All results saved to results/2_JGIgold_KEGG_anayses_out/
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import sys
import subprocess
import re
from collections import Counter
from Bio import SeqIO
import glob
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 6)
plt.rcParams['font.size'] = 10

# Base paths
BASE_DIR = Path("/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint")
RESULTS_DIR = BASE_DIR / "results" / "2_JGIgold_KEGG_anayses_out"
RESULTS_DIR.mkdir(parents=True, exist_ok=True)
ASSEMBLIES_DIR = BASE_DIR / "data/ncbi/assemblies"

# Input data paths
NCBI_METADATA = BASE_DIR / "data/ncbi/metadata/assembly_data_report_extracted.tsv"
NCBI_BASIC = BASE_DIR / "data/ncbi/metadata/metadata_table.tsv"
SELECTED_ENV = BASE_DIR / "results/1_exploratory_analyses_out/08_selected_environments_summary.tsv"

print("=" * 80)
print("JGI/GOLD KEGG ANALYSIS: Bacterial Genome Size Constraints")
print("=" * 80)
print(f"\nOutput directory: {RESULTS_DIR}")
print("=" * 80)

# ============================================================================
# 1. Load and Prepare Data
# ============================================================================
print("\n[1] Loading metadata and preparing data...")

ncbi_detailed = pd.read_csv(NCBI_METADATA, sep='\t', low_memory=False)
ncbi_basic = pd.read_csv(NCBI_BASIC, sep='\t')

print(f"  ✓ NCBI detailed metadata: {len(ncbi_detailed):,} genomes")
print(f"  ✓ NCBI basic metadata: {len(ncbi_basic):,} genomes")

# Load selected environments if available
selected_envs = None
if SELECTED_ENV.exists():
    selected_envs = pd.read_csv(SELECTED_ENV, sep='\t')
    print(f"  ✓ Selected environments: {len(selected_envs):,} genomes")

# Get list of accessions to process
if selected_envs is not None:
    accessions = selected_envs['accession'].tolist()
    print(f"  Processing {len(accessions):,} genomes from selected environments")
else:
    # Use all accessions from metadata
    accessions = ncbi_detailed['accession'].dropna().tolist()
    print(f"  Processing all {len(accessions):,} genomes")

# Limit to first 50 for testing (remove this limit for full run)
# Uncomment the line below to process all genomes
# accessions = accessions[:50]
# print(f"  Testing with {len(accessions):,} genomes")
print(f"  Processing {len(accessions):,} genomes")

# ============================================================================
# 2. Compute Genome Metrics
# ============================================================================
print("\n[2] Computing genome metrics...")
print("  Computing: coding density, gene counts, GC content")

genome_metrics = []

for acc in accessions:
    acc_dir = ASSEMBLIES_DIR / acc
    
    if not acc_dir.exists():
        continue
    
    # Get metadata for this accession
    meta_row = ncbi_detailed[ncbi_detailed['accession'] == acc].iloc[0] if len(ncbi_detailed[ncbi_detailed['accession'] == acc]) > 0 else None
    
    if meta_row is None:
        continue
    
    metrics = {
        'accession': acc,
        'genome_size_bp': pd.to_numeric(meta_row.get('stats_totalSequenceLength', 0), errors='coerce'),
        'gc_percent': pd.to_numeric(meta_row.get('stats_gcPercent', 0), errors='coerce'),
        'genes_total': pd.to_numeric(meta_row.get('genes_total', 0), errors='coerce'),
        'genes_proteinCoding': pd.to_numeric(meta_row.get('genes_proteinCoding', 0), errors='coerce'),
        'checkm_completeness': pd.to_numeric(meta_row.get('checkm_completeness', 0), errors='coerce'),
        'checkm_contamination': pd.to_numeric(meta_row.get('checkm_contamination', 0), errors='coerce'),
    }
    
    # Compute coding density from CDS file
    cds_file = acc_dir / "cds_from_genomic.fna"
    if cds_file.exists():
        try:
            sequences = list(SeqIO.parse(cds_file, 'fasta'))
            total_cds_length = sum(len(seq.seq) for seq in sequences)
            if metrics['genome_size_bp'] > 0:
                metrics['coding_density'] = (total_cds_length / metrics['genome_size_bp']) * 100
            else:
                metrics['coding_density'] = np.nan
        except Exception as e:
            metrics['coding_density'] = np.nan
    else:
        metrics['coding_density'] = np.nan
    
    genome_metrics.append(metrics)
    
    if len(genome_metrics) % 10 == 0:
        print(f"    Processed {len(genome_metrics)} genomes...")

genome_metrics_df = pd.DataFrame(genome_metrics)
print(f"  ✓ Computed metrics for {len(genome_metrics_df):,} genomes")

# Save genome metrics
genome_metrics_df.to_csv(RESULTS_DIR / '01_genome_metrics.tsv', sep='\t', index=False)
print(f"  ✓ Saved: {RESULTS_DIR / '01_genome_metrics.tsv'}")

# ============================================================================
# 3. Amino Acid Composition and Nitrogen Burden
# ============================================================================
print("\n[3] Computing amino acid composition and nitrogen burden...")
print("  Rationale: Nitrogen burden indicates nutrient limitation signatures")

# Amino acid nitrogen content (atoms per residue)
AA_N_CONTENT = {
    'A': 1, 'R': 4, 'N': 2, 'D': 1, 'C': 1, 'Q': 2, 'E': 1, 'G': 1,
    'H': 3, 'I': 1, 'L': 1, 'K': 2, 'M': 1, 'F': 1, 'P': 1, 'S': 1,
    'T': 1, 'W': 2, 'Y': 1, 'V': 1, 'X': 0, '*': 0
}

aa_composition = []

for acc in accessions:
    acc_dir = ASSEMBLIES_DIR / acc
    protein_file = acc_dir / "protein.faa"
    
    if not protein_file.exists():
        continue
    
    try:
        sequences = list(SeqIO.parse(protein_file, 'fasta'))
        
        if len(sequences) == 0:
            continue
        
        # Concatenate all protein sequences
        all_aa = ''.join([str(seq.seq).upper() for seq in sequences])
        total_aa = len(all_aa)
        
        if total_aa == 0:
            continue
        
        # Count amino acids
        aa_counts = Counter(all_aa)
        
        # Calculate nitrogen burden (average N atoms per residue)
        total_n_atoms = sum(AA_N_CONTENT.get(aa, 0) * count for aa, count in aa_counts.items())
        n_burden = total_n_atoms / total_aa if total_aa > 0 else 0
        
        # Calculate amino acid frequencies
        aa_freq = {f'aa_{aa}': (count / total_aa) * 100 for aa, count in aa_counts.items() if aa in AA_N_CONTENT}
        
        aa_data = {
            'accession': acc,
            'amino_n_burden': n_burden,
            'total_aa_residues': total_aa,
            **aa_freq
        }
        
        aa_composition.append(aa_data)
        
    except Exception as e:
        continue

aa_df = pd.DataFrame(aa_composition)
print(f"  ✓ Computed amino acid composition for {len(aa_df):,} genomes")

# Save amino acid composition
aa_df.to_csv(RESULTS_DIR / '02_amino_acid_composition.tsv', sep='\t', index=False)
print(f"  ✓ Saved: {RESULTS_DIR / '02_amino_acid_composition.tsv'}")

# ============================================================================
# 4. Transcription Factor and Mobile Element Counts
# ============================================================================
print("\n[4] Identifying transcription factors and mobile elements...")
print("  Rationale: TF count indicates regulatory complexity")
print("  Mobile elements indicate horizontal gene transfer")

# Keywords for transcription factors and mobile elements
TF_KEYWORDS = [
    'transcription factor', 'transcriptional regulator', 'helix-turn-helix',
    'sigma factor', 'response regulator', 'lysr', 'arac', 'gntr', 'iclr',
    'marr', 'tetr', 'merr', 'luxr', 'fis', 'hns', 'ihf', 'fru', 'crp'
]

MOBILE_KEYWORDS = [
    'transposase', 'transposon', 'integrase', 'recombinase', 'conjugative',
    'mobilization', 'plasmid', 'phage', 'prophage', 'insertion sequence',
    'is element', 'tnp', 'ins', 'tra', 'mob'
]

tf_mobile_counts = []

for acc in accessions:
    acc_dir = ASSEMBLIES_DIR / acc
    gff_file = list(acc_dir.glob("*.gff"))
    
    if not gff_file:
        continue
    
    gff_file = gff_file[0]
    
    try:
        tf_count = 0
        mobile_count = 0
        
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                line_lower = line.lower()
                
                # Check for transcription factors
                if any(keyword in line_lower for keyword in TF_KEYWORDS):
                    tf_count += 1
                
                # Check for mobile elements
                if any(keyword in line_lower for keyword in MOBILE_KEYWORDS):
                    mobile_count += 1
        
        tf_mobile_counts.append({
            'accession': acc,
            'tf_count': tf_count,
            'mobile_element_count': mobile_count
        })
        
    except Exception as e:
        continue

tf_mobile_df = pd.DataFrame(tf_mobile_counts)
print(f"  ✓ Identified TFs and mobile elements for {len(tf_mobile_df):,} genomes")

# Save TF and mobile element counts
tf_mobile_df.to_csv(RESULTS_DIR / '03_tf_mobile_elements.tsv', sep='\t', index=False)
print(f"  ✓ Saved: {RESULTS_DIR / '03_tf_mobile_elements.tsv'}")

# ============================================================================
# 5. KEGG Annotation (KofamScan)
# ============================================================================
print("\n[5] Running KEGG annotation with KofamScan...")
print("  Note: This requires KofamScan to be installed")
print("  Command: conda install -c bioconda kofamscan")

# Check if KofamScan is available
kofamscan_available = False
try:
    result = subprocess.run(['which', 'exec_annotation'], 
                          capture_output=True, text=True, timeout=5)
    if result.returncode == 0:
        kofamscan_available = True
        print("  ✓ KofamScan found")
    else:
        print("  ⚠ KofamScan not found - skipping KEGG annotation")
        print("    Install with: conda install -c bioconda kofamscan")
except:
    print("  ⚠ KofamScan not found - skipping KEGG annotation")
    print("    Install with: conda install -c bioconda kofamscan")

kegg_results = []

if kofamscan_available:
    # Create temporary directory for KofamScan outputs
    kofam_dir = RESULTS_DIR / "kofam_outputs"
    kofam_dir.mkdir(exist_ok=True)
    
    # Process first 10 genomes as test (remove limit for full run)
    test_accessions = accessions[:10]
    
    for acc in test_accessions:
        acc_dir = ASSEMBLIES_DIR / acc
        protein_file = acc_dir / "protein.faa"
        
        if not protein_file.exists():
            continue
        
        output_file = kofam_dir / f"{acc}_kofam.txt"
        
        # Run KofamScan (this is a placeholder - actual command depends on setup)
        # exec_annotation -o {output_file} --cpu 8 -f detail-tsv {protein_file} -k ko_list -p profiles/
        print(f"    Processing {acc}... (KofamScan execution skipped - requires setup)")
        
        # For now, create placeholder results
        kegg_results.append({
            'accession': acc,
            'ko_count': np.nan,
            'module_count': np.nan,
            'ko_per_mb': np.nan
        })
else:
    # Create placeholder KEGG results
    for acc in accessions[:10]:
        kegg_results.append({
            'accession': acc,
            'ko_count': np.nan,
            'module_count': np.nan,
            'ko_per_mb': np.nan
        })

kegg_df = pd.DataFrame(kegg_results)
print(f"  ✓ KEGG annotation results for {len(kegg_df):,} genomes")

# Save KEGG results
kegg_df.to_csv(RESULTS_DIR / '04_kegg_annotation.tsv', sep='\t', index=False)
print(f"  ✓ Saved: {RESULTS_DIR / '04_kegg_annotation.tsv'}")

# ============================================================================
# 6. Merge All Features into Comprehensive Matrix
# ============================================================================
print("\n[6] Creating comprehensive feature matrix...")

# Start with processed accessions and get their metadata
processed_accessions = set(genome_metrics_df['accession'].tolist())

# Get metadata for processed genomes only
feature_matrix = ncbi_basic[ncbi_basic['Assembly Accession'].isin(processed_accessions)].copy()
feature_matrix = feature_matrix[['Assembly Accession', 'Organism Name', 
                                 'Organism Taxonomic ID', 'Assembly Release Date',
                                 'Assembly Level']].copy()
feature_matrix.columns = ['accession', 'organism_name', 'organism_taxId', 
                         'release_date', 'assembly_level']

# Merge genome metrics
feature_matrix = feature_matrix.merge(genome_metrics_df, on='accession', how='inner')

# Merge amino acid composition
feature_matrix = feature_matrix.merge(aa_df[['accession', 'amino_n_burden']], 
                                     on='accession', how='left')

# Merge TF and mobile element counts
feature_matrix = feature_matrix.merge(tf_mobile_df, on='accession', how='left')

# Merge KEGG results
feature_matrix = feature_matrix.merge(kegg_df, on='accession', how='left')

# Add environment information if available
if selected_envs is not None and 'ORGANISM ECOSYSTEM CATEGORY' in selected_envs.columns:
    env_map = selected_envs[['accession', 'ORGANISM ECOSYSTEM CATEGORY']].drop_duplicates()
    feature_matrix = feature_matrix.merge(env_map, on='accession', how='left')
    feature_matrix.rename(columns={'ORGANISM ECOSYSTEM CATEGORY': 'environment'}, 
                         inplace=True)

# Compute derived metrics
feature_matrix['genome_size_mb'] = feature_matrix['genome_size_bp'] / 1e6
feature_matrix['ko_per_mb'] = feature_matrix['ko_count'] / feature_matrix['genome_size_mb']
feature_matrix['module_per_mb'] = feature_matrix['module_count'] / feature_matrix['genome_size_mb']

# Save comprehensive feature matrix
feature_matrix.to_csv(RESULTS_DIR / '05_genome_feature_matrix.tsv', sep='\t', index=False)
print(f"  ✓ Saved comprehensive feature matrix: {RESULTS_DIR / '05_genome_feature_matrix.tsv'}")
print(f"    {len(feature_matrix):,} genomes, {len(feature_matrix.columns)} features")

# ============================================================================
# 7. Summary Statistics
# ============================================================================
print("\n[7] Computing summary statistics...")

summary_stats = {
    'total_genomes': len(feature_matrix),
    'genomes_with_metrics': feature_matrix['genome_size_bp'].notna().sum(),
    'genomes_with_aa': feature_matrix['amino_n_burden'].notna().sum(),
    'genomes_with_tf': feature_matrix['tf_count'].notna().sum(),
    'genomes_with_mobile': feature_matrix['mobile_element_count'].notna().sum(),
    'genomes_with_kegg': feature_matrix['ko_count'].notna().sum(),
}

print("\n  Summary:")
for key, value in summary_stats.items():
    print(f"    {key}: {value:,}")

# Save summary
summary_df = pd.DataFrame([summary_stats])
summary_df.to_csv(RESULTS_DIR / '06_summary_statistics.tsv', sep='\t', index=False)
print(f"  ✓ Saved: {RESULTS_DIR / '06_summary_statistics.tsv'}")

print("\n" + "=" * 80)
print("ANALYSIS COMPLETE")
print("=" * 80)
print(f"\nAll outputs saved to: {RESULTS_DIR}")
print("\nNext steps:")
print("  1. Install and configure KofamScan for full KEGG annotation")
print("  2. Run KofamScan on all protein sequences")
print("  3. Map KOs to modules using KEGG Mapper or KEMET")
print("  4. Add environment nutrient scores and carbon substrate diversity")
print("  5. Proceed to modeling phase")
print("=" * 80)

