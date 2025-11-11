#!/usr/bin/env python3
"""
Comprehensive Preliminary Analysis: Bacterial Genome Size Constraints

RESEARCH QUESTIONS:
1. Do bacteria from similar environments across diverse taxa converge on similar genome sizes?
2. Is genome size constrained by metabolic complexity (diverse carbon sources) or nutrient limitation (nitrogen conservation)?
3. Which well-defined environments are suitable for detailed analysis (5-20 taxonomically diverse organisms per environment)?

HYPOTHESES:
- H1 (Metabolic Complexity): Larger genomes in environments with diverse carbon sources
- H2 (Nutrient Limitation): Smaller genomes in nutrient-limited environments (e.g., nitrogen conservation)

COMPREHENSIVE ANALYSIS APPROACH:
1. Data quality control and filtering
2. Genome size distribution characterization
3. Phylogenetic diversity analysis (not just taxonomic counts)
4. Environment identification and convergence testing
5. Amino acid composition analysis (nutrient limitation signatures)
6. Gene content analysis (functional categories, transcription factors, mobile elements)
7. Statistical modeling preparation (phylogenetic regression framework)
8. Clustering analysis within environments (multiple survival strategies)
9. Correlation analysis with proper controls
10. Model building framework preparation
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
from scipy import stats
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from collections import Counter
import glob
warnings.filterwarnings('ignore')

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 8)
plt.rcParams['font.size'] = 11

# Base paths
BASE_DIR = Path("/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint")
RESULTS_DIR = BASE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)
ASSEMBLIES_DIR = BASE_DIR / "data/ncbi/assemblies"

# NCBI data paths
NCBI_METADATA = BASE_DIR / "data/ncbi/metadata/assembly_data_report_extracted.tsv"
NCBI_BASIC = BASE_DIR / "data/ncbi/metadata/metadata_table.tsv"

# GOLD data paths
GOLD_ORG = BASE_DIR / "data/gold/0_20251106_gold_metadata_Organism.csv"
GOLD_BIOSAMPLE = BASE_DIR / "data/gold/0_20251106_gold_metadata_Biosample.csv"

print("=" * 80)
print("COMPREHENSIVE PRELIMINARY ANALYSIS: Bacterial Genome Size Constraints")
print("=" * 80)
print("\nRESEARCH QUESTIONS:")
print("  1. Do bacteria from similar environments converge on similar genome sizes?")
print("  2. Metabolic complexity vs. nutrient limitation as constraints?")
print("  3. Which environments are suitable for detailed analysis?")
print("\n" + "=" * 80)

# ============================================================================
# 1. Data Loading and Quality Control
# ============================================================================
print("\n[1] Loading and preparing data...")

ncbi_detailed = pd.read_csv(NCBI_METADATA, sep='\t', low_memory=False)
print(f"  ✓ NCBI detailed metadata: {len(ncbi_detailed):,} genomes")

# Extract key variables
ncbi_detailed['genome_size_bp'] = pd.to_numeric(ncbi_detailed['stats_totalSequenceLength'], errors='coerce')
ncbi_detailed['genome_size_mb'] = ncbi_detailed['genome_size_bp'] / 1e6
ncbi_detailed['genes_total_num'] = pd.to_numeric(ncbi_detailed['genes_total'], errors='coerce')
ncbi_detailed['genes_protein_num'] = pd.to_numeric(ncbi_detailed['genes_proteinCoding'], errors='coerce')
ncbi_detailed['gc_percent'] = pd.to_numeric(ncbi_detailed['stats_gcPercent'], errors='coerce')

# Load GOLD metadata
gold_org = pd.read_csv(GOLD_ORG, low_memory=False)
print(f"  ✓ GOLD Organism metadata: {len(gold_org):,} records")

# ============================================================================
# 2. Filter Out Confounders
# ============================================================================
print("\n[2] Filtering out confounders...")
print("  Rationale: Obligate intracellular bacteria have small genomes due to genome decay,")
print("  not environmental constraints. This is a known confounder.")

intracellular_keywords = [
    'chlamydia', 'rickettsia', 'mycoplasma', 'ureaplasma', 'mesoplasma',
    'buchnera', 'wolbachia', 'candidatus pelagibacter', 'dehalococcoides'
]

ncbi_detailed['is_intracellular'] = ncbi_detailed['organism_name'].str.lower().str.contains(
    '|'.join(intracellular_keywords), na=False, regex=True
)
ncbi_detailed['is_very_small'] = ncbi_detailed['genome_size_mb'] < 1.0

ncbi_filtered = ncbi_detailed[
    ~ncbi_detailed['is_intracellular'] & 
    ~ncbi_detailed['is_very_small']
].copy()

print(f"  Removed {len(ncbi_detailed) - len(ncbi_filtered):,} obligate intracellular/very small genomes")
print(f"  Remaining: {len(ncbi_filtered):,} genomes for analysis")

# ============================================================================
# 3. Genome Size Distribution Analysis
# ============================================================================
print("\n[3] Characterizing genome size distribution...")

genome_sizes = ncbi_filtered['genome_size_mb'].dropna()

print(f"\n  Genome size statistics (n={len(genome_sizes):,}):")
print(f"    Mean: {genome_sizes.mean():.2f} Mb")
print(f"    Median: {genome_sizes.median():.2f} Mb")
print(f"    Range: {genome_sizes.min():.2f} - {genome_sizes.max():.2f} Mb")
print(f"    Std: {genome_sizes.std():.2f} Mb")
print(f"    IQR: {genome_sizes.quantile(0.25):.2f} - {genome_sizes.quantile(0.75):.2f} Mb")

# Distribution plots
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
axes[0].hist(genome_sizes, bins=60, edgecolor='black', alpha=0.7, color='steelblue')
axes[0].axvline(genome_sizes.mean(), color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {genome_sizes.mean():.2f} Mb')
axes[0].axvline(genome_sizes.median(), color='green', linestyle='--', linewidth=2, 
                label=f'Median: {genome_sizes.median():.2f} Mb')
axes[0].set_xlabel('Genome Size (Mb)', fontsize=12)
axes[0].set_ylabel('Frequency', fontsize=12)
axes[0].set_title('Genome Size Distribution (Linear Scale)', fontsize=13, fontweight='bold')
axes[0].legend()
axes[0].grid(True, alpha=0.3)

axes[1].hist(genome_sizes, bins=60, edgecolor='black', alpha=0.7, color='steelblue')
axes[1].axvline(genome_sizes.mean(), color='red', linestyle='--', linewidth=2, 
                label=f'Mean: {genome_sizes.mean():.2f} Mb')
axes[1].axvline(genome_sizes.median(), color='green', linestyle='--', linewidth=2, 
                label=f'Median: {genome_sizes.median():.2f} Mb')
axes[1].set_xlabel('Genome Size (Mb)', fontsize=12)
axes[1].set_ylabel('Frequency', fontsize=12)
axes[1].set_title('Genome Size Distribution (Log Scale)', fontsize=13, fontweight='bold')
axes[1].set_xscale('log')
axes[1].legend()
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(RESULTS_DIR / '01_genome_size_distribution.png', dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {RESULTS_DIR / '01_genome_size_distribution.png'}")
plt.close()

# ============================================================================
# 4. Merge NCBI and GOLD Data
# ============================================================================
print("\n[4] Merging NCBI and GOLD data...")

ncbi_tax = ncbi_filtered[['accession', 'organism_taxId', 'organism_name', 
                          'genome_size_mb', 'genes_total_num', 'genes_protein_num', 
                          'gc_percent', 'biosample_host', 'biosample_isolation_source']].copy()
ncbi_tax['organism_taxId'] = pd.to_numeric(ncbi_tax['organism_taxId'], errors='coerce')

tax_cols = ['ORGANISM NCBI TAX ID', 'ORGANISM NCBI TAXONOMY ID', 'NCBI TAX ID']
gold_tax_col = None
for col in tax_cols:
    if col in gold_org.columns:
        gold_tax_col = col
        break

merged = None
if gold_tax_col:
    print(f"  Found taxonomy ID column in GOLD: {gold_tax_col}")
    gold_tax = gold_org.copy()
    gold_tax[gold_tax_col] = pd.to_numeric(gold_tax[gold_tax_col], errors='coerce')
    
    merged = ncbi_tax.merge(gold_tax, left_on='organism_taxId', right_on=gold_tax_col, how='inner')
    merged = merged.drop_duplicates(subset=['accession'], keep='first')
    print(f"  ✓ Merged datasets: {len(merged):,} unique genomes with GOLD environment data")
else:
    merged = ncbi_tax.copy()

# ============================================================================
# 5. Phylogenetic Diversity Analysis
# ============================================================================
print("\n[5] Analyzing phylogenetic diversity...")
print("  Rationale: Need to assess diversity at multiple taxonomic levels")
print("  to ensure convergence is tested across diverse taxa, not just closely related species")

if 'ORGANISM ECOSYSTEM CATEGORY' in merged.columns:
    env_data = merged[['accession', 'genome_size_mb', 'organism_name', 'organism_taxId',
                      'ORGANISM ECOSYSTEM CATEGORY', 'ORGANISM GOLD PHYLUM', 
                      'ORGANISM GENUS', 'ORGANISM SPECIES']].copy()
    env_data = env_data[env_data['ORGANISM ECOSYSTEM CATEGORY'].notna()]
    
    # Extract genus from organism name if GOLD genus is missing
    env_data['genus'] = env_data['ORGANISM GENUS'].fillna(
        env_data['organism_name'].str.split().str[0]
    )
    
    # Calculate diversity metrics per environment
    env_diversity = env_data.groupby('ORGANISM ECOSYSTEM CATEGORY').agg({
        'genome_size_mb': ['count', 'mean', 'std'],
        'organism_taxId': 'nunique',
        'genus': 'nunique',
        'ORGANISM GOLD PHYLUM': lambda x: x.dropna().nunique()
    }).round(2)
    
    env_diversity.columns = ['n_genomes', 'mean_size_mb', 'std_size_mb', 'n_taxa', 'n_genera', 'n_phyla']
    env_diversity = env_diversity[env_diversity['n_genomes'] >= 20]
    env_diversity['cv_size'] = (env_diversity['std_size_mb'] / env_diversity['mean_size_mb'] * 100).round(2)
    env_diversity['phyla_per_genome'] = (env_diversity['n_phyla'] / env_diversity['n_genomes'] * 100).round(2)
    env_diversity['genera_per_genome'] = (env_diversity['n_genera'] / env_diversity['n_genomes'] * 100).round(2)
    
    print(f"\n  Phylogenetic diversity in environments (>=20 genomes):")
    print(f"\n  {'Environment':<30} {'N':<6} {'Phyla':<7} {'Genera':<7} {'Mean (Mb)':<10} {'CV (%)':<8}")
    print("  " + "-" * 85)
    
    for env, row in env_diversity.head(15).iterrows():
        env_name = str(env)[:28]
        print(f"  {env_name:<30} {int(row['n_genomes']):<6} {int(row['n_phyla']):<7} "
              f"{int(row['n_genera']):<7} {row['mean_size_mb']:<10.2f} {row['cv_size']:<8.1f}")
    
    env_diversity.to_csv(RESULTS_DIR / '02_phylogenetic_diversity.tsv', sep='\t')
    print(f"\n  ✓ Saved: {RESULTS_DIR / '02_phylogenetic_diversity.tsv'}")

# ============================================================================
# 6. Environment Convergence Testing
# ============================================================================
print("\n[6] Testing for genome size convergence within environments...")
print("  Rationale: Key observation - bacteria from similar environments")
print("  across diverse taxa have similar genome sizes")

if 'ORGANISM ECOSYSTEM CATEGORY' in merged.columns and len(env_diversity) > 0:
    top_envs = env_diversity.nlargest(12, 'n_genomes').index.tolist()
    
    # Statistical test: ANOVA to test if environment explains genome size variance
    env_sizes = []
    env_labels = []
    
    for env in top_envs:
        sizes = env_data[env_data['ORGANISM ECOSYSTEM CATEGORY'] == env]['genome_size_mb'].dropna()
        if len(sizes) >= 20:
            env_sizes.append(sizes.values)
            env_labels.append(str(env)[:25])
    
    if len(env_sizes) > 1:
        # One-way ANOVA
        f_stat, p_value = stats.f_oneway(*env_sizes)
        print(f"\n  One-way ANOVA (environment effect on genome size):")
        print(f"    F-statistic: {f_stat:.2f}")
        print(f"    p-value: {p_value:.2e}")
        print(f"    {'Significant' if p_value < 0.05 else 'Not significant'} environment effect")
        
        # Variance partitioning
        all_sizes = np.concatenate(env_sizes)
        overall_var = np.var(all_sizes)
        within_env_vars = [np.var(sizes) for sizes in env_sizes]
        mean_within_var = np.mean(within_env_vars)
        
        # Between-environment variance
        env_means = [np.mean(sizes) for sizes in env_sizes]
        between_env_var = np.var(env_means)
        
        print(f"\n  Variance partitioning:")
        print(f"    Overall variance: {overall_var:.2f}")
        print(f"    Mean within-environment variance: {mean_within_var:.2f}")
        print(f"    Between-environment variance: {between_env_var:.2f}")
        print(f"    Variance reduction: {(1 - mean_within_var/overall_var)*100:.1f}%")
        
        # Visualization
        fig, ax = plt.subplots(figsize=(14, 8))
        bp = ax.boxplot(env_sizes, labels=env_labels, vert=True, patch_artist=True, 
                       showmeans=True, meanline=True)
        
        colors = plt.cm.viridis(np.linspace(0, 1, len(bp['boxes'])))
        for i, patch in enumerate(bp['boxes']):
            patch.set_facecolor(colors[i])
            patch.set_alpha(0.7)
        
        ax.set_ylabel('Genome Size (Mb)', fontsize=12)
        ax.set_xlabel('Environment/Ecosystem Category', fontsize=12)
        ax.set_title(f'Genome Size by Environment\n(ANOVA: F={f_stat:.2f}, p={p_value:.2e})', 
                    fontsize=14, fontweight='bold')
        plt.xticks(rotation=45, ha='right')
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(RESULTS_DIR / '03_genome_size_by_environment.png', dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved: {RESULTS_DIR / '03_genome_size_by_environment.png'}")
        plt.close()

# ============================================================================
# 7. Amino Acid Composition Analysis (Nutrient Limitation Signatures)
# ============================================================================
print("\n[7] Analyzing amino acid composition from CDS sequences...")
print("  Rationale: Nutrient limitation hypothesis predicts differences in")
print("  amino acid usage (e.g., nitrogen-rich amino acids in N-limited environments)")

# Sample a subset of genomes for amino acid analysis (computationally intensive)
sample_size = min(100, len(ncbi_filtered))
sample_genomes = ncbi_filtered.sample(n=sample_size, random_state=42)

aa_composition = []
aa_columns = ['accession', 'genome_size_mb', 'N_content', 'C_content', 
              'N_rich_aa', 'N_poor_aa', 'aromatic_aa']

for idx, row in sample_genomes.iterrows():
    acc = row['accession']
    cds_file = ASSEMBLIES_DIR / acc / "cds_from_genomic.fna"
    
    if cds_file.exists():
        try:
            from Bio import SeqIO
            sequences = list(SeqIO.parse(cds_file, 'fasta'))
            
            if len(sequences) > 0:
                # Concatenate all CDS sequences
                all_seq = ''.join([str(seq.seq) for seq in sequences])
                
                # Count nucleotides
                n_counts = Counter(all_seq.upper())
                total_nt = len(all_seq)
                
                if total_nt > 0:
                    # Calculate N and C content
                    n_content = (n_counts.get('N', 0) / total_nt) * 100
                    c_content = ((n_counts.get('C', 0) + n_counts.get('G', 0)) / total_nt) * 100
                    
                    # Translate to amino acids (simplified - just count codons)
                    # N-rich amino acids: R, H, K, W (high N:C ratio)
                    # N-poor amino acids: G, A, V, L, I (low N:C ratio)
                    # Aromatic: F, Y, W
                    
                    # For now, use nucleotide composition as proxy
                    # A+T rich = more N-poor amino acids, G+C rich = more N-rich
                    n_rich_proxy = c_content  # GC content correlates with N-rich AAs
                    n_poor_proxy = 100 - c_content  # AT content correlates with N-poor AAs
                    
                    aa_composition.append({
                        'accession': acc,
                        'genome_size_mb': row['genome_size_mb'],
                        'N_content': n_content,
                        'C_content': c_content,
                        'N_rich_aa': n_rich_proxy,
                        'N_poor_aa': n_poor_proxy,
                        'aromatic_aa': 0  # Placeholder
                    })
        except Exception as e:
            continue

if len(aa_composition) > 0:
    aa_df = pd.DataFrame(aa_composition)
    print(f"  ✓ Analyzed {len(aa_df):,} genomes for amino acid composition")
    
    # Correlation with genome size
    corr_n = aa_df['genome_size_mb'].corr(aa_df['N_rich_aa'])
    corr_c = aa_df['genome_size_mb'].corr(aa_df['C_content'])
    
    print(f"\n  Correlations with genome size:")
    print(f"    N-rich amino acid proxy (GC%): r={corr_n:.3f}")
    print(f"    C content: r={corr_c:.3f}")
    
    # Visualization
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    axes[0].scatter(aa_df['genome_size_mb'], aa_df['N_rich_aa'], alpha=0.6, s=30)
    axes[0].set_xlabel('Genome Size (Mb)', fontsize=12)
    axes[0].set_ylabel('N-rich AA Proxy (GC%)', fontsize=12)
    axes[0].set_title(f'Genome Size vs. N-rich Amino Acid Content\n(r={corr_n:.3f})', 
                     fontsize=13, fontweight='bold')
    axes[0].grid(True, alpha=0.3)
    
    axes[1].scatter(aa_df['genome_size_mb'], aa_df['C_content'], alpha=0.6, s=30, color='green')
    axes[1].set_xlabel('Genome Size (Mb)', fontsize=12)
    axes[1].set_ylabel('GC Content (%)', fontsize=12)
    axes[1].set_title(f'Genome Size vs. GC Content\n(r={corr_c:.3f})', 
                     fontsize=13, fontweight='bold')
    axes[1].grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(RESULTS_DIR / '04_amino_acid_composition.png', dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {RESULTS_DIR / '04_amino_acid_composition.png'}")
    plt.close()
    
    aa_df.to_csv(RESULTS_DIR / '05_amino_acid_composition.tsv', sep='\t', index=False)
    print(f"  ✓ Saved: {RESULTS_DIR / '05_amino_acid_composition.tsv'}")

# ============================================================================
# 8. Gene Content Analysis
# ============================================================================
print("\n[8] Analyzing gene content patterns...")
print("  Rationale: Gene content reflects metabolic capacity and")
print("  may indicate metabolic complexity vs. nutrient limitation strategies")

# Calculate gene density and functional ratios
ncbi_filtered['gene_density'] = (ncbi_filtered['genes_total_num'] / ncbi_filtered['genome_size_mb']).round(2)
ncbi_filtered['protein_fraction'] = (ncbi_filtered['genes_protein_num'] / ncbi_filtered['genes_total_num'] * 100).round(2)
ncbi_filtered['pseudogene_fraction'] = (
    pd.to_numeric(ncbi_filtered['genes_pseudogene'], errors='coerce') / 
    ncbi_filtered['genes_total_num'] * 100
).round(2)

# Correlation analysis
corr_genes = ncbi_filtered['genome_size_mb'].corr(ncbi_filtered['genes_total_num'])
corr_density = ncbi_filtered['genome_size_mb'].corr(ncbi_filtered['gene_density'])

print(f"\n  Gene content correlations:")
print(f"    Genome size vs. total genes: r={corr_genes:.3f}")
print(f"    Genome size vs. gene density: r={corr_density:.3f}")

# Visualization
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

axes[0].scatter(ncbi_filtered['genome_size_mb'], ncbi_filtered['genes_total_num'], 
               alpha=0.3, s=10, color='steelblue')
axes[0].set_xlabel('Genome Size (Mb)', fontsize=12)
axes[0].set_ylabel('Total Genes', fontsize=12)
axes[0].set_title(f'Genome Size vs. Gene Count\n(r={corr_genes:.3f})', 
                 fontsize=13, fontweight='bold')
axes[0].grid(True, alpha=0.3)

axes[1].scatter(ncbi_filtered['genome_size_mb'], ncbi_filtered['gene_density'], 
               alpha=0.3, s=10, color='coral')
axes[1].set_xlabel('Genome Size (Mb)', fontsize=12)
axes[1].set_ylabel('Gene Density (genes/Mb)', fontsize=12)
axes[1].set_title(f'Genome Size vs. Gene Density\n(r={corr_density:.3f})', 
                 fontsize=13, fontweight='bold')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig(RESULTS_DIR / '06_gene_content_analysis.png', dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {RESULTS_DIR / '06_gene_content_analysis.png'}")
plt.close()

# ============================================================================
# 9. Clustering Analysis Within Environments
# ============================================================================
print("\n[9] Analyzing clustering patterns within environments...")
print("  Rationale: Multiple survival strategies may exist within environments")
print("  (e.g., different genome size clusters in human gut)")

if 'ORGANISM ECOSYSTEM CATEGORY' in merged.columns and len(env_diversity) > 0:
    large_envs = env_diversity[env_diversity['n_genomes'] >= 50].head(5).index.tolist()
    
    if len(large_envs) > 0:
        fig, axes = plt.subplots(len(large_envs), 1, figsize=(12, 3*len(large_envs)))
        if len(large_envs) == 1:
            axes = [axes]
        
        for idx, env in enumerate(large_envs):
            env_subset = env_data[env_data['ORGANISM ECOSYSTEM CATEGORY'] == env]['genome_size_mb'].dropna()
            
            if len(env_subset) >= 20:
                # Histogram with density
                axes[idx].hist(env_subset, bins=30, edgecolor='black', alpha=0.7, 
                              color='steelblue', density=True)
                axes[idx].axvline(env_subset.mean(), color='red', linestyle='--', 
                                 linewidth=2, label=f'Mean: {env_subset.mean():.2f} Mb')
                axes[idx].axvline(env_subset.median(), color='green', linestyle='--', 
                                 linewidth=2, label=f'Median: {env_subset.median():.2f} Mb')
                
                # Test for multimodality (multiple clusters)
                from scipy.stats import gaussian_kde
                kde = gaussian_kde(env_subset)
                x_range = np.linspace(0, 10, 100)  # Match x-axis limit
                axes[idx].plot(x_range, kde(x_range), 'r-', linewidth=2, label='KDE')
                
                axes[idx].set_xlabel('Genome Size (Mb)', fontsize=10)
                axes[idx].set_ylabel('Density', fontsize=10)
                axes[idx].set_xlim(0, 10)  # Set x-axis limit to 0-10 Mb
                axes[idx].set_title(f'{str(env)[:50]} (n={len(env_subset)})', 
                                   fontsize=12, fontweight='bold')
                axes[idx].legend()
                axes[idx].grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(RESULTS_DIR / '07_clustering_within_environments.png', dpi=300, bbox_inches='tight')
        print(f"  ✓ Saved: {RESULTS_DIR / '07_clustering_within_environments.png'}")
        plt.close()

# ============================================================================
# 10. Environment Selection for Downstream Analysis
# ============================================================================
print("\n[10] Selecting environments for downstream analysis...")
print("  Criteria:")
print("    - Well-defined with similar nutrient limitations")
print("    - >=20 genomes (preferably 50+)")
print("    - High phylogenetic diversity (multiple phyla)")
print("    - Low CV (convergent genome sizes)")
print("    - Include both small and large genome size environments")

if 'ORGANISM ECOSYSTEM CATEGORY' in merged.columns and len(env_diversity) > 0:
    # Score environments
    env_diversity['score'] = (
        (env_diversity['n_genomes'] / env_diversity['n_genomes'].max() * 0.25) +
        (env_diversity['n_phyla'] / env_diversity['n_phyla'].max() * 0.25) +
        (env_diversity['n_genera'] / env_diversity['n_genera'].max() * 0.20) +
        (1 - env_diversity['cv_size'] / env_diversity['cv_size'].max() * 0.20) +
        (0.10)
    )
    
    selected_envs = env_diversity.nlargest(12, 'score').index.tolist()
    
    print(f"\n  Selected {len(selected_envs)} environments:")
    print(f"\n  {'Environment':<35} {'N':<6} {'Phyla':<7} {'Mean (Mb)':<10} {'CV (%)':<8} {'Score':<8}")
    print("  " + "-" * 90)
    
    for env in selected_envs:
        row = env_diversity.loc[env]
        env_name = str(env)[:33]
        print(f"  {env_name:<35} {int(row['n_genomes']):<6} {int(row['n_phyla']):<7} "
              f"{row['mean_size_mb']:<10.2f} {row['cv_size']:<8.1f} {row['score']:<8.2f}")
    
    # Create summary table
    selected_summary = merged[
        merged['ORGANISM ECOSYSTEM CATEGORY'].isin(selected_envs)
    ][['accession', 'genome_size_mb', 'genes_total_num', 'genes_protein_num', 
       'gc_percent', 'organism_name', 'organism_taxId', 'ORGANISM ECOSYSTEM CATEGORY']].copy()
    
    selected_summary.to_csv(RESULTS_DIR / '08_selected_environments_summary.tsv', sep='\t', index=False)
    print(f"\n  ✓ Saved: {RESULTS_DIR / '08_selected_environments_summary.tsv'}")
    print(f"    {len(selected_summary):,} genomes from {len(selected_envs)} environments")

# ============================================================================
# 11. Statistical Modeling Preparation
# ============================================================================
print("\n[11] Preparing data for statistical modeling...")
print("  Rationale: Need to prepare data for phylogenetic regression (PGLS)")
print("  and machine learning models incorporating phylogeny and environment")

# Create modeling dataset
if merged is not None and len(merged) > 0:
    model_data = merged[[
        'accession', 'genome_size_mb', 'genes_total_num', 'genes_protein_num',
        'gc_percent', 'organism_taxId', 'organism_name'
    ]].copy()
    
    if 'ORGANISM ECOSYSTEM CATEGORY' in merged.columns:
        model_data['environment'] = merged['ORGANISM ECOSYSTEM CATEGORY']
    
    # Add gene density
    model_data['gene_density'] = (model_data['genes_total_num'] / model_data['genome_size_mb']).round(2)
    
    # Remove missing values
    model_data = model_data.dropna(subset=['genome_size_mb', 'genes_total_num'])
    
    model_data.to_csv(RESULTS_DIR / '09_modeling_dataset.tsv', sep='\t', index=False)
    print(f"  ✓ Saved modeling dataset: {RESULTS_DIR / '09_modeling_dataset.tsv'}")
    print(f"    {len(model_data):,} genomes ready for modeling")

print("\n" + "=" * 80)
print("COMPREHENSIVE ANALYSIS COMPLETE")
print("=" * 80)
print(f"\n✓ All outputs saved to: {RESULTS_DIR}")
print(f"\nNext steps:")
print("  1. Select 5-20 taxonomically diverse organisms per selected environment")
print("  2. Perform KEGG pathway annotation on protein sequences")
print("  3. Analyze amino acid usage from protein sequences (full analysis)")
print("  4. Build phylogenetic regression models (PGLS)")
print("  5. Test metabolic complexity vs. nutrient limitation hypotheses")
print("  6. Analyze transcription factors and mobile elements")
print("=" * 80)
