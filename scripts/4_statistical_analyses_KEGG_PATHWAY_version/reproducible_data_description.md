# Reproducible Data Description for Statistical Analyses

## Overview

This document provides a comprehensive description of all datasets used in the statistical analyses for the bacterial genome constraint project. The goal is to enable complete reproducibility of the scaling law analyses (similar to van Nimwegen 2003) that examine how gene content scales with genome size across different environments and functional categories.

**Project Goal**: Model how environmental factors constrain bacterial genome size by analyzing genome size, metabolic coding capacity (via KEGG), and other genomic signatures across multiple well-defined environments.

**Analysis Framework**: Replicate and extend the scaling law analyses from van Nimwegen (2003) "Scaling laws in the functional content of genomes", but stratified by **environment** (from GOLD metadata) rather than by clade or lifestyle.

---

## Data Sources Summary

| Source | Location | Purpose | Key Linking Field |
|--------|----------|---------|------------------|
| **NCBI Genomes** | `data/ncbi/` | Primary genome sequences and annotations | `accession` (GCF_*) |
| **GOLD Metadata** | `data/gold/` | Environment/ecology annotations | `ORGANISM NCBI TAX ID` |
| **GO Annotations** | `results/3_GO_analyses/` | Functional category counts | `accession` |
| **KEGG Annotations** | `results/3_GO_analyses/2_JGIgold_KEGG_anayses_out/` | Metabolic pathway annotations | `accession` |
| **Exploratory Results** | `results/1_exploratory_analyses_out/` | Pre-processed environment summaries | `accession` |

**Total Genomes Analyzed**: 5,986 NCBI reference genomes → 3,101 with complete annotations → 3,088 with GO annotations

---

## 1. NCBI Genome Data (`data/ncbi/`)

### 1.1 Directory Structure

```
data/ncbi/
├── assemblies/                    # 5,986 genome directories (126 GB)
│   └── GCF_*/                    # One directory per assembly accession
│       ├── *_genomic.fna         # Genome sequences (FASTA)
│       ├── protein.faa           # Protein sequences (FASTA)
│       ├── cds_from_genomic.fna  # Coding sequences (FASTA)
│       ├── genomic.gff           # Gene annotations (GFF3 format)
│       └── genomic.gbff          # GenBank format
│
├── metadata/                      # Summary files and manifests (29 MB)
│   ├── assembly_manifest.txt              # List of all accessions (5,986 rows)
│   ├── metadata_table.tsv                 # Basic metadata (6 columns, 5,987 rows)
│   ├── assembly_data_report_extracted.tsv # Detailed metadata (73 columns, 5,987 rows)
│   └── dataset_catalog_extracted.tsv      # File catalog (4 columns, 35,918 rows)
│
└── bac_genomes_5988_dir/         # Original NCBI download structure (38 MB)
    └── ncbi_dataset/
        └── data/
            ├── assembly_data_report.jsonl  # Complete metadata (JSON Lines, 23 MB)
            └── dataset_catalog.json        # File catalog (JSON, 5.3 MB)
```

### 1.2 Key Metadata Files

#### `metadata/assembly_manifest.txt`
- **Format**: One accession per line (plain text)
- **Rows**: 5,986
- **Content**: List of all assembly accessions (e.g., `GCF_000006945.2`)
- **Usage**: Iterating through all genomes in processing pipelines
- **Example**:
  ```
  GCF_000006945.2
  GCF_000195955.2
  GCF_000009045.1
  ```

#### `metadata/metadata_table.tsv`
- **Format**: Tab-separated values
- **Rows**: 5,987 (header + 5,986 genomes)
- **Columns**: 6
  1. `Assembly Accession` - RefSeq accession (GCF_*)
  2. `Assembly Release Date` - YYYY-MM-DD format
  3. `Assembly Level` - Complete, Chromosome, Scaffold, etc.
  4. `Annotation Name` - Annotation version
  5. `Organism Name` - Scientific name
  6. `Organism Taxonomic ID` - NCBI taxonomy ID
- **Usage**: Basic genome-to-taxonomy mapping, filtering by date/level
- **Key for joining**: `Assembly Accession` → `accession` in other files

#### `metadata/assembly_data_report_extracted.tsv`
- **Format**: Tab-separated values
- **Rows**: 5,987 (header + 5,986 genomes)
- **Columns**: 73 (comprehensive metadata)
- **Size**: 4.8 MB
- **Key Columns for Analysis**:
  - **Identifiers**: `accession`, `organism_taxId`, `organism_name`
  - **Assembly Statistics**: `stats_totalSequenceLength` (genome size in bp), `stats_gcPercent`, `stats_numberOfContigs`
  - **Gene Counts**: `genes_total`, `genes_proteinCoding`, `genes_nonCoding`, `genes_pseudogene`
  - **Quality Metrics**: `checkm_completeness`, `checkm_contamination`
  - **Biosample Info**: `biosample_isolation_source`, `biosample_geo_loc_name`, `biosample_host`
- **Usage**: Comprehensive metadata queries, filtering by quality metrics
- **Note**: Numeric fields stored as strings (convert as needed)

### 1.3 Genome Sequence Files

Each assembly directory (`data/ncbi/assemblies/GCF_*/`) contains:

- **`*_genomic.fna`**: Genomic DNA sequences (FASTA format)
  - Variable naming: `{ACCESSION}_{ASSEMBLY_NAME}_genomic.fna`
  - Use glob pattern `*_genomic.fna` to find files
  - Used for: Genome size calculation, sequence analysis

- **`protein.faa`**: Protein sequences (FASTA format, amino acids)
  - Used for: KEGG annotation, protein domain analysis (Pfam)

- **`cds_from_genomic.fna`**: Coding DNA sequences (FASTA format)
  - Used for: Gene analysis, codon usage, amino acid composition

- **`genomic.gff`**: Gene annotation (GFF3 format)
  - Contains: Gene features, CDS, GO terms in `Ontology_term` attribute
  - Used for: GO term extraction, gene location analysis

- **`genomic.gbff`**: GenBank flat file format
  - Complete annotation with full metadata

**File Access Pattern**:
```python
import glob
acc = 'GCF_000006945.2'
genome_dir = f'data/ncbi/assemblies/{acc}/'
genome_fna = glob.glob(f'{genome_dir}/*_genomic.fna')[0]
protein_faa = f'{genome_dir}/protein.faa'
gff = f'{genome_dir}/genomic.gff'
```

---

## 2. GOLD Metadata (`data/gold/`)

### 2.1 Overview

GOLD (Genomes Online Database) provides comprehensive environment and ecology metadata for bacterial genomes. The metadata is split across multiple CSV files extracted from a master Excel file.

**Total Size**: ~1.1 GB (all files combined)
**Source**: JGI GOLD database (downloaded 2025-11-06)

### 2.2 File Structure

| File | Rows | Size | Purpose |
|------|------|------|---------|
| `0_20251106_gold_metadata_Organism.csv` | 527,569 | 167 MB | **Primary file** - Organism-level metadata with environment |
| `0_20251106_gold_metadata_Analysis Project.csv` | 449,140 | 197 MB | Analysis project metadata |
| `0_20251106_gold_metadata_Sequencing Project.csv` | 615,411 | 395 MB | Sequencing project metadata |
| `0_20251106_gold_metadata_Biosample.csv` | 233,866 | 66 MB | Biosample metadata |
| `0_20251106_gold_metadata_Study.csv` | 63,759 | 11 MB | Study metadata |
| `0_20251106_gold_metadata.xlsx` | - | 224 MB | Master Excel file (all sheets) |

### 2.3 Primary File: `0_20251106_gold_metadata_Organism.csv`

**Format**: CSV (comma-separated, quoted fields with headers)
**Rows**: 527,569 (header + 527,568 organisms)
**Columns**: 50 (exact count)
**Size**: 167 MB
**Encoding**: UTF-8 (with quoted fields containing commas)

#### Key Columns for Environment Stratification (Exact Names):

1. **`ORGANISM NCBI TAX ID`** (integer/string) - **CRITICAL JOIN KEY**: Links to NCBI genomes via `organism_taxId` column
   - Format: NCBI taxonomy ID (integer, e.g., `224914`)
   - Usage: Primary join key: `GOLD['ORGANISM NCBI TAX ID']` = `NCBI['organism_taxId']`

2. **`ORGANISM ECOSYSTEM CATEGORY`** (string) - **PRIMARY ENVIRONMENT VARIABLE**
   - **Top values** (from GOLD sample of 10,000 organisms): 
     - "Mammals: Human" (4,607), "Aquatic" (1,105), "Terrestrial" (512), "Plants" (440), "Mammals" (416), "Food production" (200), "Arthropoda: Insects" (124), "Microbial" (121), "Wastewater" (103), "Lab synthesis" (74), "Unclassified" (65), "Birds" (58), "Fish" (53), "Mollusca" (32), "Algae" (26), "Solid waste" (23), "Arthropoda: Chelicerates" (23), "Invertebrates" (18), "Industrial production" (17), "Fungi" (16)
   - This is the **main stratification variable** for scaling law analyses
   - **Note**: Values may include colons (e.g., "Mammals: Human") - handle in string matching
   - **Note**: Some values may be empty strings - filter these out before analysis

3. **`ORGANISM ECOSYSTEM`** (string) - More specific ecosystem classification
4. **`ORGANISM ECOSYSTEM TYPE`** (string) - Ecosystem subtype
5. **`ORGANISM ECOSYSTEM SUBTYPE`** (string) - Further refinement
6. **`ORGANISM SPECIFIC ECOSYSTEM`** (string) - Most specific classification

#### Additional Environment-Relevant Columns:

- **`ORGANISM HABITAT`** - Free text description (e.g., "Aquatic|Hot spring|Fresh water")
- **`ORGANISM OXYGEN REQUIREMENT`** - "Aerobe", "Anaerobe", "Facultative", "Obligate aerobe", "Obligate anaerobe"
- **`ORGANISM METABOLISM`** - Metabolic type (e.g., "Methanogen", "Cellulose degrader")
- **`ORGANISM TEMPERATURE RANGE`** - "Mesophile", "Thermophile", "Hyperthermophile", "Psychrophile"
- **`ORGANISM SALINITY`** - Salinity preference
- **`ORGANISM ISOLATION HOST NAME`** - Host organism (for host-associated)
- **`ORGANISM HOST BODY SITE`** - Body site (for host-associated)
- **`ORGANISM GEOGRAPHIC LOCATION`** - Geographic origin
- **`ORGANISM LATITUDE`**, **`ORGANISM LONGITUDE`** - Geographic coordinates

#### Other Important Columns:

- **`ORGANISM GOLD ID`** - GOLD-specific identifier (e.g., "Go0000002")
- **`ORGANISM NAME`** - Organism name (may differ from NCBI name)
- **`ORGANISM GENUS`**, **`ORGANISM SPECIES`** - Taxonomic classification
- **`ORGANISM BIOTIC RELATIONSHIPS`** - "Free living", "Host-associated", etc.

#### Usage for Joining:

```python
import pandas as pd

# Load GOLD metadata
gold = pd.read_csv('data/gold/0_20251106_gold_metadata_Organism.csv')

# Load NCBI metadata
ncbi = pd.read_csv('data/ncbi/metadata/assembly_data_report_extracted.tsv', sep='\t')

# Join on taxonomy ID
merged = ncbi.merge(
    gold,
    left_on='organism_taxId',
    right_on='ORGANISM NCBI TAX ID',
    how='inner'
)

# Filter by environment
aquatic_genomes = merged[merged['ORGANISM ECOSYSTEM CATEGORY'] == 'Aquatic']
```

**Note**: Not all NCBI genomes have GOLD metadata. The exploratory analysis found 5,770 genomes successfully merged with GOLD data.

---

## 3. GO Analysis Results (`results/3_GO_analyses/`)

### 3.1 Overview

Gene Ontology (GO) annotations extracted from GFF files for functional category analysis. These files enable scaling law analyses similar to van Nimwegen (2003).

**Total Genomes with GO Annotations**: 3,088
**Genomes Excluded**: 13 (no GO terms found)

### 3.2 Key Files

#### `presence_absence_table.txt`
- **Format**: Tab-separated values (TSV)
- **Rows**: 3,089 (header + 3,088 genomes)
- **Columns**: 3,834 (1 `Genome` column + 3,833 unique GO term columns)
- **Structure**:
  - First column: `Genome` (string, accession in GCF_* format, e.g., `GCF_000006985.1`)
  - Remaining 3,833 columns: GO term IDs in 7-digit zero-padded format (e.g., `0000009`, `0000010`, `0000014`, `0000015`)
  - Values: Integer `1` if term present, `0` if absent
- **Data Types**: 
  - `Genome` column: String
  - All GO term columns: Integer (0 or 1)
- **GO Term ID Format**: 7-digit zero-padded (e.g., `0000015` = `GO:0000015`)
  - To convert to standard GO format: Prepend `GO:` to each term ID
- **Usage**: Binary presence/absence analysis, identifying co-occurring terms, filtering genomes by functional capability
- **Example Row** (first 10 columns):
  ```
  Genome	0000009	0000010	0000014	0000015	0000023	0000027	0000028	0000030	0000034
  GCF_000006985.1	0	0	0	1	0	0	0	0	0
  ```

#### `ubiquitous_terms.txt`
- **Format**: One GO term ID per line (plain text)
- **Rows**: 334 (334 ubiquitous terms)
- **Content**: GO terms present in ≥95% of genomes
- **Usage**: Filtering for core functional categories (equivalent to van Nimwegen's 129 categories)
- **Example**:
  ```
  0000015
  0000049
  0000155
  0000156
  0000160
  ```

#### `unique_terms.txt`
- **Format**: One GO term ID per line (plain text)
- **Rows**: 3,833 (all unique GO terms found)
- **Content**: Complete list of all GO term IDs found across all genomes
- **Usage**: Reference list for mapping term IDs to names

#### `ubiquitous_counts_table.txt`
- **Format**: Tab-separated values (TSV)
- **Rows**: 3,089 (header + 3,088 genomes)
- **Columns**: 335 (1 `Genome` column + 334 ubiquitous GO term columns)
- **Structure**:
  - First column: `Genome` (string, accession in GCF_* format, index column)
  - Remaining 334 columns: GO term IDs in 7-digit format (e.g., `0000015`, `0000049`, `0000155`)
  - Values: Integer counts (non-negative) of how many times each term appears in each genome
- **Data Types**:
  - `Genome` column: String (use as index when loading)
  - All GO term columns: Integer (counts, can be 0)
- **Value Ranges**: 
  - Minimum: 0 (term absent or count = 0)
  - Maximum: Varies by category (sample shows counts up to at least 3 for first category)
  - Typical range: 0 to hundreds/thousands depending on category and genome size
- **GO Term Selection**: Only includes 334 terms present in ≥95% of genomes (from `ubiquitous_terms.txt`)
- **Usage**: Quantitative analysis of ubiquitous term abundance, **PRIMARY DATA for scaling law analyses**
- **Key for Scaling Laws**: This provides **nc(g)** values - the count of genes/domains in category c for genome g
- **Example Row** (first 10 columns):
  ```
  Genome	0000015	0000049	0000155	0000156	0000160	0000166	0000179	0000271	0000287
  GCF_000006985.1	2	0	1	4	6	19	1	2	9
  ```

#### `ubiquitous_counts_table_normalized.txt`
- **Format**: Tab-separated values (TSV)
- **Rows**: 3,089 (header + 3,088 genomes)
- **Columns**: 335 (1 `Genome` column + 334 ubiquitous GO term columns)
- **Structure**:
  - First column: `Genome` (string, accession in GCF_* format)
  - Remaining 334 columns: GO term IDs in 7-digit format (same as `ubiquitous_counts_table.txt`)
  - Values: Floating-point numbers (normalized counts)
- **Data Types**:
  - `Genome` column: String
  - All GO term columns: Float (normalized proportions)
- **Normalization Formula**: `normalized_count = raw_count / genes_proteinCoding`
  - Where `genes_proteinCoding` comes from the genome's protein-coding gene count
  - Result: Proportion of protein-coding genes annotated with each term
- **Value Ranges**: 
  - Minimum: 0.0 (term absent or count = 0)
  - Maximum: Typically < 1.0 (proportion of genes)
  - Typical range: 0.0 to ~0.01-0.1 depending on category (e.g., 0.001006 = 0.1% of genes)
- **Usage**: Size-normalized comparisons across genomes of different sizes
- **Example Row** (first 5 columns):
  ```
  Genome	0000015	0000049	0000155	0000156
  GCF_000006985.1	0.001006036217303823	0.0	0.0005030181086519115	0.002012072434607646
  ```

#### `missing_terms.txt`
- **Format**: One filename per line
- **Rows**: 13
- **Content**: Genomes with no GO annotations (excluded from analyses)

### 3.3 GO Ontology File

**Location**: `data/go/go-basic.obo`
- **Format**: OBO (Open Biomedical Ontology)
- **Size**: ~550,491 lines
- **Content**: Complete GO ontology hierarchy with term definitions, relationships, names
- **Usage**: Mapping GO term IDs to names, understanding term relationships
- **Helper Script**: `scripts/3_GO_analyses/load_go_ontology.py` provides functions to load ontology data

### 3.4 Usage for Scaling Law Analysis

```python
import pandas as pd
import numpy as np

# Load ubiquitous term counts (this is nc(g) for each category c)
counts = pd.read_csv('results/3_GO_analyses/ubiquitous_counts_table.txt', sep='\t', index_col=0)

# Load genome metadata to get total gene count n(g)
metadata = pd.read_csv('results/3_GO_analyses/2_JGIgold_KEGG_anayses_out/01_genome_metrics.tsv', sep='\t')

# For each GO category c, fit: nc(g) = βc * n(g)^αc
# where n(g) is total genes and nc(g) is count for category c

for category in counts.columns:
    n_g = metadata['genes_total'].values  # Total genes per genome
    nc_g = counts[category].values  # Category count per genome
    
    # Fit power law: log(nc_g) = log(βc) + αc * log(n_g)
    # Use linear regression on log-transformed data
    log_n = np.log(n_g)
    log_nc = np.log(nc_g + 1)  # Add 1 to avoid log(0)
    
    # Fit: log_nc = β + α * log_n
    alpha, beta = np.polyfit(log_n, log_nc, 1)
    beta = np.exp(beta)  # Convert back from log space
```

---

## 4. KEGG and Genome Feature Analysis (`results/3_GO_analyses/2_JGIgold_KEGG_anayses_out/`)

### 4.1 Overview

Comprehensive genome feature matrices combining NCBI metadata, GO annotations, KEGG annotations, transcription factors, mobile elements, and amino acid composition.

**Total Genomes**: 3,101 with complete metrics
**Genomes with KEGG**: 0 (KEGG annotation not yet completed)

### 4.2 Key Files

#### `01_genome_metrics.tsv`
- **Format**: Tab-separated values (TSV)
- **Rows**: 3,102 (header + 3,101 genomes)
- **Columns**: 8 (exact order)
  1. `accession` (string) - Assembly accession (GCF_* format)
  2. `genome_size_bp` (integer) - Genome size in base pairs
  3. `gc_percent` (float) - GC content percentage
  4. `genes_total` (integer) - Total number of genes
  5. `genes_proteinCoding` (integer) - Number of protein-coding genes
  6. `checkm_completeness` (float) - CheckM completeness score (0-100)
  7. `checkm_contamination` (float) - CheckM contamination score (typically 0-10)
  8. `coding_density` (float) - Percentage of genome that is coding sequence
- **Data Types**: 
  - Integers: `genome_size_bp`, `genes_total`, `genes_proteinCoding`
  - Floats: `gc_percent`, `checkm_completeness`, `checkm_contamination`, `coding_density`
- **Value Ranges**:
  - Genome size: 1.00 - 13.61 Mb (1,003,404 - 13,605,706 bp)
  - GC content: 22.5% - 76.0%
  - Total genes: 633 - 11,093
  - CheckM completeness: Typically 90-100% (1 missing value)
  - CheckM contamination: Typically 0-10% (363 missing values)
- **Missing Data**: 
  - `checkm_completeness`: 1 missing value
  - `checkm_contamination`: 363 missing values (empty strings, not NaN)
- **Usage**: Basic genome statistics, filtering by quality metrics
- **Key for Scaling Laws**: `genes_total` is **n(g)** - the total gene count used in scaling law fits
- **Example Row**:
  ```
  GCF_000007125.1	3294931	57	3234	2941	97.2	0.58	87.6825341714288
  ```

#### `03_tf_mobile_elements.tsv`
- **Format**: Tab-separated values (TSV)
- **Rows**: 3,102 (header + 3,101 genomes)
- **Columns**: 3 (exact order)
  1. `accession` (string) - Assembly accession (GCF_* format)
  2. `tf_count` (integer) - Number of transcription factors
  3. `mobile_element_count` (integer) - Number of mobile genetic elements
- **Data Types**: All integers, no missing values
- **Value Ranges** (from sample data):
  - Transcription factors: 28 - 924 (wide range, scales with genome size)
  - Mobile elements: 935 - 7,795 (includes transposons, IS elements, etc.)
- **Usage**: Regulatory complexity analysis, mobile element analysis
- **Note**: Transcription factor count is a key category in van Nimwegen (2003) - shows quadratic scaling (α ≈ 2)
- **Example Row**:
  ```
  GCF_000007125.1	218	3183
  ```

#### `04_kegg_annotation.tsv`
- **Format**: Tab-separated values
- **Rows**: 11 (header + 10 genomes, mostly empty)
- **Columns**: 4
  1. `accession` - Assembly accession
  2. `ko_count` - Number of KEGG orthologs (KO)
  3. `module_count` - Number of KEGG modules
  4. `ko_per_mb` - KO density per megabase
- **Status**: **NOT YET COMPLETED** - KEGG annotation in progress
- **Usage**: Metabolic pathway analysis (when complete)

#### `05_genome_feature_matrix.tsv`
- **Format**: Tab-separated values (TSV)
- **Rows**: 3,102 (header + 3,101 genomes)
- **Columns**: 21 (exact order, comprehensive feature matrix)
  1. `accession` (string) - Assembly accession (GCF_* format)
  2. `organism_name` (string) - Organism scientific name
  3. `organism_taxId` (integer) - NCBI taxonomy ID
  4. `release_date` (string) - Assembly release date (YYYY-MM-DD format)
  5. `assembly_level` (string) - Assembly completeness level (e.g., "Complete Genome")
  6. `genome_size_bp` (integer) - Genome size in base pairs
  7. `gc_percent` (float) - GC content percentage
  8. `genes_total` (integer) - Total genes (**n(g) for scaling laws**)
  9. `genes_proteinCoding` (integer) - Protein-coding genes
  10. `checkm_completeness` (float) - CheckM completeness score (0-100)
  11. `checkm_contamination` (float) - CheckM contamination score
  12. `coding_density` (float) - Coding sequence percentage
  13. `amino_n_burden` (float) - Amino acid nitrogen burden (nutrient limitation signature)
  14. `tf_count` (integer) - Transcription factor count
  15. `mobile_element_count` (integer) - Mobile element count
  16. `ko_count` (float) - KEGG ortholog count (**ALL MISSING** - 3,101 empty)
  17. `module_count` (float) - KEGG module count (**ALL MISSING** - 3,101 empty)
  18. `ko_per_mb` (float) - KO density per megabase (**ALL MISSING** - 3,101 empty)
  19. `environment` (string) - **ENVIRONMENT CATEGORY** (from GOLD, **PRIMARY STRATIFICATION VARIABLE**)
  20. `genome_size_mb` (float) - Genome size in megabases (derived from `genome_size_bp`)
  21. `module_per_mb` (float) - Module density per megabase (**ALL MISSING** - 3,101 empty)
- **Data Types**: Mixed (strings, integers, floats)
- **Missing Data Pattern**:
  - KEGG fields (`ko_count`, `module_count`, `ko_per_mb`, `module_per_mb`): **100% missing** (3,101/3,101)
  - `checkm_contamination`: 363 missing (11.7%)
  - `checkm_completeness`: 1 missing (0.03%)
  - All other fields: Complete
- **Environment Distribution** (12 environments with data):
  - Aquatic: 841 genomes (27.1%)
  - Terrestrial: 645 genomes (20.8%)
  - Mammals: Human: 528 genomes (17.0%)
  - Plants: 350 genomes (11.3%)
  - Mammals: 274 genomes (8.8%)
  - Food production: 146 genomes (4.7%)
  - Unclassified: 102 genomes (3.3%)
  - Wastewater: 93 genomes (3.0%)
  - Birds: 52 genomes (1.7%)
  - Solid waste: 24 genomes (0.8%)
  - Algae: 24 genomes (0.8%)
  - Bioreactor: 22 genomes (0.7%)
- **Usage**: **PRIMARY DATA FILE** for statistical analyses - combines all features
- **Key for Environment Stratification**: `environment` column contains GOLD ecosystem categories
- **Example Row**:
  ```
  GCF_000007125.1	Brucella melitensis bv. 1 str. 16M	224914	2001-12-28	Complete Genome	3294931	57.0	3234	2941	97.2	0.58	87.68253417142878	1.356950352751071	218	3183				Mammals	3.294931	
  ```

#### `06_summary_statistics.tsv`
- **Format**: Tab-separated values
- **Rows**: 2 (header + 1 summary row)
- **Columns**: 6
  - `total_genomes` - 3,101
  - `genomes_with_metrics` - 3,101
  - `genomes_with_aa` - 3,101 (amino acid composition)
  - `genomes_with_tf` - 3,101 (transcription factors)
  - `genomes_with_mobile` - 3,101 (mobile elements)
  - `genomes_with_kegg` - 0 (KEGG not yet completed)

### 4.3 Amino Acid Composition File

#### `02_amino_acid_composition.tsv`
- **Format**: Tab-separated values (TSV)
- **Rows**: 62,856 (header + 62,855 entries)
- **Structure**: Long format with one row per genome-amino acid combination
- **Columns**: 6 (exact order)
  1. `accession` (string) - Assembly accession (GCF_* format)
  2. `amino_n_burden` (float) - Amino acid nitrogen burden (calculated metric)
  3. `total_aa_residues` (integer) - Total amino acid residues in proteome
  4. `aa_Var1` (string) - First amino acid variable (appears to be "A" in all rows)
  5. `aa_Var2` (string) - Second amino acid variable (actual amino acid code: A, C, D, E, F, G, H, I, K, L, M, N, P, Q, R, S, T, V, W, Y)
  6. `aa_Freq` (float) - Frequency of this amino acid in the proteome (percentage)
- **Data Types**: Mixed (strings, integers, floats)
- **Structure Details**:
  - Each genome has 20 rows (one per standard amino acid)
  - `aa_Var1` appears to be a constant ("A" in all observed rows)
  - `aa_Var2` contains the actual amino acid single-letter code
  - `aa_Freq` is the percentage frequency (e.g., 11.58 = 11.58%)
- **Value Ranges** (from sample):
  - `amino_n_burden`: ~1.31 - 1.41 (nutrient limitation signature)
  - `total_aa_residues`: Varies by genome size (e.g., 877,531 for GCF_000007125.1)
  - `aa_Freq`: 0.0 - ~15% (varies by amino acid and genome)
- **Usage**: Amino acid usage analysis, nitrogen burden calculation
- **Note**: This data is used to calculate `amino_n_burden` in the feature matrix (`05_genome_feature_matrix.tsv`)
- **Example Rows** (same genome, different amino acids):
  ```
  accession	amino_n_burden	total_aa_residues	aa_Var1	aa_Var2	aa_Freq
  GCF_000007125.1	1.35695035275107	877531	A	A	11.5804455910959
  GCF_000007125.1	1.35695035275107	877531	A	C	0.839400545393838
  GCF_000007125.1	1.35695035275107	877531	A	D	5.54818006429403
  ```

---

## 5. Exploratory Analysis Results (`results/1_exploratory_analyses_out/`)

### 5.1 Overview

Pre-processed datasets from initial exploratory analyses, ready for downstream statistical modeling.

### 5.2 Key Files

#### `08_selected_environments_summary.tsv`
- **Format**: Tab-separated values (TSV)
- **Rows**: 3,102 (header + 3,101 genomes)
- **Columns**: 8 (exact order)
  1. `accession` (string) - Assembly accession (GCF_* format)
  2. `genome_size_mb` (float) - Genome size in megabases
  3. `genes_total_num` (integer) - Total gene count
  4. `genes_protein_num` (integer) - Protein-coding gene count
  5. `gc_percent` (float) - GC content percentage
  6. `organism_name` (string) - Organism scientific name
  7. `organism_taxId` (integer) - NCBI taxonomy ID
  8. `ORGANISM ECOSYSTEM CATEGORY` (string) - **PRIMARY ENVIRONMENT VARIABLE** (from GOLD)
- **Data Types**: All columns present, no missing values in key fields
- **Usage**: Environment-stratified genome size analysis
- **Key Finding**: 16 environments identified with ≥20 genomes suitable for detailed analysis
- **Example Row**:
  ```
  GCF_000007125.1	3.294931	3234	2941	57.0	Brucella melitensis bv. 1 str. 16M	224914	Mammals
  ```

#### `09_modeling_dataset.tsv`
- **Format**: Tab-separated values (TSV)
- **Rows**: 5,331 (header + 5,330 genomes)
- **Columns**: 9 (exact order)
  1. `accession` (string) - Assembly accession (GCF_* format)
  2. `genome_size_mb` (float) - Genome size in megabases
  3. `genes_total_num` (integer) - Total gene count
  4. `genes_protein_num` (integer) - Protein-coding gene count
  5. `gc_percent` (float) - GC content percentage
  6. `organism_taxId` (integer) - NCBI taxonomy ID
  7. `organism_name` (string) - Organism scientific name
  8. `environment` (string) - **ENVIRONMENT CATEGORY** (simplified/cleaned from GOLD)
  9. `gene_density` (float) - Genes per megabase (calculated: genes_total_num / genome_size_mb)
- **Data Types**: All columns present, no missing values
- **Usage**: **PRIMARY MODELING DATASET** - ready for regression/ML modeling
- **Environment Distribution**: 39 unique environments
  - Top environments: Aquatic (841), Terrestrial (645), Mammals: Human (528), Plants (350), Mammals (274), Food production (146)
  - Minimum for analysis: 20 genomes (Invertebrates, Air, Porifera, etc.)
- **Note**: More genomes than GO analysis (5,330 vs 3,088) because this includes all genomes with GOLD metadata, not just those with GO annotations
- **Example Row**:
  ```
  GCF_000007125.1	3.294931	3234	2941	57.0	224914	Brucella melitensis bv. 1 str. 16M	Mammals	981.51
  ```

---

## 6. Data Mapping and Joining Strategy

### 6.1 Primary Join Keys

| Dataset | Join Key | Format | Example |
|---------|----------|--------|---------|
| NCBI Metadata | `accession` | GCF_* | `GCF_000006945.2` |
| GOLD Metadata | `ORGANISM NCBI TAX ID` | Integer | `224914` |
| GO Counts | `Genome` (or `accession`) | GCF_* | `GCF_000006985.1` |
| Feature Matrix | `accession` | GCF_* | `GCF_000007125.1` |

### 6.2 Complete Data Pipeline

```python
import pandas as pd
import numpy as np
from scipy import stats

# Step 1: Load genome feature matrix (has environment and n(g))
features = pd.read_csv(
    'results/3_GO_analyses/2_JGIgold_KEGG_anayses_out/05_genome_feature_matrix.tsv',
    sep='\t',
    na_values=['']  # Handle empty strings as NaN
)

# Step 2: Load GO term counts for ubiquitous categories (nc(g) values)
go_counts = pd.read_csv(
    'results/3_GO_analyses/ubiquitous_counts_table.txt',
    sep='\t',
    index_col=0  # 'Genome' column becomes index
)

# Step 3: Join on accession
merged = features.merge(
    go_counts,
    left_on='accession',
    right_index=True,
    how='inner'
)

# Step 4: Filter for high quality and sufficient sample size
merged = merged[
    (merged['checkm_completeness'] > 90) & 
    (merged['checkm_contamination'] < 5)
]

# Get environments with sufficient genomes (≥20)
env_counts = merged['environment'].value_counts()
valid_envs = env_counts[env_counts >= 20].index
merged = merged[merged['environment'].isin(valid_envs)]

# Step 5: Fit scaling laws for each environment × category
results = []
environments = merged['environment'].unique()
go_categories = go_counts.columns

for env in environments:
    env_data = merged[merged['environment'] == env]
    
    # Total genes n(g) for this environment
    n_g = env_data['genes_total'].values
    
    for category in go_categories:
        # Category counts nc(g) for this environment
        nc_g = env_data[category].values
        
        # Filter out zeros for log transformation
        mask = (n_g > 0) & (nc_g > 0)
        if mask.sum() < 10:  # Need at least 10 data points
            continue
            
        log_n = np.log(n_g[mask])
        log_nc = np.log(nc_g[mask])
        
        # Fit power law: log(nc) = log(β) + α * log(n)
        slope, intercept, r_value, p_value, std_err = stats.linregress(log_n, log_nc)
        
        alpha = slope  # Scaling exponent
        beta = np.exp(intercept)  # Scaling offset
        
        results.append({
            'environment': env,
            'category': category,
            'alpha': alpha,
            'beta': beta,
            'r_squared': r_value**2,
            'p_value': p_value,
            'n_genomes': mask.sum()
        })

results_df = pd.DataFrame(results)
```

### 6.3 Environment Stratification

**Primary Environment Variable**: `environment` column in feature matrix (derived from `ORGANISM ECOSYSTEM CATEGORY` in GOLD)

**Available Environments in Feature Matrix** (3,101 genomes, 12 environments with data):
**Exact environment values** (alphabetically sorted):
1. **Algae**: 24 genomes (0.8%)
2. **Aquatic**: 841 genomes (27.1%) - Largest group
3. **Bioreactor**: 22 genomes (0.7%)
4. **Birds**: 52 genomes (1.7%)
5. **Food production**: 146 genomes (4.7%)
6. **Mammals**: 274 genomes (8.8%)
7. **Mammals: Human**: 528 genomes (17.0%)
8. **Plants**: 350 genomes (11.3%)
9. **Solid waste**: 24 genomes (0.8%)
10. **Terrestrial**: 645 genomes (20.8%)
11. **Unclassified**: 102 genomes (3.3%)
12. **Wastewater**: 93 genomes (3.0%)

**Note**: Environment names are case-sensitive and use exact spelling as shown above. The colon in "Mammals: Human" is part of the value.

**Available Environments in Modeling Dataset** (5,330 genomes, 39 unique environments):
- **All 12 from feature matrix** (see above)
- **Additional environments** (27 more):
  - Arthropoda: Insects (89 genomes)
  - Fish (44 genomes)
  - Mollusca (35 genomes)
  - Invertebrates (20 genomes)
  - Air (19 genomes)
  - Porifera (17 genomes)
  - Lab synthesis (13 genomes)
  - Built environment (12 genomes)
  - Arthropoda: Chelicerates (11 genomes)
  - Microbial (11 genomes)
  - Industrial production (11 genomes)
  - Lab culture (10 genomes)
  - Arthropoda: Crustaceans (10 genomes)
  - Biotransformation (9 genomes)
  - Fungi (5 genomes)
  - Annelida (4 genomes)
  - Tunicates (4 genomes)
  - Sewage treatment plant (3 genomes)
  - Amphibia (3 genomes)
  - Lab enrichment (3 genomes)
  - Artificial ecosystem (3 genomes)
  - Reptilia (3 genomes)
  - Bioremediation (3 genomes)
  - Animal feed production (2 genomes)
  - Laboratory developed (2 genomes)
  - Protists (2 genomes)
  - WWTP (2 genomes)

**Environment Selection Criteria** (from exploratory analysis):
- **Minimum 20 genomes per environment** for statistical power in scaling law fits
- **16+ environments meet this criterion** in modeling dataset
- **Extreme environments identified**: 
  - Smallest average: Insects (2.44 Mb avg)
  - Largest average: Terrestrial (5.73 Mb avg)
- **Recommended for analysis**: Environments with ≥20 genomes (sufficient for regression)

---

## 7. Scaling Law Analysis Framework

### 7.1 Mathematical Framework

Following van Nimwegen (2003), for each functional category c and each environment e:

**Power Law**: nc,e(g) = βc,e * n(g)^αc,e

Where:
- `nc,e(g)` = number of genes/domains in category c for genome g in environment e
- `n(g)` = total number of genes in genome g
- `αc,e` = scaling exponent (category c, environment e)
- `βc,e` = scaling offset (category c, environment e)

**Log-Linear Form**: log(nc,e(g)) = log(βc,e) + αc,e * log(n(g))

### 7.2 Data Requirements

For each environment e and category c, need:
1. **n(g)**: Total gene count → `genes_total` from feature matrix
2. **nc(g)**: Category count → `ubiquitous_counts_table.txt` columns
3. **Environment label**: `environment` column from feature matrix

### 7.3 Statistical Analysis Steps

1. **Filter by Environment**: Subset genomes by `environment` column
2. **For Each GO Category**: Extract category counts from `ubiquitous_counts_table.txt`
3. **Fit Power Law**: Use linear regression on log-transformed data
   - Dependent: log(nc(g) + 1)
   - Independent: log(n(g))
   - Slope = αc,e (exponent)
   - Intercept = log(βc,e) → βc,e = exp(intercept)
4. **Calculate Z-Statistics**: Compare exponents across environments
   - Z = (αc,e1 - αc,e2) / sqrt(σ²c,e1 + σ²c,e2)
5. **Identify Variable Categories**: Categories with Z > 2 across environments

### 7.4 Expected Results (Based on van Nimwegen 2003)

- **Most categories**: Z < 2 (statistically indistinguishable across environments)
- **Transcription regulators**: Should show quadratic scaling (α ≈ 2) across most environments
- **Metabolic processes**: Should show near-linear scaling (α ≈ 1)
- **Variable categories**: Small set (like amino acid metabolism) may show significant variation

---

## 8. Data Quality and Completeness

### 8.1 Genome Coverage

| Stage | Count | Description | Percentage of NCBI |
|------|-------|-------------|-------------------|
| **NCBI Download** | 5,986 | Complete reference genomes | 100% |
| **GOLD Merge** | 5,770 | Genomes with environment metadata | 96.4% |
| **GO Annotations** | 3,088 | Genomes with GO term annotations | 51.6% |
| **Complete Metrics** | 3,101 | Genomes with all feature metrics | 51.8% |
| **KEGG Annotations** | 0 | Not yet completed | 0% |
| **Modeling Dataset** | 5,330 | Genomes with environment + basic metrics | 89.0% |

**Key Filtering Steps**:
1. **NCBI → GOLD**: 216 genomes lost (3.6%) - no GOLD metadata available
2. **GOLD → GO**: 2,682 genomes lost (46.5%) - no GO annotations in GFF files
3. **GO → Complete Metrics**: 13 genomes lost (0.4%) - missing quality metrics

### 8.2 Missing Data Handling

**Complete Missing (100% missing)**:
- **KEGG Annotations**: All 4 KEGG-related fields (`ko_count`, `module_count`, `ko_per_mb`, `module_per_mb`) are **100% empty** (3,101/3,101 genomes) - metabolic pathway analysis pending

**Partial Missing**:
- **CheckM Contamination**: 363 missing values (11.7% of 3,101 genomes) - empty strings, not NaN
- **CheckM Completeness**: 1 missing value (0.03% of 3,101 genomes)

**Genome-Level Missing**:
- **GO Annotations**: 13 genomes missing (excluded from GO analyses) - listed in `missing_terms.txt`
- **Environment Metadata**: 216 genomes without GOLD metadata (excluded from environment-stratified analyses)

**Data Type Notes**:
- Empty fields in TSV files are **empty strings** (`""`), not NaN
- When loading with pandas, use `na_values=['']` or `keep_default_na=False` to handle properly
- GOLD CSV files use quoted empty strings - pandas handles these automatically

### 8.3 Quality Filters Applied

**Pre-Download Filters** (NCBI):
- **Assembly Level**: Complete genomes only (from NCBI filter criteria)
- **Reference Quality**: Reference genomes only (not MAGs, not multi-isolate)
- **Annotation**: Must have RefSeq and GenBank annotations
- **Release Date**: Released since 2010

**Post-Download Filters** (Applied in analyses):
- **CheckM Completeness**: >90% recommended (most genomes meet this - range: 90-100%)
- **CheckM Contamination**: <5% recommended (most genomes meet this - range: 0-10%)
- **GO Coverage**: Only genomes with GO annotations included (3,088 of 3,101)
- **Environment Assignment**: Only genomes with GOLD metadata included (5,770 of 5,986)

**Recommended Filters for Scaling Law Analysis**:
```python
# Filter for high-quality genomes
high_quality = data[
    (data['checkm_completeness'] > 90) & 
    (data['checkm_contamination'] < 5) &
    (data['genes_total'] > 0)
]

# Filter for environments with sufficient genomes
env_counts = data['environment'].value_counts()
valid_envs = env_counts[env_counts >= 20].index
filtered_data = data[data['environment'].isin(valid_envs)]
```

---

## 9. File Size Summary

| Dataset | Size | Rows | Columns |
|---------|------|------|---------|
| NCBI assemblies/ | 126 GB | 5,986 genomes | - |
| NCBI metadata/ | 29 MB | 5,986-35,918 | 4-73 |
| GOLD Organism.csv | 167 MB | 527,569 | 50 |
| GOLD (all files) | 1.1 GB | ~2.7M total | - |
| GO presence_absence | ~12 MB | 3,088 | 3,834 |
| GO ubiquitous_counts | ~1 MB | 3,088 | 334 |
| Feature matrix | ~500 KB | 3,101 | 21 |
| Modeling dataset | ~400 KB | 5,330 | 9 |

---

## 10. Reproducibility Checklist

To reproduce the scaling law analyses, ensure you have:

- [ ] **NCBI Genome Data**: `data/ncbi/assemblies/` with all 5,986 genomes
- [ ] **GOLD Metadata**: `data/gold/0_20251106_gold_metadata_Organism.csv`
- [ ] **GO Counts**: `results/3_GO_analyses/ubiquitous_counts_table.txt`
- [ ] **Feature Matrix**: `results/3_GO_analyses/2_JGIgold_KEGG_anayses_out/05_genome_feature_matrix.tsv`
- [ ] **Environment Labels**: Available in feature matrix `environment` column
- [ ] **Ubiquitous Terms List**: `results/3_GO_analyses/ubiquitous_terms.txt` (334 terms)

### Required Python Packages:
- pandas, numpy, scipy, matplotlib, seaborn
- scikit-learn (for regression)
- statsmodels (for statistical tests)

### Analysis Steps:
1. Load feature matrix and GO counts
2. Join on `accession`
3. Filter by `environment` for each environment category
4. For each environment × GO category combination:
   - Extract n(g) and nc(g)
   - Fit power law on log-transformed data
   - Calculate exponent α and offset β
5. Compare exponents across environments using Z-statistics
6. Visualize results (similar to van Nimwegen Figure 1)

---

## 11. Key Differences from van Nimwegen (2003)

| Aspect | van Nimwegen (2003) | This Analysis |
|--------|-------------------|---------------|
| **Stratification** | 13 prokaryotic clades | Environment categories (from GOLD) |
| **Genomes** | 682 | 3,088 (with GO) or 5,330 (with environment) |
| **Categories** | 129 ubiquitous GO categories | 334 ubiquitous GO categories |
| **Annotation Method** | Pfam domains → GO mapping | Direct GO terms from GFF files |
| **Lifestyle Classes** | 24 (habitat, oxygen, temperature, etc.) | Environment categories (ecosystem-based) |
| **Expected Variation** | Minimal across clades | May show more variation across environments |

---

## 12. Contact and References

**Project**: Constraints on Bacterial Genome Size (HST.508 Final Project, Fall 2025)
**Authors**: Aidan Pavao, Terry Cho

**Key Reference**: van Nimwegen, E. (2003). Scaling laws in the functional content of genomes. *Trends in Genetics*, 19(9), 479-484.

**Data Sources**:
- NCBI RefSeq: https://www.ncbi.nlm.nih.gov/refseq/
- GOLD Database: https://gold.jgi.doe.gov/
- Gene Ontology: http://geneontology.org/

---

## Appendix: Quick Reference

### Most Important Files for Scaling Law Analysis:

1. **`results/3_GO_analyses/ubiquitous_counts_table.txt`** - GO category counts **nc(g)** (3,088 genomes × 334 categories)
2. **`results/3_GO_analyses/2_JGIgold_KEGG_anayses_out/05_genome_feature_matrix.tsv`** - Total genes **n(g)** and environment labels (3,101 genomes × 21 features)
3. **`results/3_GO_analyses/ubiquitous_terms.txt`** - List of 334 category IDs to analyze

### Complete Join and Analysis Command:
```python
import pandas as pd
import numpy as np
from scipy import stats

# Load data
features = pd.read_csv(
    'results/3_GO_analyses/2_JGIgold_KEGG_anayses_out/05_genome_feature_matrix.tsv',
    sep='\t',
    na_values=['']
)
go_counts = pd.read_csv(
    'results/3_GO_analyses/ubiquitous_counts_table.txt',
    sep='\t',
    index_col=0
)

# Join on accession
data = features.merge(go_counts, left_on='accession', right_index=True, how='inner')

# Filter for quality and sample size
data = data[(data['checkm_completeness'] > 90) & (data['checkm_contamination'] < 5)]
env_counts = data['environment'].value_counts()
valid_envs = env_counts[env_counts >= 20].index
data = data[data['environment'].isin(valid_envs)]

# Now ready for scaling law analysis
# n(g) = data['genes_total']
# nc(g) = data[category] for each category in go_counts.columns
```

### Environment Filtering:
```python
# Get valid environments (≥20 genomes)
valid_envs = data['environment'].value_counts()
valid_envs = valid_envs[valid_envs >= 20].index

# Filter data
data_filtered = data[data['environment'].isin(valid_envs)]

# Iterate by environment
for env in valid_envs:
    env_data = data_filtered[data_filtered['environment'] == env]
    n_g = env_data['genes_total'].values  # Total genes
    # Fit scaling laws for each category...
```

---

**End of Document**

