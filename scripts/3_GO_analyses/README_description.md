
# GO Analysis Results Documentation

## Overview
This directory contains results from Gene Ontology (GO) annotation analyses performed on bacterial genomes. The analyses extract GO terms from GFF files, identify unique terms across all genomes, and create presence/absence and count tables.

## Scripts

### `scripts/3_GO_analyses/go_analyses.ipynb`
A Jupyter notebook that performs the following analyses:

1. **Extract GO annotations from GFF files** (Cell 1-4):
   - Reads GFF files for each genome accession
   - Extracts GO terms from the `Ontology_term` field in the attributes column
   - Saves unique GO terms for each genome to individual text files in `results/3_GO_analyses/GO_lists/`
   - Also creates files with term counts in `results/3_GO_analyses/GO_lists_with_counts/`

2. **Identify unique GO terms** (Cell 5-6):
   - Collects all GO terms from all genome files
   - Creates a sorted list of unique GO term IDs
   - Identifies genomes with no GO terms (saved to `missing_terms.txt`)

3. **Create presence/absence table** (Cell 7):
   - Generates a matrix where rows are genomes and columns are unique GO terms
   - Values are 1 (term present) or 0 (term absent) for each genome-term combination

4. **Identify ubiquitous terms** (Cell 8):
   - Finds GO terms present in ≥95% of genomes
   - Saves list to `ubiquitous_terms.txt`

5. **Generate count tables** (Cell 9-12):
   - Creates tables with actual counts (not just presence/absence) for ubiquitous terms
   - Normalizes counts by number of protein-coding genes per genome

### `scripts/3_GO_analyses/load_go_ontology.py`
A Python module with two functions:

1. **`load_full_ontology(term_list=None)`**:
   - Loads the full GO ontology from `go-basic.obo` file
   - If `term_list` is provided, only loads data for those terms
   - Returns a dictionary where keys are GO term IDs and values are dictionaries containing all ontology data for that term
   - **Note**: There is a syntax error on line 31: `len(term_list == 0)` should be `len(term_list) == 0`

2. **`load_ids_and_names(term_list=None)`**:
   - Loads GO term IDs and their names from the ontology file
   - Returns a dictionary mapping GO term IDs to their names
   - More efficient than `load_full_ontology` when only IDs and names are needed

## Result Files

### `unique_terms.txt`
- **Format**: One GO term ID per line (7-digit zero-padded format, e.g., `0000015`)
- **Content**: Sorted list of all unique GO term IDs found across all analyzed genomes
- **Usage**: Reference list of all GO terms in the dataset. Can be used to:
  - Map term IDs to names using the GO ontology
  - Filter analyses to specific terms
  - Understand the scope of functional annotations

### `ubiquitous_terms.txt`
- **Format**: One GO term ID per line (7-digit zero-padded format)
- **Content**: List of GO terms present in ≥95% of genomes (334 terms)
- **Usage**: Identifies core functional terms that are nearly universal across the bacterial genomes. Useful for:
  - Identifying essential/core functions
  - Filtering out rare terms for certain analyses
  - Understanding conserved functional elements

### `missing_terms.txt`
- **Format**: One filename per line (e.g., `GCF_000006765.1.txt`)
- **Content**: List of genome files that contained no GO terms (13 genomes)
- **Usage**: Identifies genomes that may have annotation issues or lack functional annotations. These genomes were excluded from the presence/absence table.

### `missing_genomes.txt`
- **Format**: Empty file (or one accession per line if any genomes were missing)
- **Content**: Genomes that were expected but not found in the data directory
- **Usage**: Tracks which genomes from the modeling accessions list were not available for analysis

### `presence_absence_table.txt`
- **Format**: Tab-separated values (TSV)
- **Structure**:
  - First row: Header with "Genome" followed by all unique GO term IDs (7-digit format)
  - Subsequent rows: One row per genome
    - First column: Genome accession (e.g., `GCF_000006985.1`)
    - Remaining columns: 1 if term is present, 0 if absent
- **Dimensions**: ~3,088 genomes × ~3,834 unique GO terms
- **Usage**: 
  - Binary matrix for presence/absence analysis
  - Can be loaded into R/Python for statistical analysis
  - Useful for clustering genomes by functional similarity
  - Can identify co-occurring or mutually exclusive GO terms

### `ubiquitous_counts_table.txt`
- **Format**: Tab-separated values (TSV)
- **Structure**:
  - First row: Header with "Genome" followed by ubiquitous GO term IDs
  - Subsequent rows: One row per genome
    - First column: Genome accession
    - Remaining columns: Integer counts of how many times each ubiquitous term appears in that genome
- **Dimensions**: ~3,088 genomes × 334 ubiquitous terms
- **Usage**:
  - Quantitative analysis of ubiquitous term abundance
  - Can identify genomes with unusually high/low counts of specific terms
  - Useful for understanding gene copy number variation in core functions

### `ubiquitous_counts_table_normalized.txt`
- **Format**: Tab-separated values (TSV)
- **Structure**: Same as `ubiquitous_counts_table.txt` but values are normalized
- **Normalization**: Each count is divided by the number of protein-coding genes in that genome
- **Values**: Floating-point numbers (e.g., `0.001006036217303823`)
- **Usage**:
  - Allows comparison across genomes of different sizes
  - Normalized values represent the proportion of genes annotated with each term
  - Useful for identifying genomes with unusual functional composition relative to genome size
  - Can be used for comparative genomics analyses that account for genome size differences

## How to Use These Files

### In Python:
```python
import pandas as pd

# Load presence/absence table
pa_table = pd.read_csv('presence_absence_table.txt', sep='\t', index_col=0)

# Load normalized counts
norm_counts = pd.read_csv('ubiquitous_counts_table_normalized.txt', sep='\t', index_col=0)

# Load unique terms list
with open('unique_terms.txt', 'r') as f:
    unique_terms = [line.strip() for line in f]

# Load ubiquitous terms
with open('ubiquitous_terms.txt', 'r') as f:
    ubiquitous_terms = [line.strip() for line in f]
```

### In R:
```r
# Load presence/absence table
pa_table <- read.table('presence_absence_table.txt', sep='\t', header=TRUE, row.names=1)

# Load normalized counts
norm_counts <- read.table('ubiquitous_counts_table_normalized.txt', sep='\t', header=TRUE, row.names=1)

# Load unique terms
unique_terms <- readLines('unique_terms.txt')

# Load ubiquitous terms
ubiquitous_terms <- readLines('ubiquitous_terms.txt')
```

### Mapping GO IDs to Names:
Use the `load_go_ontology.py` script or query the GO database to map 7-digit IDs (e.g., `0000015`) to full GO IDs (e.g., `GO:0000015`) and names.

## Notes
- GO term IDs are stored in 7-digit zero-padded format (e.g., `0000015` instead of `GO:0000015`)
- To convert to standard GO format, prepend `GO:` to each ID
- The analysis includes ~3,088 genomes with GO annotations
- 13 genomes had no GO terms and were excluded from tables
- 334 terms are considered "ubiquitous" (present in ≥95% of genomes)

