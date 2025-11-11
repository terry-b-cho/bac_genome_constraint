# NCBI Genome Data Download

This directory contains scripts to download and organize ~5,988 bacterial reference genomes from NCBI RefSeq. The genomes are filtered for complete assemblies, annotated sequences, and reference quality standards.

## Prerequisites

1. Install NCBI Datasets CLI:
```bash
conda create -n ncbi_datasets -c conda-forge ncbi-datasets-cli
conda activate ncbi_datasets
```

2. Ensure write permissions in:
   - `/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint/data/ncbi`
   - `/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint/data/metadata`

## Workflow

Run scripts in order:

1. **01_extract_accessions.sh** - Extract assembly accessions from TSV
   ```bash
   bash 01_extract_accessions.sh
   ```

2. **02_download_dehydrated.sh** - Download dehydrated genome package
   ```bash
   sbatch 02_download_dehydrated.sh
   ```

3. **03_rehydrate.sh** - Rehydrate (download actual sequence files)
   ```bash
   sbatch 03_rehydrate.sh
   ```

4. **04_organize.sh** - Organize files into standard structure
   ```bash
   sbatch 04_organize.sh
   ```

5. **05_qc.sh** - Quality control checks
   ```bash
   sbatch 05_qc.sh
   ```

6. **06_extract_metadata.sh** - Extract basic metadata table
   ```bash
   sbatch 06_extract_metadata.sh
   ```

7. **07_assembly_data_processing.sh** - Extract detailed metadata from JSON files
   ```bash
   sbatch 07_assembly_data_processing.sh
   ```

## Directory Structure

Complete directory tree with inline comments:

```
data/ncbi/
├── assemblies/                    # [Created by: 04_organize.sh] One directory per assembly accession (5,986 total)
│   ├── GCF_000005845.2/          # Assembly accession directory
│   │   ├── GCF_000005845.2_ASM584v2_genomic.fna  # [03_rehydrate.sh] Genomic DNA sequences (FASTA)
│   │   ├── protein.faa                            # [03_rehydrate.sh] Protein sequences (FASTA, amino acids)
│   │   ├── cds_from_genomic.fna                   # [03_rehydrate.sh] Coding DNA sequences (FASTA)
│   │   ├── genomic.gff                            # [03_rehydrate.sh] Gene annotations (GFF3 format)
│   │   ├── genomic.gbff                           # [03_rehydrate.sh] Complete annotation (GenBank format)
│   │   └── sequence_report.jsonl                  # [03_rehydrate.sh] Sequence statistics (optional, may be missing)
│   ├── GCF_000006685.1/
│   └── ... (5,986 total directories)
│
├── metadata/                      # [Created by: Multiple scripts] Summary files, manifests, and extracted metadata
│   ├── assembly_manifest.txt      # [04_organize.sh] List of all assembly accessions (one per line, 5,986 rows)
│   ├── qc_manifest.txt            # [05_qc.sh] Quality control results (OK/MISSING status for each assembly)
│   ├── metadata_table.tsv         # [06_extract_metadata.sh] Basic metadata summary (6 columns: accession, date, level, annotation, organism, taxID)
│   ├── assembly_data_report.jsonl # [02_download_dehydrated.sh, 04_organize.sh] Complete assembly metadata (JSON Lines, 23 MB, 5,986 rows)
│   ├── assembly_data_report_extracted.tsv  # [07_assembly_data_processing.sh] Extracted detailed metadata (TSV, 73 columns, 5,987 rows)
│   ├── dataset_catalog.json       # [02_download_dehydrated.sh, 04_organize.sh] File catalog (JSON, 5.3 MB)
│   └── dataset_catalog_extracted.tsv  # [07_assembly_data_processing.sh] Extracted file catalog (TSV, 4 columns, 35,918 rows)
│
├── bac_genomes_5988_dir/          # [Created by: 03_rehydrate.sh] Original extracted structure from NCBI (kept for reference)
│   ├── README.md                  # NCBI-provided documentation
│   ├── md5sum.txt                 # MD5 checksums for file integrity verification
│   └── ncbi_dataset/
│       ├── fetch.txt              # URLs for downloading files (35,917 lines)
│       └── data/
│           ├── assembly_data_report.jsonl  # Original metadata (copy)
│           └── dataset_catalog.json        # Original catalog (copy)
│
└── bac_genomes_5988.zip           # [Created by: 02_download_dehydrated.sh] Dehydrated package (5.8 MB, contains metadata and file pointers)
```

## File Descriptions

### Genome Files (`assemblies/`)

Each assembly directory contains:

- **`*_genomic.fna`**: Genomic DNA sequences in FASTA format. Naming pattern: `{ACCESSION}_{ASSEMBLY_NAME}_genomic.fna` (e.g., `GCF_000006945.2_ASM694v2_genomic.fna`). Use glob pattern `*_genomic.fna` to find files.

- **`protein.faa`**: Protein sequences in FASTA format (amino acids). Used for KEGG annotation and protein analysis.

- **`cds_from_genomic.fna`**: Coding DNA sequences (CDS) extracted from genomic sequences. Used for gene analysis and codon usage studies.

- **`genomic.gff`**: Gene annotation in GFF3 format. Contains features: genes, CDS, etc. Columns: seqid, source, type, start, end, score, strand, phase, attributes.

- **`genomic.gbff`**: GenBank flat file format with complete annotation. Contains full annotation details and submission metadata.

- **`sequence_report.jsonl`** (optional): Sequence statistics and metadata in JSON Lines format. Some assemblies may be missing this file due to download failures.

### Metadata Files (`metadata/`)

#### `assembly_manifest.txt`
- **Format**: One accession per line
- **Rows**: 5,986
- **Use**: Iterating through all assemblies in processing pipelines
- **Example**:
  ```
  GCF_000005845.2
  GCF_000006685.1
  GCF_000006945.2
  ```

#### `qc_manifest.txt`
- **Format**: Text file with QC results
- **Rows**: 5,986
- **Format**: `OK {accession}` or `MISSING {accession}` followed by missing file details
- **Use**: Identifying assemblies with missing required files
- **Status**: All 5,986 assemblies passed QC (OK=5986, MISSING=0)

#### `metadata_table.tsv`
- **Format**: Tab-separated values
- **Rows**: 5,987 (header + 5,986 genomes)
- **Columns**:
  1. `Assembly Accession` - RefSeq accession (GCF_*)
  2. `Assembly Release Date` - YYYY-MM-DD format
  3. `Assembly Level` - Complete, Chromosome, Scaffold, etc.
  4. `Annotation Name` - Annotation version (e.g., GCF_000006945.2-RS_2025_03_27)
  5. `Organism Name` - Scientific name
  6. `Organism Taxonomic ID` - NCBI taxonomy ID
- **Use**: Linking genomes to taxonomy, filtering by date/level, joining with other data

#### `assembly_data_report.jsonl`
- **Format**: JSON Lines (one JSON object per line)
- **Rows**: 5,986
- **Size**: 23 MB
- **Content**: Complete assembly metadata including taxonomy, statistics, annotation info, biosample data, CheckM quality metrics, ANI matches, etc.
- **Use**: Detailed metadata queries, programmatic access to full metadata

#### `assembly_data_report_extracted.tsv`
- **Format**: Tab-separated values
- **Rows**: 5,987 (header + 5,986 genomes)
- **Columns**: 73
- **Size**: 4.8 MB
- **Source**: Extracted from `assembly_data_report.jsonl` using `07_assembly_data_processing.sh`
- **Column Categories**:
  - **Identifiers** (4): accession, currentAccession, pairedAccession, sourceDatabase
  - **Organism Info** (4): organism_taxId, organism_name, organism_strain, organism_subspecies
  - **Assembly Info** (11): assembly_level, assembly_name, assembly_type, assembly_status, assembly_releaseDate, assembly_submitter, assembly_refseqCategory, assembly_bioprojectAccession, pairedAssembly_accession, pairedAssembly_status
  - **Biosample Info** (9): biosample_accession, biosample_strain, biosample_subSpecies, biosample_publicationDate, biosample_submissionDate, biosample_typeMaterial, biosample_isolation_source, biosample_collection_date, biosample_geo_loc_name, biosample_host
  - **Assembly Statistics** (13): stats_totalSequenceLength, stats_totalUngappedLength, stats_numberOfContigs, stats_numberOfScaffolds, stats_numberOfChromosomes, stats_contigN50, stats_contigL50, stats_scaffoldN50, stats_scaffoldL50, stats_gcPercent, stats_gcCount, stats_atgcCount, stats_genomeCoverage
  - **Annotation Info** (7): annotation_name, annotation_provider, annotation_releaseDate, annotation_method, annotation_pipeline, annotation_softwareVersion
  - **Gene Counts** (4): genes_total, genes_proteinCoding, genes_nonCoding, genes_pseudogene
  - **CheckM Quality** (7): checkm_completeness, checkm_contamination, checkm_completenessPercentile, checkm_markerSet, checkm_markerSetRank, checkm_version, checkm_speciesTaxId
  - **ANI** (12): ani_taxonomyCheckStatus, ani_matchStatus, ani_submittedOrganism, ani_submittedSpecies, ani_category, ani_bestMatch_assembly, ani_bestMatch_organism, ani_bestMatch_ani, ani_bestMatch_coverage, ani_submittedMatch_assembly, ani_submittedMatch_organism, ani_submittedMatch_ani, ani_submittedMatch_coverage
  - **Type Material** (2): typeMaterial_label, typeMaterial_displayText
- **Use**: Comprehensive metadata analysis, filtering, joining with other datasets
- **Notes**: Empty fields are empty strings. Some fields may be missing for certain assemblies. Numeric fields stored as strings (convert as needed). Dates in ISO format.

#### `dataset_catalog.json`
- **Format**: JSON
- **Size**: 5.3 MB
- **Content**: Catalog of all files in dataset with types and locations
- **Use**: Understanding file structure, verifying completeness

#### `dataset_catalog_extracted.tsv`
- **Format**: Tab-separated values
- **Rows**: 35,918 (header + 35,917 file entries)
- **Columns**: 4
- **Source**: Extracted from `dataset_catalog.json` using `07_assembly_data_processing.sh`
- **Columns**:
  1. `accession` - Assembly accession (GCF_*) or empty for report files
  2. `filePath` - Relative path to file
  3. `fileType` - File type (GENOMIC_NUCLEOTIDE_FASTA, PROTEIN_FASTA, GFF3, etc.)
  4. `uncompressedLengthBytes` - File size in bytes
- **Use**: File inventory, verifying file presence, calculating total sizes

## Usage Patterns

### Pattern 1: Iterate Through All Genomes

```python
# Read manifest
with open('data/ncbi/metadata/assembly_manifest.txt') as f:
    accessions = [line.strip() for line in f]

# Process each genome
for acc in accessions:
    genome_dir = f'data/ncbi/assemblies/{acc}/'
    genome_fna = glob.glob(f'{genome_dir}/*_genomic.fna')[0]
    protein_faa = f'{genome_dir}/protein.faa'
    # Process files...
```

### Pattern 2: Map Accession to Metadata

```python
import pandas as pd

# Load metadata table
metadata = pd.read_csv('data/ncbi/metadata/metadata_table.tsv', sep='\t')

# Map accession to organism
acc_to_organism = dict(zip(metadata['Assembly Accession'], 
                           metadata['Organism Name']))

# Get genome for specific organism
target_org = 'Salmonella enterica'
target_acc = metadata[metadata['Organism Name'].str.contains(target_org)]['Assembly Accession'].values[0]
genome_file = f'data/ncbi/assemblies/{target_acc}/*_genomic.fna'
```

### Pattern 3: Filter by Assembly Level

```python
# Filter for complete genomes only
complete_genomes = metadata[metadata['Assembly Level'] == 'Complete Genome']
complete_accessions = complete_genomes['Assembly Accession'].tolist()

# Process only complete genomes
for acc in complete_accessions:
    genome_dir = f'data/ncbi/assemblies/{acc}/'
    # Process...
```

### Pattern 4: Access Specific File Types

```python
import glob

acc = 'GCF_000006945.2'
genome_dir = f'data/ncbi/assemblies/{acc}/'

# Find files (handles variable naming)
genome_fna = glob.glob(f'{genome_dir}/*_genomic.fna')[0]
protein_faa = f'{genome_dir}/protein.faa'
cds_fna = f'{genome_dir}/cds_from_genomic.fna'
gff = glob.glob(f'{genome_dir}/genomic.gff')[0]
gbff = f'{genome_dir}/genomic.gbff'
```

### Pattern 5: Use Detailed Metadata

```python
# Load detailed extracted metadata
detailed = pd.read_csv('data/ncbi/metadata/assembly_data_report_extracted.tsv', sep='\t')

# Filter by genome size
large_genomes = detailed[detailed['stats_totalSequenceLength'].astype(float) > 5000000]

# Filter by CheckM completeness
high_quality = detailed[detailed['checkm_completeness'].astype(float) > 95.0]

# Join with basic metadata
basic = pd.read_csv('data/ncbi/metadata/metadata_table.tsv', sep='\t')
merged = detailed.merge(basic, left_on='accession', right_on='Assembly Accession', how='inner')
```

### Pattern 6: Join with External Data

```python
# Load metadata
metadata = pd.read_csv('data/ncbi/metadata/metadata_table.tsv', sep='\t')

# Join with external data (e.g., environment data)
external_data = pd.read_csv('other_data.csv')
merged = metadata.merge(external_data, 
                       left_on='Assembly Accession',
                       right_on='genome_id',
                       how='inner')

# Access genomes for specific environment
env_genomes = merged[merged['environment'] == 'human_skin']
for acc in env_genomes['Assembly Accession']:
    genome_file = f'data/ncbi/assemblies/{acc}/*_genomic.fna'
```

## Key Mappings

| Accession | Directory | Genome File | Protein File | GFF File |
|-----------|-----------|-------------|--------------|----------|
| GCF_000006945.2 | `assemblies/GCF_000006945.2/` | `*_genomic.fna` | `protein.faa` | `genomic.gff` |
| GCF_000006685.1 | `assemblies/GCF_000006685.1/` | `*_genomic.fna` | `protein.faa` | `genomic.gff` |

**Note**: Genome FASTA files have variable naming: `{ACCESSION}_{ASSEMBLY_NAME}_genomic.fna`. Use glob pattern `*_genomic.fna` to find them.

## Script Workflow Summary

1. **01_extract_accessions.sh** → Creates `data/metadata/selected_genomes.txt` (list of 5,988 accessions)

2. **02_download_dehydrated.sh** → Downloads `bac_genomes_5988.zip` (5.8 MB, contains metadata and file pointers)

3. **03_rehydrate.sh** → 
   - Unzips `bac_genomes_5988.zip` → Creates `bac_genomes_5988_dir/`
   - Downloads all sequence files → Populates `bac_genomes_5988_dir/ncbi_dataset/data/GCF_*/`

4. **04_organize.sh** → 
   - Moves assembly directories → `assemblies/`
   - Copies metadata files → `metadata/`
   - Creates `assembly_manifest.txt`

5. **05_qc.sh** → Creates `metadata/qc_manifest.txt` (verifies all required files present)

6. **06_extract_metadata.sh** → Creates `metadata/metadata_table.tsv` (basic 6-column summary)

7. **07_assembly_data_processing.sh** → Creates:
   - `metadata/assembly_data_report_extracted.tsv` (detailed 73-column metadata)
   - `metadata/dataset_catalog_extracted.tsv` (file catalog)

## File Size Summary

- `assemblies/`: **126 GB** (5,986 genomes with sequences and annotations)
- `metadata/`: **29 MB** total
  - `assembly_data_report.jsonl`: 23 MB
  - `dataset_catalog.json`: 5.3 MB
  - `assembly_data_report_extracted.tsv`: 4.8 MB
  - Other files: <1 MB each
- `bac_genomes_5988_dir/`: **38 MB** (original extracted structure, mostly metadata)
- `bac_genomes_5988.zip`: **5.8 MB** (dehydrated package only)

## Data Access Summary

- **Total assemblies**: 5,986
- **Assembly directory**: `data/ncbi/assemblies/GCF_XXXXXXXX.X/`
- **Manifest**: `data/ncbi/metadata/assembly_manifest.txt`
- **Basic metadata**: `data/ncbi/metadata/metadata_table.tsv` (6 columns)
- **Detailed metadata**: `data/ncbi/metadata/assembly_data_report_extracted.tsv` (73 columns)
- **QC status**: All passed (see `data/ncbi/metadata/qc_manifest.txt`)

## Notes

- All scripts use Slurm job scheduling (except step 1)
- Logs are written to `logs/` directory
- Rehydration can be resumed if interrupted (simply re-run `03_rehydrate.sh`)
- QC script reports missing files for debugging
- Metadata extraction scripts process JSON/JSONL files into TSV format for easier analysis
- Genome FASTA files use variable naming; always use glob patterns (`*_genomic.fna`) to find them
