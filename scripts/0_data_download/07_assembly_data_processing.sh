#!/bin/bash
#SBATCH --job-name=ncbi_process_jsonl
#SBATCH --output=logs/process_jsonl_%j.out
#SBATCH --error=logs/process_jsonl_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=20G
#SBATCH --partition=short

set -e

BASE_DIR="/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint"
DATA_DIR="${BASE_DIR}/data/ncbi"
METADATA_DIR="${DATA_DIR}/metadata"
INPUT_JSONL="${METADATA_DIR}/assembly_data_report.jsonl"
OUTPUT_TSV="${METADATA_DIR}/assembly_data_report_extracted.tsv"
INPUT_CATALOG="${METADATA_DIR}/dataset_catalog.json"
OUTPUT_CATALOG_TSV="${METADATA_DIR}/dataset_catalog_extracted.tsv"
LOG_DIR="${BASE_DIR}/scripts/0_data_download/logs"

mkdir -p "${LOG_DIR}"

if [ ! -f "${INPUT_JSONL}" ]; then
    echo "Error: ${INPUT_JSONL} not found."
    exit 1
fi

if [ ! -f "${INPUT_CATALOG}" ]; then
    echo "Error: ${INPUT_CATALOG} not found."
    exit 1
fi

echo "Processing assembly_data_report.jsonl at $(date)"
echo "Input: ${INPUT_JSONL}"
echo "Output: ${OUTPUT_TSV}"

export INPUT_JSONL="${INPUT_JSONL}"
export OUTPUT_TSV="${OUTPUT_TSV}"
export INPUT_CATALOG="${INPUT_CATALOG}"
export OUTPUT_CATALOG_TSV="${OUTPUT_CATALOG_TSV}"

python3 << 'PYTHON_SCRIPT'
import json
import csv
import os
import sys
from pathlib import Path

def safe_get(data, *keys, default=''):
    """Safely get nested dictionary values"""
    for key in keys:
        if isinstance(data, dict):
            data = data.get(key, {})
        else:
            return default
    return data if data != {} else default

def safe_get_nested(data, path, default=''):
    """Get value from nested path like 'assemblyInfo.assemblyLevel'"""
    keys = path.split('.')
    return safe_get(data, *keys, default=default)

def extract_list_value(data, path, index=0, subkey=None, default=''):
    """Extract value from list, optionally from a specific item"""
    keys = path.split('.')
    value = safe_get(data, *keys, default=[])
    if isinstance(value, list) and len(value) > index:
        item = value[index]
        if subkey and isinstance(item, dict):
            return item.get(subkey, default)
        return item if not subkey else default
    return default

def extract_attributes(data, path, name_key='name', value_key='value'):
    """Extract biosample attributes as key-value pairs"""
    keys = path.split('.')
    attrs = safe_get(data, *keys, default=[])
    if isinstance(attrs, list):
        result = {}
        for attr in attrs:
            if isinstance(attr, dict):
                name = attr.get(name_key, '')
                value = attr.get(value_key, '')
                if name:
                    result[name] = value
        return result
    return {}

input_file = Path(os.environ['INPUT_JSONL'])
output_file = Path(os.environ['OUTPUT_TSV'])

# Define all columns to extract
columns = [
    # Basic identifiers
    ('accession', lambda d: safe_get_nested(d, 'accession')),
    ('currentAccession', lambda d: safe_get_nested(d, 'currentAccession')),
    ('pairedAccession', lambda d: safe_get_nested(d, 'pairedAccession')),
    ('sourceDatabase', lambda d: safe_get_nested(d, 'sourceDatabase')),
    
    # Organism info
    ('organism_taxId', lambda d: safe_get_nested(d, 'organism.taxId')),
    ('organism_name', lambda d: safe_get_nested(d, 'organism.organismName')),
    ('organism_strain', lambda d: safe_get_nested(d, 'organism.infraspecificNames.strain', default='')),
    ('organism_subspecies', lambda d: safe_get_nested(d, 'organism.infraspecificNames.subspecies', default='')),
    
    # Assembly info
    ('assembly_level', lambda d: safe_get_nested(d, 'assemblyInfo.assemblyLevel')),
    ('assembly_name', lambda d: safe_get_nested(d, 'assemblyInfo.assemblyName')),
    ('assembly_type', lambda d: safe_get_nested(d, 'assemblyInfo.assemblyType')),
    ('assembly_status', lambda d: safe_get_nested(d, 'assemblyInfo.assemblyStatus')),
    ('assembly_releaseDate', lambda d: safe_get_nested(d, 'assemblyInfo.releaseDate')),
    ('assembly_submitter', lambda d: safe_get_nested(d, 'assemblyInfo.submitter')),
    ('assembly_refseqCategory', lambda d: safe_get_nested(d, 'assemblyInfo.refseqCategory')),
    ('assembly_bioprojectAccession', lambda d: safe_get_nested(d, 'assemblyInfo.bioprojectAccession')),
    
    # Paired assembly
    ('pairedAssembly_accession', lambda d: safe_get_nested(d, 'assemblyInfo.pairedAssembly.accession', default='')),
    ('pairedAssembly_status', lambda d: safe_get_nested(d, 'assemblyInfo.pairedAssembly.status', default='')),
    
    # Biosample
    ('biosample_accession', lambda d: safe_get_nested(d, 'assemblyInfo.biosample.accession', default='')),
    ('biosample_strain', lambda d: safe_get_nested(d, 'assemblyInfo.biosample.strain', default='')),
    ('biosample_subSpecies', lambda d: safe_get_nested(d, 'assemblyInfo.biosample.subSpecies', default='')),
    ('biosample_publicationDate', lambda d: safe_get_nested(d, 'assemblyInfo.biosample.publicationDate', default='')),
    ('biosample_submissionDate', lambda d: safe_get_nested(d, 'assemblyInfo.biosample.submissionDate', default='')),
    
    # Biosample attributes (common ones)
    ('biosample_typeMaterial', lambda d: extract_attributes(d, 'assemblyInfo.biosample.attributes').get('type-material', '')),
    ('biosample_isolation_source', lambda d: extract_attributes(d, 'assemblyInfo.biosample.attributes').get('isolation-source', '')),
    ('biosample_collection_date', lambda d: extract_attributes(d, 'assemblyInfo.biosample.attributes').get('collection-date', '')),
    ('biosample_geo_loc_name', lambda d: extract_attributes(d, 'assemblyInfo.biosample.attributes').get('geo_loc-name', '')),
    ('biosample_host', lambda d: extract_attributes(d, 'assemblyInfo.biosample.attributes').get('host', '')),
    
    # Assembly stats
    ('stats_totalSequenceLength', lambda d: safe_get_nested(d, 'assemblyStats.totalSequenceLength')),
    ('stats_totalUngappedLength', lambda d: safe_get_nested(d, 'assemblyStats.totalUngappedLength')),
    ('stats_numberOfContigs', lambda d: safe_get_nested(d, 'assemblyStats.numberOfContigs')),
    ('stats_numberOfScaffolds', lambda d: safe_get_nested(d, 'assemblyStats.numberOfScaffolds')),
    ('stats_numberOfChromosomes', lambda d: safe_get_nested(d, 'assemblyStats.totalNumberOfChromosomes')),
    ('stats_contigN50', lambda d: safe_get_nested(d, 'assemblyStats.contigN50')),
    ('stats_contigL50', lambda d: safe_get_nested(d, 'assemblyStats.contigL50')),
    ('stats_scaffoldN50', lambda d: safe_get_nested(d, 'assemblyStats.scaffoldN50')),
    ('stats_scaffoldL50', lambda d: safe_get_nested(d, 'assemblyStats.scaffoldL50')),
    ('stats_gcPercent', lambda d: safe_get_nested(d, 'assemblyStats.gcPercent')),
    ('stats_gcCount', lambda d: safe_get_nested(d, 'assemblyStats.gcCount')),
    ('stats_atgcCount', lambda d: safe_get_nested(d, 'assemblyStats.atgcCount')),
    ('stats_genomeCoverage', lambda d: safe_get_nested(d, 'assemblyStats.genomeCoverage')),
    
    # Annotation info
    ('annotation_name', lambda d: safe_get_nested(d, 'annotationInfo.name', default='')),
    ('annotation_provider', lambda d: safe_get_nested(d, 'annotationInfo.provider', default='')),
    ('annotation_releaseDate', lambda d: safe_get_nested(d, 'annotationInfo.releaseDate', default='')),
    ('annotation_method', lambda d: safe_get_nested(d, 'annotationInfo.method', default='')),
    ('annotation_pipeline', lambda d: safe_get_nested(d, 'annotationInfo.pipeline', default='')),
    ('annotation_softwareVersion', lambda d: safe_get_nested(d, 'annotationInfo.softwareVersion', default='')),
    
    # Gene counts
    ('genes_total', lambda d: safe_get_nested(d, 'annotationInfo.stats.geneCounts.total', default='')),
    ('genes_proteinCoding', lambda d: safe_get_nested(d, 'annotationInfo.stats.geneCounts.proteinCoding', default='')),
    ('genes_nonCoding', lambda d: safe_get_nested(d, 'annotationInfo.stats.geneCounts.nonCoding', default='')),
    ('genes_pseudogene', lambda d: safe_get_nested(d, 'annotationInfo.stats.geneCounts.pseudogene', default='')),
    
    # CheckM info
    ('checkm_completeness', lambda d: safe_get_nested(d, 'checkmInfo.completeness', default='')),
    ('checkm_contamination', lambda d: safe_get_nested(d, 'checkmInfo.contamination', default='')),
    ('checkm_completenessPercentile', lambda d: safe_get_nested(d, 'checkmInfo.completenessPercentile', default='')),
    ('checkm_markerSet', lambda d: safe_get_nested(d, 'checkmInfo.checkmMarkerSet', default='')),
    ('checkm_markerSetRank', lambda d: safe_get_nested(d, 'checkmInfo.checkmMarkerSetRank', default='')),
    ('checkm_version', lambda d: safe_get_nested(d, 'checkmInfo.checkmVersion', default='')),
    ('checkm_speciesTaxId', lambda d: safe_get_nested(d, 'checkmInfo.checkmSpeciesTaxId', default='')),
    
    # ANI info
    ('ani_taxonomyCheckStatus', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.taxonomyCheckStatus', default='')),
    ('ani_matchStatus', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.matchStatus', default='')),
    ('ani_submittedOrganism', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.submittedOrganism', default='')),
    ('ani_submittedSpecies', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.submittedSpecies', default='')),
    ('ani_category', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.category', default='')),
    ('ani_bestMatch_assembly', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.bestAniMatch.assembly', default='')),
    ('ani_bestMatch_organism', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.bestAniMatch.organismName', default='')),
    ('ani_bestMatch_ani', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.bestAniMatch.ani', default='')),
    ('ani_bestMatch_coverage', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.bestAniMatch.assemblyCoverage', default='')),
    ('ani_submittedMatch_assembly', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.submittedAniMatch.assembly', default='')),
    ('ani_submittedMatch_organism', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.submittedAniMatch.organismName', default='')),
    ('ani_submittedMatch_ani', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.submittedAniMatch.ani', default='')),
    ('ani_submittedMatch_coverage', lambda d: safe_get_nested(d, 'averageNucleotideIdentity.submittedAniMatch.assemblyCoverage', default='')),
    
    # Type material
    ('typeMaterial_label', lambda d: safe_get_nested(d, 'typeMaterial.typeLabel', default='')),
    ('typeMaterial_displayText', lambda d: safe_get_nested(d, 'typeMaterial.typeDisplayText', default='')),
]

# Process JSONL file
rows = []
with open(input_file, 'r') as f:
    for line_num, line in enumerate(f, 1):
        try:
            data = json.loads(line.strip())
            row = {}
            for col_name, extract_func in columns:
                try:
                    value = extract_func(data)
                    # Convert to string, handle None
                    row[col_name] = str(value) if value not in [None, {}, []] else ''
                except Exception as e:
                    row[col_name] = ''
            rows.append(row)
        except json.JSONDecodeError as e:
            print(f"Warning: Skipping line {line_num} due to JSON error: {e}", file=sys.stderr)
            continue

# Write to TSV
if rows:
    with open(output_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[col[0] for col in columns], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"Successfully processed {len(rows)} entries")
    print(f"Output written to: {output_file}")
    print(f"Total columns: {len(columns)}")
else:
    print("Error: No rows processed")
    exit(1)
PYTHON_SCRIPT

echo "Processing completed at $(date)"
echo "Output: ${OUTPUT_TSV}"
echo "Lines: $(wc -l < ${OUTPUT_TSV})"

echo ""
echo "Processing dataset_catalog.json at $(date)"
echo "Input: ${INPUT_CATALOG}"
echo "Output: ${OUTPUT_CATALOG_TSV}"

python3 << 'PYTHON_SCRIPT2'
import json
import csv
import os
from pathlib import Path

input_catalog = Path(os.environ['INPUT_CATALOG'])
output_catalog = Path(os.environ['OUTPUT_CATALOG_TSV'])

# Process catalog JSON
with open(input_catalog, 'r') as f:
    catalog_data = json.load(f)

rows = []
for assembly in catalog_data.get('assemblies', []):
    accession = assembly.get('accession', '')
    files = assembly.get('files', [])
    
    for file_info in files:
        row = {
            'accession': accession,
            'filePath': file_info.get('filePath', ''),
            'fileType': file_info.get('fileType', ''),
            'uncompressedLengthBytes': str(file_info.get('uncompressedLengthBytes', ''))
        }
        rows.append(row)

# Write to TSV
if rows:
    with open(output_catalog, 'w', newline='') as f:
        fieldnames = ['accession', 'filePath', 'fileType', 'uncompressedLengthBytes']
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    
    print(f"Successfully processed {len(rows)} file entries")
    print(f"Output written to: {output_catalog}")
    print(f"Total columns: {len(fieldnames)}")
else:
    print("Error: No rows processed")
    exit(1)
PYTHON_SCRIPT2

echo "Catalog processing completed at $(date)"
echo "Output: ${OUTPUT_CATALOG_TSV}"
echo "Lines: $(wc -l < ${OUTPUT_CATALOG_TSV})"
