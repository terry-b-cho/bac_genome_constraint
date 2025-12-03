#!/bin/bash
#SBATCH --job-name=ncbi_organize
#SBATCH --output=logs/organize_%j.out
#SBATCH --error=logs/organize_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --partition=short

set -e

BASE_DIR="/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint"
DATA_DIR="${BASE_DIR}/data/ncbi"
EXTRACT_DIR="${DATA_DIR}/bac_genomes_5988_dir"
ASSEMBLIES_DIR="${DATA_DIR}/assemblies"
METADATA_DIR="${DATA_DIR}/metadata"
LOG_DIR="${BASE_DIR}/scripts/0_data_download/logs"

mkdir -p "${LOG_DIR}"

if [ ! -d "${EXTRACT_DIR}/ncbi_dataset/data" ]; then
    echo "Error: Rehydrated data not found. Run 03_rehydrate.sh first."
    exit 1
fi

echo "Organizing files at $(date)"

mkdir -p "${ASSEMBLIES_DIR}"
mkdir -p "${METADATA_DIR}"

# Move assembly directories (exclude JSON files)
find "${EXTRACT_DIR}/ncbi_dataset/data" -maxdepth 1 -type d ! -path "${EXTRACT_DIR}/ncbi_dataset/data" -exec mv {} "${ASSEMBLIES_DIR}/" \;

# Copy metadata files
if [ -f "${EXTRACT_DIR}/ncbi_dataset/data/assembly_data_report.jsonl" ]; then
    cp "${EXTRACT_DIR}/ncbi_dataset/data/assembly_data_report.jsonl" "${METADATA_DIR}/"
fi

if [ -f "${EXTRACT_DIR}/ncbi_dataset/data/dataset_catalog.json" ]; then
    cp "${EXTRACT_DIR}/ncbi_dataset/data/dataset_catalog.json" "${METADATA_DIR}/"
fi

ls -1 "${ASSEMBLIES_DIR}" > "${METADATA_DIR}/assembly_manifest.txt"

echo "Organization completed at $(date)"
echo "Total assemblies: $(wc -l < ${METADATA_DIR}/assembly_manifest.txt)"


