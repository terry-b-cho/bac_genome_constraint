#!/bin/bash
#SBATCH --job-name=ncbi_metadata
#SBATCH --output=logs/metadata_%j.out
#SBATCH --error=logs/metadata_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=10G
#SBATCH --partition=short

set -e

BASE_DIR="/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint"
DATA_DIR="${BASE_DIR}/data/ncbi"
ZIP_FILE="${DATA_DIR}/bac_genomes_5988.zip"
METADATA_DIR="${DATA_DIR}/metadata"
OUTPUT_TSV="${METADATA_DIR}/metadata_table.tsv"
LOG_DIR="${BASE_DIR}/scripts/0_data_download/logs"

mkdir -p "${LOG_DIR}"
mkdir -p "${METADATA_DIR}"

if [ ! -f "${ZIP_FILE}" ]; then
    echo "Error: ${ZIP_FILE} not found."
    exit 1
fi

module load conda/miniforge3/24.11.3-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ncbi_datasets

echo "Extracting metadata at $(date)"

dataformat tsv genome \
  --package "${ZIP_FILE}" \
  --fields accession,assminfo-release-date,assminfo-level,annotinfo-name,organism-name,organism-tax-id \
  > "${OUTPUT_TSV}"

echo "Metadata extraction completed at $(date)"
echo "Output: ${OUTPUT_TSV}"
echo "Lines: $(wc -l < ${OUTPUT_TSV})"

