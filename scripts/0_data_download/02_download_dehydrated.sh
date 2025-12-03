#!/bin/bash
#SBATCH --job-name=ncbi_download
#SBATCH --output=logs/download_%j.out
#SBATCH --error=logs/download_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=50G
#SBATCH --partition=short

set -e

BASE_DIR="/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint"
DATA_DIR="${BASE_DIR}/data/ncbi"
METADATA_DIR="${BASE_DIR}/data/metadata"
ACCESSIONS_FILE="${METADATA_DIR}/selected_genomes.txt"
ZIP_FILE="${DATA_DIR}/bac_genomes_5988.zip"
LOG_DIR="${BASE_DIR}/scripts/0_data_download/logs"

mkdir -p "${DATA_DIR}"
mkdir -p "${LOG_DIR}"

if [ ! -f "${ACCESSIONS_FILE}" ]; then
    echo "Error: ${ACCESSIONS_FILE} not found. Run 01_extract_accessions.sh first."
    exit 1
fi

module load conda/miniforge3/24.11.3-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ncbi_datasets

cd "${DATA_DIR}"

echo "Starting dehydrated download at $(date)"
echo "Input file: ${ACCESSIONS_FILE}"
echo "Output file: ${ZIP_FILE}"

datasets download genome accession \
  --inputfile "${ACCESSIONS_FILE}" \
  --assembly-source refseq \
  --annotated \
  --assembly-level complete \
  --exclude-atypical \
  --mag exclude \
  --include genome,gff3,gbff,protein,cds,seq-report \
  --filename "${ZIP_FILE}" \
  --dehydrated

echo "Download completed at $(date)"
echo "File size: $(du -h ${ZIP_FILE} | cut -f1)"

