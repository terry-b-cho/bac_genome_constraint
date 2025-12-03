#!/bin/bash
#SBATCH --job-name=ncbi_rehydrate
#SBATCH --output=logs/rehydrate_%j.out
#SBATCH --error=logs/rehydrate_%j.err
#SBATCH --time=5-00:00:00
#SBATCH --mem=100G
#SBATCH --partition=medium

set -e

BASE_DIR="/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint"
DATA_DIR="${BASE_DIR}/data/ncbi"
ZIP_FILE="${DATA_DIR}/bac_genomes_5988.zip"
EXTRACT_DIR="${DATA_DIR}/bac_genomes_5988_dir"
LOG_DIR="${BASE_DIR}/scripts/0_data_download/logs"

mkdir -p "${LOG_DIR}"

if [ ! -f "${ZIP_FILE}" ]; then
    echo "Error: ${ZIP_FILE} not found. Run 02_download_dehydrated.sh first."
    exit 1
fi

module load conda/miniforge3/24.11.3-0
source $(conda info --base)/etc/profile.d/conda.sh
conda activate ncbi_datasets

cd "${DATA_DIR}"

if [ ! -d "${EXTRACT_DIR}" ]; then
    echo "Extracting zip file at $(date)"
    unzip -q "${ZIP_FILE}" -d "${EXTRACT_DIR}"
    echo "Extraction completed at $(date)"
fi

cd "${EXTRACT_DIR}"

echo "Starting rehydration at $(date)"
datasets rehydrate --directory .
echo "Rehydration completed at $(date)"

