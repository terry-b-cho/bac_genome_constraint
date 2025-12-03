#!/bin/bash
#SBATCH --job-name=ncbi_qc
#SBATCH --output=logs/qc_%j.out
#SBATCH --error=logs/qc_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=20G
#SBATCH --partition=short

set -e

BASE_DIR="/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint"
DATA_DIR="${BASE_DIR}/data/ncbi"
ASSEMBLIES_DIR="${DATA_DIR}/assemblies"
METADATA_DIR="${DATA_DIR}/metadata"
QC_OUTPUT="${METADATA_DIR}/qc_manifest.txt"
LOG_DIR="${BASE_DIR}/scripts/0_data_download/logs"

mkdir -p "${LOG_DIR}"
mkdir -p "${METADATA_DIR}"

if [ ! -d "${ASSEMBLIES_DIR}" ]; then
    echo "Error: Assemblies directory not found. Run 04_organize.sh first."
    exit 1
fi

echo "Running QC checks at $(date)"

ASSEMBLY_COUNT=$(ls -1 "${ASSEMBLIES_DIR}" | wc -l)
echo "Total assembly directories: ${ASSEMBLY_COUNT}"

echo "Checking required files for each assembly..."
> "${QC_OUTPUT}"

find "${ASSEMBLIES_DIR}" -maxdepth 1 -type d ! -path "${ASSEMBLIES_DIR}" | while read d; do
    acc=$(basename "$d")
    
    # Check for files with actual naming pattern (may have assembly name in filename)
    genomic_fna=$(find "$d" -name "*_genomic.fna" -type f | head -1)
    protein_faa=$(find "$d" -name "protein.faa" -type f | head -1)
    cds_fna=$(find "$d" -name "cds_from_genomic.fna" -type f | head -1)
    gff=$(find "$d" -name "genomic.gff" -o -name "*.gff3" -type f | head -1)
    gbff=$(find "$d" -name "genomic.gbff" -type f | head -1)
    
    if [ -n "$genomic_fna" ] && [ -n "$protein_faa" ] && [ -n "$cds_fna" ] && [ -n "$gff" ]; then
        echo "OK $acc" >> "${QC_OUTPUT}"
    else
        echo "MISSING $acc" >> "${QC_OUTPUT}"
        [ -z "$genomic_fna" ] && echo "  Missing: *_genomic.fna" >> "${QC_OUTPUT}"
        [ -z "$protein_faa" ] && echo "  Missing: protein.faa" >> "${QC_OUTPUT}"
        [ -z "$cds_fna" ] && echo "  Missing: cds_from_genomic.fna" >> "${QC_OUTPUT}"
        [ -z "$gff" ] && echo "  Missing: genomic.gff/gff3" >> "${QC_OUTPUT}"
    fi
done

OK_COUNT=$(grep -c "^OK" "${QC_OUTPUT}" || echo "0")
MISSING_COUNT=$(grep -c "^MISSING" "${QC_OUTPUT}" || echo "0")

echo "QC completed at $(date)"
echo "OK: ${OK_COUNT}"
echo "MISSING: ${MISSING_COUNT}"
echo "QC report: ${QC_OUTPUT}"

