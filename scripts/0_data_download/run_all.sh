#!/bin/bash
# Master script to run all download steps sequentially

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${SCRIPT_DIR}"

echo "Step 1: Extracting accessions..."
bash 01_extract_accessions.sh

echo "Step 2: Submitting download job..."
JOB1=$(sbatch --parsable 02_download_dehydrated.sh)
echo "Download job ID: ${JOB1}"

echo "Step 3: Submitting rehydration job (depends on step 2)..."
JOB2=$(sbatch --parsable --dependency=afterok:${JOB1} 03_rehydrate.sh)
echo "Rehydration job ID: ${JOB2}"

echo "Step 4: Submitting organization job (depends on step 3)..."
JOB3=$(sbatch --parsable --dependency=afterok:${JOB2} 04_organize.sh)
echo "Organization job ID: ${JOB3}"

echo "Step 5: Submitting QC job (depends on step 4)..."
JOB4=$(sbatch --parsable --dependency=afterok:${JOB3} 05_qc.sh)
echo "QC job ID: ${JOB4}"

echo "Step 6: Submitting metadata extraction job (depends on step 2)..."
JOB5=$(sbatch --parsable --dependency=afterok:${JOB1} 06_extract_metadata.sh)
echo "Metadata job ID: ${JOB5}"

echo "Step 7: Submitting assembly data processing job (depends on step 6)..."
JOB6=$(sbatch --parsable --dependency=afterok:${JOB5} 07_assembly_data_processing.sh)
echo "Assembly data processing job ID: ${JOB6}"

echo "All jobs submitted. Monitor with: squeue -u \$USER"



