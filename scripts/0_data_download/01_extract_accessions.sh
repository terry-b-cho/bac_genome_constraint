#!/bin/bash
# Extract assembly accessions from TSV file

INPUT_TSV="/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint/data/0_20241106_ncbi_ref_genome_table.tsv"
OUTPUT_TXT="/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint/data/metadata/selected_genomes.txt"
METADATA_DIR="/n/data1/joslin/icrb/kostic/terry/github/bac_genome_constraint/data/metadata"

mkdir -p "${METADATA_DIR}"

tail -n +2 "${INPUT_TSV}" | cut -f1 > "${OUTPUT_TXT}"

echo "Extracted $(wc -l < ${OUTPUT_TXT}) accessions to ${OUTPUT_TXT}"



