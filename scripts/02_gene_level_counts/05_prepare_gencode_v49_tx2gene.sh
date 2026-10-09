#!/usr/bin/env bash

# 05_prepare_gencode_v49_tx2gene.sh
# Download GENCODE v49 annotation and create a transcript-to-gene mapping.

set -euo pipefail

# Define reference directory
REFERENCE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/inputs/reference"

# Define reference files
GTF_FILE="${REFERENCE_DIR}/gencode.v49.annotation.gtf.gz"
TX2GENE_FILE="${REFERENCE_DIR}/tx2gene_gencode_v49.csv"

# Define GENCODE v49 annotation URL
GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"

# Create reference directory if needed
mkdir -p "${REFERENCE_DIR}"

# Download GENCODE v49 annotation
if [[ ! -s "${GTF_FILE}" ]]; then
    wget -c "${GTF_URL}" -O "${GTF_FILE}"
fi

# Create version-free transcript-to-gene mapping
zcat "${GTF_FILE}" | \
awk '$3 == "transcript" {
    match($0, /gene_id "([^"]+)"/, gene)
    match($0, /transcript_id "([^"]+)"/, tx)

    if (gene[1] != "" && tx[1] != "") {

        tx_id = tx[1]
        gene_id = gene[1]

        sub(/\.[0-9]+$/, "", tx_id)
        sub(/\.[0-9]+$/, "", gene_id)

        print tx_id "," gene_id
    }
}' | \
sort -u | \
{ echo "transcript_id,gene_id"; cat; } > "${TX2GENE_FILE}"

# Report mapping count
MAPPING_COUNT=$(($(wc -l < "${TX2GENE_FILE}") - 1))

echo "Created: ${TX2GENE_FILE}"
echo "Mappings: ${MAPPING_COUNT}"