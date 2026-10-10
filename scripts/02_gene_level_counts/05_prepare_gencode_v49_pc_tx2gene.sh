#!/usr/bin/env bash

# 05_prepare_gencode_v49_pc_tx2gene.sh
# Create a transcript-to-gene mapping from GENCODE v49 protein-coding transcripts.

set -euo pipefail

# Define project and reference directories
PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
REFERENCE_DIR="${PROJECT_ROOT}/inputs/reference"

# Define reference files
TRANSCRIPT_FASTA="${REFERENCE_DIR}/gencode.v49.pc_transcripts.fa.gz"
TX2GENE_FILE="${REFERENCE_DIR}/tx2gene_gencode_v49_pc.csv"

# Create reference directory if needed
mkdir -p "${REFERENCE_DIR}"

# Check transcript FASTA
if [[ ! -s "${TRANSCRIPT_FASTA}" ]]; then
    echo "ERROR: Transcript FASTA not found: ${TRANSCRIPT_FASTA}" >&2
    exit 1
fi

# Validate compressed FASTA
gzip -t "${TRANSCRIPT_FASTA}"

# Create transcript-to-gene mapping from FASTA headers
{
    printf 'transcript_id,gene_id\n'

    zcat "${TRANSCRIPT_FASTA}" |
    awk -F'|' '/^>/ {
        sub(/^>/, "", $1)

        transcript_id = $1
        gene_id = $2

        if (transcript_id != "" && gene_id != "") {
            print transcript_id "," gene_id
        }
    }' |
    sort -u
} > "${TX2GENE_FILE}"

# Validate mapping
MAPPING_COUNT=$(($(wc -l < "${TX2GENE_FILE}") - 1))

DUPLICATE_COUNT=$(
    tail -n +2 "${TX2GENE_FILE}" |
    cut -d',' -f1 |
    sort |
    uniq -d |
    wc -l
)

if [[ "${MAPPING_COUNT}" -le 0 || "${DUPLICATE_COUNT}" -ne 0 ]]; then
    echo "ERROR: Mapping validation failed." >&2
    echo "Mappings: ${MAPPING_COUNT}" >&2
    echo "Duplicate transcript IDs: ${DUPLICATE_COUNT}" >&2
    exit 1
fi

echo "Created: ${TX2GENE_FILE}"
echo "Mappings: ${MAPPING_COUNT}"
echo "Duplicate transcript IDs: ${DUPLICATE_COUNT}"
echo "Validation: OK"
