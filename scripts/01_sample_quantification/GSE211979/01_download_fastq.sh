#!/usr/bin/env bash

# Download FASTQ files using parallel-fastq-dump
# Dataset: GSE211979

# Reads accessions from inputs/discovery/GSE211979_SRR_Acc_List.txt
# Saves all FASTQ files in the root-level fastq/GSE211979/ directory

# Example:
# parallel-fastq-dump --sra-id SRR24919175 --threads 4 --outdir out/ --split-files --gzip

# Data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211979

set -euo pipefail

PROJECT_ROOT="../../../"

INPUT="${PROJECT_ROOT}/inputs/discovery/GSE211979_SRR_Acc_List.txt"
OUTPUT="${PROJECT_ROOT}/fastq/GSE211979"

mkdir -p "$OUTPUT"

while read -r ACC; do
    [[ -z "$ACC" ]] && continue
    echo "Downloading $ACC ..."
    
    parallel-fastq-dump \
        --sra-id "$ACC" \
        --threads 4 \
        --outdir "$OUTPUT" \
        --split-files \
        --gzip
done < "$INPUT"
echo "All FASTQ files are saved in: $OUTPUT"