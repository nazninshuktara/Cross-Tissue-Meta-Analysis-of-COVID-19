#!/usr/bin/env bash

# Run quality check with FastQC and MultiQC
# Dataset: GSE211979

# Input: FASTQ files in fastq/GSE211979/*.fastq.gz
# Output: QC reports in outputs/qc/GSE211979/

set -euo pipefail

PROJECT_ROOT="../../../"

INPUT="${PROJECT_ROOT}/fastq/GSE211979"
OUTPUT="${PROJECT_ROOT}/outputs/qc/GSE211979"

mkdir -p "$OUTPUT"

echo "Running FastQC..."
fastqc \
    "${INPUT}"/*.fastq.gz \
    --outdir "$OUTPUT"

echo "Running MultiQC..."
multiqc \
    "$OUTPUT" \
    --outdir "$OUTPUT"

echo "Quality check complete! Reports saved in: $OUTPUT"
