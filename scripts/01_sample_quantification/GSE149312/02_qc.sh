#!/usr/bin/env bash

# Run quality check with FastQC and MultiQC
# Dataset: GSE149312

# Input: FASTQ files in fastq/GSE149312/*.fastq.gz
# Output: QC reports in outputs/qc/GSE149312/

set -euo pipefail

PROJECT_ROOT="../../../"

INPUT="${PROJECT_ROOT}/fastq/GSE149312"
OUTPUT="${PROJECT_ROOT}/outputs/qc/GSE149312"

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
