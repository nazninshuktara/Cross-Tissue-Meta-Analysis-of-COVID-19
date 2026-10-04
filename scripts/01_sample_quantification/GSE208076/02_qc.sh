#!/usr/bin/env bash

# Run quality check with FastQC and MultiQC
# Dataset: GSE208076

# Input: FASTQ files in fastq/GSE208076/*.fastq.gz
# Output: QC reports in outputs/qc/GSE208076/

set -euo pipefail

PROJECT_ROOT="../../../"

INPUT="${PROJECT_ROOT}/fastq/GSE208076"
OUTPUT="${PROJECT_ROOT}/outputs/qc/GSE208076"

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
