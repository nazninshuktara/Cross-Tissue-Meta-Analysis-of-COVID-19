#!/usr/bin/env bash

# Run quality check with FastQC and MultiQC
# Dataset: GSE182917

# Input: FASTQ files in fastq/GSE182917/*.fastq.gz
# Output: QC reports in outputs/qc/GSE182917/

set -euo pipefail

PROJECT_ROOT="../../../"

INPUT="${PROJECT_ROOT}/fastq/GSE182917"
OUTPUT="${PROJECT_ROOT}/outputs/qc/GSE182917"

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
