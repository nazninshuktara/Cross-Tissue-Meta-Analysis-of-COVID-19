#!/usr/bin/env bash

# Run quality check with FastQC and MultiQC
# Dataset: GSE243217

# Input: FASTQ files in fastq/GSE243217/*.fastq.gz
# Output: QC reports in outputs/qc/GSE243217/

set -euo pipefail

PROJECT_ROOT="../../../"

INPUT="${PROJECT_ROOT}/fastq/GSE243217"
OUTPUT="${PROJECT_ROOT}/outputs/qc/GSE243217"

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
