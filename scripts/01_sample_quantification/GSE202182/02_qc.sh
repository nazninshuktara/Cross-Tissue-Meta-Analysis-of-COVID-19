#!/usr/bin/env bash

# Run quality check with FastQC and MultiQC
# Dataset: GSE202182

# Input: FASTQ files in fastq/GSE202182/*.fastq.gz
# Output: QC reports in outputs/qc/GSE202182/

set -euo pipefail

PROJECT_ROOT="../../../"

INPUT="${PROJECT_ROOT}/fastq/GSE202182"
OUTPUT="${PROJECT_ROOT}/outputs/qc/GSE202182"

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
