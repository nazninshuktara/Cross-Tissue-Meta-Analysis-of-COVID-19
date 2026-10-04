#!/usr/bin/env bash

# Run quality check with FastQC and MultiQC
# Dataset: GSE169241

# Input: FASTQ files in fastq/GSE169241/*.fastq.gz
# Output: QC reports in outputs/qc/GSE169241/

set -euo pipefail

PROJECT_ROOT="../../../"

INPUT="${PROJECT_ROOT}/fastq/GSE169241"
OUTPUT="${PROJECT_ROOT}/outputs/qc/GSE169241"

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
