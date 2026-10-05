#!/usr/bin/env Rscript

# 05_tximport_gene_counts.R
# Convert Salmon transcript-level quantification to gene-level counts
# using tximport and GENCODE v49 transcript-to-gene mapping.

# Load required packages
suppressPackageStartupMessages({
  library(tximport)  # Import Salmon quantification and summarize to genes
  library(readr)     # Read and write CSV files
  library(dplyr)     # Filter and manipulate metadata
  library(tibble)    # Handle data frames and gene IDs
})

# Define project paths
PROJECT_ROOT <- "../../"

METADATA_FILE <- file.path(
  PROJECT_ROOT,
  "data",
  "metadata",
  "metadata_curated_sample.csv"
)

QUANT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "salmon_out"
)

TX2GENE_FILE <- file.path(
  PROJECT_ROOT,
  "inputs",
  "reference",
  "tx2gene_gencode_v49.csv"
)

OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "data",
  "gene_counts"
)

# Create output directory
dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)