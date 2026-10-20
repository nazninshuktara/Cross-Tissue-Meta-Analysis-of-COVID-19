#!/usr/bin/env Rscript

# 02_GSE208342_run_level_DE.R
# Verify GSE208342 differential expression using individual SRR runs.

# Load required packages
suppressPackageStartupMessages({
  library(tximport)    # Summarize transcript-level estimates to gene-level
  library(DESeq2)      # RNA-seq count processing and VST transformation
  library(readr)       # Read and write CSV files
  library(dplyr)       # Data manipulation
})

# Define project paths
PROJECT_ROOT <- "."
DATASET <- "GSE208342"

# Run-level metadata.
METADATA_FILE <- file.path(
  PROJECT_ROOT, 
  "data", 
  "metadata", 
  "metadata_run_level.csv"
)

# GENCODE v49 protein-coding transcript-to-gene mapping
TX2GENE_FILE <- file.path(
  PROJECT_ROOT, 
  "inputs", 
  "reference", 
  "tx2gene_gencode_v49.csv"
)

# Salmon quant
QUANT_DIR <- file.path(
  PROJECT_ROOT, 
  "outputs", 
  "salmon_out", 
  DATASET
)

# Directory for storing differential expression analysis results.
OUTPUT_DIR <- file.path(
  PROJECT_ROOT, 
  "outputs", 
  "differential_expression"
)

# Create the output directory if it does not already exist.
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Read run-level metadata
metadata <- read_csv(METADATA_FILE, show_col_types = FALSE) %>%
  filter(dataset_id == DATASET)

# Define condition from curated disease-group labels
metadata <- metadata %>%
  mutate(
    condition = factor(
      disease_group,
      levels = c("Control", "COVID-19")
    )
  )

# Check run counts by condition
print(dplyr::count(metadata, condition))

# Check required columns
required <- c("run_id", "condition")
if (!all(required %in% names(metadata))) {
  stop(
    "Check metadata column names. Required: run_id and condition. ",
    "Available: ", paste(names(metadata), collapse = ", ")
  )
}

metadata <- metadata %>%
  mutate(
    run_id = as.character(run_id),
    condition = as.character(condition)
  ) %>%
  filter(condition %in% c("COVID-19", "Control"))

if (anyDuplicated(metadata$run_id)) {
  stop("Duplicate SRR run IDs found.")
}

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$run_id

# Locate Salmon quantification files
quant_files <- file.path(
  QUANT_DIR, metadata$run_id, "quant.sf"
)
names(quant_files) <- metadata$run_id

if (any(!file.exists(quant_files))) {
  stop(
    "Missing quant.sf files: ",
    paste(names(quant_files)[!file.exists(quant_files)], collapse = ", ")
  )
}

# Read transcript-to-gene mapping
tx2gene <- read_csv(TX2GENE_FILE, show_col_types = FALSE) %>%
  select(transcript_id, gene_id) %>%
  mutate(
    transcript_id = sub("\\.[0-9]+$", "", transcript_id)
  )

if (anyDuplicated(tx2gene$transcript_id)) {
  stop("Duplicate transcript IDs in tx2gene mapping.")
}

# Import each SRR as a separate column and summarize to genes
txi <- tximport(
  quant_files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "no",
  ignoreTxVersion = TRUE
)

# Confirm run order matches metadata
metadata <- metadata[
  match(colnames(txi$counts), metadata$run_id),
]

if (anyNA(metadata$run_id)) {
  stop("Imported runs do not match metadata.")
}

# Build DESeq2 dataset: each SRR is treated as an independent sample
metadata$condition <- factor(
  metadata$condition,
  levels = c("Control", "COVID-19")
)

rownames(metadata) <- metadata$run_id

dds <- DESeqDataSetFromMatrix(
  countData = round(txi$counts),
  colData = metadata,
  design = ~ condition
)

# Apply the same low-count filter used in the current workflow
keep <- rowSums(counts(dds) >= 10) >= ceiling(0.10 * ncol(dds))
dds <- dds[keep, ]

# Run differential expression
dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "COVID-19", "Control"),
  alpha = 0.05
)

result_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene_id") %>%
  mutate(
    significant = !is.na(padj) & padj < 0.05
  ) %>%
  arrange(padj)

# Save results and summary
write_csv(
  result_df,
  file.path(OUTPUT_DIR, "GSE208342_DE_results.csv")
)

summary_df <- tibble::tibble(
  dataset_id = DATASET,
  n_runs = ncol(dds),
  n_covid_runs = sum(metadata$condition == "COVID-19"),
  n_control_runs = sum(metadata$condition == "Control"),
  genes_tested = nrow(result_df),
  significant_genes = sum(result_df$significant, na.rm = TRUE)
)

write_csv(
  summary_df,
  file.path(OUTPUT_DIR, "GSE208342_summary.csv")
)

print(summary_df)
