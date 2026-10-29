#!/usr/bin/env Rscript

# 02_tximport_gene_counts.R
# Import Salmon quantification, summarize transcript-level estimates at the gene level,
# and aggregate technical runs for each biological sample.

# Load required packages
suppressPackageStartupMessages({
  library(tximport)    # Import Salmon quantifications and summarize transcript-level estimates to gene-level
  library(readr)       # Read metadata and write gene-level count
  library(dplyr)       # Filter, select, and organize metadata and sequencing-run information
})

# Define project paths
PROJECT_ROOT <- "."

# Curated sample-level metadata.
METADATA_FILE <- file.path(
  PROJECT_ROOT, 
  "data", 
  "metadata",
  "metadata_curated_sample.csv"
)

# GENCODE v49 protein-coding transcript-to-gene mapping
TX2GENE_FILE <- file.path(
  PROJECT_ROOT, 
  "inputs", 
  "reference",
  "tx2gene_gencode_v49.csv"
)

# Directory containing Salmon quantification results
QUANT_DIR <- file.path(
  PROJECT_ROOT, 
  "outputs", 
  "salmon_out"
)

# Directory for storing gene-level count and length tables
OUTPUT_DIR <- file.path(
  PROJECT_ROOT, 
  "outputs", 
  "gene_counts"
)

# Create the output directory if it does not already exist.
dir.create(
  OUTPUT_DIR, 
  recursive = TRUE, 
  showWarnings = FALSE
)

# Validate input files
if (!file.exists(METADATA_FILE)) {
  stop("Metadata file not found: ", METADATA_FILE)
}

if (!file.exists(TX2GENE_FILE)) {
  stop("Transcript-to-gene mapping not found: ", TX2GENE_FILE)
}

# Read metadata and select primary-analysis samples
metadata <- read_csv(
  METADATA_FILE, 
  show_col_types = FALSE
)

# Check required metadata columns
required_columns <- c(
  "dataset_id", 
  "sample_id", 
  "run_id", 
  "inclusion_status"
)

missing_columns <- setdiff(required_columns, names(metadata))

if (length(missing_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Retain included samples and standardize identifier types
metadata <- metadata %>%
  filter(inclusion_status == "included") %>%
  mutate(
    dataset_id = as.character(dataset_id),
    sample_id = as.character(sample_id),
    run_id = as.character(run_id)
)

# Ensure that included samples are available
if (nrow(metadata) == 0) {
  stop("No included biological samples found in metadata.")
}

# Check for missing identifiers
if (anyNA(metadata$dataset_id) ||
    anyNA(metadata$sample_id) ||
    anyNA(metadata$run_id)) {
  stop("Included metadata contains missing dataset, sample, or run IDs.")
}

# Each dataset/sample combination must identify one biological sample
if (anyDuplicated(metadata[c("dataset_id", "sample_id")])) {
  stop("Duplicate dataset_id + sample_id rows found in included metadata.")
}

# Expand comma-separated run IDs into one row per sequencing run
run_map <- do.call(
  rbind,
  lapply(seq_len(nrow(metadata)), function(i) {
    runs <- trimws(
      strsplit(as.character(metadata$run_id[i]), ",", fixed = TRUE)[[1]]
    )
    runs <- runs[nzchar(runs)]

    if (length(runs) == 0) {
      stop("No run IDs found for sample: ", metadata$sample_id[i])
    }

    data.frame(
      dataset_id = as.character(metadata$dataset_id[i]),
      sample_id = as.character(metadata$sample_id[i]),
      run_id = runs,
      stringsAsFactors = FALSE
    )
  })
)

# Ensure that each sequencing run is assigned only once per dataset
if (anyDuplicated(run_map[c("dataset_id", "run_id")])) {
  stop("A sequencing run is assigned more than once within a dataset.")
}

# Read and validate transcript-to-gene mapping
tx2gene <- read_csv(
  TX2GENE_FILE, 
  show_col_types = FALSE
)

if (!all(c("transcript_id", "gene_id") %in% names(tx2gene))) {
  stop("Mapping must contain transcript_id and gene_id columns.")
}

if (anyNA(tx2gene$transcript_id) || anyNA(tx2gene$gene_id)) {
  stop("Mapping contains missing transcript or gene IDs.")
}

# Remove transcript version suffixes to match Salmon quantification IDs
tx2gene$transcript_id <- sub(
  "\\.[0-9]+$",
  "",
  tx2gene$transcript_id
)

if (anyDuplicated(tx2gene$transcript_id)) {
  stop("Duplicate transcript IDs found after removing version suffixes.")
}

tx2gene <- as.data.frame(tx2gene)

# Report expected primary-analysis input size
message("Included biological samples: ", nrow(metadata))
message("Included sequencing runs: ", nrow(run_map))
message("Datasets to process: ", n_distinct(metadata$dataset_id))
message("Mapping transcript IDs: ", nrow(tx2gene))

# Process each dataset independently
for (dataset in unique(metadata$dataset_id)) {
  message("\nProcessing dataset: ", dataset)

  dataset_map <- run_map %>%
    filter(dataset_id == dataset)

  dataset_metadata <- metadata %>%
    filter(dataset_id == dataset)

  run_ids <- dataset_map$run_id

  quant_files <- file.path(
    QUANT_DIR, 
    dataset, 
    run_ids, 
    "quant.sf"
  )
  names(quant_files) <- run_ids

  # Stop rather than silently dropping missing runs
  missing_files <- quant_files[!file.exists(quant_files)]

  if (length(missing_files) > 0) {
    stop(
      dataset, ": missing quant.sf files for: ",
      paste(names(missing_files), collapse = ", ")
    )
  }

  message("  Runs imported: ", length(quant_files))
  message("  Biological samples: ", nrow(dataset_metadata))
  
  # Import Salmon quantifications and summarize transcripts to genes
  txi <- tximport(
    quant_files,
    type = "salmon",
    tx2gene = tx2gene,
    ignoreTxVersion = TRUE,
    countsFromAbundance = "no"
  )
  
  # Group technical runs by biological sample within this dataset.
  dataset_run_map <- dataset_map |>
    dplyr::select(sample_id, run_id)
  
  sample_run_groups <- split(
    dataset_run_map$run_id,
    dataset_run_map$sample_id
  )
  
  # Confirm that every run maps to an imported quantification column.
  if (!all(dataset_run_map$run_id %in% colnames(txi$counts))) {
    stop(dataset, ": run IDs do not match quantification columns.")
  }
  
  # Sum gene-level counts across technical runs
  gene_counts <- vapply(
    sample_run_groups,
    function(sample_runs) {
      rowSums(
        txi$counts[, sample_runs, drop = FALSE]
      )
    },
    FUN.VALUE = numeric(nrow(txi$counts))
  )
  rownames(gene_counts) <- rownames(txi$counts)
  
  # Round estimated counts for DESeq2
  gene_counts <- round(gene_counts)
  storage.mode(gene_counts) <- "integer"
  
  # Save dataset-level gene count matrix
  readr::write_csv(
    data.frame(
      gene_id = rownames(gene_counts),
      gene_counts,
      check.names = FALSE
    ),
    file.path(OUTPUT_DIR, paste0("gene_counts_", dataset, ".csv"))
  )
  
  message("  Gene-level count matrix: ",
          nrow(gene_counts), " genes x ",
          ncol(gene_counts), " biological samples")
}
message("\nAll datasets processed successfully.")

# Combine all datasets
# Identify dataset-specific count files
count_files <- file.path(
  OUTPUT_DIR,
  paste0("gene_counts_", unique(metadata$dataset_id), ".csv")
)

# Read all dataset-specific count matrices
count_tables <- lapply(
  count_files, 
  readr::read_csv, 
  show_col_types = FALSE
)

# Add dataset prefixes to sample column names
count_tables <- lapply(
  seq_along(count_tables),
  function(i) {
    counts <- count_tables[[i]]
    dataset <- unique(metadata$dataset_id)[i]
    
    sample_cols <- setdiff(names(counts), "gene_id")
    names(counts)[match(sample_cols, names(counts))] <-
      paste(dataset, sample_cols, sep = "_")
    
    counts
  }
)

# Combine datasets by gene ID
combined_counts <- Reduce(
  function(x, y) dplyr::full_join(x, y, by = "gene_id"),
  count_tables
)

# Replace missing values with zero for genes absent from a dataset
sample_cols <- setdiff(names(combined_counts), "gene_id")
combined_counts[sample_cols] <- lapply(
  combined_counts[sample_cols],
  function(x) {
    x[is.na(x)] <- 0
    x
  }
)

# Save the combined gene-level count matrix
readr::write_csv(
  combined_counts,
  file.path(OUTPUT_DIR, "gene_counts_all_datasets.csv")
)

message(
  "Combined gene count matrix: ",
  nrow(combined_counts), " genes x ",
  length(sample_cols), " biological samples"
)
message("\nGene-level count processing completed.")

