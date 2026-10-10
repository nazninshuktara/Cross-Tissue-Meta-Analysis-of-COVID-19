#!/usr/bin/env Rscript

# 06_tximport_gene_counts.R
# Import Salmon quantification and aggregate technical runs by biological sample.

# Load required packages
suppressPackageStartupMessages({
  library(tximport)    # Import Salmon quantifications and summarize transcript-level estimates to gene-level
  library(readr)       # Read metadata and write gene-level count and length tables
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
  "tx2gene_gencode_v49_pc.csv"
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
  "data", 
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
metadata <- read_csv(METADATA_FILE, show_col_types = FALSE)

# Check required metadata columns
required_columns <- c(
  "dataset_id", "sample_id", "run_id", "inclusion_status"
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

if (anyDuplicated(tx2gene$transcript_id)) {
  stop("Mapping contains duplicate transcript IDs.")
}

tx2gene <- as.data.frame(tx2gene)

# Remove transcript version suffixes to match Salmon quantification IDs
tx2gene$transcript_id <- sub(
  "\\.[0-9]+$",
  "",
  tx2gene$transcript_id
)

# Validate transcript IDs after removing version suffixes

if (anyDuplicated(tx2gene$transcript_id)) {
stop("Duplicate transcript IDs found after removing version suffixes.")
}

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

  # Import Salmon transcript-level estimates and summarize by gene
  txi <- tximport(
    quant_files,
    type = "salmon",
    tx2gene = tx2gene,
    countsFromAbundance = "no",
    ignoreTxVersion = TRUE
  )

  # Confirm that imported run columns match the run map
  if (!all(run_ids %in% colnames(txi$counts))) {
    stop(dataset, ": imported count columns do not match run IDs.")
  }
  
  sample_ids <- dataset_metadata$sample_id
  
  sample_counts <- matrix(
    0,
    nrow = nrow(txi$counts),
    ncol = length(sample_ids),
    dimnames = list(rownames(txi$counts), sample_ids)
  )
  
  sample_lengths <- sample_counts
  
  # Aggregate technical runs within each biological sample
  for (sample in sample_ids) {
    sample_runs <- dataset_map$run_id[
      dataset_map$sample_id == sample
    ]
    
    run_counts <- txi$counts[, sample_runs, drop = FALSE]
    run_lengths <- txi$length[, sample_runs, drop = FALSE]

    # Sum estimated counts across runs
    aggregated_counts <- rowSums(run_counts)

    # Count-weighted effective length across runs
    valid_lengths <- is.finite(run_lengths) & run_lengths > 0
    
    weighted_counts <- run_counts
    weighted_counts[!valid_lengths] <- 0

    length_values <- run_lengths
    length_values[!valid_lengths] <- 0

    numerator <- rowSums(weighted_counts * length_values)
    denominator <- rowSums(weighted_counts)

    aggregated_lengths <- rep(NA_real_, length(aggregated_counts))
    names(aggregated_lengths) <- rownames(txi$counts)

    has_counts <- denominator > 0
    aggregated_lengths[has_counts] <-
      numerator[has_counts] / denominator[has_counts]

    # For zero-count genes, use the mean valid run length
    no_counts <- !has_counts

    if (any(no_counts)) {
      fallback_lengths <- run_lengths[no_counts, , drop = FALSE]
      fallback_lengths[!is.finite(fallback_lengths) |
                         fallback_lengths <= 0] <- NA_real_
      
      fallback_means <- rowMeans(fallback_lengths, na.rm = TRUE)
      fallback_means[!is.finite(fallback_means)] <- NA_real_

      aggregated_lengths[no_counts] <- fallback_means
    }
    

    if (any(!is.finite(aggregated_lengths)) ||
        any(aggregated_lengths <= 0)) {
      stop(
        dataset, ": invalid effective lengths after aggregation for sample ",
        sample
      )
    }

    sample_counts[, sample] <- aggregated_counts
    sample_lengths[, sample] <- aggregated_lengths
  }

  # Convert gene-level estimated counts to integers for downstream DESeq2
  gene_counts <- round(sample_counts)
  storage.mode(gene_counts) <- "integer"
  
  # Define output files
  count_file <- file.path(
    OUTPUT_DIR, paste0("gene_counts_", dataset, ".csv")
  )
  
  length_file <- file.path(
    OUTPUT_DIR, paste0("gene_lengths_", dataset, ".csv")
  )
  
  # Avoid silently overwriting existing results
  if (file.exists(count_file) || file.exists(length_file)) {
    stop(
      "Output already exists for ", dataset,
      ". Move or back up existing files before rerunning: ",
      count_file, " ; ", length_file
    )
  }
  
  # Write gene-level counts and effective lengths
  counts_df <- data.frame(
    gene_id = rownames(gene_counts),
    gene_counts,
    check.names = FALSE
  )
  
  lengths_df <- data.frame(
    gene_id = rownames(sample_lengths),
    sample_lengths,
    check.names = FALSE
  )
  
  write_csv(counts_df, count_file)
  write_csv(lengths_df, length_file)
  
  message("  Gene-level counts saved: ", count_file)
  message("  Effective lengths saved: ", length_file)
  message("  Genes: ", nrow(gene_counts))
  message("  Samples: ", ncol(gene_counts))
}

message("\nAll datasets processed successfully.")

