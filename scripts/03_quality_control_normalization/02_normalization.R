#!/usr/bin/env Rscript

# 02_normalization.R
# Estimate DESeq2 size factors and generate normalized gene-level counts.

# Load required packages
suppressPackageStartupMessages({
  library(DESeq2)   # RNA-seq count processing and VST transformation
  library(dplyr)    # Data manipulation
  library(readr)    # Read and write CSV files
  library(tibble)   # Tidy data frame utilities
})

# Define project root
PROJECT_ROOT <- "."

# Input directories and files
# Gene counts
COUNT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "gene_counts"
)

# Curated sample-level metadata.
METADATA_FILE <- file.path(
  PROJECT_ROOT,
  "data",
  "metadata",
  "metadata_curated_sample.csv"
)

# Output directory
# Directory for storing gene-level count normalization results.
OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "gene_counts_normalization"
)

# Create the output directory if it does not already exist.
dir.create(
  OUTPUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

# Define discovery datasets
DATASETS <- c(
  "GSE208342",
  "GSE273149",
  "GSE243217",
  "GSE202182",
  "GSE208076",
  "GSE211979",
  "GSE169241",
  "GSE149312",
  "GSE182917"
)

# Samples excluded after quality control
QC_EXCLUDED_SAMPLES <- c(
  "GSM6896317"
)

# Read curated sample-level metadata
metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

# Keep included primary analysis samples
metadata <- metadata %>%
  filter(
    inclusion_status == "included",
    dataset_id %in% DATASETS,
    !sample_id %in% QC_EXCLUDED_SAMPLES
  )

if (anyDuplicated(metadata$sample_id)) {
  stop("Duplicate sample IDs found in metadata.")
}

if (nrow(metadata) == 0) {
  stop("No included samples found.")
}

# Report sample counts before normalization
cat("\nPrimary analysis samples by dataset:\n")

sample_summary <- metadata %>%
  dplyr::count(dataset_id, name = "n_samples")

print(sample_summary)

cat(
  "\nTotal primary analysis samples:",
  nrow(metadata),
  "\n"
)

# Initialize output containers
size_factor_list <- list()
normalization_summary_list <- list()

# Start normalization
cat("\nPerforming DESeq2 size-factor normalization...\n")

# Dataset-wise normalization
for (dataset in DATASETS) {
  
  cat("\nProcessing:", dataset, "\n")
  
  # Define gene-count input file
  count_file <- file.path(
    COUNT_DIR,
    paste0("gene_counts_", dataset, ".csv")
  )
  
  if (!file.exists(count_file)) {
    stop("Missing gene-count file: ", count_file)
  }
  
  # Read gene-count matrix
  counts_df <- read_csv(
    count_file,
    show_col_types = FALSE
  )
  names(counts_df)[1] <- "gene_id"
  
  # Get included samples for this dataset
  dataset_metadata <- metadata %>%
    filter(dataset_id == dataset)
  
  metadata_samples <- dataset_metadata$sample_id
  
  # Retain only included primary analysis samples
  counts_df <- counts_df %>%
    select(
      gene_id,
      all_of(metadata_samples)
    )
  
  # Convert counts to a numeric matrix
  count_matrix <- as.matrix(
    counts_df[, -1]
  )
  
  rownames(count_matrix) <- counts_df$gene_id
  
  # Round estimated counts for DESeq2
  count_matrix <- round(count_matrix)
  storage.mode(count_matrix) <- "integer"
  
  # Filter low-count genes
  min_samples <- max(
    1L,
    ceiling(
      0.10 * ncol(count_matrix)
    )
  )
  
  keep_genes <- rowSums(
    count_matrix >= 10
  ) >= min_samples
  
  count_matrix <- count_matrix[
    keep_genes,
    ,
    drop = FALSE
  ]
  
  if (nrow(count_matrix) == 0) {
    stop("No genes remain after filtering in ", dataset)
  }
  
  cat("Samples:", ncol(count_matrix), "\n")
  cat("Genes after filtering:", nrow(count_matrix), "\n")
  
  # DESeq2 size factor estimation
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = data.frame(
      row.names = colnames(count_matrix)
    ),
    design = ~ 1
  )
  
  # Estimate size factors
  dds <- estimateSizeFactors(dds)
  
  size_factors <- sizeFactors(dds)
  
  if (any(!is.finite(size_factors)) ||
      any(size_factors <= 0)) {
    stop("Invalid size factors estimated for ", dataset)
  }
  
  # Store sample-level size factors
  size_factor_list[[dataset]] <- tibble(
    dataset_id = dataset,
    sample_id = names(size_factors),
    size_factor = as.numeric(size_factors)
  )
  
  cat(
    "Size factors: min =",
    round(min(size_factors), 3),
    "median =",
    round(median(size_factors), 3),
    "max =",
    round(max(size_factors), 3),
    "\n"
  )
  
  # Generate normalized counts
  normalized_counts <- counts(
    dds,
    normalized = TRUE
  )
  
  normalized_counts_df <- as.data.frame(
    normalized_counts
  ) %>%
    rownames_to_column("gene_id")
  
  # Save normalized counts
  normalized_file <- file.path(
    OUTPUT_DIR,
    paste0("normalized_counts_", dataset, ".csv")
  )
  
  write_csv(
    normalized_counts_df,
    normalized_file
  )
  
  cat("Saved:", normalized_file, "\n")
  
  # Store dataset-level summary
  normalization_summary_list[[dataset]] <- tibble(
    dataset_id = dataset,
    n_samples = ncol(count_matrix),
    n_genes = nrow(count_matrix),
    min_size_factor = min(size_factors),
    median_size_factor = median(size_factors),
    max_size_factor = max(size_factors)
  )
}

# Combine size-factor tables
size_factors_all <- bind_rows(size_factor_list)

# Save size factors
write_csv(
  size_factors_all,
  file.path(OUTPUT_DIR, "size_factors.csv")
)

# Combine dataset-level normalization summaries
normalization_summary <- bind_rows(
  normalization_summary_list
)

# Save normalization summary
write_csv(
  normalization_summary,
  file.path(OUTPUT_DIR, "normalization_summary.csv")
)

# Report final summary
cat("\nNormalization completed successfully.\n\n")

cat("Datasets:", length(DATASETS), "\n")
cat("Samples normalized:", nrow(size_factors_all), "\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")
