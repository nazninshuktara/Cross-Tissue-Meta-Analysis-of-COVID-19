#!/usr/bin/env Rscript

# 09_differential_expression.R
# Perform per-dataset differential expression analysis using DESeq2.

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
# Directory for storing per-dataset differential expression analysis results.
OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "differential_expression"
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

# Samples excluded from the primary analysis after QC
QC_EXCLUDED_SAMPLES <- c("GSM6896317")

# Read curated sample-level metadata
metadata <- read_csv(
  METADATA_FILE,
  show_col_types = FALSE
) %>%
  dplyr::filter(
    inclusion_status == "included",
    dataset_id %in% DATASETS,
    !sample_id %in% QC_EXCLUDED_SAMPLES
  ) %>%
  dplyr::select(
    dataset_id,
    sample_id,
    tissue_type,
    disease_status,
    disease_group,
    disease_severity,
    inclusion_status,
    age,
    sex
)

# Validate primary metadata
if (nrow(metadata) == 0) {
  stop("No primary samples found in curated metadata.")
}

if (anyDuplicated(metadata$sample_id)) {
  stop("Duplicate sample_id values found in curated metadata.")
}

if (!all(DATASETS %in% unique(metadata$dataset_id))) {
  stop("One or more discovery datasets are missing from curated metadata.")
}

# Check disease groups
unexpected_groups <- setdiff(
  unique(metadata$disease_group),
  c("COVID-19", "Control")
)

if (length(unexpected_groups) > 0) {
  stop(
    "Unexpected disease_group values found: ",
    paste(unexpected_groups, collapse = ", ")
  )
}

# Initialize result summary
result_summary <- list()

# Process each dataset
for (dataset in DATASETS) {

  cat(
    "\nProcessing differential expression:",
    dataset,
    "\n"
  )

  # Select dataset metadata
  dataset_metadata <- metadata %>%
    filter(
      dataset_id == dataset
    )

  # Check disease groups
  group_counts <- table(
    dataset_metadata$disease_group
  )

  if (!all(
    c("COVID-19", "Control") %in%
      names(group_counts)
  )) {
    stop(
      dataset,
      ": both COVID-19 and Control groups are required."
    )
  }

  cat(
    "COVID-19:",
    group_counts["COVID-19"],
    "\n"
  )

  cat(
    "Control:",
    group_counts["Control"],
    "\n"
  )

  # Read gene-count matrix
  count_file <- file.path(
    COUNT_DIR,
    paste0(
      "gene_counts_",
      dataset,
      ".csv"
    )
  )

  if (!file.exists(count_file)) {
    stop(
      "Missing gene-count file: ",
      count_file
    )
  }

  counts_df <- read_csv(
    count_file,
    show_col_types = FALSE
  ) %>%
    rename(
      gene_id = 1
    )

  # Check duplicate genes
  if (anyDuplicated(counts_df$gene_id)) {
    stop(
      dataset,
      ": duplicate gene IDs found."
    )
  }

  # Keep only primary metadata samples
  sample_ids <- dataset_metadata$sample_id

  missing_samples <- setdiff(
    sample_ids,
    colnames(counts_df)
  )

  if (length(missing_samples) > 0) {
    stop(
      dataset,
      ": samples present in metadata but missing from counts: ",
      paste(missing_samples, collapse = ", ")
    )
  }

  counts_df <- counts_df %>%
    select(
      gene_id,
      all_of(sample_ids)
    )

  # Convert to matrix
  count_matrix <- as.matrix(
    counts_df[, -1]
  )

  rownames(count_matrix) <- counts_df$gene_id

  storage.mode(count_matrix) <- "numeric"

  # Convert to integer counts
  count_matrix <- round(
    count_matrix
  )

  storage.mode(count_matrix) <- "integer"

  # Define condition
  condition <- factor(
    dataset_metadata$disease_group,
    levels = c(
      "Control",
      "COVID-19"
    )
  )

  # Validate condition
  if (any(is.na(condition))) {
    stop(
      dataset,
      ": invalid disease_group values."
    )
  }

  # Create DESeq2 sample metadata
  col_data <- dataset_metadata %>%
    mutate(
      condition = condition
    ) %>%
    as.data.frame()

  rownames(col_data) <- col_data$sample_id

  # Ensure count matrix and metadata are aligned
  if (!identical(
    colnames(count_matrix),
    rownames(col_data)
  )) {
    stop(
      dataset,
      ": count matrix and metadata sample order do not match."
    )
  }

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = col_data,
    design = ~ condition
  )

  # Filter low-count genes
  keep <- rowSums(
    counts(dds) >= 10
  ) >= ceiling(
    0.10 * ncol(dds)
  )

  dds <- dds[
    keep,
    ]

  cat(
    "Genes after filtering:",
    nrow(dds),
    "\n"
  )

  # Run DESeq2
  dds <- DESeq(
    dds
  )

  # COVID-19 vs Control
  res <- results(
    dds,
    contrast = c(
      "condition",
      "COVID-19",
      "Control"
    ),
    alpha = 0.05
  )

  # Convert result to data frame
  result_df <- as.data.frame(
    res
  ) %>%
    rownames_to_column(
      "gene_id"
    ) %>%
    arrange(
      padj
    )

  # Add significance categories
  result_df <- result_df %>%
    mutate(
      significant = !is.na(padj) &
        padj < 0.05,
      strong_effect = !is.na(padj) &
        padj < 0.05 &
        abs(log2FoldChange) > 1
    )

  # Save complete DE results
  output_file <- file.path(
    OUTPUT_DIR,
    paste0(
      "DESeq2_",
      dataset,
      "_COVID_vs_Control.csv"
    )
  )

  write_csv(
    result_df,
    output_file
  )

  # Save significant genes
  significant_df <- result_df %>%
    filter(
      significant
    )

  significant_file <- file.path(
    OUTPUT_DIR,
    paste0(
      "DESeq2_",
      dataset,
      "_significant.csv"
    )
  )

  write_csv(
    significant_df,
    significant_file
  )

  # Save summary
  result_summary[[dataset]] <- tibble(
    dataset_id = dataset,
    n_samples = ncol(dds),
    n_covid = sum(
      condition == "COVID-19"
    ),
    n_control = sum(
      condition == "Control"
    ),
    genes_tested = nrow(result_df),
    significant_genes = sum(
      result_df$significant,
      na.rm = TRUE
    ),
    strong_effect_genes = sum(
      result_df$strong_effect,
      na.rm = TRUE
    )
  )

  cat(
    "Significant genes:",
    sum(
      result_df$significant,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Strong-effect genes:",
    sum(
      result_df$strong_effect,
      na.rm = TRUE
    ),
    "\n"
  )

  cat(
    "Saved:",
    output_file,
    "\n"
  )
}

# Combine DE summary
result_summary <- bind_rows(
  result_summary
)

write_csv(
  result_summary,
  file.path(
    OUTPUT_DIR,
    "DESeq2_summary.csv"
  )
)

# --------------------------------------------------
# GSE208342 sensitivity analysis
# Adjust for age and sex.
# --------------------------------------------------

dataset <- "GSE208342"

cat(
  "\nRunning GSE208342 sensitivity analysis.\n"
)

# Select GSE208342 metadata
sensitivity_metadata <- metadata %>%
  filter(
    dataset_id == dataset
  )

# Check required covariates
if (any(is.na(sensitivity_metadata$age))) {
  stop(
    "GSE208342 sensitivity analysis requires complete age data."
  )
}

if (any(is.na(sensitivity_metadata$sex)) ||
    any(trimws(as.character(sensitivity_metadata$sex)) == "")) {
  stop(
    "GSE208342 sensitivity analysis requires complete sex data."
  )
}

if (!is.numeric(sensitivity_metadata$age)) {
  stop(
    "GSE208342 sensitivity analysis requires numeric age data."
  )
}

# Read GSE208342 gene-level counts
count_file <- file.path(
  COUNT_DIR,
  "gene_counts_GSE208342.csv"
)

counts_df <- read_csv(
  count_file,
  show_col_types = FALSE
) %>%
  rename(
    gene_id = 1
  )

# Check sample availability in the count matrix
sample_ids <- sensitivity_metadata$sample_id

missing_samples <- setdiff(
  sample_ids,
  colnames(counts_df)
)

if (length(missing_samples) > 0) {
  stop(
    "GSE208342 sensitivity: samples missing from counts: ",
    paste(
      missing_samples,
      collapse = ", "
    )
  )
}


# Select samples in metadata order
counts_df <- counts_df %>%
  select(
    gene_id,
    all_of(sample_ids)
  )

# Convert gene counts to an integer matrix
count_matrix <- as.matrix(
  counts_df[, -1]
)

rownames(count_matrix) <- counts_df$gene_id

storage.mode(count_matrix) <- "numeric"

count_matrix <- round(
  count_matrix
)

storage.mode(count_matrix) <- "integer"

# Define sensitivity model variables
sensitivity_metadata <- sensitivity_metadata %>%
  mutate(
    # Center age to improve numerical stability
    age_centered = age - mean(age),
    condition = factor(
      disease_group,
      levels = c("Control", "COVID-19")
    ),
    sex = factor(sex)
  ) %>%
  as.data.frame()

# Set sample IDs as row names
rownames(sensitivity_metadata) <- sensitivity_metadata$sample_id

# Ensure count matrix and metadata are aligned
stopifnot(
  identical(
    colnames(count_matrix),
    rownames(sensitivity_metadata)
  )
)

# Check condition groups
if (anyNA(sensitivity_metadata$condition) ||
    any(table(sensitivity_metadata$condition) == 0)) {
  stop(
    "GSE208342 sensitivity: invalid or missing condition groups."
  )
}

# Check model matrix for rank deficiency
model_matrix <- model.matrix(
  ~ age_centered + sex + condition,
  data = sensitivity_metadata
)

if (qr(model_matrix)$rank < ncol(model_matrix)) {
  stop(
    "GSE208342 sensitivity: model matrix is not full rank."
  )
}

# Create DESeq2 object with age- and sex-adjusted design
dds_sensitivity <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sensitivity_metadata,
  design = ~ age_centered + sex + condition
)

# Filter low-count genes
keep <- rowSums(
  counts(dds_sensitivity) >= 10
) >= ceiling(
  0.10 * ncol(dds_sensitivity)
)

dds_sensitivity <- dds_sensitivity[
  keep,
]

# Run DESeq2
dds_sensitivity <- DESeq(
  dds_sensitivity
)

# Test COVID-19 versus Control
res_sensitivity <- results(
  dds_sensitivity,
  contrast = c(
    "condition",
    "COVID-19",
    "Control"
  ),
  alpha = 0.05
)

# Prepare sensitivity analysis results
sensitivity_df <- as.data.frame(
  res_sensitivity
) %>%
  rownames_to_column(
    "gene_id"
  ) %>%
  arrange(
    padj
  ) %>%
  mutate(
    significant = !is.na(padj) &
      padj < 0.05,
    strong_effect = !is.na(padj) &
      padj < 0.05 &
      abs(log2FoldChange) > 1
  )

# Save sensitivity analysis results
write_csv(
  sensitivity_df,
  file.path(
    OUTPUT_DIR,
    "DESeq2_GSE208342_COVID_vs_Control_age_sex_adjusted.csv"
  )
)

# Report significant genes
cat(
  "GSE208342 sensitivity significant genes:",
  sum(
    sensitivity_df$significant,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nDifferential expression analysis completed successfully.\n"
)

cat(
  "Output directory:",
  OUTPUT_DIR,
  "\n"
)