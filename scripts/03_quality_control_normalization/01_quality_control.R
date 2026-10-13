#!/usr/bin/env Rscript

# 01_quality_control.R
# Perform library-size QC, gene filtering, VST, global PCA, and outlier detection.

# Load required packages
suppressPackageStartupMessages({
  library(DESeq2)   # RNA-seq count processing and VST transformation
  library(dplyr)    # Data manipulation
  library(readr)    # Read and write CSV files
  library(tibble)   # Tidy data frame utilities
  library(ggplot2)  # Data Visualization
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
# Directory for storing gene-level counts QC results.
OUTPUT_DIR <- file.path(
  PROJECT_ROOT, 
  "outputs", 
  "gene_counts_qc"
)

# Directory for storing plots, tables, reports
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
REPORT_DIR <- file.path(OUTPUT_DIR, "reports")

# Create the output directory if it does not already exist.
dir.create(
  OUTPUT_DIR, 
  recursive = TRUE, 
  showWarnings = FALSE
)

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(REPORT_DIR, recursive = TRUE, showWarnings = FALSE)

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

# Define dataset-specific low-count filtering thresholds
MIN_COUNT <- 10 
MIN_SAMPLE_FRACTION <- 0.10

# Read curated sample-level metadata
metadata <- readr::read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

required_metadata_columns <- c(
  "dataset_id", 
  "sample_id", 
  "tissue_type", 
  "disease_status", 
  "disease_group", 
  "disease_severity", 
  "inclusion_status" 
)

if (!all(required_metadata_columns %in% colnames(metadata))) { 
  stop("Required metadata columns are missing.") 
}

# Keep primary-analysis samples from discovery datasets
metadata <- metadata %>%
  dplyr::filter(inclusion_status == "included") %>%
  dplyr::select(dplyr::all_of(required_metadata_columns))

# Validate curated metadata
if (nrow(metadata) == 0) {
  stop("No primary samples found in curated metadata.")
}

if (anyNA(metadata$dataset_id) || anyNA(metadata$sample_id)) { 
  stop("Missing dataset_id or sample_id in curated metadata.")
}

if (anyDuplicated(metadata[c("dataset_id", "sample_id")])) {
  stop("Duplicate dataset_id and sample_id combinations found.")
}

if (!all(DATASETS %in% unique(metadata$dataset_id))) { 
  stop("One or more discovery datasets are missing from metadata.") 
}


# Read and validate count and gene-length files
count_matrices <- list()
gene_length_matrices <- list()
input_summary <- list()

for (dataset in DATASETS) {
  
  count_file <- file.path(
    COUNT_DIR,
    paste0("gene_counts_", dataset, ".csv")
  )
  
  length_file <- file.path(
    COUNT_DIR,
    paste0("gene_lengths_", dataset, ".csv")
  )
  
  if (!file.exists(count_file)) {
    stop("Missing gene-count file: ", count_file)
  }
  
  if (!file.exists(length_file)) {
    stop("Missing gene-length file: ", length_file)
  }
  
  # Read gene counts
  counts <- readr::read_csv(
    count_file,
    show_col_types = FALSE
  )
  
  if (!"gene_id" %in% names(counts)) {
    stop(dataset, ": gene-count file has no gene_id column.")
  }
  
  if (ncol(counts) < 2 ||
      anyNA(counts$gene_id) ||
      anyDuplicated(counts$gene_id) ||
      anyDuplicated(names(counts)[-1])) {
    stop(dataset, ": invalid gene-count matrix structure.")
  }
  
  # Ensure gene_id is the first column
  if (names(counts)[1] != "gene_id") {
    stop(dataset, ": gene_id must be the first column.")
  }
  
  # Validate numeric, finite, non-negative integer counts
  count_values <- as.matrix(counts[, -1, drop = FALSE])
  
  suppressWarnings(
    storage.mode(count_values) <- "numeric"
  )
  
  if (anyNA(count_values) ||
      any(!is.finite(count_values)) ||
      any(count_values < 0) ||
      any(abs(count_values - round(count_values)) > 1e-8)) {
    stop(dataset, ": invalid gene counts detected.")
  }
  
  # Read gene lengths
  gene_lengths <- readr::read_csv(
    length_file,
    show_col_types = FALSE
  )
  
  if (ncol(gene_lengths) < 2) {
    stop(dataset, ": gene-length file has no length column.")
  }
  
  # Standardize the first column name
  names(gene_lengths)[1] <- "gene_id"
  
  if (anyNA(gene_lengths$gene_id) ||
      anyDuplicated(gene_lengths$gene_id) ||
      !setequal(counts$gene_id, gene_lengths$gene_id)) {
    stop(dataset, ": gene IDs differ between count and length files.")
  }
  
  # Validate gene-length values
  length_values <- as.matrix(
    gene_lengths[, -1, drop = FALSE]
  )
  
  suppressWarnings(
    storage.mode(length_values) <- "numeric"
  )
  
  if (anyNA(length_values) ||
      any(!is.finite(length_values)) ||
      any(length_values <= 0)) {
    stop(dataset, ": invalid gene-length values detected.")
  }
  
  # Store validated inputs
  count_matrices[[dataset]] <- counts
  gene_length_matrices[[dataset]] <- gene_lengths
  
  input_summary[[dataset]] <- tibble::tibble(
    dataset_id = dataset,
    genes = nrow(counts),
    count_matrix_samples = ncol(counts) - 1,
    gene_length_genes = nrow(gene_lengths),
    gene_length_columns = ncol(gene_lengths) - 1
  )
  
  cat(dataset, ": input validation passed\n")
}

# Validate sample IDs against primary metadata
for (dataset in DATASETS) {
  
  matrix_samples <- names(count_matrices[[dataset]])[-1]
  
  metadata_samples <- metadata %>%
    dplyr::filter(.data$dataset_id == dataset) %>%
    dplyr::pull(.data$sample_id)
  
  missing_from_counts <- setdiff(
    metadata_samples,
    matrix_samples
  )
  
  missing_from_metadata <- setdiff(
    matrix_samples,
    metadata_samples
  )
  
  if (length(missing_from_counts) > 0 ||
      length(missing_from_metadata) > 0) {
    stop(
      dataset, ": sample IDs do not match.\n",
      "Missing from counts: ",
      paste(missing_from_counts, collapse = ", "),
      "\nMissing from metadata: ",
      paste(missing_from_metadata, collapse = ", ")
    )
  }
  
  cat(
    dataset, ": sample IDs match;",
    length(metadata_samples),
    "biological samples\n"
  )
}

# Calculate library sizes and filter low-count genes
library_sizes <- list()
filtering_summary <- list()
filtered_matrices <- list()

for (dataset in DATASETS) {
  
  # Retrieve count matrix and dataset metadata
  counts_df <- count_matrices[[dataset]]
  
  dataset_metadata <- metadata %>%
    dplyr::filter(dataset_id == dataset)
  
  metadata_samples <- dataset_metadata$sample_id
  
  # Align count matrix columns with metadata sample order
  counts_df <- counts_df %>%
    dplyr::select(
      gene_id, 
      dplyr::all_of(metadata_samples)
    )
  
  count_matrix <- as.matrix(
    counts_df[, -1, drop = FALSE]
  )
  
  suppressWarnings(
    storage.mode(count_matrix) <- "numeric"
  )
  
  rownames(count_matrix) <- counts_df$gene_id
  
  # Validate count matrix
  if (anyNA(count_matrix) ||
      any(!is.finite(count_matrix)) ||
      any(count_matrix < 0) ||
      any(abs(count_matrix - round(count_matrix)) > 1e-8)) {
    stop(dataset, ": invalid count values detected.")
  }
  
  if (anyDuplicated(rownames(count_matrix))) {
    stop(dataset, ": duplicated gene IDs detected.")
  }
  
  if (any(colSums(count_matrix) <= 0)) {
    stop(dataset, ": zero-library-size sample detected.")
  }
  
  # Calculate library sizes
  library_sizes[[dataset]] <- tibble::tibble(
    dataset_id = dataset,
    sample_id = colnames(count_matrix),
    library_size = colSums(count_matrix),
    detected_genes = colSums(count_matrix > 0),
    genes_ge_10 = colSums(count_matrix >= MIN_COUNT)
  ) %>%
    dplyr::left_join(
      dataset_metadata %>%
        dplyr::select(sample_id, 
                       tissue_type, 
                       disease_group
                      ),
      by = "sample_id"
    )
  
  # Filter low-count genes within each dataset
  n_samples <- ncol(count_matrix)
  
  min_samples <- max(
    1L, 
    ceiling(MIN_SAMPLE_FRACTION * n_samples)
  )
  
  keep_genes <- rowSums(count_matrix >= MIN_COUNT) >= min_samples
  
  filtered_matrix <- count_matrix[keep_genes, , drop = FALSE]
  
  if (nrow(filtered_matrix) == 0) {
    stop(dataset, ": no genes remain after filtering.")
  }
  
  # Validate filtered matrix
  if (anyDuplicated(rownames(filtered_matrix))) {
    stop(dataset, ": duplicated gene IDs after filtering.")
  }
  
  if (anyNA(filtered_matrix) ||
      any(!is.finite(filtered_matrix)) ||
      any(filtered_matrix < 0)) {
    stop(dataset, ": invalid values in filtered matrix.")
  }
  
  filtered_matrices[[dataset]] <- filtered_matrix
  
  # Record gene-filtering summary
  filtering_summary[[dataset]] <- tibble::tibble(
    dataset_id = dataset,
    n_samples = n_samples,
    min_count = MIN_COUNT,
    min_sample_fraction = MIN_SAMPLE_FRACTION,
    min_samples_required = min_samples,
    genes_before_filtering = nrow(count_matrix),
    genes_after_filtering = nrow(filtered_matrix),
    genes_removed = nrow(count_matrix) - nrow(filtered_matrix)
  )
}

# Combine QC summaries
library_size_summary <- dplyr::bind_rows(library_sizes)
filtering_summary_df <- dplyr::bind_rows(filtering_summary)
input_summary_df <- dplyr::bind_rows(input_summary)

# Save QC summaries
readr::write_csv(
  library_size_summary,
  file.path(TABLE_DIR, "library_size_summary.csv")
)

readr::write_csv(
  filtering_summary_df,
  file.path(TABLE_DIR, "gene_filtering_summary.csv")
)

readr::write_csv(
  input_summary_df,
  file.path(TABLE_DIR, "input_validation_summary.csv")
)

print(filtering_summary_df, n = Inf)

# Identify genes retained across all datasets
common_genes <- Reduce(
  intersect,
  lapply(filtered_matrices, rownames)
)

if (length(common_genes) < 1000) {
  stop("Fewer than 1,000 common genes remain after filtering.")
}

cat("Common filtered genes:", length(common_genes), "\n")

# Build combined count matrix for exploratory PCA
combined_counts <- do.call(
  cbind,
  lapply(
    filtered_matrices,
    function(x) x[common_genes, , drop = FALSE]
  )
)

if (anyDuplicated(colnames(combined_counts))) {
  stop("Duplicate sample IDs in combined count matrix.")
}

# Align metadata with count matrix
metadata_index <- match(
  colnames(combined_counts),
  metadata$sample_id
)

if (anyNA(metadata_index)) {
  stop("One or more count-matrix samples are missing from metadata.")
}

combined_metadata <- metadata[metadata_index, , drop = FALSE]

# Verify sample order
if (!identical(
  combined_metadata$sample_id,
  colnames(combined_counts)
)) {
  stop("Metadata and count matrix sample order does not match.")
}

cat("Metadata aligned with count matrix successfully.\n")

# Create DESeq2 object for exploratory VST
coldata <- as.data.frame(combined_metadata)
rownames(coldata) <- coldata$sample_id

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = round(combined_counts),
  colData = coldata,
  design = ~ dataset_id
)

# Apply variance-stabilizing transformation
vsd <- DESeq2::vst(
  dds,
  blind = TRUE
)
vst_matrix <- SummarizedExperiment::assay(vsd)

# Save the global VST matrix
readr::write_csv(
  as.data.frame(vst_matrix) %>%
    tibble::rownames_to_column("gene_id"),
  file.path(TABLE_DIR, "global_vst_matrix.csv")
)

# Perform global PCA
pca <- prcomp(t(vst_matrix), center = TRUE, scale. = FALSE)

variance_explained <- pca$sdev^2 / sum(pca$sdev^2)

pca_variance <- tibble::tibble(
  PC = paste0("PC", seq_along(variance_explained)),
  variance_explained = variance_explained,
  variance_percent = variance_explained * 100
)

readr::write_csv(
  pca_variance,
  file.path(TABLE_DIR, "PCA_variance.csv")
)

# Attach metadata to PCA coordinates
pca_coordinates <- as.data.frame(pca$x) %>%
  tibble::rownames_to_column("sample_id") %>%
  dplyr::left_join(combined_metadata, by = "sample_id")

if (anyNA(pca_coordinates$dataset_id)) {
  stop("Missing metadata for one or more PCA samples.")
}

# Flag extreme PC1 values
pc1_sd <- sd(pca_coordinates$PC1)

if (is.finite(pc1_sd) && pc1_sd > 0) {
  pca_coordinates <- pca_coordinates %>%
    dplyr::mutate(
      PC1_z = (PC1 - mean(PC1)) / pc1_sd,
      outlier_PC1 = abs(PC1_z) > 3
    )
} else {
  pca_coordinates <- pca_coordinates %>%
    dplyr::mutate(
      PC1_z = 0,
      outlier_PC1 = FALSE
    )
}

# Save PCA coordinates and outlier table
readr::write_csv(
  pca_coordinates,
  file.path(TABLE_DIR, "PCA_coordinates.csv")
)

readr::write_csv(
  pca_coordinates %>%
    dplyr::select(
      sample_id, dataset_id, tissue_type,
      disease_status, disease_group,
      PC1, PC1_z, outlier_PC1
    ),
  file.path(TABLE_DIR, "PCA_outliers.csv")
)

# Quantify dataset-associated variation in PC1
dataset_model <- stats::lm(PC1 ~ dataset_id, data = pca_coordinates)
dataset_R2 <- summary(dataset_model)$r.squared

readr::write_csv(
  tibble::tibble(
    model = "PC1 ~ dataset_id",
    R_squared = dataset_R2
  ),
  file.path(TABLE_DIR, "dataset_associated_PC1_variation.csv")
)

# Prepare library-size summary for plotting
library_size_summary <- dplyr::bind_rows(library_sizes)

# Plot library sizes by dataset
p_library_size <- ggplot(
  library_size_summary, 
  aes( 
    x = dataset_id, 
    y = library_size, 
    fill = dataset_id 
  ) 
) +
  geom_boxplot( 
    width = 0.7, 
    alpha = 0.85, 
    outlier.shape = 21, 
    outlier.size = 1.8 
  ) +
  scale_y_log10() + 
  scale_fill_brewer(palette = "Set3") + 
  labs( 
    title = "Library Size Distribution", 
    x = "Dataset", 
    y = "Library size (log10 scale)" 
  ) +
  theme_classic(base_size = 12) + 
  theme( 
    plot.title = element_text( 
      face = "bold", hjust = 0.5 
    ), 
    axis.text.x = element_text( 
      angle = 45, 
      hjust = 1 ), 
    legend.position = "none" 
  )

ggsave( 
  file.path(PLOT_DIR, "library_sizes.pdf"), 
  plot = p_library_size, 
  width = 10, 
  height = 6 
)
  
# Plot PCA by dataset
p_pca_dataset <- ggplot(
  pca_coordinates, 
  aes( 
    x = PC1, 
    y = PC2, 
    color = dataset_id
    
  )
) + 
  geom_point( 
    size = 3, 
    alpha = 0.85 
  ) +
  scale_color_brewer(palette = "Set1") +
  labs(
    title = "Global PCA by Dataset",
    x = paste0(
      "PC1 (",
      round(pca_variance$variance_percent[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance$variance_percent[2], 2),
      "%)"
    ),
    color = "Dataset"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.position = "right"
  )

ggsave(
  file.path(PLOT_DIR, "PCA_by_dataset.pdf"),
  plot = p_pca_dataset,
  width = 10,
  height = 7
)

# Plot PCA by disease group
p_pca_disease <- ggplot(
  pca_coordinates,
  aes( 
    x = PC1, 
    y = PC2, 
    color = disease_group 
  )
) + 
  geom_point( 
    size = 3, 
    alpha = 0.85 
  ) +
  scale_color_brewer(palette = "Dark2") +
  labs( 
    title = "Global PCA by Disease Group",
    x = paste0(
      "PC1 (",
      round(pca_variance$variance_percent[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(pca_variance$variance_percent[2], 2),
      "%)"
    ), 
    color = "Disease group" 
  ) +
  theme_classic(base_size = 12) + 
  theme( plot.title = element_text( 
    face = "bold", 
    hjust = 0.5 
    ), 
    legend.position = "right" 
  )
  
ggsave( 
  file.path(PLOT_DIR, "PCA_by_disease_group.pdf"), 
  plot = p_pca_disease, 
  width = 10, 
  height = 7 
)
  
# Save high-resolution PNG versions 
ggsave(  
  file.path(PLOT_DIR, "library_sizes.png"), 
  plot = p_library_size, 
  width = 10, 
  height = 6, 
  dpi = 300 
)

ggsave( 
  file.path(PLOT_DIR, "PCA_by_dataset.png"), 
  plot = p_pca_dataset, 
  width = 10, 
  height = 7, 
  dpi = 300 
)

ggsave( 
  file.path(PLOT_DIR, "PCA_by_disease_group.png"), 
  plot = p_pca_disease, 
  width = 10, 
  height = 7, 
  dpi = 300 
)

# Print QC summary
cat("\nGene-level count QC completed successfully.\n")
cat("Datasets:", length(DATASETS), "\n")
cat("Primary biological samples:", nrow(combined_metadata), "\n")
cat("Common filtered genes:", length(common_genes), "\n")
cat("PC1 variance explained:", round(pca_variance$variance_percent[1], 2), "%\n")
cat("PC2 variance explained:", round(pca_variance$variance_percent[2], 2), "%\n")
cat("PC1 samples flagged:", sum(pca_coordinates$outlier_PC1), "\n")
cat("Dataset-associated PC1 R-squared:", round(dataset_R2, 4), "\n")
