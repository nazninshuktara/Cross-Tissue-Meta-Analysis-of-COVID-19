#!/usr/bin/env Rscript

# 01_random_effects_meta_analysis.R
# Combine per-dataset differential expression results using random-effects meta-analysis (metafor REML).
# Calculate Per-gene pooled log2FC, SE, 95% CI, z, p, I², τ², Cochran's Q test p.
# Apply BH-FDR correction and assess cross-tissue effect-direction concordance.
# Generate meta significant plots and Table.

# Load required packages
suppressPackageStartupMessages({
  library(readr)         # Read and write CSV files
  library(dplyr)         # Data manipulation
  library(tibble)        # Structured result tables
  library(metafor)       # Random-effects meta-analysis
  library(ggplot2)       # Data visualization
  library(gt)            # Publication-ready tables
  library(RColorBrewer)  # Publication-friendly color palettes
})

# Define project paths
PROJECT_ROOT <- "."

# Input directories and files
# Directory containing complete DESeq2 result files.
INPUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "differential_expression"
)

# Output directory
# Store meta-analysis results, tables, and plots.
OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "meta_analysis",
  "random_effects_meta_analysis"
)

# Directory for storing plots, tables
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")

# Create the output directory if it does not already exist.
dir.create(
  OUTPUT_DIR, 
  recursive = TRUE, 
  showWarnings = FALSE
)

dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)


# Define analysis settings
FDR_THRESHOLD <- 0.05
MIN_DATASETS <- 2

# Expected discovery datasets
EXPECTED_DATASETS <- c(
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

# Identify complete DESeq2 result files
result_files <- list.files(
  INPUT_DIR,
  pattern = "^DESeq2_GSE[0-9]+_COVID_vs_Control\\.csv$",
  full.names = TRUE
)

file_ids <- sub(
  "^DESeq2_(GSE[0-9]+)_COVID_vs_Control\\.csv$",
  "\\1",
  basename(result_files)
)

# Check that all expected datasets are present
missing_datasets <- setdiff(EXPECTED_DATASETS, file_ids)
unexpected_datasets <- setdiff(file_ids, EXPECTED_DATASETS)

if (length(missing_datasets) > 0) {
  stop(
    "Missing DESeq2 result files for: ",
    paste(missing_datasets, collapse = ", ")
  )
}

if (anyDuplicated(file_ids)) {
  stop("Duplicate dataset result files detected.")
}

if (length(unexpected_datasets) > 0) {
  stop(
    "Unexpected DESeq2 result files for: ",
    paste(unexpected_datasets, collapse = ", ")
  )
}

# Keep only the prespecified discovery datasets in analysis order
result_files <- result_files[
  match(EXPECTED_DATASETS, file_ids)
]

message("Expected datasets found: ", length(result_files))

# Validate and import per-dataset effect estimates
required_columns <- c(
  "gene_id",
  "log2FoldChange",
  "lfcSE"
)

# Initialize lists for dataset-level results and input quality control
dataset_results <- list()
input_qc <- vector("list", length(result_files))

# Process each dataset's DESeq2 result file
for (i in seq_along(result_files)) {
  
  # Define the current file path and dataset identifier
  file <- result_files[i]
  dataset_id <- EXPECTED_DATASETS[i]
  
  # Read the DESeq2 result table
  result <- read_csv(file, show_col_types = FALSE)
  
  # Check required columns
  missing_columns <- setdiff(required_columns, names(result))
  
  # Stop if any required column is missing
  if (length(missing_columns) > 0) {
    stop(
      basename(file),
      " is missing required columns: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Check gene identifiers
  # Stop if gene IDs are missing, empty, or duplicated
  if (anyNA(result$gene_id) ||
      any(trimws(as.character(result$gene_id)) == "")) {
    stop("Missing or empty gene IDs in: ", basename(file))
  }
  
  if (anyDuplicated(result$gene_id)) {
    stop("Duplicate gene IDs detected in: ", basename(file))
  }
  
  # Check numeric effect estimates
  # Stop if log2 fold changes or standard errors are not numeric
  if (!is.numeric(result$log2FoldChange) ||
      !is.numeric(result$lfcSE)) {
    stop(
      "log2FoldChange and lfcSE must be numeric in: ",
      basename(file)
    )
  }
  
  # Retain finite effects and positive standard errors
  valid <- is.finite(result$log2FoldChange) &
    is.finite(result$lfcSE) &
    result$lfcSE > 0
  
  # Record input row counts and valid/invalid estimate counts
  input_qc[[i]] <- tibble(
    dataset_id = dataset_id,
    input_rows = nrow(result),
    valid_estimates = sum(valid),
    invalid_estimates = sum(!valid)
  )
 
  # Retain required columns and attach the dataset identifier
  # Exclude non-finite effects and non-positive standard errors 
  result <- result %>%
    transmute(
      gene_id = as.character(gene_id),
      log2FoldChange = log2FoldChange,
      lfcSE = lfcSE,
      dataset_id = dataset_id
    ) %>%
    filter(
      is.finite(log2FoldChange),
      is.finite(lfcSE),
      lfcSE > 0
    )
  
  # Store the validated estimates under the dataset identifier
  dataset_results[[dataset_id]] <- result
  
  # Report the number of valid estimates retained for this dataset
  message(
    dataset_id,
    ": ",
    nrow(result),
    " valid estimates from ",
    input_qc[[i]]$input_rows,
    " input rows"
  )
}

# Combine input quality-control summaries across all datasets
input_qc_table <- bind_rows(input_qc)

# Save the input quality-control summary as a CSV file
write_csv(
  input_qc_table,
  file.path(TABLE_DIR, "meta_analysis_input_qc.csv")
)

# Combine validated effect estimates from all discovery datasets
all_estimates <- bind_rows(dataset_results)

# Stop if no valid effect estimates remain for meta-analysis
if (nrow(all_estimates) == 0) {
  stop("No valid effect estimates remain after input validation.")
}

# Count the number of datasets contributing valid estimates per gene
gene_counts <- all_estimates %>%
  dplyr::count(gene_id, name = "n_datasets")

# Retain genes with estimates from at least MIN_DATASETS datasets
eligible_genes <- gene_counts %>%
  filter(n_datasets >= MIN_DATASETS) %>%
  pull(gene_id)

# Report the number of unique and eligible genes
message("Unique genes with valid estimates: ", nrow(gene_counts))
message("Genes eligible for meta-analysis: ", length(eligible_genes))

# Stop if no genes meet the minimum dataset requirement
if (length(eligible_genes) == 0) {
  stop("No genes have valid estimates in at least two datasets.")
}

# Fit a separate random-effects meta-analysis model for each eligible gene
meta_results <- vector("list", length(eligible_genes))
names(meta_results) <- eligible_genes

# Iterate over genes eligible for meta-analysis
for (i in seq_along(eligible_genes)) {
  
  # Select the current gene
  gene <- eligible_genes[i]
  
  # Extract and order valid estimates for the current gene
  gene_data <- all_estimates %>%
    filter(gene_id == gene) %>%
    arrange(dataset_id)
  
  # Fit REML model using log2 fold changes and their standard errors
  model <- tryCatch(
    metafor::rma(
      yi = gene_data$log2FoldChange,
      sei = gene_data$lfcSE,
      method = "REML",
      test = "z"
    ),
    error = function(e) e
  )
  
  # Record model failures without stopping the entire analysis
  if (inherits(model, "error")) {
    
    # Store missing model statistics and the failure reason
    meta_results[[i]] <- tibble(
      gene_id = gene,
      n_datasets = nrow(gene_data),
      datasets = paste(gene_data$dataset_id, collapse = ";"),
      pooled_log2FC = NA_real_,
      standard_error = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      H2 = NA_real_,
      Q = NA_real_,
      Q_p_value = NA_real_,
      model_status = paste(
        "failed:",
        conditionMessage(model)
      )
    )
    # Skip to the next gene after recording the failure
    next
  }
  
  # Extract pooled effects and heterogeneity statistics
  meta_results[[i]] <- tibble(
    gene_id = gene,
    n_datasets = nrow(gene_data),
    datasets = paste(gene_data$dataset_id, collapse = ";"),
    pooled_log2FC = as.numeric(model$b[1, 1]),
    standard_error = as.numeric(model$se[1]),
    ci_lower = as.numeric(model$ci.lb[1]),
    ci_upper = as.numeric(model$ci.ub[1]),
    z_value = as.numeric(model$zval[1]),
    p_value = as.numeric(model$pval[1]),
    tau2 = as.numeric(model$tau2),
    I2 = as.numeric(model$I2),
    H2 = as.numeric(model$H2),
    Q = as.numeric(model$QE),
    Q_p_value = as.numeric(model$QEp),
    model_status = "success"
  )
  
  if (i %% 1000 == 0 || i == length(eligible_genes)) {
    message(
      "Processed ",
      i,
      " / ",
      length(eligible_genes),
      " genes"
    )
  }
}

# Combine gene-level meta-analysis results
meta_results_table <- bind_rows(meta_results)

# Summarize initial model-fitting outcomes
model_status_summary <- meta_results_table %>%
  dplyr::count(model_status, name = "n_genes")

print(model_status_summary)

# Identify genes whose initial models failed
failed_gene_ids <- meta_results_table %>%
  filter(grepl("^failed:", model_status)) %>%
  pull(gene_id)

message(
  "Genes requiring model refitting: ",
  length(failed_gene_ids)
)

# Inspect initial model-fitting failures
failed_models <- meta_results_table %>%
  filter(grepl("^failed:", model_status)) %>%
  select(gene_id, n_datasets, datasets, model_status)

print(failed_models, n = Inf)

# Retry failed models using increased iteration limits
retry_results <- lapply(failed_gene_ids, function(gene) {
  
  # Retrieve the valid estimates for the failed gene
  gene_data <- all_estimates %>%
    filter(gene_id == gene) %>%
    arrange(dataset_id)
  
  # First retry: increase the maximum number of iterations
  model <- tryCatch(
    metafor::rma(
      yi = gene_data$log2FoldChange,
      sei = gene_data$lfcSE,
      method = "REML",
      test = "z",
      control = list(maxiter = 1000)
    ),
    error = function(e) e
  )
  
  retry_method <- "maxiter_1000"
  
  # Second retry: reduce the step adjustment if the first retry fails
  if (inherits(model, "error")) {
    model <- tryCatch(
      metafor::rma(
        yi = gene_data$log2FoldChange,
        sei = gene_data$lfcSE,
        method = "REML",
        test = "z",
        control = list(stepadj = 0.5, maxiter = 1000)
      ),
      error = function(e) e
    )
    
    retry_method <- "stepadj_0.5_maxiter_1000"
  }
  
  # Preserve the failure reason if both retries fail
  if (inherits(model, "error")) {
    return(tibble(
      gene_id = gene,
      n_datasets = nrow(gene_data),
      datasets = paste(gene_data$dataset_id, collapse = ";"),
      pooled_log2FC = NA_real_,
      standard_error = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      z_value = NA_real_,
      p_value = NA_real_,
      tau2 = NA_real_,
      I2 = NA_real_,
      H2 = NA_real_,
      Q = NA_real_,
      Q_p_value = NA_real_,
      model_status = paste(
        "failed after retries:",
        conditionMessage(model)
      ),
      retry_method = retry_method
    ))
  }
  
  # Extract estimates from the successfully refitted model
  tibble(
    gene_id = gene,
    n_datasets = nrow(gene_data),
    datasets = paste(gene_data$dataset_id, collapse = ";"),
    pooled_log2FC = as.numeric(model$b[1, 1]),
    standard_error = as.numeric(model$se[1]),
    ci_lower = as.numeric(model$ci.lb[1]),
    ci_upper = as.numeric(model$ci.ub[1]),
    z_value = as.numeric(model$zval[1]),
    p_value = as.numeric(model$pval[1]),
    tau2 = as.numeric(model$tau2),
    I2 = as.numeric(model$I2),
    H2 = as.numeric(model$H2),
    Q = as.numeric(model$QE),
    Q_p_value = as.numeric(model$QEp),
    model_status = "success_after_retry",
    retry_method = retry_method
  )
})

# Combine the retry results
retry_results <- bind_rows(retry_results)

# Review the retry outcomes
print(retry_results %>% 
        dplyr::count(model_status, retry_method), n = Inf)


# Confirm that all failed genes were successfully refitted
stopifnot(
  length(failed_gene_ids) == nrow(retry_results),
  !anyDuplicated(retry_results$gene_id),
  setequal(failed_gene_ids, retry_results$gene_id),
  all(retry_results$model_status == "success_after_retry"),
  all(is.finite(retry_results$p_value))
)

if (any(retry_results$model_status != "success_after_retry") ||
    any(!is.finite(retry_results$p_value))) {
  stop("Some failed models could not be successfully refitted.")
}

# Remove the original failed rows and insert the successful retry results
meta_results_table <- meta_results_table %>%
  filter(!gene_id %in% failed_gene_ids) %>%
  bind_rows(retry_results %>% select(-retry_method))

# Verify the final model status
print(meta_results_table %>% 
        dplyr::count(model_status))

# Confirm that the expected number of genes is retained
stopifnot(
  nrow(meta_results_table) == length(eligible_genes),
  !anyDuplicated(meta_results_table$gene_id),
  all(meta_results_table$model_status == "success_after_retry" |
        meta_results_table$model_status == "success")
)

# Validate meta-analysis p-values before multiple-testing correction
valid_p <- !is.na(meta_results_table$p_value) &
  is.finite(meta_results_table$p_value) &
  meta_results_table$p_value >= 0 &
  meta_results_table$p_value <= 1

# Stop if any gene has an invalid p-value
if (!all(valid_p)) {
  stop(
    "Invalid meta-analysis p-values detected for ",
    sum(!valid_p),
    " genes."
  )
}

# Apply Benjamini-Hochberg false discovery rate correction
meta_results_table <- meta_results_table %>%
  mutate(
    padj = p.adjust(p_value, method = "BH"),
    meta_significant = padj < FDR_THRESHOLD
  )

# Save retry details 
write_csv(
  retry_results,
  file.path(TABLE_DIR, "meta_analysis_retry_results.csv")
)

# Save the complete meta-analysis results
readr::write_csv(
  meta_results_table,
  file.path(TABLE_DIR, "random_effects_meta_analysis.csv")
)

# Save meta-significant genes
meta_significant_genes <- meta_results_table %>%
  filter(meta_significant) %>%
  arrange(padj)

readr::write_csv(
  meta_significant_genes,
  file.path(TABLE_DIR, "meta_significant_genes.csv")
)

# Assess effect-direction concordance across datasets
# Concordance is the fraction of available estimates matching the pooled
# effect direction. It is descriptive and does not replace heterogeneity tests.
concordance_results <- all_estimates %>%
  filter(
    is.finite(log2FoldChange)
  ) %>%
  inner_join(
    meta_results_table %>%
      select(gene_id, 
             pooled_log2FC, 
             significant = meta_significant),
    by = "gene_id"
  ) %>%
  filter(
    is.finite(pooled_log2FC),
    pooled_log2FC != 0
  ) %>%
  group_by(gene_id, pooled_log2FC, significant) %>%
  summarise(
    n_datasets = n(),
    n_concordant = sum(
      sign(log2FoldChange) ==
        sign(pooled_log2FC[1])
    ),
    direction_concordance = n_concordant / n_datasets,
    .groups = "drop"
  ) %>%
  arrange(
    desc(significant), 
    desc(direction_concordance)
  )

# Save cross-dataset effect-direction concordance results
readr::write_csv(
  concordance_results,
  file.path(TABLE_DIR, "cross_dataset_effect_direction_concordance.csv"),
  na = ""
)

# Report the number of genes assessed
message(
  "Effect-direction concordance results saved for ",
  nrow(concordance_results),
  " genes."
)

# Create a volcano plot of pooled effects
volcano_data <- meta_results_table %>%
  filter(
    is.finite(pooled_log2FC),
    is.finite(padj)
  ) %>%
  mutate(
    neg_log10_fdr = -log10(pmax(padj, .Machine$double.xmin)),
    regulation = case_when(
      meta_significant & pooled_log2FC > 0 ~ "Upregulated",
      meta_significant & pooled_log2FC < 0 ~ "Downregulated",
      TRUE ~ "Not significant"
    )
  )

# Generate the volcano plot
volcano_plot <- ggplot(
  volcano_data,
  aes(
    x = pooled_log2FC,
    y = neg_log10_fdr,
    color = regulation
  )
) +
  geom_point(alpha = 0.65, size = 1.2) +
  scale_color_manual(
    values = c(
      "Upregulated" = "#B2182B",
      "Downregulated" = "#2166AC",
      "Not significant" = "grey70"
    ),
    breaks = c(
      "Upregulated",
      "Downregulated",
      "Not significant"
    )
  ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey35"
  ) +
  geom_hline(
    yintercept = -log10(FDR_THRESHOLD),
    linetype = "dashed",
    color = "grey35"
  ) +
  labs(
    title = "Cross-tissue random-effects meta-analysis",
    x = expression("Pooled log"[2] * " fold change"),
    y = expression(-log[10]("BH-adjusted p-value")),
    color = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold")
)

# Save the volcano plot
ggsave(
  file.path(PLOT_DIR, "meta_analysis_volcano.png"),
  plot = volcano_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# Plot direction concordance among meta-significant genes
concordance_plot_data <- concordance_results %>%
  filter(
    significant,
    is.finite(direction_concordance),
    direction_concordance >= 0,
    direction_concordance <= 1
  )

if (nrow(concordance_plot_data) > 0) {
  
  concordance_plot <- ggplot(
    concordance_plot_data,
    aes(x = direction_concordance)
  ) +
    geom_histogram(
      bins = 20,
      fill = brewer.pal(3, "Set2")[2],
      color = "white"
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, 0.2)
    ) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(
      title = "Effect-direction concordance of meta-significant genes",
      x = "Fraction of datasets agreeing with pooled effect direction",
      y = "Number of genes"
    ) +
    theme_bw(base_size = 12)

 # Save the concordance plot  
  ggsave(
    file.path(PLOT_DIR, "cross_tissue_concordance.png"),
    concordance_plot,
    width = 8,
    height = 6,
    dpi = 300
  )
}

# Create a table for top 20 meta significant genes
# Map gene IDs to gene symbols
gene_symbol_map <- readr::read_csv(
  "inputs/reference/gencode.v49_gene_id_symbol_map.csv",
  show_col_types = FALSE
) %>%
  mutate(gene_id_clean = sub("\\.[0-9]+$", "", gene_id)) %>%
  select(gene_id_clean, gene_symbol) %>%
  distinct(gene_id_clean, .keep_all = TRUE)

# Select the top 20 meta-significant genes by BH-adjusted p-value.
top20_genes <- meta_results_table %>%
  filter(
    meta_significant,
    is.finite(padj),
    is.finite(pooled_log2FC)
  ) %>%
  arrange(padj, desc(abs(pooled_log2FC)), gene_id) %>%
  slice_head(n = 20) %>%
  mutate(
    Rank = row_number(),
    gene_id_clean = sub("\\.[0-9]+$", "", as.character(gene_id))
  ) %>%
  left_join(gene_symbol_map, by = "gene_id_clean") %>%
  mutate(
    gene_label = if_else(
      !is.na(gene_symbol) & gene_symbol != "",
      gene_symbol,
      as.character(gene_id)
    )
  )

# Stop if no eligible genes are available.
if (nrow(top20_genes) == 0) {
  stop("No eligible meta-significant genes were found.")
}

# Count significant genes in each discovery dataset.
# Dataset-level significance is defined as DESeq2 padj < FDR_THRESHOLD.
dataset_significance <- bind_rows(lapply(seq_along(result_files), function(i) {
  readr::read_csv(result_files[i], show_col_types = FALSE) %>%
    transmute(
      gene_id_clean = sub("\\.[0-9]+$", "", as.character(gene_id)),
      significant = !is.na(padj) & is.finite(padj) & 
        padj < FDR_THRESHOLD
    ) %>%
    distinct(gene_id_clean, .keep_all = TRUE) %>%
    mutate(dataset_id = EXPECTED_DATASETS[i])
})) %>%
  group_by(gene_id_clean) %>%
  summarise(
    n_datasets_significant = sum(significant),
    .groups = "drop"
  )

# Assemble Table 2 using the existing ranking.
table_2_top20 <- top20_genes %>%
  left_join(
    concordance_results %>%
      select(gene_id, direction_concordance) %>%
      distinct(gene_id, .keep_all = TRUE),
    by = "gene_id"
  ) %>%
  left_join(dataset_significance, by = "gene_id_clean") %>%
  mutate(n_datasets_significant = coalesce(n_datasets_significant, 0L)) %>%
  arrange(Rank) %>%
  transmute(
    Rank,
    `Gene symbol` = gene_label,
    `Pooled log2FC` = pooled_log2FC,
    `Pooled padj` = padj,
    `95% CI lower` = ci_lower,
    `95% CI upper` = ci_upper,
    `I2 (%)` = I2,
    `Datasets tested` = n_datasets,
    `Datasets significant` = n_datasets_significant,
    `Concordance fraction` = direction_concordance
  )

stopifnot(
  nrow(table_2_top20) == nrow(top20_genes),
  identical(table_2_top20$Rank, seq_len(nrow(table_2_top20)))
)

# Format and export the table.
top20_significant_gene_table <- table_2_top20 %>%
  gt() %>%
  tab_header(
    title = md("**Table 2. Top 20 Meta-significant Genes**"),
    subtitle = paste(
      "Top 20 genes from the random-effects meta-analysis,",
      "ranked by BH-adjusted p-value."
    )
  ) %>%
  cols_label(
    Rank = "Rank",
    `Gene symbol` = "Gene symbol",
    `Pooled log2FC` = "Pooled log2FC",
    `Pooled padj` = "Pooled padj",
    `95% CI lower` = "95% CI lower",
    `95% CI upper` = "95% CI upper",
    `I2 (%)` = "I² (%)",
    `Datasets tested` = "Datasets tested",
    `Datasets significant` = "Datasets significant",
    `Concordance fraction` = "Concordance fraction"
  ) %>%
  # Format numeric columns
  fmt_number(
    columns = c(
      `Pooled log2FC`,
      `95% CI lower`,
      `95% CI upper`,
      `I2 (%)`,
      `Concordance fraction`
    ),
    decimals = 2
  ) %>%
  fmt_scientific(
    columns = `Pooled padj`,
    decimals = 2
  ) %>%
  
  # Style the table header
  tab_style(
    style = list(
      cell_fill(color = "#24476B"),
      cell_text(
        color = "white",
        weight = "bold"
      )
    ),
    locations = cells_column_labels(columns = everything())
  ) %>%
  
  # Add vertical separators between columns
  tab_style(
    style = cell_borders(
      sides = "right",
      color = "#B7C3D0",
      weight = px(1)
    ),
    locations = cells_body(columns = everything())
  ) %>%
  
  # Add horizontal separators between rows
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "#B7C3D0",
      weight = px(1)
    ),
    locations = cells_body(columns = everything())
  ) %>%
  
  # Add vertical separators to the header
  tab_style(
    style = cell_borders(
      sides = "right",
      color = "#FFFFFF",
      weight = px(1)
    ),
    locations = cells_column_labels(columns = everything())
  ) %>%
  
  # Add a bottom border to the header
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "#24476B",
      weight = px(2)
    ),
    locations = cells_column_labels(columns = everything())
  ) %>%
  
  # Add an explanatory footnote
  tab_source_note(
    source_note = paste0(
      "Genes are ranked by BH-adjusted p-value from the ",
      "random-effects meta-analysis. Pooled log2FC represents ",
      "the estimated overall effect. The 95% CI shows the ",
      "confidence interval of the pooled effect. I² indicates ",
      "between-dataset heterogeneity. Dataset-level significance ",
      "is defined as DESeq2 BH-adjusted p-value < ",
      FDR_THRESHOLD,
      ". Concordance fraction indicates agreement with the ",
      "pooled-effect direction."
    )
  ) %>%
  
  # Italicize the footnote
  tab_style(
    style = cell_text(
      style = "italic",
      size = px(11)
    ),
    locations = cells_source_notes()
  ) %>%
  
  # Set table-wide formatting
  tab_options(
    table.font.size = px(11),
    column_labels.font.weight = "bold",
    heading.title.font.size = px(16),
    heading.subtitle.font.size = px(11),
    data_row.padding = px(5),
    table.border.top.color = "#24476B",
    table.border.top.width = px(2),
    table.border.bottom.color = "#24476B",
    table.border.bottom.width = px(2)
  )

# Save the table
gt::gtsave(
  top20_significant_gene_table,
  file.path(TABLE_DIR, "top20_meta_significant_genes.html")
)

# Save the random-effects meta-analysis summary
analysis_summary <- tibble(
  metric = c(
    "datasets_expected",
    "datasets_analyzed",
    "unique_genes_with_valid_estimates",
    "genes_eligible_for_meta_analysis",
    "successful_models",
    "failed_models",
    "genes_with_valid_p_values",
    "meta_significant_genes",
    "non_significant_genes",
    "minimum_datasets_per_gene",
    "median_datasets_per_gene",
    "maximum_datasets_per_gene",
    "FDR_threshold",
    "meta_analysis_method",
    "multiple_testing_method"
  ),
  value = as.character(c(
    length(EXPECTED_DATASETS),
    length(result_files),
    nrow(gene_counts),
    length(eligible_genes),
    sum(meta_results_table$model_status %in%
          c("success", "success_after_retry")),
    sum(!meta_results_table$model_status %in%
          c("success", "success_after_retry")),
    sum(
      is.finite(meta_results_table$p_value) &
        meta_results_table$p_value >= 0 &
        meta_results_table$p_value <= 1
    ),
    sum(meta_results_table$meta_significant, na.rm = TRUE),
    sum(!meta_results_table$meta_significant, na.rm = TRUE),
    MIN_DATASETS,
    median(meta_results_table$n_datasets),
    max(meta_results_table$n_datasets),
    FDR_THRESHOLD,
    "Random-effects model (REML, z-test)",
    "Benjamini-Hochberg FDR"
  ))
)

write_csv(
  analysis_summary,
  file.path(TABLE_DIR, "meta_analysis_summary.csv")
)

# Report analysis completion and output locations
message("\nRandom-effects meta-analysis completed.")
message("Summary: ", file.path(TABLE_DIR, "meta_analysis_summary.csv"))
message("Full results: ", file.path(TABLE_DIR, "random_effects_meta_analysis.csv"))
message("Significant genes: ", file.path(TABLE_DIR, "meta_significant_genes.csv"))
message("Retry results: ", file.path(TABLE_DIR, "meta_analysis_retry_results.csv"))
message("Input QC: ", file.path(TABLE_DIR, "meta_analysis_input_qc.csv"))
message("Plots directory: ", PLOT_DIR)
message("Meta-significant genes (BH-FDR < ", FDR_THRESHOLD, "): ",sum(meta_results_table$meta_significant %in% TRUE))
message("Non-significant genes (BH-FDR >= ", FDR_THRESHOLD, "): ",sum(meta_results_table$meta_significant %in% FALSE))
