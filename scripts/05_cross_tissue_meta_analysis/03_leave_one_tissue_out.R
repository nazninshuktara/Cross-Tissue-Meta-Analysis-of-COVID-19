#!/usr/bin/env Rscript

# 03_leave_one_tissue_out.R
# Assess core-gene robustness by excluding one tissue group at a time.

suppressPackageStartupMessages({
  library(readr)      # Import data
  library(dplyr)      # Data wrangling
  library(tibble)     # Tibble utilities
  library(metafor)    # Meta-analysis
  library(ggplot2)    # Plotting
})

# Define project paths
PROJECT_ROOT <- "."

# Input directories and files
# Directory containing DESeq2 results
DE_DIR <- file.path(
  PROJECT_ROOT, 
  "outputs", 
  "differential_expression"
)

# Core/shared gene classification from meta-analysis
META_FILE <- file.path(
  PROJECT_ROOT,
  "outputs", 
  "meta_analysis",
  "core_tissue_specific_decomposition",
  "tables",
  "gene_classification.csv"
)

# Output directories
OUTPUT_DIR <- file.path(
  PROJECT_ROOT, 
  "outputs", 
  "meta_analysis",
  "leave_one_tissue_out"
)

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

# Analysis settings
FDR_THRESHOLD <- 0.05
MIN_DATASETS <- 2

# Discovery datasets
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

# Prespecified grouping: six tissue groups across nine datasets
TISSUE_GROUPS <- list(
  Nasopharyngeal = c("GSE208342"),
  Lung = c("GSE182917", "GSE208076"),
  Blood = c("GSE211979", "GSE273149", "GSE243217"),
  Kidney = c("GSE202182"),
  Heart = c("GSE169241"),
  Gut_organoids = c("GSE149312")
)

# Validate tissue-group assignment
assigned_datasets <- unlist(TISSUE_GROUPS, use.names = FALSE)

# Load gene classification
gene_classes <- read_csv(
  META_FILE,
  show_col_types = FALSE
)

# Required columns
required_core_columns <- c(
  "gene_id",
  "gene_symbol",
  "pooled_log2FC"
)

# Validate classification file
if (!all(required_core_columns %in% names(gene_classes))) {
  stop("Gene classification file is missing required columns.")
}

if (anyDuplicated(gene_classes$gene_id)) {
  stop("Duplicate gene IDs detected.")
}

# Retain core genes only
core_genes <- gene_classes %>%
  filter(component == "core") %>%
  transmute(
    gene_id,
    gene_symbol,
    original_log2FC = pooled_log2FC,
    original_direction = case_when(
      pooled_log2FC > 0 ~ "Positive",
      pooled_log2FC < 0 ~ "Negative",
      TRUE ~ "Zero_or_missing"
    )
  )

# Validate pooled effect sizes
if (any(!is.finite(core_genes$original_log2FC))) {
  stop("Core genes contain missing or non-finite pooled log2FC values.")
}

message("Core genes loaded: ", nrow(core_genes))

# Locate the expected DESeq2 result files
result_files <- list.files(
  DE_DIR,
  pattern = "^DESeq2_GSE[0-9]+_COVID_vs_Control\\.csv$",
  full.names = TRUE
)

# Extract dataset IDs
file_ids <- sub(
  "^DESeq2_(GSE[0-9]+)_COVID_vs_Control\\.csv$",
  "\\1",
  basename(result_files)
)

# Check duplicates
if (anyDuplicated(file_ids)) {
  stop("Duplicate DESeq2 result files detected.")
}

# Check missing datasets
missing_datasets <- setdiff(
  EXPECTED_DATASETS,
  file_ids
)

if (length(missing_datasets) > 0) {
  stop(
    "Missing DESeq2 result files: ",
    paste(missing_datasets, collapse = ", ")
  )
}

# Order files consistently
result_files <- result_files[
  match(EXPECTED_DATASETS, file_ids)
]
message("DESeq2 result files validated: ", length(result_files))

# Load and validate DESeq2 results
dataset_results <- vector("list", length(EXPECTED_DATASETS))
names(dataset_results) <- EXPECTED_DATASETS

for (i in seq_along(EXPECTED_DATASETS)) {
  
  # Dataset identifier and file
  dataset_id <- EXPECTED_DATASETS[i]
  result_file <- result_files[i]
  # Read DE results
  result <- read_csv(result_file, show_col_types = FALSE)
  
  # Required columns
  required_columns <- c("gene_id", "log2FoldChange", "lfcSE")
  # Validate columns
  missing_columns <- setdiff(required_columns, names(result))
  
  if (length(missing_columns) > 0) {
    stop(
      basename(result_file),
      " is missing: ",
      paste(missing_columns, collapse = ", ")
    )
  }
  
  # Keep valid effect estimates
  dataset_results[[dataset_id]] <- result %>%
    filter(
      is.finite(log2FoldChange),
      is.finite(lfcSE),
      lfcSE > 0
    ) %>%
    transmute(
      gene_id = as.character(gene_id),
      log2FoldChange = log2FoldChange,
      lfcSE = lfcSE,
      dataset_id = dataset_id
    )
  message(
    dataset_id, ": ",
    nrow(dataset_results[[dataset_id]]),
    " valid effect estimates"
  )
}

# Merge all datasets  
all_estimates <- bind_rows(dataset_results)

# Validate merged estimates
if (nrow(all_estimates) == 0) {
  stop("No valid effect estimates available.")
}

# Report total estimates
message("Total effect estimates: ", nrow(all_estimates))

# Run leave-one-tissue-out meta-analysis
run_loto <- function(excluded_tissue, excluded_datasets) {
  
  remaining <- all_estimates %>%
    filter(!dataset_id %in% excluded_datasets)
  
  if (length(unique(remaining$dataset_id)) < MIN_DATASETS) {
    stop("Too few datasets remain after excluding ", excluded_tissue)
  }
  
  gene_counts <- remaining %>%
    count(gene_id, name = "n_datasets") %>%
    filter(n_datasets >= MIN_DATASETS)
  
  eligible_genes <- intersect(
    gene_counts$gene_id,
    core_genes$gene_id
  )
  
  message(
    "\nExcluding ", excluded_tissue, " [",
    paste(excluded_datasets, collapse = ", "),
    "]; remaining datasets: ",
    length(unique(remaining$dataset_id)),
    "; eligible genes: ", length(eligible_genes)
  )
  
  results <- vector("list", length(eligible_genes))
  
  for (i in seq_along(eligible_genes)) {
    
    gene <- eligible_genes[i]
    
    gene_data <- remaining %>%
      filter(gene_id == gene) %>%
      arrange(dataset_id)
    
    model <- tryCatch(
      metafor::rma.uni(
        yi = gene_data$log2FoldChange,
        sei = gene_data$lfcSE,
        method = "REML",
        test = "z"
      ),
      error = function(e) e
    )
    
    if (inherits(model, "error")) {
      
      results[[i]] <- tibble(
        gene_id = gene,
        n_datasets = nrow(gene_data),
        datasets = paste(gene_data$dataset_id, collapse = ";"),
        pooled_log2FC = NA_real_,
        standard_error = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        p_value = NA_real_,
        tau2 = NA_real_,
        I2 = NA_real_,
        Q_p_value = NA_real_,
        model_status = "failed"
      )
      
    } else {
      
      results[[i]] <- tibble(
        gene_id = gene,
        n_datasets = nrow(gene_data),
        datasets = paste(gene_data$dataset_id, collapse = ";"),
        pooled_log2FC = as.numeric(model$beta[1]),
        standard_error = as.numeric(model$se[1]),
        ci_lower = as.numeric(model$ci.lb),
        ci_upper = as.numeric(model$ci.ub),
        p_value = as.numeric(model$pval),
        tau2 = as.numeric(model$tau2),
        I2 = as.numeric(model$I2),
        Q_p_value = as.numeric(model$QEp),
        model_status = "success"
      )
    }
  }
  
  bind_rows(results) %>%
    mutate(
      # Recalculate BH-FDR across all eligible genes in this iteration
      padj = p.adjust(p_value, method = "BH"),
      excluded_tissue = excluded_tissue,
      excluded_datasets = paste(excluded_datasets, collapse = ";"),
      remaining_datasets = length(unique(remaining$dataset_id))
    )
}

# Run six tissue-exclusion iterations
loto_results <- bind_rows(
  lapply(names(TISSUE_GROUPS), function(tissue) {
    run_loto(tissue, TISSUE_GROUPS[[tissue]])
  })
)

# Evaluate robustness of core genes
core_loto <- loto_results %>%
  inner_join(core_genes, by = "gene_id") %>%
  mutate(
    meta_significant = !is.na(padj) & padj < FDR_THRESHOLD,
    direction_preserved = case_when(
      !is.finite(pooled_log2FC) ~ NA,
      original_log2FC == 0 ~ NA,
      pooled_log2FC == 0 ~ FALSE,
      TRUE ~ sign(pooled_log2FC) == sign(original_log2FC)
    ),
    significant_and_direction_preserved =
      meta_significant & direction_preserved
  ) %>%
  arrange(excluded_tissue, gene_id)

# Summarize each tissue-exclusion iteration
tissue_summary <- core_loto %>%
  group_by(excluded_tissue, excluded_datasets) %>%
  summarise(
    n_core_genes = n(),
    n_significant = sum(meta_significant, na.rm = TRUE),
    n_direction_preserved = sum(direction_preserved, na.rm = TRUE),
    n_significant_direction_preserved =
      sum(significant_and_direction_preserved, na.rm = TRUE),
    n_model_failed = sum(model_status != "success"),
    .groups = "drop"
  ) %>%
  mutate(
    percent_significant = 100 * n_significant / n_core_genes,
    percent_significant_direction_preserved =
      100 * n_significant_direction_preserved / n_core_genes
  )

# Summarize robustness across all six iterations for each core gene
gene_summary <- core_loto %>%
  group_by(gene_id, gene_symbol, original_log2FC, original_direction) %>%
  summarise(
    n_tissue_exclusions = n(),
    n_significant = sum(meta_significant, na.rm = TRUE),
    n_direction_preserved = sum(direction_preserved, na.rm = TRUE),
    n_significant_direction_preserved =
      sum(significant_and_direction_preserved, na.rm = TRUE),
    robust_across_all_six =
      n() == length(TISSUE_GROUPS) &&
      all(significant_and_direction_preserved, na.rm = FALSE),
    .groups = "drop"
  )

message(
  "Robust core genes: ",
  sum(gene_summary$robust_across_all_six),
  " / ",
  nrow(gene_summary)
)

# Save results
write_csv(
  core_loto,
  file.path(TABLE_DIR, "loto_gene_level_results.csv")
)

write_csv(
  tissue_summary,
  file.path(TABLE_DIR, "loto_robustness_summary.csv")
)

write_csv(
  gene_summary,
  file.path(TABLE_DIR, "loto_core_gene_summary.csv")
)

# Identify non-robust core genes
non_robust_genes <- core_loto %>%
  filter(!significant_and_direction_preserved) %>%
  dplyr::select(
    gene_symbol,
    gene_id,
    excluded_tissue,
    pooled_log2FC,
    padj,
    meta_significant,
    direction_preserved
  )

print(non_robust_genes)

# Plot the number of core genes that remain significant and direction-consistent
plot_data <- tissue_summary %>%
  mutate(
    excluded_tissue = factor(
      excluded_tissue,
      levels = names(TISSUE_GROUPS)
    )
  )

p <- ggplot(
  plot_data,
  aes(
    x = excluded_tissue,
    y = n_significant_direction_preserved
  )
) +
  geom_col(fill = "#287D8E", width = 0.7) +
  geom_text(
    aes(label = paste0(
      n_significant_direction_preserved, "/", n_core_genes
    )),
    vjust = -0.4,
    size = 3.5
  ) +
  scale_y_continuous(
    limits = c(0, nrow(core_genes) * 1.12),
    breaks = pretty(c(0, nrow(core_genes)))
  ) +
  labs(
    title = "Robustness of core genes across leave-one-tissue-out analyses",
    subtitle = "Genes remaining significant (FDR < 0.05) with preserved effect direction",
    x = "Excluded tissue group",
    y = "Robust core genes"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(
  file.path(PLOT_DIR, "loto_robustness.png"),
  plot = p,
  width = 9,
  height = 5.5,
  dpi = 300
)

# Overall analysis summary
overall_summary <- tibble(
  n_core_genes = nrow(core_genes),
  n_tissue_groups = length(TISSUE_GROUPS),
  n_gene_tissue_tests = nrow(core_loto),
  n_core_genes_robust_across_all_six =
    sum(gene_summary$robust_across_all_six),
  fdr_threshold = FDR_THRESHOLD,
  meta_analysis_method = "REML",
  p_value_test = "z",
  multiple_testing_correction = "Benjamini-Hochberg"
)

write_csv(
  overall_summary,
  file.path(TABLE_DIR, "loto_analysis_summary.csv")
)

message("\nLeave-one-tissue-out analysis completed.")
message("Output directory: ", OUTPUT_DIR)
message("Core genes: ", nrow(core_genes))
message("Robust across all six tissue exclusions: ",
  sum(gene_summary$robust_across_all_six), "/", nrow(core_genes)
)