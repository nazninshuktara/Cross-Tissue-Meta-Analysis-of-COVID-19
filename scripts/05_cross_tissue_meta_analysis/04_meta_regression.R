#!/usr/bin/env Rscript

# 04_meta_regression.R
# Explore study-level moderators of cross-tissue gene effects.

suppressPackageStartupMessages({
  library(readr)      # Import data
  library(dplyr)      # Data wrangling
  library(tibble)     # Tibble utilities
  library(metafor)    # Meta-regression
})

# Define project paths
PROJECT_ROOT <- "."

# Input directories and files
# Directory containing DESeq2 results
DE_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "differential_expression")

# Core gene classification
META_FILE <- file.path(
  PROJECT_ROOT,
  "outputs",
  "meta_analysis",
  "core_tissue_specific_decomposition",
  "tables",
  "gene_classification.csv")

# Output directories
OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "meta_analysis",
  "meta_regression")

TABLE_DIR <- file.path(OUTPUT_DIR, "tables")
PLOT_DIR  <- file.path(OUTPUT_DIR, "plots")

# Create output directories
dir.create(OUTPUT_DIR, recursive = TRUE,showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)

# Analysis settings
FDR_THRESHOLD <- 0.05
MIN_DATASETS  <- 2

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
  "GSE182917")

# Tissue groups
TISSUE_GROUPS <- list(
  Nasopharyngeal = c("GSE208342"),
  Lung = c("GSE182917", "GSE208076"),
  Blood = c("GSE211979", "GSE273149", "GSE243217"),
  Kidney = c("GSE202182"),
  Heart = c("GSE169241"),
  Gut_organoids = c("GSE149312"))

# Load gene classification
gene_classes <- read_csv(
  META_FILE,
  show_col_types = FALSE)

# Required columns
required_columns <- c(
  "gene_id",
  "gene_symbol",
  "pooled_log2FC",
  "component")

# Validate classification file
if (!all(required_columns %in% names(gene_classes))) {
  stop("Gene classification file is missing required columns.")
}

# Retain core genes
core_genes <- gene_classes %>%
  filter(component == "core") %>%
  dplyr::select(
    gene_id,
    gene_symbol,
    pooled_log2FC)

core_gene_ids <- core_genes$gene_id
message("Core genes loaded: ", nrow(core_genes))

# Load study-level moderators
study_mods <- read_csv(
  "outputs/meta_analysis/tables/meta_regression/study_level_moderators.csv",
  show_col_types = FALSE)

message("Studies with moderator data: ", nrow(study_mods))

# Locate DESeq2 result files
result_files <- list.files(
  DE_DIR,
  pattern = "^DESeq2_GSE[0-9]+_COVID_vs_Control\\.csv$",
  full.names = TRUE)

file_ids <- sub(
  "^DESeq2_(GSE[0-9]+)_COVID_vs_Control\\.csv$",
  "\\1",
  basename(result_files))

result_files <- result_files[
  match(EXPECTED_DATASETS, file_ids)]

# Load DESeq2 effect estimates
dataset_results <- vector(
  "list",
  length(EXPECTED_DATASETS)
)

names(dataset_results) <- EXPECTED_DATASETS

for (i in seq_along(EXPECTED_DATASETS)) {
  
  dataset_id <- EXPECTED_DATASETS[i]
  
  result <- read_csv(
    result_files[i],
    show_col_types = FALSE
  )
  
  dataset_results[[dataset_id]] <- result %>%
    filter(
      is.finite(log2FoldChange),
      is.finite(lfcSE),
      lfcSE > 0
    ) %>%
    transmute(
      gene_id = as.character(gene_id),
      log2FoldChange,
      lfcSE,
      dataset_id
    )
}

all_estimates <- bind_rows(dataset_results)

# Keep core genes only
all_estimates <- all_estimates %>%
  filter(
    gene_id %in% core_genes$gene_id
  )

# Add study-level moderators
all_estimates <- all_estimates %>%
  inner_join(
    study_mods,
    by = "dataset_id")

message("Core-gene effect estimates: ",nrow(all_estimates))

# Fit one moderator model for one gene
fit_one_model <- function(gene, moderator, dat) {
  base_result <- tibble(
    gene_id = gene,
    moderator = moderator,
    k = nrow(dat),
    estimate = NA_real_,
    standard_error = NA_real_,
    tau2 = NA_real_,
    QM = NA_real_,
    QM_df = NA_real_,
    pvalue = NA_real_,
    warning_message = NA_character_,
    error_message = NA_character_
  )
  
  warnings_seen <- character()
  
  formula_mods <- switch(
    moderator,
    "Case/control ratio" = ~ case_control_ratio,
    "Total sample size (per 10 samples)" = ~ total_n_per10,
    "Tissue group (exploratory)" = ~ factor(tissue_group)
  )
  
  fit <- tryCatch(
    withCallingHandlers(
      metafor::rma.uni(
        yi = log2FoldChange,
        sei = lfcSE,
        mods = formula_mods,
        data = dat,
        method = "REML",
        test = "z"
      ),
      warning = function(w) {
        warnings_seen <<- c(warnings_seen, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      base_result$error_message <<- conditionMessage(e)
      NULL
    }
  )
  
  if (!is.null(fit)) {
    # Extract scalar model statistics safely.
    base_result$tau2 <- as.numeric(fit$tau2)[1]
    base_result$QM <- as.numeric(fit$QM)[1]
    base_result$QM_df <- fit$m
    base_result$pvalue <- as.numeric(fit$QMp)[1]
    
    # Extract the slope for continuous moderators.
    if (
      moderator != "Tissue group (exploratory)" &&
      length(fit$b) >= 2 &&
      length(fit$se) >= 2
    ) {
      base_result$estimate <- as.numeric(fit$b[2, 1])
      base_result$standard_error <- as.numeric(fit$se[2])
    }
  }
  
  if (length(warnings_seen) > 0) {
    base_result$warning_message <- paste(
      unique(warnings_seen),
      collapse = " | "
    )
  }
  
  base_result
}

# Fit models separately for each core gene and moderator
moderators <- c(
  "Case/control ratio",
  "Tissue group (exploratory)",
  "Total sample size (per 10 samples)"
)

results <- lapply(core_gene_ids, function(gene) {
  gene_data <- all_estimates %>%
    filter(gene_id == gene)
  
  lapply(moderators, function(moderator) {
    if (moderator == "Case/control ratio") {
      dat <- gene_data %>%
        filter(is.finite(case_control_ratio))
    } else if (moderator == "Total sample size (per 10 samples)") {
      dat <- gene_data %>%
        filter(is.finite(total_n_per10))
    } else {
      dat <- gene_data %>%
        filter(!is.na(tissue_group)) %>%
        mutate(tissue_group = factor(tissue_group))
    }
    
    fit_one_model(gene, moderator, dat)
  }) %>%
    bind_rows()
}) %>%
  bind_rows() %>%
  group_by(moderator) %>%
  mutate(
    padj = ifelse(
      is.finite(pvalue),
      p.adjust(pvalue, method = "BH"),
      NA_real_
    )
  ) %>%
  ungroup()

# Save gene-level model results
write_csv(
  results,
  file.path(TABLE_DIR, "meta_regression_gene_level.csv"))

# Summarise model status and multiple-testing results
summary_table <- results %>%
  group_by(moderator) %>%
  summarise(
    total_models = n(),
    tested = sum(is.finite(pvalue)),
    failed = sum(!is.na(error_message)),
    warnings = sum(!is.na(warning_message)),
    n_FDR_below_0.05 = sum(padj < 0.05, na.rm = TRUE),
    min_pvalue = if (any(is.finite(pvalue))) {
      min(pvalue, na.rm = TRUE)
    } else {
      NA_real_
    },
    min_FDR = if (any(is.finite(padj))) {
      min(padj, na.rm = TRUE)
    } else {
      NA_real_
    },
    .groups = "drop"
  )

write_csv(
  summary_table,
  file.path(TABLE_DIR, "meta_regression_summary.csv")
)

# Print report
cat("\nMeta-regression completed.\n\n")
print(summary_table)