#!/usr/bin/env Rscript

# 00_primary_sample_distribution.R
# Summarize and visualize primary-analysis sample distribution.

# Load required packages
suppressPackageStartupMessages({
  library(ggplot2)       # Create publication-quality figures
  library(dplyr)         # Manipulate and summarize sample metadata
  library(readr)         # Read and write CSV files
  library(tidyr)         # Reshape and complete sample-count tables
  library(scales)        # Format plot axes and tick marks
  library(gtsummary)     # Format publication-ready tables
  library(gt)            # Create Word documents
  library(flextable)
  library(officer)
})

# Define project root
PROJECT_ROOT <- "."

# Define input file
# Curated sample-level metadata.
METADATA_FILE <- file.path(
  PROJECT_ROOT,
  "data",
  "metadata",
  "metadata_curated_sample.csv"
)

# Define output directories
OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "primary_analysis_sample_distribution"
)

# Directory for storing plots and tables
PLOT_DIR <- file.path(OUTPUT_DIR, "plots")
TABLE_DIR <- file.path(OUTPUT_DIR, "tables")

# Create output directories if they do not exist
dir.create(PLOT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# Define discovery dataset order
datasets <- c(
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

# Define tissue labels for table and figure
tissue_mapping <- c(
  "GSE208342" = "Nasopharyngeal",
  "GSE273149" = "PBMC",
  "GSE243217" = "Peripheral blood",
  "GSE202182" = "Kidney",
  "GSE208076" = "Lung",
  "GSE211979" = "Blood",
  "GSE169241" = "Heart",
  "GSE149312" = "Gut organoid",
  "GSE182917" = "Lung"
)

# Define dataset-specific control descriptions
control_mapping <- c(
  "GSE208342" = "SARS-CoV-2 test-negative",
  "GSE273149" = "Non-COVID ARDS",
  "GSE243217" = "Healthy donors",
  "GSE202182" = "Healthy kidney",
  "GSE208076" = "Normal lung",
  "GSE211979" = "Healthy donors",
  "GSE169241" = "Uninfected heart",
  "GSE149312" = "Untreated organoids",
  "GSE182917" = "Non-COVID lung"
)

# Read curated sample metadata
metadata <- read_csv(
  METADATA_FILE,
  show_col_types = FALSE
)

# Validate required metadata columns
required_columns <- c(
  "dataset_id",
  "disease_group",
  "inclusion_status",
  "tissue_type"
)

missing_columns <- setdiff(
  required_columns,
  names(metadata)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

# Select samples included in the primary analysis
primary_samples <- metadata %>%
  filter(
    inclusion_status == "included",
    disease_group %in% c("COVID-19", "Control"),
    dataset_id %in% datasets
  )

if (nrow(primary_samples) == 0) {
  stop("No primary-analysis samples found. Check metadata labels.")
}

# Check tissue types represented in the primary analysis
tissue_check <- primary_samples %>%
  distinct(dataset_id, tissue_type) %>%
  arrange(dataset_id, tissue_type)

print(tissue_check, n = Inf)

# Count biological samples by dataset and disease group
sample_counts <- primary_samples %>%
  count(dataset_id, disease_group, name = "sample_count") %>%
  complete(
    dataset_id = datasets,
    disease_group = c("COVID-19", "Control"),
    fill = list(sample_count = 0)
  ) %>%
  mutate(
    dataset_id = factor(dataset_id, levels = datasets),
    disease_group = factor(
      disease_group,
      levels = c("COVID-19", "Control")
    )
  ) %>%
  arrange(dataset_id, disease_group)

# Create the Primary analysis sample distribution table
table1 <- sample_counts %>%
  mutate(
    dataset_id = as.character(dataset_id),
    disease_group = as.character(disease_group)
  ) %>%
  select(dataset_id, disease_group, sample_count) %>%
  pivot_wider(
    names_from = disease_group,
    values_from = sample_count,
    values_fill = 0
  )

# Ensure both disease-group columns exist
if (!"COVID-19" %in% names(table1)) {
  table1[["COVID-19"]] <- 0L
}

if (!"Control" %in% names(table1)) {
  table1[["Control"]] <- 0L
}

# Format the primary analysis sample distribution table 
sample_distribution_table <- table1 %>%
  mutate(
    Tissue = unname(tissue_mapping[dataset_id]),
    `Total samples (N)` = `COVID-19` + Control,
    `COVID-19 (n)` = `COVID-19`,
    `Controls (n)` = Control,
    `Control type` = unname(control_mapping[dataset_id])
  ) %>%
  select(
    `Dataset ID` = dataset_id,
    `Tissue Type` = Tissue,
    `COVID-19 (n)`,
    `Controls (n)`,
    `Total samples (N)`,
    `Control type`
  )

# Add the total row
total_row <- tibble(
  `Dataset ID` = "Total",
  `Tissue Type` = "",
  `COVID-19 (n)` = sum(sample_distribution_table$`COVID-19 (n)`),
  `Controls (n)` = sum(sample_distribution_table$`Controls (n)`),
  `Total samples (N)` = sum(sample_distribution_table$`Total samples (N)`),
  `Control type` = ""
)

sample_distribution_table <- bind_rows(sample_distribution_table, total_row)

# Validate the sample counts
if (nrow(sample_distribution_table) != length(datasets) + 1) {
  stop("Unexpected number of rows in the sample characteristics table.")
}

if (sum(sample_distribution_table$`Total samples (N)`[-nrow(sample_distribution_table)]) !=
    nrow(primary_samples)) {
  stop("Sample counts do not match the primary-analysis metadata.")
}

# Create a publication-ready 'primary analysis sample distribution' table
sample_distribution <- sample_distribution_table %>%
  gt() %>%
  tab_header(
    title = "Table 1. Primary Analysis Sample Distribution",
    subtitle = "Summary of the nine discovery datasets for primary analysis: dataset ID, tissue type, number of COVID-19(cases) and Controls, total number of primary analysis samples, and control type."
  ) %>%
  cols_label(
    `Dataset ID` = "Dataset ID",
    `Tissue Type` = "Tissue Type",
    `COVID-19 (n)` = "COVID-19 (n)",
    `Controls (n)` = "Controls (n)",
    `Total samples (N)` = "Total (N)",
    `Control type` = "Control type"
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
    locations = cells_column_labels(
      columns = everything()
    )
  ) %>%
  # Bold the total row
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      rows = `Dataset ID` == "Total"
    )
  ) %>%
  # Add vertical separators between columns
  tab_style(
    style = cell_borders(
      sides = "right",
      color = "#B7C3D0",
      weight = px(1)
    ),
    locations = cells_body(
      columns = everything()
    )
  ) %>%
  # Add horizontal separators between rows
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "#B7C3D0",
      weight = px(1)
    ),
    locations = cells_body(
      columns = everything()
    )
  ) %>%
  # Add separators to the header as well
  tab_style(
    style = cell_borders(
      sides = "right",
      color = "#FFFFFF",
      weight = px(1)
    ),
    locations = cells_column_labels(
      columns = everything()
    )
  ) %>%
  # Add a subtle background to the total row
  tab_style(
    style = cell_fill(color = "#EAF0F6"),
    locations = cells_body(
      rows = `Dataset ID` == "Total"
    )
  ) %>%
  # Add explanatory footnote
  tab_source_note(
    source_note = paste(
      "Sample counts represent samples selected for the primary analysis based on metadata curation, before quality control."
    )
  ) %>%
  # Italicize the footnote
  tab_style(
    style = cell_text(
      style = "italic",
      size = px(11)
    ),
    locations = cells_source_notes()
  )

# Save the formatted table as HTML
gtsave(
  sample_distribution,
  filename = file.path(
    TABLE_DIR,
    "primary_sample_distribution.html"
  )
)

# Create the primary analysis sample-distribution figure
plot_data <- sample_counts %>%
  mutate(
    dataset_id = as.character(dataset_id),
    tissue = unname(tissue_mapping[dataset_id]),
    disease_group = factor(
      disease_group,
      levels = c("COVID-19", "Control")
    )
  )

sample_distribution_plot <- ggplot(
  plot_data,
  aes(
    x = dataset_id,
    y = sample_count,
    fill = disease_group
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7,
    color = "black",
    linewidth = 0.4
  ) +
  geom_text(
    aes(label = ifelse(sample_count > 0, sample_count, "")),
    position = position_dodge(width = 0.8),
    vjust = -0.35,
    size = 3.2
  ) +
  scale_fill_manual(
    values = c(
      "COVID-19" = "#C0392B",
      "Control" = "#4C72B0"
    ),
    labels = c(
      "COVID-19" = "COVID-19",
      "Control" = "Control"
    )
  ) +
  scale_x_discrete(
    limits = datasets,
    labels = setNames(
      paste0(datasets, "\n", tissue_mapping[datasets]),
      datasets
    )
  ) +
  scale_y_continuous(
    breaks = scales::breaks_pretty(),
    expand = expansion(mult = c(0, 0.12))
  ) +
  labs(
    title = "Primary-Analysis sample distribution",
    x = "Discovery datasets (tissue type)",
    y = "Number of biological samples",
    fill = "Disease Group"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      vjust = 1
    ),
    legend.position = "top",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

# Save the figure as a high-resolution PNG
ggsave(
  filename = file.path(PLOT_DIR, "primary_sample_distribution.png"),
  plot = sample_distribution_plot,
  width = 12,
  height = 6.5,
  dpi = 600
)


# Print the final outputs in the console
print(sample_distribution_table)
print(sample_distribution_plot)

message("Sample-distribution table and figures created successfully.")