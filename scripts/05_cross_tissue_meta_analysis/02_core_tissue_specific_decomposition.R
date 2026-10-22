#!/usr/bin/env Rscript

# 02_core_tissue_specific_decomposition.R
# Classify meta-significant genes into core, shared and tissue-specific components
# Perform GO, KEGG and GSEA analyses.

# Load required packages
suppressPackageStartupMessages({
  # Data import, wrangling, and table manipulation
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  
  # Visualization and figure assembly
  library(ggplot2)
  library(stringr)
  library(patchwork)
  
  # Specialized plotting utilities
  library(ggridges)
  library(gridExtra)
  
  # Functional enrichment analysis
  library(clusterProfiler)
  
  # Human gene annotation databases
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  
})

# Define project paths
PROJECT_ROOT <- "."

# Input directories
# Directory containing complete DESeq2 result files.
DE_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "differential_expression"
)

# Directory containing random-effects meta-analysis results
META_DIR <- file.path(
  PROJECT_ROOT, 
  "outputs", 
  "meta_analysis",
  "random_effects_meta_analysis"
)

# Output directories for decomposition results, enrichment analyses, tables, and figures.
OUTPUT_DIR <- file.path(
  PROJECT_ROOT,
  "outputs",
  "meta_analysis",
  "core_tissue_specific_decomposition"
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


# Define the nine discovery datasets
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

# Classification criteria:
# Core: >=4 significant datasets, I² <25%, concordance >=80%
# Shared: >=2 significant datasets but not core
# Tissue-specific candidate: <2 significant datasets
FDR_THRESHOLD <- 0.05
CORE_MIN_DATASETS <- 4
CORE_MAX_I2 <- 25
CORE_MIN_CONCORDANCE <- 0.80
MIN_GO_GENES <- 5

# Define input files
META_FILE <- file.path(
  META_DIR,
  "tables",
  "random_effects_meta_analysis.csv"
)

CONCORDANCE_FILE <- file.path(
  META_DIR,
  "tables",
  "cross_dataset_effect_direction_concordance.csv"
)

# Check that required inputs exist
if (!file.exists(META_FILE)) {
  stop("Meta-analysis results not found: ", META_FILE)
}

if (!file.exists(CONCORDANCE_FILE)) {
  stop("Direction concordance results not found: ", CONCORDANCE_FILE)
}

# Locate the differential expression results for all datasets
DE_FILES <- file.path(
  DE_DIR,
  paste0("DESeq2_", DATASETS, "_COVID_vs_Control.csv")
)

# Check whether any differential expression files are missing
missing_files <- DE_FILES[!file.exists(DE_FILES)]

if (length(missing_files) > 0) {
  stop(
    "Missing differential expression files:\n",
    paste(missing_files, collapse = "\n")
  )
}

# Read the random-effects meta-analysis results
meta_results <- read_csv(
  META_FILE,
  show_col_types = FALSE
)

# Read direction concordance results
concordance_results <- read_csv(
  CONCORDANCE_FILE,
  show_col_types = FALSE
)

# Report the dimensions of the input tables
message(
  "Meta-analysis table: ",
  nrow(meta_results), " rows, ",
  ncol(meta_results), " columns"
)

message(
  "Concordance table: ",
  nrow(concordance_results), " rows, ",
  ncol(concordance_results), " columns"
)

# Validate required columns
required_meta_columns <- c(
  "gene_id",
  "n_datasets",
  "pooled_log2FC",
  "I2",
  "padj",
  "meta_significant"
)

required_concordance_columns <- c(
  "gene_id",
  "significant",
  "n_datasets",
  "n_concordant",
  "direction_concordance"
)

# Check the meta-analysis table
if (!all(required_meta_columns %in% names(meta_results))) {
  stop(
    "Meta-analysis file is missing required columns: ",
    paste(
      setdiff(required_meta_columns, names(meta_results)),
      collapse = ", "
    )
  )
}

# Check the direction concordance table
if (!all(required_concordance_columns %in% names(concordance_results))) {
  stop(
    "Concordance file is missing required columns: ",
    paste(
      setdiff(required_concordance_columns, names(concordance_results)),
      collapse = ", "
    )
  )
}

message("Input column validation passed.")

# Combine DESeq2 results from all nine datasets
de_results_list <- lapply(seq_along(DE_FILES), function(i) {
  
  read_csv(
    DE_FILES[i],
    show_col_types = FALSE
  ) %>%
    mutate(dataset = DATASETS[i])
})

de_all <- bind_rows(de_results_list)

# Report the combined results
message("Total rows: ", nrow(de_all))
message("Unique genes: ", n_distinct(de_all$gene_id))
message("Datasets: ", n_distinct(de_all$dataset))


# Summarize dataset-level evidence for all genes
gene_evidence_all <- de_all %>%
  group_by(gene_id) %>%
  summarise(
    n_de_datasets = sum(
      !is.na(padj) & is.finite(padj) & padj < FDR_THRESHOLD
    ),
    .groups = "drop"
  )

message("Genes with dataset-level evidence: ", nrow(gene_evidence_all))

# Retain the existing meta-significant gene evidence for classification
meta_significant_genes <- meta_results %>%
  dplyr::filter(meta_significant %in% TRUE) %>%
  dplyr::select(gene_id)

gene_evidence <- gene_evidence_all %>%
  semi_join(meta_significant_genes, by = "gene_id")

message("Meta-significant genes: ", nrow(gene_evidence))


# Merge dataset-level evidence with meta-analysis and concordance results
gene_evidence_merged <- gene_evidence %>%
  dplyr::left_join(
    meta_results %>%
      dplyr::select(
        gene_id,
        n_datasets,
        pooled_log2FC,
        I2,
        padj,
        meta_significant
      ),
    by = "gene_id"
  ) %>%
  dplyr::left_join(
    concordance_results %>%
      dplyr::select(
        gene_id,
        n_concordant,
        direction_concordance
      ),
    by = "gene_id"
  )

# Check merged results
print(head(gene_evidence_merged, 10))

message(
  "Genes after merging: ",
  nrow(gene_evidence_merged)
)

message(
  "Genes with missing I2: ",
  sum(is.na(gene_evidence_merged$I2))
)

message(
  "Genes with missing direction concordance: ",
  sum(is.na(gene_evidence_merged$direction_concordance))
)

# Supplementary analysis:
# Identify concordant genes that are not meta-significant
concordant_not_sig <- gene_evidence_all %>%
  left_join(
    meta_results %>%
      dplyr::select(
        gene_id, n_datasets, pooled_log2FC,
        I2, padj, meta_significant
      ),
    by = "gene_id"
  ) %>%
  left_join(
    concordance_results %>%
      dplyr::select(
        gene_id, n_concordant, direction_concordance
      ),
    by = "gene_id"
  )

# Apply classification criteria
keep <- !is.na(concordant_not_sig[["meta_significant"]]) &
  concordant_not_sig[["meta_significant"]] == FALSE &
  concordant_not_sig[["n_de_datasets"]] >= 2 &
  !is.na(concordant_not_sig[["direction_concordance"]]) &
  concordant_not_sig[["direction_concordance"]] >= CORE_MIN_CONCORDANCE

concordant_not_sig <- concordant_not_sig[keep, ] %>%
  mutate(component = "concordant_not_sig")

# Save the supplementary gene table
readr::write_csv(
  concordant_not_sig,
  file.path(TABLE_DIR, "concordant_not_meta_significant_genes.csv"))

message("Concordant but non-significant genes: ", nrow(concordant_not_sig))

# Classify meta-significant genes into core, shared and tissue specific components
gene_classification <- gene_evidence_merged %>%
  mutate(
    component = case_when(
      
      n_de_datasets >= CORE_MIN_DATASETS &
        !is.na(I2) &
        I2 < CORE_MAX_I2 &
        !is.na(direction_concordance) &
        direction_concordance >= CORE_MIN_CONCORDANCE
      ~ "core",
      
      n_de_datasets >= 2
      ~ "shared",
      
      n_de_datasets < 2
      ~ "tissue_specific_candidate"
    )
  )

# Summarize component counts
classification_summary <- gene_classification %>%
  dplyr::count(component, name = "n_genes") %>%
  dplyr::arrange(component)

print(classification_summary)

# Confirm that all meta-significant genes are retained
stopifnot(
  nrow(gene_classification) == nrow(gene_evidence_merged),
  !anyNA(gene_classification$component)
)

message(
  "Total classified genes: ",
  nrow(gene_classification)
)

# Save the results
readr::write_csv(
  gene_classification,
  file.path(TABLE_DIR, "gene_classification.csv")
)

# Summarize gene counts by component
classification_summary <- gene_classification %>%
  dplyr::count(component, name = "n_genes")

print(classification_summary)

# GO Biological Process enrichment analysis
# Prepare Ensembl gene IDs by removing version suffixes
gene_classification <- gene_classification %>%
  mutate(
    ensembl_id = sub("\\.[0-9]+$", "", gene_id)
  )
# Map Ensembl IDs to Entrez IDs.
id_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(gene_classification$ensembl_id),
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

# Attach Entrez IDs to the existing classification table.
gene_classification <- gene_classification %>%
  mutate(
    ENTREZID = unname(id_map[ensembl_id])
  )

# Check mapping success.
cat("Total genes:", nrow(gene_classification), "\n")
cat("Mapped genes:", sum(!is.na(gene_classification$ENTREZID)), "\n")
cat("Unmapped genes:", sum(is.na(gene_classification$ENTREZID)), "\n")

# Prepare the background gene universe
background_genes <- meta_results %>%
  mutate(
    ensembl_id = sub("\\.[0-9]+$", "", gene_id)
  ) %>%
  filter(!is.na(ensembl_id), ensembl_id != "") %>%
  distinct(ensembl_id)

# Map background Ensembl IDs to Entrez IDs
background_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = background_genes$ensembl_id,
  keytype = "ENSEMBL",
  column = "ENTREZID",
  multiVals = "first"
)

# Keep unique, successfully mapped Entrez IDs
gene_universe <- unique(
  na.omit(unname(background_map))
)

# Check the background mapping
cat("Background Ensembl genes:", nrow(background_genes), "\n")
cat("Mapped background Entrez IDs:", length(gene_universe), "\n")

# Prepare unique Entrez IDs for each component
core_ids <- gene_classification %>%
  filter(component == "core", !is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

shared_ids <- gene_classification %>%
  filter(component == "shared", !is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

tissue_specific_ids <- gene_classification %>%
  filter(
    component == "tissue_specific_candidate",
    !is.na(ENTREZID)
  ) %>%
  pull(ENTREZID) %>%
  unique()

# Keep only genes present in the background universe
core_ids <- intersect(core_ids, gene_universe)
shared_ids <- intersect(shared_ids, gene_universe)
tissue_specific_ids <- intersect(tissue_specific_ids, gene_universe)

# Report the final gene-set sizes
cat("Core genes:", length(core_ids), "\n")
cat("Shared genes:", length(shared_ids), "\n")
cat("Tissue-specific candidate genes:", length(tissue_specific_ids), "\n")
cat("Background genes:", length(gene_universe), "\n")

# Run GO Biological Process enrichment for each component
run_go_enrichment <- function(gene_ids, universe_ids) {
  
  # Remove missing IDs and duplicates
  gene_ids <- unique(gene_ids[!is.na(gene_ids)])
  universe_ids <- unique(universe_ids[!is.na(universe_ids)])
  
  # Restrict input genes to the background universe
  gene_ids <- intersect(gene_ids, universe_ids)
  
  # Skip enrichment if fewer than MIN_GO_GENES remain
  if (length(gene_ids) < MIN_GO_GENES) {
    message("Insufficient genes for GO enrichment.")
    return(NULL)
  }
  
  clusterProfiler::enrichGO(
    gene          = gene_ids,
    universe      = universe_ids,
    OrgDb         = org.Hs.eg.db,
    keyType       = "ENTREZID",
    ont           = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    minGSSize     = MIN_GO_GENES,
    readable      = TRUE
  )
}

# Run enrichment for all three components
enrich_core <- run_go_enrichment(core_ids, gene_universe)
enrich_shared <- run_go_enrichment(shared_ids, gene_universe)
enrich_tissue_specific <- run_go_enrichment(tissue_specific_ids, gene_universe)

# Count enriched GO terms
cat("Core enriched terms:",
    if (is.null(enrich_core)) 0 else nrow(as.data.frame(enrich_core)), "\n")

cat("Shared enriched terms:",
    if (is.null(enrich_shared)) 0 else nrow(as.data.frame(enrich_shared)), "\n")

cat("Tissue-specific candidate enriched terms:",
    if (is.null(enrich_tissue_specific)) 0 else nrow(as.data.frame(enrich_tissue_specific)), "\n")

# Combine enrichment results into one table
enrichment_results <- bind_rows(
  if (!is.null(enrich_core))
    mutate(as.data.frame(enrich_core),
           component = "core"),
  
  if (!is.null(enrich_shared))
    mutate(as.data.frame(enrich_shared),
           component = "shared"),
  
  if (!is.null(enrich_tissue_specific))
    mutate(as.data.frame(enrich_tissue_specific),
           component = "tissue_specific_candidate")
)

# Save the combined enrichment table
readr::write_csv(
  enrichment_results,
  file.path(TABLE_DIR, "go_bp_enrichment_by_component.csv")
)

# Select the top 10 enriched GO terms per component
top_go_terms <- enrichment_results |>
  dplyr::filter(!is.na(p.adjust), !is.na(Description)) |>
  dplyr::group_by(component) |>
  dplyr::slice_min(
    order_by = p.adjust,
    n = 10,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

# Create a dot plot for each component
plot_go_component <- function(data, component_name) {
  
  plot_data <- data |>
    dplyr::filter(component == component_name) |>
    dplyr::mutate(
      Description = factor(
        Description,
        levels = rev(Description[order(GeneRatio)])
      )
    )
  
  if (nrow(plot_data) == 0) {
    message("No enriched GO terms for: ", component_name)
    return(NULL)
  }
  
  ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = GeneRatio,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    ggplot2::geom_point(alpha = 0.9) +
    ggplot2::scale_color_gradient(
      low = "firebrick",
      high = "steelblue",
      trans = "reverse"
    ) +
    ggplot2::scale_y_discrete(
      labels = function(x) stringr::str_wrap(x, width = 45)
    ) +
    ggplot2::labs(
      title = paste("GO Biological Process:", component_name),
      subtitle = "Top 10 terms ranked by adjusted p-value",
      x = "Gene ratio",
      y = NULL,
      color = "Adjusted p-value",
      size = "Gene count"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10),
      panel.grid.minor = ggplot2::element_blank()
    )
}

# Generate separate plots
p_core <- plot_go_component(top_go_terms, "core")
p_shared <- plot_go_component(top_go_terms, "shared")

# Save Core GO enrichment plot
if (!is.null(p_core)) {
  ggplot2::ggsave(
    filename = file.path(PLOT_DIR, "GO_BP_Core.png"),
    plot = p_core,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

# Save Shared GO enrichment plot
if (!is.null(p_shared)) {
  ggplot2::ggsave(
    filename = file.path(PLOT_DIR, "GO_BP_Shared.png"),
    plot = p_shared,
    width = 11,
    height = 9,
    dpi = 300,
    bg = "white"
  )
}

# Combine all plots vertically
combined_go_plot <- patchwork::wrap_plots(
  p_core,
  p_shared,
  ncol = 1,
  heights = c(1, 1)
)

# Save the combined figure
ggplot2::ggsave(
  filename = file.path(PLOT_DIR, "GO_BP_All_Components.png"),
  plot = combined_go_plot,
  width = 13,
  height = 16,
  dpi = 300,
  bg = "white"
)

# Run KEGG pathway enrichment
run_kegg_enrichment <- function(gene_ids, universe_ids) {
  
  # Remove missing IDs and duplicates
  gene_ids <- unique(as.character(gene_ids[!is.na(gene_ids)]))
  universe_ids <- unique(
    as.character(universe_ids[!is.na(universe_ids)])
  )
  
  # Restrict input genes to the background universe
  gene_ids <- intersect(gene_ids, universe_ids)
  
  # Skip if too few genes remain
  if (length(gene_ids) < MIN_GO_GENES) {
    message("Insufficient genes for KEGG enrichment.")
    return(NULL)
  }
  
  clusterProfiler::enrichKEGG(
    gene          = gene_ids,
    universe      = universe_ids,
    organism      = "hsa",
    keyType       = "ncbi-geneid",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.10,
    minGSSize     = MIN_GO_GENES
  )
}

# Run KEGG enrichment for all three components
kegg_core <- run_kegg_enrichment(core_ids, gene_universe)

kegg_shared <- run_kegg_enrichment(shared_ids, gene_universe)

kegg_tissue_specific <- run_kegg_enrichment(
  tissue_specific_ids,
  gene_universe
)

# Count significant KEGG pathways
cat(
  "Core KEGG pathways:",
  if (is.null(kegg_core)) 0 else nrow(as.data.frame(kegg_core)),
  "\n"
)

cat(
  "Shared KEGG pathways:",
  if (is.null(kegg_shared)) 0 else nrow(as.data.frame(kegg_shared)),
  "\n"
)

cat(
  "Tissue-specific candidate KEGG pathways:",
  if (is.null(kegg_tissue_specific)) 0 else
    nrow(as.data.frame(kegg_tissue_specific)),
  "\n"
)

# Preview significant KEGG pathways
if (!is.null(kegg_core)) {
  head(as.data.frame(kegg_core)[
    , c("ID", "Description", "GeneRatio", "Count", "p.adjust")
  ])
}

if (!is.null(kegg_shared)) {
  head(as.data.frame(kegg_shared)[
    , c("ID", "Description", "GeneRatio", "Count", "p.adjust")
  ])
}

if (!is.null(kegg_tissue_specific)) {
  head(as.data.frame(kegg_tissue_specific)[
    , c("ID", "Description", "GeneRatio", "Count", "p.adjust")
  ])
}
# Combine non-empty KEGG enrichment results
kegg_results <- dplyr::bind_rows(
  if (!is.null(kegg_core) && nrow(as.data.frame(kegg_core)) > 0) {
    dplyr::mutate(as.data.frame(kegg_core), component = "core")
  },
  if (!is.null(kegg_shared) && nrow(as.data.frame(kegg_shared)) > 0) {
    dplyr::mutate(as.data.frame(kegg_shared), component = "shared")
  },
  if (
    !is.null(kegg_tissue_specific) &&
    nrow(as.data.frame(kegg_tissue_specific)) > 0
  ) {
    dplyr::mutate(
      as.data.frame(kegg_tissue_specific),
      component = "tissue_specific_candidate"
    )
  }
)

# Check combined results
dim(kegg_results)
head(kegg_results)

# Save KEGG enrichment results
readr::write_csv(
  kegg_results,
  file.path(
    "outputs/meta_analysis/core_tissue_specific_decomposition/tables",
    "kegg_enrichment_by_component.csv"
  )
)

# Select the top 10 KEGG pathways per component
top_kegg <- kegg_results %>%
  dplyr::filter(
    !is.na(p.adjust),
    !is.na(Description)
  ) %>%
  dplyr::group_by(component) %>%
  dplyr::slice_min(
    order_by = p.adjust,
    n = 10,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

# Create KEGG dot plot
plot_kegg_component <- function(data, component_name) {
  
  plot_data <- data %>%
    dplyr::filter(component == component_name)
  
  if (nrow(plot_data) == 0) {
    message("No KEGG pathways for: ", component_name)
    return(NULL)
  }
  
  plot_data <- plot_data %>%
    dplyr::mutate(
      Description = factor(
        Description,
        levels = rev(
          Description[order(p.adjust)]
        )
      )
    )
  
  ggplot(
    plot_data,
    aes(
      x = Count,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point(alpha = 0.9) +
    scale_color_gradient(
      low = "firebrick",
      high = "steelblue",
      trans = "reverse"
    ) +
    scale_y_discrete(
      labels = function(x)
        stringr::str_wrap(x, width = 40)
    ) +
    labs(
      title = paste(
        "KEGG pathway enrichment:",
        component_name
      ),
      subtitle = "Top pathways ranked by adjusted p-value",
      x = "Gene count",
      y = NULL,
      color = "Adjusted p-value",
      size = "Gene count"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size = 10),
      panel.grid.minor = element_blank()
    )
}

# Generate plots
p_kegg_core <-
  plot_kegg_component(
    top_kegg,
    "core"
  )

p_kegg_shared <-
  plot_kegg_component(
    top_kegg,
    "shared"
  )

p_kegg_tissue <-
  plot_kegg_component(
    top_kegg,
    "tissue_specific_candidate"
  )

# Save Plot
if (!is.null(p_kegg_core)) {
  ggsave(
    file.path(PLOT_DIR, "KEGG_Core.png"),
    p_kegg_core,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

if (!is.null(p_kegg_shared)) {
  ggsave(
    file.path(PLOT_DIR, "KEGG_Shared.png"),
    p_kegg_shared,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

if (!is.null(p_kegg_tissue)) {
  ggsave(
    file.path(PLOT_DIR,
              "KEGG_Tissue_Specific.png"),
    p_kegg_tissue,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
}

# GSEA Analysis
# Map Ensembl gene IDs to Entrez IDs
gsea_genes <- meta_results |>
  dplyr::filter(
    !is.na(gene_id),
    !is.na(z_value),
    is.finite(z_value)
  ) |>
  dplyr::mutate(
    gene_id = sub("\\..*$", "", gene_id)
  ) |>
  dplyr::distinct(gene_id, .keep_all = TRUE)

gene_mapping <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(gsea_genes$gene_id),
  keytype = "ENSEMBL",
  columns = "ENTREZID"
)

# Keep mapped genes and rank by signed Z-score
gsea_ranked <- gsea_genes |>
  dplyr::inner_join(
    gene_mapping,
    by = c("gene_id" = "ENSEMBL")
  ) |>
  dplyr::filter(!is.na(ENTREZID)) |>
  dplyr::group_by(ENTREZID) |>
  dplyr::slice_max(
    order_by = abs(z_value),
    n = 1,
    with_ties = FALSE
  ) |>
  dplyr::ungroup()

gene_list <- gsea_ranked$z_value
names(gene_list) <- gsea_ranked$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# Check ranking
cat("Ranked Entrez genes:", length(gene_list), "\n")
cat("Unique Entrez IDs:", length(unique(names(gene_list))), "\n")

# Run KEGG GSEA
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  keyType = "ncbi-geneid",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  eps = 1e-10,
  verbose = TRUE
)

# Inspect significant pathways
gsea_results <- as.data.frame(gsea_kegg)

cat("Significant KEGG GSEA pathways:", nrow(gsea_results), "\n")

head(
  gsea_results[
    ,
    intersect(
      c("ID", "Description", "setSize", "enrichmentScore",
        "NES", "pvalue", "p.adjust", "qvalues"),
      names(gsea_results)
    )
  ]
)

# Save GSEA results
readr::write_csv(
  gsea_results,
  file.path(
    "outputs/meta_analysis/core_tissue_specific_decomposition/tables",
    "gsea_kegg_results.csv"
  )
)

# Plot top significant pathways
gsea_plot <- clusterProfiler::dotplot(
  gsea_kegg,
  showCategory = min(20, nrow(gsea_results)),
  x = "NES",
  color = "p.adjust"
) +
  ggplot2::labs(
    title = "KEGG GSEA of COVID-19 Meta-analysis",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  ggplot2::theme_bw()

ggsave(
  file.path(
    PLOT_DIR,
    "GSEA_KEGG_Dotplot.png"
  ),
  plot = gsea_plot,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)

# Additional GSEA visualizations
# Ridgeplot of top enriched pathways
gsea_ridge <- enrichplot::ridgeplot(
  gsea_kegg,
  showCategory = 20
) +
  ggplot2::labs(
    title = "KEGG GSEA enrichment landscape"
  ) +
  ggplot2::theme_bw()

ggsave(
  filename = file.path(
    PLOT_DIR,
    "GSEA_KEGG_Ridgeplot.png"
  ),
  plot = gsea_ridge,
  width = 12,
  height = 10,
  dpi = 300,
  bg = "white"
)

# Top pathway enrichment curves
# Neutrophil extracellular trap formation
p_net <- enrichplot::gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04613",
  title = "Neutrophil extracellular trap formation",
  pvalue_table = TRUE
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "GSEA_NET_Formation.png"
  ),
  plot = p_net,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# NOD-like receptor signaling pathway
p_nod <- enrichplot::gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04621",
  title = "NOD-like receptor signaling pathway",
  pvalue_table = TRUE
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "GSEA_NOD_Like_Receptor.png"
  ),
  plot = p_nod,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

# Combined figure
combined_gsea_plot <- patchwork::wrap_plots(
  p_net,
  p_nod,
  ncol = 1
)

ggsave(
  filename = file.path(
    PLOT_DIR,
    "GSEA_Top_Pathways.png"
  ),
  plot = combined_gsea_plot,
  width = 12,
  height = 12,
  dpi = 300,
  bg = "white"
)