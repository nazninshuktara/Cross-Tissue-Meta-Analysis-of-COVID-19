#!/usr/bin/env python3

# 03_salmon_mapping_stats
# Dataset: all discovery datasets

"""
Summarize Salmon quantification for the samples carried forward into the primary
analysis and build Supplementary Table S3 (one row per dataset).

Inputs
  data/metadata/metadata_curated_sample.csv   curated sample metadata (run_id = comma-separated runs)
  outputs/salmon_out/<dataset>/<run>/         Salmon output (aux_info/meta_info.json, lib_format_counts.json)
  outputs/qc/<dataset>/multiqc_report_data/   MultiQC logs (tool versions)

Outputs
  outputs/mapping_stats/salmon_run_stats.tsv           run-level statistics (all quantified runs)
  outputs/mapping_stats/Supplementary_Table_S3.csv     dataset-level summary
"""

import json
import re
from pathlib import Path

import pandas as pd

# Paths (script location: scripts/02_gene_level_counts/)
ROOT = Path(__file__).resolve().parents[2]
META_FILE = ROOT / "data" / "metadata" / "metadata_curated_sample.csv"
QUANT_DIR = ROOT / "outputs" / "salmon_out"
QC_DIR = ROOT / "outputs" / "qc"
OUT_DIR = ROOT / "outputs" / "mapping_stats" 

STATS_FILE = OUT_DIR / "salmon_run_stats.tsv"
SUMMARY_FILE = OUT_DIR / "mapping_summary.csv"

OUT_DIR.mkdir(parents=True, exist_ok=True)

# Step 1. Collect run-level Salmon statistics
rows = []
for meta_json in sorted(QUANT_DIR.glob("*/*/aux_info/meta_info.json")):
    run_dir = meta_json.parent.parent
    info = json.loads(meta_json.read_text())

    lib_json = run_dir / "lib_format_counts.json"
    lib_type = None
    if lib_json.exists():
        lib_type = json.loads(lib_json.read_text()).get("expected_format")

    rows.append({
        "dataset_id": run_dir.parent.name,
        "run_id": run_dir.name,
        "salmon_version": info.get("salmon_version"),
        "num_processed": info.get("num_processed"),
        "num_mapped": info.get("num_mapped"),
        "percent_mapped": info.get("percent_mapped"),
        "lib_type": lib_type,
    })

stats = pd.DataFrame(rows)
if stats.empty:
    raise SystemExit(f"No Salmon output found under {QUANT_DIR}")
stats.to_csv(STATS_FILE, sep="\t", index=False)
print(f"Runs quantified: {len(stats)} | saved: {STATS_FILE}")

# Step 2. Select included samples and expand run IDs (one row per run)
meta = pd.read_csv(META_FILE)
inc = meta[meta["inclusion_status"] == "included"].copy()
inc["run_id"] = inc["run_id"].astype(str).str.split(",")
runs = inc[["dataset_id", "sample_id", "run_id"]].explode("run_id")
runs["run_id"] = runs["run_id"].str.strip()

merged = runs.merge(stats, on=["dataset_id", "run_id"], how="left")
print(f"Included samples: {merged['sample_id'].nunique()} | "
      f"included runs: {len(merged)}")

missing = merged[merged["num_mapped"].isna()]
if len(missing) > 0:
    raise SystemExit(
        "Included runs without Salmon output:\n"
        + missing[["dataset_id", "run_id"]].to_string(index=False)
    )

# Step 3. Sample-level mapped reads (technical-replicate runs summed)
samples = (
    merged.groupby(["dataset_id", "sample_id"], as_index=False)[["num_processed", "num_mapped"]]
    .sum()
)
samples["mapping_rate"] = 100 * samples["num_mapped"] / samples["num_processed"]
samples["mapped_M"] = samples["num_mapped"] / 1e6

# Step 4. Tool versions from MultiQC logs
def multiqc_version(dataset):
    log = QC_DIR / dataset / "multiqc_report_data" / "multiqc.log"
    if log.exists():
        match = re.search(r"This is MultiQC v([0-9][0-9.]*)", log.read_text())
        if match:
            return match.group(1)
    return "NA"


def fastqc_version(dataset):
    path = QC_DIR / dataset / "multiqc_report_data" / "multiqc_software_versions.txt"

    if not path.exists():
        return "NA"

    text = path.read_text()

    match = re.search(
        r"FastQC\s+([0-9.]+)",
        text
    )

    return match.group(1) if match else "NA"

# Step 5. Dataset-level table
def count_summary(values):
    return "; ".join(f"{k} (n={v})" for k, v in values.value_counts().items())


by_run = merged.groupby("dataset_id")
by_sample = samples.groupby("dataset_id")

table = pd.DataFrame({
    "Samples": by_run["sample_id"].nunique(),
    "Runs": by_run["run_id"].nunique(),
    "Salmon version": by_run["salmon_version"].agg(count_summary),
    "Library type (detected)": by_run["lib_type"].agg(count_summary),
    "Mapped reads, median (M)": by_sample["mapped_M"].median().round(2),
    "Mapped reads, min (M)": by_sample["mapped_M"].min().round(3),
    "Mapping rate, median (%)": by_sample["mapping_rate"].median().round(1),
    "Mapping rate, min-max (%)": by_sample["mapping_rate"].agg(
        lambda x: f"{x.min():.1f}-{x.max():.1f}"
    ),
}).reset_index().rename(columns={"dataset_id": "Dataset"})

table.insert(1, "Tissue", table["Dataset"].map(inc.groupby("dataset_id")["tissue_type"].first()))
table["FastQC version"] = table["Dataset"].map(fastqc_version)
table["MultiQC version"] = table["Dataset"].map(multiqc_version)

# Step 6. Save
table.to_csv(OUT_DIR / "mapping_summary.csv", index=False)
print(f"Saved: {OUT_DIR / 'mapping_summary.csv'}")
print(table.to_string(index=False))
