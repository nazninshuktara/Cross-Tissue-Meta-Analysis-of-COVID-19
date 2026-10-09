#!/usr/bin/env python3

"""
Metadata curation for the COVID-19 cross-tissue meta-analysis.

Input:
    inputs/discovery/*_SraRunTable.csv

Output:
    data/metadata/metadata_run_level.csv
    data/metadata/metadata_sample_level.csv
    data/metadata/metadata_curated_sample.csv
    data/metadata/discovery_dataset_summary.csv

Run-level:
    One row per SRA run.

Sample-level:
    One row per biological sample/GSM.
    Technical replicate SRRs are collapsed into one sample.

Curated sample-level:
    Sample-level metadata with inclusion/exclusion status.

Dataset-specific notes:
    GSE273149: Day 1 (baseline) samples only.
    GSE243217: sepsis samples were excluded from the primary analysis.
    GSE202182: 4 hantavirus and 10 non-COVID acute tubular injury samples were excluded from the primary analysis.
    GSE182917: lung samples only. Heart, kidney, liver and spleen samples
               are COVID-19 cases without matched controls, so they are
               excluded from the primary analysis.
"""

import os
import re
import pandas as pd

# Define input and output directories.
INPUT_DIR = "inputs/discovery"
OUTPUT_DIR = "data/metadata"
# Use a consistent label for missing metadata.
NA = "NA"

# Discovery datasets used in the analysis
DATASETS = [
    "GSE208342",
    "GSE273149",
    "GSE243217",
    "GSE202182",
    "GSE208076",
    "GSE211979",
    "GSE169241",
    "GSE149312",
    "GSE182917",
]
# Shared metadata columns
BASE_COLUMNS = [
    "dataset_id",
    "run_id",
    "sample_id",
    "tissue_type",
    "disease_status",
    "disease_group",
    "disease_severity",
    "treatment",
    "outcome",
    "timepoint",
    "age",
    "sex",
    "comorbidities",
    "vaccination_status",
    "infection_variant",
    "sequencing_platform",
    "library_type",
    "data_source",
    "metadata_confidence",
]
# Additional columns for curated metadata
CURATED_COLUMNS = BASE_COLUMNS + [
    "inclusion_status",
    "exclusion_reason"
]


def get_value(row, column, default=NA):
    """Return a cleaned metadata value."""
    if column not in row.index:
        return default

    value = row[column]

    if pd.isna(value):
        return default

    value = str(value).strip()

    if value.lower() in {
        "",
        "missing",
        "nan",
        "none",
        "not reported",
        "na",
        "n/a",
        "null",
    }:
        return default

    return value


def get_first_value(row, columns, default=NA):
    """Return the first available metadata value."""
    for column in columns:
        value = get_value(row, column)
        if value != NA:
            return value

    return default


def normalize_age(value):
    """Normalize age values."""
    if value == NA:
        return NA

    value = str(value).strip()

    match = re.match(
        r"^(\d+(?:\.\d+)?)\s*(?:Y|y|years?|yrs?)?$",
        value,
    )

    if match:
        return match.group(1)

    return value


def normalize_sex(value):
    """Normalize sex labels."""
    if value == NA:
        return NA

    value = str(value).strip().lower()

    if value in ("male", "m"):
        return "Male"

    if value in ("female", "f"):
        return "Female"

    return value


def normalize_tissue(gse, row):
    """Normalize tissue and cell-type labels."""
    tissue = get_first_value(
        row,
        (
            "tissue",
            "tissue/cell_type",
            "source_name",
        ),
    )

    if tissue == NA:
        return NA

    tissue_lower = tissue.lower()

    if "nasophary" in tissue_lower:
        return "Nasopharyngeal"

    if "pbmc" in tissue_lower:
        return "PBMC"

    if "peripheral blood" in tissue_lower:
        return "Peripheral blood"

    if tissue_lower == "blood":
        return "Blood"

    if "kidney" in tissue_lower:
        return "Kidney"

    if "lung" in tissue_lower:
        return "Lung"

    if (
        "gut" in tissue_lower
        or "intestinal" in tissue_lower
        or "organoid" in tissue_lower
    ):
        return "Gut organoids"

    if gse == "GSE169241":
        if tissue_lower == "human heart":
            return "Heart"

        if tissue_lower == "human es-derived cardiomyocyte":
            return "Human ES-derived cardiomyocyte"

        if tissue_lower == "human es-derived macrophage":
            return "Human ES-derived macrophage"

    return tissue

def normalize_disease(gse, row):
    """Keep original disease status labels unchanged."""

    if gse == "GSE169241":
        return get_first_value(
            row,
            (
                "source_name",
                "virus",
                "disease_state",
                "disease",
                "disease_status",
                "condition",
                "status",
            ),
        )

    return get_first_value(
        row,
        (
            "disease_state",
            "disease",
            "disease_status",
            "condition",
            "status",
        ),
    )

def assign_disease_group(dataset_id, disease_status):
    """Assign standardized disease categories."""

    disease = str(disease_status).strip().lower()

    if dataset_id == "GSE208342":
        if disease == "covid-19":
            return "COVID-19"
        if disease == "negative control":
            return "Control"

    elif dataset_id == "GSE273149":
        if disease in {
            "covid_nonsurvivor",
            "covid_survivor",
        }:
            return "COVID-19"
        if disease == "ards_con":
            return "Control"

    elif dataset_id == "GSE243217":
        if disease == "covid-19":
            return "COVID-19"
        if disease == "healthy donor":
            return "Control"
        if "sepsis" in disease:
            return "Sepsis"

    elif dataset_id == "GSE202182":
        if disease == "covid-19":
            return "COVID-19"
        if disease == "healthy kidney":
            return "Control"
        if disease == "non-covid-19 ati kidneys":
            return "Non-COVID-19 ATI"
        if "hantavirus" in disease:
            return "Hantavirus"

    elif dataset_id == "GSE208076":
        if disease == "covid-19":
            return "COVID-19"
        if disease == "normal":
            return "Control"

    elif dataset_id == "GSE211979":
        if disease in {
            "covid19 icu",
            "covid19 nonicu",
        }:
            return "COVID-19"
        
        if disease == "healthy":
            return "Control"

    elif dataset_id == "GSE169241":
        if disease == "sars-cov-2 infected_heart":
            return "COVID-19"
        if disease == "uninfected_heart":
            return "Control"

    elif dataset_id == "GSE149312":
        if disease == "sars-cov-2 infected":
            return "COVID-19"
        if disease == "sars-cov infected":
            return "SARS-COV"
        if disease == "control (untreated)":
            return "Control"

    elif dataset_id == "GSE182917":
        if disease in {
            "covid-19",
            "covid-19 (severe pneumonia)",
            "covid-19 (critical)",
        }:
            return "COVID-19"

        if disease in {
            "lung squamous cell carcinoma",
            "pneumonia",
            "pulmonary bronchiectasis",
        }:
            return "Control"

    return NA

def normalize_severity(gse, row):
    """Extract and standardize disease severity."""

    # GSE208342
    if gse == "GSE208342":
        return get_value(row, "level_of_severity")


    # GSE182917
    if gse == "GSE182917":
        disease_state = get_value(row, "disease_state")
        
        if disease_state != NA:
            disease_lower = disease_state.lower()
            
            if "critical" in disease_lower:
                return "Critical"

            if "severe pneumonia" in disease_lower:
                return "Severe"

            if "severe" in disease_lower:
                return "Severe"
            
            return NA


    # GSE211979
    if gse == "GSE211979":
        disease = get_value(row, "disease_state")
        
        if disease == NA:
            return NA

        disease_lower = disease.lower()

        if disease_lower == "covid19 icu":
            return "ICU"

        if disease_lower == "covid19 nonicu":
            return "Non-ICU"

        return NA


    # Other datasets with explicit severity columns
    severity = get_first_value(
        row,
        (
            "disease_severity",
            "severity",
            "severity_status",
        ),
    )


    if severity != NA:
        severity_lower = severity.lower()

        if "critical" in severity_lower:
            return "Critical"

        if "severe" in severity_lower:
            return "Severe"

        if "moderate" in severity_lower:
            return "Moderate"

        if "mild" in severity_lower:
            return "Mild"

        if "asymptomatic" in severity_lower:
            return "Asymptomatic"


    return severity


def get_timepoint(gse, row):
    """Extract study timepoint information."""
    return get_first_value(
        row,
        (
            "timepoint",
            "time_point",
            "time",
            "time_post_infection",
            "days_post_infection",
            "collection_time",
        ),
    )


def get_outcome(gse, row):
    """Extract clinical outcome information."""
    if gse == "GSE273149":
        disease_status = normalize_disease(gse, row)

        if disease_status != NA:
            disease_lower = disease_status.lower()

            if "nonsurvivor" in disease_lower:
                return "died"

            if "survivor" in disease_lower:
                return "survived"

    return get_first_value(
        row,
        (
            "outcome",
            "clinical_outcome",
            "patient_outcome",
            "survival",
        ),
    )


def get_treatment(row):
    """Extract treatment information."""
    return get_first_value(
        row,
        (
            "treatment",
            "treatment_ch1",
            "treatment_group",
            "drug",
            "therapy",
            "stimulus",
        ),
    )


def get_age(row):
    """Extract and normalize age."""
    age = get_first_value(
        row,
        (
            "AGE",
            "age",
            "age_ch1",
            "age_years",
            "patient_age",
            "donor_age",
        ),
    )

    return normalize_age(age)


def get_sex(row):
    """Extract and normalize sex."""
    sex = get_first_value(
        row,
        (
            "sex",
            "sex_ch1",
            "sex_assigned_at_birth",
            "gender",
            "gender_ch1",
        ),
    )

    return normalize_sex(sex)


def get_comorbidities(row):
    """Extract comorbidity information."""
    return get_first_value(
        row,
        (
            "comorbidities",
            "comorbidity",
            "comorbidities_ch1",
            "co-morbidities",
            "medical_history",
        ),
    )


def get_vaccination_status(row):
    """Extract vaccination information."""
    return get_first_value(
        row,
        (
            "vaccination_status",
            "vaccination",
            "vaccination_status_ch1",
            "vaccinated",
            "vaccine",
        ),
    )


def get_infection_variant(row):
    """Extract infection variant information."""
    return get_first_value(
        row,
        (
            "infection_variant",
            "variant",
            "virus_variant",
            "strain",
            "viral_strain",
        ),
    )


def get_platform(row):
    """Extract sequencing platform information."""
    return get_first_value(
        row,
        (
            "Platform",
            "platform",
            "sequencing_platform",
            "instrument_model",
        ),
    )


def get_library_type(row):
    """Combine library strategy, layout, and selection."""
    values = []

    for columns in [
        ("LibraryStrategy", "library_strategy"),
        ("LibraryLayout", "library_layout"),
        ("LibrarySelection", "library_selection"),
    ]:
        value = get_first_value(row, columns)

        if value != NA:
            values.append(value)

    if not values:
        return NA

    return "; ".join(dict.fromkeys(values))


def map_row(gse, row):
    """Map one SRA run to the harmonized metadata schema."""
    sample_id = get_value(row, "Sample Name")

    if sample_id == NA:
        sample_id = get_value(row, "BioSample")

    run_id = get_value(row, "Run")

    disease_status = normalize_disease(gse, row)
    
    disease_group = assign_disease_group(gse, disease_status)

    disease_severity = normalize_severity(gse, row)

    metadata_confidence = "high"

    if gse == "GSE182917":
        metadata_confidence = "medium"

    return {
        "dataset_id": gse,
        "run_id": run_id,
        "sample_id": sample_id,
        "tissue_type": normalize_tissue(gse, row),
        "disease_status": disease_status,
        "disease_group": disease_group,
        "disease_severity": disease_severity,
        "treatment": get_treatment(row),
        "outcome": get_outcome(gse, row),
        "timepoint": get_timepoint(gse, row),
        "age": get_age(row),
        "sex": get_sex(row),
        "comorbidities": get_comorbidities(row),
        "vaccination_status": get_vaccination_status(row),
        "infection_variant": get_infection_variant(row),
        "sequencing_platform": get_platform(row),
        "library_type": get_library_type(row),
        "data_source": "NCBI-GEO",
        "metadata_confidence": metadata_confidence,
    }


def reconcile_values(values):
    """Reconcile metadata values across technical replicates."""
    cleaned = []

    for value in values:
        if pd.isna(value):
            continue

        value = str(value).strip()

        if value in ("", NA):
            continue

        if value not in cleaned:
            cleaned.append(value)

    if not cleaned:
        return NA

    return "; ".join(cleaned)


def collapse_dataset(run_metadata):
    """Collapse technical replicate SRRs into biological samples."""
    collapsed_records = []

    grouped = run_metadata.groupby(
        ["dataset_id", "sample_id"],
        sort=False,
        dropna=False,
        observed=True,
    )

    for (dataset_id, sample_id), group in grouped:
        record = {
            "dataset_id": dataset_id,
            "sample_id": sample_id,
        }

        run_ids = group["run_id"].tolist()

        record["run_id"] = ",".join(run_ids)

        for column in BASE_COLUMNS:
            if column in {
                "dataset_id",
                "run_id",
                "sample_id",
            }:
                continue

            record[column] = reconcile_values(
                group[column].tolist()
            )

        collapsed_records.append(record)

    return pd.DataFrame(
        collapsed_records
    )[BASE_COLUMNS]


def apply_analysis_criteria(sample_metadata):
    """Apply dataset-specific primary-analysis criteria."""
    curated = sample_metadata.copy()

    curated["inclusion_status"] = "included"
    curated["exclusion_reason"] = NA


    for index, row in curated.iterrows():
        gse = row["dataset_id"]
        disease = str(row["disease_group"]).lower()
        timepoint = str(row["timepoint"]).lower()
        tissue = str(row["tissue_type"]).lower()

        if gse == "GSE208342":
            if disease in {
                "covid-19",
                "control",
            }:
                curated.at[index, "inclusion_status"] = "included"
            else:
                curated.at[index, "inclusion_status"] = "excluded"
                curated.at[
                    index,
                    "exclusion_reason",
                ] = "Not COVID-19 or negative control"

        elif gse == "GSE273149":
            is_day_1 = bool(
                re.fullmatch(
                    r"(day\s*0?1|d\s*0?1|baseline)",
                    timepoint,
                )
            )

            if (
                is_day_1
                and disease in {
                    "covid-19",
                    "control",
                }
            ):
                curated.at[index, "inclusion_status"] = "included"
            else:
                curated.at[index, "inclusion_status"] = "excluded"

                if not is_day_1:
                    curated.at[
                        index,
                        "exclusion_reason",
                    ] = "Non-Day 1 sample"
                else:
                    curated.at[
                        index,
                        "exclusion_reason",
                    ] = "Not COVID-19 or ARDS_CON"

        elif gse == "GSE243217":
            if "sepsis" in disease:
                curated.at[
                    index,
                    "inclusion_status",
                ] = "specificity_analysis"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = "Sepsis sample excluded from primary analysis"

            elif disease in {
                "covid-19",
                "control",
            }:
                curated.at[index, "inclusion_status"] = "included"

            else:
                curated.at[index, "inclusion_status"] = "excluded"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = "Not COVID-19 or healthy donor"

        elif gse == "GSE202182":
            if "hantavirus" in disease:
                curated.at[
                    index,
                    "inclusion_status",
                ] = "excluded"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = "Hantavirus infection"

            elif (
                "non-covid" in disease
                or "non covid" in disease
                or "acute tubular injury" in disease
            ):
                curated.at[
                    index,
                    "inclusion_status",
                ] = "excluded"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = "Non-COVID acute tubular injury"

            elif disease in {
                "covid-19",
                "control",
            }:
                curated.at[index, "inclusion_status"] = "included"

            else:
                curated.at[index, "inclusion_status"] = "excluded"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = "Not COVID-19 or healthy kidney"

        elif gse == "GSE208076":
            if disease in {
                "covid-19",
                "control",
            }:
                curated.at[index, "inclusion_status"] = "included"
            else:
                curated.at[index, "inclusion_status"] = "excluded"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = "Not COVID-19 or normal lung"

        elif gse == "GSE211979":
            if disease in {
                "covid-19",
                "control",
            }:
                curated.at[index, "inclusion_status"] = "included"
            else:
                curated.at[index, "inclusion_status"] = "excluded"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = "Not COVID-19 or healthy donor"

        elif gse == "GSE169241":
            if (
                tissue == "heart"
                and disease in {
                    "covid-19",
                    "control",
                }
            ):
                curated.at[index, "inclusion_status"] = "included"
            else:
                curated.at[index, "inclusion_status"] = "excluded"

                if tissue != "heart":
                    reason = (
                        "ES-derived model rather than human heart"
                    )
                else:
                    reason = "Not COVID-19 or control heart"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = reason

        elif gse == "GSE149312":
            if disease in {
                "covid-19",
                "control",
            }:
                curated.at[index, "inclusion_status"] = "included"

            elif "sars-cov" in disease:
                curated.at[
                    index,
                    "inclusion_status",
                ] = "excluded"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = "SARS-CoV comparator"

            else:
                curated.at[index, "inclusion_status"] = "excluded"

                curated.at[
                    index,
                    "exclusion_reason",
                ] = "Not COVID-19"

        elif gse == "GSE182917":
            # Lung only: heart, kidney, liver and spleen samples are
            # COVID-19 cases with no matched controls, so organ and
            # disease are confounded and DE cannot be estimated.
            is_cov_ctrl = disease in {
                "covid-19",
                "control",
            }
 
            if tissue == "lung" and is_cov_ctrl:
                curated.at[index, "inclusion_status"] = "included"
            else:
                curated.at[index, "inclusion_status"] = "excluded"
 
                if is_cov_ctrl and tissue != "lung":
                    reason = "Non-lung organ without matched control"
                else:
                    reason = "Not COVID-19 or control lung"
 
                curated.at[
                    index,
                    "exclusion_reason",
                ] = reason
 
    return curated[CURATED_COLUMNS]

def validate_metadata(
    run_metadata,
    sample_metadata,
    curated_metadata,
    input_run_counts,
):
    """Validate metadata structure and run accounting."""

    observed_run_counts = (
        run_metadata["dataset_id"]
        .value_counts()
        .reindex(DATASETS, fill_value=0)
        .to_dict()
    )

    assert observed_run_counts == input_run_counts, (
        "Run-level counts do not match input files."
    )

    unique_samples = (
        sample_metadata[
            [
                "dataset_id",
                "sample_id",
            ]
        ]
        .drop_duplicates()
        .shape[0]
    )

    assert len(sample_metadata) == unique_samples
    assert len(curated_metadata) == len(sample_metadata)

    assert list(run_metadata.columns) == BASE_COLUMNS
    assert list(sample_metadata.columns) == BASE_COLUMNS
    assert list(curated_metadata.columns) == CURATED_COLUMNS

    assert (
        run_metadata["dataset_id"]
        .drop_duplicates()
        .tolist()
        == DATASETS
    )

    assert not run_metadata.isna().any().any()
    assert not sample_metadata.isna().any().any()
    assert not curated_metadata.isna().any().any()

    missing_run_disease = run_metadata[
        run_metadata["disease_group"] == NA
    ][
        [
            "dataset_id",
            "run_id",
            "sample_id",
            "disease_status",
            "disease_group",
        ]
    ]

    if not missing_run_disease.empty:
        print("\nRuns with missing disease_group:")
        print(missing_run_disease.to_string(index=False))

    missing_sample_disease = sample_metadata[
        sample_metadata["disease_group"] == NA
    ][
        [
            "dataset_id",
            "sample_id",
            "disease_status",
            "disease_group",
        ]
    ]

    if not missing_sample_disease.empty:
        print("\nSamples with missing disease_group:")
        print(missing_sample_disease.to_string(index=False))

    assert run_metadata["run_id"].is_unique
    assert (run_metadata["run_id"] != NA).all()

    assert (run_metadata["sample_id"] != NA).all()
    assert (sample_metadata["sample_id"] != NA).all()

    assert not sample_metadata[
        [
            "dataset_id",
            "sample_id",
        ]
    ].duplicated().any()

    run_ids = run_metadata["run_id"].tolist()

    collapsed_run_ids = []

    for value in sample_metadata["run_id"]:
        collapsed_run_ids.extend(value.split(","))

    assert len(run_ids) == len(collapsed_run_ids)
    assert set(run_ids) == set(collapsed_run_ids)
    assert len(set(collapsed_run_ids)) == len(collapsed_run_ids)

    valid_status = {
        "included",
        "excluded",
        "specificity_analysis",
    }

    assert set(
        curated_metadata["inclusion_status"]
    ).issubset(valid_status)

    excluded = (
        curated_metadata["inclusion_status"]
        == "excluded"
    )

    assert (
        curated_metadata.loc[
            excluded,
            "exclusion_reason",
        ]
        != NA
    ).all()

    specificity = (
        curated_metadata["inclusion_status"]
        == "specificity_analysis"
    )

    assert (
        curated_metadata.loc[
            specificity,
            "exclusion_reason",
        ]
        != NA
    ).all()

    included = (
        curated_metadata["inclusion_status"]
        == "included"
    )

    assert (
        curated_metadata.loc[
            included,
            "exclusion_reason",
        ]
        == NA
    ).all()

    included_missing_group = curated_metadata[
        (
            curated_metadata["inclusion_status"] == "included"
        ) & (
            curated_metadata["disease_group"] == NA
        )
    ]

    assert included_missing_group.empty, (
        "Included samples have missing disease_group."
    )

    print(f"Run-level rows: {len(run_metadata)}")
    print(f"Sample-level rows: {len(sample_metadata)}")
    print(f"Curated rows: {len(curated_metadata)}")
    print("Column validation: passed")
    print("Input run-count validation: passed")
    print("NA validation: passed")
    print("Run accounting: passed")
    print("Sample uniqueness: passed")
    print("Inclusion/exclusion validation: passed")


def main():
    """Build, validate, and save metadata tables."""

    os.makedirs(
        OUTPUT_DIR,
        exist_ok=True,
    )

    run_frames = []
    input_run_counts = {}

    for gse in DATASETS:
        input_file = os.path.join(
            INPUT_DIR,
            f"{gse}_SraRunTable.csv",
        )

        if not os.path.exists(input_file):
            raise FileNotFoundError(
                f"Metadata file not found: {input_file}"
            )

        df = pd.read_csv(
            input_file,
            dtype=str,
            keep_default_na=False,
        )

        input_run_counts[gse] = len(df)

        rows = [
            map_row(gse, row)
            for _, row in df.iterrows()
        ]

        run_metadata = pd.DataFrame(
            rows
        )[BASE_COLUMNS]

        run_frames.append(run_metadata)

    metadata_run = pd.concat(
        run_frames,
        ignore_index=True,
    )

    metadata_run["dataset_id"] = pd.Categorical(
        metadata_run["dataset_id"],
        categories=DATASETS,
        ordered=True,
    )

    metadata_run = (
        metadata_run
        .sort_values(
            [
                "dataset_id",
                "sample_id",
                "run_id",
            ]
        )
        .reset_index(drop=True)
    )

    metadata_run["dataset_id"] = (
        metadata_run["dataset_id"].astype(str)
    )

    metadata_sample = collapse_dataset(
        metadata_run
    )

    metadata_sample["dataset_id"] = pd.Categorical(
        metadata_sample["dataset_id"],
        categories=DATASETS,
        ordered=True,
    )

    metadata_sample = (
        metadata_sample
        .sort_values(
            [
                "dataset_id",
                "sample_id",
            ]
        )
        .reset_index(drop=True)
    )

    metadata_sample["dataset_id"] = (
        metadata_sample["dataset_id"].astype(str)
    )

    metadata_curated = apply_analysis_criteria(
        metadata_sample
    )

    validate_metadata(
        metadata_run,
        metadata_sample,
        metadata_curated,
        input_run_counts,
    )

    run_output = os.path.join(
        OUTPUT_DIR,
        "metadata_run_level.csv",
    )

    sample_output = os.path.join(
        OUTPUT_DIR,
        "metadata_sample_level.csv",
    )

    curated_output = os.path.join(
        OUTPUT_DIR,
        "metadata_curated_sample.csv",
    )

    summary_output = os.path.join(
    OUTPUT_DIR,
    "discovery_dataset_summary.csv",
    )

    metadata_run.to_csv(
        run_output,
        index=False,
    )

    metadata_sample.to_csv(
        sample_output,
        index=False,
    )

    metadata_curated.to_csv(
        curated_output,
        index=False,
    )

    discovery_dataset_summary = []

    for gse in DATASETS:
        run_count = (
            metadata_run["dataset_id"] == gse
        ).sum()

        sample_count = (
            metadata_sample["dataset_id"] == gse
        ).sum()

        curated_data = metadata_curated[
            metadata_curated["dataset_id"] == gse
        ]
        
        primary_data = curated_data[
            curated_data["inclusion_status"] == "included"
        ]

        discovery_dataset_summary.append(
            {
                "dataset_id": gse,
                "run_count": run_count,
                "sample_count": sample_count,
                "included_count": (
                    curated_data["inclusion_status"]
                    == "included"
                ).sum(),
                "excluded_count": (
                    curated_data["inclusion_status"]
                    == "excluded"
                ).sum(),
                "specificity_count": (
                    curated_data["inclusion_status"]
                    == "specificity_analysis"
                ).sum(),
                "covid_count": (
                    primary_data["disease_group"]
                    == "COVID-19"
                ).sum(),
                "control_count": (
                    primary_data["disease_group"]
                    == "Control"
                ).sum(),
            }
        )
        
    discovery_dataset_summary = pd.DataFrame(
        discovery_dataset_summary
    )

    total_row = {
        "dataset_id": "TOTAL",
        "run_count": discovery_dataset_summary["run_count"].sum(),
        "sample_count": discovery_dataset_summary["sample_count"].sum(),
        "included_count": discovery_dataset_summary["included_count"].sum(),
        "excluded_count": discovery_dataset_summary["excluded_count"].sum(),
        "specificity_count": discovery_dataset_summary["specificity_count"].sum(),
        "covid_count": discovery_dataset_summary["covid_count"].sum(),
        "control_count": discovery_dataset_summary["control_count"].sum(),
    }

    discovery_dataset_summary = pd.concat(
        [
            discovery_dataset_summary,
            pd.DataFrame([total_row]),
        ],
        ignore_index=True,
    )

    discovery_dataset_summary.to_csv(
        summary_output,
        index=False,
    )


    print(
        "\nPer-dataset run-level -> "
        "sample-level -> curated counts:"
    )

    for gse in DATASETS:
        run_count = (
            metadata_run["dataset_id"] == gse
        ).sum()

        sample_count = (
            metadata_sample["dataset_id"] == gse
        ).sum()

        curated_count = (
            metadata_curated["dataset_id"] == gse
        ).sum()

        print(
            f"{gse}: "
            f"{run_count} -> "
            f"{sample_count} -> "
            f"{curated_count}"
        )

    print("\nCurated inclusion status:")

    print(
        metadata_curated[
            "inclusion_status"
        ].value_counts()
    )

    primary = metadata_curated[
        metadata_curated["inclusion_status"]
        == "included"
    ]

    print(
        "\nPrimary analysis disease-status counts:"
    )

    print(
        primary[
            "disease_status"
        ].value_counts()
    )

    print(
        "\nPer-dataset primary disease-status counts:"
    )

    for gse in DATASETS:
        dataset = primary[
            primary["dataset_id"] == gse
        ]

        counts = (
            dataset["disease_status"]
            .value_counts()
            .to_dict()
        )

        print(
            f"{gse}: "
            f"{counts}"
        )

    print(
        "\nPer-dataset COVID-19 and control counts:"
    )
    
    for gse in DATASETS:
        subset = metadata_curated[
            (
                metadata_curated["dataset_id"]
                == gse
            )
            & (
                metadata_curated["inclusion_status"]
                == "included"
            )
        ]

        covid_count = (
            subset["disease_group"]
            == "COVID-19"
        ).sum()
        
        control_count = (
            subset["disease_group"]
            == "Control"
        ).sum()

        print(
            f"{gse}: "
            f"COVID-19 = {covid_count}, "
            f"Control = {control_count}"
        )

    print("\nOutput files:")
    print(run_output)
    print(sample_output)
    print(curated_output)
    print(summary_output)
    print("\nMetadata curation completed.")


if __name__ == "__main__":
    main()