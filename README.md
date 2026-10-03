# Cross-Tissue-Meta-Analysis-of-COVID-19
A two-tier **tissue-independent immune transcriptional signature of COVID-19** identified through cross-tissue random-effects meta-analysis of bulk RNA-seq datasets.

---

## Project Overview
COVID-19 affects multiple tissues and biological systems, making it challenging to distinguish disease-associated transcriptional responses from tissue-specific background variation. This project investigates whether COVID-19 induces a **tissue-independent transcriptional immune response** that can be consistently detected across different human tissues and biological sample types.

We perform a **cross-tissue random-effects meta-analysis of bulk RNA-seq datasets** from COVID-19 cases and controls to identify a robust tissue-independent immune core signature.

The resulting signature is further evaluated through comparison with influenza, cell-type deconvolution, diagnostic classification, host genetics, tissue-specific expression, protein–protein interaction analysis, transcription-factor activity, perturbational signatures, single-cell references, and independent cohort validation.

The study ultimately defines a **two-tier core structure** consisting of:

- COVID-19-specific transcriptional components
- Shared respiratory-inflammatory components

---

## Study Rationale
COVID-19 produces transcriptional changes across multiple tissues and biological compartments. However, tissue-specific expression patterns and differences in sample composition can make it difficult to identify signals that are consistently associated with disease. A cross-tissue meta-analysis provides a framework for identifying transcriptional changes that are reproducible across heterogeneous datasets while accounting for between-study variation.

This study therefore combines dataset-level differential-expression analysis with random-effects meta-analysis to identify a **tissue-independent immune core**, followed by influenza comparison and independent validation.

---

## Datasets
### Discovery Datasets
Nine GEO datasets spanning eight tissue types were included in the discovery meta-analysis, comprising 307 biological samples (155 COVID-19 and 66 controls).

| **GEO accession** | **Tissue type** | **Sample size** |
|---|---|---|
| GSE208342 | Nasopharyngeal | 20 |
| GSE273149 | PBMC | 49 |
| GSE243217 | Peripheral blood | 72 |
| GSE202182 | Kidney | 37 |
| GSE208076 | Lung | 10 |
| GSE211979 | PBMC | 50 |
| GSE169241 | Heart | 23 |
| GSE149312 | Gut organoids | 22 |
| GSE182917 | Liver, Kidney, Lung, Heart, Spleen | 24 |

> **Data processing:** Raw sequencing data were obtained from the Sequence Read Archive (SRA). Technical replicates were identified from study metadata and collapsed at the appropriate level before differential-expression analysis. Dataset-specific differential expression was performed independently before cross-tissue meta-analysis.

### Independent Validation Datasets
6 independent GEO datasets are used to evaluate the reproducibility of the core signature:

| GEO accession |
|---|
| GSE261002 |
| GSE245922 |
| GSE251849 |
| GSE233557 |
| GSE244488 |
| GSE236841 |

### Influenza Comparator
The influenza dataset(GSE162632) is used to distinguish COVID-19-associated transcriptional responses from responses shared with another respiratory viral infection.

- **89 influenza cases**
- **89 mock controls**
- PBMC samples

---

## External Resources

The project integrates several publicly available resources:

- **GTEx v8** — tissue-specific gene expression
- **COVID-19 Host Genetics Initiative, Freeze 7** — host genetic associations
- **CZ CELLxGENE Census** — single-cell expression reference
- **LINCS L1000 / CLUE** — perturbational signatures and compound analysis
- **STRING** — protein–protein interaction analysis
- **DoRothEA** — transcription-factor/regulatory activity

---

## Research Questions

1. Can a tissue-independent COVID-19 immune core be identified?
2. Which genes are consistently altered across tissues?
3. Which components are COVID-19-specific versus shared with influenza?
4. Is the core signature reproducible in independent cohorts?
5. Can the signature support COVID-19 classification?
6. What biological, genetic, cellular, and regulatory mechanisms support the signature?

---

## Hypotheses

- COVID-19 contains a reproducible tissue-independent transcriptional immune response.
- A subset of genes will remain consistent across tissues and datasets.
- The core will contain both COVID-19-specific and shared respiratory-inflammatory components.
- The signature will show independent support from cell-type, genetic, tissue-expression, and single-cell analyses.

---

# Project Pipeline
The analysis is organized as a numbered and reproducible computational pipeline.

```text
                 ┌─────────────────────────┐
                 │   Public GEO / SRA Data │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ Data Retrieval & QC     │
                 │ GEO / SRA / metadata    │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ Salmon quant.sf files   │
                 │Transcript quantification│
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │       tximport          │
                 │ Transcript → Gene       │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ Gene-level expression   │
                 │ matrix                  │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ Technical Replicate     │
                 │ Collapsing              │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ QC / Filtering /        │
                 │ Normalization            │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │       DESeq2            │
                 │ Per-dataset DE analysis │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ Random-effects          │
                 │ Meta-analysis (metafor) │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ 1,532 Meta-significant  │
                 │ Genes                   │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ Tissue-independent      │
                 │ Core Signature           │
                 └────────────┬────────────┘
                              │
             ┌────────────────┼─────────────────┐
             │                │                 │
             ▼                ▼                 ▼
      Leave-one-tissue   Meta-regression   Publication Bias
          analysis                           / Bayesian
             │                │                 │
             └────────────────┼─────────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ Influenza Comparison    │
                 │ GSE162632               │
                 └────────────┬────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │      Two-tier Core      │
                 │                         │
                 │ COVID-specific          │
                 │ Shared respiratory      │
                 │ inflammatory            │
                 └────────────┬────────────┘
                              │
             ┌────────────────┼────────────────────┐
             │                │                    │
             ▼                ▼                    ▼
      Cell-type          Diagnostic           Independent
      Deconvolution      Classifier            Validation
             │                │                    │
             │                ▼                    │
             │          Elastic Net                │
             │                │                    │
             └────────────────┼────────────────────┘
                              │
                              ▼
                 ┌─────────────────────────┐
                 │ Orthogonal Validation   │
                 ├─────────────────────────┤
                 │ GWAS                    │
                 │ GTEx                    │
                 │ PPI / STRING            │
                 │ DoRothEA                │
                 │ LINCS / CLUE            │
                 │ CELLxGENE Census        │
                 └─────────────────────────┘
```
---

## Project Structure
```
Cross-Tissue-Meta-Analysis-of-COVID-19/
│
├── pixi.toml
├── pixi.lock
├── README.md
├── LICENSE
├── CITATION.cff
│
├── scripts/
│   ├── 00_*.R
│   ├── 01_download_geo_sra.sh
│   ├── 02_*.R
│   ├── 03_*.R
│   ├── 04_*.R
│   ├── 05_*.R
│   ├── 06_meta_analysis.R
│   ├── 07_*.R
│   ├── 08_*.R
│   ├── ...
│   ├── 18_*.R
│   ├── 19_*.R
│   ├── 20_*.R
│   │
│   └── SCRIPT_MANIFEST.md
│
├── data/
│   ├── metadata/
│   │   ├── discovery/
│   │   ├── validation/
│   │   └── influenza/
│   │
│   ├── raw/
│   │   ├── sra/
│   │   └── fastq/
│   │
│   └── processed/
│       ├── salmon/
│       ├── tximport/
│       ├── collapsed/
│       └── counts/
│
├── outputs/
│   ├── qc/
│   └── salmon_out/
│
├── results/
│   ├── tables/
│   │   ├── differential_expression/
│   │   ├── meta_analysis/
│   │   ├── core_signature/
│   │   ├── validation/
│   │   └── classifier/
│   │
│   └── figures/
│
├── inputs/
│   ├── tables/
│   └── figures/
│
└── docs/
```
---

## Reproducibility
The computational environment is managed using Pixi.

The repository records the software environment through:

- `pixi.toml`
- `pixi.lock`

The analysis is implemented as numbered scripts in `scripts/`, with metadata, intermediate data, quality-control outputs, results, and supplementary materials organized according to the project structure.

---
