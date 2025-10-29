# Example—Bulk RNA-seq (DESeq2)

This repository provides a simple example of performing **differential gene expression analysis** between **Alzheimer’s disease (AD) cases and controls** using the DESeq2 package in R.

---

## Overview

The included R script (`scripts/Example-DESeq2.R`) demonstrates how to:

- Import and merge assay, biospecimen, and clinical metadata
- Match RNA-seq expression counts to sample metadata
- Filter and clean metadata for complete cases
- Build a DESeq2 dataset with key covariates
- Run differential expression for AD vs. control
- Export results as ranked tables

---

## Inputs

- **Raw counts matrix** (genes × samples, integer counts)
- **Metadata files**, including:
  - Assay metadata (sample-level RNA-seq info)
  - Biospecimen metadata (links specimenID ↔ individualID)
  - Clinical metadata (AD status, age at death, sex, RIN, etc.)

> Replace file paths in the script with your own local dataset paths.

---

## Output

- Differential expression results table (`.csv` or `.tsv`) with:
  - `gene`, `log2FoldChange`, `pvalue`, and `padj`
- Results are sorted by adjusted p-value (FDR-corrected)

---

## Requirements

- **R ≥ 4.1**
- R packages: `DESeq2`, `dplyr`, `data.table`, `readr`

---

## Example usage

```bash
Rscript scripts/Example-DESeq2.R
