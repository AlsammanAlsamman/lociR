# lociR

A small R workflow for working with GWAS genomic loci across multiple sources (FUMA outputs, reported loci tables, and GWAS summary statistics), including:

- Reading multiple FUMA loci files (per ancestry/cohort)
- Merging overlapping/nearby loci across sources
- Aligning merged loci to a set of previously reported loci
- Exporting a visual Excel workbook for review and *round-trip* manual refinement

This repo is organized as a set of functions in `functions/` plus example driver scripts.

---

## Quick start

1) Open R in the project folder.

2) Install/load required packages:

- See `libraries.R` (auto-installs `openxlsx2` if missing).

3) Run the example workflow:

- See `main.R` for an end-to-end example.

> Note (Windows): `main.R` currently uses `setwd("/s/.../")` (a Linux path). Update that line to your local path or remove it and run from the project root.

---

## Key scripts

- `main.R`
  - Example end-to-end workflow: read FUMA → merge loci → align to reported loci → read GWAS (optional) → export Excel.

- `calculate_loci.R`
  - A separate loci-identification workflow that uses a PLINK reference panel to produce `GenomicRiskLoci.txt` (see `prompt/whatThisisDoing.md` for a detailed description).

---

## Functions (high level)

All functions live in `functions/` and are sourced by `functions/load_functions.R`.

Commonly used:

- `read_fuma()`
  - Reads one or more FUMA loci files from `inputs/`.

- `merge_fuma_loci()`
  - Merges loci across multiple FUMA sources.

- `read_reported_loci()`
  - Reads a table of previously reported loci (literature/curation).

- `align_with_reported()`
  - Aligns merged loci to reported loci.

- `read_gwas()`
  - Reads GWAS summary statistic file(s), standardizes columns, and filters by p-value threshold.
  - Supports an on-disk cache via `rdata_file` so repeated runs can skip re-reading and re-filtering GWAS TSV files.

- `export_aligned_to_excel_roundtrip()`
  - Writes a workbook designed for a human-in-the-loop refinement workflow:
    - `Aligned_Loci`: main sheet for visual review + manual edits
    - `Legend` / `Sources`: helpful reference
    - `Meta`, `Parameters`, `Metadata`: machine-facing sheets used to safely round-trip edits back to R

---

## Round-trip Excel workflow (stars editor)

The intended workflow is documented in:

- `prompt/README_refined_loci_workflow.md`

Conceptually:

1) Export an Excel workbook.
2) In Excel, place `*` in the `tuning` row bins to indicate refined interval(s).
3) Provide names in the `naming` row (semicolon-separated for multiple intervals).
4) Import the edited workbook back into R (importer is planned/iterating).

---

## Inputs / Outputs

- `inputs/`
  - FUMA loci files (e.g., `Hisp.txt`) and reported loci tables (e.g., `mehdi2024.txt`).

- `inputs/GWAS/`
  - GWAS summary statistics tables (TSV).

- `outputs/`
  - Generated results (e.g., Excel exports, locus tables).

- `plots/`
  - Plot outputs (if using plotting helpers).

---

## Example usage (outline)

```r
source("libraries.R")
source("functions/load_functions.R")
load_all_functions()

# Read FUMA loci
fuma_list <- read_fuma(c(
  "inputs/CHI.txt",
  "inputs/EUR.txt",
  "inputs/Hisp.txt"
))

# Merge and align
merged <- merge_fuma_loci(fuma_list)
reported <- read_reported_loci("inputs/mehdi2024.txt")
aligned <- align_with_reported(merged, reported)

# Optional GWAS (cached)
gwas_list <- read_gwas(
  c("inputs/GWAS/CHI.tsv", "inputs/GWAS/EUR.tsv", "inputs/GWAS/Hisp.tsv"),
  name = c("CHI", "EUR", "Hisp"),
  pvalue_threshold = 5e-4,
  rdata_file = "outputs/gwas_cache.RData"
)

# Export workbook for manual refinement
export_aligned_to_excel_roundtrip(
  aligned,
  output_file = "outputs/aligned_roundtrip.xlsx",
  gwas_list = gwas_list,
  flank_bp = 20000,
  n_bins = 20
)
```

---

## Notes

- Excel export runtime is usually dominated by per-cell styling/writes (`openxlsx2`) and per-row GWAS bin aggregation when `gwas_list` is provided.
- This repository is currently script-first (not packaged as an installed R package).
