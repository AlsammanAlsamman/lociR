#!/usr/bin/env bash
# Extract a SNP subset from the AMR 1000G reference panel using PLINK.
#
# Usage:
#   ./extract_amr_subset.sh <gwas_tsv> <plink_bfile_prefix> <output_prefix>
#
# Example:
#   ./extract_amr_subset.sh Hisp_filtered_p0.1.tsv \
#       /s/nath-lab/alsamman/__Data__/1000GenomeData_phase3_2013_PEL/plinkformat/AMR \
#       amr_subsets/Hisp_AMR_subset

set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: $0 <gwas_tsv> <plink_bfile_prefix> <output_prefix>" >&2
    exit 1
fi

GWAS_TSV=$1
PLINK_PREFIX=$2
OUT_PREFIX=$3

if [[ ! -f "$GWAS_TSV" ]]; then
    echo "Error: GWAS file '$GWAS_TSV' not found." >&2
    exit 1
fi

for ext in bed bim fam; do
    if [[ ! -f "${PLINK_PREFIX}.${ext}" ]]; then
        echo "Error: PLINK reference '${PLINK_PREFIX}.${ext}' not found." >&2
        exit 1
    fi
done

if ! command -v plink >/dev/null 2>&1; then
    echo "Error: plink not found in PATH." >&2
    exit 1
fi

mkdir -p "$(dirname "$OUT_PREFIX")"

SNP_LIST=$(mktemp)
trap 'rm -f "$SNP_LIST"' EXIT

awk -F '\t' 'NR > 1 && $3 != "" {print $3}' "$GWAS_TSV" | sort -u > "$SNP_LIST"

if [[ ! -s "$SNP_LIST" ]]; then
    echo "Error: no SNP IDs extracted from '$GWAS_TSV'." >&2
    exit 1
fi

echo "Extracting $(wc -l < "$SNP_LIST") SNPs from $PLINK_PREFIX ..."

plink \
    --bfile "$PLINK_PREFIX" \
    --extract "$SNP_LIST" \
    --make-bed \
    --out "$OUT_PREFIX" \
    --allow-extra-chr

echo "Subset written to ${OUT_PREFIX}.bed/.bim/.fam"
