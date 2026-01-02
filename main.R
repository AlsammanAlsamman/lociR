setwd("/s/nath-lab/alsamman/____MyCodes____/lociR")
# Main script for lociR package
# Load libraries
source("libraries.R")

# Load all functions
source("functions/load_functions.R")
load_all_functions("functions")

# Note: the locus wizard exporter now uses helper templates in
# functions/locus_wizard_templates.R (loaded by load_all_functions()).

# Main workflow
fuma <- read_fuma("inputs/Hisp.txt")
print(fuma)      # Quick overview
summary(fuma)    # Detailed statistics
fuma$data        # Access the data frame directly

# Single file (returns FUMA object)
fuma <- read_fuma("inputs/Hisp.txt")

# Multiple files (returns named list of FUMA objects)
fuma_list <- read_fuma(c(
  "inputs/CHI.txt",
  "inputs/EUR.txt",
  "inputs/Hisp.txt",
  "inputs/HispCHI.txt",
  "inputs/HispEUR.txt",
  "inputs/HispEURCHI.txt"
))


# Read FUMA files
# fuma_list <- read_fuma(c("inputs/Hisp.txt", "inputs/HispCHI.txt", "inputs/HispEUR.txt"))

# Merge loci
merged <- merge_fuma_loci(fuma_list)

# Export to Excel with visualization
#export_merged_to_excel(merged, fuma_list, output_file = "merged_loci.xlsx")

# Merge overlapping loci only
merged <- merge_fuma_loci(fuma_list)

# Merge loci within 100kb
#merged <- merge_fuma_loci(fuma_list, max_distance = 100000)

# View results
#print_merged_loci(merged, n = 5)

# Save all plots to folder (PDF)
#export_merged_plots(merged, output_folder = "plots/merged_loci")

# PNG format
#export_merged_plots(merged, output_folder = "plots/merged_loci", format = "png", dpi = 300)

# SVG format for high quality
#export_merged_plots(merged, output_folder = "plots/merged_loci", format = "svg")

# Custom size
#export_merged_plots(merged, output_folder = "plots", width = 12, height = 5)

# Read reported loci
mehdi <- read_reported_loci("inputs/mehdi2024.txt")
mehdi

# Align with reported loci
aligned <- align_with_reported(merged, mehdi)
summary(aligned)



# Multiple files
gwas_list <- read_gwas(
  c("inputs/GWAS/CHI.tsv", "inputs/GWAS/EUR.tsv", "inputs/GWAS/Hisp.tsv","inputs/GWAS/HispCHI.tsv","inputs/GWAS/HispEUR.tsv","inputs/GWAS/HispEURCHI.tsv"),
  name = c("CHI", "EUR", "Hisp","HispCHI","HispEUR","HispEURCHI"),
  pvalue_threshold = 5e-4,
  rdata_file = "outputs/gwas_cache.RData"
)


# Export aligned loci to Excel
#export_aligned_to_excel(aligned, "aligned_loci.xlsx", gwas_list = gwas_list)

# export_aligned_to_excel_roundtrip(aligned, "outputs/aligned_roundtrip.xlsx", gwas_list = gwas_list, flank_bp = 20000, n_bins = 20)
# 
# export_aligned_to_html_tsv(
#   aligned = aligned,
#   out_prefix = "outputs/aligned_review",
#   gwas_list = gwas_list,
#   flank_bp = 20000
# )
# 
# export_aligned_to_html_editable(
#   aligned_loci = aligned,
#   output_html = "outputs/aligned_editable.html",
#   flank_bp = 20000,
#   n_bins = 20
# )


export_aligned_to_html_locus_wizard(
  aligned_loci = aligned,
  output_html = "outputs/aligned_locus_wizard.html",
  flank_bp = 20000,
  n_bins = 20,
  gwas_list = gwas_list,
  window_bp = 20000
)
