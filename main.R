setwd("/s/nath-lab/alsamman/____MyCodes____/lociR")
# Main script for lociR package
# Load libraries
source("libraries.R")

# Load all functions
source("functions/load_functions.R")
load_all_functions("functions")

# Note: the locus wizard exporter now uses helper templates in
# functions/locus_wizard_templates.R (loaded by load_all_functions()).

# # Main workflow
# fuma <- read_fuma("inputs/Hisp.txt")
# print(fuma)      # Quick overview
# summary(fuma)    # Detailed statistics
# fuma$data        # Access the data frame directly
# 
# # Single file (returns FUMA object)
# fuma <- read_fuma("inputs/Hisp.txt")

# Multiple files (returns named list of FUMA objects)
fuma_list <- read_fuma(c(
  "inputs/Asianq2_eur_hisp.txt",
  "inputs/Asianq2_eur.txt",
  "inputs/disc.txt",
  "inputs/Asianq2.txt",
  "inputs/EUR.txt",
  "inputs/Hisp.txt"
))


# Read FUMA files
# fuma_list <- read_fuma(c("inputs/Hisp.txt", "inputs/HispCHI.txt", "inputs/HispEUR.txt"))

# Merge loci
merged <- merge_fuma_loci(fuma_list)
merged
# Export to Excel with visualization
#export_merged_to_excel(merged, fuma_list, output_file = "merged_loci.xlsx")

# Merge overlapping loci only
#merged <- merge_fuma_loci(fuma_list)

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
aligned


# Multiple files
gwas_list <- read_gwas(
  c(
    "inputs/GWAS/Asianq2_eur_hisp.tsv",
    "inputs/GWAS/Asianq2_eur.tsv",
    "inputs/GWAS/asianq2.tsv",
    "inputs/GWAS/disc.tsv",
    "inputs/GWAS/EUR.tsv",
    "inputs/GWAS/hisp.tsv"
  ),
  name = c("Asianq2_eur_hisp", "Asianq2_eur", "asianq2","disc","EUR", "hisp"),
  pvalue_threshold = 5e-4,
  rdata_file = "outputs/gwas_cache.RData"
)


# Export aligned loci to Excel
#export_aligned_to_excel(aligned, "aligned_loci_all.xlsx", gwas_list = gwas_list)

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


# export_aligned_to_html_locus_wizard(
#   aligned_loci = aligned,
#   output_html = "outputs/aligned_locus_wizard.html",
#   flank_bp = 20000,
#   n_bins = 20,
#   gwas_list = gwas_list,
#   window_bp = 20000
# )

export_aligned_summary_excel(
  aligned_loci = aligned,
  output_file = "outputs/aligned_loci_summary_4.xlsx",
  gwas_list = gwas_list,
  annotation_file = "inputs/NCBI37.3.gene.loc"
)

export_merged_loci_gwas_table_excel(
  aligned_loci = aligned,
  gwas_list = gwas_list,
  output_file = "outputs/merged_loci_gwas_table.xlsx",
  annotation_file = "inputs/NCBI37.3.gene.loc"
)



source("libraries.R")

# Load all functions
source("functions/load_functions.R")
load_all_functions("functions")
# Example: one loci file + multiple GWAS files (read one-by-one with fread)
export_loci_gwas_target_followup_excel(
  loci_file = "inputs/myTargetLoci.txt",
  gwas_files = c(
    "inputs/GWAS/asianq2.tsv",
    "inputs/GWAS/disc.tsv",
    "inputs/GWAS/TWN.hg19.txt",
    "inputs/GWAS/hisp.tsv",
    "inputs/GWAS/EUR.tsv"
  ),
  dataset_names = c("asianq2", "disc", "TWN","Hisp","EUR"),
  target_datasets = c("asianq2", "disc", "TWN"),
  validation_datasets = c("TWN"),
  validation_p_threshold = 0.05,
  output_file = "outputs/loci_gwas_target_followup_Final.xlsx"
)

#Example: FUMA loci + one GWAS file -> lead SNP per locus (lowest p-value)
export_fuma_lead_snps_excel(
  fuma_file = "inputs/disc.txt",
  gwas_file = "raw_data/Chb8_standardized.tsv",
  output_file = "outputs/disc_lead_snps_from_gwas.xlsx"
)




source("libraries.R")

# Load all functions
source("functions/load_functions.R")
load_all_functions("functions")
# NewDisc + meta + individual CHB report (uses cached subsets when available)
export_newdisc_meta_individual_report_excel(
  loci_file = "inputs/newDisc.txt",
  meta_file = "inputs/GWAS/disc.tsv",
  individual_files = c(
    "inputs/GWAS/Chb/Chb1f.tsv", "inputs/GWAS/Chb/Chb2f.tsv",
    "inputs/GWAS/Chb/Chb3f.tsv", "inputs/GWAS/Chb/Chb4f.tsv",
    "inputs/GWAS/Chb/Chb1m.tsv", "inputs/GWAS/Chb/Chb2m.tsv",
    "inputs/GWAS/Chb/Chb3m.tsv", "inputs/GWAS/Chb/Chb4m.tsv"
  ),
  output_file = "outputs/newdisc_meta_individual_report.xlsx",
  subset_root_dir = "outputs/gwas_subsets_cache",
  p_threshold = 5e-5,
  annotation_file = "inputs/NCBI37.3.gene.loc",
  reported_file = "inputs/mehdi2024.txt",
  reported_max_distance_bp = 250000,
  individual_names = NULL,
  force_rebuild_subset = FALSE,
  summary_sheet = "Summary",
  report_sheet = "Locus_Report"
)



