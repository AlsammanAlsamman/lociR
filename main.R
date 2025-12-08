setwd("/s/nath-lab/alsamman/____MyCodes____/lociR")
# Main script for lociR package
# Load libraries
source("libraries.R")

# Load all functions
source("functions/load_functions.R")
load_all_functions()

# Main workflow
fuma <- read_fuma("inputs/Hisp.txt")
print(fuma)      # Quick overview
summary(fuma)    # Detailed statistics
fuma$data        # Access the data frame directly

# Single file (returns FUMA object)
fuma <- read_fuma("inputs/Hisp.txt")

# Multiple files (returns named list of FUMA objects)
fuma_list <- read_fuma(c(
  "inputs/Hisp.txt",
  "inputs/HispCHI.txt",
  "inputs/HispEUR.txt",
  "inputs/HispEURCHI.txt"
))


# Read FUMA files
fuma_list <- read_fuma(c("inputs/Hisp.txt", "inputs/HispCHI.txt", "inputs/HispEUR.txt"))

# Merge loci
merged <- merge_fuma_loci(fuma_list)

# Export to Excel with visualization
export_merged_to_excel(merged, fuma_list, output_file = "merged_loci.xlsx")

# Merge overlapping loci only
merged <- merge_fuma_loci(fuma_list)

# Merge loci within 100kb
merged <- merge_fuma_loci(fuma_list, max_distance = 100000)

# View results
print_merged_loci(merged, n = 5)

# Save all plots to folder (PDF)
export_merged_plots(merged, output_folder = "plots/merged_loci")

# PNG format
export_merged_plots(merged, output_folder = "plots/merged_loci", format = "png", dpi = 300)

# SVG format for high quality
export_merged_plots(merged, output_folder = "plots/merged_loci", format = "svg")

# Custom size
export_merged_plots(merged, output_folder = "plots", width = 12, height = 5)

# Read reported loci
mehdi <- read_reported_loci("inputs/mehdi2024.txt")
mehdi

# Refine with reported loci
refined <- refine_with_reported(merged, mehdi)
summary(refined)

# Export refined loci to Excel
export_refined_to_excel(refined, "refined_loci.xlsx")
