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



# Merge overlapping loci only
merged <- merge_fuma_loci(fuma_list)

# Merge loci within 100kb
merged <- merge_fuma_loci(fuma_list, max_distance = 100000)

# View results
print_merged_loci(merged, n = 5)

