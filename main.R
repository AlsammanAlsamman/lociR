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

# Access individual files
fuma_list$Hisp.txt
fuma_list$HispCHI.txt

# Or use list.files to get all FUMA files
all_fuma <- read_fuma(list.files("inputs", pattern = "^Hisp.*\\.txt$", full.names = TRUE))

