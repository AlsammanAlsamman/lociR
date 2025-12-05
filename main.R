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

