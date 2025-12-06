# Required libraries for lociR

# Genomic ranges manipulation
library(GenomicRanges)
library(IRanges)

# Visualization
library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)

# Excel export: install if missing and load
if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    install.packages("openxlsx2")
}
library(openxlsx2)
