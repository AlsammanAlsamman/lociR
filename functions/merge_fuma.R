#' Merge overlapping or nearby FUMA loci across multiple files
#'
#' @param fuma_list List of FUMA objects or single FUMA object
#' @param max_distance Integer. Maximum distance (bp) between loci to merge. Default: 0 (only overlapping)
#' @param report Logical. Print merge report. Default: TRUE
#' @return Data frame with merged loci and their source information
#' @export
merge_fuma_loci <- function(fuma_list, max_distance = 0, report = TRUE) {
  
  # Convert single FUMA to list
  if ("FUMA" %in% class(fuma_list)) {
    fuma_list <- list(fuma_list)
    names(fuma_list) <- fuma_list[[1]]$metadata$file_name
  }
  
  # Check if list is empty
  if (length(fuma_list) == 0) {
    stop("fuma_list is empty")
  }
  
  # Combine all loci from all FUMA files
  all_loci <- do.call(rbind, lapply(names(fuma_list), function(source_name) {
    fuma <- fuma_list[[source_name]]
    data <- fuma$data
    
    # Create data frame with necessary columns
    df <- data.frame(
      CHR = data$CHR,
      START = data$START,
      END = data$END,
      SOURCE = source_name,
      LOCUS_ID = if ("LOCUS" %in% colnames(data)) data$LOCUS else seq_len(nrow(data)),
      stringsAsFactors = FALSE
    )
    
    # Add original locus identifier
    df$ORIGINAL_LOCUS <- sprintf("%s(%s:%d:%d)", 
                                  df$LOCUS_ID, 
                                  df$SOURCE, 
                                  df$START, 
                                  df$END)
    
    return(df)
  }))
  
  if (nrow(all_loci) == 0) {
    stop("No loci found in FUMA files")
  }
  
  # Sort by chromosome and start position
  all_loci <- all_loci[order(all_loci$CHR, all_loci$START), ]
  
  message(sprintf("\nMerging %d loci from %d source(s) with max_distance = %d bp...\n", 
                  nrow(all_loci), length(fuma_list), max_distance))
  
  # Create GRanges object
  gr <- GenomicRanges::GRanges(
    seqnames = all_loci$CHR,
    ranges = IRanges::IRanges(
      start = all_loci$START,
      end = all_loci$END
    )
  )
  
  # Add metadata
  GenomicRanges::mcols(gr)$ORIGINAL_LOCUS <- all_loci$ORIGINAL_LOCUS
  GenomicRanges::mcols(gr)$SOURCE <- all_loci$SOURCE
  GenomicRanges::mcols(gr)$LOCUS_ID <- all_loci$LOCUS_ID
  
  # Reduce/merge overlapping and nearby regions
  if (max_distance > 0) {
    # Add distance to ranges before merging
    gr_reduced <- GenomicRanges::reduce(gr, 
                                        min.gapwidth = max_distance + 1,
                                        with.revmap = TRUE)
  } else {
    # Only merge overlapping
    gr_reduced <- GenomicRanges::reduce(gr, with.revmap = TRUE)
  }
  
  # Build merged loci data frame
  merged_loci <- data.frame(
    MERGED_LOCUS_ID = seq_along(gr_reduced),
    CHR = as.character(GenomicRanges::seqnames(gr_reduced)),
    START = GenomicRanges::start(gr_reduced),
    END = GenomicRanges::end(gr_reduced),
    WIDTH = GenomicRanges::width(gr_reduced),
    stringsAsFactors = FALSE
  )
  
  # Map back to original loci
  revmap <- GenomicRanges::mcols(gr_reduced)$revmap
  
  merged_loci$N_ORIGINAL_LOCI <- lengths(revmap)
  merged_loci$ORIGINAL_LOCI <- sapply(revmap, function(idx) {
    paste(all_loci$ORIGINAL_LOCUS[idx], collapse = ";")
  })
  
  merged_loci$SOURCES <- sapply(revmap, function(idx) {
    paste(unique(all_loci$SOURCE[idx]), collapse = ";")
  })
  
  # Generate report
  if (report) {
    cat("\n=== FUMA Loci Merge Report ===\n")
    cat(sprintf("Input loci: %d\n", nrow(all_loci)))
    cat(sprintf("Merged loci: %d\n", nrow(merged_loci)))
    cat(sprintf("Reduction: %d loci (%.1f%%)\n", 
                nrow(all_loci) - nrow(merged_loci),
                100 * (nrow(all_loci) - nrow(merged_loci)) / nrow(all_loci)))
    
    cat("\nLoci per chromosome:\n")
    print(table(merged_loci$CHR))
    
    cat("\nMerge statistics:\n")
    cat(sprintf("  Single locus (no merge): %d\n", sum(merged_loci$N_ORIGINAL_LOCI == 1)))
    cat(sprintf("  Merged from 2+ loci: %d\n", sum(merged_loci$N_ORIGINAL_LOCI > 1)))
    cat(sprintf("  Max loci merged: %d\n", max(merged_loci$N_ORIGINAL_LOCI)))
    
    cat("\nMerged locus size (bp):\n")
    print(summary(merged_loci$WIDTH))
    
    # Show examples of merged loci
    multi_merge <- merged_loci[merged_loci$N_ORIGINAL_LOCI > 1, ]
    if (nrow(multi_merge) > 0) {
      cat("\nExample merged loci (first 3):\n")
      for (i in seq_len(min(3, nrow(multi_merge)))) {
        cat(sprintf("\n  Merged Locus %d: chr%s:%d-%d (%d bp)\n", 
                    multi_merge$MERGED_LOCUS_ID[i],
                    multi_merge$CHR[i],
                    multi_merge$START[i],
                    multi_merge$END[i],
                    multi_merge$WIDTH[i]))
        cat(sprintf("    Combines %d loci from: %s\n", 
                    multi_merge$N_ORIGINAL_LOCI[i],
                    multi_merge$SOURCES[i]))
        
        # Show original loci (limit output)
        orig_loci <- strsplit(multi_merge$ORIGINAL_LOCI[i], ";")[[1]]
        if (length(orig_loci) > 5) {
          cat(sprintf("    Original loci: %s ... (and %d more)\n", 
                      paste(orig_loci[1:5], collapse = "; "),
                      length(orig_loci) - 5))
        } else {
          cat(sprintf("    Original loci: %s\n", paste(orig_loci, collapse = "; ")))
        }
      }
    }
    
    cat("\n==============================\n")
  }
  
  return(merged_loci)
}


#' Print method for merged loci
#' @export
print_merged_loci <- function(merged_loci, n = 10) {
  cat(sprintf("Merged FUMA Loci (%d total)\n", nrow(merged_loci)))
  cat("===============================\n\n")
  
  n_show <- min(n, nrow(merged_loci))
  
  for (i in seq_len(n_show)) {
    cat(sprintf("Locus %d: chr%s:%d-%d (%s bp)\n", 
                merged_loci$MERGED_LOCUS_ID[i],
                merged_loci$CHR[i],
                merged_loci$START[i],
                merged_loci$END[i],
                format(merged_loci$WIDTH[i], big.mark = ",")))
    
    cat(sprintf("  Sources: %s\n", merged_loci$SOURCES[i]))
    cat(sprintf("  Original loci (%d): %s\n\n", 
                merged_loci$N_ORIGINAL_LOCI[i],
                merged_loci$ORIGINAL_LOCI[i]))
  }
  
  if (nrow(merged_loci) > n) {
    cat(sprintf("... and %d more loci\n", nrow(merged_loci) - n))
  }
}
