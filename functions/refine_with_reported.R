#' Refine Merged Loci Using Reported Loci
#'
#' Uses previously reported loci as guides to refine and split merged loci.
#' For each overlap, boundaries are expanded (min start, max end), and merged
#' loci are split based on reported loci positions.
#'
#' @param merged_loci data.frame with merged loci (output from merge_fuma_loci)
#' @param reported_loci data.frame or list of ReportedLoci objects
#' @param merge_overlapping_reported Logical, whether to merge overlapping reported 
#'        loci within the same merged region before splitting (default: TRUE)
#' @param keep_unmatched Logical, whether to keep merged loci with no overlapping 
#'        reported loci (default: TRUE)
#'
#' @return data.frame with refined loci including:
#'   - REFINED_LOCUS_ID: New locus identifier
#'   - CHR, START, END, WIDTH: Genomic coordinates
#'   - MERGED_LOCUS_ID: Original merged locus ID
#'   - REPORTED_LOCUS: Reported locus name(s) that guided the split
#'   - REPORTED_SOURCE: Source of the reported locus
#'   - N_ORIGINAL_LOCI: Number of original FUMA loci
#'   - ORIGINAL_LOCI: List of original locus IDs
#'   - SOURCES: Sources contributing to this region
#'   - REFINEMENT_TYPE: How this locus was created (split, expanded, unchanged)
#'
#' @export
refine_with_reported <- function(merged_loci,
                                  reported_loci,
                                  merge_overlapping_reported = TRUE,
                                  keep_unmatched = TRUE) {
  
  # Load required library
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges package is required. Install with: BiocManager::install('GenomicRanges')")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("IRanges package is required. Install with: BiocManager::install('IRanges')")
  }
  
  # Handle list of reported loci (combine them)
  if (is.list(reported_loci) && !is.data.frame(reported_loci)) {
    message("Combining ", length(reported_loci), " reported loci sources...")
    reported_loci <- do.call(rbind, reported_loci)
  }
  
  # Validate required columns
  required_merged <- c("CHR", "START", "END", "MERGED_LOCUS_ID")
  required_reported <- c("CHR", "START", "END", "LOCUS")
  
  if (!all(required_merged %in% colnames(merged_loci))) {
    stop("merged_loci missing required columns: ", 
         paste(setdiff(required_merged, colnames(merged_loci)), collapse = ", "))
  }
  
  if (!all(required_reported %in% colnames(reported_loci))) {
    stop("reported_loci missing required columns: ", 
         paste(setdiff(required_reported, colnames(reported_loci)), collapse = ", "))
  }
  
  message("Refining ", nrow(merged_loci), " merged loci using ", 
          nrow(reported_loci), " reported loci...")
  
  # Create GRanges objects
  merged_gr <- GenomicRanges::GRanges(
    seqnames = merged_loci$CHR,
    ranges = IRanges::IRanges(start = merged_loci$START, end = merged_loci$END)
  )
  
  reported_gr <- GenomicRanges::GRanges(
    seqnames = reported_loci$CHR,
    ranges = IRanges::IRanges(start = reported_loci$START, end = reported_loci$END)
  )
  
  # Find overlaps
  overlaps <- GenomicRanges::findOverlaps(merged_gr, reported_gr)
  
  message("Found ", length(unique(S4Vectors::queryHits(overlaps))), 
          " merged loci with overlapping reported loci")
  
  # Initialize results list
  refined_list <- list()
  refined_counter <- 1
  
  # Process each merged locus
  for (i in seq_len(nrow(merged_loci))) {
    
    # Find which reported loci overlap this merged locus
    reported_indices <- S4Vectors::subjectHits(overlaps)[S4Vectors::queryHits(overlaps) == i]
    
    if (length(reported_indices) == 0) {
      # No overlapping reported loci
      if (keep_unmatched) {
        refined_list[[refined_counter]] <- data.frame(
          REFINED_LOCUS_ID = paste0("REFINED_", refined_counter),
          CHR = merged_loci$CHR[i],
          START = merged_loci$START[i],
          END = merged_loci$END[i],
          WIDTH = merged_loci$END[i] - merged_loci$START[i] + 1,
          MERGED_LOCUS_ID = merged_loci$MERGED_LOCUS_ID[i],
          REPORTED_LOCUS = NA,
          REPORTED_SOURCE = NA,
          N_ORIGINAL_LOCI = merged_loci$N_ORIGINAL_LOCI[i],
          ORIGINAL_LOCI = merged_loci$ORIGINAL_LOCI[i],
          SOURCES = merged_loci$SOURCES[i],
          REFINEMENT_TYPE = "unchanged",
          stringsAsFactors = FALSE
        )
        refined_counter <- refined_counter + 1
      }
      next
    }
    
    # Get overlapping reported loci
    overlapping_reported <- reported_loci[reported_indices, ]
    
    # Merge overlapping reported loci if requested
    if (merge_overlapping_reported && nrow(overlapping_reported) > 1) {
      overlapping_reported <- merge_overlapping_reported_loci(overlapping_reported)
    }
    
    # Sort by start position
    overlapping_reported <- overlapping_reported[order(overlapping_reported$START), ]
    
    # Create refined loci based on reported guides
    for (j in seq_len(nrow(overlapping_reported))) {
      
      # Use reported locus boundaries, expanded only to include overlap with merged
      # Take the minimum start between reported and merged START
      refined_start <- min(overlapping_reported$START[j], merged_loci$START[i])
      # Take the maximum end between reported and merged END  
      refined_end <- max(overlapping_reported$END[j], merged_loci$END[i])
      
      # Constrain to the merged region boundaries (don't extend beyond merged)
      refined_start <- max(refined_start, merged_loci$START[i])
      refined_end <- min(refined_end, merged_loci$END[i])
      
      # Actually, use the reported boundaries as-is if they're within merged
      # Otherwise clip them to merged boundaries
      if (overlapping_reported$START[j] >= merged_loci$START[i] && 
          overlapping_reported$END[j] <= merged_loci$END[i]) {
        # Reported is fully within merged - use reported boundaries
        refined_start <- overlapping_reported$START[j]
        refined_end <- overlapping_reported$END[j]
        ref_type <- if (nrow(overlapping_reported) > 1) "split" else "aligned"
      } else {
        # Reported extends beyond merged - take the overlap region
        refined_start <- max(overlapping_reported$START[j], merged_loci$START[i])
        refined_end <- min(overlapping_reported$END[j], merged_loci$END[i])
        ref_type <- "expanded"
      }
      
      refined_list[[refined_counter]] <- data.frame(
        REFINED_LOCUS_ID = paste0("REFINED_", refined_counter),
        CHR = merged_loci$CHR[i],
        START = refined_start,
        END = refined_end,
        WIDTH = refined_end - refined_start + 1,
        MERGED_LOCUS_ID = merged_loci$MERGED_LOCUS_ID[i],
        REPORTED_LOCUS = overlapping_reported$LOCUS[j],
        REPORTED_SOURCE = if ("SOURCE_FILE" %in% colnames(overlapping_reported)) {
          overlapping_reported$SOURCE_FILE[j]
        } else {
          NA
        },
        N_ORIGINAL_LOCI = merged_loci$N_ORIGINAL_LOCI[i],
        ORIGINAL_LOCI = merged_loci$ORIGINAL_LOCI[i],
        SOURCES = merged_loci$SOURCES[i],
        REFINEMENT_TYPE = ref_type,
        stringsAsFactors = FALSE
      )
      
      refined_counter <- refined_counter + 1
    }
  }
  
  # Combine results
  refined_result <- do.call(rbind, refined_list)
  
  # Add class
  class(refined_result) <- c("RefinedLoci", "data.frame")
  
  # Print summary
  message("\n=== Refinement Summary ===")
  message("Total refined loci: ", nrow(refined_result))
  message("  - Unchanged: ", sum(refined_result$REFINEMENT_TYPE == "unchanged", na.rm = TRUE))
  message("  - Split: ", sum(refined_result$REFINEMENT_TYPE == "split", na.rm = TRUE))
  message("  - Expanded: ", sum(refined_result$REFINEMENT_TYPE == "expanded", na.rm = TRUE))
  message("  - Aligned: ", sum(refined_result$REFINEMENT_TYPE == "aligned", na.rm = TRUE))
  
  return(refined_result)
}


#' Merge Overlapping Reported Loci
#'
#' Internal function to merge overlapping reported loci within the same merged region
#'
#' @param reported_loci data.frame with reported loci (must be on same chromosome)
#'
#' @return data.frame with merged reported loci
#'
#' @keywords internal
merge_overlapping_reported_loci <- function(reported_loci) {
  
  # All should be on same chromosome
  if (length(unique(reported_loci$CHR)) > 1) {
    stop("All reported loci must be on the same chromosome for merging")
  }
  
  # Create GRanges
  reported_gr <- GenomicRanges::GRanges(
    seqnames = reported_loci$CHR,
    ranges = IRanges::IRanges(start = reported_loci$START, end = reported_loci$END)
  )
  
  # Reduce to merge overlapping
  reduced_gr <- GenomicRanges::reduce(reported_gr, with.revmap = TRUE)
  
  # Extract revmap
  revmap <- S4Vectors::mcols(reduced_gr)$revmap
  
  # Build merged data frame
  merged_reported <- data.frame(
    CHR = as.character(GenomicRanges::seqnames(reduced_gr)),
    START = GenomicRanges::start(reduced_gr),
    END = GenomicRanges::end(reduced_gr),
    stringsAsFactors = FALSE
  )
  
  # Combine locus names
  merged_reported$LOCUS <- sapply(revmap, function(indices) {
    paste(reported_loci$LOCUS[indices], collapse = ";")
  })
  
  # Combine source files if available
  if ("SOURCE_FILE" %in% colnames(reported_loci)) {
    merged_reported$SOURCE_FILE <- sapply(revmap, function(indices) {
      paste(unique(reported_loci$SOURCE_FILE[indices]), collapse = ";")
    })
  }
  
  return(merged_reported)
}


#' Print Method for RefinedLoci
#'
#' @param x RefinedLoci object
#' @param ... Additional arguments
#'
#' @export
print.RefinedLoci <- function(x, ...) {
  cat("Refined Loci Data\n")
  cat("=================\n")
  cat("Number of refined loci:", nrow(x), "\n")
  cat("Chromosomes:", length(unique(x$CHR)), "unique\n")
  cat("Refinement types:\n")
  cat("  - Unchanged:", sum(x$REFINEMENT_TYPE == "unchanged", na.rm = TRUE), "\n")
  cat("  - Split:", sum(x$REFINEMENT_TYPE == "split", na.rm = TRUE), "\n")
  cat("  - Expanded:", sum(x$REFINEMENT_TYPE == "expanded", na.rm = TRUE), "\n")
  cat("  - Aligned:", sum(x$REFINEMENT_TYPE == "aligned", na.rm = TRUE), "\n")
  
  cat("\nColumns:", paste(colnames(x), collapse = ", "), "\n")
  cat("\nFirst few loci:\n")
  print.data.frame(head(x, 10))
  
  invisible(x)
}


#' Summary Method for RefinedLoci
#'
#' @param object RefinedLoci object
#' @param ... Additional arguments
#'
#' @export
summary.RefinedLoci <- function(object, ...) {
  cat("Refined Loci Summary\n")
  cat("====================\n")
  cat("Total refined loci:", nrow(object), "\n")
  cat("Chromosomes:", paste(sort(unique(object$CHR)), collapse = ", "), "\n")
  
  cat("\nRefinement breakdown:\n")
  ref_table <- table(object$REFINEMENT_TYPE)
  for (ref_type in names(ref_table)) {
    cat("  ", ref_type, ": ", ref_table[ref_type], " (", 
        round(100 * ref_table[ref_type] / nrow(object), 1), "%)\n", sep = "")
  }
  
  cat("\nLocus width statistics (bp):\n")
  cat("  Min:", min(object$WIDTH, na.rm = TRUE), "\n")
  cat("  Mean:", round(mean(object$WIDTH, na.rm = TRUE), 2), "\n")
  cat("  Median:", median(object$WIDTH, na.rm = TRUE), "\n")
  cat("  Max:", max(object$WIDTH, na.rm = TRUE), "\n")
  
  # Count matched vs unmatched
  n_matched <- sum(!is.na(object$REPORTED_LOCUS))
  cat("\nMatched with reported loci:", n_matched, "(", 
      round(100 * n_matched / nrow(object), 1), "%)\n")
  
  invisible(object)
}
