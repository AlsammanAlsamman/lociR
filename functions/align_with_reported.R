#' Align Merged Loci Using Reported Loci
#'
#' Aligns merged loci with previously reported loci without altering locus boundaries
#' through splitting or expansion. For each merged locus the function records the
#' reported loci that overlap (or fall within a specified proximity) and returns a
#' single aligned summary row.
#'
#' @param merged_loci data.frame with merged loci (output from merge_fuma_loci).
#'        Must contain columns `CHR`, `START`, `END`, and `MERGED_LOCUS_ID`.
#' @param reported_loci data.frame or list of ReportedLoci objects. Must contain
#'        columns `CHR`, `START`, `END`, and `LOCUS`.
#' @param max_distance integer. Maximum distance (bp) to consider a reported locus
#'        as aligned with a merged locus. Default: 0 (direct overlaps only).
#' @param keep_unmatched logical, whether to keep merged loci that have no aligned
#'        reported loci. Default: TRUE.
#'
#' @return data.frame with class `AlignedLoci`. Each row corresponds to a merged
#'         locus and includes alignment metadata:
#'   - ALIGNED_LOCUS_ID: Row identifier
#'   - CHR: Chromosome
#'   - MERGED_START / MERGED_END / MERGED_WIDTH: Original merged boundaries
#'   - ALIGNED_START / ALIGNED_END / ALIGNED_WIDTH: Overlap range if available
#'   - REPORTED_START / REPORTED_END / REPORTED_LOCI / REPORTED_SOURCE: Reported summary
#'   - ALIGNMENT_STATUS: `aligned`, `proximal`, or `unmatched`
#'   - N_MATCHED_REPORTED: Number of reported loci associated with the merged locus
#'   - DISTANCE_TO_REPORTED: Minimum distance to a reported locus when proximal
#'   - MERGED_LOCUS_ID, N_ORIGINAL_LOCI, ORIGINAL_LOCI, SOURCES: Carried metadata
#'
#' @export
align_with_reported <- function(merged_loci,
                                reported_loci,
                                max_distance = 0,
                                keep_unmatched = TRUE) {
  # Dependency checks (required for interval arithmetic)
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges package is required. Install with: BiocManager::install('GenomicRanges')")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("IRanges package is required. Install with: BiocManager::install('IRanges')")
  }

  # Helper to safely extract column values
  get_col_value <- function(df, row_idx, col_name, default = NA) {
    if (col_name %in% colnames(df)) {
      return(df[[col_name]][row_idx])
    }
    default
  }

  # Handle list input of reported loci (combine)
  if (is.list(reported_loci) && !is.data.frame(reported_loci)) {
    message("Combining ", length(reported_loci), " reported loci sources ...")
    reported_loci <- do.call(rbind, lapply(reported_loci, function(x) {
      if (!is.data.frame(x)) {
        stop("All elements in reported_loci list must be data frames")
      }
      x
    }))
  }

  # Validate inputs
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

  merged_loci$CHR <- as.character(merged_loci$CHR)
  reported_loci$CHR <- as.character(reported_loci$CHR)

  message("Aligning ", nrow(merged_loci), " merged loci with ",
          nrow(reported_loci), " reported loci ...")
  message("Maximum distance threshold: ", format(max_distance, big.mark = ","), " bp")

  # Construct GRanges
  merged_gr <- GenomicRanges::GRanges(
    seqnames = merged_loci$CHR,
    ranges = IRanges::IRanges(start = merged_loci$START, end = merged_loci$END)
  )
  reported_gr <- GenomicRanges::GRanges(
    seqnames = reported_loci$CHR,
    ranges = IRanges::IRanges(start = reported_loci$START, end = reported_loci$END)
  )

  overlaps <- GenomicRanges::findOverlaps(merged_gr, reported_gr, maxgap = max(0, max_distance))

  message("Found ", length(unique(S4Vectors::queryHits(overlaps))),
          " merged loci with aligned or proximal reported loci")

  aligned_rows <- vector("list", length = nrow(merged_loci))
  aligned_counter <- 1L

  for (i in seq_len(nrow(merged_loci))) {
    merged_start <- merged_loci$START[i]
    merged_end <- merged_loci$END[i]

    subject_idx <- S4Vectors::subjectHits(overlaps)[S4Vectors::queryHits(overlaps) == i]

    if (length(subject_idx) == 0) {
      if (!keep_unmatched) {
        next
      }

      aligned_rows[[aligned_counter]] <- data.frame(
        ALIGNED_LOCUS_ID = paste0("ALIGNED_", aligned_counter),
        CHR = merged_loci$CHR[i],
        MERGED_START = merged_start,
        MERGED_END = merged_end,
        MERGED_WIDTH = merged_end - merged_start + 1,
        ALIGNED_START = NA_integer_,
        ALIGNED_END = NA_integer_,
        ALIGNED_WIDTH = NA_integer_,
        REPORTED_START = NA_integer_,
        REPORTED_END = NA_integer_,
        REPORTED_LOCI = NA_character_,
        REPORTED_SOURCE = NA_character_,
        ALIGNMENT_STATUS = "unmatched",
        N_MATCHED_REPORTED = 0L,
        DISTANCE_TO_REPORTED = NA_integer_,
        MERGED_LOCUS_ID = merged_loci$MERGED_LOCUS_ID[i],
        N_ORIGINAL_LOCI = get_col_value(merged_loci, i, "N_ORIGINAL_LOCI", NA_integer_),
        ORIGINAL_LOCI = get_col_value(merged_loci, i, "ORIGINAL_LOCI", NA_character_),
        SOURCES = get_col_value(merged_loci, i, "SOURCES", NA_character_),
        stringsAsFactors = FALSE
      )
      aligned_counter <- aligned_counter + 1L
      next
    }

    reported_subset <- reported_loci[subject_idx, , drop = FALSE]

    reported_start <- min(reported_subset$START, na.rm = TRUE)
    reported_end <- max(reported_subset$END, na.rm = TRUE)

    overlap_start <- max(merged_start, reported_start)
    overlap_end <- min(merged_end, reported_end)
    overlap_width <- overlap_end - overlap_start + 1

    if (overlap_width > 0) {
      alignment_status <- "aligned"
      aligned_start <- overlap_start
      aligned_end <- overlap_end
      aligned_width <- overlap_width
      distance_to_reported <- 0L
    } else {
      alignment_status <- "proximal"
      aligned_start <- NA_integer_
      aligned_end <- NA_integer_
      aligned_width <- NA_integer_

      # Compute minimum gap between merged locus and reported loci
      gaps_left <- merged_start - reported_subset$END
      gaps_left <- gaps_left[gaps_left > 0]
      gaps_right <- reported_subset$START - merged_end
      gaps_right <- gaps_right[gaps_right > 0]
      candidate_gaps <- c(gaps_left, gaps_right)
      distance_to_reported <- if (length(candidate_gaps) == 0) 0L else min(candidate_gaps)
    }

    reported_names <- paste(unique(reported_subset$LOCUS), collapse = ";")
    reported_sources <- if ("SOURCE_FILE" %in% colnames(reported_subset)) {
      paste(unique(reported_subset$SOURCE_FILE), collapse = ";")
    } else {
      NA_character_
    }

    aligned_rows[[aligned_counter]] <- data.frame(
      ALIGNED_LOCUS_ID = paste0("ALIGNED_", aligned_counter),
      CHR = merged_loci$CHR[i],
      MERGED_START = merged_start,
      MERGED_END = merged_end,
      MERGED_WIDTH = merged_end - merged_start + 1,
      ALIGNED_START = aligned_start,
      ALIGNED_END = aligned_end,
      ALIGNED_WIDTH = aligned_width,
      REPORTED_START = reported_start,
      REPORTED_END = reported_end,
      REPORTED_LOCI = reported_names,
      REPORTED_SOURCE = reported_sources,
      ALIGNMENT_STATUS = alignment_status,
      N_MATCHED_REPORTED = nrow(reported_subset),
      DISTANCE_TO_REPORTED = distance_to_reported,
      MERGED_LOCUS_ID = merged_loci$MERGED_LOCUS_ID[i],
      N_ORIGINAL_LOCI = get_col_value(merged_loci, i, "N_ORIGINAL_LOCI", NA_integer_),
      ORIGINAL_LOCI = get_col_value(merged_loci, i, "ORIGINAL_LOCI", NA_character_),
      SOURCES = get_col_value(merged_loci, i, "SOURCES", NA_character_),
      stringsAsFactors = FALSE
    )

    aligned_counter <- aligned_counter + 1L
  }

  aligned_rows <- aligned_rows[seq_len(aligned_counter - 1L)]
  if (length(aligned_rows) == 0) {
    message("No aligned loci returned (keep_unmatched = FALSE and no matches found)")
    return(data.frame())
  }

  aligned_result <- do.call(rbind, aligned_rows)
  class(aligned_result) <- c("AlignedLoci", "data.frame")

  message("\n=== Alignment Summary ===")
  message("Total aligned loci: ", nrow(aligned_result))
  status_counts <- table(aligned_result$ALIGNMENT_STATUS, useNA = "ifany")
  for (status in names(status_counts)) {
    message("  - ", status, ": ", status_counts[[status]])
  }

  aligned_result
}

#' Print method for AlignedLoci
#'
#' @param x AlignedLoci object
#' @param ... Additional arguments
#'
#' @export
print.AlignedLoci <- function(x, ...) {
  cat("Aligned Loci Data\n")
  cat("=================\n")
  cat("Number of loci:", nrow(x), "\n")
  cat("Chromosomes:", length(unique(x$CHR)), "unique\n")
  cat("Alignment statuses:\n")
  status_counts <- table(x$ALIGNMENT_STATUS, useNA = "ifany")
  for (status in names(status_counts)) {
    cat("  - ", status, ": ", status_counts[[status]], "\n", sep = "")
  }

  cat("\nColumns:", paste(colnames(x), collapse = ", "), "\n")
  cat("\nFirst few loci:\n")
  print.data.frame(head(x, 10))

  invisible(x)
}

#' Summary method for AlignedLoci
#'
#' @param object AlignedLoci object
#' @param ... Additional arguments
#'
#' @export
summary.AlignedLoci <- function(object, ...) {
  cat("Aligned Loci Summary\n")
  cat("=====================\n")
  cat("Total loci:", nrow(object), "\n")
  cat("Chromosomes:", paste(sort(unique(object$CHR)), collapse = ", "), "\n")

  status_counts <- table(object$ALIGNMENT_STATUS, useNA = "ifany")
  cat("\nAlignment breakdown:\n")
  for (status in names(status_counts)) {
    cat("  ", status, ": ", status_counts[[status]], " (",
        round(100 * status_counts[[status]] / nrow(object), 1), "%)\n", sep = "")
  }

  aligned_rows <- object$ALIGNMENT_STATUS == "aligned" & !is.na(object$ALIGNED_WIDTH)
  if (any(aligned_rows)) {
    cat("\nAligned overlap width (bp):\n")
    cat("  Min:", min(object$ALIGNED_WIDTH[aligned_rows], na.rm = TRUE), "\n")
    cat("  Mean:", round(mean(object$ALIGNED_WIDTH[aligned_rows], na.rm = TRUE), 2), "\n")
    cat("  Median:", median(object$ALIGNED_WIDTH[aligned_rows], na.rm = TRUE), "\n")
    cat("  Max:", max(object$ALIGNED_WIDTH[aligned_rows], na.rm = TRUE), "\n")
  }

  proximal_rows <- object$ALIGNMENT_STATUS == "proximal"
  if (any(proximal_rows)) {
    cat("\nProximal distance to reported loci (bp):\n")
    cat("  Min:", min(object$DISTANCE_TO_REPORTED[proximal_rows], na.rm = TRUE), "\n")
    cat("  Mean:", round(mean(object$DISTANCE_TO_REPORTED[proximal_rows], na.rm = TRUE), 2), "\n")
    cat("  Median:", median(object$DISTANCE_TO_REPORTED[proximal_rows], na.rm = TRUE), "\n")
    cat("  Max:", max(object$DISTANCE_TO_REPORTED[proximal_rows], na.rm = TRUE), "\n")
  }

  invisible(object)
}
