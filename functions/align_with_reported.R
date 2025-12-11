#' Align Merged Loci Using Reported Loci
#'
#' Aligns merged loci with previously reported loci without altering boundaries.
#' For each merged locus the function records overlapping (or nearby) reported
#' loci and preserves the original merged interval. No splitting or extension of
#' merged loci is performed.
#'
#' @param merged_loci data.frame produced by merge_fuma_loci(). Must contain
#'        columns `CHR`, `START`, `END`, and `MERGED_LOCUS_ID`.
#' @param reported_loci data.frame or list of ReportedLoci objects. Must contain
#'        columns `CHR`, `START`, `END`, and `LOCUS`.
#' @param max_distance integer. Maximum distance (bp) to consider a reported locus
#'        as aligned with a merged locus. Default: 0 (direct overlaps only).
#' @param keep_unmatched logical. Keep merged loci with no nearby reported loci.
#'        Default: TRUE.
#'
#' @return data.frame with class `AlignedLoci`. Each row corresponds to a merged
#'         locus paired with either an aligned/proximal reported locus or an
#'         unmatched record when no reported loci were found. Columns include:
#'   - ALIGNED_LOCUS_ID: Row identifier
#'   - ROW_TYPE: `aligned`, `proximal`, or `unmatched`
#'   - CHR, START, END, WIDTH: Reported locus coordinates (or NA when unmatched)
#'   - MERGED_LOCUS_ID, MERGED_START, MERGED_END, MERGED_WIDTH: Original merge
#'   - REPORTED_LOCUS, REPORTED_SOURCE: Reported locus metadata
#'   - ALIGNMENT_STATUS: Same as ROW_TYPE for convenience
#'   - DISTANCE_TO_MERGED: Gap distance when proximal (bp)
#'   - N_ORIGINAL_LOCI, ORIGINAL_LOCI, SOURCES: Metadata inherited from merge
#'
#' @export
align_with_reported <- function(merged_loci,
                                reported_loci,
                                max_distance = 0,
                                keep_unmatched = TRUE) {
  if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
    stop("GenomicRanges package is required. Install with: BiocManager::install('GenomicRanges')")
  }
  if (!requireNamespace("IRanges", quietly = TRUE)) {
    stop("IRanges package is required. Install with: BiocManager::install('IRanges')")
  }

  # Helper to safely access optional columns
  fetch_col <- function(df, idx, col_name, default = NA) {
    if (col_name %in% colnames(df)) {
      return(df[[col_name]][idx])
    }
    default
  }

  # Combine list input of reported loci
  if (is.list(reported_loci) && !is.data.frame(reported_loci)) {
    message("Combining ", length(reported_loci), " reported loci sources...")
    reported_loci <- do.call(rbind, lapply(reported_loci, function(x) {
      if (!is.data.frame(x)) {
        stop("All elements in reported_loci list must be data frames")
      }
      x
    }))
  }

  required_merged <- c("CHR", "START", "END", "MERGED_LOCUS_ID")
  required_reported <- c("CHR", "START", "END", "LOCUS")

  if (!all(required_merged %in% colnames(merged_loci))) {
    stop("merged_loci missing required columns: ",
         paste(setdiff(required_merged, colnames(merged_loci)), collapse = ", "))
  }
  if (nrow(reported_loci) > 0 && !all(required_reported %in% colnames(reported_loci))) {
    stop("reported_loci missing required columns: ",
         paste(setdiff(required_reported, colnames(reported_loci)), collapse = ", "))
  }

  merged_loci <- as.data.frame(merged_loci, stringsAsFactors = FALSE)
  reported_loci <- as.data.frame(reported_loci, stringsAsFactors = FALSE)

  merged_loci$CHR <- as.character(merged_loci$CHR)
  merged_loci$START <- as.numeric(merged_loci$START)
  merged_loci$END <- as.numeric(merged_loci$END)
  if (!"WIDTH" %in% colnames(merged_loci)) {
    merged_loci$WIDTH <- merged_loci$END - merged_loci$START + 1
  }

  if (nrow(reported_loci) > 0) {
    reported_loci$CHR <- as.character(reported_loci$CHR)
    reported_loci$START <- as.numeric(reported_loci$START)
    reported_loci$END <- as.numeric(reported_loci$END)
  }

  message("Aligning ", nrow(merged_loci), " merged loci with ",
          nrow(reported_loci), " reported loci ...")
  message("Maximum distance threshold: ", format(max_distance, big.mark = ","), " bp")

  merged_gr <- GenomicRanges::GRanges(
    seqnames = merged_loci$CHR,
    ranges = IRanges::IRanges(start = merged_loci$START, end = merged_loci$END)
  )

  if (nrow(reported_loci) > 0) {
    reported_gr <- GenomicRanges::GRanges(
      seqnames = reported_loci$CHR,
      ranges = IRanges::IRanges(start = reported_loci$START, end = reported_loci$END)
    )
    overlaps <- GenomicRanges::findOverlaps(merged_gr,
                                            reported_gr,
                                            maxgap = max(0, max_distance))
    query_hits <- S4Vectors::queryHits(overlaps)
    subject_hits <- S4Vectors::subjectHits(overlaps)
  } else {
    query_hits <- integer(0)
    subject_hits <- integer(0)
  }

  aligned_rows <- list()
  row_counter <- 1L
  matched_flags <- rep(FALSE, nrow(merged_loci))

  for (i in seq_len(nrow(merged_loci))) {
    merged_row <- merged_loci[i, ]
    merged_start <- merged_row$START
    merged_end <- merged_row$END

    if (length(query_hits) > 0) {
      reported_idx <- subject_hits[query_hits == i]
    } else {
      reported_idx <- integer(0)
    }

    if (length(reported_idx) == 0) {
      if (keep_unmatched) {
        aligned_rows[[row_counter]] <- data.frame(
          ALIGNED_LOCUS_ID = paste0("ALIGNED_", row_counter),
          ROW_TYPE = "unmatched",
          CHR = merged_row$CHR,
          START = NA_real_,
          END = NA_real_,
          WIDTH = NA_real_,
          MERGED_LOCUS_ID = merged_row$MERGED_LOCUS_ID,
          MERGED_START = merged_start,
          MERGED_END = merged_end,
          MERGED_WIDTH = merged_row$WIDTH,
          REPORTED_LOCUS = NA_character_,
          REPORTED_SOURCE = NA_character_,
          ALIGNMENT_STATUS = "unmatched",
          DISTANCE_TO_MERGED = NA_real_,
          N_ORIGINAL_LOCI = fetch_col(merged_loci, i, "N_ORIGINAL_LOCI", NA_real_),
          ORIGINAL_LOCI = fetch_col(merged_loci, i, "ORIGINAL_LOCI", NA_character_),
          SOURCES = fetch_col(merged_loci, i, "SOURCES", NA_character_),
          stringsAsFactors = FALSE
        )
        row_counter <- row_counter + 1L
      }
      next
    }

    matched_flags[i] <- TRUE

    for (idx in reported_idx) {
      rep_row <- reported_loci[idx, , drop = FALSE]
      rep_start <- rep_row$START
      rep_end <- rep_row$END
      rep_width <- if (is.finite(rep_start) && is.finite(rep_end)) rep_end - rep_start + 1 else NA_real_

      overlap_start <- max(merged_start, rep_start, na.rm = TRUE)
      overlap_end <- min(merged_end, rep_end, na.rm = TRUE)
      overlap_width <- overlap_end - overlap_start + 1

      has_overlap <- is.finite(overlap_width) && overlap_width > 0
      row_type <- if (has_overlap) "aligned" else "proximal"

      if (has_overlap) {
        distance_to_merged <- 0
      } else {
        if (!is.finite(rep_start) || !is.finite(rep_end)) {
          distance_to_merged <- NA_real_
        } else if (rep_end < merged_start) {
          distance_to_merged <- merged_start - rep_end
        } else if (rep_start > merged_end) {
          distance_to_merged <- rep_start - merged_end
        } else {
          distance_to_merged <- 0
        }
      }

      reported_source <- if ("SOURCE_FILE" %in% colnames(rep_row)) {
        rep_row$SOURCE_FILE
      } else if ("SOURCE" %in% colnames(rep_row)) {
        rep_row$SOURCE
      } else {
        NA_character_
      }

      aligned_rows[[row_counter]] <- data.frame(
        ALIGNED_LOCUS_ID = paste0("ALIGNED_", row_counter),
        ROW_TYPE = row_type,
        CHR = merged_row$CHR,
        START = rep_start,
        END = rep_end,
        WIDTH = rep_width,
        MERGED_LOCUS_ID = merged_row$MERGED_LOCUS_ID,
        MERGED_START = merged_start,
        MERGED_END = merged_end,
        MERGED_WIDTH = merged_row$WIDTH,
        REPORTED_LOCUS = rep_row$LOCUS,
        REPORTED_SOURCE = reported_source,
        ALIGNMENT_STATUS = row_type,
        DISTANCE_TO_MERGED = distance_to_merged,
        N_ORIGINAL_LOCI = fetch_col(merged_loci, i, "N_ORIGINAL_LOCI", NA_real_),
        ORIGINAL_LOCI = fetch_col(merged_loci, i, "ORIGINAL_LOCI", NA_character_),
        SOURCES = fetch_col(merged_loci, i, "SOURCES", NA_character_),
        stringsAsFactors = FALSE
      )
      row_counter <- row_counter + 1L
    }
  }

  if (length(aligned_rows) == 0) {
    message("No aligned loci returned (keep_unmatched = FALSE and no reported loci matched)")
    return(data.frame())
  }

  aligned_result <- do.call(rbind, aligned_rows)
  aligned_result$ROW_TYPE <- as.character(aligned_result$ROW_TYPE)
  aligned_result$ALIGNMENT_STATUS <- as.character(aligned_result$ALIGNMENT_STATUS)
  aligned_result$MERGED_LOCUS_ID <- as.character(aligned_result$MERGED_LOCUS_ID)

  class(aligned_result) <- c("AlignedLoci", "data.frame")

  message("\n=== Alignment Summary ===")
  message("Merged loci considered: ", nrow(merged_loci))
  message("Merged loci with reported matches: ", sum(matched_flags))
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
  cat("Rows:", nrow(x), "\n")
  cat("Unique merged loci:", length(unique(x$MERGED_LOCUS_ID)), "\n")
  cat("Chromosomes:", length(unique(x$CHR)), "\n")

  status_counts <- table(x$ALIGNMENT_STATUS, useNA = "ifany")
  cat("Alignment statuses:\n")
  for (status in names(status_counts)) {
    cat("  - ", status, ": ", status_counts[[status]], "\n", sep = "")
  }

  cat("\nColumns:", paste(colnames(x), collapse = ", "), "\n")
  cat("\nFirst few rows:\n")
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
  cat("Rows:", nrow(object), "\n")
  cat("Unique merged loci:", length(unique(object$MERGED_LOCUS_ID)), "\n")
  cat("Chromosomes:", paste(sort(unique(object$CHR)), collapse = ", "), "\n")

  status_counts <- table(object$ALIGNMENT_STATUS, useNA = "ifany")
  cat("\nAlignment breakdown:\n")
  for (status in names(status_counts)) {
    cat("  ", status, ": ", status_counts[[status]], " (",
        round(100 * status_counts[[status]] / nrow(object), 1), "%)\n", sep = "")
  }

  aligned_rows <- object$ALIGNMENT_STATUS == "aligned" & is.finite(object$WIDTH)
  if (any(aligned_rows)) {
    cat("\nAligned reported widths (bp):\n")
    cat("  Min:", min(object$WIDTH[aligned_rows], na.rm = TRUE), "\n")
    cat("  Mean:", round(mean(object$WIDTH[aligned_rows], na.rm = TRUE), 2), "\n")
    cat("  Median:", stats::median(object$WIDTH[aligned_rows], na.rm = TRUE), "\n")
    cat("  Max:", max(object$WIDTH[aligned_rows], na.rm = TRUE), "\n")
  }

  proximal_rows <- object$ALIGNMENT_STATUS == "proximal" & is.finite(object$DISTANCE_TO_MERGED)
  if (any(proximal_rows)) {
    cat("\nProximal distances (bp):\n")
    cat("  Min:", min(object$DISTANCE_TO_MERGED[proximal_rows], na.rm = TRUE), "\n")
    cat("  Mean:", round(mean(object$DISTANCE_TO_MERGED[proximal_rows], na.rm = TRUE), 2), "\n")
    cat("  Median:", stats::median(object$DISTANCE_TO_MERGED[proximal_rows], na.rm = TRUE), "\n")
    cat("  Max:", max(object$DISTANCE_TO_MERGED[proximal_rows], na.rm = TRUE), "\n")
  }

  invisible(object)
}
