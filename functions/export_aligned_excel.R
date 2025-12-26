#' Export aligned loci to Excel with visual alignment summary
#'
#' Produces an Excel workbook that juxtaposes the original FUMA loci, the merged
#' locus interval, and the reported loci that align (or fall within the
#' `max_distance` window) for each merged locus. The visualization mirrors the
#' previous refinement report but keeps merged boundaries intact—no splitting or
#' extension beyond the observed merged locus span. Optionally annotates each
#' FUMA block with the maximum `-log10(p)` observed in matching GWAS summary
#' statistics.
#'
#' @param aligned_loci data.frame returned by `align_with_reported()`
#' @param output_file character. Output Excel file path (required)
#' @param gwas_list optional named list of `GWAS` objects returned by
#'   `read_gwas()`. Names must match the FUMA source names used in the alignment
#'   (e.g., "Hisp", "HispCHI"). When supplied, the exporter computes the
#'   maximum `-log10(p)` for each FUMA locus and writes it inside that timeline
#'   block. Leave `NULL` (default) to skip the annotation.
#' @export
export_aligned_to_excel <- function(aligned_loci, output_file, gwas_list = NULL) {
  if (missing(output_file)) {
    stop("output_file is required. Please provide an Excel file path (e.g., 'aligned_loci.xlsx')")
  }
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("openxlsx2 package is required. Install with: install.packages('openxlsx2')")
  }
  if (!is.data.frame(aligned_loci) || !inherits(aligned_loci, "AlignedLoci")) {
    stop("aligned_loci must be the data frame produced by align_with_reported()")
  }
  if (nrow(aligned_loci) == 0) {
    stop("aligned_loci is empty. Run align_with_reported() with keep_unmatched = TRUE or provide reported loci.")
  }

  `%||%` <- function(x, y) if (!is.null(x)) x else y
  sanitize_source_name <- function(x) {
    if (is.null(x)) {
      return(NA_character_)
    }
    x <- trimws(as.character(x))
    x <- sub("\\.[^.]*$", "", x)
    x
  }

  gwas_lookup <- list()
  if (!is.null(gwas_list)) {
    if (inherits(gwas_list, "GWAS")) {
      inferred_name <- gwas_list$name %||% gwas_list$metadata$file_name
      gwas_list <- structure(list(gwas_list), names = if (nzchar(inferred_name)) inferred_name else "GWAS")
    }
    if (!is.list(gwas_list) || length(gwas_list) == 0) {
      stop("gwas_list must be NULL, a GWAS object, or a named list of GWAS objects returned by read_gwas().")
    }
    if (is.null(names(gwas_list)) || any(!nzchar(trimws(names(gwas_list))))) {
      stop("gwas_list must be a *named* list so sources can be matched.")
    }

    source_names <- sanitize_source_name(names(gwas_list))
    if (any(!nzchar(source_names))) {
      stop("All GWAS datasets must have non-empty names after sanitizing (remove extensions if needed).")
    }
    keep_idx <- !duplicated(source_names)
    if (any(!keep_idx)) {
      warning("Duplicate GWAS source names detected after sanitizing; keeping first instance of each.")
    }
    gwas_list <- gwas_list[keep_idx]
    source_names <- source_names[keep_idx]

    gwas_lookup <- lapply(gwas_list, function(obj) {
      if (!is.list(obj) || is.null(obj$data)) {
        stop("Each element of gwas_list must be a GWAS object produced by read_gwas().")
      }
      data <- obj$data
      required_cols <- c("CHR", "POS", "P")
      missing_cols <- setdiff(required_cols, colnames(data))
      if (length(missing_cols) > 0) {
        stop(sprintf("GWAS dataset '%s' is missing required columns: %s",
                     obj$name %||% obj$metadata$file_name,
                     paste(missing_cols, collapse = ", ")))
      }
      data$CHR <- as.character(data$CHR)
      data$POS <- as.numeric(data$POS)
      data$P <- as.numeric(data$P)
      data$NEG_LOG10_P <- -log10(pmax(data$P, .Machine$double.xmin))
      data <- data[is.finite(data$POS) & is.finite(data$NEG_LOG10_P),
                   c("CHR", "POS", "P", "NEG_LOG10_P"), drop = FALSE]
      data
    })
    names(gwas_lookup) <- source_names
  }

  compute_max_neg_log10 <- function(source_name, chr, start_pos, end_pos) {
    if (length(gwas_lookup) == 0) {
      return(NA_real_)
    }
    if (is.na(source_name) || !nzchar(source_name)) {
      return(NA_real_)
    }
    key <- sanitize_source_name(source_name)
    if (!key %in% names(gwas_lookup)) {
      return(NA_real_)
    }
    gwas_df <- gwas_lookup[[key]]
    if (!is.finite(start_pos) || !is.finite(end_pos) || end_pos < start_pos) {
      return(NA_real_)
    }
    chr_str <- as.character(chr)
    chr_df <- gwas_df[gwas_df$CHR == chr_str, , drop = FALSE]
    if (nrow(chr_df) == 0) {
      return(NA_real_)
    }
    chr_df <- chr_df[chr_df$POS >= start_pos & chr_df$POS <= end_pos, , drop = FALSE]
    if (nrow(chr_df) == 0) {
      return(NA_real_)
    }
    max_val <- max(chr_df$NEG_LOG10_P, na.rm = TRUE)
    if (!is.finite(max_val)) {
      return(NA_real_)
    }
    max_val
  }

  aligned_loci$MERGED_LOCUS_ID <- as.character(aligned_loci$MERGED_LOCUS_ID)
  aligned_loci$ORIGINAL_LOCI <- as.character(aligned_loci$ORIGINAL_LOCI)
  aligned_loci$SOURCES <- as.character(aligned_loci$SOURCES)

  all_sources <- unique(trimws(unlist(strsplit(aligned_loci$SOURCES, ";"))))
  all_sources <- all_sources[nzchar(all_sources)]
  source_colors_hex <- character(0)
  if (length(all_sources) > 0) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("RColorBrewer package is required to colour FUMA sources. Install with: install.packages('RColorBrewer')")
    }
    palette_colors <- if (length(all_sources) <= 8) {
      RColorBrewer::brewer.pal(max(3, length(all_sources)), "Set2")
    } else {
      grDevices::rainbow(length(all_sources))
    }
    source_colors_hex <- vapply(seq_along(all_sources), function(i) {
      rgb_vals <- grDevices::col2rgb(palette_colors[i])
      sprintf("FF%02X%02X%02X", rgb_vals[1], rgb_vals[2], rgb_vals[3])
    }, character(1L))
    names(source_colors_hex) <- all_sources
  }

  status_colors <- list(
    source = "FFFFD700",
    merged = "FF2F4F4F",
    aligned = "FF90EE90",
    proximal = "FFFFA500",
    unmatched = "FFD3D3D3",
    tuning = "FFADD8E6"
  )

  excel_rows <- list()
  merge_ranges <- list()
  unique_merged <- unique(aligned_loci$MERGED_LOCUS_ID)

  for (merged_id in unique_merged) {
    block <- aligned_loci[aligned_loci$MERGED_LOCUS_ID == merged_id, , drop = FALSE]
    block <- block[order(factor(block$ALIGNMENT_STATUS, levels = c("aligned", "proximal", "unmatched"), ordered = TRUE),
                         block$START,
                         na.last = TRUE), ]

    merged_chr <- block$CHR[1]
    merged_start <- block$MERGED_START[1]
    merged_end <- block$MERGED_END[1]
    merged_width <- merged_end - merged_start
    if (!is.finite(merged_width) || merged_width <= 0) {
      merged_width <- 1
    }

    start_row_idx <- length(excel_rows) + 1

    original_loci_string <- block$ORIGINAL_LOCI[1]
    if (!is.na(original_loci_string) && nzchar(original_loci_string)) {
      for (locus_key in trimws(strsplit(original_loci_string, ";")[[1]])) {
        matches <- regexec("([^(]+)\\(([^:]+):([^:]+):([^)]+)\\)", locus_key)
        tokens <- regmatches(locus_key, matches)[[1]]
        if (length(tokens) == 5) {
          locus_id <- tokens[2]
          source_name <- tokens[3]
          locus_start <- suppressWarnings(as.numeric(tokens[4]))
          locus_end <- suppressWarnings(as.numeric(tokens[5]))
          locus_width <- if (is.finite(locus_start) && is.finite(locus_end)) locus_end - locus_start + 1 else NA_real_
          max_neg_log10 <- compute_max_neg_log10(source_name, merged_chr, locus_start, locus_end)

          excel_rows[[length(excel_rows) + 1]] <- data.frame(
            MERGED_LOCUS_ID = merged_id,
            ALIGNMENT_ID = locus_id,
            CHR = merged_chr,
            START = locus_start,
            END = locus_end,
            WIDTH = locus_width,
            REPORTED_LOCUS = source_name,
            REPORTED_SOURCE = "FUMA",
            ALIGNMENT_STATUS = "source",
            DISTANCE_TO_MERGED = NA_real_,
            N_ORIGINAL_LOCI = block$N_ORIGINAL_LOCI[1],
            SOURCES = block$SOURCES[1],
            MAX_NEG_LOG10_P = max_neg_log10,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    excel_rows[[length(excel_rows) + 1]] <- data.frame(
      MERGED_LOCUS_ID = merged_id,
      ALIGNMENT_ID = paste0(merged_id, "_MERGED"),
      CHR = merged_chr,
      START = merged_start,
      END = merged_end,
      WIDTH = merged_width,
      REPORTED_LOCUS = "MERGED",
      REPORTED_SOURCE = "MERGED",
      ALIGNMENT_STATUS = "merged",
      DISTANCE_TO_MERGED = NA_real_,
      N_ORIGINAL_LOCI = block$N_ORIGINAL_LOCI[1],
      SOURCES = block$SOURCES[1],
      MAX_NEG_LOG10_P = NA_real_,
      stringsAsFactors = FALSE
    )

    for (r in seq_len(nrow(block))) {
      row <- block[r, ]
      row_start <- row$START
      row_end <- row$END
      row_width <- row$WIDTH

      if (!is.finite(row_start) || !is.finite(row_end) || row$ALIGNMENT_STATUS == "unmatched") {
        row_start <- merged_start
        row_end <- merged_end
        row_width <- merged_width
      }

      excel_rows[[length(excel_rows) + 1]] <- data.frame(
        MERGED_LOCUS_ID = merged_id,
        ALIGNMENT_ID = row$ALIGNED_LOCUS_ID,
        CHR = merged_chr,
        START = row_start,
        END = row_end,
        WIDTH = row_width,
        REPORTED_LOCUS = if (!is.na(row$REPORTED_LOCUS)) row$REPORTED_LOCUS else "N/A",
        REPORTED_SOURCE = if (!is.na(row$REPORTED_SOURCE)) row$REPORTED_SOURCE else "Reported",
        ALIGNMENT_STATUS = row$ALIGNMENT_STATUS,
        DISTANCE_TO_MERGED = row$DISTANCE_TO_MERGED,
        N_ORIGINAL_LOCI = row$N_ORIGINAL_LOCI,
        SOURCES = row$SOURCES,
        MAX_NEG_LOG10_P = NA_real_,
        stringsAsFactors = FALSE
      )
    }

    excel_rows[[length(excel_rows) + 1]] <- data.frame(
      MERGED_LOCUS_ID = merged_id,
      ALIGNMENT_ID = paste0(merged_id, "_TUNING"),
      CHR = merged_chr,
      START = NA_real_,
      END = NA_real_,
      WIDTH = NA_real_,
      REPORTED_LOCUS = paste0(merged_id, " tuning"),
      REPORTED_SOURCE = "Manual",
      ALIGNMENT_STATUS = "tuning",
      DISTANCE_TO_MERGED = NA_real_,
      N_ORIGINAL_LOCI = block$N_ORIGINAL_LOCI[1],
      SOURCES = block$SOURCES[1],
      MAX_NEG_LOG10_P = NA_real_,
      stringsAsFactors = FALSE
    )

    end_row_idx <- length(excel_rows)
    merge_ranges[[length(merge_ranges) + 1]] <- list(
      start = start_row_idx,
      end = end_row_idx,
      merged_start = merged_start,
      merged_end = merged_end
    )
  }

  final_data <- do.call(rbind, excel_rows)

  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("Aligned_Loci")
  wb$add_data(sheet = "Aligned_Loci", x = final_data, start_row = 1, start_col = 1)

  n_vis_cols <- 20
  vis_start_col <- ncol(final_data) + 1

  compute_row_bin_values <- function(source_name, chr, row_start, row_end,
                                     merged_start, merged_end) {
    if (length(gwas_lookup) == 0) {
      return(rep(NA_real_, n_vis_cols))
    }
    if (is.na(source_name) || !nzchar(source_name)) {
      return(rep(NA_real_, n_vis_cols))
    }
    key <- sanitize_source_name(source_name)
    if (!key %in% names(gwas_lookup)) {
      return(rep(NA_real_, n_vis_cols))
    }
    if (!is.finite(row_start) || !is.finite(row_end) || row_end < row_start) {
      return(rep(NA_real_, n_vis_cols))
    }
    chr_df <- gwas_lookup[[key]]
    chr_str <- as.character(chr)
    chr_df <- chr_df[chr_df$CHR == chr_str, , drop = FALSE]
    if (nrow(chr_df) == 0) {
      return(rep(NA_real_, n_vis_cols))
    }
    merged_width <- merged_end - merged_start
    if (!is.finite(merged_width) || merged_width <= 0) {
      merged_width <- 1
    }
    bin_step <- merged_width / n_vis_cols
    bin_vals <- rep(NA_real_, n_vis_cols)

    chr_df <- chr_df[chr_df$POS >= row_start & chr_df$POS <= row_end, , drop = FALSE]
    if (nrow(chr_df) == 0) {
      return(bin_vals)
    }

    bin_indices <- floor((chr_df$POS - merged_start) / bin_step) + 1
    bin_indices <- pmin(pmax(bin_indices, 1L), n_vis_cols)
    for (i in seq_len(nrow(chr_df))) {
      idx <- bin_indices[i]
      bin_vals[idx] <- max(bin_vals[idx], chr_df$NEG_LOG10_P[i], na.rm = TRUE)
    }

    bin_vals
  }

  for (block_idx in seq_along(merge_ranges)) {
    range_info <- merge_ranges[[block_idx]]
    merged_start <- range_info$merged_start
    merged_end <- range_info$merged_end
    merged_width <- merged_end - merged_start
    if (!is.finite(merged_width) || merged_width <= 0) {
      merged_width <- 1
    }

    start_excel_row <- range_info$start + 1
    end_excel_row <- range_info$end + 1
    base_fill <- if (block_idx %% 2 == 1) "FFEEEEEE" else "FFFFFFFF"

    for (offset in 0:(end_excel_row - start_excel_row)) {
      excel_row <- start_excel_row + offset
      data_row <- range_info$start + offset
      if (data_row > nrow(final_data)) {
        next
      }

      status <- final_data$ALIGNMENT_STATUS[data_row]
      source_name <- final_data$REPORTED_LOCUS[data_row]
      row_chr <- final_data$CHR[data_row]
      start_pos <- final_data$START[data_row]
      end_pos <- final_data$END[data_row]

      row_range <- paste0("A", excel_row, ":", openxlsx2::int2col(ncol(final_data)), excel_row)
      if (!is.na(status) && status %in% status_colors && status != "source") {
        wb$add_fill(sheet = "Aligned_Loci", dims = row_range,
                   color = openxlsx2::wb_color(hex = status_colors[[status]]))
        if (status == "merged") {
          wb$add_font(sheet = "Aligned_Loci", dims = row_range, bold = TRUE)
        }
      } else if (base_fill != "FFFFFFFF") {
        wb$add_fill(sheet = "Aligned_Loci", dims = row_range,
                   color = openxlsx2::wb_color(hex = base_fill))
      }
      if (status == "merged") {
        wb$add_font(sheet = "Aligned_Loci", dims = row_range, bold = TRUE)
      }

      if (!is.finite(start_pos) || !is.finite(end_pos)) {
        next
      }

      rel_start <- (start_pos - merged_start) / merged_width
      rel_end <- (end_pos - merged_start) / merged_width
      rel_start <- max(0, min(1, rel_start))
      rel_end <- max(0, min(1, rel_end))

      scaled_start <- vis_start_col + floor(rel_start * n_vis_cols)
      scaled_end <- vis_start_col + ceiling(rel_end * n_vis_cols) - 1
      if (scaled_start > scaled_end) {
        scaled_end <- scaled_start
      }
      scaled_start <- max(scaled_start, vis_start_col)
      scaled_start <- min(scaled_start, vis_start_col + n_vis_cols - 1)
      scaled_end <- min(scaled_end, vis_start_col + n_vis_cols - 1)

      fill_hex <- if (status == "source" && source_name %in% names(source_colors_hex)) {
        source_colors_hex[[source_name]]
      } else if (status %in% names(status_colors)) {
        status_colors[[status]]
      } else {
        status_colors$aligned
      }

      row_is_fuma <- identical(status, "source") && identical(final_data$REPORTED_SOURCE[data_row], "FUMA")
      row_bin_values <- NULL
      row_bin_max <- NA_real_
      if (row_is_fuma) {
        row_bin_values <- compute_row_bin_values(
          source_name = source_name,
          chr = row_chr,
          row_start = start_pos,
          row_end = end_pos,
          merged_start = merged_start,
          merged_end = merged_end
        )
        if (!all(is.na(row_bin_values))) {
          row_bin_max <- max(row_bin_values, na.rm = TRUE)
        }
      }

      target_cols <- scaled_start:scaled_end
      for (col_idx in target_cols) {
        cell <- paste0(openxlsx2::int2col(col_idx), excel_row)
        wb$add_fill(sheet = "Aligned_Loci", dims = cell,
                    color = openxlsx2::wb_color(hex = fill_hex))
        if (row_is_fuma && !is.null(row_bin_values)) {
          bin_idx <- col_idx - vis_start_col + 1
          if (bin_idx >= 1 && bin_idx <= length(row_bin_values)) {
            cell_value <- row_bin_values[bin_idx]
            if (!is.na(cell_value)) {
              val_text <- as.character(round(cell_value))
              wb$add_data(sheet = "Aligned_Loci", start_row = excel_row, start_col = col_idx,
                          x = val_text, col_names = FALSE)
              wb$add_cell_style(sheet = "Aligned_Loci",
                                 dims = cell,
                                 horizontal = "center", vertical = "center")
              if (!is.na(row_bin_max) && abs(cell_value - row_bin_max) <= 1e-6) {
                wb$add_font(sheet = "Aligned_Loci", dims = cell,
                            color = openxlsx2::wb_color(hex = "FFFF0000"), bold = TRUE)
              }
            }
          }
        }
      }
    }
  }

  wb$set_col_widths(sheet = "Aligned_Loci", cols = 1:ncol(final_data), widths = "auto")
  wb$set_col_widths(sheet = "Aligned_Loci",
                    cols = vis_start_col:(vis_start_col + n_vis_cols - 1), widths = 2)

  wb$add_worksheet("Legend")
  legend_data <- data.frame(
    ALIGNMENT_STATUS = c("source", "merged", "aligned", "proximal", "unmatched", "tuning"),
    DESCRIPTION = c(
      "Original FUMA loci that contributed to the merged locus",
      "Merged locus boundaries (unchanged)",
      "Reported loci overlapping the merged locus",
      "Reported loci within the specified max_distance window",
      "Merged loci without matching reported loci",
      "Blank row reserved for manual boundary tuning"
    ),
    stringsAsFactors = FALSE
  )
  wb$add_data(sheet = "Legend", x = legend_data, start_row = 1, start_col = 1)

  for (i in seq_len(nrow(legend_data))) {
    status <- legend_data$ALIGNMENT_STATUS[i]
    if (status %in% names(status_colors)) {
      cell <- paste0("A", i + 1)
      wb$add_fill(sheet = "Legend", dims = cell,
                 color = openxlsx2::wb_color(hex = status_colors[[status]]))
    }
  }

  if (length(source_colors_hex) > 0) {
    wb$add_worksheet("Sources")
    source_data <- data.frame(
      SOURCE = names(source_colors_hex),
      stringsAsFactors = FALSE
    )
    wb$add_data(sheet = "Sources", x = source_data, start_row = 1, start_col = 1)
    for (i in seq_len(nrow(source_data))) {
      cell <- paste0("A", i + 1)
      wb$add_fill(sheet = "Sources", dims = cell,
                 color = openxlsx2::wb_color(hex = source_colors_hex[[source_data$SOURCE[i]]]))
    }
    wb$set_col_widths(sheet = "Sources", cols = 1, widths = "auto")
  }

  wb$save(output_file)

  message(sprintf("Excel file saved to: %s", output_file))
  message(sprintf("  Merged loci reported: %d", length(unique_merged)))
  message(sprintf("  Alignment rows exported: %d", nrow(final_data)))

  invisible(output_file)
}

#' Export aligned loci to Excel for round-trip manual refinement
#'
#' Creates an Excel workbook similar to `export_aligned_to_excel()` but adds
#' per-merged-locus manual refinement rows and machine-readable metadata so the
#' workbook can be edited and imported back into R.
#'
#' Differences vs `export_aligned_to_excel()`:
#' - Adds a TUNING row (place `*` in bins) and a NAMING row (names separated by `;`)
#'   for each `MERGED_LOCUS_ID` block.
#' - Adds `Meta` sheet (mapping bins -> genomic coordinates) and `Parameters` sheet.
#' - Uses a display window of MERGED_START/END ± `flank_bp` so stars can extend
#'   beyond the merged bounds.
#'
#' @param aligned_loci data.frame returned by `align_with_reported()`
#' @param output_file character. Output Excel file path (required)
#' @param gwas_list optional named list of `GWAS` objects returned by `read_gwas()`.
#' @param flank_bp numeric. Number of base-pairs to extend the display window on
#'   both sides of the merged locus. Default: 20000.
#' @param n_bins integer. Number of visualization bins/columns. Default: 20.
#' @export
export_aligned_to_excel_roundtrip <- function(aligned_loci,
                                              output_file,
                                              gwas_list = NULL,
                                              flank_bp = 20000,
                                              n_bins = 20) {
  if (missing(output_file)) {
    stop("output_file is required. Please provide an Excel file path (e.g., 'aligned_loci_roundtrip.xlsx')")
  }
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("openxlsx2 package is required. Install with: install.packages('openxlsx2')")
  }
  if (!is.data.frame(aligned_loci) || !inherits(aligned_loci, "AlignedLoci")) {
    stop("aligned_loci must be the data frame produced by align_with_reported()")
  }
  if (nrow(aligned_loci) == 0) {
    stop("aligned_loci is empty. Run align_with_reported() with keep_unmatched = TRUE or provide reported loci.")
  }
  flank_bp <- suppressWarnings(as.numeric(flank_bp))
  if (!is.finite(flank_bp) || flank_bp < 0) {
    stop("flank_bp must be a non-negative number (base-pairs).")
  }
  n_bins <- suppressWarnings(as.integer(n_bins))
  if (!is.finite(n_bins) || n_bins < 5) {
    stop("n_bins must be an integer >= 5")
  }

  protect_sheet_best_effort <- function(wb, sheet) {
    # openxlsx2 has had a few naming differences across versions; try a few.
    try(wb$protect_worksheet(sheet = sheet), silent = TRUE)
    try(wb$worksheet_protect(sheet = sheet), silent = TRUE)
    try(wb$sheet_protect(sheet = sheet), silent = TRUE)
    invisible(NULL)
  }

  `%||%` <- function(x, y) if (!is.null(x)) x else y
  sanitize_source_name <- function(x) {
    if (is.null(x)) {
      return(NA_character_)
    }
    x <- trimws(as.character(x))
    x <- sub("\\.[^.]*$", "", x)
    x
  }

  gwas_lookup <- list()
  if (!is.null(gwas_list)) {
    if (inherits(gwas_list, "GWAS")) {
      inferred_name <- gwas_list$name %||% gwas_list$metadata$file_name
      gwas_list <- structure(list(gwas_list), names = if (nzchar(inferred_name)) inferred_name else "GWAS")
    }
    if (!is.list(gwas_list) || length(gwas_list) == 0) {
      stop("gwas_list must be NULL, a GWAS object, or a named list of GWAS objects returned by read_gwas().")
    }
    if (is.null(names(gwas_list)) || any(!nzchar(trimws(names(gwas_list))))) {
      stop("gwas_list must be a *named* list so sources can be matched.")
    }

    source_names <- sanitize_source_name(names(gwas_list))
    if (any(!nzchar(source_names))) {
      stop("All GWAS datasets must have non-empty names after sanitizing (remove extensions if needed).")
    }
    keep_idx <- !duplicated(source_names)
    if (any(!keep_idx)) {
      warning("Duplicate GWAS source names detected after sanitizing; keeping first instance of each.")
    }
    gwas_list <- gwas_list[keep_idx]
    source_names <- source_names[keep_idx]

    gwas_lookup <- lapply(gwas_list, function(obj) {
      if (!is.list(obj) || is.null(obj$data)) {
        stop("Each element of gwas_list must be a GWAS object produced by read_gwas().")
      }
      data <- obj$data
      required_cols <- c("CHR", "POS", "P")
      missing_cols <- setdiff(required_cols, colnames(data))
      if (length(missing_cols) > 0) {
        stop(sprintf("GWAS dataset '%s' is missing required columns: %s",
                     obj$name %||% obj$metadata$file_name,
                     paste(missing_cols, collapse = ", ")))
      }
      data$CHR <- as.character(data$CHR)
      data$POS <- as.numeric(data$POS)
      data$P <- as.numeric(data$P)
      data$NEG_LOG10_P <- -log10(pmax(data$P, .Machine$double.xmin))
      data <- data[is.finite(data$POS) & is.finite(data$NEG_LOG10_P),
                   c("CHR", "POS", "P", "NEG_LOG10_P"), drop = FALSE]
      data
    })
    names(gwas_lookup) <- source_names
  }

  compute_max_neg_log10 <- function(source_name, chr, start_pos, end_pos) {
    if (length(gwas_lookup) == 0) {
      return(NA_real_)
    }
    if (is.na(source_name) || !nzchar(source_name)) {
      return(NA_real_)
    }
    key <- sanitize_source_name(source_name)
    if (!key %in% names(gwas_lookup)) {
      return(NA_real_)
    }
    gwas_df <- gwas_lookup[[key]]
    if (!is.finite(start_pos) || !is.finite(end_pos) || end_pos < start_pos) {
      return(NA_real_)
    }
    chr_str <- as.character(chr)
    chr_df <- gwas_df[gwas_df$CHR == chr_str, , drop = FALSE]
    if (nrow(chr_df) == 0) {
      return(NA_real_)
    }
    chr_df <- chr_df[chr_df$POS >= start_pos & chr_df$POS <= end_pos, , drop = FALSE]
    if (nrow(chr_df) == 0) {
      return(NA_real_)
    }
    max_val <- max(chr_df$NEG_LOG10_P, na.rm = TRUE)
    if (!is.finite(max_val)) {
      return(NA_real_)
    }
    max_val
  }

  aligned_loci$MERGED_LOCUS_ID <- as.character(aligned_loci$MERGED_LOCUS_ID)
  aligned_loci$ORIGINAL_LOCI <- as.character(aligned_loci$ORIGINAL_LOCI)
  aligned_loci$SOURCES <- as.character(aligned_loci$SOURCES)

  all_sources <- unique(trimws(unlist(strsplit(aligned_loci$SOURCES, ";"))))
  all_sources <- all_sources[nzchar(all_sources)]
  source_colors_hex <- character(0)
  if (length(all_sources) > 0) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("RColorBrewer package is required to colour FUMA sources. Install with: install.packages('RColorBrewer')")
    }
    palette_colors <- if (length(all_sources) <= 8) {
      RColorBrewer::brewer.pal(max(3, length(all_sources)), "Set2")
    } else {
      grDevices::rainbow(length(all_sources))
    }
    source_colors_hex <- vapply(seq_along(all_sources), function(i) {
      rgb_vals <- grDevices::col2rgb(palette_colors[i])
      sprintf("FF%02X%02X%02X", rgb_vals[1], rgb_vals[2], rgb_vals[3])
    }, character(1L))
    names(source_colors_hex) <- all_sources
  }

  status_colors <- list(
    source = "FFFFD700",
    merged = "FF2F4F4F",
    aligned = "FF90EE90",
    proximal = "FFFFA500",
    unmatched = "FFD3D3D3",
    tuning = "FFADD8E6",
    naming = "FFADD8E6"
  )

  excel_rows <- list()
  merge_ranges <- list()
  meta_rows <- list()
  unique_merged <- unique(aligned_loci$MERGED_LOCUS_ID)

  for (merged_id in unique_merged) {
    block <- aligned_loci[aligned_loci$MERGED_LOCUS_ID == merged_id, , drop = FALSE]
    block <- block[order(factor(block$ALIGNMENT_STATUS, levels = c("aligned", "proximal", "unmatched"), ordered = TRUE),
                         block$START,
                         na.last = TRUE), ]

    merged_chr <- block$CHR[1]
    merged_start <- block$MERGED_START[1]
    merged_end <- block$MERGED_END[1]
    merged_width <- merged_end - merged_start
    if (!is.finite(merged_width) || merged_width <= 0) {
      merged_width <- 1
    }

    display_start <- suppressWarnings(as.numeric(merged_start) - flank_bp)
    display_end <- suppressWarnings(as.numeric(merged_end) + flank_bp)
    if (!is.finite(display_start) || !is.finite(display_end) || display_end <= display_start) {
      display_start <- merged_start
      display_end <- merged_end
    }
    if (is.finite(display_start) && display_start < 0) {
      display_start <- 0
    }

    start_row_idx <- length(excel_rows) + 1

    original_loci_string <- block$ORIGINAL_LOCI[1]
    if (!is.na(original_loci_string) && nzchar(original_loci_string)) {
      for (locus_key in trimws(strsplit(original_loci_string, ";")[[1]])) {
        matches <- regexec("([^(]+)\\(([^:]+):([^:]+):([^)]+)\\)", locus_key)
        tokens <- regmatches(locus_key, matches)[[1]]
        if (length(tokens) == 5) {
          locus_id <- tokens[2]
          source_name <- tokens[3]
          locus_start <- suppressWarnings(as.numeric(tokens[4]))
          locus_end <- suppressWarnings(as.numeric(tokens[5]))
          locus_width <- if (is.finite(locus_start) && is.finite(locus_end)) locus_end - locus_start + 1 else NA_real_
          max_neg_log10 <- compute_max_neg_log10(source_name, merged_chr, locus_start, locus_end)

          excel_rows[[length(excel_rows) + 1]] <- data.frame(
            MERGED_LOCUS_ID = merged_id,
            ALIGNMENT_ID = locus_id,
            CHR = merged_chr,
            START = locus_start,
            END = locus_end,
            WIDTH = locus_width,
            REPORTED_LOCUS = source_name,
            REPORTED_SOURCE = "FUMA",
            ALIGNMENT_STATUS = "source",
            DISTANCE_TO_MERGED = NA_real_,
            N_ORIGINAL_LOCI = block$N_ORIGINAL_LOCI[1],
            SOURCES = block$SOURCES[1],
            MAX_NEG_LOG10_P = max_neg_log10,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    excel_rows[[length(excel_rows) + 1]] <- data.frame(
      MERGED_LOCUS_ID = merged_id,
      ALIGNMENT_ID = paste0(merged_id, "_MERGED"),
      CHR = merged_chr,
      START = merged_start,
      END = merged_end,
      WIDTH = merged_width,
      REPORTED_LOCUS = "MERGED",
      REPORTED_SOURCE = "MERGED",
      ALIGNMENT_STATUS = "merged",
      DISTANCE_TO_MERGED = NA_real_,
      N_ORIGINAL_LOCI = block$N_ORIGINAL_LOCI[1],
      SOURCES = block$SOURCES[1],
      MAX_NEG_LOG10_P = NA_real_,
      stringsAsFactors = FALSE
    )

    for (r in seq_len(nrow(block))) {
      row <- block[r, ]
      row_start <- row$START
      row_end <- row$END
      row_width <- row$WIDTH

      if (!is.finite(row_start) || !is.finite(row_end) || row$ALIGNMENT_STATUS == "unmatched") {
        row_start <- merged_start
        row_end <- merged_end
        row_width <- merged_width
      }

      excel_rows[[length(excel_rows) + 1]] <- data.frame(
        MERGED_LOCUS_ID = merged_id,
        ALIGNMENT_ID = row$ALIGNED_LOCUS_ID,
        CHR = merged_chr,
        START = row_start,
        END = row_end,
        WIDTH = row_width,
        REPORTED_LOCUS = if (!is.na(row$REPORTED_LOCUS)) row$REPORTED_LOCUS else "N/A",
        REPORTED_SOURCE = if (!is.na(row$REPORTED_SOURCE)) row$REPORTED_SOURCE else "Reported",
        ALIGNMENT_STATUS = row$ALIGNMENT_STATUS,
        DISTANCE_TO_MERGED = row$DISTANCE_TO_MERGED,
        N_ORIGINAL_LOCI = row$N_ORIGINAL_LOCI,
        SOURCES = row$SOURCES,
        MAX_NEG_LOG10_P = NA_real_,
        stringsAsFactors = FALSE
      )
    }

    excel_rows[[length(excel_rows) + 1]] <- data.frame(
      MERGED_LOCUS_ID = merged_id,
      ALIGNMENT_ID = paste0(merged_id, "_TUNING"),
      CHR = merged_chr,
      START = NA_real_,
      END = NA_real_,
      WIDTH = NA_real_,
      REPORTED_LOCUS = paste0(merged_id, " tuning"),
      REPORTED_SOURCE = "Manual",
      ALIGNMENT_STATUS = "tuning",
      DISTANCE_TO_MERGED = NA_real_,
      N_ORIGINAL_LOCI = block$N_ORIGINAL_LOCI[1],
      SOURCES = block$SOURCES[1],
      MAX_NEG_LOG10_P = NA_real_,
      stringsAsFactors = FALSE
    )

    excel_rows[[length(excel_rows) + 1]] <- data.frame(
      MERGED_LOCUS_ID = merged_id,
      ALIGNMENT_ID = paste0(merged_id, "_NAMING"),
      CHR = merged_chr,
      START = NA_real_,
      END = NA_real_,
      WIDTH = NA_real_,
      REPORTED_LOCUS = paste0(merged_id, " names"),
      REPORTED_SOURCE = "Manual",
      ALIGNMENT_STATUS = "naming",
      DISTANCE_TO_MERGED = NA_real_,
      N_ORIGINAL_LOCI = block$N_ORIGINAL_LOCI[1],
      SOURCES = block$SOURCES[1],
      MAX_NEG_LOG10_P = NA_real_,
      stringsAsFactors = FALSE
    )

    end_row_idx <- length(excel_rows)
    merge_ranges[[length(merge_ranges) + 1]] <- list(
      start = start_row_idx,
      end = end_row_idx,
      display_start = display_start,
      display_end = display_end
    )

    meta_rows[[length(meta_rows) + 1]] <- data.frame(
      MERGED_LOCUS_ID = merged_id,
      CHR = merged_chr,
      MERGED_START = merged_start,
      MERGED_END = merged_end,
      DISPLAY_START = display_start,
      DISPLAY_END = display_end,
      N_BINS = n_bins,
      stringsAsFactors = FALSE
    )
  }

  final_data <- do.call(rbind, excel_rows)
  meta_data <- do.call(rbind, meta_rows)

  n_vis_cols <- n_bins
  vis_start_col <- ncol(final_data) + 1
  meta_data$VIS_COL_START <- vis_start_col
  meta_data$VIS_COL_END <- vis_start_col + n_vis_cols - 1

  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("Aligned_Loci")
  wb$add_data(sheet = "Aligned_Loci", x = final_data, start_row = 1, start_col = 1)

  compute_row_bin_values <- function(source_name, chr, row_start, row_end,
                                     display_start, display_end) {
    if (length(gwas_lookup) == 0) {
      return(rep(NA_real_, n_vis_cols))
    }
    if (is.na(source_name) || !nzchar(source_name)) {
      return(rep(NA_real_, n_vis_cols))
    }
    key <- sanitize_source_name(source_name)
    if (!key %in% names(gwas_lookup)) {
      return(rep(NA_real_, n_vis_cols))
    }
    if (!is.finite(row_start) || !is.finite(row_end) || row_end < row_start) {
      return(rep(NA_real_, n_vis_cols))
    }
    chr_df <- gwas_lookup[[key]]
    chr_str <- as.character(chr)
    chr_df <- chr_df[chr_df$CHR == chr_str, , drop = FALSE]
    if (nrow(chr_df) == 0) {
      return(rep(NA_real_, n_vis_cols))
    }

    display_width <- display_end - display_start
    if (!is.finite(display_width) || display_width <= 0) {
      display_width <- 1
    }
    bin_step <- display_width / n_vis_cols
    bin_vals <- rep(NA_real_, n_vis_cols)

    chr_df <- chr_df[chr_df$POS >= row_start & chr_df$POS <= row_end, , drop = FALSE]
    if (nrow(chr_df) == 0) {
      return(bin_vals)
    }

    bin_indices <- floor((chr_df$POS - display_start) / bin_step) + 1
    bin_indices <- pmin(pmax(bin_indices, 1L), n_vis_cols)
    for (i in seq_len(nrow(chr_df))) {
      idx <- bin_indices[i]
      bin_vals[idx] <- max(bin_vals[idx], chr_df$NEG_LOG10_P[i], na.rm = TRUE)
    }

    bin_vals
  }

  for (block_idx in seq_along(merge_ranges)) {
    range_info <- merge_ranges[[block_idx]]
    display_start <- range_info$display_start
    display_end <- range_info$display_end
    display_width <- display_end - display_start
    if (!is.finite(display_width) || display_width <= 0) {
      display_width <- 1
    }

    start_excel_row <- range_info$start + 1
    end_excel_row <- range_info$end + 1
    base_fill <- if (block_idx %% 2 == 1) "FFEEEEEE" else "FFFFFFFF"

    for (offset in 0:(end_excel_row - start_excel_row)) {
      excel_row <- start_excel_row + offset
      data_row <- range_info$start + offset
      if (data_row > nrow(final_data)) {
        next
      }

      status <- final_data$ALIGNMENT_STATUS[data_row]
      source_name <- final_data$REPORTED_LOCUS[data_row]
      row_chr <- final_data$CHR[data_row]
      start_pos <- final_data$START[data_row]
      end_pos <- final_data$END[data_row]

      row_range <- paste0("A", excel_row, ":", openxlsx2::int2col(ncol(final_data)), excel_row)
      if (!is.na(status) && status %in% status_colors && status != "source") {
        wb$add_fill(sheet = "Aligned_Loci", dims = row_range,
                    color = openxlsx2::wb_color(hex = status_colors[[status]]))
        if (status == "merged") {
          wb$add_font(sheet = "Aligned_Loci", dims = row_range, bold = TRUE)
        }
      } else if (base_fill != "FFFFFFFF") {
        wb$add_fill(sheet = "Aligned_Loci", dims = row_range,
                    color = openxlsx2::wb_color(hex = base_fill))
      }
      if (status == "merged") {
        wb$add_font(sheet = "Aligned_Loci", dims = row_range, bold = TRUE)
      }

      if (!is.finite(start_pos) || !is.finite(end_pos)) {
        next
      }

      rel_start <- (start_pos - display_start) / display_width
      rel_end <- (end_pos - display_start) / display_width
      rel_start <- max(0, min(1, rel_start))
      rel_end <- max(0, min(1, rel_end))

      scaled_start <- vis_start_col + floor(rel_start * n_vis_cols)
      scaled_end <- vis_start_col + ceiling(rel_end * n_vis_cols) - 1
      if (scaled_start > scaled_end) {
        scaled_end <- scaled_start
      }
      scaled_start <- max(scaled_start, vis_start_col)
      scaled_start <- min(scaled_start, vis_start_col + n_vis_cols - 1)
      scaled_end <- min(scaled_end, vis_start_col + n_vis_cols - 1)

      fill_hex <- if (status == "source" && source_name %in% names(source_colors_hex)) {
        source_colors_hex[[source_name]]
      } else if (status %in% names(status_colors)) {
        status_colors[[status]]
      } else {
        status_colors$aligned
      }

      row_is_fuma <- identical(status, "source") && identical(final_data$REPORTED_SOURCE[data_row], "FUMA")
      row_bin_values <- NULL
      row_bin_max <- NA_real_
      if (row_is_fuma) {
        row_bin_values <- compute_row_bin_values(
          source_name = source_name,
          chr = row_chr,
          row_start = start_pos,
          row_end = end_pos,
          display_start = display_start,
          display_end = display_end
        )
        if (!all(is.na(row_bin_values))) {
          row_bin_max <- max(row_bin_values, na.rm = TRUE)
        }
      }

      target_cols <- scaled_start:scaled_end
      for (col_idx in target_cols) {
        cell <- paste0(openxlsx2::int2col(col_idx), excel_row)
        wb$add_fill(sheet = "Aligned_Loci", dims = cell,
                    color = openxlsx2::wb_color(hex = fill_hex))
        if (row_is_fuma && !is.null(row_bin_values)) {
          bin_idx <- col_idx - vis_start_col + 1
          if (bin_idx >= 1 && bin_idx <= length(row_bin_values)) {
            cell_value <- row_bin_values[bin_idx]
            if (!is.na(cell_value)) {
              val_text <- as.character(round(cell_value))
              wb$add_data(sheet = "Aligned_Loci", start_row = excel_row, start_col = col_idx,
                          x = val_text, col_names = FALSE)
              wb$add_cell_style(sheet = "Aligned_Loci",
                                dims = cell,
                                horizontal = "center", vertical = "center")
              if (!is.na(row_bin_max) && abs(cell_value - row_bin_max) <= 1e-6) {
                wb$add_font(sheet = "Aligned_Loci", dims = cell,
                            color = openxlsx2::wb_color(hex = "FFFF0000"), bold = TRUE)
              }
            }
          }
        }
      }
    }
  }

  wb$set_col_widths(sheet = "Aligned_Loci", cols = 1:ncol(final_data), widths = "auto")
  wb$set_col_widths(sheet = "Aligned_Loci",
                    cols = vis_start_col:(vis_start_col + n_vis_cols - 1), widths = 2)

  wb$add_worksheet("Parameters")
  params <- data.frame(
    PARAMETER = c("FLANK_BP", "N_BINS", "STAR_SYMBOL", "DEFAULT_NO_STARS_BEHAVIOR"),
    VALUE = c(as.character(flank_bp), as.character(n_bins), "*", "keep_merged_as_is"),
    stringsAsFactors = FALSE
  )
  wb$add_data(sheet = "Parameters", x = params, start_row = 1, start_col = 1)
  wb$set_col_widths(sheet = "Parameters", cols = 1:2, widths = "auto")

  wb$add_worksheet("Metadata")
  export_time <- format(Sys.time(), tz = "UTC", usetz = TRUE)
  gwas_names <- if (!is.null(gwas_list)) {
    if (inherits(gwas_list, "GWAS")) {
      gwas_list$name %||% gwas_list$metadata$file_name
    } else {
      paste(names(gwas_list), collapse = ";")
    }
  } else {
    ""
  }
  metadata <- data.frame(
    KEY = c(
      "EXPORT_UTC",
      "OUTPUT_FILE",
      "N_MERGED_LOCI",
      "N_ALIGNMENT_ROWS",
      "FLANK_BP",
      "N_BINS",
      "GWAS_SOURCES",
      "R_VERSION",
      "PLATFORM"
    ),
    VALUE = c(
      export_time,
      normalizePath(output_file, winslash = "/", mustWork = FALSE),
      as.character(length(unique_merged)),
      as.character(nrow(final_data)),
      as.character(flank_bp),
      as.character(n_bins),
      as.character(gwas_names),
      R.version.string,
      paste(R.version$platform, R.version$arch)
    ),
    stringsAsFactors = FALSE
  )
  wb$add_data(sheet = "Metadata", x = metadata, start_row = 1, start_col = 1)
  wb$set_col_widths(sheet = "Metadata", cols = 1:2, widths = "auto")

  wb$add_worksheet("Meta")
  wb$add_data(sheet = "Meta", x = meta_data, start_row = 1, start_col = 1)
  wb$set_col_widths(sheet = "Meta", cols = 1:ncol(meta_data), widths = "auto")

  # Protect machine-readable sheets so accidental edits are harder.
  protect_sheet_best_effort(wb, "Metadata")
  protect_sheet_best_effort(wb, "Parameters")
  protect_sheet_best_effort(wb, "Meta")

  # Best-effort hide Meta (different openxlsx2 versions name this differently)
  try(wb$set_sheet_visibility(sheet = "Meta", value = "hidden"), silent = TRUE)
  try(wb$sheet_visibility(sheet = "Meta", value = "hidden"), silent = TRUE)
  try(wb$sheet_hide(sheet = "Meta"), silent = TRUE)

  wb$add_worksheet("Legend")
  legend_data <- data.frame(
    ALIGNMENT_STATUS = c("source", "merged", "aligned", "proximal", "unmatched", "tuning", "naming"),
    DESCRIPTION = c(
      "Original FUMA loci that contributed to the merged locus",
      "Merged locus boundaries (unchanged)",
      "Reported loci overlapping the merged locus",
      "Reported loci within the specified max_distance window",
      "Merged loci without matching reported loci",
      "Manual row: place '*' under bins to define final loci",
      "Manual row: enter final locus name(s) separated by ';'"
    ),
    stringsAsFactors = FALSE
  )
  wb$add_data(sheet = "Legend", x = legend_data, start_row = 1, start_col = 1)
  for (i in seq_len(nrow(legend_data))) {
    status <- legend_data$ALIGNMENT_STATUS[i]
    if (status %in% names(status_colors)) {
      cell <- paste0("A", i + 1)
      wb$add_fill(sheet = "Legend", dims = cell,
                  color = openxlsx2::wb_color(hex = status_colors[[status]]))
    }
  }
  wb$set_col_widths(sheet = "Legend", cols = 1:2, widths = "auto")

  if (length(source_colors_hex) > 0) {
    wb$add_worksheet("Sources")
    source_data <- data.frame(
      SOURCE = names(source_colors_hex),
      stringsAsFactors = FALSE
    )
    wb$add_data(sheet = "Sources", x = source_data, start_row = 1, start_col = 1)
    for (i in seq_len(nrow(source_data))) {
      cell <- paste0("A", i + 1)
      wb$add_fill(sheet = "Sources", dims = cell,
                  color = openxlsx2::wb_color(hex = source_colors_hex[[source_data$SOURCE[i]]]))
    }
    wb$set_col_widths(sheet = "Sources", cols = 1, widths = "auto")
  }

  wb$save(output_file)

  message(sprintf("Excel file saved to: %s", output_file))
  message(sprintf("  Merged loci reported: %d", length(unique_merged)))
  message(sprintf("  Alignment rows exported: %d", nrow(final_data)))
  message(sprintf("  Display flank (bp): %s", format(flank_bp, scientific = FALSE)))
  message(sprintf("  Visual bins: %d", n_bins))

  invisible(output_file)
}
