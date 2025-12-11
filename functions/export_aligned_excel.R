#' Export aligned loci to Excel with visual alignment summary
#'
#' Produces an Excel workbook that juxtaposes the original FUMA loci, the merged
#' locus interval, and the reported loci that align (or fall within the
#' `max_distance` window) for each merged locus. The visualization mirrors the
#' previous refinement report but keeps merged boundaries intact (no splitting or
#' extension).
#'
#' @param aligned_loci data.frame returned by `align_with_reported()`
#' @param output_file character. Output Excel file path (required)
#' @export
export_aligned_to_excel <- function(aligned_loci, output_file) {
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

  aligned_loci$MERGED_LOCUS_ID <- as.character(aligned_loci$MERGED_LOCUS_ID)
  aligned_loci$ORIGINAL_LOCI <- as.character(aligned_loci$ORIGINAL_LOCI)
  aligned_loci$SOURCES <- as.character(aligned_loci$SOURCES)

  all_sources <- unique(unlist(strsplit(aligned_loci$SOURCES, ";")))
  all_sources <- sort(unique(trimws(all_sources[nzchar(all_sources)])))
  if (length(all_sources) > 0) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("RColorBrewer package is required to color FUMA sources. Install with: install.packages('RColorBrewer')")
    }
    if (length(all_sources) <= 8) {
      palette_colors <- RColorBrewer::brewer.pal(max(3, length(all_sources)), "Set2")
    } else {
      palette_colors <- grDevices::rainbow(length(all_sources))
    }
    source_colors <- setNames(palette_colors[seq_along(all_sources)], all_sources)
    source_colors_hex <- vapply(source_colors, function(col) {
      rgb_vals <- grDevices::col2rgb(col)
      sprintf("FF%02X%02X%02X", rgb_vals[1], rgb_vals[2], rgb_vals[3])
    }, character(1))
  } else {
    source_colors_hex <- character(0)
  }

  status_colors <- list(
    source = "FFFFD700",
    merged = "FF2F4F4F",
    aligned = "FF90EE90",
    proximal = "FFFFA500",
    unmatched = "FFD3D3D3"
  )

  excel_rows <- list()
  merge_ranges <- list()
  unique_merged <- unique(aligned_loci$MERGED_LOCUS_ID)

  for (merged_id in unique_merged) {
    subset_rows <- aligned_loci[aligned_loci$MERGED_LOCUS_ID == merged_id, , drop = FALSE]
    subset_rows <- subset_rows[order(factor(subset_rows$ALIGNMENT_STATUS,
                                            levels = c("aligned", "proximal", "unmatched"), ordered = TRUE),
                                     subset_rows$START,
                                     na.last = TRUE), ]

    merged_chr <- subset_rows$CHR[1]
    merged_start <- subset_rows$MERGED_START[1]
    merged_end <- subset_rows$MERGED_END[1]
    merged_width <- subset_rows$MERGED_WIDTH[1]

    start_row_idx <- length(excel_rows) + 1
    first_label <- TRUE

    original_loci_string <- subset_rows$ORIGINAL_LOCI[1]
    if (!is.na(original_loci_string) && nzchar(original_loci_string)) {
      original_keys <- strsplit(original_loci_string, ";")[[1]]
      original_keys <- trimws(original_keys[nzchar(original_keys)])

      for (locus_key in original_keys) {
        matches <- regmatches(locus_key, regexec("([^(]+)\\(([^:]+):([^:]+):([^)]+)\\)", locus_key))
        if (length(matches[[1]]) == 5) {
          locus_id <- matches[[1]][2]
          source_name <- matches[[1]][3]
          locus_start <- as.numeric(matches[[1]][4])
          locus_end <- as.numeric(matches[[1]][5])

          excel_rows[[length(excel_rows) + 1]] <- data.frame(
            MERGED_LOCUS_ID = if (first_label) merged_id else "",
            ALIGNMENT_ID = locus_id,
            CHR = merged_chr,
            START = locus_start,
            END = locus_end,
            WIDTH = locus_end - locus_start + 1,
            REPORTED_LOCUS = source_name,
            REPORTED_SOURCE = "FUMA",
            ALIGNMENT_STATUS = "source",
            DISTANCE_TO_MERGED = NA_real_,
            N_ORIGINAL_LOCI = subset_rows$N_ORIGINAL_LOCI[1],
            SOURCES = subset_rows$SOURCES[1],
            stringsAsFactors = FALSE
          )
          first_label <- FALSE
        }
      }
    }

    if (first_label) {
      excel_rows[[length(excel_rows) + 1]] <- data.frame(
        MERGED_LOCUS_ID = merged_id,
        ALIGNMENT_ID = paste0(merged_id, "_PLACEHOLDER"),
        CHR = merged_chr,
        START = merged_start,
        END = merged_end,
        WIDTH = merged_width,
        REPORTED_LOCUS = "FUMA",
        REPORTED_SOURCE = "FUMA",
        ALIGNMENT_STATUS = "source",
        DISTANCE_TO_MERGED = NA_real_,
        N_ORIGINAL_LOCI = subset_rows$N_ORIGINAL_LOCI[1],
        SOURCES = subset_rows$SOURCES[1],
        stringsAsFactors = FALSE
      )
      first_label <- FALSE
    }

    excel_rows[[length(excel_rows) + 1]] <- data.frame(
      MERGED_LOCUS_ID = "",
      ALIGNMENT_ID = paste0(merged_id, "_MERGED"),
      CHR = merged_chr,
      START = merged_start,
      END = merged_end,
      WIDTH = merged_width,
      REPORTED_LOCUS = "MERGED",
      REPORTED_SOURCE = "MERGED",
      ALIGNMENT_STATUS = "merged",
      DISTANCE_TO_MERGED = NA_real_,
      N_ORIGINAL_LOCI = subset_rows$N_ORIGINAL_LOCI[1],
      SOURCES = subset_rows$SOURCES[1],
      stringsAsFactors = FALSE
    )

    for (r in seq_len(nrow(subset_rows))) {
      row <- subset_rows[r, ]
      status <- row$ALIGNMENT_STATUS
      disp_start <- row$START
      disp_end <- row$END
      disp_width <- row$WIDTH

      if (status == "unmatched") {
        disp_start <- merged_start
        disp_end <- merged_end
        disp_width <- merged_width
      }

      excel_rows[[length(excel_rows) + 1]] <- data.frame(
        MERGED_LOCUS_ID = "",
        ALIGNMENT_ID = row$ALIGNED_LOCUS_ID,
        CHR = merged_chr,
        START = disp_start,
        END = disp_end,
        WIDTH = disp_width,
        REPORTED_LOCUS = if (!is.na(row$REPORTED_LOCUS)) row$REPORTED_LOCUS else "N/A",
        REPORTED_SOURCE = if (!is.na(row$REPORTED_SOURCE)) row$REPORTED_SOURCE else "Reported",
        ALIGNMENT_STATUS = status,
        DISTANCE_TO_MERGED = row$DISTANCE_TO_MERGED,
        N_ORIGINAL_LOCI = row$N_ORIGINAL_LOCI,
        SOURCES = row$SOURCES,
        stringsAsFactors = FALSE
      )
    }

    end_row_idx <- length(excel_rows)
    merge_ranges[[length(merge_ranges) + 1]] <- list(
      start = start_row_idx,
      end = end_row_idx,
      merged_id = merged_id,
      merged_start = merged_start,
      merged_end = merged_end,
      merged_width = merged_width
    )
  }

  final_data <- do.call(rbind, excel_rows)

  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("Aligned_Loci")
  wb$add_data(sheet = "Aligned_Loci", x = final_data, start_row = 1, start_col = 1)

  for (range_info in merge_ranges) {
    start_excel_row <- range_info$start + 1
    end_excel_row <- range_info$end + 1

    if (end_excel_row > start_excel_row) {
      dims <- paste0("A", start_excel_row, ":A", end_excel_row)
      wb$merge_cells(sheet = "Aligned_Loci", dims = dims)
      wb$add_cell_style(sheet = "Aligned_Loci", dims = paste0("A", start_excel_row),
                        horizontal = "center", vertical = "center")
    }
  }

  n_vis_cols <- 20
  vis_start_col <- ncol(final_data) + 1

  for (merge_info in merge_ranges) {
    start_excel_row <- merge_info$start + 1
    end_excel_row <- merge_info$end + 1
    merged_width <- merge_info$merged_width

    for (offset in 0:(end_excel_row - start_excel_row)) {
      excel_row <- start_excel_row + offset
      data_row <- merge_info$start + offset
      if (data_row > nrow(final_data)) {
        next
      }

      status <- final_data$ALIGNMENT_STATUS[data_row]
      source_name <- final_data$REPORTED_LOCUS[data_row]
      start_pos <- final_data$START[data_row]
      end_pos <- final_data$END[data_row]

      if (status == "merged") {
        row_range <- paste0("A", excel_row, ":", openxlsx2::int2col(ncol(final_data)), excel_row)
        wb$add_fill(sheet = "Aligned_Loci", dims = row_range,
                   color = openxlsx2::wb_color(hex = status_colors$merged))
        wb$add_font(sheet = "Aligned_Loci", dims = row_range, bold = TRUE)

        merged_rel_start <- 0
        merged_rel_end <- 1
        if (is.finite(merged_width) && merged_width > 0) {
          merged_rel_start <- (merge_info$merged_start - merge_info$merged_start) / merged_width
          merged_rel_end <- (merge_info$merged_end - merge_info$merged_start) / merged_width
        }

        scaled_start <- vis_start_col + floor(merged_rel_start * n_vis_cols)
        scaled_end <- vis_start_col + ceiling(merged_rel_end * n_vis_cols) - 1
        if (scaled_start > scaled_end) {
          scaled_end <- scaled_start
        }
        scaled_start <- max(scaled_start, vis_start_col)
        scaled_end <- min(scaled_end, vis_start_col + n_vis_cols - 1)

        for (col_idx in scaled_start:scaled_end) {
          cell_ref <- paste0(openxlsx2::int2col(col_idx), excel_row)
          wb$add_fill(sheet = "Aligned_Loci", dims = cell_ref,
                     color = openxlsx2::wb_color(hex = status_colors$merged))
        }
        next
      }

      if (status == "unmatched") {
        row_range <- paste0("A", excel_row, ":", openxlsx2::int2col(ncol(final_data)), excel_row)
        wb$add_fill(sheet = "Aligned_Loci", dims = row_range,
                   color = openxlsx2::wb_color(hex = status_colors$unmatched))
        next
      }

      if (!is.finite(start_pos) || !is.finite(end_pos) || !is.finite(merged_width) || merged_width <= 0) {
        next
      }

      rel_start <- (start_pos - merge_info$merged_start) / merged_width
      rel_end <- (end_pos - merge_info$merged_start) / merged_width
      rel_start <- max(0, min(1, rel_start))
      rel_end <- max(0, min(1, rel_end))

      scaled_start <- vis_start_col + floor(rel_start * n_vis_cols)
      scaled_end <- vis_start_col + ceiling(rel_end * n_vis_cols) - 1
      if (scaled_start > scaled_end) {
        scaled_end <- scaled_start
      }
      scaled_start <- max(scaled_start, vis_start_col)
      scaled_end <- min(scaled_end, vis_start_col + n_vis_cols - 1)

      if (status == "source") {
        if (source_name %in% names(source_colors_hex)) {
          fill_hex <- source_colors_hex[[source_name]]
        } else {
          fill_hex <- status_colors$source
        }
      } else if (status %in% names(status_colors)) {
        fill_hex <- status_colors[[status]]
      } else {
        fill_hex <- status_colors$aligned
      }

      for (col_idx in scaled_start:scaled_end) {
        cell_address <- paste0(openxlsx2::int2col(col_idx), excel_row)
        wb$add_fill(sheet = "Aligned_Loci", dims = cell_address,
                   color = openxlsx2::wb_color(hex = fill_hex))
      }
    }
  }

  wb$set_col_widths(sheet = "Aligned_Loci", cols = 1:ncol(final_data), widths = "auto")
  wb$set_col_widths(sheet = "Aligned_Loci",
                    cols = vis_start_col:(vis_start_col + n_vis_cols - 1), widths = 2)

  wb$add_worksheet("Legend")
  legend_data <- data.frame(
    ALIGNMENT_STATUS = c("source", "merged", "aligned", "proximal", "unmatched"),
    DESCRIPTION = c(
      "Original FUMA loci that contributed to the merged locus",
      "Merged locus boundaries (unchanged)",
      "Reported loci overlapping the merged locus",
      "Reported loci within the specified max_distance window",
      "Merged loci without matching reported loci"
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
