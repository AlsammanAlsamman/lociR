#' Export aligned loci to Excel
#'
#' Generate an Excel report summarizing how merged loci align with reported loci.
#' Each merged locus is presented on a single row with alignment metadata, without
#' splitting or extending the genomic intervals.
#'
#' @param aligned_loci data.frame produced by `align_with_reported()`
#' @param output_file character. Output Excel file path (required)
#' @export
export_aligned_to_excel <- function(aligned_loci, output_file) {
  if (missing(output_file)) {
    stop("output_file is required. Please provide an Excel file path (e.g., 'aligned_loci.xlsx')")
  }

  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("openxlsx2 package is required. Install with: install.packages('openxlsx2')")
  }

  if (!is.data.frame(aligned_loci)) {
    stop("aligned_loci must be a data frame returned by align_with_reported()")
  }

  col_order <- c(
    "MERGED_LOCUS_ID", "ALIGNED_LOCUS_ID", "CHR",
    "MERGED_START", "MERGED_END", "MERGED_WIDTH",
    "ALIGNED_START", "ALIGNED_END", "ALIGNED_WIDTH",
    "REPORTED_START", "REPORTED_END", "REPORTED_LOCI", "REPORTED_SOURCE",
    "ALIGNMENT_STATUS", "N_MATCHED_REPORTED", "DISTANCE_TO_REPORTED",
    "N_ORIGINAL_LOCI", "ORIGINAL_LOCI", "SOURCES"
  )

  available_cols <- intersect(col_order, colnames(aligned_loci))
  extra_cols <- setdiff(colnames(aligned_loci), available_cols)
  data_to_write <- aligned_loci[, c(available_cols, extra_cols), drop = FALSE]

  message(sprintf("Creating Excel file for %d aligned loci ...", nrow(data_to_write)))

  status_palette <- list(
    aligned = "FF9ACD32",   # Yellow green
    proximal = "FFFFA500",  # Orange
    unmatched = "FFB0C4DE"  # Light steel blue
  )

  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("Aligned_Loci")
  wb$add_data(sheet = "Aligned_Loci", x = data_to_write, start_row = 1, start_col = 1)

  if (nrow(data_to_write) > 0) {
    last_col_letter <- openxlsx2::int2col(ncol(data_to_write))

    for (row_idx in seq_len(nrow(data_to_write))) {
      status <- data_to_write$ALIGNMENT_STATUS[row_idx]
      if (!is.null(status) && !is.na(status) && status %in% names(status_palette)) {
        row_range <- paste0("A", row_idx + 1, ":", last_col_letter, row_idx + 1)
        wb$add_fill(
          sheet = "Aligned_Loci",
          dims = row_range,
          color = openxlsx2::wb_color(hex = status_palette[[status]])
        )
      }
    }
  }

  wb$set_col_widths(sheet = "Aligned_Loci", cols = 1:ncol(data_to_write), widths = "auto")

  wb$add_worksheet("Legend")
  legend_data <- data.frame(
    ALIGNMENT_STATUS = names(status_palette),
    DESCRIPTION = c(
      "Reported loci overlap the merged locus",
      "Reported loci are within the max_distance window",
      "No reported loci matched"
    ),
    stringsAsFactors = FALSE
  )
  wb$add_data(sheet = "Legend", x = legend_data, start_row = 1, start_col = 1)

  for (i in seq_len(nrow(legend_data))) {
    status <- legend_data$ALIGNMENT_STATUS[i]
    cell_address <- paste0("A", i + 1)
    wb$add_fill(
      sheet = "Legend",
      dims = cell_address,
      color = openxlsx2::wb_color(hex = status_palette[[status]])
    )
  }
  wb$set_col_widths(sheet = "Legend", cols = 1:2, widths = "auto")

  wb$save(output_file)

  message(sprintf("Excel file saved to: %s", output_file))
  invisible(output_file)
}
