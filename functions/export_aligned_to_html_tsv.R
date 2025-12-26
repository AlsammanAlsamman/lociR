#' Export aligned loci to TSV + HTML (round-trip via TSV)
#'
#' This exporter replaces the Excel-based review workflow with a simpler hybrid:
#'
#' - A TSV file is the single source of truth for manual edits.
#' - A static HTML file is generated only for visualization/review.
#'
#' The TSV contains one row per aligned locus record (as produced by
#' `align_with_reported()`), plus user-editable refinement columns:
#' `REFINED_START`, `REFINED_END`, `LOCUS_NAME`, and `KEEP`.
#'
#' The HTML renders the TSV as an interactive table (DT) grouped by
#' `MERGED_LOCUS_ID`, with row coloring by `ALIGNMENT_STATUS`. If `gwas_list`
#' is provided, the HTML also includes per-merged-locus plots showing
#' GWAS points in the display window.
#'
#' @param aligned_loci data.frame returned by `align_with_reported()`
#' @param output_prefix,out_prefix character. Prefix used to write `<prefix>.tsv` and
#'   `<prefix>.html`.
#' @param gwas_list optional named list of `GWAS` objects returned by `read_gwas()`.
#'   Used only for optional plots.
#' @param flank_bp numeric. Base-pairs to extend the display window beyond the
#'   merged locus boundaries on both sides. Default: 20000.
#' @param n_bins integer. Number of visualization bins across the display window.
#'   Saved into the TSV for reproducibility. Default: 20.
#'
#' @return Named list with paths: `list(tsv = <path>, html = <path>)`.
#' @export
export_aligned_to_html_tsv <- function(aligned_loci,
                                       output_prefix = NULL,
                                       out_prefix = NULL,
                                       gwas_list = NULL,
                                       flank_bp = 20000,
                                       n_bins = 20) {
  if (is.null(output_prefix) || !nzchar(trimws(output_prefix))) {
    output_prefix <- out_prefix
  }
  if (is.null(output_prefix) || !nzchar(trimws(output_prefix))) {
    stop("output_prefix/out_prefix is required (e.g., 'outputs/aligned_review').")
  }
  if (!is.data.frame(aligned_loci) || !inherits(aligned_loci, "AlignedLoci")) {
    stop("aligned_loci must be the data frame produced by align_with_reported() (class 'AlignedLoci').")
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

  if (!requireNamespace("readr", quietly = TRUE)) {
    stop("readr package is required. Install with: install.packages('readr')")
  }
  if (!requireNamespace("DT", quietly = TRUE)) {
    stop("DT package is required. Install with: install.packages('DT')")
  }
  if (!requireNamespace("htmltools", quietly = TRUE)) {
    stop("htmltools package is required. Install with: install.packages('htmltools')")
  }

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  required_cols <- c(
    "MERGED_LOCUS_ID",
    "CHR",
    "MERGED_START",
    "MERGED_END",
    "ALIGNED_LOCUS_ID",
    "START",
    "END",
    "ALIGNMENT_STATUS",
    "REPORTED_LOCUS",
    "REPORTED_SOURCE",
    "DISTANCE_TO_MERGED",
    "N_ORIGINAL_LOCI",
    "SOURCES"
  )
  missing_cols <- setdiff(required_cols, colnames(aligned_loci))
  if (length(missing_cols) > 0) {
    stop("aligned_loci is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  tsv_path <- paste0(output_prefix, ".tsv")
  html_path <- paste0(output_prefix, ".html")

  out_dir <- dirname(tsv_path)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  # Build round-trip TSV (authoritative editable artifact)
  export_df <- aligned_loci[, required_cols, drop = FALSE]
  export_df$CHR <- as.character(export_df$CHR)
  export_df$MERGED_START <- suppressWarnings(as.numeric(export_df$MERGED_START))
  export_df$MERGED_END <- suppressWarnings(as.numeric(export_df$MERGED_END))
  export_df$START <- suppressWarnings(as.numeric(export_df$START))
  export_df$END <- suppressWarnings(as.numeric(export_df$END))

  export_df$DISPLAY_START <- pmax(0, export_df$MERGED_START - flank_bp)
  export_df$DISPLAY_END <- export_df$MERGED_END + flank_bp

  # User-editable refinement columns
  export_df$REFINED_START <- NA_real_
  export_df$REFINED_END <- NA_real_
  export_df$LOCUS_NAME <- NA_character_
  export_df$KEEP <- TRUE

  # Reorder exactly as requested
  export_df <- export_df[, c(
    "MERGED_LOCUS_ID",
    "CHR",
    "MERGED_START",
    "MERGED_END",
    "DISPLAY_START",
    "DISPLAY_END",
    "ALIGNED_LOCUS_ID",
    "START",
    "END",
    "ALIGNMENT_STATUS",
    "REPORTED_LOCUS",
    "REPORTED_SOURCE",
    "DISTANCE_TO_MERGED",
    "N_ORIGINAL_LOCI",
    "SOURCES",
    "REFINED_START",
    "REFINED_END",
    "LOCUS_NAME",
    "KEEP"
  ), drop = FALSE]

  readr::write_tsv(export_df, file = tsv_path, na = "NA")

  # --- HTML generation (read-only visualization)
  # Load the TSV we just wrote (TSV is the source of truth)
  tsv_df <- readr::read_tsv(tsv_path, show_col_types = FALSE, progress = FALSE)

  # Row colors by status
  status_colors <- c(
    aligned = "#d7f5d7",
    proximal = "#fff0cc",
    unmatched = "#f2f2f2",
    source = "#fff7b3",
    merged = "#d9d9d9"
  )

  # Ensure group column is first so rowGroup can use dataSrc = 0
  if (!identical(names(tsv_df)[1], "MERGED_LOCUS_ID")) {
    tsv_df <- tsv_df[, c("MERGED_LOCUS_ID", setdiff(names(tsv_df), "MERGED_LOCUS_ID"))]
  }

  dt <- DT::datatable(
    tsv_df,
    rownames = FALSE,
    extensions = c("RowGroup"),
    options = list(
      pageLength = 50,
      scrollX = TRUE,
      rowGroup = list(dataSrc = 0),
      order = list(list(0, "asc"))
    )
  )

  dt <- DT::formatStyle(
    dt,
    "ALIGNMENT_STATUS",
    target = "row",
    backgroundColor = DT::styleEqual(names(status_colors), unname(status_colors))
  )

  # Optional per-locus plots (static SVG), using gwas_list if provided
  plot_section <- NULL
  if (!is.null(gwas_list)) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("ggplot2 package is required for plotting. Install with: install.packages('ggplot2')")
    }

    sanitize_source_name <- function(x) {
      if (is.null(x)) return(NA_character_)
      x <- trimws(as.character(x))
      sub("\\.[^.]*$", "", x)
    }

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
    keep_idx <- !duplicated(source_names)
    gwas_list <- gwas_list[keep_idx]
    names(gwas_list) <- source_names[keep_idx]

    gwas_lookup <- lapply(gwas_list, function(obj) {
      if (!is.list(obj) || is.null(obj$data)) {
        stop("Each element of gwas_list must be a GWAS object produced by read_gwas().")
      }
      data <- obj$data
      required_gwas_cols <- c("CHR", "POS", "P")
      miss <- setdiff(required_gwas_cols, colnames(data))
      if (length(miss) > 0) {
        stop("GWAS dataset is missing required columns: ", paste(miss, collapse = ", "))
      }
      data$CHR <- as.character(data$CHR)
      data$POS <- as.numeric(data$POS)
      data$P <- as.numeric(data$P)
      data$NEG_LOG10_P <- -log10(pmax(data$P, .Machine$double.xmin))
      data <- data[is.finite(data$POS) & is.finite(data$NEG_LOG10_P), c("CHR", "POS", "P", "NEG_LOG10_P"), drop = FALSE]
      data
    })

    merged_ids <- unique(tsv_df$MERGED_LOCUS_ID)
    max_plots <- 30L
    plot_ids <- merged_ids[seq_len(min(length(merged_ids), max_plots))]

    plot_tags <- list(
      htmltools::tags$h2("Optional GWAS plots"),
      htmltools::tags$p(
        "Each plot shows GWAS points in the display window (",
        flank_bp,
        " bp flank). Vertical dashed lines mark MERGED_START/END."
      )
    )

    if (length(merged_ids) > max_plots) {
      plot_tags <- c(plot_tags, list(htmltools::tags$p(
        sprintf("Showing first %d merged loci (of %d total) to keep the HTML size manageable.",
                max_plots, length(merged_ids))
      )))
    }

    make_svg_tag <- function(p) {
      tmp <- tempfile(fileext = ".svg")
      grDevices::svg(filename = tmp, width = 9, height = 3)
      print(p)
      grDevices::dev.off()
      svg_txt <- paste(readLines(tmp, warn = FALSE), collapse = "\n")
      unlink(tmp)
      htmltools::HTML(svg_txt)
    }

    for (merged_id in plot_ids) {
      block <- tsv_df[tsv_df$MERGED_LOCUS_ID == merged_id, , drop = FALSE]
      chr <- as.character(block$CHR[1])
      merged_start <- suppressWarnings(as.numeric(block$MERGED_START[1]))
      merged_end <- suppressWarnings(as.numeric(block$MERGED_END[1]))
      display_start <- suppressWarnings(as.numeric(block$DISPLAY_START[1]))
      display_end <- suppressWarnings(as.numeric(block$DISPLAY_END[1]))

      plot_df <- NULL
      for (src in names(gwas_lookup)) {
        gdf <- gwas_lookup[[src]]
        gdf <- gdf[gdf$CHR == chr & gdf$POS >= display_start & gdf$POS <= display_end, , drop = FALSE]
        if (nrow(gdf) == 0) {
          next
        }
        gdf$SOURCE <- src
        plot_df <- rbind(plot_df, gdf)
      }

      if (is.null(plot_df) || nrow(plot_df) == 0) {
        plot_tags <- c(plot_tags, list(
          htmltools::tags$h3(paste0("MERGED_LOCUS_ID: ", merged_id)),
          htmltools::tags$p("No GWAS points found in display window for the provided gwas_list.")
        ))
        next
      }

      p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = POS, y = NEG_LOG10_P, color = SOURCE)) +
        ggplot2::geom_point(alpha = 0.6, size = 0.8) +
        ggplot2::geom_vline(xintercept = merged_start, linetype = "dashed", color = "black") +
        ggplot2::geom_vline(xintercept = merged_end, linetype = "dashed", color = "black") +
        ggplot2::labs(
          title = paste0("MERGED_LOCUS_ID: ", merged_id, " (chr", chr, ")"),
          x = "Position (bp)",
          y = "-log10(P)"
        ) +
        ggplot2::theme_minimal(base_size = 11) +
        ggplot2::theme(legend.position = "bottom")

      plot_tags <- c(plot_tags, list(
        htmltools::tags$h3(paste0("MERGED_LOCUS_ID: ", merged_id)),
        make_svg_tag(p)
      ))
    }

    plot_section <- htmltools::tagList(plot_tags)
  }

  page <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title("Aligned loci review (TSV + HTML)"),
      htmltools::tags$style(htmltools::HTML(
        "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif; margin: 18px; }\n",
        "h1,h2,h3 { margin-top: 1.2em; }\n",
        ".small { color: #555; }\n"
      ))
    ),
    htmltools::tags$body(
      htmltools::tags$h1("Aligned loci review (TSV + HTML)"),
      htmltools::tags$p(class = "small",
        "TSV is the ONLY editable artifact. HTML is read-only and can be regenerated at any time."),
      htmltools::tags$ul(
        htmltools::tags$li(paste0("TSV: ", normalizePath(tsv_path, winslash = "/", mustWork = FALSE))),
        htmltools::tags$li(paste0("HTML: ", normalizePath(html_path, winslash = "/", mustWork = FALSE))),
        htmltools::tags$li(paste0("Merged loci: ", length(unique(tsv_df$MERGED_LOCUS_ID)))),
        htmltools::tags$li(paste0("Rows: ", nrow(tsv_df))),
        htmltools::tags$li(paste0("Display flank (bp): ", flank_bp, "; bins: ", n_bins))
      ),
      htmltools::tags$h2("Table"),
      dt,
      plot_section
    )
  )

  htmltools::save_html(page, file = html_path)

  message(sprintf("TSV saved to: %s", tsv_path))
  message(sprintf("HTML saved to: %s", html_path))

  invisible(list(tsv = tsv_path, html = html_path))
}
