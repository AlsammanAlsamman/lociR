#' Export aligned loci to an interactive locus-at-a-time HTML editor (wizard)
#'
#' Generates a single self-contained HTML file (offline) that shows one
#' MERGED_LOCUS_ID at a time with enhanced visualizations.
#'
#' @param aligned_loci data.frame returned by `align_with_reported()` (class 'AlignedLoci')
#' @param output_html character. Output HTML filepath (e.g., 'outputs/loci_wizard.html')
#' @param flank_bp numeric. Base-pairs to extend display window for visualization.
#' @param n_bins integer. Stored in exported TSV for reproducibility.
#' @param gwas_list Optional list of GWAS objects to include p-value tracks
#' @param window_bp Window size for p-value binning
#' @param show_genes logical. Whether to show gene annotations if available
#' @param gene_padding numeric. Padding around genes for display
#'
#' @return Named list with path: `list(html = <path>)`.
#' @export
export_aligned_to_html_locus_wizard <- function(aligned_loci,
                                                output_html,
                                                flank_bp = 20000,
                                                n_bins = 20,
                                                gwas_list = NULL,
                                                window_bp = 20000,
                                                show_genes = TRUE,
                                                gene_padding = 5000) {
  
  # Input validation (keep original code)
  if (missing(output_html) || !nzchar(trimws(output_html))) {
    stop("output_html is required (e.g., 'outputs/loci_wizard.html').")
  }
  if (!is.data.frame(aligned_loci) || !inherits(aligned_loci, "AlignedLoci")) {
    stop("aligned_loci must be the data frame produced by align_with_reported() (class 'AlignedLoci').")
  }
  if (nrow(aligned_loci) == 0) {
    stop("aligned_loci is empty.")
  }
  
  flank_bp <- suppressWarnings(as.numeric(flank_bp))
  if (!is.finite(flank_bp) || flank_bp < 0) {
    stop("flank_bp must be a non-negative number (base-pairs).")
  }
  n_bins <- suppressWarnings(as.integer(n_bins))
  if (!is.finite(n_bins) || n_bins < 5) {
    stop("n_bins must be an integer >= 5")
  }
  
  if (!requireNamespace("htmltools", quietly = TRUE)) {
    stop("htmltools package is required. Install with: install.packages('htmltools')")
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite package is required. Install with: install.packages('jsonlite')")
  }
  
  window_bp <- suppressWarnings(as.integer(window_bp))
  if (!is.finite(window_bp) || window_bp < 1000) {
    stop("window_bp must be an integer >= 1000 (base-pairs).")
  }
  
  required_cols <- c(
    "MERGED_LOCUS_ID",
    "CHR",
    "MERGED_START",
    "MERGED_END",
    "ALIGNED_LOCUS_ID",
    "START",
    "END",
    "ALIGNMENT_STATUS",
    "N_ORIGINAL_LOCI",
    "SOURCES"
  )
  missing_cols <- setdiff(required_cols, colnames(aligned_loci))
  if (length(missing_cols) > 0) {
    stop("aligned_loci is missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  df <- aligned_loci
  df$DISPLAY_START <- pmax(0, df$MERGED_START - flank_bp)
  df$DISPLAY_END <- df$MERGED_END + flank_bp
  
  # ORIGINAL_LOCI is produced by merge_fuma_loci(); used to colour by GWAS/FUMA source.
  if (!"ORIGINAL_LOCI" %in% colnames(df)) {
    df$ORIGINAL_LOCI <- NA_character_
  }
  
  # Ensure editable columns exist (per row); wizard exports locus-level TSV.
  if (!"REFINED_START" %in% colnames(df)) df$REFINED_START <- df$MERGED_START
  if (!"REFINED_END" %in% colnames(df)) df$REFINED_END <- df$MERGED_END
  if (!"LOCUS_NAME" %in% colnames(df)) df$LOCUS_NAME <- ""
  if (!"KEEP" %in% colnames(df)) df$KEEP <- TRUE
  
  # Sort for stable ordering
  df <- df[order(df$CHR, df$MERGED_START, df$MERGED_END, df$MERGED_LOCUS_ID), , drop = FALSE]
  
  # Build a locus-level summary table (one row per MERGED_LOCUS_ID)
  loci_ids <- unique(df$MERGED_LOCUS_ID)
  loci_rows <- vector("list", length(loci_ids))
  
  for (i in seq_along(loci_ids)) {
    id <- loci_ids[[i]]
    sub <- df[df$MERGED_LOCUS_ID == id, , drop = FALSE]
    loci_rows[[i]] <- data.frame(
      MERGED_LOCUS_ID = id,
      CHR = sub$CHR[[1]],
      MERGED_START = sub$MERGED_START[[1]],
      MERGED_END = sub$MERGED_END[[1]],
      DISPLAY_START = sub$DISPLAY_START[[1]],
      DISPLAY_END = sub$DISPLAY_END[[1]],
      REFINED_START = sub$REFINED_START[[1]],
      REFINED_END = sub$REFINED_END[[1]],
      LOCUS_NAME = sub$LOCUS_NAME[[1]],
      KEEP = isTRUE(sub$KEEP[[1]]),
      N_ALIGNED_ROWS = nrow(sub),
      N_ORIGINAL_LOCI = sub$N_ORIGINAL_LOCI[[1]],
      SOURCES = sub$SOURCES[[1]],
      FLANK_BP = flank_bp,
      N_BINS = n_bins,
      SPLITS = "",
      stringsAsFactors = FALSE
    )
  }
  
  loci_df <- do.call(rbind, loci_rows)
  
  # Parse ORIGINAL_LOCI into a per-locus list of intervals.
  parse_original_loci <- function(original_loci_string) {
    if (is.na(original_loci_string) || !nzchar(original_loci_string)) {
      return(data.frame())
    }
    parts <- trimws(strsplit(original_loci_string, ";", fixed = TRUE)[[1]])
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) {
      return(data.frame())
    }
    rows <- vector("list", length(parts))
    keep <- 0L
    for (j in seq_along(parts)) {
      locus_str <- parts[[j]]
      pattern <- "^(.+)\\((.+):(\\d+):(\\d+)\\)$"
      if (!grepl(pattern, locus_str)) {
        next
      }
      locus_id <- sub(pattern, "\\1", locus_str)
      source <- sub(pattern, "\\2", locus_str)
      start <- suppressWarnings(as.numeric(sub(pattern, "\\3", locus_str)))
      end <- suppressWarnings(as.numeric(sub(pattern, "\\4", locus_str)))
      if (!is.finite(start) || !is.finite(end)) {
        next
      }
      keep <- keep + 1L
      rows[[keep]] <- data.frame(
        LOCUS_ID = locus_id,
        SOURCE = source,
        START = start,
        END = end,
        stringsAsFactors = FALSE
      )
    }
    if (keep == 0L) {
      return(data.frame())
    }
    do.call(rbind, rows[seq_len(keep)])
  }
  
  orig_rows <- list()
  for (i in seq_len(nrow(loci_df))) {
    merged_id <- loci_df$MERGED_LOCUS_ID[[i]]
    block <- df[df$MERGED_LOCUS_ID == merged_id, , drop = FALSE]
    original_loci_string <- block$ORIGINAL_LOCI[[1]]
    parsed <- parse_original_loci(original_loci_string)
    if (nrow(parsed) == 0) {
      next
    }
    parsed$MERGED_LOCUS_ID <- merged_id
    orig_rows[[length(orig_rows) + 1L]] <- parsed
  }
  orig_df <- if (length(orig_rows) > 0) do.call(rbind, orig_rows) else data.frame()
  
  # Colour palette for original loci sources.
  unique_sources <- character(0)
  if (nrow(orig_df) > 0 && "SOURCE" %in% colnames(orig_df)) {
    unique_sources <- unique(as.character(orig_df$SOURCE))
    unique_sources <- unique_sources[nzchar(unique_sources)]
  }
  source_colors <- list()
  if (length(unique_sources) > 0) {
    if (exists("hcl.colors", where = asNamespace("grDevices"), inherits = FALSE)) {
      cols <- grDevices::hcl.colors(length(unique_sources), palette = "Dark 3")
    } else {
      cols <- grDevices::rainbow(length(unique_sources))
    }
    source_colors <- stats::setNames(as.list(cols), unique_sources)
  }
  
  # Optional: build per-window max -log10(p) track(s) using gwas_list.
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
                   c("CHR", "POS", "NEG_LOG10_P"), drop = FALSE]
      data
    })
    names(gwas_lookup) <- source_names
  }
  
  pvals_df <- data.frame()
  if (length(gwas_lookup) > 0 && nrow(loci_df) > 0) {
    track_sources <- if (length(unique_sources) > 0) unique_sources else names(gwas_lookup)
    track_sources <- sanitize_source_name(track_sources)
    track_sources <- intersect(track_sources, names(gwas_lookup))
    
    pval_rows <- list()
    for (i in seq_len(nrow(loci_df))) {
      merged_id <- loci_df$MERGED_LOCUS_ID[[i]]
      chr <- as.character(loci_df$CHR[[i]])
      ds <- suppressWarnings(as.numeric(loci_df$DISPLAY_START[[i]]))
      de <- suppressWarnings(as.numeric(loci_df$DISPLAY_END[[i]]))
      if (!is.finite(ds) || !is.finite(de) || de < ds) {
        next
      }
      win_starts <- seq(ds, de, by = window_bp)
      win_ends <- pmin(win_starts + window_bp - 1, de)
      for (src in track_sources) {
        gdf <- gwas_lookup[[src]]
        chr_df <- gdf[gdf$CHR == chr, , drop = FALSE]
        if (nrow(chr_df) == 0) {
          next
        }
        chr_df <- chr_df[chr_df$POS >= ds & chr_df$POS <= de, , drop = FALSE]
        if (nrow(chr_df) == 0) {
          next
        }
        for (w in seq_along(win_starts)) {
          ws <- win_starts[[w]]
          we <- win_ends[[w]]
          sub <- chr_df[chr_df$POS >= ws & chr_df$POS <= we, , drop = FALSE]
          maxp <- if (nrow(sub) > 0) max(sub$NEG_LOG10_P, na.rm = TRUE) else NA_real_
          if (!is.finite(maxp)) {
            next
          }
          pval_rows[[length(pval_rows) + 1L]] <- data.frame(
            MERGED_LOCUS_ID = merged_id,
            SOURCE = src,
            WIN_START = ws,
            WIN_END = we,
            MAX_NEG_LOG10_P = maxp,
            stringsAsFactors = FALSE
          )
        }
      }
    }
    pvals_df <- if (length(pval_rows) > 0) do.call(rbind, pval_rows) else data.frame()
  }
  
  # Build segments for reported loci
  seg_keep <- c("MERGED_LOCUS_ID", "ALIGNED_LOCUS_ID", "START", "END", "ALIGNMENT_STATUS")
  if ("REPORTED_SOURCE" %in% colnames(df)) {
    seg_keep <- c(seg_keep, "REPORTED_SOURCE")
  }
  seg_df <- df[, seg_keep, drop = FALSE]
  seg_df$START <- suppressWarnings(as.numeric(seg_df$START))
  seg_df$END <- suppressWarnings(as.numeric(seg_df$END))
  
  # Create JSON data for JavaScript
  loci_js <- jsonlite::toJSON(loci_df, dataframe = "rows", na = "null", auto_unbox = TRUE)
  seg_js <- jsonlite::toJSON(seg_df, dataframe = "rows", na = "null", auto_unbox = TRUE)
  orig_js <- jsonlite::toJSON(orig_df, dataframe = "rows", na = "null", auto_unbox = TRUE)
  colors_js <- jsonlite::toJSON(source_colors, auto_unbox = TRUE)
  pvals_js <- jsonlite::toJSON(pvals_df, dataframe = "rows", na = "null", auto_unbox = TRUE)
  
  # FIX: Define download_name here
  download_name <- paste0(tools::file_path_sans_ext(basename(output_html)), ".tsv")
  download_name_js <- jsonlite::toJSON(download_name, auto_unbox = TRUE)

  # Load large JS/CSS blobs from helper file to keep this function readable.
  if (!exists("locus_wizard_build_js", mode = "function") || !exists("locus_wizard_build_css", mode = "function")) {
    tpl_path <- file.path("functions", "locus_wizard_templates.R")
    if (file.exists(tpl_path)) {
      source(tpl_path)
    }
  }
  if (!exists("locus_wizard_build_js", mode = "function") || !exists("locus_wizard_build_css", mode = "function")) {
    stop("Missing locus wizard template helpers. Source 'functions/locus_wizard_templates.R' or run load_all_functions().")
  }

  js <- locus_wizard_build_js(
    download_name_js = download_name_js,
    loci_js = loci_js,
    seg_js = seg_js,
    orig_js = orig_js,
    colors_js = colors_js,
    pvals_js = pvals_js
  )
  css <- htmltools::HTML(locus_wizard_build_css())
  
  # Create HTML page
  page <- htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title("Enhanced Locus Wizard"),
      htmltools::tags$style(css)
    ),
    htmltools::tags$body(
      htmltools::tags$h1("Enhanced Locus Wizard"),
      htmltools::tags$pre(
        id = "js_error",
        "(JavaScript error output will appear here.)"
      ),
      htmltools::tags$p(class = "small",
        "Edit REFINED_START/REFINED_END/LOCUS_NAME/KEEP. ",
        "Drag on the plot to select a region and click 'Add split'. ",
        "Use 'Save (Download TSV)' to export edits."
      ),
      htmltools::tags$div(class = "controls",
        htmltools::tags$button(id = "prev", class = "btn", "Prev"),
        htmltools::tags$button(id = "next", class = "btn", "Next"),
        htmltools::tags$button(id = "save", class = "btn", "Save (Download TSV)"),
        htmltools::tags$span(id = "counter", class = "small"),
        htmltools::tags$span(id = "dirty", class = "badge badge-ok", "Saved")
      ),
      htmltools::tags$div(class = "locus-info-grid",
        htmltools::tags$div(class = "info-item",
          htmltools::tags$span(class = "info-label", "Chromosome: "),
          htmltools::tags$span(id = "chr", class = "info-value")
        ),
        htmltools::tags$div(class = "info-item",
          htmltools::tags$span(class = "info-label", "Locus name: "),
          htmltools::tags$span(id = "locus_name_hdr", class = "info-value")
        )
      ),
      htmltools::tags$div(class = "panel left",
        htmltools::tags$svg(id = "viz", viewBox = "0 0 900 220"),
        htmltools::tags$div(style = "margin-top: 10px;",
          htmltools::tags$button(id = "add_split", class = "btn", "Add split from selection"),
          htmltools::tags$span(class = "small", " (drag on visualization to select region)")
        )
      ),
      htmltools::tags$div(class = "row",
        htmltools::tags$div(class = "panel left",
          htmltools::tags$h3("Edits"),
          htmltools::tags$label("REFINED_START"),
          htmltools::tags$input(id = "refined_start", type = "number", step = "1"),
          htmltools::tags$label("REFINED_END"),
          htmltools::tags$input(id = "refined_end", type = "number", step = "1"),
          htmltools::tags$label("LOCUS_NAME"),
          htmltools::tags$input(id = "locus_name", type = "text"),
          htmltools::tags$label("KEEP"),
          htmltools::tags$input(id = "keep", type = "checkbox")
        ),
        htmltools::tags$div(class = "panel right",
          htmltools::tags$h3("Splits"),
          htmltools::tags$ul(id = "splits"),
          htmltools::tags$p(class = "small",
            "Splits are saved as start-end pairs in the exported TSV column 'SPLITS'."
          )
        )
      ),
      htmltools::tags$script(htmltools::HTML(js))
    )
  )
  
  # Save the HTML file
  out_dir <- dirname(output_html)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  htmltools::save_html(page, file = output_html)
  
  message(sprintf("Enhanced locus wizard HTML saved to: %s", output_html))
  invisible(list(html = output_html))
}