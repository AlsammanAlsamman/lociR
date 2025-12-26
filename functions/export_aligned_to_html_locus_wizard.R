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
  
  # Simplified JavaScript with the enhanced visualization
  js <- paste0(
    "(function(){\n",
    "  // Data\n",
    "  var DOWNLOAD_NAME = ", download_name_js, ";\n",
    "  var LOCI = ", loci_js, ";\n",
    "  var SEGS = ", seg_js, ";\n",
    "  var ORIGS = ", orig_js, ";\n",
    "  var COLORS = ", colors_js, ";\n",
    "  var PVALS = ", pvals_js, ";\n",
    "  var idx = 0;\n",
    "  var dirty = false;\n",
    "  var selection = null;\n",
    "  \n",
    "  // Utility functions\n",
    "  function byId(id) { return document.getElementById(id); }\n",
    "  function showError(msg) {\n",
    "    var box = byId('js_error');\n",
    "    if (box) { box.style.display = 'block'; box.textContent = String(msg); }\n",
    "  }\n",
    "  function escTSV(v) {\n",
    "    if (v === null || v === undefined) return '';\n",
    "    var s = String(v);\n",
    "    s = s.replace(/\\t/g, ' ');\n",
    "    s = s.replace(/\\r/g, ' ');\n",
    "    s = s.replace(/\\n/g, ' ');\n",
    "    return s;\n",
    "  }\n",
    "  function isInteger(n) { return typeof n === 'number' && isFinite(n) && Math.floor(n) === n; }\n",
    "  function splitsToString(arr) {\n",
    "    if (!arr || arr.length === 0) return '';\n",
    "    var out = [];\n",
    "    for (var i = 0; i < arr.length; i++) { out.push(arr[i].start + '-' + arr[i].end); }\n",
    "    return out.join(';');\n",
    "  }\n",
    "  function parseSplits(s) {\n",
    "    var out = [];\n",
    "    if (!s) return out;\n",
    "    var parts = String(s).split(';');\n",
    "    for (var i = 0; i < parts.length; i++) {\n",
    "      var t = parts[i];\n",
    "      if (!t) continue;\n",
    "      t = t.replace(/^\\s+|\\s+$/g, '');\n",
    "      if (t.length === 0) continue;\n",
    "      var p = t.split('-');\n",
    "      if (p.length < 2) continue;\n",
    "      out.push({ start: Number(p[0]), end: Number(p[1]) });\n",
    "    }\n",
    "    return out;\n",
    "  }\n",
    "  function setDirty(v) { \n",
    "    dirty = v; \n",
    "    var elem = byId('dirty');\n",
    "    elem.textContent = dirty ? 'Unsaved changes' : 'Saved'; \n",
    "    elem.className = dirty ? 'badge badge-warn' : 'badge badge-ok'; \n",
    "  }\n",
    "  function current() { return LOCI[idx]; }\n",
    "  function locusSegments(id) {\n",
    "    var out = [];\n",
    "    for (var i = 0; i < SEGS.length; i++) { \n",
    "      if (SEGS[i].MERGED_LOCUS_ID === id) out.push(SEGS[i]); \n",
    "    }\n",
    "    return out;\n",
    "  }\n",
    "  function locusOriginal(id) {\n",
    "    var out = [];\n",
    "    for (var i = 0; i < ORIGS.length; i++) { \n",
    "      if (ORIGS[i].MERGED_LOCUS_ID === id) out.push(ORIGS[i]); \n",
    "    }\n",
    "    return out;\n",
    "  }\n",
    "  function locusPvals(id) {\n",
    "    var out = [];\n",
    "    for (var i = 0; i < PVALS.length; i++) { \n",
    "      if (PVALS[i].MERGED_LOCUS_ID === id) out.push(PVALS[i]); \n",
    "    }\n",
    "    return out;\n",
    "  }\n",
    "  function clamp(x, min, max) { return Math.max(min, Math.min(max, x)); }\n",
    "  \n",
    "  // Enhanced coordinate functions\n",
    "  function bpToX(bp, ds, de, w) { \n",
    "    var t = (bp - ds) / (de - ds); \n",
    "    return 40 + t * (w - 80); \n",
    "  }\n",
    "  function xToBp(x, ds, de, w) { \n",
    "    var t = (x - 40) / (w - 80); \n",
    "    return Math.round(ds + clamp(t, 0, 1) * (de - ds)); \n",
    "  }\n",
    "  \n",
    "  // Main render function\n",
    "  function render() {\n",
    "    var L = current();\n",
    "    byId('counter').textContent = (idx + 1) + ' / ' + LOCI.length;\n",
    "    byId('merged_id').textContent = L.MERGED_LOCUS_ID;\n",
    "    byId('chr').textContent = L.CHR;\n",
    "    byId('sources').textContent = L.SOURCES || '';\n",
    "    byId('rows').textContent = L.N_ALIGNED_ROWS;\n",
    "    byId('merged_start').textContent = formatNumber(L.MERGED_START);\n",
    "    byId('merged_end').textContent = formatNumber(L.MERGED_END);\n",
    "    byId('refined_start').value = (L.REFINED_START === null || L.REFINED_START === undefined) ? '' : L.REFINED_START;\n",
    "    byId('refined_end').value = (L.REFINED_END === null || L.REFINED_END === undefined) ? '' : L.REFINED_END;\n",
    "    byId('locus_name').value = L.LOCUS_NAME || '';\n",
    "    byId('keep').checked = !!L.KEEP;\n",
    "    var splits = parseSplits(L.SPLITS);\n",
    "    renderSplitsList(splits);\n",
    "    drawEnhancedLocus(L, splits);\n",
    "    setDirty(false);\n",
    "  }\n",
    "  \n",
    "  function formatNumber(n) {\n",
    "    if (n === null || n === undefined) return '';\n",
    "    return Number(n).toLocaleString();\n",
    "  }\n",
    "  \n",
    "  function renderSplitsList(splits) {\n",
    "    var ul = byId('splits');\n",
    "    ul.innerHTML = '';\n",
    "    if (!splits || splits.length === 0) { \n",
    "      ul.innerHTML = '<li class=\"small\">(none)</li>'; \n",
    "      return; \n",
    "    }\n",
    "    for (var i = 0; i < splits.length; i++) {\n",
    "      var s = splits[i];\n",
    "      var li = document.createElement('li');\n",
    "      li.textContent = formatNumber(s.start) + ' - ' + formatNumber(s.end);\n",
    "      var btn = document.createElement('button');\n",
    "      btn.textContent = 'remove';\n",
    "      btn.className = 'btn btn-small';\n",
    "      (function(idxToRemove) {\n",
    "        btn.addEventListener('click', function() {\n",
    "          var L = current();\n",
    "          var arr = parseSplits(L.SPLITS);\n",
    "          arr.splice(idxToRemove, 1);\n",
    "          L.SPLITS = splitsToString(arr);\n",
    "          setDirty(true);\n",
    "          renderSplitsList(arr);\n",
    "          drawEnhancedLocus(L, arr);\n",
    "        });\n",
    "      })(i);\n",
    "      li.appendChild(document.createTextNode(' '));\n",
    "      li.appendChild(btn);\n",
    "      ul.appendChild(li);\n",
    "    }\n",
    "  }\n",
    "  \n",
    "  function updateFromInputs() {\n",
    "    var L = current();\n",
    "    L.REFINED_START = Number(byId('refined_start').value);\n",
    "    L.REFINED_END = Number(byId('refined_end').value);\n",
    "    L.LOCUS_NAME = byId('locus_name').value || '';\n",
    "    L.KEEP = !!byId('keep').checked;\n",
    "    setDirty(true);\n",
    "    drawEnhancedLocus(L, parseSplits(L.SPLITS));\n",
    "  }\n",
    "  \n",
    "  function validateInputs() {\n",
    "    var s = Number(byId('refined_start').value);\n",
    "    var e = Number(byId('refined_end').value);\n",
    "    var ok = isInteger(s) && isInteger(e) && s <= e;\n",
    "    byId('refined_start').classList.toggle('invalid', !ok);\n",
    "    byId('refined_end').classList.toggle('invalid', !ok);\n",
    "    return ok;\n",
    "  }\n",
    "  \n",
    "  function downloadTSV() {\n",
    "    var cols = ['MERGED_LOCUS_ID', 'CHR', 'MERGED_START', 'MERGED_END', 'REFINED_START', 'REFINED_END', 'LOCUS_NAME', 'KEEP', 'SPLITS', 'DISPLAY_START', 'DISPLAY_END', 'FLANK_BP', 'N_BINS', 'N_ALIGNED_ROWS', 'N_ORIGINAL_LOCI', 'SOURCES'];\n",
    "    var lines = [];\n",
    "    lines.push(cols.join('\\t'));\n",
    "    for (var i = 0; i < LOCI.length; i++) {\n",
    "      var L = LOCI[i];\n",
    "      var row = [];\n",
    "      for (var j = 0; j < cols.length; j++) { row.push(escTSV(L[cols[j]])); }\n",
    "      lines.push(row.join('\\t'));\n",
    "    }\n",
    "    var blob = new Blob([lines.join('\\n')], { type: 'text/tab-separated-values;charset=utf-8' });\n",
    "    var url = URL.createObjectURL(blob);\n",
    "    var a = document.createElement('a');\n",
    "    a.href = url;\n",
    "    a.download = DOWNLOAD_NAME;\n",
    "    document.body.appendChild(a);\n",
    "    a.click();\n",
    "    a.remove();\n",
    "    URL.revokeObjectURL(url);\n",
    "    setDirty(false);\n",
    "  }\n",
    "  \n",
    "  function maybeNext(delta) {\n",
    "    if (dirty) {\n",
    "      var ok = confirm('You have unsaved changes. Click OK to continue without saving, or Cancel to stay and save.');\n",
    "      if (!ok) return;\n",
    "    }\n",
    "    idx = clamp(idx + delta, 0, LOCI.length - 1);\n",
    "    selection = null;\n",
    "    render();\n",
    "  }\n",
    "  \n",
    "  // Enhanced drawing function\n",
    "  function drawEnhancedLocus(L, splits) {\n",
    "    var svg = byId('viz');\n",
    "    var w = svg.viewBox.baseVal.width || 900;\n",
    "    var h = svg.viewBox.baseVal.height || 220;\n",
    "    var ds = Number(L.DISPLAY_START);\n",
    "    var de = Number(L.DISPLAY_END);\n",
    "    \n",
    "    // Clear SVG\n",
    "    while (svg.firstChild) svg.removeChild(svg.firstChild);\n",
    "    var ns = 'http://www.w3.org/2000/svg';\n",
    "    \n",
    "    // Helper to create elements\n",
    "    function createElement(type, attrs) {\n",
    "      var el = document.createElementNS(ns, type);\n",
    "      for (var key in attrs) {\n",
    "        if (attrs.hasOwnProperty(key)) {\n",
    "          el.setAttribute(key, attrs[key]);\n",
    "        }\n",
    "      }\n",
    "      return el;\n",
    "    }\n",
    "    \n",
    "    // Title\n",
    "    var title = createElement('text', {\n",
    "      x: '40',\n",
    " y: '20',\n",
    "      'class': 'viz-title',\n",
    "      'font-size': '14px',\n",
    "      'font-weight': 'bold'\n",
    "    });\n",
    "    title.textContent = 'chr' + L.CHR + ': ' + formatNumber(ds) + ' - ' + formatNumber(de);\n",
    "    svg.appendChild(title);\n",
    "    \n",
    "    // Ruler\n",
    "    var rulerY = 40;\n",
    "    var rulerLine = createElement('line', {\n",
    "      x1: '40',\n",
    "      y1: rulerY,\n",
    "      x2: (w - 40).toString(),\n",
    "      y2: rulerY,\n",
    "      'class': 'ruler-line',\n",
    "      'stroke': '#333',\n",
    "      'stroke-width': '1'\n",
    "    });\n",
    "    svg.appendChild(rulerLine);\n",
    "    \n",
    "    // Add tick marks\n",
    "    var span = de - ds;\n",
    "    var tickInterval = 10000;\n",
    "    if (span > 100000) tickInterval = 50000;\n",
    "    if (span > 500000) tickInterval = 100000;\n",
    "    \n",
    "    var startTick = Math.ceil(ds / tickInterval) * tickInterval;\n",
    "    for (var pos = startTick; pos <= de; pos += tickInterval) {\n",
    "      var x = bpToX(pos, ds, de, w);\n",
    "      if (x < 40 || x > w - 40) continue;\n",
    "      \n",
    "      var tick = createElement('line', {\n",
    "        x1: x.toString(),\n",
    "        y1: (rulerY - 5).toString(),\n",
    "        x2: x.toString(),\n",
    "        y2: (rulerY + 5).toString(),\n",
    "        'stroke': '#666',\n",
    "        'stroke-width': '1'\n",
    "      });\n",
    "      svg.appendChild(tick);\n",
    "      \n",
    "      var label = createElement('text', {\n",
    "        x: x.toString(),\n",
    "        y: (rulerY - 8).toString(),\n",
    "        'class': 'ruler-label',\n",
    "        'text-anchor': 'middle',\n",
    "        'font-size': '10px',\n",
    "        'fill': '#666'\n",
    "      });\n",
    "      label.textContent = (pos / 1000).toFixed(0) + 'kb';\n",
    "      svg.appendChild(label);\n",
    "    }\n",
    "    \n",
    "    var trackY = 60;\n",
    "    var trackHeight = 25;\n",
    "    var trackSpacing = 30;\n",
    "    \n",
    "    // P-value track\n",
    "    var pvals = locusPvals(L.MERGED_LOCUS_ID);\n",
    "    if (pvals.length > 0) {\n",
    "      drawPvalueTrack(svg, pvals, ds, de, w, trackY, trackHeight);\n",
    "      trackY += trackSpacing;\n",
    "    }\n",
    "    \n",
    "    // Original loci track\n",
    "    var origs = locusOriginal(L.MERGED_LOCUS_ID);\n",
    "    if (origs.length > 0) {\n",
    "      drawOriginalTrack(svg, origs, ds, de, w, trackY, trackHeight);\n",
    "      trackY += trackSpacing;\n",
    "    }\n",
    "    \n",
    "    // Merged/Refined track\n",
    "    drawMergedTrack(svg, L, ds, de, w, trackY, trackHeight);\n",
    "    trackY += trackSpacing;\n",
    "    \n",
    "    // Segments track\n",
    "    var segs = locusSegments(L.MERGED_LOCUS_ID);\n",
    "    if (segs.length > 0) {\n",
    "      drawSegmentsTrack(svg, segs, ds, de, w, trackY, trackHeight);\n",
    "    }\n",
    "    \n",
    "    // Splits overlay\n",
    "    if (splits && splits.length > 0) {\n",
    "      drawSplitsOverlay(svg, splits, ds, de, w, 80, trackHeight * 3);\n",
    "    }\n",
    "    \n",
    "    // Selection overlay\n",
    "    if (selection) {\n",
    "      drawSelectionOverlay(svg, selection, ds, de, w, 80, trackHeight * 3);\n",
    "    }\n",
    "  }\n",
    "  \n",
    "  function drawPvalueTrack(svg, pvals, ds, de, w, y, height) {\n",
    "    var ns = 'http://www.w3.org/2000/svg';\n",
    "    \n",
    "    // Track background\n",
    "    var bg = document.createElementNS(ns, 'rect');\n",
    "    bg.setAttribute('x', '40');\n",
    "    bg.setAttribute('y', y.toString());\n",
    "    bg.setAttribute('width', (w - 80).toString());\n",
    "    bg.setAttribute('height', height.toString());\n",
    "    bg.setAttribute('fill', '#f8f9fa');\n",
    "    bg.setAttribute('stroke', '#e9ecef');\n",
    "    bg.setAttribute('stroke-width', '1');\n",
    "    svg.appendChild(bg);\n",
    "    \n",
    "    // Find max p-value\n",
    "    var maxP = 0;\n",
    "    for (var i = 0; i < pvals.length; i++) {\n",
    "      if (pvals[i].MAX_NEG_LOG10_P > maxP) maxP = pvals[i].MAX_NEG_LOG10_P;\n",
    "    }\n",
    "    maxP = Math.max(maxP, 1);\n",
    "    \n",
    "    // Draw bars\n",
    "    for (var i = 0; i < pvals.length; i++) {\n",
    "      var p = pvals[i];\n",
    "      var x1 = bpToX(p.WIN_START, ds, de, w);\n",
    "      var x2 = bpToX(p.WIN_END, ds, de, w);\n",
    "      var barWidth = Math.max(1, x2 - x1);\n",
    "      \n",
    "      var barHeight = (p.MAX_NEG_LOG10_P / maxP) * (height - 4);\n",
    "      \n",
    "      var bar = document.createElementNS(ns, 'rect');\n",
    "      bar.setAttribute('x', x1.toString());\n",
    "      bar.setAttribute('y', (y + height - barHeight).toString());\n",
    "      bar.setAttribute('width', barWidth.toString());\n",
    "      bar.setAttribute('height', barHeight.toString());\n",
    "      bar.setAttribute('fill', '#e53e3e');\n",
    "      bar.setAttribute('opacity', '0.7');\n",
    "      svg.appendChild(bar);\n",
    "    }\n",
    "    \n",
    "    // Label\n",
    "    var label = document.createElementNS(ns, 'text');\n",
    "    label.setAttribute('x', '20');\n",
    "    label.setAttribute('y', (y + height/2 + 4).toString());\n",
    "    label.setAttribute('font-size', '11px');\n",
    "    label.setAttribute('fill', '#6c757d');\n",
    "    label.setAttribute('text-anchor', 'end');\n",
    "    label.textContent = 'P-values';\n",
    "    svg.appendChild(label);\n",
    "  }\n",
    "  \n",
    "  function drawOriginalTrack(svg, origs, ds, de, w, y, height) {\n",
    "    var ns = 'http://www.w3.org/2000/svg';\n",
    "    \n",
    "    // Track background\n",
    "    var bg = document.createElementNS(ns, 'rect');\n",
    "    bg.setAttribute('x', '40');\n",
    "    bg.setAttribute('y', y.toString());\n",
    "    bg.setAttribute('width', (w - 80).toString());\n",
    "    bg.setAttribute('height', height.toString());\n",
    "    bg.setAttribute('fill', '#f8f9fa');\n",
    "    bg.setAttribute('stroke', '#e9ecef');\n",
    "    bg.setAttribute('stroke-width', '1');\n",
    "    svg.appendChild(bg);\n",
    "    \n",
    "    // Group by source\n",
    "    var sources = {};\n",
    "    for (var i = 0; i < origs.length; i++) {\n",
    "      var o = origs[i];\n",
    "      var source = o.SOURCE || 'unknown';\n",
    "      if (!sources[source]) sources[source] = [];\n",
    "      sources[source].push(o);\n",
    "    }\n",
    "    \n",
    "    var sourceKeys = Object.keys(sources);\n",
    "    var rowHeight = height / Math.max(1, sourceKeys.length);\n",
    "    \n",
    "    sourceKeys.forEach(function(source, idx) {\n",
    "      var sourceY = y + idx * rowHeight + 2;\n",
    "      var sourceColor = COLORS[source] || '#718096';\n",
    "      \n",
    "      sources[source].forEach(function(o) {\n",
    "        var x1 = bpToX(o.START, ds, de, w);\n",
    "        var x2 = bpToX(o.END, ds, de, w);\n",
    "        var width = Math.max(1, x2 - x1);\n",
    "        \n",
    "        var rect = document.createElementNS(ns, 'rect');\n",
    "        rect.setAttribute('x', x1.toString());\n",
    "        rect.setAttribute('y', sourceY.toString());\n",
    "        rect.setAttribute('width', width.toString());\n",
    "        rect.setAttribute('height', (rowHeight - 4).toString());\n",
    "        rect.setAttribute('fill', sourceColor);\n",
    "        rect.setAttribute('opacity', '0.8');\n",
    "        svg.appendChild(rect);\n",
    "      });\n",
    "    });\n",
    "    \n",
    "    // Label\n",
    "    var label = document.createElementNS(ns, 'text');\n",
    "    label.setAttribute('x', '20');\n",
    "    label.setAttribute('y', (y + height/2 + 4).toString());\n",
    "    label.setAttribute('font-size', '11px');\n",
    "    label.setAttribute('fill', '#6c757d');\n",
    "    label.setAttribute('text-anchor', 'end');\n",
    "    label.textContent = 'Original';\n",
    "    svg.appendChild(label);\n",
    "  }\n",
    "  \n",
    "  function drawMergedTrack(svg, L, ds, de, w, y, height) {\n",
    "    var ns = 'http://www.w3.org/2000/svg';\n",
    "    \n",
    "    // Track background\n",
    "    var bg = document.createElementNS(ns, 'rect');\n",
    "    bg.setAttribute('x', '40');\n",
    "    bg.setAttribute('y', y.toString());\n",
    "    bg.setAttribute('width', (w - 80).toString());\n",
    "    bg.setAttribute('height', height.toString());\n",
    "    bg.setAttribute('fill', '#f8f9fa');\n",
    "    bg.setAttribute('stroke', '#e9ecef');\n",
    "    bg.setAttribute('stroke-width', '1');\n",
    "    svg.appendChild(bg);\n",
    "    \n",
    "    // Merged locus\n",
    "    var mergedX1 = bpToX(L.MERGED_START, ds, de, w);\n",
    "    var mergedX2 = bpToX(L.MERGED_END, ds, de, w);\n",
    "    var mergedRect = document.createElementNS(ns, 'rect');\n",
    "    mergedRect.setAttribute('x', mergedX1.toString());\n",
    "    mergedRect.setAttribute('y', (y + 5).toString());\n",
    "    mergedRect.setAttribute('width', Math.max(1, mergedX2 - mergedX1).toString());\n",
    "    mergedRect.setAttribute('height', (height - 10).toString());\n",
    "    mergedRect.setAttribute('fill', '#4299e1');\n",
    "    mergedRect.setAttribute('opacity', '0.3');\n",
    "    mergedRect.setAttribute('stroke', '#3182ce');\n",
    "    mergedRect.setAttribute('stroke-width', '1');\n",
    "    svg.appendChild(mergedRect);\n",
    "    \n",
    "    // Refined locus\n",
    "    if (L.REFINED_START && L.REFINED_END) {\n",
    "      var refinedX1 = bpToX(L.REFINED_START, ds, de, w);\n",
    "      var refinedX2 = bpToX(L.REFINED_END, ds, de, w);\n",
    "      var refinedRect = document.createElementNS(ns, 'rect');\n",
    "      refinedRect.setAttribute('x', refinedX1.toString());\n",
    "      refinedRect.setAttribute('y', (y + 2).toString());\n",
    "      refinedRect.setAttribute('width', Math.max(1, refinedX2 - refinedX1).toString());\n",
    "      refinedRect.setAttribute('height', (height - 4).toString());\n",
    "      refinedRect.setAttribute('fill', '#38a169');\n",
    "      refinedRect.setAttribute('opacity', '0.6');\n",
    "      refinedRect.setAttribute('stroke', '#2f855a');\n",
    "      refinedRect.setAttribute('stroke-width', '2');\n",
    "      svg.appendChild(refinedRect);\n",
    "    }\n",
    "    \n",
    "    // Label\n",
    "    var label = document.createElementNS(ns, 'text');\n",
    "    label.setAttribute('x', '20');\n",
    "    label.setAttribute('y', (y + height/2 + 4).toString());\n",
    "    label.setAttribute('font-size', '11px');\n",
    "    label.setAttribute('fill', '#6c757d');\n",
    "    label.setAttribute('text-anchor', 'end');\n",
    "    label.textContent = 'Merged/Refined';\n",
    "    svg.appendChild(label);\n",
    "  }\n",
    "  \n",
    "  function drawSegmentsTrack(svg, segs, ds, de, w, y, height) {\n",
    "    var ns = 'http://www.w3.org/2000/svg';\n",
    "    \n",
    "    // Track background\n",
    "    var bg = document.createElementNS(ns, 'rect');\n",
    "    bg.setAttribute('x', '40');\n",
    "    bg.setAttribute('y', y.toString());\n",
    "    bg.setAttribute('width', (w - 80).toString());\n",
    "    bg.setAttribute('height', height.toString());\n",
    "    bg.setAttribute('fill', '#f8f9fa');\n",
    "    bg.setAttribute('stroke', '#e9ecef');\n",
    "    bg.setAttribute('stroke-width', '1');\n",
    "    svg.appendChild(bg);\n",
    "    \n",
    "    // Draw segments\n",
    "    var maxSegs = Math.min(segs.length, 30);\n",
    "    var segHeight = height - 10;\n",
    "    \n",
    "    for (var i = 0; i < maxSegs; i++) {\n",
    "      var s = segs[i];\n",
    "      if (!s.START || !s.END) continue;\n",
    "      \n",
    "      var x1 = bpToX(s.START, ds, de, w);\n",
    "      var x2 = bpToX(s.END, ds, de, w);\n",
    "      var segWidth = Math.max(1, x2 - x1);\n",
    "      \n",
    "      var rect = document.createElementNS(ns, 'rect');\n",
    "      rect.setAttribute('x', x1.toString());\n",
    "      rect.setAttribute('y', (y + 5).toString());\n",
    "      rect.setAttribute('width', segWidth.toString());\n",
    "      rect.setAttribute('height', segHeight.toString());\n",
    "      rect.setAttribute('fill', '#444');\n",
    "      rect.setAttribute('opacity', '0.55');\n",
    "      svg.appendChild(rect);\n",
    "    }\n",
    "    \n",
    "    // Label\n",
    "    var label = document.createElementNS(ns, 'text');\n",
    "    label.setAttribute('x', '20');\n",
    "    label.setAttribute('y', (y + height/2 + 4).toString());\n",
    "    label.setAttribute('font-size', '11px');\n",
    "    label.setAttribute('fill', '#6c757d');\n",
    "    label.setAttribute('text-anchor', 'end');\n",
    "    label.textContent = 'Segments';\n",
    "    svg.appendChild(label);\n",
    "  }\n",
    "  \n",
    "  function drawSplitsOverlay(svg, splits, ds, de, w, y, height) {\n",
    "    var ns = 'http://www.w3.org/2000/svg';\n",
    "    \n",
    "    for (var i = 0; i < splits.length; i++) {\n",
    "      var sp = splits[i];\n",
    "      var x1 = bpToX(sp.start, ds, de, w);\n",
    "      var x2 = bpToX(sp.end, ds, de, w);\n",
    "      var width = Math.max(1, x2 - x1);\n",
    "      \n",
    "      var rect = document.createElementNS(ns, 'rect');\n",
    "      rect.setAttribute('x', x1.toString());\n",
    "      rect.setAttribute('y', y.toString());\n",
    "      rect.setAttribute('width', width.toString());\n",
    "      rect.setAttribute('height', height.toString());\n",
    "      rect.setAttribute('fill', '#fef3c7');\n",
    "      rect.setAttribute('opacity', '0.5');\n",
    "      rect.setAttribute('stroke', '#f59e0b');\n",
    "      rect.setAttribute('stroke-width', '1');\n",
    "      svg.appendChild(rect);\n",
    "    }\n",
    "  }\n",
    "  \n",
    "  function drawSelectionOverlay(svg, selection, ds, de, w, y, height) {\n",
    "    var ns = 'http://www.w3.org/2000/svg';\n",
    "    var a = Math.min(selection.start, selection.end);\n",
    "    var b = Math.max(selection.start, selection.end);\n",
    "    \n",
    "    var x1 = bpToX(a, ds, de, w);\n",
    "    var x2 = bpToX(b, ds, de, w);\n",
    "    var width = Math.max(1, x2 - x1);\n",
    "    \n",
    "    var rect = document.createElementNS(ns, 'rect');\n",
    "    rect.setAttribute('x', x1.toString());\n",
    "    rect.setAttribute('y', y.toString());\n",
    "    rect.setAttribute('width', width.toString());\n",
    "    rect.setAttribute('height', height.toString());\n",
    "    rect.setAttribute('fill', '#c7d2fe');\n",
    "    rect.setAttribute('opacity', '0.5');\n",
    "    svg.appendChild(rect);\n",
    "  }\n",
    "  \n",
    "  // Setup mouse interactions for selection\n",
    "  function setupVizMouse() {\n",
    "    var svg = byId('viz');\n",
    "    var dragging = false;\n",
    "    var x0 = null;\n",
    "    \n",
    "    svg.addEventListener('mousedown', function(ev) {\n",
    "      var L = current();\n",
    "      var w = svg.viewBox.baseVal.width || 900;\n",
    "      var ds = Number(L.DISPLAY_START);\n",
    "      var de = Number(L.DISPLAY_END);\n",
    "      dragging = true;\n",
    "      x0 = ev.offsetX;\n",
    "      var bp0 = xToBp(x0, ds, de, w);\n",
    "      selection = { start: bp0, end: bp0 };\n",
    "      drawEnhancedLocus(L, parseSplits(L.SPLITS));\n",
    "    });\n",
    "    \n",
    "    svg.addEventListener('mousemove', function(ev) {\n",
    "      if (!dragging || !selection) return;\n",
    "      var L = current();\n",
    "      var w = svg.viewBox.baseVal.width || 900;\n",
    "      var ds = Number(L.DISPLAY_START);\n",
    "      var de = Number(L.DISPLAY_END);\n",
    "      var bp = xToBp(ev.offsetX, ds, de, w);\n",
    "      selection.end = bp;\n",
    "      drawEnhancedLocus(L, parseSplits(L.SPLITS));\n",
    "    });\n",
    "    \n",
    "    window.addEventListener('mouseup', function() {\n",
    "      dragging = false;\n",
    "    });\n",
    "  }\n",
    "  \n",
    "  function addSplitFromSelection() {\n",
    "    if (!selection) { \n",
    "      alert('Drag on the locus visualization to select a region first.'); \n",
    "      return; \n",
    "    }\n",
    "    var a = Math.min(selection.start, selection.end);\n",
    "    var b = Math.max(selection.start, selection.end);\n",
    "    var L = current();\n",
    "    var arr = parseSplits(L.SPLITS);\n",
    "    arr.push({ start: a, end: b });\n",
    "    arr.sort(function(x, y) { return x.start - y.start; });\n",
    "    var merged = [];\n",
    "    for (var i = 0; i < arr.length; i++) {\n",
    "      var s = arr[i];\n",
    "      if (merged.length === 0) { \n",
    "        merged.push({ start: s.start, end: s.end }); \n",
    "        continue; \n",
    "      }\n",
    "      var last = merged[merged.length - 1];\n",
    "      if (s.start <= last.end) { \n",
    "        last.end = Math.max(last.end, s.end); \n",
    "      } else { \n",
    "        merged.push({ start: s.start, end: s.end }); \n",
    "      }\n",
    "    }\n",
    "    L.SPLITS = splitsToString(merged);\n",
    "    selection = null;\n",
    "    setDirty(true);\n",
    "    renderSplitsList(merged);\n",
    "    drawEnhancedLocus(L, merged);\n",
    "  }\n",
    "  \n",
    "  // Initialize event listeners\n",
    "  function initEvents() {\n",
    "    byId('prev').addEventListener('click', function() { maybeNext(-1); });\n",
    "    byId('next').addEventListener('click', function() { maybeNext(1); });\n",
    "    byId('save').addEventListener('click', downloadTSV);\n",
    "    byId('add_split').addEventListener('click', addSplitFromSelection);\n",
    "    \n",
    "    var inputs = ['refined_start', 'refined_end', 'locus_name', 'keep'];\n",
    "    inputs.forEach(function(id) {\n",
    "      var el = byId(id);\n",
    "      el.addEventListener('input', function() { \n",
    "        updateFromInputs(); \n",
    "        validateInputs(); \n",
    "      });\n",
    "      el.addEventListener('change', function() { \n",
    "        updateFromInputs(); \n",
    "        validateInputs(); \n",
    "      });\n",
    "    });\n",
    "    \n",
    "    setupVizMouse();\n",
    "    render();\n",
    "  }\n",
    "  \n",
    "  // Start everything\n",
    "  try {\n",
    "    initEvents();\n",
    "  } catch (e) {\n",
    "    showError(e && e.stack ? e.stack : e);\n",
    "  }\n",
    "})();\n"
  )
  
  # Enhanced CSS
  css <- htmltools::HTML("
    body { 
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif; 
      margin: 20px; 
      line-height: 1.5;
    }
    
    h1 { 
      color: #2d3748; 
      margin-bottom: 10px; 
      font-size: 24px;
    }
    
    .row { 
      display: flex; 
      gap: 20px; 
      align-items: flex-start; 
      margin-top: 20px;
    }
    
    .panel { 
      border: 1px solid #e2e8f0; 
      padding: 16px; 
      border-radius: 8px; 
      background: white;
      box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
    }
    
    .left { 
      flex: 1 1 auto; 
      min-width: 600px; 
    }
    
    .right { 
      width: 320px; 
      flex-shrink: 0;
    }
    
    .small { 
      color: #718096; 
      font-size: 12px; 
    }
    
    .controls { 
      display: flex; 
      gap: 10px; 
      align-items: center; 
      margin-bottom: 20px;
      padding: 12px;
      background: #f7fafc;
      border-radius: 6px;
    }
    
    .btn { 
      padding: 8px 16px; 
      cursor: pointer; 
      background: #4299e1;
      color: white;
      border: none;
      border-radius: 4px;
      font-size: 14px;
      font-weight: 500;
      transition: background-color 0.2s;
    }
    
    .btn:hover { 
      background: #3182ce; 
    }
    
    .btn-small { 
      padding: 2px 8px; 
      font-size: 12px; 
      cursor: pointer; 
      background: #e2e8f0;
      color: #4a5568;
      border: 1px solid #cbd5e0;
      border-radius: 3px;
    }
    
    .btn-small:hover {
      background: #cbd5e0;
    }
    
    .badge { 
      padding: 4px 10px; 
      border-radius: 12px; 
      font-size: 12px; 
      font-weight: 500;
      margin-left: auto;
    }
    
    .badge-ok { 
      background: #c6f6d5; 
      color: #22543d; 
    }
    
    .badge-warn { 
      background: #fed7d7; 
      color: #742a2a; 
    }
    
    label { 
      display: block; 
      margin-top: 12px; 
      font-weight: 600; 
      color: #4a5568;
      font-size: 14px;
    }
    
    input[type='text'], 
    input[type='number'] { 
      width: 100%; 
      box-sizing: border-box; 
      padding: 8px; 
      border: 1px solid #cbd5e0;
      border-radius: 4px;
      font-size: 14px;
      margin-top: 4px;
    }
    
    input.invalid { 
      outline: 2px solid #fc8181; 
      border-color: #fc8181;
    }
    
    input[type='checkbox'] {
      margin-top: 12px;
      width: 18px;
      height: 18px;
    }
    
    ul { 
      padding-left: 20px; 
      margin: 8px 0;
    }
    
    li {
      margin: 4px 0;
      font-size: 13px;
    }
    
    #viz { 
      width: 100%; 
      height: 240px; 
      border: 1px solid #e2e8f0; 
      border-radius: 6px; 
      background: white;
      margin: 10px 0;
    }
    
    .viz-title {
      font-weight: bold;
      fill: #2d3748;
    }
    
    .ruler-line {
      stroke: #4a5568;
      stroke-width: 1;
    }
    
    .ruler-label {
      font-size: 10px;
      fill: #718096;
    }
    
    .track-label {
      font-size: 11px;
      fill: #718096;
    }
    
    #js_error {
      display: none; 
      white-space: pre-wrap; 
      background: #fed7d7; 
      border: 1px solid #fc8181; 
      padding: 12px; 
      border-radius: 6px; 
      margin: 10px 0;
      color: #742a2a;
      font-size: 13px;
    }
    
    .locus-info-grid {
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
      gap: 12px;
      margin: 15px 0;
      padding: 12px;
      background: #f7fafc;
      border-radius: 6px;
    }
    
    .info-item {
      font-size: 13px;
    }
    
    .info-label {
      font-weight: 600;
      color: #4a5568;
    }
    
    .info-value {
      color: #2d3748;
      font-family: 'Monaco', 'Menlo', 'Ubuntu Mono', monospace;
    }
  ")
  
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
          htmltools::tags$span(class = "info-label", "MERGED_LOCUS_ID: "),
          htmltools::tags$span(id = "merged_id", class = "info-value")
        ),
        htmltools::tags$div(class = "info-item",
          htmltools::tags$span(class = "info-label", "Chromosome: "),
          htmltools::tags$span(id = "chr", class = "info-value")
        ),
        htmltools::tags$div(class = "info-item",
          htmltools::tags$span(class = "info-label", "Rows: "),
          htmltools::tags$span(id = "rows", class = "info-value")
        ),
        htmltools::tags$div(class = "info-item",
          htmltools::tags$span(class = "info-label", "Sources: "),
          htmltools::tags$span(id = "sources", class = "info-value")
        ),
        htmltools::tags$div(class = "info-item",
          htmltools::tags$span(class = "info-label", "Merged: "),
          htmltools::tags$span(id = "merged_start", class = "info-value"),
          " - ",
          htmltools::tags$span(id = "merged_end", class = "info-value")
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