#' Export aligned loci to an interactive locus-at-a-time HTML editor (wizard)
#'
#' Generates a single self-contained HTML file (offline) that shows one
#' MERGED_LOCUS_ID at a time. The user can:
#' - Navigate Prev/Next through loci
#' - Edit REFINED_START/REFINED_END, LOCUS_NAME, KEEP
#' - Drag-select a region on the locus visualization and add it as a "split"
#' - Click "Save (Download TSV)" to download the current edits as TSV
#'
#' Note: A static HTML file cannot write back to disk automatically. "Save"
#' means "download a TSV".
#'
#' @param aligned_loci data.frame returned by `align_with_reported()` (class 'AlignedLoci')
#' @param output_html character. Output HTML filepath (e.g., 'outputs/loci_wizard.html')
#' @param flank_bp numeric. Base-pairs to extend display window for visualization.
#' @param n_bins integer. Stored in exported TSV for reproducibility.
#'
#' @return Named list with path: `list(html = <path>)`.
#' @export
export_aligned_to_html_locus_wizard <- function(aligned_loci,
                                                output_html,
                                                flank_bp = 20000,
                                                n_bins = 20,
                                                gwas_list = NULL,
                                                window_bp = 20000) {
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
  # Keep first occurrence values for merged coordinates and defaults.
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
      # locus_id(source:start:end)
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
    # Limit to sources seen in original loci if available; otherwise include all GWAS names.
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
        # Pre-filter to locus span for speed.
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

  # We'll embed data as proper JSON to avoid JS parse issues.

  # Build segments for reported loci (optional extra context track).
  seg_keep <- c("MERGED_LOCUS_ID", "ALIGNED_LOCUS_ID", "START", "END", "ALIGNMENT_STATUS")
  if ("REPORTED_SOURCE" %in% colnames(df)) {
    seg_keep <- c(seg_keep, "REPORTED_SOURCE")
  }
  seg_df <- df[, seg_keep, drop = FALSE]
  # Coerce START/END numeric
  seg_df$START <- suppressWarnings(as.numeric(seg_df$START))
  seg_df$END <- suppressWarnings(as.numeric(seg_df$END))

  loci_js <- jsonlite::toJSON(loci_df, dataframe = "rows", na = "null", auto_unbox = TRUE)
  seg_js <- jsonlite::toJSON(seg_df, dataframe = "rows", na = "null", auto_unbox = TRUE)
  orig_js <- jsonlite::toJSON(orig_df, dataframe = "rows", na = "null", auto_unbox = TRUE)
  colors_js <- jsonlite::toJSON(source_colors, auto_unbox = TRUE)
  pvals_js <- jsonlite::toJSON(pvals_df, dataframe = "rows", na = "null", auto_unbox = TRUE)

  download_name <- paste0(tools::file_path_sans_ext(basename(output_html)), ".tsv")
  download_name_js <- jsonlite::toJSON(download_name, auto_unbox = TRUE)

  js <- paste0(
    "(function(){\n",
    "  function byId(id){ return document.getElementById(id); }\n",
    "  function showError(msg){\n",
    "    try {\n",
    "      var box = byId('js_error');\n",
    "      if(box){ box.style.display='block'; box.textContent = String(msg); }\n",
    "    } catch(e) {}\n",
    "  }\n",
    "  try {\n",
    "    var DOWNLOAD_NAME = ", download_name_js, ";\n",
    "    var N_BINS = ", as.integer(n_bins), ";\n",
    "    var FLANK_BP = ", as.numeric(flank_bp), ";\n",
    "    var WINDOW_BP = ", as.integer(window_bp), ";\n",
    "    var LOCI = ", loci_js, ";\n",
    "    var SEGS = ", seg_js, ";\n",
    "    var ORIGS = ", orig_js, ";\n",
    "    var COLORS = ", colors_js, ";\n",
    "    var PVALS = ", pvals_js, ";\n",
    "    var idx = 0;\n",
    "    var dirty = false;\n",
    "    var selection = null;\n",
    "    function escTSV(v){\n",
    "      if(v===null||v===undefined) return '';\n",
    "      var s = String(v);\n",
    "      // Keep TSV one-line-per-row by removing tabs/newlines.\n",
    "      s = s.replace(/\\t/g,' ');\n",
    "      s = s.replace(/\\r/g,' ');\n",
    "      s = s.replace(/\\n/g,' ');\n",
    "      return s;\n",
    "    }\n",
    "    function isInteger(n){ return typeof n==='number' && isFinite(n) && Math.floor(n)===n; }\n",
    "    function splitsToString(arr){\n",
    "      if(!arr || arr.length===0) return '';\n",
    "      var out=[];\n",
    "      for(var i=0;i<arr.length;i++){ out.push(arr[i].start+'-'+arr[i].end); }\n",
    "      return out.join(';');\n",
    "    }\n",
    "    function parseSplits(s){\n",
    "      var out=[];\n",
    "      if(!s) return out;\n",
    "      var parts = String(s).split(';');\n",
    "      for(var i=0;i<parts.length;i++){\n",
    "        var t = parts[i];\n",
    "        if(!t) continue;\n",
    "        t = t.replace(/^\\s+|\\s+$/g, '');\n",
    "        if(t.length===0) continue;\n",
    "        var p = t.split('-');\n",
    "        if(p.length<2) continue;\n",
    "        out.push({start:Number(p[0]), end:Number(p[1])});\n",
    "      }\n",
    "      return out;\n",
    "    }\n",
    "    function setDirty(v){ dirty=v; byId('dirty').textContent = dirty ? 'Unsaved changes' : 'Saved'; byId('dirty').className = dirty ? 'badge badge-warn' : 'badge badge-ok'; }\n",
    "    function current(){ return LOCI[idx]; }\n",
    "    function locusSegments(id){\n",
    "      var out=[];\n",
    "      for(var i=0;i<SEGS.length;i++){ if(SEGS[i].MERGED_LOCUS_ID===id) out.push(SEGS[i]); }\n",
    "      return out;\n",
    "    }\n",
    "    function locusOriginal(id){\n",
    "      var out=[];\n",
    "      for(var i=0;i<ORIGS.length;i++){ if(ORIGS[i].MERGED_LOCUS_ID===id) out.push(ORIGS[i]); }\n",
    "      return out;\n",
    "    }\n",
    "    function locusPvals(id){\n",
    "      var out=[];\n",
    "      for(var i=0;i<PVALS.length;i++){ if(PVALS[i].MERGED_LOCUS_ID===id) out.push(PVALS[i]); }\n",
    "      return out;\n",
    "    }\n",
    "function clamp(x,min,max){ return Math.max(min, Math.min(max, x)); }\n",
    "function bpToX(bp, ds, de, w){ const t=(bp-ds)/(de-ds); return 40 + t*(w-80); }\n",
    "function xToBp(x, ds, de, w){ const t=(x-40)/(w-80); return Math.round(ds + clamp(t,0,1)*(de-ds)); }\n",
    "function render(){\n",
    "  const L = current();\n",
    "  byId('counter').textContent = (idx+1)+' / '+LOCI.length;\n",
    "  byId('merged_id').textContent = L.MERGED_LOCUS_ID;\n",
    "  byId('chr').textContent = L.CHR;\n",
    "  byId('sources').textContent = L.SOURCES || '';\n",
    "  byId('rows').textContent = L.N_ALIGNED_ROWS;\n",
    "  byId('merged_start').textContent = L.MERGED_START;\n",
    "  byId('merged_end').textContent = L.MERGED_END;\n",
    "  byId('refined_start').value = (L.REFINED_START===null || L.REFINED_START===undefined) ? '' : L.REFINED_START;\n",
    "  byId('refined_end').value = (L.REFINED_END===null || L.REFINED_END===undefined) ? '' : L.REFINED_END;\n",
    "  byId('locus_name').value = L.LOCUS_NAME || '';\n",
    "  byId('keep').checked = !!L.KEEP;\n",
    "  const splits = parseSplits(L.SPLITS);\n",
    "  renderSplitsList(splits);\n",
    "  drawLocus(L, splits);\n",
    "  setDirty(false);\n",
    "}\n",
    "function renderSplitsList(splits){\n",
    "  const ul = byId('splits');\n",
    "  ul.innerHTML='';\n",
    "  if(!splits || splits.length===0){ ul.innerHTML='<li class=small>(none)</li>'; return; }\n",
    "  for(var i=0;i<splits.length;i++){\n",
    "    var s = splits[i];\n",
    "    const li=document.createElement('li');\n",
    "    li.textContent = s.start+' - '+s.end;\n",
    "    const btn=document.createElement('button');\n",
    "    btn.textContent='remove';\n",
    "    btn.className='btn btn-small';\n",
    "    (function(idxToRemove){\n",
    "      btn.addEventListener('click', function(){\n",
    "      const L=current();\n",
    "      const arr=parseSplits(L.SPLITS);\n",
    "      arr.splice(idxToRemove,1);\n",
    "      L.SPLITS=splitsToString(arr);\n",
    "      setDirty(true);\n",
    "      renderSplitsList(arr);\n",
    "      drawLocus(L, arr);\n",
    "      });\n",
    "    })(i);\n",
    "    li.appendChild(document.createTextNode(' '));\n",
    "    li.appendChild(btn);\n",
    "    ul.appendChild(li);\n",
    "  }\n",
    "}\n",
    "function updateFromInputs(){\n",
    "  const L=current();\n",
    "  L.REFINED_START = Number(byId('refined_start').value);\n",
    "  L.REFINED_END = Number(byId('refined_end').value);\n",
    "  L.LOCUS_NAME = byId('locus_name').value || '';\n",
    "  L.KEEP = !!byId('keep').checked;\n",
    "  setDirty(true);\n",
    "  drawLocus(L, parseSplits(L.SPLITS));\n",
    "}\n",
    "function validateInputs(){\n",
    "  var s=Number(byId('refined_start').value);\n",
    "  var e=Number(byId('refined_end').value);\n",
    "  var ok = isInteger(s) && isInteger(e) && s<=e;\n",
    "  byId('refined_start').classList.toggle('invalid', !ok);\n",
    "  byId('refined_end').classList.toggle('invalid', !ok);\n",
    "  return ok;\n",
    "}\n",
    "function downloadTSV(){\n",
    "  // locus-level TSV\n",
    "  const cols = ['MERGED_LOCUS_ID','CHR','MERGED_START','MERGED_END','REFINED_START','REFINED_END','LOCUS_NAME','KEEP','SPLITS','DISPLAY_START','DISPLAY_END','FLANK_BP','N_BINS','N_ALIGNED_ROWS','N_ORIGINAL_LOCI','SOURCES'];\n",
    "  const lines=[];\n",
    "  lines.push(cols.join('\t'));\n",
    "  for(var i=0;i<LOCI.length;i++){\n",
    "    var L = LOCI[i];\n",
    "    var row=[];\n",
    "    for(var j=0;j<cols.length;j++){ row.push(escTSV(L[cols[j]])); }\n",
    "    lines.push(row.join('\\t'));\n",
    "  }\n",
    "  const blob = new Blob([lines.join('\\n')], {type:'text/tab-separated-values;charset=utf-8'});\n",
    "  const url = URL.createObjectURL(blob);\n",
    "  const a = document.createElement('a');\n",
    "  a.href = url;\n",
    "  a.download = DOWNLOAD_NAME;\n",
    "  document.body.appendChild(a);\n",
    "  a.click();\n",
    "  a.remove();\n",
    "  URL.revokeObjectURL(url);\n",
    "  setDirty(false);\n",
    "}\n",
    "function maybeNext(delta){\n",
    "  if(dirty){\n",
    "    const ok = confirm('You have unsaved changes. Click OK to continue without saving, or Cancel to stay and save.');\n",
    "    if(!ok) return;\n",
    "  }\n",
    "  idx = clamp(idx + delta, 0, LOCI.length-1);\n",
    "  selection = null;\n",
    "  render();\n",
    "}\n",
    "function drawLocus(L, splits){\n",
    "  const svg = byId('viz');\n",
    "  const w = svg.viewBox.baseVal.width || 900;\n",
    "  const h = svg.viewBox.baseVal.height || 180;\n",
    "  const ds = Number(L.DISPLAY_START);\n",
    "  const de = Number(L.DISPLAY_END);\n",
    "  while(svg.firstChild) svg.removeChild(svg.firstChild);\n",
    "  const ns='http://www.w3.org/2000/svg';\n",
    "  function line(x1,y1,x2,y2,cls){ const el=document.createElementNS(ns,'line'); el.setAttribute('x1',x1); el.setAttribute('y1',y1); el.setAttribute('x2',x2); el.setAttribute('y2',y2); if(cls) el.setAttribute('class',cls); svg.appendChild(el); return el; }\n",
    "  function rect(x,y,w,h,cls){ const el=document.createElementNS(ns,'rect'); el.setAttribute('x',x); el.setAttribute('y',y); el.setAttribute('width',w); el.setAttribute('height',h); if(cls) el.setAttribute('class',cls); svg.appendChild(el); return el; }\n",
    "  function text(x,y,t,cls){ const el=document.createElementNS(ns,'text'); el.setAttribute('x',x); el.setAttribute('y',y); el.textContent=t; if(cls) el.setAttribute('class',cls); svg.appendChild(el); return el; }\n",
    "  // title\n",
    "  text(40, 16, 'chr'+L.CHR+': '+ds+' - '+de, 'label');\n",
    "  // p-value track (max -log10(p) per window)\n",
    "  var pv = locusPvals(L.MERGED_LOCUS_ID);\n",
    "  if(pv && pv.length>0){\n",
    "    var maxY=0;\n",
    "    for(var i=0;i<pv.length;i++){ if(pv[i].MAX_NEG_LOG10_P>maxY) maxY=pv[i].MAX_NEG_LOG10_P; }\n",
    "    maxY = Math.max(1, maxY);\n",
    "    for(var i=0;i<pv.length;i++){\n",
    "      var rec=pv[i];\n",
    "      var x1=bpToX(Number(rec.WIN_START), ds, de, w);\n",
    "      var x2=bpToX(Number(rec.WIN_END), ds, de, w);\n",
    "      var v=Number(rec.MAX_NEG_LOG10_P);\n",
    "      if(!isFinite(v)) continue;\n",
    "      var barH = Math.round(30 * (v/maxY));\n",
    "      var yTop = 52 - barH;\n",
    "      var col = (COLORS && rec.SOURCE && COLORS[rec.SOURCE]) ? COLORS[rec.SOURCE] : '#444';\n",
    "      var r=rect(x1, yTop, Math.max(1, x2-x1), barH, null);\n",
    "      r.setAttribute('fill', col);\n",
    "      r.setAttribute('opacity', '0.35');\n",
    "    }\n",
    "    text(40, 52, 'max -log10(p) per '+WINDOW_BP+'bp window', 'small');\n",
    "  }\n",
    "  // axis for coordinate context\n",
    "  line(40, 128, w-40, 128, 'axis');\n",

    "  // original loci (colored by source)\n",
    "  var orig = locusOriginal(L.MERGED_LOCUS_ID);\n",
    "  if(orig && orig.length>0){\n",
    "    var maxOrig = Math.min(orig.length, 24);\n",
    "    for(var i=0;i<maxOrig;i++){\n",
    "      var o=orig[i];\n",
    "      var s=Number(o.START);\n",
    "      var e=Number(o.END);\n",
    "      if(!isFinite(s) || !isFinite(e)) continue;\n",
    "      var x1=bpToX(s, ds, de, w);\n",
    "      var x2=bpToX(e, ds, de, w);\n",
    "      var y=64 + (i%8)*7;\n",
    "      var col = (COLORS && o.SOURCE && COLORS[o.SOURCE]) ? COLORS[o.SOURCE] : '#888';\n",
    "      var rr=rect(x1, y, Math.max(1, x2-x1), 5, null);\n",
    "      rr.setAttribute('fill', col);\n",
    "      rr.setAttribute('opacity', '0.9');\n",
    "    }\n",
    "  }\n",
    "  // merged locus\n",
    "  const mx1 = bpToX(Number(L.MERGED_START), ds, de, w);\n",
    "  const mx2 = bpToX(Number(L.MERGED_END), ds, de, w);\n",
    "  rect(mx1, 110, Math.max(1, mx2-mx1), 40, 'merged');\n",
    "  // refined locus\n",
    "  if(validateInputs()){\n",
    "    const rx1 = bpToX(Number(L.REFINED_START), ds, de, w);\n",
    "    const rx2 = bpToX(Number(L.REFINED_END), ds, de, w);\n",
    "    rect(rx1, 100, Math.max(1, rx2-rx1), 8, 'refined');\n",
    "  }\n",
    "  // reported loci segments as small bars\n",
    "  const segs = locusSegments(L.MERGED_LOCUS_ID);\n",
    "  const maxBars = Math.min(segs.length, 30);\n",
    "  for(var i=0;i<maxBars;i++){\n",
    "    var s=segs[i];\n",
    "    if(s.START===null || s.END===null) continue;\n",
    "    var x1 = bpToX(Number(s.START), ds, de, w);\n",
    "    var x2 = bpToX(Number(s.END), ds, de, w);\n",
    "    rect(x1, 154 + (i%3)*6, Math.max(1, x2-x1), 4, 'seg');\n",
    "  }\n",
    "  // splits overlay\n",
    "  var _splits = splits || [];\n",
    "  for(var i=0;i<_splits.length;i++){\n",
    "    var sp = _splits[i];\n",
    "    const sx1 = bpToX(Number(sp.start), ds, de, w);\n",
    "    const sx2 = bpToX(Number(sp.end), ds, de, w);\n",
    "    rect(sx1, 68, Math.max(1, sx2-sx1), 44, 'split');\n",
    "  }\n",
    "  // selection overlay (if present)\n",
    "  if(selection){\n",
    "    const sx1 = bpToX(selection.start, ds, de, w);\n",
    "    const sx2 = bpToX(selection.end, ds, de, w);\n",
    "    rect(Math.min(sx1,sx2), 56, Math.abs(sx2-sx1), 120, 'selection');\n",
    "  }\n",
    "}\n",
    "function setupVizMouse(){\n",
    "  const svg=byId('viz');\n",
    "  let dragging=false;\n",
    "  let x0=null;\n",
    "  svg.addEventListener('mousedown', function(ev){\n",
    "    const L=current();\n",
    "    const w=svg.viewBox.baseVal.width||900;\n",
    "    const ds=Number(L.DISPLAY_START);\n",
    "    const de=Number(L.DISPLAY_END);\n",
    "    dragging=true;\n",
    "    x0 = ev.offsetX;\n",
    "    const bp0 = xToBp(x0, ds, de, w);\n",
    "    selection = {start: bp0, end: bp0};\n",
    "    drawLocus(L, parseSplits(L.SPLITS));\n",
    "  });\n",
    "  svg.addEventListener('mousemove', function(ev){\n",
    "    if(!dragging || !selection) return;\n",
    "    const L=current();\n",
    "    const w=svg.viewBox.baseVal.width||900;\n",
    "    const ds=Number(L.DISPLAY_START);\n",
    "    const de=Number(L.DISPLAY_END);\n",
    "    const bp = xToBp(ev.offsetX, ds, de, w);\n",
    "    selection.end = bp;\n",
    "    drawLocus(L, parseSplits(L.SPLITS));\n",
    "  });\n",
    "  window.addEventListener('mouseup', function(){ dragging=false; });\n",
    "}\n",
    "function addSplitFromSelection(){\n",
    "  if(!selection) { alert('Drag on the locus visualization to select a region first.'); return; }\n",
    "  const a = Math.min(selection.start, selection.end);\n",
    "  const b = Math.max(selection.start, selection.end);\n",
    "  const L=current();\n",
    "  const arr=parseSplits(L.SPLITS);\n",
    "  arr.push({start:a, end:b});\n",
    "  // sort and merge overlaps lightly\n",
    "  arr.sort(function(x,y){ return x.start - y.start; });\n",
    "  const merged=[];\n",
    "  for(var i=0;i<arr.length;i++){\n",
    "    var s = arr[i];\n",
    "    if(merged.length===0) { merged.push({start:s.start,end:s.end}); continue; }\n",
    "    var last=merged[merged.length-1];\n",
    "    if(s.start <= last.end){ last.end = Math.max(last.end, s.end); } else { merged.push({start:s.start,end:s.end}); }\n",
    "  }\n",
    "  L.SPLITS=splitsToString(merged);\n",
    "  selection=null;\n",
    "  setDirty(true);\n",
    "  renderSplitsList(merged);\n",
    "  drawLocus(L, merged);\n",
    "}\n",
    "byId('prev').addEventListener('click', function(){ maybeNext(-1); });\n",
    "byId('next').addEventListener('click', function(){ maybeNext(1); });\n",
    "byId('save').addEventListener('click', downloadTSV);\n",
    "byId('add_split').addEventListener('click', addSplitFromSelection);\n",
    "var ids=['refined_start','refined_end','locus_name','keep'];\n",
    "for(var i=0;i<ids.length;i++){\n",
    "  var el=byId(ids[i]);\n",
    "  (function(_el){\n",
    "    _el.addEventListener('input', function(){ updateFromInputs(); validateInputs(); });\n",
    "    _el.addEventListener('change', function(){ updateFromInputs(); validateInputs(); });\n",
    "  })(el);\n",
    "}\n",
    "setupVizMouse();\n",
    "render();\n",
    "  } catch(e) { showError(e && e.stack ? e.stack : e); }\n",
    "})();\n"
  )

  # Make the embedded JS older-browser friendly (ES5). This avoids silent failures
  # in older embedded viewers that don't support const/let.
  js <- gsub("\\bconst\\b", "var", js)
  js <- gsub("\\blet\\b", "var", js)

  # Let save_html create the outer <html> wrapper; provide head/body only.
  page <- htmltools::tagList(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title("Locus wizard (editable HTML)"),
      htmltools::tags$style(htmltools::HTML(
        "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif; margin: 18px; }\n",
        ".row { display: flex; gap: 18px; align-items: flex-start; }\n",
        ".panel { border: 1px solid #ddd; padding: 12px; border-radius: 6px; }\n",
        ".left { flex: 1 1 auto; min-width: 520px; }\n",
        ".right { width: 360px; }\n",
        ".small { color: #555; font-size: 12px; }\n",
        ".controls { display:flex; gap: 8px; align-items: center; margin-bottom: 10px; }\n",
        ".btn { padding: 8px 10px; cursor: pointer; }\n",
        ".btn-small { padding: 2px 6px; font-size: 12px; cursor: pointer; }\n",
        ".badge { padding: 2px 8px; border-radius: 10px; font-size: 12px; }\n",
        ".badge-ok { background: #e6ffed; color: #0a5; border: 1px solid #bfecc9; }\n",
        ".badge-warn { background: #fff4e5; color: #a60; border: 1px solid #ffd7a6; }\n",
        "label { display:block; margin-top: 8px; font-weight: 600; }\n",
        "input[type=text], input[type=number] { width: 100%; box-sizing: border-box; padding: 6px; }\n",
        "input.invalid { outline: 2px solid #cc0000; }\n",
        "ul { padding-left: 18px; }\n",
        "svg { width: 100%; height: 220px; border: 1px solid #ddd; border-radius: 6px; background: white; }\n",
        ".axis { stroke: #333; stroke-width: 2; }\n",
        ".label { fill: #111; font-size: 12px; }\n",
        ".small { fill: #555; font-size: 11px; }\n",
        ".merged { fill: #e6f2ff; stroke: #339; stroke-width: 1; }\n",
        ".refined { fill: #2b6cb0; opacity: 0.9; }\n",
        ".seg { fill: #444; opacity: 0.55; }\n",
        ".split { fill: #fef3c7; stroke: #f59e0b; stroke-width: 1; opacity: 0.8; }\n",
        ".selection { fill: #c7d2fe; opacity: 0.5; }\n"
      ))
    ),
    htmltools::tags$body(
      htmltools::tags$h1("Locus wizard (one locus at a time)"),
      htmltools::tags$pre(
        id = "js_error",
        style = "display:none; white-space: pre-wrap; background:#fff5f5; border:1px solid #feb2b2; padding:10px; border-radius:6px; margin:10px 0;",
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
      htmltools::tags$div(class = "row",
        htmltools::tags$div(class = "panel left",
          htmltools::tags$div(class = "small",
            htmltools::tags$b("MERGED_LOCUS_ID: "), htmltools::tags$span(id = "merged_id"),
            " | ", htmltools::tags$b("chr"), htmltools::tags$span(id = "chr"),
            " | ", htmltools::tags$b("rows: "), htmltools::tags$span(id = "rows"),
            " | ", htmltools::tags$b("sources: "), htmltools::tags$span(id = "sources")
          ),
          htmltools::tags$svg(id = "viz", viewBox = "0 0 900 180"),
          htmltools::tags$div(style = "margin-top:10px;",
            htmltools::tags$button(id = "add_split", class = "btn", "Add split from selection"),
            htmltools::tags$span(class = "small", " (drag to select region on plot)")
          ),
          htmltools::tags$div(class = "small", style = "margin-top:8px;",
            "Merged: ", htmltools::tags$span(id = "merged_start"), " - ", htmltools::tags$span(id = "merged_end")
          )
        ),
        htmltools::tags$div(class = "panel right",
          htmltools::tags$h3("Edits"),
          htmltools::tags$label("REFINED_START"),
          htmltools::tags$input(id = "refined_start", type = "number", step = "1"),
          htmltools::tags$label("REFINED_END"),
          htmltools::tags$input(id = "refined_end", type = "number", step = "1"),
          htmltools::tags$label("LOCUS_NAME"),
          htmltools::tags$input(id = "locus_name", type = "text"),
          htmltools::tags$label("KEEP"),
          htmltools::tags$input(id = "keep", type = "checkbox"),
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

  out_dir <- dirname(output_html)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  htmltools::save_html(page, file = output_html)

  message(sprintf("Locus wizard HTML saved to: %s", output_html))
  invisible(list(html = output_html))
}
