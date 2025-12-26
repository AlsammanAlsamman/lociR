#' Export aligned loci to an editable, JS-enabled HTML table
#'
#' This exporter generates a static HTML file that renders an interactive table
#' in the browser and allows editing a small set of columns. Because the HTML
#' is static, edits cannot be written back to disk automatically; instead, the
#' page provides a "Download TSV" button to save your edits.
#'
#' Recommended workflow:
#' 1) Generate the HTML.
#' 2) Open it in a browser, edit `REFINED_START`, `REFINED_END`, `LOCUS_NAME`, `KEEP`.
#' 3) Click "Download TSV" and use the downloaded TSV as your authoritative artifact.
#'
#' @param aligned_loci data.frame returned by `align_with_reported()` (class 'AlignedLoci')
#' @param output_html character. Output HTML filepath (e.g., 'outputs/aligned_editable.html')
#' @param flank_bp numeric. Base-pairs to extend the display window beyond the merged locus.
#' @param n_bins integer. Number of bins (saved into the table for reproducibility).
#'
#' @return Named list with path: `list(html = <path>)`.
#' @export
export_aligned_to_html_editable <- function(aligned_loci,
                                            output_html,
                                            flank_bp = 20000,
                                            n_bins = 20) {
  if (missing(output_html) || !nzchar(trimws(output_html))) {
    stop("output_html is required (e.g., 'outputs/aligned_editable.html').")
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

  df <- aligned_loci

  # Derived columns used for review
  df$DISPLAY_START <- pmax(0, df$MERGED_START - flank_bp)
  df$DISPLAY_END <- df$MERGED_END + flank_bp
  df$FLANK_BP <- flank_bp
  df$N_BINS <- n_bins

  # Editable columns (default to merged boundaries)
  if (!"REFINED_START" %in% colnames(df)) df$REFINED_START <- df$MERGED_START
  if (!"REFINED_END" %in% colnames(df)) df$REFINED_END <- df$MERGED_END
  if (!"LOCUS_NAME" %in% colnames(df)) df$LOCUS_NAME <- ""
  if (!"KEEP" %in% colnames(df)) df$KEEP <- TRUE

  # Keep a stable, explicit column order (editables near the front)
  out_cols <- c(
    "MERGED_LOCUS_ID",
    "CHR",
    "MERGED_START",
    "MERGED_END",
    "REFINED_START",
    "REFINED_END",
    "LOCUS_NAME",
    "KEEP",
    "DISPLAY_START",
    "DISPLAY_END",
    "FLANK_BP",
    "N_BINS",
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
  out_cols <- out_cols[out_cols %in% colnames(df)]
  df <- df[, out_cols, drop = FALSE]

  # Build a self-contained HTML table with embedded data (offline; no CDN)
  df <- df[order(df$MERGED_LOCUS_ID, df$ALIGNED_LOCUS_ID), , drop = FALSE]

  esc_txt <- function(x) {
    x <- ifelse(is.na(x), "", as.character(x))
    x
  }

  editable_int_cols <- c("REFINED_START", "REFINED_END")
  editable_text_cols <- c("LOCUS_NAME")
  editable_bool_cols <- c("KEEP")

  ths <- lapply(names(df), function(nm) htmltools::tags$th(nm))

  make_td <- function(field, value) {
    v <- esc_txt(value)
    if (field %in% editable_int_cols) {
      return(htmltools::tags$td(`data-field` = field,
                               htmltools::tags$input(
                                 type = "number",
                                 step = "1",
                                 value = v,
                                 class = "edit edit-int",
                                 `data-field` = field
                               )))
    }
    if (field %in% editable_text_cols) {
      return(htmltools::tags$td(`data-field` = field,
                               htmltools::tags$input(
                                 type = "text",
                                 value = v,
                                 class = "edit edit-text",
                                 `data-field` = field
                               )))
    }
    if (field %in% editable_bool_cols) {
      checked <- isTRUE(value)
      return(htmltools::tags$td(`data-field` = field,
                               htmltools::tags$input(
                                 type = "checkbox",
                                 class = "edit edit-bool",
                                 `data-field` = field,
                                 checked = if (checked) "checked" else NULL
                               )))
    }
    htmltools::tags$td(`data-field` = field, htmltools::tags$span(v))
  }

  rows <- vector("list", length = 0)
  last_merged <- NULL
  ncol_tbl <- ncol(df)

  for (i in seq_len(nrow(df))) {
    merged_id <- df$MERGED_LOCUS_ID[[i]]
    if (!identical(merged_id, last_merged)) {
      rows <- c(rows, list(
        htmltools::tags$tr(
          class = "group-row",
          htmltools::tags$td(colspan = ncol_tbl, paste0("MERGED_LOCUS_ID: ", esc_txt(merged_id)))
        )
      ))
      last_merged <- merged_id
    }

    status <- tolower(esc_txt(df$ALIGNMENT_STATUS[[i]]))
    row_class <- ""
    if (grepl("matched", status)) row_class <- "row-matched"
    if (grepl("unmatched", status)) row_class <- "row-unmatched"
    if (grepl("overlap", status)) row_class <- "row-overlap"

    tds <- Map(make_td, names(df), df[i, , drop = FALSE])
    rows <- c(rows, list(htmltools::tags$tr(class = row_class, tds)))
  }

  table_tag <- htmltools::tags$table(
    id = "loci-table",
    class = "loci-table",
    htmltools::tags$thead(htmltools::tags$tr(ths)),
    htmltools::tags$tbody(rows)
  )

  download_name <- paste0(tools::file_path_sans_ext(basename(output_html)), ".tsv")

  js <- paste0(
    "const DOWNLOAD_NAME = ", dQuote(download_name), ";\n",
    "function escTSV(v){\n",
    "  if(v===null || v===undefined) return '';\n",
    "  return String(v).replace(/\t/g,' ').replace(/\r?\n/g,' ');\n",
    "}\n",
    "function getCellValue(td){\n",
    "  const inp = td.querySelector('input');\n",
    "  if(!inp) return td.textContent.trim();\n",
    "  if(inp.type==='checkbox') return inp.checked ? 'TRUE' : 'FALSE';\n",
    "  return inp.value;\n",
    "}\n",
    "function validate(){\n",
    "  let bad = 0;\n",
    "  const rows = document.querySelectorAll('#loci-table tbody tr');\n",
    "  rows.forEach(r=>{\n",
    "    if(r.classList.contains('group-row')) return;\n",
    "    const s = r.querySelector('td[data-field=\"REFINED_START\"] input');\n",
    "    const e = r.querySelector('td[data-field=\"REFINED_END\"] input');\n",
    "    if(!s || !e) return;\n",
    "    const sv = s.value.trim();\n",
    "    const ev = e.value.trim();\n",
    "    const si = Number(sv);\n",
    "    const ei = Number(ev);\n",
    "    const ok = Number.isInteger(si) && Number.isInteger(ei) && si <= ei;\n",
    "    s.classList.toggle('invalid', !ok);\n",
    "    e.classList.toggle('invalid', !ok);\n",
    "    if(!ok) bad++;\n",
    "  });\n",
    "  return bad;\n",
    "}\n",
    "function downloadTSV(){\n",
    "  const bad = validate();\n",
    "  if(bad>0){\n",
    "    alert('There are '+bad+' row(s) with invalid REFINED_START/REFINED_END (must be integers and start<=end). Fix them before using the TSV.');\n",
    "  }\n",
    "  const header = Array.from(document.querySelectorAll('#loci-table thead th')).map(th=>th.textContent.trim());\n",
    "  const lines = [];\n",
    "  lines.push(header.join('\t'));\n",
    "  const rows = document.querySelectorAll('#loci-table tbody tr');\n",
    "  rows.forEach(r=>{\n",
    "    if(r.classList.contains('group-row')) return;\n",
    "    const tds = r.querySelectorAll('td');\n",
    "    const vals = Array.from(tds).map(td=>escTSV(getCellValue(td)));\n",
    "    lines.push(vals.join('\t'));\n",
    "  });\n",
    "  const blob = new Blob([lines.join('\\n')], {type:'text/tab-separated-values;charset=utf-8'});\n",
    "  const url = URL.createObjectURL(blob);\n",
    "  const a = document.createElement('a');\n",
    "  a.href = url;\n",
    "  a.download = DOWNLOAD_NAME;\n",
    "  document.body.appendChild(a);\n",
    "  a.click();\n",
    "  a.remove();\n",
    "  URL.revokeObjectURL(url);\n",
    "}\n",
    "document.getElementById('download-tsv').addEventListener('click', downloadTSV);\n",
    "document.querySelectorAll('#loci-table input').forEach(inp=>{\n",
    "  inp.addEventListener('input', ()=>{ validate(); });\n",
    "  inp.addEventListener('change', ()=>{ validate(); });\n",
    "});\n",
    "validate();\n"
  )

  page <- htmltools::tags$html(
    htmltools::tags$head(
      htmltools::tags$meta(charset = "utf-8"),
      htmltools::tags$title("Aligned loci (editable HTML)"),
      htmltools::tags$style(htmltools::HTML(
        "body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Arial, sans-serif; margin: 18px; }\n",
        ".small { color: #555; }\n",
        ".loci-table { border-collapse: collapse; width: 100%; margin-top: 12px; }\n",
        ".loci-table th, .loci-table td { border: 1px solid #ddd; padding: 6px 8px; vertical-align: top; }\n",
        ".loci-table th { position: sticky; top: 0; background: #fafafa; z-index: 2; }\n",
        ".loci-table input.edit { width: 100%; box-sizing: border-box; }\n",
        ".loci-table input.edit-bool { width: auto; }\n",
        ".group-row td { background: #f2f2f2; font-weight: 600; }\n",
        ".row-matched { background: #f0fff4; }\n",
        ".row-unmatched { background: #fff5f5; }\n",
        ".row-overlap { background: #f7fafc; }\n",
        "input.invalid { outline: 2px solid #cc0000; }\n",
        "button { padding: 8px 12px; cursor: pointer; }\n"
      ))
    ),
    htmltools::tags$body(
      htmltools::tags$h1("Aligned loci (editable HTML)"),
      htmltools::tags$p(class = "small",
        "This is an interactive table. Edit ONLY: REFINED_START, REFINED_END, LOCUS_NAME, KEEP. ",
        "Edits are not written back to disk automatically; use Download TSV to save."),
      htmltools::tags$ul(
        htmltools::tags$li(paste0("Rows: ", nrow(df))),
        htmltools::tags$li(paste0("Merged loci: ", length(unique(df$MERGED_LOCUS_ID)))),
        htmltools::tags$li(paste0("Display flank (bp): ", flank_bp, "; bins: ", n_bins))
      ),
      htmltools::tags$button(id = "download-tsv", "Download TSV"),
      table_tag,
      htmltools::tags$script(htmltools::HTML(js))
    )
  )

  out_dir <- dirname(output_html)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  htmltools::save_html(page, file = output_html)

  message(sprintf("Editable HTML saved to: %s", output_html))
  invisible(list(html = output_html))
}
