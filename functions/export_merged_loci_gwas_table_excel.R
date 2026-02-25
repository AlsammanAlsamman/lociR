#' Export SNP-level merged GWAS table for merged loci
#'
#' Builds an SNP-level table by collecting all GWAS SNPs that fall inside the
#' merged loci intervals, merging across datasets by chromosome/position.
#' Output is written to Excel with:
#' - locus metadata and reported locus
#' - gene assignment (inside gene, otherwise nearest gene)
#' - distance from SNP to assigned gene (0 if inside)
#' - per-dataset locus index columns (green fill when non-empty)
#' - per-dataset P/BETA/OR/SE columns
#' - p-value cells highlighted light red when p < 5e-8
#' Missing values are exported as "-".
#'
#' @param aligned_loci data.frame returned by `align_with_reported()`
#' @param gwas_list named list of `GWAS` objects returned by `read_gwas()`
#' @param output_file character. Output Excel filepath.
#' @param annotation_file character. Gene annotation file path.
#'   Default: "inputs/NCBI37.3.gene.loc".
#' @param apply_styles logical. If TRUE, apply green/red cell background styling.
#'   Default: FALSE for faster export on large tables.
#'
#' @return Invisible output file path.
#' @export
export_merged_loci_gwas_table_excel <- function(aligned_loci,
                                                gwas_list,
                                                output_file,
                                                annotation_file = "inputs/NCBI37.3.gene.loc",
                                                apply_styles = FALSE) {
  start_time <- Sys.time()
  if (missing(output_file) || !nzchar(trimws(output_file))) {
    stop("output_file is required (e.g., 'outputs/merged_loci_gwas_table.xlsx').")
  }
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("openxlsx2 package is required. Install with: install.packages('openxlsx2')")
  }
  if (!is.data.frame(aligned_loci) || !inherits(aligned_loci, "AlignedLoci")) {
    stop("aligned_loci must be the data frame produced by align_with_reported() (class 'AlignedLoci').")
  }
  if (nrow(aligned_loci) == 0) {
    stop("aligned_loci is empty.")
  }
  if (missing(gwas_list) || is.null(gwas_list) || length(gwas_list) == 0) {
    stop("gwas_list is required and must contain at least one GWAS dataset.")
  }

  `%||%` <- function(x, y) if (!is.null(x)) x else y

  normalize_chr <- function(chr) {
    out <- trimws(as.character(chr))
    out <- gsub("^chr", "", out, ignore.case = TRUE)

    is_numeric_like <- grepl("^[0-9]+$", out)
    numeric_part <- suppressWarnings(as.integer(out))

    out[is_numeric_like & !is.na(numeric_part) & numeric_part == 23L] <- "X"
    out[is_numeric_like & !is.na(numeric_part) & numeric_part == 24L] <- "Y"
    out[is_numeric_like & !is.na(numeric_part) & numeric_part %in% c(25L, 26L)] <- "MT"

    out[toupper(out) %in% c("X")] <- "X"
    out[toupper(out) %in% c("Y")] <- "Y"
    out[toupper(out) %in% c("M", "MT", "CHRM", "MTR")] <- "MT"

    out
  }

  sanitize_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

  format_num <- function(x, scientific = FALSE, digits = 4) {
    if (!is.finite(x)) return("-")
    if (scientific) return(format(x, scientific = TRUE, digits = digits))
    format(round(x, digits), trim = TRUE, scientific = FALSE)
  }

  split_sources <- function(x) {
    if (is.na(x) || !nzchar(trimws(x))) return(character(0))
    vals <- trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
    unique(vals[nzchar(vals)])
  }

  has_source <- function(source_vec, dataset_name) {
    if (length(source_vec) == 0 || !nzchar(dataset_name)) return(FALSE)
    any(tolower(source_vec) == tolower(dataset_name))
  }

  read_gene_annotation <- function(path) {
    if (!file.exists(path)) {
      warning("Annotation file not found: ", path)
      return(data.frame(CHR = character(0), START = numeric(0), END = numeric(0), GENE = character(0)))
    }

    ann <- tryCatch(
      data.table::fread(path, header = FALSE, fill = TRUE, data.table = FALSE),
      error = function(e) NULL
    )

    if (is.null(ann) || nrow(ann) == 0 || ncol(ann) < 4) {
      warning("Could not parse annotation file: ", path)
      return(data.frame(CHR = character(0), START = numeric(0), END = numeric(0), GENE = character(0)))
    }

    if (ncol(ann) >= 6) {
      chr_col <- ann[[2]]
      start_col <- ann[[3]]
      end_col <- ann[[4]]
      gene_col <- ann[[6]]
    } else {
      chr_col <- ann[[2]]
      start_col <- ann[[3]]
      end_col <- ann[[4]]
      gene_col <- ann[[1]]
    }

    out <- data.frame(
      CHR = normalize_chr(chr_col),
      START = suppressWarnings(as.numeric(start_col)),
      END = suppressWarnings(as.numeric(end_col)),
      GENE = as.character(gene_col),
      stringsAsFactors = FALSE
    )

    out <- out[is.finite(out$START) & is.finite(out$END) & nzchar(trimws(out$GENE)), , drop = FALSE]
    out
  }

  assign_gene <- function(chr, pos, ann_by_chr) {
    if (!nzchar(chr) || !is.finite(pos)) {
      return(list(gene = "-", distance = "-"))
    }
    chr_ann <- ann_by_chr[[as.character(chr)]]
    if (is.null(chr_ann) || nrow(chr_ann) == 0) {
      return(list(gene = "-", distance = "-"))
    }

    hit <- chr_ann$START <= pos & chr_ann$END >= pos
    hit_genes <- unique(chr_ann$GENE[hit])
    hit_genes <- hit_genes[nzchar(trimws(hit_genes))]
    if (length(hit_genes) > 0) {
      return(list(gene = paste(hit_genes, collapse = "; "), distance = "0"))
    }

    dist_to_gene <- ifelse(
      pos < chr_ann$START,
      chr_ann$START - pos,
      ifelse(pos > chr_ann$END, pos - chr_ann$END, 0)
    )
    min_dist <- suppressWarnings(min(dist_to_gene, na.rm = TRUE))
    if (!is.finite(min_dist)) {
      return(list(gene = "-", distance = "-"))
    }

    nearest_genes <- unique(chr_ann$GENE[is.finite(dist_to_gene) & dist_to_gene == min_dist])
    nearest_genes <- nearest_genes[nzchar(trimws(nearest_genes))]
    list(
      gene = if (length(nearest_genes) > 0) paste(nearest_genes, collapse = "; ") else "-",
      distance = as.character(as.integer(round(min_dist)))
    )
  }

  add_fill_runs <- function(wb, sheet, col_idx, row_idx_vec, color_obj) {
    if (length(row_idx_vec) == 0 || !is.finite(col_idx)) return(invisible(NULL))
    row_idx_vec <- sort(unique(as.integer(row_idx_vec)))
    if (length(row_idx_vec) == 0) return(invisible(NULL))

    run_start <- row_idx_vec[[1]]
    run_prev <- row_idx_vec[[1]]

    flush_run <- function(s, e) {
      col_letter <- openxlsx2::int2col(col_idx)
      dims <- if (s == e) {
        paste0(col_letter, s + 1L)
      } else {
        paste0(col_letter, s + 1L, ":", col_letter, e + 1L)
      }
      wb$add_fill(sheet = sheet, dims = dims, color = color_obj)
    }

    if (length(row_idx_vec) == 1) {
      flush_run(run_start, run_prev)
      return(invisible(NULL))
    }

    for (k in 2:length(row_idx_vec)) {
      cur <- row_idx_vec[[k]]
      if (cur == run_prev + 1L) {
        run_prev <- cur
      } else {
        flush_run(run_start, run_prev)
        run_start <- cur
        run_prev <- cur
      }
    }
    flush_run(run_start, run_prev)
    invisible(NULL)
  }

  # Validate required aligned columns
  req_aligned <- c("MERGED_LOCUS_ID", "CHR", "MERGED_START", "MERGED_END", "SOURCES")
  missing_aligned <- setdiff(req_aligned, colnames(aligned_loci))
  if (length(missing_aligned) > 0) {
    stop("aligned_loci is missing required columns: ", paste(missing_aligned, collapse = ", "))
  }

  # Build merged-locus metadata (one row per merged locus)
  df <- as.data.frame(aligned_loci, stringsAsFactors = FALSE)
  df$MERGED_LOCUS_ID <- as.character(df$MERGED_LOCUS_ID)
  df$CHR <- normalize_chr(df$CHR)
  df$MERGED_START <- suppressWarnings(as.numeric(df$MERGED_START))
  df$MERGED_END <- suppressWarnings(as.numeric(df$MERGED_END))
  if (!"REPORTED_LOCUS" %in% colnames(df)) df$REPORTED_LOCUS <- NA_character_

  merged_ids <- unique(df$MERGED_LOCUS_ID)
  merged_info <- do.call(rbind, lapply(merged_ids, function(mid) {
    sub <- df[df$MERGED_LOCUS_ID == mid, , drop = FALSE]
    rep_vals <- unique(trimws(as.character(sub$REPORTED_LOCUS)))
    rep_vals <- rep_vals[!is.na(rep_vals) & nzchar(rep_vals)]
    rep_label <- if (length(rep_vals) > 0) paste(rep_vals, collapse = "; ") else "Novel"

    data.frame(
      MERGED_LOCUS_ID = mid,
      CHR = as.character(sub$CHR[[1]]),
      MERGED_START = as.numeric(sub$MERGED_START[[1]]),
      MERGED_END = as.numeric(sub$MERGED_END[[1]]),
      REPORTED_LOCUS = rep_label,
      SOURCES = as.character(sub$SOURCES[[1]]),
      stringsAsFactors = FALSE
    )
  }))

  # Normalize and validate GWAS list
  if (inherits(gwas_list, "GWAS")) {
    gwas_list <- list(gwas_list)
    names(gwas_list) <- gwas_list[[1]]$name %||% "GWAS"
  }
  if (is.null(names(gwas_list)) || any(!nzchar(trimws(names(gwas_list))))) {
    stop("gwas_list must be a named list. Each dataset needs a unique name.")
  }

  dataset_names <- names(gwas_list)
  gwas_lookup <- list()
  gwas_by_chr <- list()
  message(sprintf("Preparing GWAS datasets (%d)...", length(dataset_names)))
  for (nm in dataset_names) {
    obj <- gwas_list[[nm]]
    dat <- if (is.list(obj) && !is.null(obj$data)) obj$data else NULL
    if (is.null(dat) || !is.data.frame(dat)) {
      gwas_lookup[[nm]] <- data.frame()
      gwas_by_chr[[nm]] <- list()
      next
    }

    dat <- as.data.frame(dat, stringsAsFactors = FALSE)
    req_gwas <- c("CHR", "POS", "P")
    if (!all(req_gwas %in% colnames(dat))) {
      gwas_lookup[[nm]] <- data.frame()
      gwas_by_chr[[nm]] <- list()
      next
    }

    dat$CHR <- normalize_chr(dat$CHR)
    dat$POS <- suppressWarnings(as.numeric(dat$POS))
    dat$P <- suppressWarnings(as.numeric(dat$P))
    if ("BETA" %in% colnames(dat)) dat$BETA <- suppressWarnings(as.numeric(dat$BETA))
    if ("OR" %in% colnames(dat)) dat$OR <- suppressWarnings(as.numeric(dat$OR))
    if ("SE" %in% colnames(dat)) dat$SE <- suppressWarnings(as.numeric(dat$SE))

    gwas_lookup[[nm]] <- dat
    gwas_by_chr[[nm]] <- split(dat, dat$CHR)
  }

  # Determine sequential locus index in each dataset (present loci only)
  source_lists <- lapply(merged_info$SOURCES, split_sources)
  locus_index_map <- setNames(vector("list", length(dataset_names)), dataset_names)
  for (ds in dataset_names) {
    present <- vapply(source_lists, has_source, logical(1), dataset_name = ds)
    idx_vec <- rep(NA_integer_, nrow(merged_info))
    n_present <- sum(present)
    if (n_present > 0) idx_vec[present] <- seq_len(n_present)
    names(idx_vec) <- merged_info$MERGED_LOCUS_ID
    locus_index_map[[ds]] <- idx_vec
  }

  source_cols <- paste0("SOURCE_IDX_", sanitize_name(dataset_names))
  p_cols <- paste0("P_", sanitize_name(dataset_names))
  beta_cols <- paste0("BETA_", sanitize_name(dataset_names))
  or_cols <- paste0("OR_", sanitize_name(dataset_names))
  se_cols <- paste0("SE_", sanitize_name(dataset_names))

  col_order <- c(
    "CHR", "POS", "MERGED_LOCUS_ID", "MERGED_START", "MERGED_END", "REPORTED_LOCUS",
    "GENE", "DISTANCE_TO_GENE",
    source_cols,
    p_cols,
    beta_cols,
    or_cols,
    se_cols
  )

  ann <- read_gene_annotation(annotation_file)
  ann_by_chr <- split(ann, ann$CHR)
  message(sprintf("Building merged SNP table across %d loci...", nrow(merged_info)))

  # Build SNP-level merged rows
  out_rows <- list()
  row_counter <- 1L

  for (i in seq_len(nrow(merged_info))) {
    mid <- merged_info$MERGED_LOCUS_ID[[i]]
    chr <- merged_info$CHR[[i]]
    st <- merged_info$MERGED_START[[i]]
    en <- merged_info$MERGED_END[[i]]
    rep_locus <- merged_info$REPORTED_LOCUS[[i]]

    if (i %% 50 == 0 || i == nrow(merged_info)) {
      message(sprintf("  Locus progress: %d/%d", i, nrow(merged_info)))
    }

    ds_slices <- list()
    ds_best_map <- list()
    for (ds in dataset_names) {
      chr_bucket <- gwas_by_chr[[ds]][[as.character(chr)]]
      if (is.null(chr_bucket) || nrow(chr_bucket) == 0) {
        ds_slices[[ds]] <- data.frame()
        ds_best_map[[ds]] <- data.frame()
        next
      }

      keep <- is.finite(chr_bucket$POS) & chr_bucket$POS >= st & chr_bucket$POS <= en
      ds_dat <- chr_bucket[keep, , drop = FALSE]
      ds_slices[[ds]] <- ds_dat

      if (nrow(ds_dat) == 0) {
        ds_best_map[[ds]] <- data.frame()
        next
      }

      ord <- order(ds_dat$POS, ifelse(is.finite(ds_dat$P), ds_dat$P, Inf))
      ds_sorted <- ds_dat[ord, , drop = FALSE]
      ds_best_map[[ds]] <- ds_sorted[!duplicated(ds_sorted$POS), , drop = FALSE]
    }

    pos_union <- sort(unique(unlist(lapply(ds_slices, function(d) if (is.data.frame(d) && nrow(d) > 0) d$POS else numeric(0)))))
    pos_union <- pos_union[is.finite(pos_union)]
    if (length(pos_union) == 0) next

    for (pos in pos_union) {
      row <- as.list(setNames(rep("-", length(col_order)), col_order))

      row$CHR <- as.character(chr)
      row$POS <- as.character(as.integer(round(pos)))
      row$MERGED_LOCUS_ID <- as.character(mid)
      row$MERGED_START <- as.character(as.integer(round(st)))
      row$MERGED_END <- as.character(as.integer(round(en)))
      row$REPORTED_LOCUS <- as.character(rep_locus)

      gene_info <- assign_gene(chr = chr, pos = pos, ann_by_chr = ann_by_chr)
      row$GENE <- gene_info$gene
      row$DISTANCE_TO_GENE <- gene_info$distance

      for (j in seq_along(dataset_names)) {
        ds <- dataset_names[[j]]
        ds_tag <- sanitize_name(ds)

        idx_val <- locus_index_map[[ds]][[as.character(mid)]]
        row[[paste0("SOURCE_IDX_", ds_tag)]] <- if (is.finite(idx_val)) as.character(idx_val) else "-"

        best_tbl <- ds_best_map[[ds]]
        if (is.null(best_tbl) || nrow(best_tbl) == 0) next
        k <- match(pos, best_tbl$POS)
        if (is.na(k)) next
        best <- best_tbl[k, , drop = FALSE]

        p_val <- if ("P" %in% colnames(best)) best$P[[1]] else NA_real_
        beta_val <- if ("BETA" %in% colnames(best)) best$BETA[[1]] else NA_real_
        or_val <- if ("OR" %in% colnames(best)) best$OR[[1]] else NA_real_
        se_val <- if ("SE" %in% colnames(best)) best$SE[[1]] else NA_real_

        row[[paste0("P_", ds_tag)]] <- format_num(p_val, scientific = TRUE, digits = 4)
        row[[paste0("BETA_", ds_tag)]] <- format_num(beta_val, scientific = FALSE, digits = 6)
        row[[paste0("OR_", ds_tag)]] <- format_num(or_val, scientific = FALSE, digits = 6)
        row[[paste0("SE_", ds_tag)]] <- format_num(se_val, scientific = FALSE, digits = 6)
      }

      out_rows[[row_counter]] <- as.data.frame(row, stringsAsFactors = FALSE)
      row_counter <- row_counter + 1L
    }
  }

  if (length(out_rows) == 0) {
    stop("No GWAS SNPs found inside merged loci for the provided datasets.")
  }

  out_table <- do.call(rbind, out_rows)
  out_table <- out_table[, col_order, drop = FALSE]

  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("GWAS_Merged")
  wb$add_data(sheet = "GWAS_Merged", x = out_table, start_row = 1, start_col = 1)
  if (isTRUE(apply_styles)) {
    message("Applying Excel styling...")

    # Green fill for non-empty source index cells
    green_fill <- openxlsx2::wb_color(hex = "FF92D050")
    for (src_col in source_cols) {
      col_idx <- which(colnames(out_table) == src_col)
      if (length(col_idx) != 1) next
      rows_hit <- which(out_table[[src_col]] != "-" & nzchar(out_table[[src_col]]))
      if (length(rows_hit) == 0) next
      add_fill_runs(wb, "GWAS_Merged", col_idx, rows_hit, green_fill)
    }

    # Light red fill for significant p-values (< 5e-8)
    red_fill <- openxlsx2::wb_color(hex = "FFFFC7CE")
    for (p_col in p_cols) {
      col_idx <- which(colnames(out_table) == p_col)
      if (length(col_idx) != 1) next
      pvals <- suppressWarnings(as.numeric(out_table[[p_col]]))
      rows_hit <- which(is.finite(pvals) & pvals < 5e-8)
      if (length(rows_hit) == 0) next
      add_fill_runs(wb, "GWAS_Merged", col_idx, rows_hit, red_fill)
    }
  }

  wb$set_col_widths(sheet = "GWAS_Merged", cols = 1:ncol(out_table), widths = "auto")
  wb$save(output_file)

  message(sprintf("Merged GWAS SNP table saved: %s", output_file))
  message(sprintf("  Rows (SNP-locus pairs): %d", nrow(out_table)))
  message(sprintf("  Datasets: %d", length(dataset_names)))
  message(sprintf("  Elapsed: %.1f sec", as.numeric(difftime(Sys.time(), start_time, units = "secs"))))

  invisible(output_file)
}
