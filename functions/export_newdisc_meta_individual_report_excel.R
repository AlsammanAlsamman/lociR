#' Export NewDisc meta + individual follow-up report to Excel
#'
#' Workflow:
#' 1) Read loci (`newDisc`-style) and GWAS files with `data.table::fread`.
#' 2) Subset each GWAS file to loci intervals (CHR/START/END) and cache subsets.
#' 3) If cached subset files already exist, reuse them to save time.
#' 4) For each locus:
#'    - count meta SNPs with p-value below threshold,
#'    - pick lead SNP from meta (lowest p) and include full meta row,
#'    - add closest gene to the lead SNP position,
#'    - lookup the same SNP position in each individual GWAS and report p + OR
#'      converted from beta (`exp(beta)`),
#'    - annotate gene (inside-gene, else nearest gene),
#'    - compute distance to previous/next locus,
#'    - assign nearest reported locus when boundary distance <= threshold.
#'
#' @param loci_file character. Loci file path (default: inputs/newDisc.txt).
#' @param meta_file character. Meta GWAS file path (default: inputs/GWAS/disc.tsv).
#' @param individual_files character vector. Individual GWAS file paths.
#' @param output_file character. Output Excel file.
#' @param subset_root_dir character. Root directory for cached subsets.
#' @param p_threshold numeric. P-value threshold for per-locus SNP count.
#' @param annotation_file character. Gene annotation file path.
#' @param reported_file character. Reported loci file path (mehdi-style).
#' @param reported_max_distance_bp numeric. Max boundary distance for matching.
#' @param individual_names optional character vector names for individual files.
#' @param force_rebuild_subset logical. Rebuild cached subsets even if present.
#' @param summary_sheet character. Excel summary sheet name.
#' @param report_sheet character. Excel detailed report sheet name.
#'
#' @return Invisible output file path.
#' @export
export_newdisc_meta_individual_report_excel <- function(
    loci_file = "inputs/newDisc.txt",
    meta_file = "inputs/GWAS/disc.tsv",
    individual_files = c(
      "inputs/GWAS/Chb/Chb1f.tsv", "inputs/GWAS/Chb/Chb2f.tsv",
      "inputs/GWAS/Chb/Chb3f.tsv", "inputs/GWAS/Chb/Chb4f.tsv",
      "inputs/GWAS/Chb/Chb1m.tsv", "inputs/GWAS/Chb/Chb2m.tsv",
      "inputs/GWAS/Chb/Chb3m.tsv", "inputs/GWAS/Chb/Chb4m.tsv"
    ),
    output_file = "outputs/newdisc_meta_individual_report.xlsx",
    subset_root_dir = "outputs/gwas_subsets_cache",
    p_threshold = 5e-5,
    annotation_file = "inputs/NCBI37.3.gene.loc",
    reported_file = "inputs/mehdi2024.txt",
    reported_max_distance_bp = 250000,
    individual_names = NULL,
    force_rebuild_subset = FALSE,
    summary_sheet = "Summary",
    report_sheet = "Locus_Report") {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package is required. Install with: install.packages('data.table')")
  }
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("openxlsx2 package is required. Install with: install.packages('openxlsx2')")
  }

  check_file <- function(path, label) {
    if (!nzchar(trimws(path))) stop(label, " is required.")
    if (!file.exists(path)) stop(label, " not found: ", path)
  }

  check_file(loci_file, "loci_file")
  check_file(meta_file, "meta_file")
  if (length(individual_files) == 0) stop("individual_files must contain at least one file.")
  for (fp in individual_files) check_file(fp, "individual_file")
  check_file(annotation_file, "annotation_file")
  check_file(reported_file, "reported_file")
  if (!is.finite(p_threshold) || p_threshold <= 0 || p_threshold >= 1) {
    stop("p_threshold must be between 0 and 1 (exclusive).")
  }
  if (!is.finite(reported_max_distance_bp) || reported_max_distance_bp < 0) {
    stop("reported_max_distance_bp must be >= 0.")
  }

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

  pick_first_matching_col <- function(df, aliases) {
    nms <- colnames(df)
    idx <- match(tolower(aliases), tolower(nms))
    idx <- idx[!is.na(idx)]
    if (length(idx) == 0) return(NULL)
    nms[[idx[[1]]]]
  }

  make_sanitized_names <- function(files, names_in = NULL) {
    if (is.null(names_in)) names_in <- tools::file_path_sans_ext(basename(files))
    if (length(names_in) != length(files)) {
      stop("individual_names length must match individual_files length.")
    }
    out <- gsub("[^A-Za-z0-9_]+", "_", names_in)
    out <- make.unique(out, sep = "_")
    out
  }

  distance_between_intervals <- function(s1, e1, s2, e2) {
    if (!is.finite(s1) || !is.finite(e1) || !is.finite(s2) || !is.finite(e2)) return(NA_real_)
    if (e1 < s2) return(s2 - e1)
    if (e2 < s1) return(s1 - e2)
    0
  }

  format_num <- function(x, digits = 6, scientific = FALSE) {
    if (!is.finite(x)) return(NA_character_)
    format(round(x, digits), scientific = scientific, trim = TRUE)
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

  # ---------------------------
  # Read loci
  # ---------------------------
  message("Reading loci file with fread: ", loci_file)
  loci_raw <- data.table::fread(loci_file, header = TRUE, data.table = FALSE, showProgress = FALSE)
  loci_raw <- as.data.frame(loci_raw, stringsAsFactors = FALSE)

  loci_col_id <- pick_first_matching_col(loci_raw, c("GenomicLocus", "genomiclocus", "locus", "locus_id"))
  loci_col_chr <- pick_first_matching_col(loci_raw, c("chr", "chrom", "chromosome", "CHR", "CHROM"))
  loci_col_start <- pick_first_matching_col(loci_raw, c("start", "START", "locusstart", "locus_start"))
  loci_col_end <- pick_first_matching_col(loci_raw, c("end", "END", "locusend", "locus_end"))

  if (is.null(loci_col_chr) || is.null(loci_col_start) || is.null(loci_col_end)) {
    stop("loci_file must include chr/start/end columns.")
  }

  loci <- data.frame(
    LOCUS_ID = if (is.null(loci_col_id)) seq_len(nrow(loci_raw)) else loci_raw[[loci_col_id]],
    CHR = normalize_chr(loci_raw[[loci_col_chr]]),
    START = suppressWarnings(as.numeric(loci_raw[[loci_col_start]])),
    END = suppressWarnings(as.numeric(loci_raw[[loci_col_end]])),
    stringsAsFactors = FALSE
  )

  loci <- loci[is.finite(loci$START) & is.finite(loci$END) & nzchar(loci$CHR), , drop = FALSE]
  if (nrow(loci) == 0) stop("No valid loci after parsing CHR/START/END.")

  swap_idx <- which(loci$START > loci$END)
  if (length(swap_idx) > 0) {
    tmp <- loci$START[swap_idx]
    loci$START[swap_idx] <- loci$END[swap_idx]
    loci$END[swap_idx] <- tmp
  }
  loci$LOCUS_ID <- as.character(loci$LOCUS_ID)

  # Distances to previous/next locus (within chromosome)
  loci$PREV_LOCUS_ID <- NA_character_
  loci$NEXT_LOCUS_ID <- NA_character_
  loci$DIST_TO_PREV_BP <- NA_real_
  loci$DIST_TO_NEXT_BP <- NA_real_

  for (chr in unique(loci$CHR)) {
    idx <- which(loci$CHR == chr)
    ord <- idx[order(loci$START[idx], loci$END[idx])]
    if (length(ord) == 1) next

    for (k in seq_along(ord)) {
      i <- ord[[k]]
      if (k > 1) {
        j_prev <- ord[[k - 1]]
        loci$PREV_LOCUS_ID[[i]] <- loci$LOCUS_ID[[j_prev]]
        loci$DIST_TO_PREV_BP[[i]] <- distance_between_intervals(loci$START[[i]], loci$END[[i]], loci$START[[j_prev]], loci$END[[j_prev]])
      }
      if (k < length(ord)) {
        j_next <- ord[[k + 1]]
        loci$NEXT_LOCUS_ID[[i]] <- loci$LOCUS_ID[[j_next]]
        loci$DIST_TO_NEXT_BP[[i]] <- distance_between_intervals(loci$START[[i]], loci$END[[i]], loci$START[[j_next]], loci$END[[j_next]])
      }
    }
  }

  # ---------------------------
  # Gene annotation
  # ---------------------------
  ann_raw <- data.table::fread(annotation_file, header = FALSE, data.table = FALSE, showProgress = FALSE)
  ann_raw <- as.data.frame(ann_raw, stringsAsFactors = FALSE)
  if (ncol(ann_raw) < 4) stop("annotation_file must contain at least 4 columns.")

  ann <- data.frame(
    GENE = as.character(if (ncol(ann_raw) >= 6) ann_raw[[6]] else ann_raw[[1]]),
    CHR = normalize_chr(ann_raw[[2]]),
    START = suppressWarnings(as.numeric(ann_raw[[3]])),
    END = suppressWarnings(as.numeric(ann_raw[[4]])),
    stringsAsFactors = FALSE
  )
  ann <- ann[is.finite(ann$START) & is.finite(ann$END) & nzchar(ann$GENE), , drop = FALSE]
  ann_by_chr <- split(ann, ann$CHR)

  assign_gene <- function(chr, s, e, ann_chr) {
    chr_ann <- ann_chr[[as.character(chr)]]
    if (is.null(chr_ann) || nrow(chr_ann) == 0) {
      return(list(gene_all = NA_character_, closest_gene = NA_character_, dist_bp = NA_real_))
    }

    locus_center <- (s + e) / 2

    overlap <- chr_ann$START <= e & chr_ann$END >= s
    if (any(overlap)) {
      ov <- chr_ann[overlap, , drop = FALSE]
      ov <- ov[nzchar(trimws(ov$GENE)), , drop = FALSE]
      if (nrow(ov) == 0) return(list(gene_all = NA_character_, closest_gene = NA_character_, dist_bp = 0))

      all_genes <- unique(as.character(ov$GENE))
      gene_centers <- (ov$START + ov$END) / 2
      k <- which.min(abs(gene_centers - locus_center))
      return(list(gene_all = paste(all_genes, collapse = ";"), closest_gene = as.character(ov$GENE[[k]]), dist_bp = 0))
    }

    dist_vec <- ifelse(
      e < chr_ann$START,
      chr_ann$START - e,
      ifelse(s > chr_ann$END, s - chr_ann$END, 0)
    )
    min_dist <- suppressWarnings(min(dist_vec, na.rm = TRUE))
    if (!is.finite(min_dist)) return(list(gene_all = NA_character_, closest_gene = NA_character_, dist_bp = NA_real_))

    cand <- chr_ann[is.finite(dist_vec) & dist_vec == min_dist, , drop = FALSE]
    cand <- cand[nzchar(trimws(cand$GENE)), , drop = FALSE]
    if (nrow(cand) == 0) return(list(gene_all = NA_character_, closest_gene = NA_character_, dist_bp = min_dist))

    all_genes <- unique(as.character(cand$GENE))
    gene_centers <- (cand$START + cand$END) / 2
    k <- which.min(abs(gene_centers - locus_center))
    list(gene_all = paste(all_genes, collapse = ";"), closest_gene = as.character(cand$GENE[[k]]), dist_bp = min_dist)
  }

  loci$GENE <- NA_character_
  loci$CLOSEST_GENE <- NA_character_
  loci$GENE_DISTANCE_BP <- NA_real_
  for (i in seq_len(nrow(loci))) {
    g <- assign_gene(loci$CHR[[i]], loci$START[[i]], loci$END[[i]], ann_by_chr)
    loci$GENE[[i]] <- g$gene_all
    loci$CLOSEST_GENE[[i]] <- g$closest_gene
    loci$GENE_DISTANCE_BP[[i]] <- g$dist_bp
  }

  # ---------------------------
  # Reported-locus matching
  # ---------------------------
  rep_raw <- data.table::fread(reported_file, header = TRUE, data.table = FALSE, showProgress = FALSE)
  rep_raw <- as.data.frame(rep_raw, stringsAsFactors = FALSE)

  rep_col_chr <- pick_first_matching_col(rep_raw, c("chr", "chrom", "chromosome"))
  rep_col_start <- pick_first_matching_col(rep_raw, c("locusstart", "start", "locus_start"))
  rep_col_end <- pick_first_matching_col(rep_raw, c("locusend", "end", "locus_end"))
  rep_col_name <- pick_first_matching_col(rep_raw, c("locusname", "reported_locus", "name"))
  rep_col_gene <- pick_first_matching_col(rep_raw, c("bestguessgene", "gene", "best_gene"))

  if (is.null(rep_col_chr) || is.null(rep_col_start) || is.null(rep_col_end)) {
    stop("reported_file must include chr/locusstart/locusend columns.")
  }

  reported <- data.frame(
    CHR = normalize_chr(rep_raw[[rep_col_chr]]),
    START = suppressWarnings(as.numeric(rep_raw[[rep_col_start]])),
    END = suppressWarnings(as.numeric(rep_raw[[rep_col_end]])),
    REPORTED_LOCUS_NAME = if (is.null(rep_col_name)) NA_character_ else as.character(rep_raw[[rep_col_name]]),
    REPORTED_GENE = if (is.null(rep_col_gene)) NA_character_ else as.character(rep_raw[[rep_col_gene]]),
    stringsAsFactors = FALSE
  )
  reported <- reported[is.finite(reported$START) & is.finite(reported$END) & nzchar(reported$CHR), , drop = FALSE]
  reported_by_chr <- split(reported, reported$CHR)

  match_reported <- function(chr, s, e, rep_chr, max_bp) {
    rr <- rep_chr[[as.character(chr)]]
    if (is.null(rr) || nrow(rr) == 0) {
      return(list(name = NA_character_, gene = NA_character_, matched = FALSE, diff_bp = NA_real_))
    }

    diffs <- pmin(
      abs(s - rr$START),
      abs(e - rr$END),
      abs(s - rr$END),
      abs(e - rr$START)
    )
    j <- which.min(diffs)
    if (length(j) == 0 || !is.finite(diffs[[j]]) || diffs[[j]] > max_bp) {
      return(list(name = NA_character_, gene = NA_character_, matched = FALSE, diff_bp = if (length(j) == 0) NA_real_ else diffs[[j]]))
    }

    list(
      name = rr$REPORTED_LOCUS_NAME[[j]],
      gene = rr$REPORTED_GENE[[j]],
      matched = TRUE,
      diff_bp = diffs[[j]]
    )
  }

  loci$REPORTED_LOCUS <- NA_character_
  loci$REPORTED_GENE <- NA_character_
  loci$REPORTED_MATCHED <- FALSE
  loci$REPORTED_BOUNDARY_DIFF_BP <- NA_real_

  for (i in seq_len(nrow(loci))) {
    m <- match_reported(loci$CHR[[i]], loci$START[[i]], loci$END[[i]], reported_by_chr, reported_max_distance_bp)
    loci$REPORTED_LOCUS[[i]] <- m$name
    loci$REPORTED_GENE[[i]] <- m$gene
    loci$REPORTED_MATCHED[[i]] <- isTRUE(m$matched)
    loci$REPORTED_BOUNDARY_DIFF_BP[[i]] <- m$diff_bp
  }

  # ---------------------------
  # Build subset-cache paths
  # ---------------------------
  all_files <- c(meta_file, individual_files)
  all_names <- c("meta", make_sanitized_names(individual_files, individual_names))

  cache_folder <- file.path(
    subset_root_dir,
    paste0(tools::file_path_sans_ext(basename(meta_file)), "__", tools::file_path_sans_ext(basename(loci_file)))
  )
  if (!dir.exists(cache_folder)) dir.create(cache_folder, recursive = TRUE, showWarnings = FALSE)

  cache_files <- setNames(file.path(cache_folder, paste0(all_names, ".subset.tsv")), all_names)

  subset_gwas_to_loci <- function(file_path, out_path) {
    message("Reading GWAS with fread: ", file_path)
    raw <- data.table::fread(file_path, header = TRUE, data.table = FALSE, showProgress = FALSE)
    raw <- as.data.frame(raw, stringsAsFactors = FALSE)

    c_chr <- pick_first_matching_col(raw, c("chr", "chrom", "chromosome", "CHR", "CHROM"))
    c_pos <- pick_first_matching_col(raw, c("pos", "position", "bp", "POS", "BP", "base_pair_location"))
    c_p <- pick_first_matching_col(raw, c("p", "pval", "p_value", "pvalue", "P", "PVAL", "P_VALUE", "p.value", "P-value"))
    c_beta <- pick_first_matching_col(raw, c("beta", "BETA", "effect"))

    if (is.null(c_chr) || is.null(c_pos) || is.null(c_p)) {
      stop("GWAS file missing CHR/POS/P columns after mapping: ", file_path)
    }

    raw$CHR_STD <- normalize_chr(raw[[c_chr]])
    raw$POS_STD <- suppressWarnings(as.numeric(raw[[c_pos]]))
    raw$P_STD <- suppressWarnings(as.numeric(raw[[c_p]]))
    raw$BETA_STD <- if (is.null(c_beta)) NA_real_ else suppressWarnings(as.numeric(raw[[c_beta]]))

    raw <- raw[is.finite(raw$POS_STD) & is.finite(raw$P_STD) & nzchar(raw$CHR_STD), , drop = FALSE]

    loci_idx <- split(loci, loci$CHR)
    keep <- logical(nrow(raw))
    for (chr in unique(raw$CHR_STD)) {
      idx <- which(raw$CHR_STD == chr)
      ranges <- loci_idx[[chr]]
      if (is.null(ranges) || nrow(ranges) == 0) next
      pos <- raw$POS_STD[idx]
      in_any <- rep(FALSE, length(pos))
      for (k in seq_len(nrow(ranges))) {
        in_any <- in_any | (pos >= ranges$START[[k]] & pos <= ranges$END[[k]])
      }
      keep[idx] <- in_any
    }
    sub <- raw[keep, , drop = FALSE]
    data.table::fwrite(sub, file = out_path, sep = "\t")
    sub
  }

  read_subset_or_build <- function(name_key, file_path) {
    out_path <- cache_files[[name_key]]
    if (file.exists(out_path) && !isTRUE(force_rebuild_subset)) {
      message("Found cached subset: ", out_path)
      return(as.data.frame(data.table::fread(out_path, header = TRUE, data.table = FALSE, showProgress = FALSE), stringsAsFactors = FALSE))
    }

    message("Building subset for: ", name_key)
    subset_gwas_to_loci(file_path, out_path)
  }

  subsets <- list()
  subsets[["meta"]] <- read_subset_or_build("meta", meta_file)
  for (i in seq_along(individual_files)) {
    nm <- all_names[[i + 1]]
    subsets[[nm]] <- read_subset_or_build(nm, individual_files[[i]])
  }

  # Prepare meta full-row output columns
  meta_subset <- subsets[["meta"]]
  meta_output_colnames <- make.unique(paste0("META_", colnames(meta_subset)), sep = "_")

  # ---------------------------
  # Build report table
  # ---------------------------
  report_rows <- vector("list", nrow(loci))

  for (i in seq_len(nrow(loci))) {
    chr <- loci$CHR[[i]]
    st <- loci$START[[i]]
    en <- loci$END[[i]]
    locus_id <- loci$LOCUS_ID[[i]]

    meta_in_locus <- meta_subset[
      meta_subset$CHR_STD == chr & meta_subset$POS_STD >= st & meta_subset$POS_STD <= en,
      , drop = FALSE
    ]

    n_lt <- if (nrow(meta_in_locus) == 0) 0L else sum(is.finite(meta_in_locus$P_STD) & meta_in_locus$P_STD < p_threshold)

    base <- data.frame(
      LOCUS_ID = locus_id,
      CHR = chr,
      START = as.integer(round(st)),
      END = as.integer(round(en)),
      GENE = loci$GENE[[i]],
      CLOSEST_GENE = loci$CLOSEST_GENE[[i]],
      GENE_DISTANCE_BP = if (is.finite(loci$GENE_DISTANCE_BP[[i]])) as.integer(round(loci$GENE_DISTANCE_BP[[i]])) else NA_integer_,
      DIST_TO_PREV_BP = if (is.finite(loci$DIST_TO_PREV_BP[[i]])) as.integer(round(loci$DIST_TO_PREV_BP[[i]])) else NA_integer_,
      DIST_TO_NEXT_BP = if (is.finite(loci$DIST_TO_NEXT_BP[[i]])) as.integer(round(loci$DIST_TO_NEXT_BP[[i]])) else NA_integer_,
      PREV_LOCUS_ID = loci$PREV_LOCUS_ID[[i]],
      NEXT_LOCUS_ID = loci$NEXT_LOCUS_ID[[i]],
      REPORTED_LOCUS = loci$REPORTED_LOCUS[[i]],
      REPORTED_GENE = loci$REPORTED_GENE[[i]],
      REPORTED_MATCHED = loci$REPORTED_MATCHED[[i]],
      REPORTED_BOUNDARY_DIFF_BP = if (is.finite(loci$REPORTED_BOUNDARY_DIFF_BP[[i]])) as.integer(round(loci$REPORTED_BOUNDARY_DIFF_BP[[i]])) else NA_integer_,
      N_META_SNPS = nrow(meta_in_locus),
      N_META_SNPS_P_LT_THRESHOLD = n_lt,
      stringsAsFactors = FALSE
    )

    if (nrow(meta_in_locus) == 0) {
      meta_vals <- as.list(rep(NA, length(meta_output_colnames)))
      names(meta_vals) <- meta_output_colnames
      row_df <- cbind(base, as.data.frame(meta_vals, stringsAsFactors = FALSE))
      row_df$LEAD_SNP_CLOSEST_GENE <- NA_character_
      row_df$LEAD_SNP_GENE_DISTANCE_BP <- NA_integer_

      for (nm in all_names[all_names != "meta"]) {
        row_df[[paste0(nm, "_P")]] <- NA_real_
        row_df[[paste0(nm, "_OR_FROM_BETA")]] <- NA_real_
      }

      report_rows[[i]] <- row_df
      next
    }

    lead_idx <- which.min(meta_in_locus$P_STD)
    lead_meta <- meta_in_locus[lead_idx, , drop = FALSE]
    lead_pos <- lead_meta$POS_STD[[1]]
    lead_gene <- assign_gene(chr, lead_pos, lead_pos, ann_by_chr)

    colnames(lead_meta) <- meta_output_colnames
    row_df <- cbind(base, lead_meta)
    row_df$LEAD_SNP_CLOSEST_GENE <- lead_gene$gene
    row_df$LEAD_SNP_GENE_DISTANCE_BP <- if (is.finite(lead_gene$dist_bp)) as.integer(round(lead_gene$dist_bp)) else NA_integer_

    for (nm in all_names[all_names != "meta"]) {
      ds <- subsets[[nm]]
      hit <- ds[ds$CHR_STD == chr & ds$POS_STD == lead_pos, , drop = FALSE]
      if (nrow(hit) == 0) {
        row_df[[paste0(nm, "_P")]] <- NA_real_
        row_df[[paste0(nm, "_OR_FROM_BETA")]] <- NA_real_
      } else {
        row_df[[paste0(nm, "_P")]] <- suppressWarnings(as.numeric(hit$P_STD[[1]]))
        beta_val <- suppressWarnings(as.numeric(hit$BETA_STD[[1]]))
        row_df[[paste0(nm, "_OR_FROM_BETA")]] <- if (is.finite(beta_val)) exp(beta_val) else NA_real_
      }
    }

    report_rows[[i]] <- row_df
  }

  report_df <- do.call(rbind, report_rows)

  # Format p-value columns in scientific notation for Excel readability
  pval_cols <- which(
    grepl("_P$|_P_STD$|^META_.*p($|_)", colnames(report_df), ignore.case = TRUE)
  )
  if (length(pval_cols) > 0) {
    for (j in pval_cols) {
      v <- suppressWarnings(as.numeric(report_df[[j]]))
      report_df[[j]] <- ifelse(is.finite(v), format(v, scientific = TRUE, digits = 4), NA_character_)
    }
  }

  summary_df <- data.frame(
    METRIC = c(
      "n_loci",
      "meta_file",
      "loci_file",
      "p_threshold",
      "total_meta_snps_in_loci",
      "total_meta_snps_p_lt_threshold",
      "reported_match_count",
      "subset_cache_folder"
    ),
    VALUE = c(
      nrow(loci),
      meta_file,
      loci_file,
      format_num(p_threshold, digits = 8, scientific = TRUE),
      sum(report_df$N_META_SNPS, na.rm = TRUE),
      sum(report_df$N_META_SNPS_P_LT_THRESHOLD, na.rm = TRUE),
      sum(report_df$REPORTED_MATCHED, na.rm = TRUE),
      cache_folder
    ),
    stringsAsFactors = FALSE
  )

  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  summary_sheet <- substr(ifelse(nzchar(trimws(summary_sheet)), trimws(summary_sheet), "Summary"), 1, 31)
  report_sheet <- substr(ifelse(nzchar(trimws(report_sheet)), trimws(report_sheet), "Locus_Report"), 1, 31)

  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet(summary_sheet)
  wb$add_data(sheet = summary_sheet, x = summary_df, start_row = 1, start_col = 1)
  wb$set_col_widths(sheet = summary_sheet, cols = 1:ncol(summary_df), widths = "auto")

  wb$add_worksheet(report_sheet)
  wb$add_data(sheet = report_sheet, x = report_df, start_row = 1, start_col = 1)

  # OR-based cell highlighting: >1 light red, <1 light blue
  red_fill <- openxlsx2::wb_color(hex = "FFFFC7CE")
  blue_fill <- openxlsx2::wb_color(hex = "FFDDEBF7")

  or_cols <- which(grepl("_OR_FROM_BETA$|^META_.*(^|_)OR($|_)", colnames(report_df), ignore.case = TRUE))
  if (length(or_cols) > 0) {
    for (j in or_cols) {
      v <- suppressWarnings(as.numeric(report_df[[j]]))
      rows_red <- which(is.finite(v) & v > 1)
      rows_blue <- which(is.finite(v) & v < 1)
      if (length(rows_red) > 0) add_fill_runs(wb, report_sheet, j, rows_red, red_fill)
      if (length(rows_blue) > 0) add_fill_runs(wb, report_sheet, j, rows_blue, blue_fill)
    }
  }

  wb$set_col_widths(sheet = report_sheet, cols = 1:ncol(report_df), widths = "auto")

  wb$save(output_file)

  message("Report saved: ", output_file)
  message("  Loci: ", nrow(loci))
  message("  Cache folder: ", cache_folder)

  invisible(output_file)
}
