#' Export target-loci GWAS follow-up workbook
#'
#' Reads one loci file and multiple GWAS files, then creates an Excel workbook
#' with two sheets:
#' 1) `Top_SNP_Per_Locus`: one row per locus, with the lowest-p SNP in each
#'    dataset formatted as `p|or (chr:pos:nea:ea)`.
#' 2) `Target_Followup`: for selected target datasets, picks each target
#'    dataset's lowest-p SNP per locus and reports the exact same SNP location
#'    in all datasets (same format).
#'
#' If a dataset has no SNP for a locus (or at the target location), the cell is
#' set to `-`.
#' GWAS rows are filtered immediately to keep only SNPs falling inside the
#' provided loci intervals, reducing memory usage and runtime.
#'
#' OR handling:
#' - Uses `OR` when available.
#' - Otherwise computes OR as `exp(BETA)` when `BETA` is available.
#'
#' Cell highlight:
#' - Cells with p-value `< 5e-8` are filled light red on both sheets.
#'
#' @param loci_file character. Path to loci table (must include CHR, START, END;
#'   `MERGED_LOCUS_ID` optional).
#' @param gwas_files character vector. Paths to GWAS files.
#' @param output_file character. Output Excel filepath.
#' @param dataset_names optional character vector of dataset names (same length
#'   as `gwas_files`). Defaults to filename stems.
#' @param target_datasets optional character vector subset of dataset names to
#'   use as anchors in the `Target_Followup` sheet. Default: all datasets.
#' @param validation_datasets optional character vector of validation datasets
#'   shared for all target datasets.
#'   - `NULL` (default): for each target dataset, use all other datasets.
#' @param validation_p_threshold numeric. Validation p-value threshold applied
#'   to validation dataset(s) at the same SNP position. At least one validation
#'   dataset must pass this threshold. Default: `5e-3`.
#'   Validation rule is:
#'   target SNP p-value < 5e-8 AND at least one validator dataset at the same
#'   SNP position has p-value < `validation_p_threshold` =>
#'   `validation = "validated"`;
#'   otherwise `validation = "invalid"`.
#'
#' @return Invisible output file path.
#' @export
export_loci_gwas_target_followup_excel <- function(loci_file,
                                                   gwas_files,
                                                   output_file,
                                                   dataset_names = NULL,
                                                   target_datasets = NULL,
                                                   validation_datasets = NULL,
                                                   validation_p_threshold = 5e-3) {
  if (missing(loci_file) || !nzchar(trimws(loci_file))) {
    stop("loci_file is required.")
  }
  if (!file.exists(loci_file)) {
    stop("loci_file not found: ", loci_file)
  }
  if (missing(gwas_files) || length(gwas_files) == 0) {
    stop("gwas_files is required and must contain at least one file.")
  }
  if (any(!file.exists(gwas_files))) {
    missing_files <- gwas_files[!file.exists(gwas_files)]
    stop("GWAS file(s) not found: ", paste(missing_files, collapse = ", "))
  }
  if (missing(output_file) || !nzchar(trimws(output_file))) {
    stop("output_file is required (e.g., 'outputs/loci_gwas_target_followup.xlsx').")
  }
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("openxlsx2 package is required. Install with: install.packages('openxlsx2')")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package is required. Install with: install.packages('data.table')")
  }
  validation_p_threshold <- suppressWarnings(as.numeric(validation_p_threshold))
  if (!is.finite(validation_p_threshold) || validation_p_threshold <= 0 || validation_p_threshold >= 1) {
    stop("validation_p_threshold must be a numeric value between 0 and 1 (exclusive).")
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

  sanitize_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

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

  format_p <- function(p) {
    if (!is.finite(p)) return("-")
    format(p, scientific = TRUE, digits = 4)
  }

  format_or <- function(or_value) {
    if (!is.finite(or_value)) return("-")
    format(round(or_value, 6), trim = TRUE, scientific = FALSE)
  }

  compute_or <- function(or_value, beta_value) {
    if (is.finite(or_value)) return(or_value)
    if (is.finite(beta_value)) return(exp(beta_value))
    NA_real_
  }

  format_snp <- function(chr, pos, p, or_value, nea, ea) {
    if (!is.finite(p) || !is.finite(pos) || !nzchar(chr)) return("-")
    nea_txt <- if (is.na(nea) || !nzchar(trimws(as.character(nea)))) "-" else as.character(nea)
    ea_txt <- if (is.na(ea) || !nzchar(trimws(as.character(ea)))) "-" else as.character(ea)
    paste0(
      format_p(p), "|", format_or(or_value),
      " (", chr, ":", as.integer(round(pos)), ":", nea_txt, ":", ea_txt, ")"
    )
  }

  standardize_gwas_local <- function(data) {
    mapping <- list(
      CHR = c("chr", "chrom", "chromosome", "CHR", "CHROM", "Chromosome"),
      POS = c("pos", "position", "bp", "POS", "BP", "Position", "base_pair_location"),
      P = c("p", "pval", "p_value", "pvalue", "P", "PVAL", "P_VALUE", "Pvalue", "p.value", "p.value.na", "P-value"),
      EFFECT_ALLELE = c("effectallele", "effect_allele", "allele1", "a1", "A1", "ea", "EA", "Allele1"),
      OTHER_ALLELE = c("otherallele", "other_allele", "allele2", "a2", "A2", "oa", "OA", "Allele2", "ref", "nea", "NEA", "non_effect_allele"),
      BETA = c("beta", "Beta", "BETA", "effect", "Effect"),
      OR = c("or", "OR", "odds_ratio", "OddsRatio")
    )

    original_cols <- colnames(data)
    new_cols <- original_cols
    mapped_indices <- logical(length(original_cols))

    for (i in seq_along(original_cols)) {
      col <- original_cols[i]
      for (std_name in names(mapping)) {
        if (tolower(col) %in% tolower(mapping[[std_name]])) {
          new_cols[i] <- std_name
          mapped_indices[i] <- TRUE
          break
        }
      }
    }

    data <- data[, mapped_indices, drop = FALSE]
    new_cols <- new_cols[mapped_indices]
    colnames(data) <- new_cols

    if (anyDuplicated(colnames(data))) {
      keep <- !duplicated(colnames(data))
      data <- data[, keep, drop = FALSE]
    }

    data
  }

  standardize_gwas <- function(data) {
    # Use local mapping on purpose: it explicitly supports dotted aliases
    # such as p.value / p.value.NA that appear in some cohort exports.
    standardize_gwas_local(data)
  }

  build_chr_interval_index <- function(loci_df) {
    chr_index <- list()
    chr_vals <- unique(as.character(loci_df$CHR))

    for (chr in chr_vals) {
      sub <- loci_df[loci_df$CHR == chr, c("START", "END"), drop = FALSE]
      sub <- sub[is.finite(sub$START) & is.finite(sub$END), , drop = FALSE]
      if (nrow(sub) == 0) next

      sub <- sub[order(sub$START, sub$END), , drop = FALSE]

      merged_start <- numeric(0)
      merged_end <- numeric(0)
      cur_start <- sub$START[[1]]
      cur_end <- sub$END[[1]]

      if (nrow(sub) > 1) {
        for (k in 2:nrow(sub)) {
          s <- sub$START[[k]]
          e <- sub$END[[k]]
          if (s <= cur_end) {
            cur_end <- max(cur_end, e)
          } else {
            merged_start <- c(merged_start, cur_start)
            merged_end <- c(merged_end, cur_end)
            cur_start <- s
            cur_end <- e
          }
        }
      }

      merged_start <- c(merged_start, cur_start)
      merged_end <- c(merged_end, cur_end)

      chr_index[[chr]] <- data.frame(
        START = merged_start,
        END = merged_end,
        stringsAsFactors = FALSE
      )
    }

    chr_index
  }

  keep_positions_in_loci <- function(chr, pos, chr_index) {
    ranges <- chr_index[[as.character(chr)]]
    keep <- rep(FALSE, length(pos))
    if (is.null(ranges) || nrow(ranges) == 0 || length(pos) == 0) return(keep)

    finite <- is.finite(pos)
    if (!any(finite)) return(keep)

    starts <- ranges$START
    ends <- ranges$END

    idx <- findInterval(pos[finite], starts)
    valid <- idx >= 1L & idx <= length(ends)

    keep_finite <- rep(FALSE, sum(finite))
    keep_finite[valid] <- pos[finite][valid] <= ends[idx[valid]]
    keep[finite] <- keep_finite
    keep
  }

  pick_top_snp <- function(chr_lookup, chr, start_pos, end_pos) {
    if (is.null(chr_lookup) || length(chr_lookup) == 0) return(NULL)
    bucket <- chr_lookup[[as.character(chr)]]
    if (is.null(bucket) || nrow(bucket) == 0) return(NULL)

    pos <- bucket$POS
    if (length(pos) == 0) return(NULL)

    left <- findInterval(start_pos, pos)
    start_idx <- if (left == 0L) 1L else if (pos[[left]] < start_pos) left + 1L else left
    end_idx <- findInterval(end_pos, pos)

    if (end_idx < 1L || start_idx > end_idx || start_idx > length(pos)) return(NULL)
    sub <- bucket[start_idx:end_idx, , drop = FALSE]
    if (nrow(sub) == 0) return(NULL)
    sub[which.min(sub$P), , drop = FALSE]
  }

  get_at_position <- function(chr_lookup, chr, pos) {
    if (is.null(chr_lookup) || length(chr_lookup) == 0 || !is.finite(pos)) return(NULL)
    bucket <- chr_lookup[[as.character(chr)]]
    if (is.null(bucket) || nrow(bucket) == 0) return(NULL)
    m <- match(pos, bucket$POS)
    if (is.na(m)) return(NULL)
    bucket[m, , drop = FALSE]
  }

  message("Reading loci file: ", loci_file)
  loci <- data.table::fread(loci_file, header = TRUE, sep = "\t", data.table = FALSE, showProgress = FALSE)
  loci <- as.data.frame(loci, stringsAsFactors = FALSE)

  req_loci <- c("CHR", "START", "END")
  miss_loci <- setdiff(req_loci, colnames(loci))
  if (length(miss_loci) > 0) {
    stop("loci_file is missing required columns: ", paste(miss_loci, collapse = ", "))
  }

  if (!"MERGED_LOCUS_ID" %in% colnames(loci)) {
    loci$MERGED_LOCUS_ID <- seq_len(nrow(loci))
  }

  loci$MERGED_LOCUS_ID <- as.character(loci$MERGED_LOCUS_ID)
  loci$CHR <- normalize_chr(loci$CHR)
  loci$START <- suppressWarnings(as.numeric(loci$START))
  loci$END <- suppressWarnings(as.numeric(loci$END))
  loci <- loci[is.finite(loci$START) & is.finite(loci$END) & loci$START <= loci$END, , drop = FALSE]
  if (nrow(loci) == 0) {
    stop("No valid loci rows after parsing CHR/START/END.")
  }
  loci_chr_index <- build_chr_interval_index(loci)

  if (is.null(dataset_names)) {
    dataset_names <- tools::file_path_sans_ext(basename(gwas_files))
  }
  if (length(dataset_names) != length(gwas_files)) {
    stop("dataset_names length must match gwas_files length.")
  }
  if (any(!nzchar(trimws(dataset_names)))) {
    stop("dataset_names cannot contain empty names.")
  }
  if (any(duplicated(dataset_names))) {
    stop("dataset_names must be unique.")
  }

  if (is.null(target_datasets)) {
    target_datasets <- dataset_names
  }
  if (!all(target_datasets %in% dataset_names)) {
    stop("All target_datasets must be present in dataset_names.")
  }

  # Build per-target validator mapping from one shared dataset list
  validation_map <- setNames(vector("list", length(target_datasets)), target_datasets)
  if (is.null(validation_datasets)) {
    for (target_ds in target_datasets) {
      validation_map[[target_ds]] <- setdiff(dataset_names, target_ds)
    }
  } else {
    if (!is.character(validation_datasets)) {
      stop("validation_datasets must be a character vector (shared for all target datasets) or NULL.")
    }
    shared_validators <- unique(validation_datasets[nzchar(trimws(validation_datasets))])
    bad <- setdiff(shared_validators, dataset_names)
    if (length(bad) > 0) {
      stop("validation_datasets contains unknown dataset(s): ", paste(bad, collapse = ", "))
    }
    for (target_ds in target_datasets) {
      validation_map[[target_ds]] <- setdiff(shared_validators, target_ds)
    }
  }

  message("Reading GWAS files one-by-one...")
  gwas_lookup <- list()
  for (i in seq_along(gwas_files)) {
    nm <- dataset_names[[i]]
    fp <- gwas_files[[i]]
    message(sprintf("  [%d/%d] %s", i, length(gwas_files), basename(fp)))

    dt_raw <- data.table::fread(fp, header = TRUE, sep = "\t", data.table = FALSE, showProgress = FALSE)
    dt_raw <- as.data.frame(dt_raw, stringsAsFactors = FALSE)
    raw_names_lc <- tolower(colnames(dt_raw))
    has_explicit_ea <- any(raw_names_lc %in% c("ea", "effect_allele", "effectallele"))
    has_explicit_nea <- any(raw_names_lc %in% c("nea", "other_allele", "otherallele", "non_effect_allele"))
    has_allele1 <- "allele1" %in% raw_names_lc
    has_allele2 <- "allele2" %in% raw_names_lc

    dt <- standardize_gwas(dt_raw)

    req_gwas <- c("CHR", "POS", "P")
    if (!all(req_gwas %in% colnames(dt))) {
      warning("Skipping dataset '", nm, "' due to missing required columns CHR/POS/P after standardization.")
      gwas_lookup[[nm]] <- list()
      next
    }

    # User convention: when only Allele1/Allele2 are present, Allele2 is effect allele.
    if (!has_explicit_ea && !has_explicit_nea && has_allele1 && has_allele2) {
      idx_a1 <- match("allele1", raw_names_lc)
      idx_a2 <- match("allele2", raw_names_lc)
      dt$EFFECT_ALLELE <- as.character(dt_raw[[idx_a2]])
      dt$OTHER_ALLELE <- as.character(dt_raw[[idx_a1]])
    }

    dt$CHR <- normalize_chr(dt$CHR)
    dt$POS <- suppressWarnings(as.numeric(dt$POS))
    dt$P <- suppressWarnings(as.numeric(dt$P))
    if (!"OTHER_ALLELE" %in% colnames(dt)) dt$OTHER_ALLELE <- NA_character_
    if (!"EFFECT_ALLELE" %in% colnames(dt)) dt$EFFECT_ALLELE <- NA_character_
    if (!"OR" %in% colnames(dt)) dt$OR <- NA_real_
    if (!"BETA" %in% colnames(dt)) dt$BETA <- NA_real_
    dt$OR <- suppressWarnings(as.numeric(dt$OR))
    dt$BETA <- suppressWarnings(as.numeric(dt$BETA))

    keep_loc <- logical(nrow(dt))
    if (nrow(dt) > 0) {
      chr_vals <- unique(dt$CHR)
      for (chr in chr_vals) {
        idx_chr <- which(dt$CHR == chr)
        keep_loc[idx_chr] <- keep_positions_in_loci(chr, dt$POS[idx_chr], loci_chr_index)
      }
    }

    keep <- is.finite(dt$POS) & is.finite(dt$P) & keep_loc
    dt <- dt[keep, c("CHR", "POS", "P", "OR", "BETA", "OTHER_ALLELE", "EFFECT_ALLELE"), drop = FALSE]

    if (nrow(dt) > 0) {
      ord <- order(dt$CHR, dt$POS, dt$P)
      dt <- dt[ord, , drop = FALSE]
      dt <- dt[!duplicated(paste(dt$CHR, dt$POS, sep = "::")), , drop = FALSE]
    }
    gwas_lookup[[nm]] <- if (nrow(dt) > 0) split(dt, dt$CHR) else list()
  }

  # Sheet 1: top SNP per locus in each dataset
  sheet1_rows <- vector("list", nrow(loci))
  sheet1_pvals <- list()
  for (ds in dataset_names) sheet1_pvals[[ds]] <- rep(NA_real_, nrow(loci))

  for (i in seq_len(nrow(loci))) {
    chr <- loci$CHR[[i]]
    st <- loci$START[[i]]
    en <- loci$END[[i]]

    row <- list(
      MERGED_LOCUS_ID = loci$MERGED_LOCUS_ID[[i]],
      CHR = chr,
      START = as.integer(round(st)),
      END = as.integer(round(en))
    )

    for (ds in dataset_names) {
      top <- pick_top_snp(gwas_lookup[[ds]], chr, st, en)
      if (is.null(top) || nrow(top) == 0) {
        row[[paste0("TOP_", sanitize_name(ds))]] <- "-"
        sheet1_pvals[[ds]][[i]] <- NA_real_
      } else {
        p <- top$P[[1]]
        or_val <- compute_or(top$OR[[1]], top$BETA[[1]])
        nea <- top$OTHER_ALLELE[[1]]
        ea <- top$EFFECT_ALLELE[[1]]
        row[[paste0("TOP_", sanitize_name(ds))]] <- format_snp(chr, top$POS[[1]], p, or_val, nea, ea)
        sheet1_pvals[[ds]][[i]] <- p
      }
    }
    sheet1_rows[[i]] <- as.data.frame(row, stringsAsFactors = FALSE)
  }
  sheet1 <- do.call(rbind, sheet1_rows)

  # Sheet 2: target follow-up (anchor dataset top SNP, then exact location in all datasets)
  sheet2_rows <- list()
  sheet2_target_p <- numeric(0)
  sheet2_follow_p <- list()
  for (ds in dataset_names) sheet2_follow_p[[ds]] <- numeric(0)

  row_counter <- 1L
  for (i in seq_len(nrow(loci))) {
    chr <- loci$CHR[[i]]
    st <- loci$START[[i]]
    en <- loci$END[[i]]

    for (target_ds in target_datasets) {
      top_target <- pick_top_snp(gwas_lookup[[target_ds]], chr, st, en)
      target_pos <- if (!is.null(top_target) && nrow(top_target) > 0) top_target$POS[[1]] else NA_real_
      validators <- validation_map[[target_ds]]
      pvals_at_target <- setNames(as.list(rep(NA_real_, length(dataset_names))), dataset_names)

      row <- list(
        MERGED_LOCUS_ID = loci$MERGED_LOCUS_ID[[i]],
        CHR = chr,
        START = as.integer(round(st)),
        END = as.integer(round(en)),
        TARGET_DATASET = target_ds,
        TARGET_POS = if (is.finite(target_pos)) as.integer(round(target_pos)) else NA_integer_
      )

      if (is.finite(target_pos)) {
        p_t <- top_target$P[[1]]
        or_t <- compute_or(top_target$OR[[1]], top_target$BETA[[1]])
        nea_t <- top_target$OTHER_ALLELE[[1]]
        ea_t <- top_target$EFFECT_ALLELE[[1]]
        row$TARGET_SNP <- format_snp(chr, target_pos, p_t, or_t, nea_t, ea_t)
        sheet2_target_p[row_counter] <- p_t
        pvals_at_target[[target_ds]] <- p_t
      } else {
        row$TARGET_SNP <- "-"
        sheet2_target_p[row_counter] <- NA_real_
      }

      for (ds in dataset_names) {
        val <- if (is.finite(target_pos)) get_at_position(gwas_lookup[[ds]], chr, target_pos) else NULL
        col_nm <- paste0("AT_TARGET_", sanitize_name(ds))

        if (is.null(val) || nrow(val) == 0) {
          row[[col_nm]] <- "-"
          sheet2_follow_p[[ds]][row_counter] <- NA_real_
          pvals_at_target[[ds]] <- NA_real_
        } else {
          p <- val$P[[1]]
          or_val <- compute_or(val$OR[[1]], val$BETA[[1]])
          nea <- val$OTHER_ALLELE[[1]]
          ea <- val$EFFECT_ALLELE[[1]]
          row[[col_nm]] <- format_snp(chr, target_pos, p, or_val, nea, ea)
          sheet2_follow_p[[ds]][row_counter] <- p
          pvals_at_target[[ds]] <- p
        }
      }

      target_p <- suppressWarnings(as.numeric(sheet2_target_p[row_counter]))
      validator_p <- suppressWarnings(as.numeric(unlist(pvals_at_target[validators], use.names = FALSE)))
      is_validated <- is.finite(target_p) &&
        target_p < 5e-8 &&
        any(is.finite(validator_p) & validator_p < validation_p_threshold)
      row$validation <- if (is_validated) "validated" else "invalid"

      sheet2_rows[[row_counter]] <- as.data.frame(row, stringsAsFactors = FALSE)
      row_counter <- row_counter + 1L
    }
  }
  sheet2 <- do.call(rbind, sheet2_rows)

  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("Top_SNP_Per_Locus")
  wb$add_data(sheet = "Top_SNP_Per_Locus", x = sheet1, start_row = 1, start_col = 1)

  # Light-red fill for significant p-values (< 5e-8)
  red_fill <- openxlsx2::wb_color(hex = "FFFFC7CE")

  for (ds in dataset_names) {
    col_name <- paste0("TOP_", sanitize_name(ds))
    col_idx <- which(colnames(sheet1) == col_name)
    if (length(col_idx) == 1) {
      rows_hit <- which(is.finite(sheet1_pvals[[ds]]) & sheet1_pvals[[ds]] < 5e-8)
      if (length(rows_hit) > 0) add_fill_runs(wb, "Top_SNP_Per_Locus", col_idx, rows_hit, red_fill)
    }
  }

  # One follow-up sheet per target dataset
  sheet_name_counts <- integer(0)
  build_target_sheet_name <- function(target_name) {
    base <- paste0("Followup_", sanitize_name(target_name))
    base <- substr(base, 1, 31)

    prev_n <- if (base %in% names(sheet_name_counts)) sheet_name_counts[[base]] else 0L
    n <- prev_n + 1L
    sheet_name_counts[[base]] <<- n
    if (n == 1L) return(base)

    suffix <- paste0("_", n)
    max_base <- 31 - nchar(suffix)
    paste0(substr(base, 1, max_base), suffix)
  }

  for (target_ds in target_datasets) {
    target_rows <- which(sheet2$TARGET_DATASET == target_ds)
    if (length(target_rows) == 0) next

    target_sheet <- sheet2[target_rows, , drop = FALSE]
    sheet_name <- build_target_sheet_name(target_ds)
    wb$add_worksheet(sheet_name)
    wb$add_data(sheet = sheet_name, x = target_sheet, start_row = 1, start_col = 1)

    # Highlight target SNP column for this target dataset
    target_col_idx <- which(colnames(target_sheet) == "TARGET_SNP")
    if (length(target_col_idx) == 1) {
      p_target <- sheet2_target_p[target_rows]
      rows_hit <- which(is.finite(p_target) & p_target < 5e-8)
      if (length(rows_hit) > 0) add_fill_runs(wb, sheet_name, target_col_idx, rows_hit, red_fill)
    }

    # Highlight each dataset's AT_TARGET column in this target sheet
    for (ds in dataset_names) {
      col_name <- paste0("AT_TARGET_", sanitize_name(ds))
      col_idx <- which(colnames(target_sheet) == col_name)
      if (length(col_idx) == 1) {
        pvec <- sheet2_follow_p[[ds]][target_rows]
        rows_hit <- which(is.finite(pvec) & pvec < 5e-8)
        if (length(rows_hit) > 0) add_fill_runs(wb, sheet_name, col_idx, rows_hit, red_fill)
      }
    }

    wb$set_col_widths(sheet = sheet_name, cols = 1:ncol(target_sheet), widths = "auto")
  }

  wb$set_col_widths(sheet = "Top_SNP_Per_Locus", cols = 1:ncol(sheet1), widths = "auto")
  wb$save(output_file)

  message(sprintf("Workbook saved: %s", output_file))
  message(sprintf("  Loci rows: %d", nrow(sheet1)))
  message(sprintf("  Follow-up rows (total): %d", nrow(sheet2)))
  message(sprintf("  Target sheets: %d", length(target_datasets)))
  message(sprintf("  GWAS datasets: %d", length(dataset_names)))

  invisible(output_file)
}
