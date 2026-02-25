#' Export lead SNP per FUMA locus to Excel
#'
#' Reads one FUMA loci file and one GWAS file, subsets GWAS rows to each FUMA
#' locus interval, and extracts the lead SNP (lowest p-value) per locus.
#' Also annotates each locus with the nearest neighboring locus on the same
#' chromosome and the inter-locus distance in base pairs.
#' The output also includes all original GWAS columns for the selected lead SNP
#' in each locus, prefixed as `GWAS_<column_name>`.
#'
#' Expected FUMA columns include: `GenomicLocus`, `chr`, `start`, `end`.
#' `GenomicLocus` is optional (row index is used when missing).
#'
#' Expected GWAS columns include chromosome, position, and p-value; common
#' aliases are automatically recognized (e.g., `chrom`, `pos`, `p`).
#'
#' @param fuma_file character. Path to FUMA loci file.
#' @param gwas_file character. Path to GWAS file.
#' @param output_file character. Output Excel filepath.
#' @param sheet_name character. Excel sheet name. Default: "Lead_SNPs".
#'
#' @return Invisible output file path.
#' @export
export_fuma_lead_snps_excel <- function(fuma_file,
                                        gwas_file,
                                        output_file,
                                        sheet_name = "Lead_SNPs") {
  if (missing(fuma_file) || !nzchar(trimws(fuma_file))) {
    stop("fuma_file is required.")
  }
  if (!file.exists(fuma_file)) {
    stop("fuma_file not found: ", fuma_file)
  }
  if (missing(gwas_file) || !nzchar(trimws(gwas_file))) {
    stop("gwas_file is required.")
  }
  if (!file.exists(gwas_file)) {
    stop("gwas_file not found: ", gwas_file)
  }
  if (missing(output_file) || !nzchar(trimws(output_file))) {
    stop("output_file is required (e.g., 'outputs/fuma_lead_snps.xlsx').")
  }
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("data.table package is required. Install with: install.packages('data.table')")
  }
  if (!requireNamespace("openxlsx2", quietly = TRUE)) {
    stop("openxlsx2 package is required. Install with: install.packages('openxlsx2')")
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

  distance_between_intervals <- function(s1, e1, s2, e2) {
    if (!is.finite(s1) || !is.finite(e1) || !is.finite(s2) || !is.finite(e2)) return(NA_real_)
    if (e1 < s2) return(s2 - e1)
    if (e2 < s1) return(s1 - e2)
    0
  }

  message("Reading FUMA loci file: ", fuma_file)
  fuma <- data.table::fread(fuma_file, header = TRUE, sep = "\t", data.table = FALSE, showProgress = FALSE)
  fuma <- as.data.frame(fuma, stringsAsFactors = FALSE)

  col_genomic_locus <- pick_first_matching_col(fuma, c("GenomicLocus", "genomiclocus", "locus", "locus_id"))
  col_chr <- pick_first_matching_col(fuma, c("chr", "chrom", "chromosome", "CHR", "CHROM"))
  col_start <- pick_first_matching_col(fuma, c("start", "START", "locus_start"))
  col_end <- pick_first_matching_col(fuma, c("end", "END", "locus_end"))

  if (is.null(col_chr) || is.null(col_start) || is.null(col_end)) {
    stop("FUMA file must include chromosome/start/end columns (e.g., chr, start, end).")
  }

  loci <- data.frame(
    GenomicLocus = if (is.null(col_genomic_locus)) seq_len(nrow(fuma)) else fuma[[col_genomic_locus]],
    CHR = normalize_chr(fuma[[col_chr]]),
    START = suppressWarnings(as.numeric(fuma[[col_start]])),
    END = suppressWarnings(as.numeric(fuma[[col_end]])),
    stringsAsFactors = FALSE
  )

  loci <- loci[is.finite(loci$START) & is.finite(loci$END) & nzchar(loci$CHR), , drop = FALSE]
  if (nrow(loci) == 0) {
    stop("No valid loci rows after parsing FUMA CHR/START/END.")
  }

  loci$GenomicLocus <- as.character(loci$GenomicLocus)
  loci$START <- as.numeric(loci$START)
  loci$END <- as.numeric(loci$END)

  swap_idx <- which(loci$START > loci$END)
  if (length(swap_idx) > 0) {
    tmp <- loci$START[swap_idx]
    loci$START[swap_idx] <- loci$END[swap_idx]
    loci$END[swap_idx] <- tmp
  }

  message("Reading GWAS file: ", gwas_file)
  gwas <- data.table::fread(gwas_file, header = TRUE, data.table = FALSE, showProgress = FALSE)
  gwas <- as.data.frame(gwas, stringsAsFactors = FALSE)
  extra_out_colnames <- make.unique(paste0("GWAS_", colnames(gwas)), sep = "_")

  col_g_chr <- pick_first_matching_col(gwas, c("chr", "chrom", "chromosome", "CHR", "CHROM", "chrom"))
  col_g_pos <- pick_first_matching_col(gwas, c("pos", "position", "bp", "POS", "BP", "base_pair_location"))
  col_g_p <- pick_first_matching_col(gwas, c("p", "pval", "p_value", "pvalue", "P", "PVAL", "P_VALUE", "p.value", "p.value.na", "P-value"))
  col_g_snp <- pick_first_matching_col(gwas, c("snp", "snpid", "rsid", "markername", "variant", "variant_id"))
  col_g_ea <- pick_first_matching_col(gwas, c("ea", "effect_allele", "effectallele", "a1", "allele1"))
  col_g_nea <- pick_first_matching_col(gwas, c("nea", "other_allele", "otherallele", "a2", "allele2", "ref", "non_effect_allele"))
  col_g_beta <- pick_first_matching_col(gwas, c("beta", "BETA", "effect"))
  col_g_se <- pick_first_matching_col(gwas, c("se", "SE", "stderr", "standard_error"))
  col_g_or <- pick_first_matching_col(gwas, c("or", "OR", "odds_ratio"))

  if (is.null(col_g_chr) || is.null(col_g_pos) || is.null(col_g_p)) {
    stop("GWAS file must include chromosome/position/p-value columns.")
  }

  gw <- data.frame(
    RAW_ROW = seq_len(nrow(gwas)),
    CHR = normalize_chr(gwas[[col_g_chr]]),
    POS = suppressWarnings(as.numeric(gwas[[col_g_pos]])),
    P = suppressWarnings(as.numeric(gwas[[col_g_p]])),
    SNP = if (is.null(col_g_snp)) NA_character_ else as.character(gwas[[col_g_snp]]),
    EA = if (is.null(col_g_ea)) NA_character_ else as.character(gwas[[col_g_ea]]),
    NEA = if (is.null(col_g_nea)) NA_character_ else as.character(gwas[[col_g_nea]]),
    BETA = if (is.null(col_g_beta)) NA_real_ else suppressWarnings(as.numeric(gwas[[col_g_beta]])),
    SE = if (is.null(col_g_se)) NA_real_ else suppressWarnings(as.numeric(gwas[[col_g_se]])),
    OR = if (is.null(col_g_or)) NA_real_ else suppressWarnings(as.numeric(gwas[[col_g_or]])),
    stringsAsFactors = FALSE
  )

  gw <- gw[is.finite(gw$POS) & is.finite(gw$P) & nzchar(gw$CHR), , drop = FALSE]
  if (nrow(gw) == 0) {
    stop("No valid GWAS rows after parsing CHR/POS/P.")
  }

  # Keep only SNPs in any FUMA locus interval
  locus_index <- split(loci, loci$CHR)
  keep <- logical(nrow(gw))

  for (chr in unique(gw$CHR)) {
    idx <- which(gw$CHR == chr)
    ranges <- locus_index[[chr]]
    if (is.null(ranges) || nrow(ranges) == 0) next

    pos <- gw$POS[idx]
    in_any <- rep(FALSE, length(pos))
    for (k in seq_len(nrow(ranges))) {
      in_any <- in_any | (pos >= ranges$START[[k]] & pos <= ranges$END[[k]])
    }
    keep[idx] <- in_any
  }

  gw <- gw[keep, , drop = FALSE]

  # Prepare nearest-locus annotations (same chromosome)
  nearest_locus_id <- rep(NA_character_, nrow(loci))
  nearest_locus_dist <- rep(NA_real_, nrow(loci))

  for (chr in unique(loci$CHR)) {
    idx_chr <- which(loci$CHR == chr)
    if (length(idx_chr) == 1) {
      nearest_locus_id[idx_chr] <- NA_character_
      nearest_locus_dist[idx_chr] <- NA_real_
      next
    }

    for (i_local in seq_along(idx_chr)) {
      i <- idx_chr[[i_local]]
      others <- setdiff(idx_chr, i)
      dvec <- vapply(others, function(j) {
        distance_between_intervals(loci$START[[i]], loci$END[[i]], loci$START[[j]], loci$END[[j]])
      }, numeric(1))

      if (all(!is.finite(dvec))) next
      j_best <- others[[which.min(dvec)]]
      nearest_locus_id[[i]] <- loci$GenomicLocus[[j_best]]
      nearest_locus_dist[[i]] <- dvec[[which.min(dvec)]]
    }
  }

  # Lead SNP per locus
  rows_out <- vector("list", nrow(loci))
  for (i in seq_len(nrow(loci))) {
    chr <- loci$CHR[[i]]
    st <- loci$START[[i]]
    en <- loci$END[[i]]

    sub <- gw[gw$CHR == chr & gw$POS >= st & gw$POS <= en, , drop = FALSE]

    if (nrow(sub) == 0) {
      base_row <- data.frame(
        GenomicLocus = loci$GenomicLocus[[i]],
        CHR = chr,
        START = as.integer(round(st)),
        END = as.integer(round(en)),
        LEAD_POS = NA_integer_,
        LEAD_P = NA_real_,
        LEAD_SNP = NA_character_,
        EA = NA_character_,
        NEA = NA_character_,
        BETA = NA_real_,
        SE = NA_real_,
        OR = NA_real_,
        N_SNPS_IN_LOCUS = 0L,
        NEAREST_LOCUS = nearest_locus_id[[i]],
        DIST_TO_NEAREST_LOCUS_BP = if (is.finite(nearest_locus_dist[[i]])) as.integer(round(nearest_locus_dist[[i]])) else NA_integer_,
        stringsAsFactors = FALSE
      )
      extra_values <- as.list(rep(NA, length(extra_out_colnames)))
      names(extra_values) <- extra_out_colnames
      rows_out[[i]] <- cbind(base_row, as.data.frame(extra_values, stringsAsFactors = FALSE))
      next
    }

    lead_idx <- which.min(sub$P)
    lead <- sub[lead_idx, , drop = FALSE]
    lead_raw <- gwas[lead$RAW_ROW[[1]], , drop = FALSE]
    colnames(lead_raw) <- extra_out_colnames

    base_row <- data.frame(
      GenomicLocus = loci$GenomicLocus[[i]],
      CHR = chr,
      START = as.integer(round(st)),
      END = as.integer(round(en)),
      LEAD_POS = as.integer(round(lead$POS[[1]])),
      LEAD_P = lead$P[[1]],
      LEAD_SNP = lead$SNP[[1]],
      EA = lead$EA[[1]],
      NEA = lead$NEA[[1]],
      BETA = lead$BETA[[1]],
      SE = lead$SE[[1]],
      OR = lead$OR[[1]],
      N_SNPS_IN_LOCUS = nrow(sub),
      NEAREST_LOCUS = nearest_locus_id[[i]],
      DIST_TO_NEAREST_LOCUS_BP = if (is.finite(nearest_locus_dist[[i]])) as.integer(round(nearest_locus_dist[[i]])) else NA_integer_,
      stringsAsFactors = FALSE
    )

    rows_out[[i]] <- cbind(base_row, lead_raw)
  }

  out <- do.call(rbind, rows_out)

  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  sheet_name <- substr(trimws(as.character(sheet_name)), 1, 31)
  if (!nzchar(sheet_name)) sheet_name <- "Lead_SNPs"

  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet(sheet_name)
  wb$add_data(sheet = sheet_name, x = out, start_row = 1, start_col = 1)
  wb$set_col_widths(sheet = sheet_name, cols = 1:ncol(out), widths = "auto")
  wb$save(output_file)

  message("Lead SNP workbook saved: ", output_file)
  message("  Loci rows: ", nrow(out))
  message("  GWAS rows retained in loci: ", nrow(gw))

  invisible(output_file)
}
