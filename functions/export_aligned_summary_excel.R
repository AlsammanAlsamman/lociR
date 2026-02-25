#' Export aligned loci summary workbook
#'
#' Creates an Excel workbook with:
#' 1) A per-dataset summary sheet (total/shared/unique merged loci)
#' 2) A per-locus table including reported locus (or "Novel"), top SNP per GWAS
#'    dataset formatted as "pvalue (rsid)", and genes containing those top SNPs.
#'
#' @param aligned_loci data.frame returned by `align_with_reported()` (class `AlignedLoci`)
#' @param output_file character. Output Excel path (required)
#' @param gwas_list optional named list of `GWAS` objects returned by `read_gwas()`
#' @param annotation_file character. Gene annotation file path.
#'   Default: "inputs/NCBI37.3.gene.loc"
#'
#' @return Invisible output file path.
#' @export
export_aligned_summary_excel <- function(aligned_loci,
                                         output_file,
                                         gwas_list = NULL,
                                         annotation_file = "inputs/NCBI37.3.gene.loc") {
  if (missing(output_file) || !nzchar(trimws(output_file))) {
    stop("output_file is required (e.g., 'outputs/aligned_loci_summary.xlsx').")
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

    # Expected file layout here is typically:
    # gene, chr, start, end, strand, gene
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

  parse_variant_alleles <- function(row_chr, row_pos, candidate_ids) {
    if (length(candidate_ids) == 0) return(list(nea = NA_character_, ea = NA_character_))

    chr_txt <- if (is.na(row_chr)) NA_character_ else as.character(row_chr)
    pos_txt <- if (!is.finite(row_pos)) NA_character_ else as.character(as.integer(round(row_pos)))

    for (idv in candidate_ids) {
      if (is.na(idv) || !nzchar(trimws(as.character(idv)))) next
      s <- as.character(idv)
      # Expected common formats:
      #  chr:pos:ref:alt
      #  chr_pos_ref_alt
      #  chr-pos-ref-alt
      parts <- unlist(strsplit(s, "[:_\\-]"))
      if (length(parts) < 4) next

      id_chr <- normalize_chr(parts[[1]])
      id_pos <- suppressWarnings(as.integer(parts[[2]]))
      id_ref <- parts[[3]]
      id_alt <- parts[[4]]

      if (!nzchar(id_ref) || !nzchar(id_alt)) next

      chr_ok <- is.na(chr_txt) || identical(as.character(id_chr), as.character(chr_txt))
      pos_ok <- is.na(pos_txt) || (!is.na(id_pos) && identical(as.character(id_pos), as.character(pos_txt)))

      if (chr_ok && pos_ok) {
        return(list(nea = id_ref, ea = id_alt))
      }
    }

    list(nea = NA_character_, ea = NA_character_)
  }

  format_top_snp <- function(p, chr, pos, nea, ea) {
    if (!is.finite(p)) return(NA_character_)
    chr_txt <- if (is.na(chr) || !nzchar(trimws(as.character(chr)))) "NA" else as.character(chr)
    pos_txt <- if (!is.finite(pos)) "NA" else as.character(as.integer(round(pos)))
    nea_txt <- if (is.na(nea) || !nzchar(trimws(as.character(nea)))) "?" else as.character(nea)
    ea_txt <- if (is.na(ea) || !nzchar(trimws(as.character(ea)))) "?" else as.character(ea)
    snp_txt <- paste(chr_txt, pos_txt, nea_txt, ea_txt, sep = ":")
    paste0(format(p, scientific = TRUE, digits = 3), " (", snp_txt, ")")
  }

  pick_top_snp <- function(gwas_df, chr, start_pos, end_pos) {
    if (is.null(gwas_df) || nrow(gwas_df) == 0) {
      return(list(text = NA_character_, chr = chr, pos = NA_real_))
    }

    keep <- gwas_df$CHR == chr &
      is.finite(gwas_df$POS) &
      gwas_df$POS >= start_pos &
      gwas_df$POS <= end_pos &
      is.finite(gwas_df$P)

    sub <- gwas_df[keep, , drop = FALSE]
    if (nrow(sub) == 0) {
      return(list(text = NA_character_, chr = chr, pos = NA_real_))
    }

    best_idx <- which.min(sub$P)
    best <- sub[best_idx, , drop = FALSE]

    nea <- if ("OTHER_ALLELE" %in% colnames(best)) best$OTHER_ALLELE[[1]] else NA_character_
    ea <- if ("EFFECT_ALLELE" %in% colnames(best)) best$EFFECT_ALLELE[[1]] else NA_character_

    if ((is.na(nea) || !nzchar(trimws(as.character(nea)))) ||
        (is.na(ea) || !nzchar(trimws(as.character(ea))))) {
      id_candidates <- character(0)
      for (id_col in c("RSID", "SNP", "MARKER", "VARIANT_ID", "ID")) {
        if (id_col %in% colnames(best)) {
          id_candidates <- c(id_candidates, as.character(best[[id_col]][[1]]))
        }
      }
      parsed <- parse_variant_alleles(best$CHR[[1]], best$POS[[1]], id_candidates)
      if (is.na(nea) || !nzchar(trimws(as.character(nea)))) nea <- parsed$nea
      if (is.na(ea) || !nzchar(trimws(as.character(ea)))) ea <- parsed$ea
    }
    text <- format_top_snp(
      p = best$P[[1]],
      chr = best$CHR[[1]],
      pos = best$POS[[1]],
      nea = nea,
      ea = ea
    )

    list(
      text = text,
      chr = as.character(best$CHR[[1]]),
      pos = as.numeric(best$POS[[1]])
    )
  }

  # Required columns check
  required_cols <- c("MERGED_LOCUS_ID", "CHR", "MERGED_START", "MERGED_END", "SOURCES")
  missing_cols <- setdiff(required_cols, colnames(aligned_loci))
  if (length(missing_cols) > 0) {
    stop("aligned_loci is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  df <- aligned_loci
  df$MERGED_LOCUS_ID <- as.character(df$MERGED_LOCUS_ID)
  df$CHR <- normalize_chr(df$CHR)
  df$MERGED_START <- suppressWarnings(as.numeric(df$MERGED_START))
  df$MERGED_END <- suppressWarnings(as.numeric(df$MERGED_END))

  if ("REPORTED_LOCUS" %in% colnames(df)) {
    df$REPORTED_LOCUS <- as.character(df$REPORTED_LOCUS)
  } else {
    df$REPORTED_LOCUS <- NA_character_
  }

  merged_keys <- unique(df$MERGED_LOCUS_ID)

  merged_table <- do.call(rbind, lapply(merged_keys, function(mid) {
    sub <- df[df$MERGED_LOCUS_ID == mid, , drop = FALSE]

    reported_vals <- unique(trimws(sub$REPORTED_LOCUS))
    reported_vals <- reported_vals[!is.na(reported_vals) & nzchar(reported_vals)]
    reported_label <- if (length(reported_vals) > 0) paste(reported_vals, collapse = "; ") else "Novel"

    data.frame(
      MERGED_LOCUS_ID = mid,
      CHR = sub$CHR[[1]],
      START = sub$MERGED_START[[1]],
      END = sub$MERGED_END[[1]],
      REPORTED_LOCUS = reported_label,
      SOURCES = as.character(sub$SOURCES[[1]]),
      stringsAsFactors = FALSE
    )
  }))

  # Build GWAS lookup
  gwas_lookup <- list()
  dataset_names <- character(0)

  if (!is.null(gwas_list)) {
    if (inherits(gwas_list, "GWAS")) {
      gwas_list <- list(gwas_list)
      names(gwas_list) <- gwas_list[[1]]$name %||% "GWAS"
    }

    if (is.list(gwas_list) && length(gwas_list) > 0) {
      dataset_names <- names(gwas_list)

      for (nm in dataset_names) {
        obj <- gwas_list[[nm]]
        dat <- if (is.list(obj) && !is.null(obj$data)) obj$data else NULL
        if (!is.null(dat) && is.data.frame(dat)) {
          dat <- as.data.frame(dat, stringsAsFactors = FALSE)
          if (all(c("CHR", "POS", "P") %in% colnames(dat))) {
            dat$CHR <- normalize_chr(dat$CHR)
            dat$POS <- suppressWarnings(as.numeric(dat$POS))
            dat$P <- suppressWarnings(as.numeric(dat$P))
            gwas_lookup[[nm]] <- dat
          }
        }
      }
    }
  }

  if (length(dataset_names) == 0) {
    source_pool <- unique(unlist(lapply(merged_table$SOURCES, split_sources)))
    dataset_names <- source_pool[nzchar(source_pool)]
  }

  # Dataset summary: total/shared/unique merged loci by source membership
  source_lists <- lapply(merged_table$SOURCES, split_sources)

  dataset_summary <- do.call(rbind, lapply(dataset_names, function(ds) {
    present <- vapply(source_lists, has_source, logical(1), dataset_name = ds)
    n_sources_each <- vapply(source_lists, length, numeric(1))

    total_loci <- sum(present)
    shared_loci <- sum(present & n_sources_each > 1)
    unique_loci <- sum(present & n_sources_each == 1)

    data.frame(
      DATASET = ds,
      TOTAL_LOCI = total_loci,
      SHARED_LOCI = shared_loci,
      UNIQUE_LOCI = unique_loci,
      stringsAsFactors = FALSE
    )
  }))

  sanitize_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

  # Replace combined SOURCES with one column per dataset.
  # Each dataset column gets a sequential locus index where present (1, 2, 3, ...);
  # absent loci stay empty (NA).
  source_index_cols <- character(0)
  for (ds in dataset_names) {
    col_name <- paste0("SOURCE_IDX_", sanitize_name(ds))
    present <- vapply(source_lists, has_source, logical(1), dataset_name = ds)

    idx_vec <- rep(NA_integer_, nrow(merged_table))
    n_present <- sum(present)
    if (n_present > 0) {
      idx_vec[present] <- seq_len(n_present)
    }

    merged_table[[col_name]] <- idx_vec
    source_index_cols <- c(source_index_cols, col_name)
  }

  # Add top SNP columns
  top_snp_positions <- vector("list", length = nrow(merged_table))
  for (i in seq_len(nrow(merged_table))) {
    top_snp_positions[[i]] <- data.frame(
      CHR = character(0),
      POS = numeric(0),
      stringsAsFactors = FALSE
    )
  }

  for (ds in dataset_names) {
    col_name <- paste0("TOP_SNP_", sanitize_name(ds))
    merged_table[[col_name]] <- NA_character_

    ds_data <- NULL
    if (length(gwas_lookup) > 0) {
      exact_idx <- which(names(gwas_lookup) == ds)
      if (length(exact_idx) == 1) {
        ds_data <- gwas_lookup[[exact_idx]]
      } else {
        ci_idx <- which(tolower(names(gwas_lookup)) == tolower(ds))
        if (length(ci_idx) >= 1) ds_data <- gwas_lookup[[ci_idx[1]]]
      }
    }

    for (i in seq_len(nrow(merged_table))) {
      top <- pick_top_snp(
        gwas_df = ds_data,
        chr = merged_table$CHR[[i]],
        start_pos = merged_table$START[[i]],
        end_pos = merged_table$END[[i]]
      )

      merged_table[[col_name]][[i]] <- top$text

      if (is.finite(top$pos)) {
        top_snp_positions[[i]] <- rbind(
          top_snp_positions[[i]],
          data.frame(CHR = top$chr, POS = top$pos, stringsAsFactors = FALSE)
        )
      }
    }
  }

  # Gene mapping for top SNPs
  ann <- read_gene_annotation(annotation_file)
  ann_by_chr <- split(ann, ann$CHR)

  find_genes_for_snp <- function(chr, pos) {
    if (!is.finite(pos) || !nzchar(chr)) return(character(0))
    chr_ann <- ann_by_chr[[as.character(chr)]]
    if (is.null(chr_ann) || nrow(chr_ann) == 0) return(character(0))

    hit <- chr_ann$START <= pos & chr_ann$END >= pos
    inside_genes <- unique(chr_ann$GENE[hit])
    inside_genes <- inside_genes[nzchar(trimws(inside_genes))]
    if (length(inside_genes) > 0) {
      return(inside_genes)
    }

    # If SNP is not inside a gene, return nearest gene(s)
    dist_to_gene <- ifelse(
      pos < chr_ann$START,
      chr_ann$START - pos,
      ifelse(pos > chr_ann$END, pos - chr_ann$END, 0)
    )
    min_dist <- suppressWarnings(min(dist_to_gene, na.rm = TRUE))
    if (!is.finite(min_dist)) return(character(0))

    nearest_genes <- unique(chr_ann$GENE[is.finite(dist_to_gene) & dist_to_gene == min_dist])
    nearest_genes[nzchar(trimws(nearest_genes))]
  }

  genes_col <- character(nrow(merged_table))
  for (i in seq_len(nrow(merged_table))) {
    pos_df <- top_snp_positions[[i]]
    if (nrow(pos_df) == 0) {
      genes_col[[i]] <- NA_character_
      next
    }

    genes <- unique(unlist(mapply(
      find_genes_for_snp,
      chr = as.character(pos_df$CHR),
      pos = as.numeric(pos_df$POS),
      SIMPLIFY = FALSE
    )))

    genes <- genes[nzchar(trimws(genes))]
    genes_col[[i]] <- if (length(genes) > 0) paste(genes, collapse = "; ") else NA_character_
  }

  merged_table$GENES_FROM_TOP_SNPS <- genes_col

  # Keep requested key fields first
  first_cols <- c("MERGED_LOCUS_ID", "CHR", "START", "END", "REPORTED_LOCUS")
  source_cols <- source_index_cols[source_index_cols %in% colnames(merged_table)]
  top_cols <- grep("^TOP_SNP_", colnames(merged_table), value = TRUE)
  last_cols <- c("GENES_FROM_TOP_SNPS")
  keep_cols <- c(first_cols, source_cols, top_cols, last_cols)
  keep_cols <- keep_cols[keep_cols %in% colnames(merged_table)]
  merged_table <- merged_table[, keep_cols, drop = FALSE]

  # Workbook output
  out_dir <- dirname(output_file)
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  }

  wb <- openxlsx2::wb_workbook()

  wb$add_worksheet("Summary")
  wb$add_data(sheet = "Summary", x = dataset_summary, start_row = 1, start_col = 1)
  wb$set_col_widths(sheet = "Summary", cols = 1:ncol(dataset_summary), widths = "auto")

  wb$add_worksheet("Loci_Table")
  wb$add_data(sheet = "Loci_Table", x = merged_table, start_row = 1, start_col = 1)

  # Style top SNP cells: light red when p-value < 5e-8
  top_snp_cols_in_output <- grep("^TOP_SNP_", colnames(merged_table), value = TRUE)
  if (length(top_snp_cols_in_output) > 0) {
    red_fill <- openxlsx2::wb_color(hex = "FFFFC7CE")

    extract_pvalue <- function(x) {
      if (is.na(x) || !nzchar(trimws(as.character(x)))) return(NA_real_)
      txt <- as.character(x)
      p_txt <- sub("\\s*\\(.*$", "", txt)
      suppressWarnings(as.numeric(p_txt))
    }

    for (top_col in top_snp_cols_in_output) {
      col_idx <- which(colnames(merged_table) == top_col)
      vals <- merged_table[[top_col]]
      pvals <- vapply(vals, extract_pvalue, numeric(1))
      hit_rows <- which(is.finite(pvals) & pvals < 5e-8)

      if (length(hit_rows) > 0) {
        for (row_idx in hit_rows) {
          cell_addr <- paste0(openxlsx2::int2col(col_idx), row_idx + 1L)
          wb$add_fill(sheet = "Loci_Table", dims = cell_addr, color = red_fill)
        }
      }
    }
  }

  # Style per-dataset source index cells: green if non-empty, no fill if empty.
  source_cols_in_output <- grep("^SOURCE_IDX_", colnames(merged_table), value = TRUE)
  if (length(source_cols_in_output) > 0) {
    green_fill <- openxlsx2::wb_color(hex = "FF92D050")
    for (src_col in source_cols_in_output) {
      col_idx <- which(colnames(merged_table) == src_col)
      non_empty_rows <- which(!is.na(merged_table[[src_col]]))
      if (length(non_empty_rows) > 0) {
        for (row_idx in non_empty_rows) {
          cell_addr <- paste0(openxlsx2::int2col(col_idx), row_idx + 1L)
          wb$add_fill(sheet = "Loci_Table", dims = cell_addr, color = green_fill)
        }
      }
    }
  }

  wb$set_col_widths(sheet = "Loci_Table", cols = 1:ncol(merged_table), widths = "auto")

  wb$save(output_file)

  message(sprintf("Aligned loci summary workbook saved: %s", output_file))
  message(sprintf("  Summary datasets: %d", nrow(dataset_summary)))
  message(sprintf("  Loci rows: %d", nrow(merged_table)))

  invisible(output_file)
}
