#' Read GWAS summary statistics file(s)
#'
#' Reads GWAS summary statistics with automatic column standardization,
#' p-value filtering, and selection of important columns.
#'
#' @param file_path Character string or vector. Path(s) to GWAS file(s)
#' @param name Character string or vector. Name(s) for the GWAS dataset(s). 
#'        Used for linking and identification. If NULL, uses filename.
#' @param pvalue_threshold Numeric. Maximum p-value threshold for filtering SNPs. Default: 5e-4
#' @param sep Character. Field separator. Default: "\t" (tab-separated)
#' @param header Logical. Whether file has header. Default: TRUE
#' @return If single file: Object of class 'GWAS'. If multiple files: Named list of GWAS objects
#' @export
read_gwas <- function(file_path, 
                      name = NULL,
                      pvalue_threshold = 5e-4, 
                      sep = "\t",
                      header = TRUE) {
  
  # Handle multiple files
  if (length(file_path) > 1) {
    message(sprintf("Reading %d GWAS files...\n", length(file_path)))
    
    # Set default names if not provided
    if (is.null(name)) {
      name <- tools::file_path_sans_ext(basename(file_path))
    } else if (length(name) != length(file_path)) {
      stop("Number of names must match number of file paths")
    }
    
    gwas_list <- lapply(seq_along(file_path), function(i) {
      message(sprintf("\n[%d/%d] Processing: %s", i, length(file_path), basename(file_path[i])))
      read_gwas_single(file_path[i], name[i], pvalue_threshold, sep, header)
    })
    
    names(gwas_list) <- name
    
    message("\n=== All GWAS files loaded successfully ===\n")
    return(gwas_list)
  }
  
  # Single file
  if (is.null(name)) {
    name <- tools::file_path_sans_ext(basename(file_path))
  }
  
  return(read_gwas_single(file_path, name, pvalue_threshold, sep, header))
}


#' Read a single GWAS summary statistics file
#'
#' @param file_path Character string. Path to GWAS file
#' @param name Character string. Name for the GWAS dataset
#' @param pvalue_threshold Numeric. Maximum p-value threshold
#' @param sep Character. Field separator
#' @param header Logical. Whether file has header
#' @return Object of class 'GWAS'
#' @keywords internal
read_gwas_single <- function(file_path, name, pvalue_threshold, sep, header) {
  
  # Check file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read the file using fread for faster reading
  message("Reading file...")
  gwas_data <- data.table::fread(file_path, 
                                 header = header, 
                                 sep = sep,
                                 data.table = FALSE)
  
  message(sprintf("  Initial SNPs: %s", format(nrow(gwas_data), big.mark = ",")))
  
  # Standardize column names
  message("\n--- Standardizing column names ---")
  gwas_data <- standardize_gwas_columns(gwas_data)
  message("--- Standardization complete ---\n")
  
  # Validate required columns
  required_cols <- c("CHR", "POS", "P")
  missing_cols <- setdiff(required_cols, colnames(gwas_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Filter by p-value threshold
  message(sprintf("Filtering SNPs with p-value <= %s...", format(pvalue_threshold, scientific = TRUE)))
  initial_count <- nrow(gwas_data)
  gwas_data <- gwas_data[gwas_data$P <= pvalue_threshold, ]
  filtered_count <- nrow(gwas_data)
  removed_count <- initial_count - filtered_count
  
  message(sprintf("  Retained: %s SNPs (%.2f%%)", 
                  format(filtered_count, big.mark = ","),
                  100 * filtered_count / initial_count))
  message(sprintf("  Removed:  %s SNPs (%.2f%%)", 
                  format(removed_count, big.mark = ","),
                  100 * removed_count / initial_count))
  
  if (filtered_count == 0) {
    warning("No SNPs remain after p-value filtering!")
  }
  
  # Select and order important columns
  important_cols <- c("CHR", "POS", "RSID", "EFFECT_ALLELE", "OTHER_ALLELE", 
                     "BETA", "SE", "OR", "P")
  
  # Keep only columns that exist
  keep_cols <- intersect(important_cols, colnames(gwas_data))
  gwas_data <- gwas_data[, keep_cols, drop = FALSE]
  
  message(sprintf("\nRetained columns: %s", paste(keep_cols, collapse = ", ")))
  
  # Create GWAS object
  gwas_obj <- structure(
    list(
      name = name,
      data = gwas_data,
      metadata = list(
        file_path = file_path,
        file_name = basename(file_path),
        n_snps = nrow(gwas_data),
        pvalue_threshold = pvalue_threshold,
        date_read = Sys.time(),
        chromosomes = sort(unique(gwas_data$CHR)),
        p_range = if(nrow(gwas_data) > 0) c(min(gwas_data$P), max(gwas_data$P)) else c(NA, NA)
      )
    ),
    class = "GWAS"
  )
  
  message(sprintf("\n=== GWAS dataset '%s' loaded successfully ===\n", name))
  
  return(gwas_obj)
}


#' Standardize GWAS column names
#'
#' @param data Data frame with GWAS data
#' @return Data frame with standardized column names
#' @keywords internal
standardize_gwas_columns <- function(data) {
  
  # Define column mapping for GWAS-specific columns
  gwas_mapping <- list(
    CHR = c("chr", "chrom", "chromosome", "CHR", "CHROM", "Chromosome"),
    POS = c("pos", "position", "bp", "POS", "BP", "Position", "base_pair_location"),
    RSID = c("rsid", "snp", "rsID", "rs", "SNP", "RSID", "MarkerName", "variant_id", "marker"),
    P = c("p", "pval", "p_value", "pvalue", "P", "PVAL", "P_VALUE", "Pvalue", "P-value"),
    EFFECT_ALLELE = c("effectallele", "effect_allele", "allele1", "a1", "A1", "ea", "EA", "Allele1"),
    OTHER_ALLELE = c("otherallele", "other_allele", "allele2", "a2", "A2", "oa", "OA", "Allele2", "ref"),
    BETA = c("beta", "Beta", "BETA", "effect", "Effect"),
    SE = c("se", "SE", "std_err", "stderr", "standard_error", "StdErr"),
    OR = c("or", "OR", "odds_ratio", "OddsRatio")
  )
  
  # Get current column names
  original_cols <- colnames(data)
  new_cols <- original_cols
  mapped_indices <- logical(length(original_cols))
  
  # Map columns (case-insensitive)
  for (i in seq_along(original_cols)) {
    col <- original_cols[i]
    
    # Find matching standard name
    for (std_name in names(gwas_mapping)) {
      if (tolower(col) %in% tolower(gwas_mapping[[std_name]])) {
        new_cols[i] <- std_name
        mapped_indices[i] <- TRUE
        message(sprintf("  Mapped '%s' -> '%s'", col, std_name))
        break
      }
    }
  }
  
  # Handle unmapped columns
  unmapped <- original_cols[!mapped_indices]
  
  if (length(unmapped) > 0) {
    message(sprintf("  Note: Unmapped columns will be removed: %s", 
                    paste(unmapped, collapse = ", ")))
  }
  
  # Check for duplicate column names after mapping
  if (any(duplicated(new_cols[mapped_indices]))) {
    duplicates <- new_cols[mapped_indices][duplicated(new_cols[mapped_indices])]
    warning("Duplicate standard names detected: ", paste(unique(duplicates), collapse = ", "),
            "\nKeeping first occurrence only")
    
    # Keep only first occurrence of each name
    keep_idx <- !duplicated(new_cols) & mapped_indices
    data <- data[, keep_idx, drop = FALSE]
    new_cols <- new_cols[keep_idx]
  } else {
    # Keep only mapped columns
    data <- data[, mapped_indices, drop = FALSE]
    new_cols <- new_cols[mapped_indices]
  }
  
  # Apply new column names
  colnames(data) <- new_cols
  
  # Convert CHR to character (remove 'chr' prefix if present)
  if ("CHR" %in% colnames(data)) {
    data$CHR <- as.character(data$CHR)
    data$CHR <- gsub("^chr", "", data$CHR, ignore.case = TRUE)
  }
  
  # Convert numeric columns
  numeric_cols <- c("POS", "BETA", "SE", "OR", "P")
  for (col in intersect(numeric_cols, colnames(data))) {
    data[[col]] <- as.numeric(data[[col]])
  }
  
  return(data)
}


#' Print method for GWAS objects
#' @export
print.GWAS <- function(x, ...) {
  cat("GWAS Summary Statistics\n")
  cat("=======================\n")
  cat("Name:", x$name, "\n")
  cat("File:", x$metadata$file_name, "\n")
  cat("Number of SNPs:", format(x$metadata$n_snps, big.mark = ","), "\n")
  cat("P-value threshold:", format(x$metadata$pvalue_threshold, scientific = TRUE), "\n")
  cat("Chromosomes:", paste(x$metadata$chromosomes, collapse = ", "), "\n")
  if (!is.na(x$metadata$p_range[1])) {
    cat("P-value range:", format(x$metadata$p_range[1], scientific = TRUE), 
        "to", format(x$metadata$p_range[2], scientific = TRUE), "\n")
  }
  cat("Date read:", format(x$metadata$date_read, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Columns:", paste(colnames(x$data), collapse = ", "), "\n")
  cat("\nFirst few SNPs:\n")
  print(head(x$data, 3))
}


#' Summary method for GWAS objects
#' @export
summary.GWAS <- function(object, ...) {
  cat("GWAS Summary Statistics - Detailed Summary\n")
  cat("==========================================\n")
  cat("Dataset name:", object$name, "\n")
  cat("Total SNPs:", format(object$metadata$n_snps, big.mark = ","), "\n")
  cat("Chromosomes:", length(object$metadata$chromosomes), "\n")
  
  cat("\nSNPs per chromosome:\n")
  chr_table <- table(object$data$CHR)
  print(chr_table)
  
  if ("P" %in% colnames(object$data)) {
    cat("\nP-value distribution:\n")
    cat("  Min:   ", format(min(object$data$P), scientific = TRUE), "\n")
    cat("  Q1:    ", format(quantile(object$data$P, 0.25), scientific = TRUE), "\n")
    cat("  Median:", format(median(object$data$P), scientific = TRUE), "\n")
    cat("  Q3:    ", format(quantile(object$data$P, 0.75), scientific = TRUE), "\n")
    cat("  Max:   ", format(max(object$data$P), scientific = TRUE), "\n")
    
    # Count genome-wide significant SNPs
    gw_sig <- sum(object$data$P < 5e-8)
    cat(sprintf("\nGenome-wide significant SNPs (p < 5e-8): %s (%.2f%%)\n", 
                format(gw_sig, big.mark = ","),
                100 * gw_sig / nrow(object$data)))
  }
  
  if ("BETA" %in% colnames(object$data)) {
    cat("\nEffect size (BETA) statistics:\n")
    print(summary(object$data$BETA))
  }
  
  if ("OR" %in% colnames(object$data)) {
    cat("\nOdds Ratio (OR) statistics:\n")
    print(summary(object$data$OR))
  }
  
  invisible(object)
}
