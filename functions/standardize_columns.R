#' Standardize column names in genomic data tables
#'
#' @param data Data frame with genomic data
#' @param column_map List. Custom column mapping to add/override defaults
#' @param strict Logical. If TRUE, stop on unmapped columns. If FALSE, keep original names
#' @param interactive Logical. If TRUE, prompt user for unmapped columns
#' @return Data frame with standardized column names
#' @export
standardize_columns <- function(data, 
                                column_map = NULL, 
                                strict = FALSE,
                                interactive = FALSE) {
  
  # Default column mapping: standard_name -> alternative names
  default_map <- list(
    CHR = c("chr", "chrom", "chromosome", "CHR", "CHROM", "Chromosome"),
    POS = c("pos", "position", "bp", "POS", "BP", "Position", "base_pair_location"),
    RSID = c("rsid", "snp", "rsID", "rs", "SNP", "RSID", "rsID", "MarkerName", "variant_id"),
    P = c("p", "pval", "p_value", "pvalue", "P", "PVAL", "P_VALUE", "Pvalue"),
    START = c("start", "Start", "START", "start_pos", "start_position"),
    END = c("end", "End", "END", "end_pos", "end_position"),
    LOCUS = c("locus", "Locus", "LOCUS", "GenomicLocus", "locus_id"),
    UNIQID = c("uniqid", "uniqID", "unique_id", "UniqueID"),
    NSNPS = c("nsnps", "nSNPs", "n_snps", "num_snps", "number_snps"),
    LEADSNPS = c("leadsnps", "LeadSNPs", "lead_snps", "lead_variants"),
    BETA = c("beta", "Beta", "BETA", "effect", "Effect"),
    SE = c("se", "SE", "std_err", "stderr", "standard_error"),
    OR = c("or", "OR", "odds_ratio", "OddsRatio"),
    EAF = c("eaf", "EAF", "maf", "MAF", "freq", "allele_freq"),
    A1 = c("a1", "A1", "effect_allele", "EA", "ea", "allele1"),
    A2 = c("a2", "A2", "other_allele", "OA", "oa", "allele2", "ref"),
    N = c("n", "N", "sample_size", "nsamples", "n_samples")
  )
  
  # Merge with custom mapping if provided
  if (!is.null(column_map)) {
    for (std_name in names(column_map)) {
      if (std_name %in% names(default_map)) {
        # Add to existing mapping
        default_map[[std_name]] <- unique(c(default_map[[std_name]], column_map[[std_name]]))
      } else {
        # New mapping
        default_map[[std_name]] <- column_map[[std_name]]
      }
    }
  }
  
  # Get current column names
  original_cols <- colnames(data)
  new_cols <- original_cols
  mapped_indices <- logical(length(original_cols))
  
  # Map columns
  for (i in seq_along(original_cols)) {
    col <- original_cols[i]
    
    # Find matching standard name
    for (std_name in names(default_map)) {
      if (col %in% default_map[[std_name]]) {
        new_cols[i] <- std_name
        mapped_indices[i] <- TRUE
        break
      }
    }
  }
  
  # Handle unmapped columns
  unmapped <- original_cols[!mapped_indices]
  
  if (length(unmapped) > 0) {
    # Always warn about unmapped columns
    warning("Unmapped column(s) found and will keep original name(s): ", 
            paste(unmapped, collapse = ", "),
            "\nConsider adding to column_map if standardization is needed.",
            call. = FALSE)
    
    if (interactive) {
      cat("\nUnmapped columns found:", paste(unmapped, collapse = ", "), "\n")
      cat("Options:\n")
      cat("  1. Keep original names\n")
      cat("  2. Add to existing mapping\n")
      cat("  3. Skip these columns\n")
      
      choice <- readline(prompt = "Choose option (1-3): ")
      
      if (choice == "2") {
        for (unmapped_col in unmapped) {
          cat("\nColumn:", unmapped_col, "\n")
          cat("Standard names:", paste(names(default_map), collapse = ", "), "\n")
          std_name <- readline(prompt = "Map to (or press Enter to keep original): ")
          
          if (std_name != "" && std_name %in% names(default_map)) {
            idx <- which(original_cols == unmapped_col)
            new_cols[idx] <- std_name
            mapped_indices[idx] <- TRUE
          }
        }
      } else if (choice == "3") {
        # Remove unmapped columns
    } else if (strict) {
      stop("Unmapped columns found: ", paste(unmapped, collapse = ", "), 
           "\nAdd them to column_map or set strict=FALSE")
    }
    # If not interactive and not strict, warning was already issued above stop("Unmapped columns found: ", paste(unmapped, collapse = ", "), 
           "\nAdd them to column_map or set strict=FALSE")
    } else {
      # Keep original names for unmapped columns (already set in new_cols)
      message("Keeping original names for unmapped columns: ", 
              paste(unmapped, collapse = ", "))
    }
  }
  
  # Check for duplicate column names after mapping
  if (any(duplicated(new_cols))) {
    duplicates <- new_cols[duplicated(new_cols)]
    warning("Duplicate standard names detected: ", paste(unique(duplicates), collapse = ", "),
            "\nKeeping first occurrence only")
    
    # Keep only first occurrence of each name
    keep_idx <- !duplicated(new_cols)
    data <- data[, keep_idx, drop = FALSE]
    new_cols <- new_cols[keep_idx]
    original_cols <- original_cols[keep_idx]
  }
  
  # Apply new column names
  colnames(data) <- new_cols
  
  # Print mapping summary
  changes <- original_cols != new_cols
  if (any(changes)) {
    cat("\nColumn mapping applied:\n")
    for (i in which(changes)) {
      cat(sprintf("  %s -> %s\n", original_cols[i], new_cols[i]))
    }
  }
  
  return(data)
}


#' Get current column mapping
#'
#' @return List of standard column names and their alternatives
#' @export
get_column_mapping <- function() {
  list(
    CHR = c("chr", "chrom", "chromosome", "CHR", "CHROM", "Chromosome"),
    POS = c("pos", "position", "bp", "POS", "BP", "Position", "base_pair_location"),
    RSID = c("rsid", "snp", "rsID", "rs", "SNP", "RSID", "rsID", "MarkerName", "variant_id"),
    P = c("p", "pval", "p_value", "pvalue", "P", "PVAL", "P_VALUE", "Pvalue"),
    START = c("start", "Start", "START", "start_pos", "start_position"),
    END = c("end", "End", "END", "end_pos", "end_position"),
    LOCUS = c("locus", "Locus", "LOCUS", "GenomicLocus", "locus_id"),
    UNIQID = c("uniqid", "uniqID", "unique_id", "UniqueID"),
    NSNPS = c("nsnps", "nSNPs", "n_snps", "num_snps", "number_snps"),
    LEADSNPS = c("leadsnps", "LeadSNPs", "lead_snps", "lead_variants"),
    BETA = c("beta", "Beta", "BETA", "effect", "Effect"),
    SE = c("se", "SE", "std_err", "stderr", "standard_error"),
    OR = c("or", "OR", "odds_ratio", "OddsRatio"),
    EAF = c("eaf", "EAF", "maf", "MAF", "freq", "allele_freq"),
    A1 = c("a1", "A1", "effect_allele", "EA", "ea", "allele1"),
    A2 = c("a2", "A2", "other_allele", "OA", "oa", "allele2", "ref"),
    N = c("n", "N", "sample_size", "nsamples", "n_samples")
  )
}


#' Add custom column mapping
#'
#' @param standard_name Character. The standard column name to use
#' @param alternatives Character vector. Alternative names for this column
#' @return List with the new mapping
#' @export
add_column_mapping <- function(standard_name, alternatives) {
  mapping <- list()
  mapping[[standard_name]] <- alternatives
  
  cat("New mapping created:\n")
  cat(sprintf("  %s <- %s\n", standard_name, paste(alternatives, collapse = ", ")))
  
  return(mapping)
}
