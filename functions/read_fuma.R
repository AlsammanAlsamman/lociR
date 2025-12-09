#' Read FUMA genomic loci file(s)
#'
#' @param file_path Character string or vector. Path(s) to FUMA output file(s)
#' @param strip_extension Character vector. File extensions to remove from names. Default: c(".txt", ".tsv")
#' @return If single file: Object of class 'FUMA'. If multiple files: List of FUMA objects
#' @export
read_fuma <- function(file_path, strip_extension = c(".txt", ".tsv")) {
  
  # Handle multiple files
  if (length(file_path) > 1) {
    message(sprintf("Reading %d FUMA files...\n", length(file_path)))
    
    fuma_list <- lapply(seq_along(file_path), function(i) {
      message(sprintf("\n[%d/%d] Processing: %s", i, length(file_path), basename(file_path[i])))
      read_fuma_single(file_path[i])
    })
    
    # Name the list elements by file names (without extensions)
    file_names <- basename(file_path)
    
    # Remove specified extensions
    if (!is.null(strip_extension) && length(strip_extension) > 0) {
      for (ext in strip_extension) {
        # Escape special regex characters in extension
        ext_pattern <- paste0(gsub("([.|()\\^{}+$*?])", "\\\\\\1", ext), "$")
        file_names <- sub(ext_pattern, "", file_names, ignore.case = TRUE)
      }
    }
    
    names(fuma_list) <- file_names
    
    message("\n=== All files loaded successfully ===\n")
    return(fuma_list)
  }
  
  # Single file
  return(read_fuma_single(file_path))
}


#' Read a single FUMA genomic loci file
#'
#' @param file_path Character string. Path to FUMA output file
#' @return Object of class 'FUMA' containing the loci data and metadata
#' @keywords internal
read_fuma_single <- function(file_path) {
  
  # Check file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read the file using fread for faster reading
  fuma_data <- data.table::fread(file_path, 
                                 header = TRUE, 
                                 sep = "\t",
                                 data.table = FALSE)
  
  # Standardize column names
  message("\n--- Standardizing column names ---")
  fuma_data <- standardize_columns(fuma_data, strict = FALSE)
  message("--- Standardization complete ---\n")
  
  # Validate required columns (after standardization)
  required_cols <- c("LOCUS", "CHR", "POS", "P", "START", "END", "RSID")
  missing_cols <- setdiff(required_cols, colnames(fuma_data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Create FUMA object
  fuma_obj <- structure(
    list(
      data = fuma_data,
      metadata = list(
        file_path = file_path,
        file_name = basename(file_path),
        n_loci = nrow(fuma_data),
        date_read = Sys.time(),
        chromosomes = unique(fuma_data$CHR)
      )
    ),
    class = "FUMA"
  )
  
  return(fuma_obj)
}


#' Print method for FUMA objects
#' @export
print.FUMA <- function(x, ...) {
  cat("FUMA Genomic Loci Object\n")
  cat("========================\n")
  cat("File:", x$metadata$file_name, "\n")
  cat("Number of loci:", x$metadata$n_loci, "\n")
  cat("Chromosomes:", paste(sort(x$metadata$chromosomes), collapse = ", "), "\n")
  cat("Date read:", format(x$metadata$date_read, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("(Column names standardized)\n")
  cat("\nFirst few rows:\n")
  print(head(x$data, 3))
}


#' Summary method for FUMA objects
#' @export
summary.FUMA <- function(object, ...) {
  cat("FUMA Genomic Loci Summary\n")
  cat("=========================\n")
  cat("Total loci:", object$metadata$n_loci, "\n")
  cat("Chromosomes:", length(object$metadata$chromosomes), "\n")
  # Summary statistics
  cat("\nP-value range:", format(min(object$data$P), scientific = TRUE), 
      "to", format(max(object$data$P), scientific = TRUE), "\n")
  
  cat("Loci size (bp):\n")
  loci_size <- object$data$END - object$data$START
  print(summary(loci_size))
  
  if ("NSNPS" %in% colnames(object$data)) {
    cat("\nSNPs per locus:\n")
    print(summary(object$data$NSNPS))
  }
  
  cat("\nLoci per chromosome:\n")
  print(table(object$data$CHR))
}
