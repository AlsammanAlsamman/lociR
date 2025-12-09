#' Read Previously Reported Loci Files
#'
#' Reads files containing previously reported genomic loci from literature or databases.
#' Handles single or multiple files and returns a standardized data structure.
#'
#' @param file_paths Character vector of file paths to read
#' @param file_type Type of file format (default: "tab")
#' @param header Logical, whether file has header (default: TRUE)
#' @param custom_mapping Optional named list for custom column mapping
#' @param strip_extension Logical, whether to strip file extensions from SOURCE_FILE and list names (default: TRUE)
#'
#' @return If single file: data.frame with class "ReportedLoci"
#'         If multiple files: named list of data.frames with class "ReportedLoci"
#'
#' @export
read_reported_loci <- function(file_paths, 
                                file_type = "tab",
                                header = TRUE,
                                custom_mapping = NULL,
                                strip_extension = TRUE) {
  
  # Handle single vs multiple files
  if (length(file_paths) == 1) {
    return(read_reported_loci_single(file_paths, file_type, header, custom_mapping, strip_extension))
  } else {
    # Read multiple files
    result <- lapply(file_paths, function(f) {
      read_reported_loci_single(f, file_type, header, custom_mapping, strip_extension)
    })
    
    # Name the list based on file names
    if (strip_extension) {
      names(result) <- gsub("\\.[^.]*$", "", basename(file_paths))
    } else {
      names(result) <- basename(file_paths)
    }
    
    return(result)
  }
}


#' Read Single Previously Reported Loci File
#'
#' Internal function to read a single reported loci file
#'
#' @param file_path Path to the file
#' @param file_type Type of file format
#' @param header Whether file has header
#' @param custom_mapping Optional custom column mapping
#' @param strip_extension Whether to strip file extension from SOURCE_FILE
#'
#' @return data.frame with class "ReportedLoci"
#'
#' @keywords internal
read_reported_loci_single <- function(file_path, 
                                      file_type = "tab",
                                      header = TRUE,
                                      custom_mapping = NULL,
                                      strip_extension = TRUE) {
  
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  
  # Read the file based on type using fread for faster reading
  if (file_type == "tab") {
    data <- data.table::fread(file_path, 
                             header = header, 
                             sep = "\t", 
                             data.table = FALSE,
                             fill = TRUE)
  } else {
    stop("Unsupported file type: ", file_type)
  }
  
  message("Read ", nrow(data), " reported loci from: ", basename(file_path))
  
  # Standardize column names
  data <- standardize_reported_columns(data, custom_mapping)
  
  # Validate required columns
  required_cols <- c("CHR", "START", "END", "LOCUS")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns after standardization: ", 
         paste(missing_cols, collapse = ", "),
         "\nAvailable columns: ", paste(colnames(data), collapse = ", "))
  }
  
  # Convert chromosome to character and remove 'chr' prefix if present
  data$CHR <- as.character(data$CHR)
  data$CHR <- gsub("^chr", "", data$CHR, ignore.case = TRUE)
  
  # Convert positions to numeric
  data$START <- as.numeric(data$START)
  data$END <- as.numeric(data$END)
  
  # Calculate width if not present
  if (!"WIDTH" %in% colnames(data)) {
    data$WIDTH <- data$END - data$START + 1
  }
  
  # Add source file name (strip extension by default)
  if (strip_extension) {
    data$SOURCE_FILE <- gsub("\\.[^.]*$", "", basename(file_path))
  } else {
    data$SOURCE_FILE <- basename(file_path)
  }
  
  # Set class
  class(data) <- c("ReportedLoci", "data.frame")
  
  message("Successfully standardized ", nrow(data), " reported loci")
  
  return(data)
}


#' Standardize Reported Loci Column Names
#'
#' Standardizes column names for reported loci data to ensure consistency
#'
#' @param data data.frame with reported loci
#' @param custom_mapping Optional named list for custom column mapping
#'
#' @return data.frame with standardized column names
#'
#' @keywords internal
standardize_reported_columns <- function(data, custom_mapping = NULL) {
  
  # Define default column mapping for reported loci
  default_mapping <- list(
    CHR = c("chr", "chrom", "chromosome", "seqnames", "chr_name"),
    START = c("start", "locusstart", "locus_start", "region_start", "pos_start", "begin"),
    END = c("end", "locusend", "locus_end", "region_end", "pos_end", "stop"),
    LOCUS = c("locus", "locusname", "locus_name", "region", "region_name", "id", "locus_id"),
    GENE = c("gene", "bestguessgene", "best_guess_gene", "gene_name", "top_gene", "target_gene"),
    RSID = c("rsid", "rs_id", "snp", "snp_id", "lead_snp", "index_snp", "likelycausalsnp"),
    PVALUE = c("p", "pval", "pvalue", "p_value", "p.value")
  )
  
  # Merge with custom mapping if provided
  if (!is.null(custom_mapping)) {
    for (target in names(custom_mapping)) {
      default_mapping[[target]] <- c(custom_mapping[[target]], default_mapping[[target]])
    }
  }
  
  # Store original column names
  original_cols <- colnames(data)
  new_cols <- original_cols
  
  # Track which columns were mapped
  mapped_cols <- character(0)
  
  # Apply mapping
  for (target_name in names(default_mapping)) {
    possible_names <- default_mapping[[target_name]]
    
    # Find matching column (case-insensitive)
    match_idx <- which(tolower(original_cols) %in% tolower(possible_names))
    
    if (length(match_idx) > 0) {
      # Use first match
      match_idx <- match_idx[1]
      new_cols[match_idx] <- target_name
      mapped_cols <- c(mapped_cols, original_cols[match_idx])
      
      message("  Mapped '", original_cols[match_idx], "' -> '", target_name, "'")
    }
  }
  
  # Warn about unmapped columns
  unmapped <- setdiff(original_cols, mapped_cols)
  if (length(unmapped) > 0) {
    message("  Note: Unmapped columns retained as-is: ", paste(unmapped, collapse = ", "))
  }
  
  # Apply new column names
  colnames(data) <- new_cols
  
  return(data)
}


#' Print Method for ReportedLoci
#'
#' @param x ReportedLoci object
#' @param ... Additional arguments
#'
#' @export
print.ReportedLoci <- function(x, ...) {
  cat("Reported Loci Data\n")
  cat("==================\n")
  cat("Number of loci:", nrow(x), "\n")
  cat("Chromosomes:", length(unique(x$CHR)), "unique\n")
  
  if ("SOURCE_FILE" %in% colnames(x)) {
    cat("Source:", unique(x$SOURCE_FILE), "\n")
  }
  
  cat("\nColumns:", paste(colnames(x), collapse = ", "), "\n")
  cat("\nFirst few loci:\n")
  print.data.frame(head(x, 10))
  
  invisible(x)
}


#' Summary Method for ReportedLoci
#'
#' @param object ReportedLoci object
#' @param ... Additional arguments
#'
#' @export
summary.ReportedLoci <- function(object, ...) {
  cat("Reported Loci Summary\n")
  cat("=====================\n")
  cat("Total loci:", nrow(object), "\n")
  cat("Chromosomes:", paste(sort(unique(object$CHR)), collapse = ", "), "\n")
  
  if ("WIDTH" %in% colnames(object)) {
    cat("\nLocus width statistics (bp):\n")
    cat("  Min:", min(object$WIDTH, na.rm = TRUE), "\n")
    cat("  Mean:", round(mean(object$WIDTH, na.rm = TRUE), 2), "\n")
    cat("  Median:", median(object$WIDTH, na.rm = TRUE), "\n")
    cat("  Max:", max(object$WIDTH, na.rm = TRUE), "\n")
  }
  
  if ("GENE" %in% colnames(object)) {
    n_with_gene <- sum(!is.na(object$GENE) & object$GENE != "")
    cat("\nLoci with gene annotation:", n_with_gene, "(", 
        round(100 * n_with_gene / nrow(object), 1), "%)\n")
  }
  
  if ("SOURCE_FILE" %in% colnames(object)) {
    cat("\nSource file:", unique(object$SOURCE_FILE), "\n")
  }
  
  invisible(object)
}
