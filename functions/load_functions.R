#' Load all functions from the functions directory
#'
#' @param functions_dir Character string. Path to functions directory
#' @return Invisible NULL. Functions are loaded into the global environment
#' @export
load_all_functions <- function(functions_dir = "functions") {
  
  # Get all R files in the functions directory
  r_files <- list.files(functions_dir, 
                        pattern = "\\.R$", 
                        full.names = TRUE,
                        ignore.case = TRUE)
  
  if (length(r_files) == 0) {
    warning("No R files found in: ", functions_dir)
    return(invisible(NULL))
  }
  
  # Source each file
  for (file in r_files) {
    cat("Loading:", basename(file), "\n")
    source(file)
  }
  
  cat("\nLoaded", length(r_files), "function file(s)\n")
  
  return(invisible(NULL))
}
