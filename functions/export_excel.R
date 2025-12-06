#' Export merged FUMA loci to Excel with visual representation
#'
#' @param merged_loci Data frame. Output from merge_fuma_loci()
#' @param fuma_list List of FUMA objects. Original FUMA data to get detailed information
#' @param output_file Character. Output Excel file path (required)
#' @param color_palette Character vector. Colors for sources. Default: RColorBrewer Set2
#' @export
export_merged_to_excel <- function(merged_loci, fuma_list, output_file) {
  
  if (missing(output_file)) {
    stop("output_file is required. Please provide an Excel file path (e.g., 'merged_loci.xlsx')")
  }
  
  if (missing(fuma_list)) {
    stop("fuma_list is required to extract detailed locus information")
  }
  
  # Convert single FUMA to list
  if ("FUMA" %in% class(fuma_list)) {
    fuma_list <- list(fuma_list)
    names(fuma_list) <- fuma_list[[1]]$metadata$file_name
  }
  
  message(sprintf("Creating Excel file for %d merged loci...", nrow(merged_loci)))
  
  # Create a mapping of all original loci with full details
  all_loci_details <- do.call(rbind, lapply(names(fuma_list), function(source_name) {
    fuma <- fuma_list[[source_name]]
    data <- fuma$data
    
    data$SOURCE <- source_name
    data$LOCUS_KEY <- paste0(data$LOCUS, "(", source_name, ":", data$START, ":", data$END, ")")
    
    return(data)
  }))
  
  # Prepare data for Excel
  excel_data <- list()
  current_row <- 1
  
  # Define color palette
  unique_sources <- unique(unlist(strsplit(merged_loci$SOURCES, ";")))
  n_sources <- length(unique_sources)
  
  if (n_sources <= 8) {
    colors <- RColorBrewer::brewer.pal(max(3, n_sources), "Set2")
  } else {
    colors <- rainbow(n_sources)
  }
  source_colors <- setNames(colors[seq_len(n_sources)], unique_sources)
  
  # Track row positions for merging cells
  merge_ranges <- list()
  
  # Process each merged locus
  for (i in seq_len(nrow(merged_loci))) {
    row <- merged_loci[i, ]
    
    # Parse original loci
    orig_loci_keys <- strsplit(row$ORIGINAL_LOCI, ";")[[1]]
    
    # Store starting row for this merged locus
    start_row <- length(excel_data) + 1
    
    # Extract details for each original locus
    for (j in seq_along(orig_loci_keys)) {
      locus_key <- orig_loci_keys[j]
      orig_detail <- all_loci_details[all_loci_details$LOCUS_KEY == locus_key, ]
      
      if (nrow(orig_detail) > 0) {
        orig_detail <- orig_detail[1, ]  # Take first match
        
        excel_row <- data.frame(
          MERGED_LOCUS_ID = if(j == 1) paste0("MERGED_", row$MERGED_LOCUS_ID) else "",
          SOURCE = orig_detail$SOURCE,
          LOCUS = orig_detail$LOCUS,
          CHR = orig_detail$CHR,
          START = orig_detail$START,
          END = orig_detail$END,
          WIDTH = orig_detail$END - orig_detail$START,
          stringsAsFactors = FALSE
        )
        
        excel_data[[length(excel_data) + 1]] <- excel_row
      }
    }
    
    # Add MERGED row
    end_row <- length(excel_data) + 1
    
    # Store merge range for first column (MERGED_LOCUS_ID)
    if (end_row > start_row) {
      merge_ranges[[length(merge_ranges) + 1]] <- c(start_row, end_row)
    }
    
    merged_row <- data.frame(
      MERGED_LOCUS_ID = "",
      SOURCE = "MERGED",
      LOCUS = paste0("chr", row$CHR, "_merged"),
      CHR = row$CHR,
      START = row$START,
      END = row$END,
      WIDTH = row$WIDTH,
      stringsAsFactors = FALSE
    )
    
    excel_data[[length(excel_data) + 1]] <- merged_row
  }
  
  # Combine all data
  final_data <- do.call(rbind, excel_data)
  
  # Pre-compute which merged locus each row belongs to
  row_to_merged_id <- character(nrow(final_data))
  current_id <- ""
  for (i in seq_len(nrow(final_data))) {
    if (!is.na(final_data$MERGED_LOCUS_ID[i]) && final_data$MERGED_LOCUS_ID[i] != "") {
      current_id <- final_data$MERGED_LOCUS_ID[i]
    }
    row_to_merged_id[i] <- current_id
  }
  
  # Create workbook
  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("Merged_Loci")
  
  # Write data
  wb$add_data(sheet = "Merged_Loci", x = final_data, start_row = 1, start_col = 1)
  
  # Merge cells and apply alternating background colors
  for (i in seq_along(merge_ranges)) {
    merge_range <- merge_ranges[[i]]
    start_excel_row <- merge_range[1] + 1  # +1 for header
    end_excel_row <- merge_range[2] + 1
    
    if (end_excel_row > start_excel_row) {
      merge_dims <- paste0("A", start_excel_row, ":A", end_excel_row)
      wb$merge_cells(sheet = "Merged_Loci", dims = merge_dims)
      
      # Center align the merged cell
      wb$add_cell_style(sheet = "Merged_Loci", 
                       dims = paste0("A", start_excel_row),
                       horizontal = "center",
                       vertical = "center")
    }
    
    # Apply alternating background color for entire locus group
    # Odd groups: darker gray, Even groups: white
    if (i %% 2 == 1) {
      # Darker gray background for odd groups
      for (row in start_excel_row:end_excel_row) {
        row_range <- paste0("A", row, ":", openxlsx2::int2col(ncol(final_data)), row)
        wb$add_fill(sheet = "Merged_Loci", dims = row_range, 
                   color = openxlsx2::wb_color(hex = "FFE0E0E0"))
      }
    }
    # Even groups remain white (default)
  }
  
  # Number of visualization columns
  n_vis_cols <- 20
  vis_start_col <- ncol(final_data) + 1
  
  # Create a mapping of merged locus boundaries
  merged_boundaries <- list()
  for (i in seq_len(nrow(merged_loci))) {
    merged_id <- paste0("MERGED_", merged_loci$MERGED_LOCUS_ID[i])
    merged_boundaries[[merged_id]] <- c(merged_loci$START[i], merged_loci$END[i])
  }
  
  # Apply formatting and colors
  for (row_idx in seq_len(nrow(final_data))) {
    excel_row <- row_idx + 1  # +1 for header row
    
    source <- final_data$SOURCE[row_idx]
    
    # Get merged ID from pre-computed mapping
    merged_id <- row_to_merged_id[row_idx]
    
    # Get boundaries for this merged locus
    merged_start <- NA
    merged_end <- NA
    if (merged_id != "" && merged_id %in% names(merged_boundaries)) {
      merged_start <- merged_boundaries[[merged_id]][1]
      merged_end <- merged_boundaries[[merged_id]][2]
    }
    
    if (!is.na(source) && source != "") {
      # Determine color
      if (source == "MERGED") {
        fill_color <- "FF2F4F4F"  # Dark gray
      } else if (source %in% names(source_colors)) {
        # Convert R color to hex
        fill_color <- paste0("FF", substr(rgb(col2rgb(source_colors[source])[1,1]/255,
                                              col2rgb(source_colors[source])[2,1]/255,
                                              col2rgb(source_colors[source])[3,1]/255), 2, 7))
      } else {
        fill_color <- "FFD3D3D3"  # Light gray
      }
      
      # Get start and end positions for this row
      start_pos <- final_data$START[row_idx]
      end_pos <- final_data$END[row_idx]
      
      if (!is.na(start_pos) && !is.na(end_pos) && !is.na(merged_start) && !is.na(merged_end)) {
        # Calculate positions relative to merged region
        merged_width <- merged_end - merged_start
        
        if (merged_width > 0) {
          # Scale positions to visualization columns (relative to merged region)
          rel_start <- (start_pos - merged_start) / merged_width
          rel_end <- (end_pos - merged_start) / merged_width
          
          # Clamp to [0, 1] range
          rel_start <- max(0, min(1, rel_start))
          rel_end <- max(0, min(1, rel_end))
          
          # Convert to column indices
          scaled_start <- vis_start_col + floor(rel_start * n_vis_cols)
          scaled_end <- vis_start_col + ceiling(rel_end * n_vis_cols) - 1
          
  # Auto-size columns
  wb$set_col_widths(sheet = "Merged_Loci", cols = 1:ncol(final_data), widths = "auto")
  
  # Set visualization columns to narrower width (20 columns total)
  vis_cols <- vis_start_col:(vis_start_col + n_vis_cols - 1)
  wb$set_col_widths(sheet = "Merged_Loci", cols = vis_cols, widths = 2)
          for (col_idx in scaled_start:min(scaled_end, vis_start_col + n_vis_cols - 1)) {
            cell_ref <- openxlsx2::int2col(col_idx)
            cell_address <- paste0(cell_ref, excel_row)
            
            wb$add_fill(sheet = "Merged_Loci", dims = cell_address, color = openxlsx2::wb_color(hex = fill_color))
          }
        }
      }
      
      # Highlight MERGED rows with bold
      if (source == "MERGED") {
        row_range <- paste0("A", excel_row, ":", openxlsx2::int2col(ncol(final_data)), excel_row)
        wb$add_font(sheet = "Merged_Loci", dims = row_range, bold = TRUE)
      }
    }
  }
  
  # Auto-size columns
  wb$set_col_widths(sheet = "Merged_Loci", cols = 1:ncol(final_data), widths = "auto")
  
  # Set visualization columns to narrower width
  vis_cols <- (ncol(final_data) + 1):(ncol(final_data) + 100)
  wb$set_col_widths(sheet = "Merged_Loci", cols = vis_cols, widths = 1)
  
  # Save workbook
  wb$save(output_file)
  
  message(sprintf("Excel file saved to: %s", output_file))
  message(sprintf("  Total rows: %d", nrow(final_data)))
  message(sprintf("  Merged loci: %d", nrow(merged_loci)))
  
  return(invisible(output_file))
}
