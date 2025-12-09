#' Export refined loci to Excel with visual representation
#'
#' Creates an Excel report showing how merged loci were refined/split using reported loci.
#' Visualizes the splitting with colored cells and groups by original merged locus.
#'
#' @param refined_loci Data frame. Output from refine_with_reported()
#' @param output_file Character. Output Excel file path (required)
#' @param color_palette Character vector. Colors for refinement types. Default: custom palette
#' @export
export_refined_to_excel <- function(refined_loci, output_file) {
  
  if (missing(output_file)) {
    stop("output_file is required. Please provide an Excel file path (e.g., 'refined_loci.xlsx')")
  }
  
  message(sprintf("Creating Excel file for %d refined loci from %d merged loci...", 
                  nrow(refined_loci), 
                  length(unique(refined_loci$MERGED_LOCUS_ID))))
  
  # Prepare data for Excel - group by MERGED_LOCUS_ID
  excel_data <- list()
  merge_ranges <- list()
  
  # Get unique merged loci
  unique_merged <- unique(refined_loci$MERGED_LOCUS_ID)
  
  # Define colors for refinement types
  refinement_colors <- list(
    "split" = "FF87CEEB",      # Sky blue
    "expanded" = "FF90EE90",   # Light green
    "aligned" = "FFFFA07A",    # Light salmon
    "unchanged" = "FFD3D3D3",  # Light gray
    "source" = "FFFFD700",     # Gold for source loci
    "merged_summary" = "FF2F4F4F",  # Dark gray for merged summary rows
    "merged_row" = "FF2F4F4F"  # Dark gray for merged summary rows (alias)
  )
  
  # Process each merged locus
  for (merged_id in unique_merged) {
    # Get all refined loci for this merged locus
    refined_subset <- refined_loci[refined_loci$MERGED_LOCUS_ID == merged_id, ]
    refined_subset <- refined_subset[order(refined_subset$START), ]
    
    # Store starting row for this merged locus
    start_row <- length(excel_data) + 1
    
    # Parse and add original source loci first
    original_loci_keys <- strsplit(refined_subset$ORIGINAL_LOCI[1], ";")[[1]]
    
    for (locus_key in original_loci_keys) {
      # Parse: "locus_id(source:start:end)"
      matches <- regmatches(locus_key, regexec("([^(]+)\\(([^:]+):([^:]+):([^)]+)\\)", locus_key))
      
      if (length(matches[[1]]) == 5) {
        locus_id <- matches[[1]][2]
        source_name <- matches[[1]][3]
        locus_start <- as.numeric(matches[[1]][4])
        locus_end <- as.numeric(matches[[1]][5])
        
        excel_row <- data.frame(
          MERGED_LOCUS_ID = if(length(excel_data) == start_row - 1) merged_id else "",
          REFINED_LOCUS_ID = locus_id,
          CHR = refined_subset$CHR[1],
          START = locus_start,
          END = locus_end,
          WIDTH = locus_end - locus_start + 1,
          REPORTED_LOCUS = source_name,
          REPORTED_SOURCE = "FUMA",
          REFINEMENT_TYPE = "source",
          N_ORIGINAL_LOCI = refined_subset$N_ORIGINAL_LOCI[1],
          SOURCES = refined_subset$SOURCES[1],
          stringsAsFactors = FALSE
        )
        
        excel_data[[length(excel_data) + 1]] <- excel_row
      }
    }
    
    # Calculate MERGED boundaries from all rows added so far (including source loci)
    all_starts <- sapply(excel_data[start_row:length(excel_data)], function(x) x$START)
    all_ends <- sapply(excel_data[start_row:length(excel_data)], function(x) x$END)
    merged_start <- min(all_starts, na.rm = TRUE)
    merged_end <- max(all_ends, na.rm = TRUE)
    merged_width <- merged_end - merged_start + 1
    
    # Add MERGED summary row in the middle
    
    merged_summary_row <- data.frame(
      MERGED_LOCUS_ID = "",
      REFINED_LOCUS_ID = paste0(merged_id, "_MERGED"),
      CHR = refined_subset$CHR[1],
      START = merged_start,
      END = merged_end,
      WIDTH = merged_width,
      REPORTED_LOCUS = "MERGED",
      REPORTED_SOURCE = "MERGED",
      REFINEMENT_TYPE = "merged_summary",
      N_ORIGINAL_LOCI = refined_subset$N_ORIGINAL_LOCI[1],
      SOURCES = refined_subset$SOURCES[1],
      stringsAsFactors = FALSE
    )
    
    excel_data[[length(excel_data) + 1]] <- merged_summary_row
    
    # Add each refined locus (split by reported)
    for (i in seq_len(nrow(refined_subset))) {
      row <- refined_subset[i, ]
      
      excel_row <- data.frame(
        MERGED_LOCUS_ID = "",
        REFINED_LOCUS_ID = row$REFINED_LOCUS_ID,
        CHR = row$CHR,
        START = row$START,
        END = row$END,
        WIDTH = row$WIDTH,
        REPORTED_LOCUS = if(!is.na(row$REPORTED_LOCUS)) row$REPORTED_LOCUS else "N/A",
        REPORTED_SOURCE = if(!is.na(row$REPORTED_SOURCE)) row$REPORTED_SOURCE else "N/A",
        REFINEMENT_TYPE = row$REFINEMENT_TYPE,
        N_ORIGINAL_LOCI = row$N_ORIGINAL_LOCI,
        SOURCES = row$SOURCES,
        stringsAsFactors = FALSE
      )
      
      excel_data[[length(excel_data) + 1]] <- excel_row
    }
    
    # Calculate ending row for this merged locus group
    end_row <- length(excel_data)
    
    # Store merge range for first column (MERGED_LOCUS_ID)
    if (end_row > start_row) {
      merge_ranges[[length(merge_ranges) + 1]] <- list(
        start = start_row,
        end = end_row,
        merged_id = merged_id,
        merged_start = merged_start,
        merged_end = merged_end,
        merged_width = merged_width
      )
    }
  }
  
  # Combine all data
  final_data <- do.call(rbind, excel_data)
  
  # Create workbook
  wb <- openxlsx2::wb_workbook()
  wb$add_worksheet("Refined_Loci")
  
  # Write data
  wb$add_data(sheet = "Refined_Loci", x = final_data, start_row = 1, start_col = 1)
  
  # Number of visualization columns
  n_vis_cols <- 20
  vis_start_col <- ncol(final_data) + 1
  
  # Process each merged locus group
  for (merge_idx in seq_along(merge_ranges)) {
    merge_info <- merge_ranges[[merge_idx]]
    start_excel_row <- merge_info$start + 1  # +1 for header
    end_excel_row <- merge_info$end + 1
    
    # Merge cells for MERGED_LOCUS_ID column
    if (end_excel_row > start_excel_row) {
      merge_dims <- paste0("A", start_excel_row, ":A", end_excel_row)
      wb$merge_cells(sheet = "Refined_Loci", dims = merge_dims)
      
      # Center align the merged cell
      wb$add_cell_style(sheet = "Refined_Loci", 
                       dims = paste0("A", start_excel_row),
                       horizontal = "center",
                       vertical = "center")
    }
    
    # Note: Row background colors are now applied per-row based on refinement type
    # (see below in the visualization loop)
    
    # Add visualization for each refined locus in this group
    # All loci in the group use the SAME merged region scale
    for (row_offset in 0:(end_excel_row - start_excel_row)) {
      excel_row <- start_excel_row + row_offset
      data_row <- merge_info$start + row_offset
      
      if (data_row > nrow(final_data)) break
      
      refinement_type <- final_data$REFINEMENT_TYPE[data_row]
      start_pos <- final_data$START[data_row]
      end_pos <- final_data$END[data_row]
      
      # Determine color based on refinement type
      if (refinement_type == "merged_summary") {
        fill_color <- refinement_colors$merged_row
        
        # Make summary row bold
        row_range <- paste0("A", excel_row, ":", openxlsx2::int2col(ncol(final_data)), excel_row)
        wb$add_font(sheet = "Refined_Loci", dims = row_range, bold = TRUE)
        
        # For summary row, fill entire visualization span
        for (col_idx in vis_start_col:(vis_start_col + n_vis_cols - 1)) {
          cell_ref <- openxlsx2::int2col(col_idx)
          cell_address <- paste0(cell_ref, excel_row)
          wb$add_fill(sheet = "Refined_Loci", dims = cell_address, 
                     color = openxlsx2::wb_color(hex = fill_color))
        }
      } else {
        # For refined loci, use the color for that refinement type
        if (refinement_type %in% names(refinement_colors)) {
          fill_color <- refinement_colors[[refinement_type]]
        } else {
          fill_color <- "FFD3D3D3"  # Default light gray
        }
        
        # Calculate visualization positions relative to THE SAME merged region for ALL rows
        merged_width <- merge_info$merged_width
        
        if (!is.na(start_pos) && !is.na(end_pos) && merged_width > 0) {
          # Scale positions to visualization columns (ALL relative to same merged region)
          rel_start <- (start_pos - merge_info$merged_start) / merged_width
          rel_end <- (end_pos - merge_info$merged_start) / merged_width
          
          # Clamp to [0, 1] range
          rel_start <- max(0, min(1, rel_start))
          rel_end <- max(0, min(1, rel_end))
          
          # Convert to column indices using proper rounding
          scaled_start <- vis_start_col + round(rel_start * n_vis_cols)
          scaled_end <- vis_start_col + round(rel_end * n_vis_cols)
          
          # Only fill if the region is large enough to span at least one column
          # For very small regions (like 1bp in a large merged region), they won't show
          if (scaled_end > scaled_start) {
            # Fill visualization cells
            for (col_idx in scaled_start:min(scaled_end - 1, vis_start_col + n_vis_cols - 1)) {
              cell_ref <- openxlsx2::int2col(col_idx)
              cell_address <- paste0(cell_ref, excel_row)
              
              wb$add_fill(sheet = "Refined_Loci", dims = cell_address, 
                         color = openxlsx2::wb_color(hex = fill_color))
            }
          }
          # If scaled_start == scaled_end, the region is too small to visualize, so skip
        }
      }
    }
  }
  
  # Auto-size data columns
  wb$set_col_widths(sheet = "Refined_Loci", cols = 1:ncol(final_data), widths = "auto")
  
  # Set visualization columns to narrower width (20 columns total)
  vis_cols <- vis_start_col:(vis_start_col + n_vis_cols - 1)
  wb$set_col_widths(sheet = "Refined_Loci", cols = vis_cols, widths = 2)
  
  # Add a legend worksheet
  wb$add_worksheet("Legend")
  
  legend_data <- data.frame(
    REFINEMENT_TYPE = c("source", "merged_summary", "split", "expanded", "aligned", "unchanged"),
    DESCRIPTION = c(
      "Original FUMA locus from source (before merging)",
      "Merged locus boundaries (after merging sources)",
      "Merged locus split into multiple regions based on reported loci",
      "Locus boundaries expanded to encompass reported locus",
      "Locus aligned with reported locus (no boundary changes)",
      "Merged locus unchanged (no overlapping reported loci)"
    ),
    stringsAsFactors = FALSE
  )
  
  wb$add_data(sheet = "Legend", x = legend_data, start_row = 1, start_col = 1)
  
  # Add color to legend
  for (i in seq_len(nrow(legend_data))) {
    type <- legend_data$REFINEMENT_TYPE[i]
    if (type %in% names(refinement_colors)) {
      cell_address <- paste0("A", i + 1)  # +1 for header
      wb$add_fill(sheet = "Legend", dims = cell_address, 
                 color = openxlsx2::wb_color(hex = refinement_colors[[type]]))
    }
  }
  
  wb$set_col_widths(sheet = "Legend", cols = 1:2, widths = "auto")
  
  # Save workbook
  wb$save(output_file)
  
  message(sprintf("Excel file saved to: %s", output_file))
  message(sprintf("  Total refined loci: %d", nrow(final_data)))
  message(sprintf("  Merged loci groups: %d", length(merge_ranges)))
  message(sprintf("  Average splits per merged locus: %.2f", 
                  nrow(refined_loci) / length(unique(refined_loci$MERGED_LOCUS_ID))))
  
  return(invisible(output_file))
}
