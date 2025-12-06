#' Plot merged FUMA loci showing original loci and merged boundaries
#'
#' @param merged_loci Data frame. Output from merge_fuma_loci()
#' @param locus_id Integer or vector. Which merged locus/loci to plot. Default: all (max 10)
#' @param max_plot Integer. Maximum number of loci to plot. Default: 10
#' @return ggplot object or list of ggplot objects
#' @export
plot_merged_loci <- function(merged_loci, locus_id = NULL, max_plot = 10) {
  
  # Determine which loci to plot
  if (is.null(locus_id)) {
    locus_id <- seq_len(min(max_plot, nrow(merged_loci)))
  } else {
    # Validate locus_id
    invalid <- !locus_id %in% merged_loci$MERGED_LOCUS_ID
    if (any(invalid)) {
      stop("Invalid locus_id(s): ", paste(locus_id[invalid], collapse = ", "))
    }
  }
  
  # Filter to selected loci
  plot_data <- merged_loci[merged_loci$MERGED_LOCUS_ID %in% locus_id, ]
  
  # Create plots
  plot_list <- lapply(seq_len(nrow(plot_data)), function(i) {
    row <- plot_data[i, ]
    
    # Parse original loci
    orig_loci_str <- strsplit(row$ORIGINAL_LOCI, ";")[[1]]
    
    # Extract information from each original locus
    orig_data <- do.call(rbind, lapply(seq_along(orig_loci_str), function(j) {
      locus_str <- orig_loci_str[j]
      
      # Parse: locus_id(source:start:end)
      pattern <- "^(.+)\\((.+):(\\d+):(\\d+)\\)$"
      
      if (grepl(pattern, locus_str)) {
        locus_name <- sub(pattern, "\\1", locus_str)
        source <- sub(pattern, "\\2", locus_str)
        start <- as.numeric(sub(pattern, "\\3", locus_str))
        end <- as.numeric(sub(pattern, "\\4", locus_str))
        
        data.frame(
          level = j,
          locus_name = locus_name,
          source = source,
          start = start,
          end = end,
          label = paste0(source, "\n", locus_name),
          stringsAsFactors = FALSE
        )
      } else {
        NULL
      }
    }))
    
    if (is.null(orig_data) || nrow(orig_data) == 0) {
      warning("Could not parse original loci for merged locus ", row$MERGED_LOCUS_ID)
      return(NULL)
    }
    
    # Create merged locus data
    merged_data <- data.frame(
      level = 0,
      locus_name = "MERGED",
      source = "MERGED",
      start = row$START,
      end = row$END,
      label = sprintf("Merged\nchr%s:%s-%s", 
                     row$CHR, 
                     format(row$START, big.mark = ","),
                     format(row$END, big.mark = ",")),
      stringsAsFactors = FALSE
    )
    
    # Combine data
    all_data <- rbind(merged_data, orig_data)
    all_data$midpoint <- (all_data$start + all_data$end) / 2
    
    # Get unique sources for coloring
    unique_sources <- unique(orig_data$source)
    n_sources <- length(unique_sources)
    
    # Color palette
    if (n_sources <= 8) {
      colors <- RColorBrewer::brewer.pal(max(3, n_sources), "Set2")
    } else {
      colors <- rainbow(n_sources)
    }
    
    source_colors <- setNames(colors[seq_len(n_sources)], unique_sources)
    source_colors["MERGED"] <- "gray30"
    
    # Create plot
    p <- ggplot2::ggplot(all_data, ggplot2::aes(xmin = start, xmax = end, 
                                                  ymin = level - 0.4, ymax = level + 0.4,
                                                  fill = source)) +
      ggplot2::geom_rect(color = "black", size = 0.5) +
      ggplot2::geom_text(ggplot2::aes(x = midpoint, y = level, label = label),
                        size = 3, vjust = 0.5) +
      ggplot2::scale_fill_manual(values = source_colors) +
      ggplot2::scale_y_continuous(breaks = all_data$level,
                                   labels = c("MERGED", paste0("Original ", seq_len(nrow(orig_data)))),
                                   expand = ggplot2::expansion(mult = 0.1)) +
      ggplot2::scale_x_continuous(labels = scales::comma,
                                   expand = ggplot2::expansion(mult = 0.02)) +
      ggplot2::labs(
        title = sprintf("Merged Locus %d: chr%s:%s-%s (%s bp)",
                       row$MERGED_LOCUS_ID,
                       row$CHR,
                       format(row$START, big.mark = ","),
                       format(row$END, big.mark = ","),
                       format(row$WIDTH, big.mark = ",")),
        subtitle = sprintf("Combines %d original loci from %d source(s)",
                          row$N_ORIGINAL_LOCI,
                          length(unique_sources)),
        x = "Genomic Position (bp)",
        y = "",
        fill = "Source"
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(size = 10),
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        plot.subtitle = ggplot2::element_text(size = 10, color = "gray40")
      )
    
    return(p)
  })
  
  # Remove NULL plots
  plot_list <- Filter(Negate(is.null), plot_list)
  
  # Return single plot or list
  if (length(plot_list) == 1) {
    return(plot_list[[1]])
  } else {
    return(plot_list)
  }
}


#' Plot all merged loci in a multi-panel figure
#'
#' @param merged_loci Data frame. Output from merge_fuma_loci()
#' @param ncol Integer. Number of columns in the grid. Default: 2
#' @param max_plot Integer. Maximum number of loci to plot. Default: 10
#' @return Combined ggplot object
#' @export
plot_all_merged_loci <- function(merged_loci, ncol = 2, max_plot = 10) {
  
  n_plot <- min(max_plot, nrow(merged_loci))
  locus_ids <- merged_loci$MERGED_LOCUS_ID[seq_len(n_plot)]
  
  plot_list <- plot_merged_loci(merged_loci, locus_id = locus_ids)
  
  if (!is.list(plot_list)) {
    plot_list <- list(plot_list)
  }
  
  # Combine plots
  combined <- gridExtra::grid.arrange(grobs = plot_list, ncol = ncol)
  
  return(combined)
}


#' Export merged loci plots to folder (one plot per locus)
#'
#' @param merged_loci Data frame. Output from merge_fuma_loci()
#' @param output_folder Character. Output folder path (required)
#' @param width Numeric. Plot width in inches. Default: 10
#' @param height Numeric. Plot height in inches. Default: 4
#' @param format Character. File format (e.g., "pdf", "png", "svg"). Default: "pdf"
#' @param dpi Numeric. Resolution for raster formats. Default: 300
#' @export
export_merged_plots <- function(merged_loci, 
                                output_folder,
                                width = 10, 
                                height = 4,
                                format = "pdf",
                                dpi = 300) {
  
  if (missing(output_folder)) {
    stop("output_folder is required. Please provide a folder path (e.g., 'plots/merged_loci')")
  }
  
  # Create output folder if it doesn't exist
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
    message(sprintf("Created output folder: %s", output_folder))
  }
  
  message(sprintf("Creating individual plots for %d merged loci...", nrow(merged_loci)))
  
  # Create plots for each locus
  saved_files <- character()
  
  for (i in seq_len(nrow(merged_loci))) {
    locus_id <- merged_loci$MERGED_LOCUS_ID[i]
    chr <- merged_loci$CHR[i]
    
    # Create plot
    p <- plot_merged_loci(merged_loci, locus_id = locus_id)
    
    # Generate filename
    filename <- sprintf("merged_locus_%03d_chr%s.%s", locus_id, chr, format)
    filepath <- file.path(output_folder, filename)
    
    # Save plot
    ggplot2::ggsave(
      filename = filepath,
      plot = p,
      width = width,
      height = height,
      units = "in",
      dpi = dpi
    )
    
    saved_files <- c(saved_files, filepath)
    
    if (i %% 10 == 0 || i == nrow(merged_loci)) {
      message(sprintf("  Progress: %d/%d plots saved", i, nrow(merged_loci)))
    }
  }
  
  message(sprintf("\nAll %d plots saved to: %s", length(saved_files), output_folder))
  
  return(invisible(saved_files))
}
