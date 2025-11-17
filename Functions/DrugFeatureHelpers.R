# Helper functions for DrugFeature module

#' Create Drug Comparison Plot
#' @param data Data frame with drug sensitivity data
#' @param comparison_var Name of the variable to compare by
#' @param value_column Name of the column containing drug sensitivity values
#' @param value_label Label for the y-axis
#' @param num_bins Number of bins for continuous variables (default: 4)
#' @param show_groups_boxplot Whether to show boxplot for groups (default: TRUE)
#' @return A ggplot object
createDrugComparisonPlot <- function(data, comparison_var, value_column, 
                                      value_label = "Drug Sensitivity",
                                      num_bins = 4, show_groups_boxplot = TRUE) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("Package 'ggpubr' is required for this function")
  }
  
  # Create a copy to avoid modifying original data
  plot_data <- data
  
  # Check if comparison variable is numeric
  if (is.numeric(plot_data[[comparison_var]])) {
    # Bin continuous variables
    plot_data$group <- cut(plot_data[[comparison_var]], 
                           breaks = num_bins, 
                           include.lowest = TRUE,
                           labels = paste0("Group", 1:num_bins))
    
    x_var <- "group"
    x_label <- paste(comparison_var, "(binned)")
  } else {
    # Use categorical variable as-is
    x_var <- comparison_var
    x_label <- comparison_var
    plot_data$group <- as.factor(plot_data[[comparison_var]])
  }
  
  # Remove NA values in the grouping variable
  plot_data <- plot_data[!is.na(plot_data$group), ]
  
  if (nrow(plot_data) == 0) {
    return(ggplot2::ggplot() + 
             ggplot2::annotate("text", x = 0.5, y = 0.5, 
                              label = "No data available after filtering") + 
             ggplot2::theme_void())
  }
  
  # Create base plot
  if (show_groups_boxplot) {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "group", y = value_column)) +
      ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, fill = "lightblue") +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 2)
  } else {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes_string(x = "group", y = value_column)) +
      ggplot2::geom_violin(alpha = 0.7, fill = "lightblue") +
      ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 2)
  }
  
  # Add statistical comparisons
  p <- p +
    ggpubr::stat_compare_means(method = "kruskal.test", 
                               label.y = max(plot_data[[value_column]], na.rm = TRUE) * 1.1) +
    ggplot2::labs(x = x_label, y = value_label, 
                 title = paste("Drug Sensitivity by", comparison_var)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 14)
    )
  
  return(p)
}

