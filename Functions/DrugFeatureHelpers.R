# Helper functions for DrugFeature module

#' Create Drug Comparison Plot
#' @param data Data frame with drug sensitivity data
#' @param comparison_var Name of the variable to compare by
#' @param value_column Name of the column containing drug sensitivity values
#' @param value_label Label for the y-axis
#' @param num_bins Number of bins for continuous variables (default: 4)
#' @param show_groups_boxplot Whether to show boxplot for groups (default: TRUE)
#' @param title Plot title (optional, default: NULL for auto-generated)
#' @return A ggplot object
#' @description This function delegates to the appropriate plotting function from DROMA_R package:
#'   - For continuous variables: plotContinuousComparison (scatter) or plotContinuousGroups (boxplot)
#'   - For categorical variables: plotCategoryComparison (boxplot)
createDrugComparisonPlot <- function(data, comparison_var, value_column, 
                                      value_label = "Drug Sensitivity",
                                      num_bins = 4, show_groups_boxplot = TRUE,
                                      title = NULL) {
  
  # Check required package availability
  if (!requireNamespace("DROMA.R", quietly = TRUE)) {
    stop("Package 'DROMA.R' is required for this function. Please install it first.")
  }
  
  # Check if comparison variable exists
  if (!comparison_var %in% names(data)) {
    stop(paste("Variable", comparison_var, "not found in data"))
  }
  
  # Check if value column exists
  if (!value_column %in% names(data)) {
    stop(paste("Value column", value_column, "not found in data"))
  }
  
  # Determine if comparison variable is numeric or categorical
  if (is.numeric(data[[comparison_var]])) {
    # Continuous variable
    if (show_groups_boxplot) {
      # Use plotContinuousGroups for grouped boxplot
      p <- DROMA.R::plotContinuousGroups(
        data = data,
        cont_column = comparison_var,
        value_column = value_column,
        value_label = value_label,
        num_bins = num_bins,
        title = title
      )
    } else {
      # Use plotContinuousComparison for scatter plot with correlation
      p <- DROMA.R::plotContinuousComparison(
        data = data,
        cont_column = comparison_var,
        value_column = value_column,
        value_label = value_label,
        title = title
      )
    }
  } else {
    # Categorical variable - use plotCategoryComparison
    p <- DROMA.R::plotCategoryComparison(
      data = data,
      category_column = comparison_var,
      value_column = value_column,
      value_label = value_label,
      title = title
    )
  }
  
  return(p)
}

