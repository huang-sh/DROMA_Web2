#' Create Plot with Common Axes
#'
#' @param plot A ggplot object
#' @param x_title Title for x-axis
#' @param y_title Title for y-axis
#' @return A ggplot object with updated axes
#' @export
createPlotWithCommonAxes <- function(plot, x_title, y_title) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for this function")
  }
  
  plot + ggplot2::labs(x = x_title, y = y_title) +
    ggplot2::theme(
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      axis.text = ggplot2::element_text(size = 10)
    )
}
