#' Box Plot
#'
#' This function generates a box plot from a single cell metadata dataframe
#'
#' @param data df containing sample metadata (i.e. seurat@meta.data)
#' @param y_col dat column to place on y axis
#' @param group_by dat column to group cells by on x axis
#' @param y_lab label for y axis
#' @param x_lab label for x axis
#' @return ggplot object
#' @export
BoxPlot <- function(data, y_col, group_by, y_lab, x_lab){
  ggplot(dat, aes(x = .data[[group_by]], y = .data[[y_col]], fill = .data[[group_by]])) +
    geom_boxplot() +
    xlab(x_lab) +
    ylab(y_lab) +
    theme(
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      strip.background = element_rect(colour = "white", fill = "white"),
      strip.text = element_text(size = 18),
      axis.text = element_text(size = 16),
      axis.title = element_text(size = 18))
}
