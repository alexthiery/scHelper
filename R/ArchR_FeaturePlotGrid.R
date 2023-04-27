#' ArchR function: Plots grid of featureplots
#'
#' This function makes a grid of featureplots
#'
#' @param ArchR ArchR object
#' @param matrix which matrix to plot (default = "GeneScoreMatrix")
#' @param feature_list which features to plot
#' @return grid of featureplots
#' @export
ArchR_FeaturePlotGrid <- function(ArchR, matrix = "GeneScoreMatrix", feature_list) {
  p <- plotEmbedding(ArchR, colorBy = matrix, name = feature_list, 
                     plotAs = "points", size = 1.8, baseSize = 0, labelSize = 8, legendSize = 10)
  p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6.5) +
      theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
      theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
      )
  })
  do.call(cowplot::plot_grid, c(list(ncol = 4),p2))
}