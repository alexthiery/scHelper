
#' Elbow Cutoff
#'
#' This function unbiasedly calculates where the principal components start to elbow.
#' 
#' First we take the larger value of the point where the principal components only contribute 5% of standard deviation and the point where the principal components cumulatively contribute 90% of the standard deviation.
#' Next we take the point where the percent change in variation between the consecutive PCs is less than 0.1%.
#' 
#' The smaller out of these two values is determined at the elbow cutoff
#'
#' @param seurat_object Seurat object
#' @param return either 'pc_cutoff' or 'plot' based on whether the user would like to return a ggplot object or the principle component cutoff point
#' @return ggplot object OR principle component cutoff
#' @export
ElbowCutoff <- function(seurat_object, return = 'pc_cutoff'){
  # Determine percent of variation associated with each PC
  pct <- seurat_object[["pca"]]@stdev / sum(seurat_object[["pca"]]@stdev) * 100
  
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  
  # Determine the difference between variation of PC and subsequent PC
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  
  
  if(return == 'pc_cutoff'){
    return(pcs)
  } else if (return == 'plot'){
    # Create a dataframe with values
    plot_df <- data.frame(pct = pct, 
                          cumu = cumu, 
                          rank = 1:length(pct))
    
    # Elbow plot to visualize 
    return(ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
      geom_text() + 
      geom_vline(xintercept = 90, color = "grey") + 
      geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
      ylab('Percent SD') +
      xlab('Cumulative SD') +
      theme_bw())
  } else {
    stop('return must be one of "pc_cutoff" or "plot"')
  }
}

