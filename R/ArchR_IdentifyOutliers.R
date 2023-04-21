#' ArchR function: Identify oulier clusters
#'
#' This function identifies clusters of cells based on QC metrics of interest. Groups are classed as outliers if the median for a given metric falls outside the provided quantiles.
#'
#' @param ArcgR ArchR object
#' @param group_by column to group cells by
#' @param metrics array of metrics to test for outliers
#' @param intersect_metrics boolean value for whether to intersect the outliers for each metric
#' @param quantiles Percent quantiles to plot on bar plots. Must be an array (length 2) with lower and upper percent quantiles
#' @return array of clusters which are outliers
#' @export
ArchR_IdentifyOutliers <- function(ArchR, group_by = 'clusters', metrics, intersect_metrics = TRUE, quantiles){
  outlier <- list()
  if(!length(quantiles) == 2){
    stop('quantiles must be an array of length == 2')
  }
  for(metric in metrics){
    min = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[1])
    max = quantile(getCellColData(ArchR, select = metrics)[,1], probs = quantiles[2])
    
    outlier[[metric]] <- as.tibble(getCellColData(ArchR)) %>%
      group_by((!!as.symbol(group_by))) %>%
      summarise(median = median((!!as.symbol(metric)))) %>%
      filter(median > max | median < min) %>%
      pull(!!as.symbol(group_by))
  }
  
  if(intersect_metrics){
    if(length(Reduce(intersect, outlier)) == 0){
      cat('No outliers detected!')
    } else {
      return(Reduce(intersect, outlier))
    }
  } else{
    if(length(as.character(unique(unlist(outlier)))) == 0){
      cat('No outliers detected!')
    } else {
      return(as.character(unique(unlist(outlier))))
    }
  }
}