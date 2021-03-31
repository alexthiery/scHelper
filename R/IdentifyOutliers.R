#' Identify oulier clusters
#'
#' This function identifies clusters of cells based on QC metrics of interest. Groups are classed as outliers if the median for a given metric falls outside quartiles 1 and 3.
#'
#' @param dat df containing sample metadata (i.e. seurat@meta.data)
#' @param group_by column to group cells by
#' @param metrics array of metrics to test for outliers
#' @param intersect_metrics boolean value for whether to intersect the outliers for each metric
#' @return array of clusters which are outliers
#' @export
IdentifyOutliers <- function(dat, group_by, metrics, intersect_metrics = FALSE){
  outlier <- list()
  for(metric in metrics){
    min = quantile(dat[[metric]])[[2]]
    max = quantile(dat[[metric]])[[4]]
    
    outlier[[metric]] <- dat %>%
      group_by((!!as.symbol(group_by))) %>%
      summarise(median = median((!!as.symbol(metric)))) %>%
      filter(median > max | median < min) %>%
      pull(!!as.symbol(group_by))
  }
  
  ifelse(intersect_metrics, return(Reduce(intersect, outlier)), return(as.character(unique(unlist(outlier)))))
}