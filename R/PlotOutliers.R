#' PlotOutliers
#' 
#' Plot box plots for elements in seurat_obj@meta.data with quartile cutoffs
#' 
#' @param seurat_obj Seurat object
#' @param y_elements array of column names in seurat_object@meta.data to be plotted
#' @param group_by seurat_object@meta.data column containing cell clusters
#' @param x_lab label for x axis
#' @param quantiles Percent quantiles to plot on bar plots. Can be either NULL, a single element, or an array (length 2) with lower and upper percent quantiles
#' @return list of box plots for each y_element
#' @export
PlotOutliers <- function(seurat_obj, y_elements, group_by = 'seurat_clusters', x_lab = "Cluster ID", quantiles){
  lapply(y_elements, function(x){
    p = BoxPlot(dat = seurat_obj@meta.data, y_col = x, group_by = group_by, y_lab = x, x_lab = x_lab)
    
    if(length(quantiles) > 2){
      stop("quantiles must be either, NULL, a single element or an array of length 2")
      } else {
        if(length(quantiles) > 0) {p = p + geom_hline(yintercept = quantile(seurat_obj@meta.data[[x]], probs = quantiles[1]), linetype = "dashed", color = "red")}
        if(length(quantiles) > 1) {p = p + geom_hline(yintercept = quantile(seurat_obj@meta.data[[x]], probs = quantiles[2]), linetype = "dashed", color = "red")}
        }
    return(p)
    }
  )
}