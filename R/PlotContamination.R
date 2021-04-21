#' PlotContamination
#' 
#' Calculate average expression for groups of genes across Seurat clusters and plot BoxPlots.
#' 
#' @param seurat_obj Seurat object
#' @param gene_list Named list of genes passed to seurat::AddModuleScore
#' @param group_by seurat_object@meta.data column containing cell clusters
#' @param x_lab label for x axis
#' @param quantiles Percent quantile to plot on bar plots. Can be either a single element or an array (length 2) with lower and upper percent quantiles
#' @param ... Extra arguments passed to grid.arrange
#' @return plot output from grid.arrange
#' @export
PlotContamination <- function(seurat_obj, gene_list, group_by = 'seurat_clusters', x_lab = "Cluster ID", quantiles = c(0.05, 0.95), ...){
  if(!all(names(gene_list) %in% colnames(seurat_obj@meta.data))){
    seurat_obj <- AverageGeneModules(seurat_obj = seurat_obj, gene_list = gene_list)
  }
  grid.arrange(grobs = PlotOutliers(seurat_obj = seurat_obj, y_elements = names(gene_list), group_by = group_by, x_lab = x_lab, quantiles = quantiles), ...)
}
