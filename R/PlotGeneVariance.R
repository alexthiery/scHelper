#' PlotGeneVariance
#'
#' This function plots gene variance as a function of average expression for groups of cells.
#'
#' @param seurat_obj Seurat object
#' @param group_by column to group cells by
#' @param assay assay containing expression data. Default set to 'RNA'.
#' @return plot
#' @export

PlotGeneVariance <- function(seurat_obj = seurat_data, group_by = 'seurat_clusters', assay = 'RNA'){
  DefaultAssay(seurat_data) <- "RNA"
  seurat_split <- SplitObject(seurat_data, split.by = group_by)
  
  plot_data <- data.frame()
  for(i in names(seurat_split)){
    assay <- GetAssayData(object = seurat_split[[i]])
    gene_variance <- apply(assay, 1, var)
    median_expression <- apply(assay, 1, median)
    plot_data <- rbind(plot_data, data.frame(gene_variance, median_expression, group = i))
  }
  
  
  ggplot(plot_data, aes(x = median_expression, y = gene_variance, colour = group)) +
    geom_smooth(method = 'gam', se = FALSE) +
    ylab('gene variance') +
    xlab('median gene expression')
}
