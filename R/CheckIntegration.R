#' Check Integration
#'
#' This function plots UMAPs after integration, allowing you to see that cells from the same population but sequenced as part of separate batches, overlap after integration.
#'
#' @param seurat_object Seurat object
#' @param group_by seurat_object@meta.data column to group cells by - this should correspond to the variable which identifies shared populations across batches
#' @param nrow number of rows to pass to grid.arrange
#' @param shuffle randomly plot cells in order to prevent hiding certain populations when overplotting
#' @param seed set seed for plotting consistency when shuffle is TRUE
#' @return plot generated by grid arrange
#' @export
CheckIntegration <- function(seurat_object, group_by = 'orig.ident', pt.size = 0.2, xlim = c(-15,15), ylim = c(-15,15), nrow = 1, shuffle = TRUE, seed = 123){
  plots = list()
  plots[['all cells']] <- DimPlot(seurat_object, group_by = group_by, pt.size = 0.2, shuffle = shuffle, seed = seed) + xlim(xlim) + ylim(ylim) + labs(title = 'all cells')
  
  # set colours to the same as when plotting all cells
  colours <- factor(seurat_object[[group_by, drop = T]])
  
  stages = unique(gsub('_.*', '', seurat_object[[group_by, drop=T]]))
  
  for(i in stages){
    dat <- filter(seurat_object@meta.data, grepl(pattern = i, orig.ident))
    if(length(unique(dat[[group_by]])) > 1){
      plots[[i]] <- DimPlot(seurat_object, cells = rownames(dat), group_by = group_by, pt.size = 0.2, shuffle = shuffle, seed = seed,
                            cols = setNames(ggplotColours(n = length(levels(colours))), levels(colours))) + xlim(xlim) + ylim(ylim) + labs(title = i)
    }
  }
  do.call("grid.arrange", c(plots, nrow=nrow))
}