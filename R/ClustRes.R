#' Plot multiple cluster resolutions
#'
#' This function plots UMAPs for a series of nine cluster resolutions, alongside a clustree which reveals the relationship between clusters at different resolutions
#'
#' @param seurat_object Seurat object
#' @param by value by which resolution is incremented
#' @param starting_res lowest resolution to cluster cells by
#' @param prefix prefix used to select metadata columns with cluster identities following clustering
#' @return plot generated by grid extra
#' @export
ClustRes <- function(seurat_object, by = 0.1, starting_res = 0, prefix = "RNA_snn_res."){
  plots <- list()
  resolutions <- c(seq(starting_res, starting_res+9*by, by=by))
  if(length(seurat_object@reductions) == 0){
    stop("Carry out dimensionality reduction (PCA) before clustering")
  }
  seurat_object@meta.data <- seurat_object@meta.data[,!grepl(prefix, colnames(seurat_object@meta.data))]
  seurat_object <- FindClusters(seurat_object, resolution = resolutions, verbose = F)
  plots[["clustree"]] <- clustree(seurat_object@meta.data, prefix = prefix)
  for(res in resolutions[2:length(resolutions)]){
    plots[[paste(res)]] <- DimPlot(seurat_object, group.by =  paste0(prefix, res)) +
      ggtitle(paste("resolution = ", res))
  }
  lay <- rbind(c(1,1,1,2,3,4),
               c(1,1,1,5,6,7),
               c(1,1,1,8,9,10))
  plots2 <- gridExtra::arrangeGrob(grobs = plots, layout_matrix = lay)
  return(gridExtra::grid.arrange(plots2))
}
