#' @export
ClustStagePlot <- function(data, cluster.col = "seurat_clusters", stage.col = "orig.ident"){
  clust_plot <- DimPlot(data, group.by = cluster.col) + 
    ggtitle(paste("Clusters")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  stage_plot <- DimPlot(data, group.by =  stage.col) + 
    ggtitle("Developmental Stage") +
    theme(plot.title = element_text(hjust = 0.5))
  
  gridExtra::grid.arrange(stage_plot, clust_plot, nrow = 1)
}
