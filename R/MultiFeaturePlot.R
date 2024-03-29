#' Plot FeaturePlots for a list of genes
#'
#' This function plots UMAPs for a series of nine cluster resolutions, alongside a clustree which reveals the relationship between clusters at different resolutions
#'
#' @param seurat_object Seurat object
#' @param gene_list genes to plot DimPlots for
#' @param plot_stage boolean for whether to plot DimPlot for developmental stage
#' @param stage_col seurat_object@meta.data column used to plot developmental stage
#' @param plot_clusters boolean for whether to plot DimPlot for clusters
#' @param cluster_col seurat_object@meta.data column used to plot clusters
#' @param plot_celltype boolean for whether to plot DimPlot for celltype
#' @param celltype_col seurat_object@meta.data column used to plot celltype
#' @param label plot title
#' @return plot generated by grid extra
#' @export
MultiFeaturePlot <- function(seurat_object, gene_list, plot_stage = FALSE, stage_col = "orig.ident", plot_clusters = TRUE, cluster_col = "seurat_clusters", plot_celltype = FALSE, celltype_col = NA, n_col = 4, label = "UMAP plots for GOI on normalised filtered data"){
  
  plots <- lapply(gene_list, function(x) FeaturePlot(seurat_object, features = x))
  
  if(plot_celltype == T){
    plots <- c(list(DimPlot(seurat_object, group.by = celltype_col, label = T) + 
                      ggtitle(paste("Cell Types")) +
                      theme(plot.title = element_text(hjust = 0.5)) +
                      NoLegend()), plots)
  }
  
  if(plot_clusters == T){
    plots <- c(list(DimPlot(seurat_object, group.by = cluster_col) + 
                      ggtitle(paste("Clusters")) +
                      theme(plot.title = element_text(hjust = 0.5))), plots)
  }
  
  if(plot_stage == T){
    plots <- c(list(DimPlot(seurat_object, group.by =  stage_col) + 
                         ggtitle("Developmental Stage") +
                         theme(plot.title = element_text(hjust = 0.5))), plots)
  }
  
  print(gridExtra::grid.arrange(grobs = plots, ncol = n_col, top = textGrob(label = label, gp=gpar(fontsize=20, font = 2))))
}