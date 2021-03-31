#' Generate a QC plot from Seurat metadata
#'
#' This function generates a series of box plots as quality control for a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param group_by seurat_obj@meta.data column to group cells by on x axis
#' @param y_elements seurat_obj@meta.data column to plot QC metrics for
#' @param y_lab array of labels for y axis reflecting y_elements
#' @param x_lab label for x axis
#' @param plot_quantiles boolean for whether to plot horizontal quantile cutoff lines
#' @return QC plots
#' @export
QCPlot <- function (seurat_obj,
          group_by = "seurat_clusters",
          y_elements = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
          y_lab = c("UMI Count", "Gene Count", "% MT"),
          x_lab = "Cluster ID",
          plot_quantiles = FALSE) 
{
  plots <- list()
  
  if(plot_quantiles){
    plots$a = BoxPlot(dat = seurat_obj@meta.data, y_col = y_elements[1], 
                      group_by = group_by, y_lab = y_lab[1], x_lab = x_lab) +
      geom_hline(yintercept=quantile(seurat_obj@meta.data[[y_elements[1]]])[[2]], linetype="dashed", color = "red") +
      geom_hline(yintercept=quantile(seurat_obj@meta.data[[y_elements[1]]])[[4]], linetype="dashed", color = "red")
    
    plots$b = BoxPlot(dat = seurat_obj@meta.data, y_col = y_elements[2], 
                      group_by = group_by, y_lab = y_lab[2], x_lab = x_lab) +
      geom_hline(yintercept=quantile(seurat_obj@meta.data[[y_elements[2]]])[[2]], linetype="dashed", color = "red") +
      geom_hline(yintercept=quantile(seurat_obj@meta.data[[y_elements[2]]])[[4]], linetype="dashed", color = "red")
    
    plots$c = BoxPlot(dat = seurat_obj@meta.data, y_col = y_elements[3], 
                      group_by = group_by, y_lab = y_lab[3], x_lab = x_lab) +
      geom_hline(yintercept=quantile(seurat_obj@meta.data[[y_elements[3]]])[[2]], linetype="dashed", color = "red") +
      geom_hline(yintercept=quantile(seurat_obj@meta.data[[y_elements[3]]])[[4]], linetype="dashed", color = "red")
    
  } else {
    plots$a = BoxPlot(dat = seurat_obj@meta.data, y_col = y_elements[1], 
                      group_by = group_by, y_lab = y_lab[1], x_lab = x_lab)
    plots$b = BoxPlot(dat = seurat_obj@meta.data, y_col = y_elements[2], 
                      group_by = group_by, y_lab = y_lab[2], x_lab = x_lab)
    plots$c = BoxPlot(dat = seurat_obj@meta.data, y_col = y_elements[3], 
                      group_by = group_by, y_lab = y_lab[3], x_lab = x_lab)
  }

  grid.arrange(grobs = plots, ncol = 3)
}