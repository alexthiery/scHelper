#' Generate a QC plot from Seurat metadata
#'
#' This function generates a series of box plots as quality control for a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param group_by seurat_obj@meta.data column to group cells by on x axis
#' @param y_elements seurat_obj@meta.data column to plot QC metrics for
#' @param y_lab array of labels for y axis reflecting y_elements
#' @param x_lab label for x axis
#' @return QC plots
#' @export
QCPlot <- function(seurat_obj, group_by = "seurat_clusters", y_elements = c("nCount_RNA", "nFeature_RNA", "percent.mt") ,
                    y_lab = c("UMI Count", "Gene Count", "% MT"), x_lab = "Cluster ID"){
  plots <- list()
  plots$a = box.plot(dat = seurat_obj@meta.data, y_col = y_elements[1], group_by = group_by, y_lab = y_lab[1], x_lab = x_lab)
  plots$b = box.plot(dat = seurat_obj@meta.data, y_col = y_elements[2], group_by = group_by, y_lab = y_lab[2], x_lab = x_lab)
  plots$c = box.plot(dat = seurat_obj@meta.data, y_col = y_elements[3], group_by = group_by, y_lab = y_lab[3], x_lab = x_lab)
  
  grid.arrange(grobs = plots, ncol = 3)
}