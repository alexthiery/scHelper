#' Generate a QC plot from Seurat metadata
#'
#' This function generates a series of box plots as quality control for a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param y_elements seurat_obj@meta.data column to plot QC metrics for
#' @param group_by seurat_obj@meta.data column to group cells by on x axis
#' @param stage seurat_obj@meta.data column to stack cell count bars by
#' @param y_lab array of labels for y axis reflecting y_elements
#' @param x_lab label for x axis
#' @param quantiles Percent quantiles to plot on bar plots. Can be either NULL, a single element, or an array (length 2) with lower and upper percent
#' @param ... Extra arguments to be passed to grid.arrange
#' @return plot output from grid.arrange
#' @export
QCPlot <- function(seurat_obj,
                   y_elements = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                   group_by = 'seurat_clusters',
                   y_lab = c("UMI Count", "Gene Count", "% MT"),
                   x_lab = "Cluster ID",
                   stage = "orig.ident",
                   quantiles = c(0.25, 0.75),
                   ...){
  
  p = PlotOutliers(seurat_obj = seurat_obj, y_elements = y_elements, group_by = group_by, x_lab = x_lab, quantiles = quantiles)

  bar_data <- seurat_obj@meta.data %>%
    rownames_to_column('cell_name') %>%
    dplyr::select(c(cell_name, !!sym(group_by), !!sym(stage))) %>%
    group_by(!!sym(group_by)) %>%
    count(!!sym(stage), .drop = FALSE)

  p[['bar_plot']] <- ggplot(bar_data, aes(fill=!!sym(stage), y=n, x=!!sym(group_by))) +
    geom_bar(position="stack", stat="identity") +
    ylab('Cell Count') + xlab('Clusters') + theme(legend.position = c(0.5, 0.9),
                                                  legend.direction="horizontal",
                                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                  panel.background = element_blank(), axis.line = element_line(colour = "black"),
                                                  legend.background = element_rect(fill =  alpha("white", 0)),
                                                  strip.text = element_text(size = 18), axis.text = element_text(size = 16),
                                                  axis.title = element_text(size = 18))
  
  grid.arrange(grobs = p, ...)
}