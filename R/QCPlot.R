#' Generate a QC plot from Seurat metadata
#'
#' This function generates a series of box plots as quality control for a Seurat object
#'
#' @param seurat_obj Seurat object
#' @param group_by seurat_obj@meta.data column to group cells by on x axis
#' @param y_elements seurat_obj@meta.data column to plot QC metrics for
#' @param stage seurat_obj@meta.data column to stack cell count bars by
#' @param y_lab array of labels for y axis reflecting y_elements
#' @param x_lab label for x axis
#' @param plot_quantiles boolean for whether to plot horizontal quantile cutoff lines
#' @return QC plots
#' @export
QCPlot <- function (seurat_obj,
                    group_by = "seurat_clusters",
                    y_elements = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
                    stage = "orig.ident",
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
  
  bar_data <- seurat_obj@meta.data %>%
    rownames_to_column('cell_name') %>%
    dplyr::select(c(cell_name, !!sym(group_by), !!sym(stage))) %>%
    group_by(!!sym(group_by)) %>%
    count(!!sym(stage), .drop = FALSE)
  
  bar_plot = ggplot(bar_data, aes(fill=!!sym(stage), y=n, x=!!sym(group_by))) + 
    geom_bar(position="stack", stat="identity") +
    ylab('Cell Count') + xlab('Clusters') + theme(legend.position = c(0.5, 0.9),
                                                  legend.direction="horizontal",
                                                  panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                  panel.background = element_blank(), axis.line = element_line(colour = "black"), 
                                                  strip.background = element_rect(colour = "white", fill = "white"), 
                                                  strip.text = element_text(size = 18), axis.text = element_text(size = 16), 
                                                  axis.title = element_text(size = 18))
  
  plots$d = bar_plot
  
  grid.arrange(grobs = plots, ncol = 2)
}