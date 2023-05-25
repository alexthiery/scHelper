#' Function to create cell type annotation for Complex peak module heatmap
#'
#' @param plot_data Made by PrepPeakModuleHeatmap function
#' @param scHelper_cell_type_colors vector of colors for scHelper_cell_type annotation
#' @return Heatmap annotation for scHelper_cell_type
#' @export
 
CreateCellTypeAnnotation <- function(plot_data, scHelper_cell_type_colors){
  return(
    HeatmapAnnotation(scHelper_cell_type = anno_simple(x = as.character(plot_data$col_ann$scHelper_cell_type),
                                                       col = scHelper_cell_type_colors, height = unit(0.5, "cm")), show_annotation_name = FALSE,
                      labels = anno_mark(at = cumsum(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths) - floor(rle(as.character(plot_data$col_ann$scHelper_cell_type))$lengths/2),
                                         labels = rle(as.character(plot_data$col_ann$scHelper_cell_type))$values,
                                         which = "column", side = 'bottom',
                                         labels_gp = gpar(fontsize = 10), lines_gp = gpar(lwd=2)))
  )
}