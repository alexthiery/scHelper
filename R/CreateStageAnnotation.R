#' Function to create stage annotation for Complex peak module heatmap
#'
#' @param plot_data Made by PrepPeakModuleHeatmap function
#' @param scHelper_cell_type_colors vector of colors for stage annotation
#' @return Heatmap annotation for stage
#' @export
 
CreateStageAnnotation <- function(plot_data, stage_colours){
  return(
    HeatmapAnnotation(stage = anno_block(gp = gpar(fill = stage_colours),
                                         labels = levels(plot_data$col_ann$stage),
                                         labels_gp = gpar(col = "white", fontsize = 20, fontface='bold')),
                      simple_anno_size = unit(1, "cm"),
                      annotation_label = "stage", gp = gpar(fontsize = 20))
  )
}