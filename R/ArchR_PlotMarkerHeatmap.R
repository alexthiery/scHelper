#' ArchR function: Plots marker heatmap
#'
#' This function takes matrix from seMarker and plots heatmap
#'
#' @param mat data matrix, created by: Log2norm(extract_means_from_se(seMarker))
#' @param pal color palatte for heatmap cells
#' @param labelRows boolean, whether to label rows
#' @param fontSizeRows font size for row labels
#' @param labelCols boolean, whether to label cols
#' @param fontSizeCols font size for col labels
#' @param cluster_columns boolean, whether to cluster cols
#' @param cluster_rows boolean, whether to cluster rows
#' @return Heatmap of marker features
#' @export
ArchR_PlotMarkerHeatmap <- function(mat, pal = viridis::magma(100), 
                           labelRows = FALSE, fontSizeRows = 12,
                           labelCols = TRUE, fontSizeCols = 12,
                           cluster_columns = TRUE, cluster_rows = TRUE) {
  
  # scale each feature independently and add min/max limits
  limits <- c(-2, 2) # could make this user-defined
  mat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), 
               `/`)
  mat[mat > max(limits)] <- max(limits)
  mat[mat < min(limits)] <- min(limits)
  
  # colours - set default if NULL
  if (is.null(pal) == TRUE) {
    pal <- paletteContinuous(set = "solarExtra", n = 100)
  }
  
  # # order rows by eucladian distance
  # dist_mat <- dist(mat, method = 'euclidean')
  # hclust_avg <- hclust(dist_mat, method = 'average')
  # ordered_features <- hclust_avg$labels[c(hclust_avg$order)]
  # mat <- mat[match(ordered_features, rownames(mat)), ]
  # 
  # # order columns by eucladian distance
  # dist_mat <- dist(t(mat), method = 'euclidean')
  # hclust_avg <- hclust(dist_mat, method = 'average')
  # ordered_cell_groups <- hclust_avg$labels[c(hclust_avg$order)]
  # mat <- mat[ , match(ordered_cell_groups, colnames(mat))]
  
  Heatmap(
    matrix = mat,
    col = pal,
    heatmap_legend_param = list(title = "z-scores"),
    #top_annotation = topAnno, 
    # add raster stuff?
    
    #Column Options
    cluster_columns = cluster_columns,
    show_column_names = labelCols,
    column_names_gp = gpar(fontsize = fontSizeCols),
    column_names_max_height = unit(100, "mm"),
    #column_split = colData$stage,
    
    #Row Options
    cluster_rows = cluster_rows,
    show_row_names = labelRows,
    row_names_gp = gpar(fontsize = fontSizeRows)
    #row_split = row_split_params
  )
  
}