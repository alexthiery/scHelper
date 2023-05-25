#' Function for prep data to plot peak module heatmap using Complex Heatmap
#'
#' @param peak_normalised_matrix matrix of peak normalised data, rows are peaks and columns are cells
#' @param cell_metadata data frame of cell metadata, rows are cells and columns are metadata
#' @param col_order vector of column names to order cells by
#' @param custom_order_column column name to order cells by custom ordering
#' @param custom_order vector of cell names to order cells by custom ordering
#' @param hclust_SEACells logical, whether to hclust SEACells
#' @param hclust_SEACells_within_groups logical, whether to hclust SEACells within SEACell groups
#' @param peak_modules list of peak modules, each element is a vector of peak IDs
#' @param peak_row_annotation logical, whether to add peak row annotation
#' @param scale_data logical, whether to scale data
#' @param log_path path to save logs (hclust dendograms), if NULL no logs will be saved
#' @return plot_data list containing plot data, column annotation and row annotation
#' @export
 
PrepPeakModuleHeatmap <- function (peak_normalised_matrix, cell_metadata, 
                                   col_order, custom_order_column = NULL, custom_order = NULL, 
                                   hclust_SEACells = FALSE, hclust_SEACells_within_groups = TRUE,
                                   peak_modules, peak_row_annotation = TRUE,
                                   scale_data = TRUE,
                                   log_path = NULL)
{
  
  ### Cell-level ordering and annotations ###
  
  # Initiate column anndata
  col_ann <- cell_metadata %>% mutate_if(is.character, as.factor)
  
  # If 'custom_order' is set use this to reorder cells
  if (!is.null(custom_order)) {
    if (!setequal(custom_order, unique(col_ann[[custom_order_column]]))) {
      stop("custom_order factors missing from custom_order_column \n\n")
    }
    col_ann[[custom_order_column]] <- factor(col_ann[[custom_order_column]], levels = custom_order)
    col_ann <- col_ann[order(col_ann[[custom_order_column]]),]
  }
  
  # If 'col_order' is use these columns to order cells
  if (!is.null(col_order)) {
    col_ann <- col_ann[do.call("order", c(col_ann[col_order], list(decreasing = FALSE))), , drop = FALSE]
  }
  
  # Optionally hclust SEACells
  if (hclust_SEACells == TRUE) {
    
    # Hclust SEACells within SEACell groups eg within scHelper_cell_type groups, which is specified by last col_order value
    if (hclust_SEACells_within_groups == TRUE) {
      cell_groups <- split(col_ann, col_ann[[tail(col_order, n=1)]])
      CellGroups_ordered_SEACells <- c()
      for (i in names(cell_groups)) {
        mat <- peak_normalised_matrix[rownames(cell_groups[[i]]), ]
        dist_mat <- dist(mat, method = "euclidean")
        hclust_avg <- hclust(dist_mat, method = "average")
        if (!is.null(log_path)){
          dir.create(log_path, recursive = T)
          png(paste0(log_path, i, '_SEACells_dendogram.png'),  width = 60, height = 40, units = 'cm', res = 400)
          plot(hclust_avg, main = paste0("SEACells dendogram for ", i))
          graphics.off()
        }
        ordered_SEACells <- hclust_avg$labels[c(hclust_avg$order)]
        CellGroups_ordered_SEACells[[i]] <- ordered_SEACells
      }
      col_ann <- col_ann[order(match(rownames(col_ann), unlist(CellGroups_ordered_SEACells))), , drop = FALSE]
      
      # Hclust SEACells across all SEACells, then if there is a col_order secondarily order by the last col_order value
    } else {
      dist_mat <- dist(peak_normalised_matrix, method = "euclidean")
      hclust_avg <- hclust(dist_mat, method = "average")
      if (!is.null(col_order)){ hclust_avg <- with(col_ann, reorder(hclust_avg, as.numeric(col_ann[[tail(col_order, n=1)]])))}
      if (!is.null(log_path)){
        dir.create(log_path, recursive = T)
        png(paste0(log_path, 'SEACells_dendogram.png'), width = 120, height = 40, units = 'cm', res = 400)
        plot(hclust_avg, main = "SEACells dendogram")
        graphics.off()
        }
    ordered_SEACells <- hclust_avg$labels[c(hclust_avg$order)]
    col_ann <- col_ann[order(match(rownames(col_ann), ordered_SEACells)), , drop = FALSE]
  }
  
  ### Peak-level ordering and annotations ###
  
  # Optionally annotate peaks by their modules
  if (peak_row_annotation == TRUE) {
    row_ann <- stack(peak_modules) %>% dplyr::rename(`Peak Modules` = ind) %>%
      column_to_rownames("values")
  } else {
    row_ann <- NA
  }
  
  ### Prepare data for plotting ###
  
  # Order matrix by row and column annotation orders
  plot_data <- t(peak_normalised_matrix)[unlist(peak_modules), rownames(col_ann)]
  
  # Optionally scale
  if (scale_data) {
    cat("Scaling data \n")
    plot_data <- t(scale(t(plot_data)))
    plot_data <- replace(plot_data, plot_data >= 2, 2)
    plot_data <- replace(plot_data, plot_data <= -2, -2)
  }
  
  ### Output plotting data and annotations ###
  
  output <- list(plot_data = plot_data,
                 row_ann = row_ann,
                 col_ann = col_ann)
  return(output)
  
  }
}