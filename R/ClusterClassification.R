#' ClusterClassification
#'
#' This function retrieves cluster classifications using cell module score columns in seurat object. Clusters are classiffied as a given cell type if the median for a given cell type falls outside the provided quantile. 
#'
#' @param seurat_obj Seurat object
#' @param group_by column to group cells by
#' @param cell_types column names for cell type module scores found in seurat_obj@meta.data
#' @param quantile Percent quantile to plot on bar plots. Must be a single value representing the upper percent quantile
#' @param annotation_name Name of column to be annotated with corresponding cell types
#' @return array of clusters which are outliers
#' @export

ClusterClassification <- function(seurat_obj = seurat_data, group_by = "seurat_clusters", cell_type_markers = cell_type_markers, 
                                  quantile = 0.8, annotation_name = "scHelper_cell_type", plot_path = "scHelper_log/", ...) 
{
  dir.create(plot_path, showWarnings = FALSE)
  
  if (!length(quantile) == 1) {
    stop("quantile must be a single numeric value")
  }

  iterate = TRUE
  iteration = 1
  metrics = names(cell_type_markers)
  temp_seurat <- seurat_obj
  cell_type_df = data.frame()
  
  while(iterate == TRUE){
    # Calculate average module expression for contamination gene list
    temp_seurat <- AverageGeneModules(seurat_obj = temp_seurat, gene_list = cell_type_markers)
    
    # select clusters and cell type annotations which are higher than percentile threshold
    celltype_df <- temp_seurat@meta.data %>%
      dplyr::select(c(!!as.symbol(group_by), all_of(metrics))) %>%
      reshape2::melt() %>% 
      mutate(variable = as.character(variable)) %>% 
      # Set threshold based on predefined percentile
      group_by(variable, .drop = FALSE) %>%
      mutate(threshold = quantile(value, probs = quantile)) %>%
      group_by(!!as.symbol(group_by), variable, .drop = FALSE) %>%
      mutate(median = median(value)) %>%
      # Calculate the percentile for the median value of each cluster
      group_by(variable, .drop = FALSE) %>%
      mutate(median_percentile = ecdf(value)(median)) %>%
      select(!value) %>%
      distinct(.keep_all = TRUE)
    
    # select clusters and cell type annotations which are higher than percentile threshold
    threshold_subset <- celltype_df %>%
      # Filter clusters where the median is lower than the threshold value
      filter(median >= threshold) %>%
      # If cluster passes the threshold - assign to the cell-type with highest median percentile
      group_by(!!as.symbol(group_by)) %>%
      top_n(1, median_percentile) %>%
      dplyr::select(!!as.symbol(group_by), variable)

    png(paste0(plot_path, "classification_iteration_", iteration, ".png"), width = 40, height = 30, units = "cm", res = 200)
    seurat_obj@misc[[annotation_name]] <- append(seurat_obj@misc[[annotation_name]], list( iteration = grid.arrange(grobs = PlotOutliers(seurat_obj = temp_seurat, 
                                                                                                                                         y_elements = metrics, group_by = group_by, 
                                                                                                                                         quantiles = quantile), ...)))
    graphics.off()

    # if no cells are assigned, assign remaining clusters to highest percentile
    if(nrow(threshold_subset) == 0){
      iterate = FALSE
      
      # Filter clusters which dont pass threshold
      below_threshold <- celltype_df %>%
        filter(median < threshold) %>%
        group_by(!!as.symbol(group_by)) %>%
        top_n(1, median_percentile) %>%
        dplyr::select(!!as.symbol(group_by), variable)
        
      cat('The following clusters did not pass the percentile threshold. Clusters will be assigned to the cell type with the highest percentile. \nClusters:', as.character(unique(below_threshold[[group_by]])))
      
      threshold_subset <- rbind(threshold_subset, below_threshold)
      
    } else {
      # remove assigned cells from seurat object in order to re-assign remaining cells
      temp_seurat <- subset(temp_seurat, subset = !!as.symbol(group_by) %in% threshold_subset[[group_by]], invert=TRUE)
      
      # If only one cluster remains - assign as cannot calculate percentile with only one group
      remaining_clusters <- filter(celltype_df, !(!!as.symbol(group_by)) %in% threshold_subset[[group_by]])
      
      if(length(unique(remaining_clusters[[group_by]])) == 1){
        iterate = FALSE
        cat('The following clusters did not pass the percentile threshold. Clusters will be assigned to the cell type with the highest percentile. \nClusters:', as.character(unique(remaining_clusters[[group_by]])))
        
        threshold_subset <- remaining_clusters %>%
          ungroup() %>%
          top_n(1, median_percentile) %>%
          dplyr::select(!!as.symbol(group_by), variable) %>%
          rbind(threshold_subset)
      }
    }
    
    # iteratively add cluster assignments
    cell_type_df <- rbind(cell_type_df, threshold_subset)
    iteration = iteration + 1
  }
  
  if (nrow(cell_type_df) == 0) {
    stop("No cell types detected!")
  }
  
  # add cell_type to seurat metadata
  seurat_obj@meta.data[[annotation_name]] <- unlist(apply(seurat_obj@meta.data, 
                                                          1, function(x) ifelse(x[[group_by]] %in% cell_type_df[[group_by]], 
                                                                                cell_type_df[cell_type_df[[group_by]] == x[[group_by]], "variable"], NA)))
  return(seurat_obj)
}


