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

ClusterClassification <- function (seurat_obj, group_by = "seurat_clusters", cell_type_markers, 
                                   quantile, annotation_name = "scHelper_cell_type", plot_path = "scHelper_log/", fast = FALSE, ...) 
{
  dir.create(plot_path, showWarnings = FALSE)
  
  if (!length(quantile) == 1) {
    stop("quantile must be a single numeric value")
  }
  # if (!all(metrics %in% colnames(seurat_obj@meta.data))) {
  #   missing_metrics <- paste(metrics[!metrics %in% colnames(seurat_obj@meta.data)], 
  #                            collapse = ", ")
  #   stop("The following cell types are missing from the seurat metadata: ", 
  #        missing_metrics, ".\nPlease check that scHelper::AverageGeneModules has beeen ran for the missing cell types.")
  # }
  # 
  iterate = TRUE
  iteration = 1
  metrics = names(cell_type_markers)
  temp_seurat <- seurat_obj
  cell_type_df = data.frame()
  
  if(fast == TRUE){
    # Calculate average module expression for contamination gene list
    temp_seurat <- AverageGeneModules(seurat_obj = temp_seurat, gene_list = cell_type_markers)
  }
  
  while(iterate == TRUE){
    if(fast == FALSE){
      # Calculate average module expression for contamination gene list
      temp_seurat <- AverageGeneModules(seurat_obj = temp_seurat, gene_list = cell_type_markers)
    }
    
    # select clusters and cell type annotations which are higher than percentile threshold
    threshold_subset <- temp_seurat@meta.data %>%
      dplyr::select(c(!!as.symbol(group_by), all_of(metrics))) %>%
      reshape2::melt() %>%
      mutate(variable = as.character(variable)) %>%
      group_by(variable) %>%
      mutate(threshold = quantile(value, probs = quantile)) %>%
      group_by(!!as.symbol(group_by), variable) %>%
      mutate(median = median(value)) %>%
      filter(median >= threshold) %>%
      select(!value) %>%
      distinct(.keep_all = TRUE) %>%
      group_by(variable, .drop=FALSE) %>%
      add_count(name="tally") %>%
      group_by(!!as.symbol(group_by), .drop=FALSE) %>%
      add_count(name="tally")
    
    
    temp_cell_type_df <- threshold_subset %>%
      filter(tally == 1) %>%
      dplyr::select(!!as.symbol(group_by), variable)
    
    # subset cell types which are assigned to only one cluster (remaining clusters will be re-assigned in loop)
    # if a remaining cluster has two cell-type identities, throw warning and assign to highest normalised median
    temp_cell_type_df <- threshold_subset %>%
      filter(tally != 1) %>%
      group_by(!!as.symbol(group_by)) %>%
      mutate(scaled_med = median/threshold) %>%
      top_n(1, scaled_med) %>%
      dplyr::select(!!as.symbol(group_by), variable) %>%
      rbind(temp_cell_type_df)
    
    
    png(paste0(plot_path, "classification_iteration_", iteration, ".png"), width = 40, height = 30, units = "cm", res = 200)
    seurat_obj@misc[[annotation_name]] <- append(seurat_obj@misc[[annotation_name]], list( iteration = grid.arrange(grobs = PlotOutliers(seurat_obj = temp_seurat, 
                                                                                                                                         y_elements = metrics, group_by = group_by, 
                                                                                                                                         quantiles = quantile), ...)))
    graphics.off()
    
    # remove assigned cells from seurat object in order to re-assign remaining cells
    temp_seurat <- subset(temp_seurat, subset = !!as.symbol(group_by) %in% temp_cell_type_df[[group_by]], invert=TRUE)
    
    iteration = iteration + 1
    iterate = ifelse(nrow(temp_cell_type_df) == 0, FALSE, TRUE)
    
    # iteratively add cluster assignments
    cell_type_df <- rbind(cell_type_df, temp_cell_type_df)
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