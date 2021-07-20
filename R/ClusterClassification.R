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

ClusterClassification <- function (seurat_obj, group_by = "seurat_clusters", metrics, 
                  quantile, annotation_name = "scHelper_cell_type") 
{
  if (!length(quantile) == 1) {
    stop("quantile must be a single numeric value")
  }
  if (!all(metrics %in% colnames(seurat_obj@meta.data))) {
    missing_metrics <- paste(metrics[!metrics %in% colnames(seurat_obj@meta.data)], 
                             collapse = ", ")
    stop("The following cell types are missing from the seurat metadata: ", 
         missing_metrics, ".\nPlease check that scHelper::AverageGeneModules has beeen ran for the missing cell types.")
  }
  
  tic = TRUE
  metadata <- seurat_obj@meta.data
  cell_type_df = data.frame()
  
  while(tic == TRUE){
    # select clusters and cell type annotations which are higher than percentile threshold
    threshold_subset <- metadata %>%
      dplyr::select(c(!!as.symbol(group_by), metrics)) %>%
      reshape2::melt() %>%
      mutate(variable = as.character(variable)) %>%
      group_by(variable) %>%
      mutate(threshold = quantile(value, probs = 0.9)) %>%
      group_by(!!as.symbol(group_by), variable) %>%
      mutate(median = median(value)) %>%
      filter(median >= threshold) %>%
      select(!value) %>%
      distinct(.keep_all = TRUE) %>%
      group_by(variable, .drop=FALSE) %>%
      add_count(name="tally")
    
    # subset cell types which are assigned to only one cluster (remaining clusters will be re-assigned in loop)
    # if a remaining cluster has two cell-type identities, throw warning and assign to highest normalised median
    unique_cell_type <- threshold_subset %>%
      filter(tally == 1) %>%
      group_by(!!as.symbol(group_by)) %>%
      mutate(scaled_med = median/threshold) %>%
      top_n(1, scaled_med) %>%
      dplyr::select(!!as.symbol(group_by), variable)
    
    tic = ifelse(all(threshold_subset$tally == 1), FALSE, TRUE)
    
    # remove unique cell types from metadata in order to re-assign remaining cells
    metadata <- filter(metadata, !(!!as.symbol(group_by) %in% unique_cell_type[[eval(group_by)]]))
    
    # iteratively add cluster assignments
    cell_type_df <- rbind(cell_type_df, unique_cell_type)
  }
  
  if (nrow(cell_type_df) == 0) {
    stop("No cell types detected!")
  }
  
  # NO LONGER USED
  # duplicated_clusters <- cell_type_df[[group_by]][duplicated(cell_type_df[[group_by]])]
  # if (length(duplicated_clusters) > 1) {
  #   warning("Clusters: ", paste(as.character(duplicated_clusters), 
  #                               collapse = ", "), " are annotated to more than one cell type. These labels will be concatenated in the final annotation.")
  # }
  # else if (length(duplicated_clusters) > 0) {
  #   warning("Cluster ", as.character(duplicated_clusters), 
  #           " is annotated to more than one cell type. These labels will be concatenated in the final annotation.")
  # }
  # cell_type_df <- aggregate(data = cell_type_df, tissue ~ ., FUN = paste, collapse = "/")
  
  # add cell_type to seurat metadata
  seurat_obj@meta.data[[annotation_name]] <- unlist(apply(seurat_obj@meta.data, 
                                                   1, function(x) ifelse(x[[group_by]] %in% cell_type_df[[group_by]], 
                                                                         cell_type_df[cell_type_df[[group_by]] == x[[group_by]], "variable"], NA)))
  return(seurat_obj)
}












