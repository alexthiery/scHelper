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

ClusterClassification <- function(seurat_obj, group_by = "seurat_clusters", metrics, quantile, annotation_name = "scHelper_cell_type") 
{
  if (!length(quantile) == 1) {
    stop("quantile must be a single numeric value")
  }
  
  if (!all(metrics %in% colnames(seurat_obj@meta.data))) {
    missing_metrics <- paste(metrics[!metrics %in% colnames(seurat_obj@meta.data)], collapse = ", ")
    stop("The following cell types are missing from the seurat metadata: ", missing_metrics, ".\nPlease check that scHelper::AverageGeneModules has beeen ran for the missing cell types.")
  }
  
  df <- data.frame()
  for (metric in metrics) {
    max = quantile(seurat_obj@meta.data[[metric]], probs = quantile)
    df <- seurat_obj@meta.data %>% group_by((!!as.symbol(group_by))) %>% 
      summarise(median = median((!!as.symbol(metric)))) %>% 
      filter(median > max) %>% mutate(tissue = metric) %>%
      dplyr::select(!!as.symbol(group_by), tissue) %>% rbind(df)
  }
  
  if (nrow(df) == 0) {
    stop("No cell types detected!")
  }
  
  # check for duplicated clusters
  duplicated_clusters <- df[[group_by]][duplicated(df[[group_by]])]
  
  if (length(duplicated_clusters) > 1) {
    warning("Clusters: ", paste(as.character(duplicated_clusters), collapse = ", "), " are annotated to more than one cell type. These labels will be concatenated in the final annotation.")
  } else if (length(duplicated_clusters) > 0){
    warning("Cluster ", as.character(duplicated_clusters), " is annotated to more than one cell type. These labels will be concatenated in the final annotation.")
  }
  
  # concatenate tissue labels for clusters with ambiguous annotations
  df <- aggregate(data=df, tissue~.,FUN=paste,collapse="/")
  
  # add column to metadata with resulting cell type annotations
  seurat_obj@meta.data[[annotation_name]] <- apply(seurat_obj@meta.data, 1, function(x) ifelse(x[[group_by]] %in% df[[group_by]],
                                                                                                             df[df[[group_by]] == x[[group_by]], 'tissue'],
                                                                                                             NA))
  
  return(seurat_obj)
}
