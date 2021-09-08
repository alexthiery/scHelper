#' ClusterClassification
#'
#' This function retrieves cluster classifications using cell module score columns in seurat object. Clusters are classiffied as a given cell type if the median for a given cell type falls outside the provided quantile. 
#'
#' @param seurat_obj Seurat object
#' @param group_by column to group cells by
#' @param cell_types column names for cell type module scores found in seurat_obj@meta.data
#' @param quantile Percent quantile to plot on bar plots. Must be a single value representing the upper percent quantile
#' @param annotation_name Name of column to be annotated with corresponding cell types
#' @param force_assign Boolean for whether to force assign cell identity to clusters which do not pass the filtering threshold
#' @param plot_path Directory to plot log files
#' @param ... Additional options to be passed to grid.arrange
#' @return array of clusters which are outliers
#' @export

ClusterClassification <- function(seurat_obj = seurat_data, group_by = "seurat_clusters", cell_type_markers = cell_type_markers, 
                                  quantile = 0.8, assay = 'RNA', annotation_name = "scHelper_cell_type", force_assign = FALSE, plot_path = "scHelper_log/", ...) 
{
  dir.create(plot_path, recursive = TRUE, showWarnings = FALSE)
  
  if (!length(quantile) == 1) {
    stop("\nquantile must be a single numeric value\n")
  }
  
  iterate = TRUE
  iteration = 1
  metrics = names(cell_type_markers)

  DefaultAssay(seurat_obj) <- assay

  temp_seurat <- seurat_obj

  cell_type_df = data.frame()
  
  while(iterate == TRUE){
    # Calculate average module expression for contamination gene list
    temp_seurat <- AverageGeneModules(seurat_obj = temp_seurat, gene_list = cell_type_markers)
    
    png(paste0(plot_path, "classification_iteration_", iteration, ".png"), width = 40, height = 30, units = "cm", res = 200)
    seurat_obj@misc[[annotation_name]] <- append(seurat_obj@misc[[annotation_name]], list( iteration = grid.arrange(grobs = PlotOutliers(seurat_obj = temp_seurat, 
                                                                                                                                         y_elements = metrics, group_by = group_by, 
                                                                                                                                         quantiles = quantile), ...)))
    graphics.off()
    
    # select clusters and cell type annotations which are higher than percentile threshold
    celltype_df = temp_seurat@meta.data %>%
      dplyr::select(c(!!as.symbol(group_by), all_of(metrics))) %>%
      reshape2::melt() %>%
      mutate(variable = as.character(variable)) 
    
    # Calculate threshold module scores for each gene module
    celltype_df = celltype_df %>%
      group_by(variable, .drop = FALSE) %>%
      mutate(threshold = quantile(value, probs = quantile))
    
    # Add median gene module scores for each cell cluster
    celltype_df = celltype_df %>%
      group_by(!!as.symbol(group_by), variable, .drop = FALSE) %>%
      mutate(median = median(value)) %>%
      group_by(variable, .drop = FALSE) %>%
      mutate(median_percentile = ecdf(value)(median))
    
    # Drop individual cell gene module scores and remove duplicate rows
    celltype_df = celltype_df %>%
      select(!value) %>%
      distinct(.keep_all = TRUE)
    
    # Select top classification for each cluster
    celltype_df = celltype_df %>%
      group_by(!!as.symbol(group_by)) %>% 
      top_n(1, median_percentile)
    
    # Select clusters which pass the threshold subset
    threshold_subset = filter(celltype_df, median >= threshold)
    
    # Select remaining (unclassified) clusters
    remaining_clusters <- filter(celltype_df, !(!!as.symbol(group_by)) %in% threshold_subset[[group_by]])
    
    
    # iteratively add cluster assignments
    if(nrow(threshold_subset) != 0){
      cell_type_df <- rbind(cell_type_df, threshold_subset)
      
      if(nrow(remaining_clusters) == 0){
        iterate = FALSE
        cat("\nAll clusters classified!\n")
        
      } else if (nrow(remaining_clusters) <= 1) {
        iterate = FALSE
        
        if(!force_assign){
          remaining_clusters$variable <- NA
          cat("\nThe following clusters did not pass the percentile threshold. Clusters will be assigned N/A. \nClusters:", 
              as.character(unique(remaining_clusters[[group_by]])), '\n')
        } else {
          cat("\nThe following clusters did not pass the percentile threshold. Clusters will be assigned to the cell type with the highest percentile. \nClusters:", 
              as.character(unique(remaining_clusters[[group_by]])), '\n')
        }
        cell_type_df <- rbind(cell_type_df, remaining_clusters)
        
      } else{
        iteration = iteration + 1
        # Remove assigned cell types from temp_seurat object for next round of iteration
        temp_seurat <- subset(temp_seurat, subset = !!as.symbol(group_by) %in% threshold_subset[[group_by]], invert=TRUE)
      }
      
    } else if(nrow(threshold_subset) == 0){
      iterate = FALSE
      
      if(!force_assign){
        remaining_clusters$variable <- NA
        cat("\nThe following clusters did not pass the percentile threshold. Clusters will be assigned N/A. \nClusters:", 
            as.character(unique(remaining_clusters[[group_by]])), '\n')
      } else {
        cat("\nThe following clusters did not pass the percentile threshold. Clusters will be assigned to the cell type with the highest percentile. \nClusters:", 
            as.character(unique(remaining_clusters[[group_by]])), '\n')
      }
      cell_type_df <- rbind(cell_type_df, remaining_clusters)
    }
  }
    
    # check if any clusters are annotated to multiple cell states (annotations share same percentile)
    if(any(cell_type_df %>% count(seurat_clusters) %>% pull(n) > 1)){
      duplicate_clusters <- cell_type_df %>%
        count(seurat_clusters) %>%
        filter(!n == 1) %>%
        pull(seurat_clusters) %>%
        as.character()
      
      cat("\nFollowing clusters are equally annotated to multiple cell states. Annotated cell states will be concatenated. \nClusters:", duplicate_clusters, "\n")
      
      cell_type_df = cell_type_df %>% group_by(seurat_clusters) %>% 
        mutate(variable = paste0(variable, collapse = "/")) %>% 
        distinct(!!as.symbol(group_by), variable, median_percentile)

    }
    
    # add cell_type to seurat metadata
    seurat_obj@meta.data[[annotation_name]] <- unlist(apply(seurat_obj@meta.data, 
                                                            1, function(x) ifelse(x[[group_by]] %in% cell_type_df[[group_by]], 
                                                                                  cell_type_df[cell_type_df[[group_by]] == x[[group_by]], "variable"], NA)))
    return(seurat_obj)
  }

