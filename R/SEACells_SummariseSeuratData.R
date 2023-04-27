#' SEACells function: Summarises data from seurat object across SEACells
#'
#' This function takes a seurat object and summarises the data from data_slot by SEACells grouping
#'
#' @param seurat seurat object
#' @param data_slot which assay in seurat object to summarise, default is "counts"
#' @param category slot in seurat metadata by which to group cells for summarising, default is "SEACell"
#' @return dataframe with summarised data across groupings
#' @export
SEACells_SummariseSeuratData <- function(seurat, data_slot = "counts", category = "SEACell"){
  
  # extract data into dataframe
  df <- GetAssayData(object = seurat, slot = data_slot)
  df <- as.data.frame(t(as.data.frame(df)))
  
  # convert cell ids to category ids
  category_ids <- select(seurat@meta.data, category)[,1]
  df <- df %>% mutate(category = category_ids)
  
  # aggregate df based on category
  df_summarised <- aggregate(. ~ category, df, sum)
  
  # format df so can be added back to seurat object
  df_summarised <- t(column_to_rownames(df_summarised, var = "category"))
  
  return(df_summarised)
}