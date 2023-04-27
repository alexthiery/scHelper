#' SEACells function: Summarises gene score counts across SEACells
#'
#' This function takes ATAC gene score matrix and adds up gene score counts across SEACells
#'
#' @param matrix exported gene score matrix from SEACells
#' @param dictionary dataframe with two columns: 'index' (single cell ID) and 'SEACell' (SEACell ID)
#' @return dataframe with summarised gene counts which can be used to make a seurat object of SEACells
#' @export
SEACells_SummariseGeneScoreMatrix <- function(matrix, dictionary){
  
  # turn gene score matrix into numeric dataframe
  df <- as.data.frame(t(as.data.frame(matrix)))
  df <- df[-1, ] # remove row with gene names
  df2 <- mutate_all(df, function(x) as.numeric(as.character(x)))
  
  # add dictionary of cell ids to metacell ids to dataframe
  df3 <- df2 %>% mutate(category = dictionary$SEACell)
  
  # aggregate df based on metacell ids
  df_summarised <- aggregate(. ~ category, df3, sum)
  
  # format df so can be added to seurat object
  df_summarised <- t(column_to_rownames(df_summarised, var = "category"))
  
  return(df_summarised)
}