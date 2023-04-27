#' SEACells function: Frequency of metacells in each category
#'
#' Function to take seurat object and a category and make a freq table of how frequently metacells found in each category
#'
#' @param seurat seurat object
#' @param metacell_slot slot in seurat metadata by which to group cells for summarising, default is "SEACell"
#' @param category slot in seurat metadata by which to count metacell frequencies, default is "stage"
#' @param calc_proportions boolean as to whether to calculate proportions of metacells in each category
#' @return dataframe of either frequency of metacells in each category or proportions of metacells in each category
#' @export
SEACells_MetacellFrequencies <- function(seurat, metacell_slot = "SEACell", category = "stage", calc_proportions = FALSE){
  
  df <- data.frame(FetchData(object = seurat, vars = c(metacell_slot, category)))
  colnames(df) <- c("Metacell", "Category")
  
  df <- df %>%
    group_by(Metacell, Category) %>% 
    dplyr::summarize(count = n()) %>%
    mutate(Metacell = str_split(Metacell, "-", simplify = TRUE)[ , 2]) %>%
    mutate(Metacell = as.numeric(Metacell)) %>% 
    arrange(Metacell) %>% 
    mutate(count = as.numeric(count))
  
  print(paste0("Number of metacells: ", length(unique(freq$Metacell))))
  print(paste0("Number of categories: ", length(unique(freq$Category))))

  if (calc_proportions){

    # calculate total cell counts per metacell
    totals_df <- aggregate(count ~ Metacell, data=df, FUN=sum)
    totals <- totals_df$count
    print(length(totals))
  
    # calculate proportions per metacell
    prop_table <- df %>%
      mutate(prop = count/totals[Metacell+1])

    df <- prop_table
  }

  return(df)
  
}