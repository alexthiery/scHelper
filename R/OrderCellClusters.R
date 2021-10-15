#' Automatic ordering of cell clusters based on cluster composition
#'
#' This function orders cell clusters by a column in the seurat metadata. i.e. sort seurat clusters by their developmental stage, with clusters shared between stages ordered in between those stages
#'
#' @param seurat_object Seurat object
#' @param col_to_sort seurat_object@meta.data column to be re-arranged (i.e. clusters). Only rearrange cluster if at least 10% of cells in cluster are from a different identity
#' @param sort_by array containing factor levels passed to col_to_sort (i.e. developmental stage)
#' @return array of cluster IDs automatically ordered based on their composition
#' @export
OrderCellClusters = function(seurat_object, col_to_sort, sort_by){

  # if there are any character cols in metadata convert to factor
  seurat_object@meta.data[sapply(seurat_object@meta.data, is.character)] <- lapply(seurat_object@meta.data[sapply(seurat_object@meta.data, is.character)], as.factor)
  
  dat <- seurat_object@meta.data %>%
    group_by(!!as.symbol(sort_by)) %>%
    count(!!as.symbol(col_to_sort)) %>%
    arrange(!!as.symbol(col_to_sort)) %>%
    group_by(!!as.symbol(col_to_sort)) %>%
    mutate(n = n/sum(n))
  
  
  top2 <- dat %>%
    top_n(n = 2) %>%
    arrange(!!as.symbol(sort_by), desc(n)) %>%
    filter(n > 0.1)
  
  
  top1 <- dat %>%
    top_n(n = 1) %>%
    arrange(!!as.symbol(sort_by), desc(n)) %>%
    distinct(!!as.symbol(col_to_sort), .keep_all = TRUE)
  
  for(i in 1:nrow(top1)){
    
    tempnext = unique(as.integer(top2[[sort_by]][as.character(top2[[sort_by]]) == as.character(top1[[sort_by]])[i]]) + 1)
    tempprev = unique(as.integer(top2[[sort_by]][as.character(top2[[sort_by]]) == as.character(top1[[sort_by]])[i]]) - 1)
    
    if(sum(as.integer(top2[[sort_by]]) == tempnext) == 0 & sum(as.integer(top2[[sort_by]]) == tempprev) == 0){
      next
    } else if (sum(as.integer(top2[[sort_by]]) == tempnext) != 0){
      nextfactor = top2[as.integer(top2[[sort_by]]) == tempnext,]
      
      if(top1[[col_to_sort]][i] %in% nextfactor[[2]]){
        top1[i,"n"] <- top1[i,"n"] - 1
      } else {next}
      
      
    } else if (sum(as.integer(top2[[sort_by]]) == tempprev) != 0){
      prevfactor = top2[as.integer(top2[[sort_by]]) == tempprev,]
      
      if(top1[[col_to_sort]][i] %in% prevfactor[[2]]){
        top1[i,"n"] <- (1- top1[i,"n"]) + 1
      } else {next}
      
    } else {stop("error - factor present in more than two levels")}
  }
  
  top1$n = top1$n + (rev(order(levels(top1[[sort_by]])))[top1[[sort_by]]] * 2)
  
  top1 = top1 %>%
    arrange(!!as.symbol(sort_by), desc(n))
  
  return(top1[[col_to_sort]])
}



