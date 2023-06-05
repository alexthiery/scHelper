#' Extract module contents
#'
#' This function extracts the ids of peak or gene modules in a vector so they can be plotted
#'
#' @param modules_df data frame of modules, each row name is module name and each column is a gene or peak ID
#' @param module_names vector of module names to extract, must match rownames of modules_df
#' @return vector of peak or gene IDs
#' @export
ExtractModuleContents <- function(modules_df, module_names){
  if (sum(module_names %in% rownames(modules_df)) < length(module_names)){stop("Module names must match rownames of module df!")}
  ids <- c()
  for (module_name in module_names){
    temp_ids <- sub(' ', '', as.character(as.vector(modules_df[rownames(modules_df) %in% module_name, ])) )
    temp_ids <- temp_ids[nzchar(temp_ids)]
    temp_ids <- temp_ids[!is.na(temp_ids)]
    ids <- c(ids, temp_ids)
  }
  return(ids)
}
