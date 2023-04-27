#' ArchR function: Subset ArchR
#'
#' This function subsets ArchR specified groups
#'
#' @param ArchR ArchR object
#' @param meta_col1 first column to split cells by
#' @param meta_col2 second column to split cells by
#' @param groups1 first groups to subset by
#' @param groups2 second groups to subset by
#' @param invert1 invert first groups
#' @param invert2 invert second groups
#' @param invert invert all groups
#' @return subsetted ArchR object
#' @export
ArchR_Subset <- function(ArchR, meta_col1, meta_col2, groups1, groups2, invert1, invert2, invert = FALSE){
  
  print(paste0("Filtering on first meta_col: ", meta_col1))
  if (invert1 == FALSE){
    idxPass_1 <- which(ArchR@cellColData[,meta_col1] %in% groups1)
  } else { idxPass_1 <- which(!(ArchR@cellColData[,meta_col1] %in% groups1)) }
  idxPass <- idxPass_1
  
  if (is.null(meta_col2) == FALSE){
    print(paste0("Filtering on second meta_col: ", meta_col2))
    if (invert2 == FALSE){
      idxPass_2 <- which(ArchR@cellColData[,meta_col2] %in% groups2) 
    } else { idxPass_2 <- which(!(ArchR@cellColData[,meta_col2] %in% groups2)) }
    idxPass <- idxPass_1[(idxPass_1 %in% idxPass_2)] # take intersect of pass 1 and pass 2
  }
  
  cellsPass <- ArchR$cellNames[idxPass]
  ArchR_subset <- ArchR[cellsPass, ]
}
