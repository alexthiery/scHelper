#' ArchR function: Split by stage
#'
#' This function splits ArchR object by stage
#'
#' @param ArchR ArchR object
#' @return list of ArchR objects, split by stage
#' @export
split_ArchR_by_stage <- function(ArchR){
  split_names <- unique(ArchR$stage)
  split_ArchR <- c()
  for (i in split_names){
    idxSample <- BiocGenerics::which(ArchR$stage %in% i)
    cellsSample <- ArchR$cellNames[idxSample]
    split_ArchR <- c(split_ArchR, ArchR[cellsSample, ])
  }
  return(split_ArchR)
}