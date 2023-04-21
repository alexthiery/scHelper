#' ArchR function: Extract mean matrix from se object
#'
#' This function takes seMarker and extracts the mean matrix for heatmap plotting
#'
#' @param seMarker the se object created `ArchR::getMarkerFeatures`
#' @param Log2norm boolean to log2 normalize the matrix + scale for cluster size
#' @param scaleTo scale factor to multiply by
#' @return means matrix for plotting
#' @export
ArchR_extract_means_from_se <- function(seMarker, Log2norm = TRUE, scaleTo = 10^4) {
  mat <- as.data.frame(SummarizedExperiment::assays(seMarker)[["Mean"]])
  rownames(mat) <- rowData(seMarker)$unique_id

  if (Log2norm) {
    # normalising means for depth of cluster
    mat <- log2(t(t(mat)/colSums(mat)) * scaleTo + 1) 
  }
  
  return(mat)
}