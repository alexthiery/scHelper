#' ArchR function: Extract IDs for plotting
#'
#' This function extracts the IDs to be plotted, either by cut off or cut off + top n features
#'
#' @param seMarker the se object created `ArchR::getMarkerFeatures`
#' @param cutOff string to filter features, eg by FDR or Log2FC
#' @param top_n boolean to extract top n features
#' @param n if top_n is TRUE, extract top n features
#' @param group_name name of grouping to be used when extracting top n features
#' @return means matrix for plotting
#' @export
ArchR_ExtractIds <- function(seMarker, cutOff = "FDR <= 1 & Log2FC >= 0", top_n = TRUE, n = 10, group_name = "clusters") {
  
  markerList <- getMarkers(seMarker, cutOff = cutOff) # extract features that pass threshold
  
  df <- data.frame() # merged all features into a df
  for (i in 1:length(names(markerList))) {
    print(i)
    df_i <- as.data.frame(markerList[i])
    df <- rbind(df, df_i)
  }
  
  if (top_n == FALSE){
    ids <- df$unique_id
  } else {
    df <- df %>%
      group_by(group_name) %>%
      top_n(n, Log2FC) %>%
      dplyr::arrange(Log2FC, .by_group = TRUE)
    ids <- unique(df$unique_id)
  }
  
  return(ids)
}