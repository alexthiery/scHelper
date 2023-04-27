#' ArchR function: Add unique id to se object
#'
#' This function takes seMarker and adds a unique_id column
#'
#' @param seMarker the se object created `ArchR::getMarkerFeatures`
#' @param ArchR ArchR object
#' @param matrix_type which matrix to use ("GeneScoreMatrix" or "PeakMatrix")
#' @return seMarker with unique_id column
#' @export
ArchRAddUniqueIdsToSe <- function(seMarker, ArchR, matrix_type) {
  
  if (matrix_type == "PeakMatrix") {
    tmp_peaks = data.frame(ArchR@peakSet)
    tmp_diff_peaks = data.frame(rowData(seMarker))
    diff_peaks_join_peakset = left_join(tmp_diff_peaks, tmp_peaks, 
                                        by = c("seqnames" = "seqnames", "start" = "start", "end" = "end"))
    diff_peaks_join_peakset$gene_name = paste(diff_peaks_join_peakset$nearestGene, diff_peaks_join_peakset$distToTSS,sep="_")
    diff_peaks_join_peakset$unique_id = paste0(diff_peaks_join_peakset$seqnames, ":", diff_peaks_join_peakset$start, "-", diff_peaks_join_peakset$end)
    rowData(seMarker) = diff_peaks_join_peakset
  }
  
  if (matrix_type == "GeneScoreMatrix") {
    rowData <- as.data.frame(rowData(seMarker))
    
    duplicated_gene_names <- rowData$name[duplicated(rowData$name)]
    duplicated_genes <- rowData[which(rowData$name %in% duplicated_gene_names), ]
    duplicated_genes <- duplicated_genes %>% group_by(name) %>% 
      arrange(name) %>% mutate(unique_id = paste0(name, "-", rowid(name)))
    join_df = left_join(rowData, duplicated_genes,
                        by = c("seqnames" = "seqnames", "start" = "start", "end" = "end", "name" = "name", "idx" = "idx", "strand" = "strand"))
    
    join_df <- join_df %>% mutate(unique_id = coalesce(unique_id, name))
    
    rowData(seMarker) = join_df
  }
  
  return(seMarker)
}