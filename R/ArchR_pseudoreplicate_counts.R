#' ArchR function: Plot table of cell counts in pseudoreplicates
#'
#' This function prints a table of how many cells are in each pseudoreplicate and how many samples and clusters are they made of
#'
#' @param ArchR ArchR object
#' @param pseudo_replicates object created by Archr::addGroupCoverages
#' @param group_by ArchR metadata to group by (e.g. "Sample")
#' @return table of cell counts for each pseudoreplicate across groups
#' @export
ArchR_PseudoreplicateCounts <- function(ArchR = ArchR, pseudo_replicates, group_by = "Sample") {
  
  unlisted <- unlist(pseudo_replicates, recursive=FALSE)
  print(paste0("Number of pseudoreplicates: ", length(unlisted)))
  group_cell_counts <- data.frame()
  
  # iterate through each pseudoreplicate
  for (i in c(1:length(unlisted))) {
    #print(paste0("Pseudoreplicate number: ", i))
    group_name <- names(unlisted[i])
    cell_IDs <- unlisted[i][[1]]
    ArchR_pseudo_replicate <- ArchR[cell_IDs, ]
    
    # add up contributions of each group to pseudoreplicates
    group_cell_count <- as.data.frame(table(getCellColData(ArchR_pseudo_replicate, select = group_by))) %>%
      pivot_wider(names_from = Var1, values_from = Freq) %>% 
      add_column(pseudo_replicate_ID = group_name)
    group_cell_counts <- rbind.fill(group_cell_counts, group_cell_count)
    
  }
  
  # format table
  group_cell_counts[is.na(group_cell_counts)] <- 0
  group_cell_counts <- group_cell_counts %>% relocate(pseudo_replicate_ID)
  
  grid.arrange(tableGrob(group_cell_counts, rows=NULL, theme = ttheme_minimal()))
  
}