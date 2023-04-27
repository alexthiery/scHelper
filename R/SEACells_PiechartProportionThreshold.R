#' SEACells function: Purity of metacells plotted on piechart
#'
#' Function to plot piechart of how many metacells pass threshold for proportion of cells coming from one label
#'
#' @param prop_table dataframe of proportions of metacells in each category, output from scHelper::calculate_metacell_frequencies
#' @param threshold threshold for proportion of cells coming from one label
#' @return pie chart showing proportions of cells that pass threshold purity
#' @export
SEACells_PiechartProportionThreshold <- function(prop_table, threshold = 0.5){
  
  # filter cells to only include those that pass threshold
  passed_cells <- prop_table %>% filter(prop > threshold)
  # number of cells that pass threshold:
  passed_cells <- length(unique(passed_cells$Metacell))
  # number of cells that didn't pass threshold:
  failed_cells <- length(unique(prop_table$Metacell)) - passed_cells
  # plot piechart
  slices <- c(passed_cells, failed_cells)
  lbls <- c(paste0("Passed: ", passed_cells), paste0("Didn't pass: ", failed_cells))
  
  return(pie(slices, labels = lbls, main = paste0("Number of cells that passed threshold (", threshold, ") of label proportions")))
  
}