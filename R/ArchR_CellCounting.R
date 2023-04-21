#' ArchR function: Plot table of cell counts from different groups
#'
#' This function counts how many cells from each cluster/sample are assigned the same label/cluster
#'
#' @param ArchR ArchR object
#' @param group1 first column to split cells by
#' @param group2 second column to split cells by
#' @param print_table boolean to print table or return it
#' @param scHelper_cell_type_order order in which to return cell types
#' @return table of cell counts, or plots table directly
#' @export
ArchR_cell_counting <- function(ArchR = ArchR, group1 = "clusters", group2 = "stage", print_table = TRUE, scHelper_cell_type_order = scHelper_cell_type_order) {
  
  group1_data <- getCellColData(ArchR, select = group1)[,1]
  group1_cell_counts <- as.data.frame(table(group1_data))
  colnames(group1_cell_counts) <- c("ID", "Total_count")
  
  group2_cell_counts <- data.frame()
  group2_data <- getCellColData(ArchR, select = group2)[,1]
  data_group1 <- getCellColData(ArchR, select = group1)[,1]
  for (i in unique(group1_data)) {
    cells <- ArchR$cellNames[BiocGenerics::which(data_group1 == i)]
    if (length(cells) > 1){
      ArchR_subset <- ArchR[cells, ]
      data_group2 <- getCellColData(ArchR_subset, select = group2)[,1]
      group2_cell_counts_i <- as.data.frame(table(data_group2)) %>%
        tidyr::pivot_wider(names_from = data_group2, values_from = Freq) %>% 
        tibble::add_column(ID = !!i)
      group2_cell_counts <- rbind.fill(group2_cell_counts, group2_cell_counts_i)
    }
  }
  
  cell_counts <- merge(group1_cell_counts, group2_cell_counts) %>%
    column_to_rownames(., var = "ID")
  totals <- cell_counts$Total_count
  cell_counts <- cell_counts[, -1]
  cell_counts[is.na(cell_counts)] <- 0
  
  # Order rows and columns - group 1 = order rows
  if (group1 %in% c("clusters", "stage")) {
    cell_counts <- cell_counts[ mixedsort(rownames(cell_counts)) , ] }
  if (group1 == "scHelper_cell_type_old") {
    if (is.null(scHelper_cell_type_order)){
      print("scHelper_cell_type_order not specified!")
    } else {
      order <- scHelper_cell_type_order[scHelper_cell_type_order %in% rownames(cell_counts)]
      cell_counts <- cell_counts[ order , ] }
  }
  
  # Order rows and columns - group 2 = order columns
  if (group2 %in% c("clusters", "stage")){
    cell_counts <- cell_counts[ , mixedsort(colnames(cell_counts)) ]
  }
  if (group2 == "scHelper_cell_type_old") {
    if (is.null(scHelper_cell_type_order)){
      print("scHelper_cell_type_order not specified!")
    } else {
      order <- scHelper_cell_type_order[scHelper_cell_type_order %in% colnames(cell_counts)]
      cell_counts <- cell_counts[ , order ] }
  }
  
  # either return table or print it
  if (print_table == FALSE){
    return(cell_counts)
  } else {
    cell_counts <- cell_counts %>% mutate(., Total = totals)
    grid.arrange(tableGrob(cell_counts, theme = ttheme_minimal()))
  }
}
