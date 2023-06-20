#' Cell Labels Pie Charts
#'
#' This function generates a pie chart for each cell label showing the proportion of cells in each sample
#'
#' @param counts df containing cell labels as columns and samples as rows
#' @param cols vector of colours to use for each cell label
#' @param scale boolean value to specify whether to scale the data by the number of cells in each sample
#' @return ggplot object
#' @export
CellLabelPieCharts <- function(counts, cols, scale = FALSE) {
  
  # remove any columns that have no cells in
  if (0 %in% colSums(counts)){
    counts <- counts[,-(which(colSums(counts)==0))] }
  
  # scale by number of cells in each row
  if (scale == TRUE){
    count_data <- t(apply(counts, 1, function(x) x/sum(x)))
  } else { 
    count_data <- counts }
  
  # make piecharts
  plots <- list()
  for (i in colnames(counts)){
    print(i)
    # calculate totals to add to plot titles
    raw_data <- data.frame(group = rownames(counts), value = (counts[,i]))
    total <- sum(raw_data$value)
    # extract either scaled or raw data for plotting
    data <- data.frame(group = rownames(count_data), value = (count_data[,i]))
    plot <- ggplot(data, aes(x="", y=value, fill=group)) +
      geom_bar(stat="identity", width=1, color="white") +
      coord_polar("y", start=0) +
      theme_void() + scale_fill_manual(values= cols) +
      ggtitle(paste0(i, " (cells: ", total, ")")) +
      theme(legend.position="none", plot.title = element_text(size = 20, hjust = 0.5, vjust = 0))
    plots[[i]] <- plot
  }
  do.call(grid.arrange,plots)
}
