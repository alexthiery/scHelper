#' Function for plotting gene module pheatmap from Seurat object
#'
#' col_ann_order can be specified to change the order in which the column annotations appear on the heatmap
#' @param seurat_obj Seurat object
#' @param metadata a string or array of colnames in seurat_obj@meta.data used to annotate the cells.
#' @param col_order cell order will be prioritised by the order of elements in the array. By default cells are ordered by the first element of the metadata variable.
#' @param custom_order array specifying cell group ordering. If set then custom_order_column must also be set.
#' @param custom_order_column metadata column from which custom_order is derived.
#' @param assay seurat assay to use for plotting.
#' @param slot seurat slot to use for plotting.
#' @param gene_modules list of gene modules to plot.
#' @param selected_genes optional array of genes to filter gene modules. Modules with no selected_genes present will be removed.
#' @param hide_annotation array of elements to hide from column annotation.
#' @param gaps_row boolean specifying whether to add whitespace between gene modules.
#' @param gaps_col metadata column by which to add whitespace.
#' @param gm_row_annotation boolean specifying whether to add gene module row annotations.
#' @param cell_subset array of cells to subset for plotting. If specified, data will be re-scaled.
#' @param use_seurat_colours boolean specifying whether to use default seurat colours. If FALSE, colours from colour_scheme will be used instead.
#' @param colour_scheme array or colour schemes from RColourBrewer to be used for annotating columns.
#' @param col_ann_order can be specified to change the order of the column annotations.
#' @param show_colnames default parameter passed to pheatmap.
#' @param show_rownames default parameter passed to pheatmap.
#' @param cluster_rows default parameter passed to pheatmap.
#' @param cluster_cols default parameter passed to pheatmap.
#' @param order_genes boolean specifying whether to order genes within gms based on hclustering
#' @param annotation_names_row default parameter passed to pheatmap.
#' @param ... Extra arguments to be passed to pheatmap.
#' @return plot output from pheatmap.
#' @export
GeneModulePheatmap <- function (seurat_obj, metadata, col_order = metadata[1], custom_order = NULL, 
                     custom_order_column = NULL, assay = "RNA", slot = "scale.data", 
                     gene_modules, selected_genes = NULL, hide_annotation = NULL, 
                     gaps_row = TRUE, gaps_col = NULL, gm_row_annotation = TRUE, 
                     cell_subset = NULL, use_seurat_colours = TRUE, 
                     colour_scheme = c("PRGn", "RdYlBu", "Greys"), col_ann_order = rev(metadata),
                     show_colnames = FALSE, show_rownames = TRUE, cluster_rows = FALSE, cluster_cols = FALSE,
                     order_genes = TRUE, annotation_names_row = FALSE, ...) 
{
  # Subset cells if necessary
  if (!is.null(cell_subset)) {
    seurat_obj <- subset(seurat_obj, cells = cell_subset)
  }
  
  ### Column formatting ###
  
  # Select metadata columns and convert to factor if character
  seurat_obj@meta.data <- seurat_obj@meta.data[,metadata, drop=FALSE] %>% mutate_if(is.character, as.factor)
  # Drop all unused levels
  seurat_obj@meta.data[] <- lapply(seurat_obj@meta.data, function(x) if(is.factor(x)) factor(x) else x)
  
  # Set col_ann
  col_ann <- seurat_obj@meta.data
  
  # Order cells based on col_order
  col_ann <- col_ann[do.call("order", c(col_ann[col_order], list(decreasing = FALSE))), , drop = FALSE]
  
  # Order cells based on custom_order
  if (!is.null(custom_order)) {
    if (is.null(custom_order_column)) {
      "custom_order column must be specified \n"
    }
    if (!setequal(custom_order, unique(col_ann[[custom_order_column]]))) {
      stop("custom_order factors missing from custom_order_column \n\n")
    }
    levels(col_ann[[custom_order_column]]) <- custom_order
    col_ann <- col_ann[order(col_ann[[custom_order_column]]), , drop = FALSE]
  }
  
  # Insert whitespace between columns
  if (!is.null(gaps_col)) {
    ifelse(class(gaps_col) != "character", stop("gaps_col must be a metadata column name"), gaps_col <- cumsum(rle(as.vector(col_ann[[gaps_col]]))[["lengths"]]))
  }
  
  # Hide column annotations
  if (!is.null(hide_annotation)) {
    col_ann[, hide_annotation] <- NULL
  }
  
  # Use seurat colours by default -> otherwise use colour scheme set in colour_scheme
  if (use_seurat_colours == FALSE) {
    ann_colours <- lapply(1:ncol(col_ann), function(x) setNames(colorRampPalette(brewer.pal(9,colour_scheme[x])[2:9])(length(unique(col_ann[, x]))), unique(col_ann[, x])))
  } else {
    ann_colours <- lapply(colnames(col_ann), function(x) {
      temp <- setNames(ggPlotColours(n = length(levels(col_ann[, x]))), levels(col_ann[, x]))
      temp[match(levels(col_ann[[x]]), names(temp))]
    })
  }
  
  names(ann_colours) <- colnames(col_ann)
  
  ### Row formatting ###
  
  # Subset gene modules based on selected genes and rename if modules have no names
  if (!is.null(selected_genes)) {
    selected_GM <- SubsetGeneModules(gm = gene_modules, selected_genes = selected_genes)
  } else {
    if (is.null(names(gene_modules))) {
      names(gene_modules) <- paste0("GM: ", 1:length(gene_modules))
    }
    selected_GM <- gene_modules
  }
  
  # Set row annotations
  if (gm_row_annotation == TRUE) {
    row_ann <- stack(selected_GM) %>% rename("Gene Modules" = ind) %>% column_to_rownames('values')
  } else {
    row_ann <- NA
  }
  
  # Insert whitespace between rows
  if (gaps_row == TRUE) {
    row_ann <- droplevels(row_ann)
    gaps_row = cumsum(summary(row_ann[["Gene Modules"]], maxsum = max(lengths(lapply(row_ann, unique)))))
  } else {
    gaps_row = NULL
  }
  
  # Set row annotation colours
  ann_colours[["Gene Modules"]] <- setNames(colorRampPalette(brewer.pal(9, "Paired"))(length(unique(row_ann[['Gene Modules']]))), unique(row_ann[['Gene Modules']]))
  
  # Order genes within gene modules
  if (order_genes == TRUE) {
    GMs_ordered_genes <- c()
    for (i in names(selected_GM)){
      dist_mat <- dist(seurat_obj@assays$RNA@scale.data[selected_GM[[i]], ], method = 'euclidean')
      hclust_avg <- hclust(dist_mat, method = 'average')
      ordered_genes <- hclust_avg$labels[c(hclust_avg$order)]
      GMs_ordered_genes[[i]] <- ordered_genes}
    selected_GM <- GMs_ordered_genes
  }
  
  ### Prepare data ###
  plot_data <- t(as.matrix(x = GetAssayData(object = seurat_obj, assay = assay, slot = slot)[unlist(selected_GM), rownames(col_ann), drop = FALSE]))
  
  # Scale data if cells have been subset
  if (!is.null(cell_subset)) {
    cat("rescaling data as cells have been subset \n")
    plot_data <- t(scale(t(plot_data)))
  }

  plot_data <- replace(plot_data, plot_data >= 2, 2)
  plot_data <- replace(plot_data, plot_data <= -2, -2)
  
  return(pheatmap(t(plot_data), color = PurpleAndYellow(),
                 annotation_col = col_ann[, rev(col_ann_order), drop = FALSE], annotation_row = row_ann, annotation_colors = ann_colours,
                 cluster_rows = cluster_rows, cluster_cols = cluster_cols, show_colnames = show_colnames, show_rownames = show_rownames, 
                 gaps_col = gaps_col, gaps_row = gaps_row, annotation_names_row = annotation_names_row, ...))
}