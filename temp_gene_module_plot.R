temp_gene_module_plot

GM.plot <- function (data, metadata, col_order = metadata[1], custom_order = NULL, 
          custom_order_column = NULL, assay = "RNA", slot = "scale.data", 
          gene_modules, selected_genes = NULL, main = "", hide_annotation = NULL, 
          show_rownames = TRUE, annotation_colors = NA, hclust_rows = FALSE, 
          hclust_cols = FALSE, gaps_row = TRUE, gaps_col = NULL, gm_row_annotation = TRUE, 
          cell_subset = NULL, treeheight_row = 0, use_seurat_colours = TRUE, 
          colour_scheme = c("PRGn", "RdYlBu", "Greys"), col_ann_order = rev(metadata), 
          ...) 
{
  # Subset cells if necessary
  if (!is.null(cell_subset)) {
    data <- subset(data, cells = cell_subset)
  }
  
  # Select metadata columns and convert to factor if character
  data@meta.data <- data@meta.data[,metadata] %>% mutate_if(is.character,as.factor)
  # Drop all unused levels
  data@meta.data[] <- lapply(data@meta.data, function(x) if(is.factor(x)) factor(x) else x)
  
  # Set HM.col
  HM.col <- data@meta.data
  
  # Order cells based on col_order
  HM.col <- HM.col[do.call("order", c(HM.col[col_order], list(decreasing = FALSE))), , drop = FALSE]
  
  # Order cells based on custom_order
  if (!is.null(custom_order)) {
    if (is.null(custom_order_column)) {
      "custom_order column must be specified \n"
    }
    if (!setequal(custom_order, unique(HM.col[[custom_order_column]]))) {
      stop("custom_order factors missing from custom_order_column \n\n")
    }
    levels(HM.col[[custom_order_column]]) <- custom_order
    HM.col <- HM.col[order(HM.col[[custom_order_column]]), , drop = FALSE]
  }
  
  # Subset gene modules based on selected genes
  if (!is.null(selected_genes)) {
    selected_GM <- SubsetGeneModules(gm = gene_modules, selected_genes = selected_genes)
  } else {
    if (is.null(names(gene_modules))) {
      names(gene_modules) <- paste0("GM:", 1:length(gene_modules))
    }
    selected_GM <- gene_modules
  }
  
  # Set row annotations
  if (gm_row_annotation == TRUE) {
    stack(selected_GM) %>% rename("Gene Modules" = ind) %>% column_to_rownames('values')
  } else {
    row_ann = NA
  }

  ##UP TO HERE

  
  # Insert whitespace between rows
  if (gaps_row == TRUE) {
    row_ann <- droplevels(row_ann)
    gaps_row = cumsum(summary(as.factor(row_ann[["Gene Modules"]]), 
                              maxsum = max(lengths(lapply(row_ann, unique)))))
  } else {
    gaps_row = NULL
  }
  
  # Insert whitespace between columns
  if (!is.null(gaps_col)) {
    if (class(gaps_col) != "character") {
      stop("gaps_col must be a metadata column name")
    } else {
      gaps_col = cumsum(rle(as.vector(HM.col[[gaps_col]]))[["lengths"]])
    }
  }
  
  # Hide annotations
  if (!is.null(hide_annotation)) {
    HM.col[, hide_annotation] <- NULL
  }
  
  if (use_seurat_colours == FALSE) {
    ann_colours <- list()
    for (tic in 1:ncol(HM.col)) {
      ann_colours[[colnames(HM.col[tic])]] <- setNames(colorRampPalette(brewer.pal(9, 
                                                                                   colour_scheme[tic])[2:9])(length(unique(HM.col[, 
                                                                                                                                  tic]))), unique(HM.col[, tic]))
    }
  }
  else {
    ann_colours <- list()
    for (column in colnames(HM.col)) {
      ann_colours[[column]] <- setNames(ggplotColours(n = length(levels(HM.col[, column]))), levels(HM.col[, column]))
      
      ann_colours[[column]] <- ann_colours[[column]][match(levels(HM.col[[column]]), names(ann_colours[[column]]))]
    }
  }
  ann_colours[["Gene Modules"]] <- setNames(colorRampPalette(brewer.pal(9, 
                                                                        "Paired"))(length(unique(row_ann$`Gene Modules`))), unique(row_ann$`Gene Modules`))
  new.dat <- t(as.matrix(x = GetAssayData(object = data, assay = assay, 
                                          slot = slot)[unlist(selected_GM), rownames(HM.col), drop = FALSE]))
  if (!is.null(cell_subset)) {
    cat("rescaling data as cells have been subset \n")
    new.dat <- t(scale(t(new.dat)))
  }
  else {
  }
  new.dat <- replace(new.dat, new.dat >= 2, 2)
  new.dat <- replace(new.dat, new.dat <= -2, -2)
  print(pheatmap(t(new.dat), color = PurpleAndYellow(), cluster_rows = hclust_rows, 
                 cluster_cols = hclust_cols, show_colnames = FALSE, annotation_col = HM.col[, 
                                                                                            rev(col_ann_order), drop = FALSE], gaps_col = gaps_col, 
                 gaps_row = gaps_row, main = main, show_rownames = show_rownames, 
                 annotation_row = row_ann, annotation_colors = ann_colours, 
                 treeheight_row = treeheight_row, ...))
}
