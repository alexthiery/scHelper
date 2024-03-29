% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GeneModulePheatmap.R
\name{GeneModulePheatmap}
\alias{GeneModulePheatmap}
\title{Function for plotting gene module pheatmap from Seurat object}
\usage{
GeneModulePheatmap(
  seurat_obj,
  metadata,
  col_order = metadata[1],
  custom_order = NULL,
  custom_order_column = NULL,
  assay = "RNA",
  slot = "scale.data",
  gene_modules,
  selected_genes = NULL,
  hide_annotation = NULL,
  gaps_row = TRUE,
  gaps_col = NULL,
  gm_row_annotation = TRUE,
  cell_subset = NULL,
  use_seurat_colours = TRUE,
  colour_scheme = c("PRGn", "RdYlBu", "Greys"),
  col_ann_order = rev(metadata),
  show_colnames = FALSE,
  show_rownames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  order_genes = TRUE,
  annotation_names_row = FALSE,
  ...
)
}
\arguments{
\item{seurat_obj}{Seurat object}

\item{metadata}{a string or array of colnames in seurat_obj@meta.data used to annotate the cells.}

\item{col_order}{cell order will be prioritised by the order of elements in the array. By default cells are ordered by the first element of the metadata variable.}

\item{custom_order}{array specifying cell group ordering. If set then custom_order_column must also be set.}

\item{custom_order_column}{metadata column from which custom_order is derived.}

\item{assay}{seurat assay to use for plotting.}

\item{slot}{seurat slot to use for plotting.}

\item{gene_modules}{list of gene modules to plot.}

\item{selected_genes}{optional array of genes to filter gene modules. Modules with no selected_genes present will be removed.}

\item{hide_annotation}{array of elements to hide from column annotation.}

\item{gaps_row}{boolean specifying whether to add whitespace between gene modules.}

\item{gaps_col}{metadata column by which to add whitespace.}

\item{gm_row_annotation}{boolean specifying whether to add gene module row annotations.}

\item{cell_subset}{array of cells to subset for plotting. If specified, data will be re-scaled.}

\item{use_seurat_colours}{boolean specifying whether to use default seurat colours. If FALSE, colours from colour_scheme will be used instead.}

\item{colour_scheme}{array or colour schemes from RColourBrewer to be used for annotating columns.}

\item{col_ann_order}{can be specified to change the order of the column annotations.}

\item{show_colnames}{default parameter passed to pheatmap.}

\item{show_rownames}{default parameter passed to pheatmap.}

\item{cluster_rows}{default parameter passed to pheatmap.}

\item{cluster_cols}{default parameter passed to pheatmap.}

\item{order_genes}{boolean specifying whether to order genes within gms based on hclustering}

\item{annotation_names_row}{default parameter passed to pheatmap.}

\item{...}{Extra arguments to be passed to pheatmap.}
}
\value{
plot output from pheatmap.
}
\description{
col_ann_order can be specified to change the order in which the column annotations appear on the heatmap
}
