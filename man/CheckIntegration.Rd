% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CheckIntegration.R
\name{CheckIntegration}
\alias{CheckIntegration}
\title{Check Integration}
\usage{
CheckIntegration(
  seurat_object,
  group_by = "orig.ident",
  split_chr = "-",
  pt.size = 0.2,
  xlim = c(-15, 15),
  ylim = c(-15, 15),
  nrow = 1,
  shuffle = TRUE,
  seed = 123,
  title = "all cells"
)
}
\arguments{
\item{seurat_object}{Seurat object}

\item{group_by}{seurat_object@meta.data column to group cells by - this should correspond to the variable which identifies shared populations across batches}

\item{split_chr}{character or string used to split group_by variable by batch (i.e. for 'hh6-1' you should specify '-')}

\item{nrow}{number of rows to pass to grid.arrange}

\item{shuffle}{randomly plot cells in order to prevent hiding certain populations when overplotting}

\item{seed}{set seed for plotting consistency when shuffle is TRUE}

\item{title}{title for plot displaying all cells}
}
\value{
plot generated by grid arrange
}
\description{
This function plots UMAPs after integration, allowing you to see that cells from the same population but sequenced as part of separate batches, overlap after integration.
}
