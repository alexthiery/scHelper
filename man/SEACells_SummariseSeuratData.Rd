% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SEACells_SummariseSeuratData.R
\name{SEACells_SummariseSeuratData}
\alias{SEACells_SummariseSeuratData}
\title{SEACells function: Summarises data from seurat object across SEACells}
\usage{
SEACells_SummariseSeuratData(
  seurat,
  data_slot = "counts",
  category = "SEACell"
)
}
\arguments{
\item{seurat}{seurat object}

\item{data_slot}{which assay in seurat object to summarise, default is "counts"}

\item{category}{slot in seurat metadata by which to group cells for summarising, default is "SEACell"}
}
\value{
dataframe with summarised data across groupings
}
\description{
This function takes a seurat object and summarises the data from data_slot by SEACells grouping
}
