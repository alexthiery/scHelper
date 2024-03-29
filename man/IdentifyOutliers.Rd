% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IdentifyOutliers.R
\name{IdentifyOutliers}
\alias{IdentifyOutliers}
\title{Identify oulier clusters}
\usage{
IdentifyOutliers(
  seurat_obj,
  group_by = "seurat_clusters",
  metrics,
  intersect_metrics = TRUE,
  quantiles
)
}
\arguments{
\item{seurat_obj}{Seurat object}

\item{group_by}{column to group cells by}

\item{metrics}{array of metrics to test for outliers}

\item{intersect_metrics}{boolean value for whether to intersect the outliers for each metric}

\item{quantiles}{Percent quantiles to plot on bar plots. Must be an array (length 2) with lower and upper percent quantiles}
}
\value{
array of clusters which are outliers
}
\description{
This function identifies clusters of cells based on QC metrics of interest. Groups are classed as outliers if the median for a given metric falls outside the provided quantiles.
}
