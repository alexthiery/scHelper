% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArchR_ClustRes.R
\name{ArchR_ClustRes}
\alias{ArchR_ClustRes}
\title{ArchR function: Plot multiple cluster resolutions}
\usage{
ArchR_ClustRes(ArchR, by = 0.1, starting_res = 0)
}
\arguments{
\item{ArchR}{ArchR object}

\item{by}{value by which resolution is incremented}

\item{starting_res}{lowest resolution to cluster cells by}
}
\value{
plot generated by grid extra
}
\description{
This function plots UMAPs for a series of nine cluster resolutions, alongside a clustree which reveals the relationship between clusters at different resolutions
}
