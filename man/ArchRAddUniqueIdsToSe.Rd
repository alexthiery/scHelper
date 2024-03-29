% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArchR_AddUniqueIdsToSe.R
\name{ArchRAddUniqueIdsToSe}
\alias{ArchRAddUniqueIdsToSe}
\title{ArchR function: Add unique id to se object}
\usage{
ArchRAddUniqueIdsToSe(seMarker, ArchR, matrix_type)
}
\arguments{
\item{seMarker}{the se object created \code{ArchR::getMarkerFeatures}}

\item{ArchR}{ArchR object}

\item{matrix_type}{which matrix to use ("GeneScoreMatrix" or "PeakMatrix")}
}
\value{
seMarker with unique_id column
}
\description{
This function takes seMarker and adds a unique_id column
}
