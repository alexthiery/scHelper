% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ArchR_subset.R
\name{ArchR_Subset}
\alias{ArchR_Subset}
\title{ArchR function: Subset ArchR}
\usage{
ArchR_Subset(
  ArchR,
  meta_col1,
  meta_col2,
  groups1,
  groups2,
  invert1,
  invert2,
  invert = FALSE
)
}
\arguments{
\item{ArchR}{ArchR object}

\item{meta_col1}{first column to split cells by}

\item{meta_col2}{second column to split cells by}

\item{groups1}{first groups to subset by}

\item{groups2}{second groups to subset by}

\item{invert1}{invert first groups}

\item{invert2}{invert second groups}

\item{invert}{invert all groups}
}
\value{
subsetted ArchR object
}
\description{
This function subsets ArchR specified groups
}
