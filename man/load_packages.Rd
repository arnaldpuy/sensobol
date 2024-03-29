% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_packages.R
\name{load_packages}
\alias{load_packages}
\title{Load (and install) R packages.}
\usage{
load_packages(x)
}
\arguments{
\item{x}{A character vector with the name of the packages to load.}
}
\description{
The function loads R packages. If the packages are not already
in the local system, the function also downloads, installs and loads them.
}
\examples{
# Load packages:
\dontrun{load_packages(c("tidyverse", "data.table"))}
}
