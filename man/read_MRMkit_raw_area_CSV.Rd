% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_msdata.R
\name{read_MRMkit_raw_area_CSV}
\alias{read_MRMkit_raw_area_CSV}
\title{Read and convert an Agilent MassHunter Quant CSV result file}
\usage{
read_MRMkit_raw_area_CSV(filename, use_mrmkit_normdata = FALSE, silent = FALSE)
}
\arguments{
\item{filename}{File path of MassHunter Quant CSV file}

\item{use_mrmkit_normdata}{use raw or MRMkit-normalized data}

\item{silent}{Suppress messages}
}
\value{
A tibble in the long format
}
\description{
Read and convert an Agilent MassHunter Quant CSV result file
}
