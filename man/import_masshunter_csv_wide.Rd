% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_msdata.R
\name{import_masshunter_csv_wide}
\alias{import_masshunter_csv_wide}
\title{Reads a wide CSV file with Feature Intensities}
\usage{
import_masshunter_csv_wide(file, field, silent = FALSE)
}
\arguments{
\item{file}{File name and path of the MassHunter Quant CSV file}

\item{field}{Peak parameter (e.g. Area, RT)}

\item{silent}{Suppress messages}
}
\value{
A tibble in the long format
}
\description{
Reads a wide CSV file with Feature Intensities
}
