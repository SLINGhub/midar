% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_metadata.R
\name{import_msorganizer_xlm}
\alias{import_msorganizer_xlm}
\title{Imports metadata provided by an MSOrganizer EXCEL template}
\usage{
import_msorganizer_xlm(filename, trim_ws = TRUE)
}
\arguments{
\item{filename}{File path of the MSOrganizer EXCEL template (*.xlm)}

\item{trim_ws}{Trim all white spaces and double spaces}
}
\value{
A list of tibbles with different metadata
}
\description{
Imports metadata provided by an MSOrganizer EXCEL template
}
