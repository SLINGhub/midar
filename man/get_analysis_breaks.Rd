% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-managment.R
\name{get_analysis_breaks}
\alias{get_analysis_breaks}
\title{Get the number of analysis breaks in the analysis}
\usage{
get_analysis_breaks(data, break_duration_minutes)
}
\arguments{
\item{data}{A \code{MidarExperiment} object}

\item{break_duration_minutes}{A numeric value specifying the minimum duration
(in minutes) between two consecutive analyses that qualifies as an interruption.}
}
\value{
An integer with the number of interruptions, or \code{NA_integer_} if the dataset is empty.
}
\description{
Counts the number of interruptions in the analysis, where an interruption is
defined as a time gap between consecutive acquisition timestamps that
exceeds a given threshold (\code{break_duration_minutes}).
}
