% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\docType{data}
\name{quant_lcms_dataset}
\alias{quant_lcms_dataset}
\title{LC-MS Dataset with External Calibration Curve and Metadata}
\format{
A \code{MidarExperiment} object with the following data and metadata:
\describe{
\item{dataset_orig}{Original data (peak datas).}
\item{dataset}{Annotated data}
\item{annot_analyses}{Analysis-level metadata}
\item{annot_features}{Feature-level annotations}
\item{annot_istds}{ISTD concentrations}
\item{annot_qcconcentrations}{Calibrant (\code{CAL}) and QC concentrations}
}
}
\usage{
quant_lcms_dataset
}
\description{
This demo dataset is included for use in function examples and user testing.
It is a subset of an LC-MS analysis of plasma steroids, containing an external
calibration curve for each analyte, QC samples with known concentrations,
and unknown samples.
}
\keyword{datasets}
