% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{MidarExperiment-class}
\alias{MidarExperiment-class}
\title{S4 Class Representing the MIDAR Dataset}
\description{
The \code{MidarExperiment} object is the core data structure utilized within the MiDAR workflow, encapsulating all relevant experimental data and metadata.
It also includes processing results, details of the applied processing steps, and the current status of the data.
}
\section{Slots}{

\describe{
\item{\code{title}}{Title of the experiment}

\item{\code{analysis_type}}{Analysis type, one of "lipidomics", "metabolomics", "externalcalib", "others"}

\item{\code{feature_intensity_var}}{Feature variable used as default for calculations}

\item{\code{dataset_orig}}{Original imported analysis data. Required fields:}

\item{\code{dataset}}{Processed analysis data. Required fields:}

\item{\code{dataset_filtered}}{Processed analysis data. Required fields:}

\item{\code{annot_analyses}}{Annotation of analyses/runs}

\item{\code{annot_features}}{Annotation of measured features.}

\item{\code{annot_istds}}{Annotation of Internal Standard concs.}

\item{\code{annot_responsecurves}}{Annotation of response curves (RQC). Required fields}

\item{\code{annot_qcconcentrations}}{Annotation of calibration curves. Required fields}

\item{\code{annot_studysamples}}{Annotation of study samples. Required fields:}

\item{\code{annot_batches}}{Annotation of batches. Required fields:}

\item{\code{metrics_qc}}{QC information for each measured feature}

\item{\code{metrics_calibration}}{Calibration metrics calculated from external calibration curves for each measured feature}

\item{\code{parameters_processing}}{Values of parameters used for the different processing steps}

\item{\code{status_processing}}{Status within the data processing workflow}

\item{\code{is_istd_normalized}}{Flag if data has been ISTD normalized}

\item{\code{is_quantitated}}{Flag if data has been quantitated using ISTD and sample amount}

\item{\code{is_filtered}}{Flag if data has been filtered based on QC parameters}

\item{\code{is_isotope_corr}}{Flag if one or more features have been isotope corrected}

\item{\code{has_outliers_tech}}{Flag if data has technical analysis/sample outliers}

\item{\code{analyses_excluded}}{Analyses excluded from processing, plots and reporting, unless explicitly requested}

\item{\code{features_excluded}}{Features excluded from processing, plots and reporting, unless explicitly requested}

\item{\code{var_drift_corrected}}{List indicating which variables are drift corrected}

\item{\code{var_batch_corrected}}{List indicating which variables are batch corrected}
}}

