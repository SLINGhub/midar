% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correct-drift-batch.R
\name{correct_drift_gam}
\alias{correct_drift_gam}
\title{Drift Correction by Generalized Additive Model (GAM) Smoothing}
\usage{
correct_drift_gam(
  data = NULL,
  variable,
  ref_qc_types,
  batch_wise = TRUE,
  ignore_istd = TRUE,
  replace_previous = TRUE,
  bs = "ps",
  k = -1,
  sp = NULL,
  log_transform_internal = TRUE,
  conditional_correction = FALSE,
  recalc_trend_after = FALSE,
  feature_list = NULL,
  cv_diff_threshold = 0,
  use_original_if_fail = FALSE
)
}
\arguments{
\item{data}{MidarExperiment object}

\item{variable}{The variable to be corrected for drift effects. Must be one of "intensity", "norm_intensity", or "conc"}

\item{ref_qc_types}{QC types used for drift correction}

\item{batch_wise}{Logical. Apply the correction to each batch separately (\code{TRUE}, default) or across all batches (\code{FALSE}).}

\item{ignore_istd}{Logical. Exclude internal standards (ISTDs) from correction if \code{TRUE}.}

\item{replace_previous}{Logical. Replace existing correction (\code{TRUE}, default) or layer on top of it (\code{FALSE}).}

\item{bs}{Basis type for the spline: \code{"ps"} (penalized spline, default) or \code{"tp"} (thin plate spline).}

\item{k}{Number of basis functions (default: \code{-1}, automatically chosen by GAM).}

\item{sp}{Smoothing parameter (\code{NULL} by default, estimated automatically).}

\item{log_transform_internal}{Log transform the data for correction when \code{TRUE} (the default). Note: log transformation is solely applied internally for smoothing, results will not be log-transformed.}

\item{conditional_correction}{Determines whether drift correction should be applied to all features unconditionally (\code{TRUE}) or conditionally, based on sample CV change.}

\item{recalc_trend_after}{Recalculate trend post-drift correction for \code{plot_qc_runscatter()}. This will double calculation time.}

\item{feature_list}{Subset the features for correction whose names match the specified text using regular expression. Default is \code{NULL}.}

\item{cv_diff_threshold}{Maximum allowable change in CV ratio before and after smoothing for correction to be applied.}

\item{use_original_if_fail}{Determines the action when smoothing fails or results in invalid values for a feature. If \code{FALSE} (default), the result for each feature will \code{NA} for all batches, if \code{TRUE}, the original data is kept.}
}
\value{
MidarExperiment object
}
\description{
This function corrects for run-order drifts within or across batches using
Generalized Additive Models (GAMs). The correction uses penalized splines,
with automatic selection of smoothing parameters based on cross-validation or
penalized likelihood. It is typically based on QC (Quality Control) samples
measured at specific intervals throughout the run sequence.
The correction can be applied to "intensity",
"norm_intensity", or "conc" data.

It is recommended to visually inspect the correction using the
\code{\link[=plot_runscatter]{plot_runscatter()}} function. Set the argument
\code{recalc_trend_after = TRUE} so that the trends after correction are also
available for plotting. For further details, refer to the description
of \code{\link[=plot_runscatter]{plot_runscatter()}}.

The corrections can be applied on a batch-by-batch basis (\code{batch_wise = TRUE}, default) or across all batches (\code{batch_wise = FALSE}). Existing
corrections are either replaced (\code{replace_previous = TRUE}) or added on top
of them (\code{replace_previous = FALSE}).

Furthermore, drift correction can be applied unconditionally
(\code{conditional_correction = FALSE}) or conditionally, based on whether the
sample CV change before and after correction is below a defined
threshold (\code{cv_diff_threshold}). This conditional correction is assessed
independently for each batch if \code{batch_wise = TRUE}, where the median of
the CV changes across the batch is compared with the threshold.

\strong{Note}: The function outputs a message indicating the median CV change
and the mean absolute CV before and after correction for all samples.
However, these metrics are experimental and should not be used as
definitive criteria for correction (see Details below).

This cubic spline method is implemented using the base R function
\code{\link[stats:splinefun]{stats::spline()}}.
}
\details{
In the output message, the median CV change is computed as the median of CV changes for all
features in global correction or for
features where the correction passed the defined CV difference treshold in
case of conditional correction  (\code{conditional_correction = FALSE}).
For batch-wise correction, the change is calculated per batch, with the final median CV
change being the median of these batch medians across features.

This smoothing is based on Generalized Additive Models (GAM) using penalized splines, implemented via \code{mgcv::gam()}.
}
\seealso{
\code{\link[mgcv:gam]{mgcv::gam()}}
}
