#' Plasma Lipidomics Dataset with Metadata
#'
#' This demo dataset is included for use in function examples and user testing..
#' It is a small, preprocessed subset of a plasma lipidomics dataset,
#' containing raw peak areas and analytical metadata. The original dataset
#' was published in Tan et al., ATVB, 2022.
#'
#' @format A `MidarExperiment` object with the following data and metadata:
#' \describe{
#' \item{dataset_orig}{A tibble containing the original peak data.}
#' \item{dataset}{A tibble with annotated lipidomics data.}
#' \item{annot_analyses}{Analysis-level metadata}
#' \item{annot_features}{Feature-level metadata}
#' \item{annot_batches}{Batch annotations.}
#' \item{annot_istds}{ISTD concentrations}
#' \item{annot_responsecurves}{Response curves metadata}
#' }
"lipidomics_dataset"

#' LC-MS Dataset with External Calibration Curve and Metadata
#'
#' This demo dataset is included for use in function examples and user testing.
#' It is a subset of an LC-MS analysis of plasma steroids, containing an external
#' calibration curve for each analyte, QC samples with known concentrations,
#' and unknown samples.
#'
#' @format A `MidarExperiment` object with the following data and metadata:
#' \describe{
#'   \item{dataset_orig}{Original data (peak datas).}
#'   \item{dataset}{Annotated data}
#'   \item{annot_analyses}{Analysis-level metadata}
#'   \item{annot_features}{Feature-level annotations}
#'   \item{annot_istds}{ISTD concentrations}
#'   \item{annot_qcconcentrations}{Calibrant (`CAL`) and QC concentrations}
#' }
"quant_lcms_dataset"

