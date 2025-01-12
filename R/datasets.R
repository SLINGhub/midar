#' Plasma Lipidomics Dataset with Metadata
#'
#' This is an demo dataset used in function usage examples in this package.
#' It is a small subset of a prprocessed dataset
#' containing raw peak areas and analytical metadata from a plasma lipidomics analysis
#' published in Tan et al, ATVB, 2022. The dataset is used to demonstrate
#'
#' @format A `MidarExperiment` object:
#' \describe{
#'   \item{dataset_orig}{A tibble with the original peak data.}
#'   \item{dataset}{A tibble with annotated data.}
#'   \item{annot_analyses}{A tibble with analyses annotations.}
#'   \item{annot_features}{A tibble with features annotations.}
#'   \item{annot_batches}{A tibble with batch annotations.}
#'   ...
#' }
"lipidomics_dataset"
