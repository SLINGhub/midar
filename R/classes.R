#' S4 Class Representing the MIDAR Dataset
#'
#' @description
#'
#' The `MidarExperiment` object is the core data structure utilized within the MiDAR workflow, encapsulating all relevant experimental data and metadata.
#' It also includes processing results, details of the applied processing steps, and the current status of the data.
#'
#'
#' @docType class
#'
#' @slot title Title of the experiment
#' @slot analysis_type Analysis type, one of "lipidomics", "metabolomics", "externalcalib", "others"
#' @slot feature_intensity_var Feature variable used as default for calculations
#' @slot dataset_orig Original imported analysis data. Required fields:
#' @slot dataset Processed analysis data. Required fields:
#' @slot dataset_filtered Processed analysis data. Required fields:
#' @slot annot_analyses Annotation of analyses/runs
#' @slot annot_features Annotation of measured features.
#' @slot annot_istds Annotation of Internal Standard concs.
#' @slot annot_responsecurves Annotation of  Response curves (RQC). Required fields
#' @slot annot_studysamples Annotation of study samples. Required fields:
#' @slot annot_batches Annotation of batches. Required fields:
#' @slot metrics_qc QC information for each measured feature
#' @slot parameters_processing Values of parameters used for the different processing steps
#' @slot status_processing Status within the data processing workflow
#' @slot is_istd_normalized Flag if data has been ISTD normalized
#' @slot is_quantitated Flag if data has been quantitated using ISTD and sample amount
#' @slot is_filtered Flag if data has been filtered based on QC parameters
#' @slot is_isotope_corr Flag if one or more features have been isotope corrected
#' @slot has_outliers_tech Flag if data has technical analysis/sample outliers
#' @slot analyses_excluded Flag if outliers were excluded in the QC-filtered dataset
#' @slot var_drift_corrected List indicating which variables are drift corrected
#' @slot var_batch_corrected List indicating which variables are batch corrected

#' @include midar-global-definitions.R
#' @export

setClass("MidarExperiment",
  slots = c(
    title = "character",
    analysis_type = "character",
    feature_intensity_var = "character",
    dataset_orig = "tbl_df",
    dataset = "tbl_df",
    dataset_filtered = "tbl_df",
    annot_analyses = "tbl_df",
    annot_features = "tbl_df",
    annot_istds = "tbl_df",
    annot_responsecurves = "tbl_df",
    annot_studysamples = "tbl_df",
    annot_batches = "tbl_df",
    metrics_qc = "tbl_df",
    parameters_processing = "tbl_df",
    status_processing = "character",
    is_istd_normalized = "logical",
    is_quantitated = "logical",
    is_filtered = "logical",
    has_outliers_tech = "logical",
    is_isotope_corr = "logical",
    analyses_excluded = "vector",
    var_drift_corrected = "vector",
    var_batch_corrected = "vector"
  ),
  prototype = list(
    title = "",
    analysis_type = "",
    feature_intensity_var = "",
    dataset_orig = pkg.env$table_templates$annot_analyses_template,
    dataset = pkg.env$table_templates$annot_analyses_template,
    dataset_filtered = pkg.env$table_templates$annot_analyses_template,
    annot_analyses = pkg.env$table_templates$annot_analyses_template,
    annot_features = pkg.env$table_templates$annot_features_template,
    annot_istds = pkg.env$table_templates$annot_istds_template,
    annot_responsecurves = pkg.env$table_templates$annot_responsecurves_template,
    annot_studysamples = tibble::tibble(),
    annot_batches = tibble::tibble(),
    metrics_qc = tibble::tibble(),
    parameters_processing = pkg.env$table_templates$parameters_processing_template,
    status_processing = "No Data",
    is_isotope_corr = FALSE,
    is_istd_normalized = FALSE,
    is_quantitated = FALSE,
    var_drift_corrected = c(feature_intensity = FALSE, feature_norm_intensity = FALSE, feature_conc = FALSE),
    var_batch_corrected = c(feature_intensity = FALSE, feature_norm_intensity = FALSE, feature_conc = FALSE),
    is_filtered = FALSE,
    has_outliers_tech = FALSE,
    analyses_excluded = NA
  )
)



#' Constructor for the MidarExperiment object.
#' @importFrom methods new
#' @param title Title of experiment
#' @param analysis_type Analysis type, one of "lipidomics", "metabolomics", "externalcalib", "others"
#' @return `MidarExperiment` object
#' @export
MidarExperiment <- function(title = "", analysis_type = "") {
  methods::new("MidarExperiment",
               title = title,
               analysis_type = analysis_type)

}

#' Set analysis type
#' @description
#' Set the analysis type, i.e. "lipidomics", "metabolomics, "quantitative"
#' @param x MidarExperiment object
#' @return A character string
#' @export
setGeneric("analysis_type", function(x) standardGeneric("analysis_type"))


#' Get analysis type
#' @description
#' Get the analysis type defined for the MidarExperiment object
#' @param x MidarExperiment object
#' @param value Analysis type, one of "lipidomics", "metabolomics, "quantitative"
#' @export
setGeneric("analysis_type<-", function(x, value) standardGeneric("analysis_type<-"))

#' Get `analysis_type`
#' @param x MidarExperiment object
#' @return A character string
#' @export
setMethod("analysis_type", "MidarExperiment", function(x) x@analysis_type)

#' Set `analysis_type`
#' @param x MidarExperiment object
#' @param value Analysis type, one of "lipidomics", "metabolomics, "quantitative"
#' @export
setMethod("analysis_type<-", "MidarExperiment", function(x, value) {
  x@analysis_type <- value
  x
})

# TODODO: ALIGN WITH ASSERTION and DEFINE WHERE WHEN TO RUN THIS
check_integrity <- function(data = NULL, excl_unmatched_analyses) {
  check_data(data)
  if (nrow(data@dataset_orig) > 0 & nrow(data@annot_analyses) > 0) {
    d_xy <- length(setdiff(data@dataset_orig$analysis_id |> unique(), data@annot_analyses$analysis_id))
    d_yx <- length(setdiff(data@annot_analyses$analysis_id, data@dataset_orig$analysis_id |> unique()))
    if (d_xy > 0) {
      if (d_xy == length(data@dataset_orig$analysis_id |> unique())) cli::cli_abort("Error: None of the measurements/samples have matching metadata . Please check data and metadata files.")
      if (!excl_unmatched_analyses) {
        # if (d_xy < 50) {
        #   writeLines(glue::glue(""))
        #   cli::cli_abort(call. = FALSE, glue::glue("No metadata present for {d_xy} of {data@dataset_orig$analysis_id |> unique() |> length()} analyses/samples: {paste0(setdiff(data@dataset_orig$analysis_id |> unique(), data@annot_analyses$analysis_id), collapse = ", ")}"))
        # } else {
        #   cli::cli_abort(call. = FALSE, glue::glue("{d_xy} of {data@dataset_orig$analysis_id |> unique() |> length()} measurements have no matching metadata."))
        # }
      } else {
        #cli_alert_warning(col_yellow(glue::glue("Note: {d_xy} of {data@dataset_orig$analysis_id |> unique() |> length()} measurements without matching metadata were excluded.")))
      }
    } else if (d_yx > 0) {
      if (d_yx < 50) {
        writeLines(glue::glue("Following {d_yx} analysis/samples present in measurement data are not defined in the metadata:  {paste0(setdiff(data@annot_analyses$analysis_id, data@dataset_orig$analysis_id |> unique()), collapse = ", ")}"))
      } else {
        writeLines(glue::glue("{d_yx} analysis/samples present in measurement data are not defined in the metadata (too many to show)"))
      }
      cli::cli_abort(glue::glue(""))
    } else {
      data@status_processing <- "DataMetadataLoaded"
      TRUE
    }
  }
}

#' @noRd
get_status_flag <- function(x) {
  ifelse(x, {cli::col_green(cli::symbol$tick)}, {cli::col_red(cli::symbol$cross)})
  }




#' Getter for specific slots of an MidarExperiments object
#'
#' $ syntax can be used to as a shortcut for getting specific variables and results from a MidarExperiment object
#' @return Value with a variable or a tibble
#' @param x MidarExperiment object
#' @param name MidarExperiment slot
#' @examples
#' mexp <- MidarExperiment()
#' mexp$analysis_type
#' mexp$annot_analyses
#' @importFrom methods slot
#' @export
setMethod(
  f = "$",
  signature = c("MidarExperiment"),
  definition = function(x, name) {
    # check for other struct slots
    valid <- c("title", "analysis_type", "dataset", "annot_analyses", "annot_features", "annot_istds", "metrics_qc", "annot_batches", "dataset_filtered", "is_istd_normalized", "var_drift_corrected", "var_batch_corrected")
    if (!name %in% valid) cli::cli_abort('"', name, '" is not valid for this object: ', class(x)[1])
    methods::slot(x, name)
  }
)


#' Check integrity of MidarExperiment data object
#'
#' @description
#' Helper function that checks the structure and contents of
#' a MidarExperiment object
#'
#' @param data MidarExperiment object
#' @return silent on success, prints abort message on fail
#' @noRd
check_data <- function(data = NULL){
    # fail if data is NULL
  if(is.null(data)) {
      cli::cli_div(theme = list(span.emph = list(color = "#e81744")))
      cli::cli_abort(c("x" = "`data` cannot be {.emph NULL}, please use a {.emph MidarExperiment}"))
    }
  if(is(data)[1] != 'MidarExperiment') {
    cli::cli_div(theme = list(span.emph = list(color = "#e81744")))
    cli::cli_abort(c("x" = "`data` must be a {.emph MidarExperiment}"))
  }
}

setMethod("show", "MidarExperiment", function(object) {

  cli::cli_par()
  cli::cli_h1(is(object)[[1]])
  cli::cli_text(cli::col_blue("Title: {.strong {object@title}}"))
  cli::cli_end()

  cli::cli_par()
  cli::cli_text(cli::col_blue("Processing status: {.strong {object@status_processing}}"))
  cli::cli_end()

  cli::cli_h2("Annotated Raw Data")
  cli::cli_ul(id = "A")
  cli::cli_li("Analyses: {length(unique(object@dataset$analysis_id))}")
  cli::cli_li("Features: {length(unique(object@dataset$feature_id))}")
  cli::cli_li("Raw signal used for processing: `{object@feature_intensity_var}`")
  cli::cli_end(id = "A")

  cli::cli_h2("Metadata")
  cli::cli_ul(id ="B")
  cli::cli_li("Sample annotation: {.strong {get_status_flag(nrow(object@annot_analyses) > 0)}}")
  cli::cli_li("Feature annotation: {.strong {get_status_flag(nrow(object@annot_features) > 0)}}")
  cli::cli_li("Internal standard annotation: {.strong {get_status_flag(nrow(object@annot_istds) > 0)}}")
  cli::cli_li("Response curve annotation:  {.strong {get_status_flag(nrow(object@annot_responsecurves) > 0)}}")
  cli::cli_li("Study samples annotation:  {.strong {get_status_flag(nrow(object@annot_studysamples) > 0)}}")
  cli::cli_end(id ="B")

  cli::cli_h2("Processing Status")
  cli::cli_ul(id = "C")
  cli::cli_li("Isotope corrected: {get_status_flag(object@is_isotope_corr)}")
  cli::cli_li("ISTD normalized: {get_status_flag(object@is_istd_normalized)}")
  cli::cli_li("ISTD quantitated: {get_status_flag(object@is_quantitated)}")

  get_corr_var <- function(vars) {
    vars_names <- names(vars)
    if (length(vars_names) == 0)
      return(cli::col_red(cli::symbol$cross))  # Red cross if the vector is empty
    else
      return(glue("`", stringr::str_flatten_comma(vars_names, last = " and ", na.rm = TRUE), "`"))
  }


  cli::cli_li("Drift corrected variables:  {get_corr_var(object@var_drift_corrected[object@var_drift_corrected])}")
  cli::cli_li("Batch corrected variables:  {get_corr_var(object@var_batch_corrected[object@var_batch_corrected])}")
  cli::cli_li("Feature filtering applied:  {get_status_flag(object@is_filtered)}")
  cli::cli_end(id = "C")
  cli::cli_h2("Exclusion of Analyses and Features")
  cli::cli_ul(id = "D")

  if (length(object@analyses_excluded) == 0)
    str <- cli::col_red(cli::symbol$cross)  # Red cross if the vector is empty
  else
    str <- glue::glue_collapse(object@analyses_excluded, sep = ", ", width = 80, last = ", and ")

  cli::cli_li("Analyses manually excluded (`analysis_id`): {col_red(str)}")
  cli::cli_end(id = "D")
})
