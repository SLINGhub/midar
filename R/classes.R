pkg.env <- new.env()

# Data structure templates
pkg.env$dataset_templates <- list(
  dataset_orig_template = tibble::tibble(
    "raw_data_filename" = character(),
    "acquisition_time_stamp" = as.Date(character()),
    "feature_intensity" = numeric(),
    "feature_norm_intensity" = numeric(),
    "feature_conc" = numeric()
  ),

  dataset_template = tibble::tibble(
    "run_id" = integer(),
    "analysis_id" = character(),
    "raw_data_filename" = character(),
    "acquisition_time_stamp" = as.Date(character()),
    "qc_type" = factor(),
    "replicate_no" = integer(),
    "batch_id" = character(),
    "valid_analysis" = logical(),
    "outlier_technical" = logical(),
    "outlier_technical_note" = character(),
    "feature_name" = character(),
    "feature_class" = character(),
    "is_istd" = logical(),
    "valid_integration" = logical(),
    "feature_intensity" = numeric(),
    "feature_norm_intensity" = numeric(),
    "feature_conc" = numeric()
  ),

  annot_analyses_template = tibble::tibble(
    "analysis_id" = character(),
    "raw_data_filename"= character(),
    "sample_id" = character(),
    "qc_type" = factor(),
    "replicate_no" = integer(),
    "method_id" = character(),
    "batch_id" = character(),
    "specimen" = character(),
    "sample_amount" = numeric(),
    "sample_amount_unit" = character(),
    "istd_volume" = numeric(),
    "valid_analysis" = logical(),
    "outlier_technical" = logical(),
    "outlier_technical_note" = character(),
    "remarks" = character()
  ),
  annot_features_template = tibble::tibble(
    "feature_name" = character(),
    "feature_class" = character(),
    "is_istd" = logical(),
    "norm_istd_feature_name" = character(),
    "quant_istd_feature_name" = character(),
    "is_quantifier" = logical(),
    "valid_integration" = logical(),
    "feature_response_factor" = numeric(),
    "remarks" = character()
  ),
  annot_istd_template = tibble::tibble(
    "norm_istd_feature_name" = character(),
    "quant_istd_feature_name" = character(),
    "istd_conc_nmolar" = numeric()
  ),
  annot_responsecurves_template = tibble::tibble(
    "analysis_id" = character(),
    "rqc_series_id" = character(),
    "relative_sample_amount" = numeric(),
    "injection_volume" = numeric()
  ),
  annot_batch_info_template = tibble::tibble(
    "batch_id" = character(),
    "batch_no" = numeric(),
    "id_batch_start" = numeric(),
    "id_batch_end" = numeric()
  ),

  parameters_processing_template = tibble::tibble(
    "parameter_name" = character()
  )


)

pkg.env$qc_type_annotation <- list(
  qc_type_levels = c("SBLK", "TBLK", "UBLK", "PQC", "TQC", "BQC", "RQC", "EQC", "NIST",
                     "LTR", "PBLK", "SPL", "SST", "MBLK"),

  qc_type_col = c(
    "SBLK" = "#1854f9",
    "TBLK" = "#db0202",
    "UBLK" = "#de21de",
    "BQC" = "#db0202",
    "TQC" = "#1854f9",
    "PQC" = "#f99f18",
    "RQC" = "#96a4ff",
    "EQC" = "#513c3c",
    "NIST" = "#002e6b",
    "LTR" = "#880391",
    "PBLK" = "#08c105",
    "SPL" = "#899ea3",
    "SST" = "#bafc03",
    "MBLK" = "black"),

  qc_type_fillcol = c(
    "SBLK" = "#f891ff",
    "TBLK" = "#fffb03",
    "UBLK" = "#c1bd04",
    "BQC" = "#db0202",
    "TQC" = "#1854f9",
    "PQC" = "#f99f18",
    "RQC" = "#688ff9",
    "EQC" = "NA",
    "NIST" = "#cce2ff",
    "LTR" = "#880391",
    "PBLK" = "#08c105",
    "SPL" = "NA",
    "SST" = "#aaaeaf",
    "MBLK" = "black"),

  qc_type_shape = c(
    "SBLK" = 23,
    "TBLK" = 23,
    "UBLK" = 23,
    "BQC" = 16,
    "TQC" = 25,
    "PQC" = 25,
    "RQC" = 6,
    "EQC" = 24,
    "NIST" = 23,
    "LTR" = 23,
    "PBLK" = 23,
    "SPL" = 21,
    "SST" = 10,
    "MBLK" = 10)
)

# allow S4 to see the class tbl_df
setOldClass(c("tbl_df", "tbl", "data.frame"))

#' S4 Class Representing the MIDAR Dataset
#'
#' @description
#' Core class to store raw and processed data, annotations and metadata.
#'
#' @details
#' Provides also functions to validate data integrity, to perfom normalization and quantitation based on internal standards and to calculate the QC parameters of the experiment
#'
#' @docType class
#'
#' @slot analysis_type Analysis type, one of "lipidomics", "metabolomics", "externalcalib", "others"
#' @slot dataset_orig Original imported analysis data. Required fields:
#' @slot dataset Processed analysis data. Required fields:
#' @slot dataset_filtered Processed analysis data. Required fields:
#' @slot annot_analyses Annotation of analyses/runs
#' @slot annot_features Annotation of measured features.
#' @slot annot_istd Annotation of Internal Standard concs.
#' @slot annot_responsecurves Annotation of  Response curves (RQC). Required fields
#' @slot annot_studysamples Annotation of study samples. Required fields:
#' @slot annot_batches Annotation of batches. Required fields:
#' @slot metrics_qc QC information for each measured feature
#' @slot parameters_processing Values of parameters used for the different processing steps
#' @slot status_processing Status within the data processing workflow
#' @slot is_istd_normalized Flag if data has been ISTD normalized
#' @slot is_quantitated Flag if data has been quantitated using ISTD and sample amount
#' @slot is_drift_corrected Flag if data has been drift corrected
#' @slot is_batch_corrected Flag if data has been batch corrected
#' @slot is_isotope_corr Flag if one or more features have been isotope corrected
#' @slot has_outliers_tech Flag if data has technical analysis/sample outliers
#' @slot excl_outliers_tech Flag if outliers were excluded in the QC-filtered dataset
#'
#' @export
#'
#' @importFrom utils tail
#' @importFrom tibble tibble

setClass("MidarExperiment",
         slots = c(
           analysis_type = "character",
           dataset_orig = "tbl_df",
           dataset = "tbl_df",
           dataset_filtered = "tbl_df",
           annot_analyses = "tbl_df",
           annot_features = "tbl_df",
           annot_istd = "tbl_df",
           annot_responsecurves= "tbl_df",
           annot_studysamples = "tbl_df",
           annot_batches = "tbl_df",
           metrics_qc = "tbl_df",
           parameters_processing = "tbl_df",
           status_processing = "character",
           is_istd_normalized = "logical",
           is_quantitated = "logical",
           is_drift_corrected = "logical",
           is_batch_corrected = "logical",
           has_outliers_tech = "logical",
           is_isotope_corr = "logical",
           excl_outliers_tech = "logical"
         ),
         prototype = list(
           analysis_type = "",
           dataset_orig = pkg.env$dataset_templates$annot_analyses_template,
           dataset = pkg.env$dataset_templates$annot_analyses_template,
           dataset_filtered = pkg.env$dataset_templates$annot_analyses_template,
           annot_analyses = pkg.env$dataset_templates$annot_analyses_template,
           annot_features = pkg.env$dataset_templates$annot_features_template,
           annot_istd = pkg.env$dataset_templates$annot_istd_template,
           annot_responsecurves= pkg.env$dataset_templates$annot_responsecurves_template,
           annot_studysamples = tibble::tibble(),
           annot_batches = tibble::tibble(),
           metrics_qc = tibble::tibble(),
           parameters_processing = pkg.env$dataset_templates$parameters_processing_template,
           status_processing = "No Data",
           is_istd_normalized = FALSE,
           is_quantitated = FALSE,
           is_drift_corrected = FALSE,
           is_batch_corrected = FALSE,
           has_outliers_tech = FALSE,
           is_isotope_corr = FALSE,
           excl_outliers_tech = FALSE
         )
)

#' Constructor for the MidarExperiment object.
#' @importFrom methods new
#' @param analysis_type Analysis type, one of "lipidomics", "metabolomics", "externalcalib", "others"
#' @return `MidarExperiment` object
#' @export
MidarExperiment <- function(analysis_type = "") {
  methods::new("MidarExperiment", analysis_type = analysis_type)
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
setMethod("analysis_type", "MidarExperiment", function(x) x@analysis_type)

#' Set `analysis_type`
#' @param x MidarExperiment object
#' @param value Analysis type, one of "lipidomics", "metabolomics, "quantitative"
setMethod("analysis_type<-", "MidarExperiment", function(x, value) {
  x@analysis_type <- value
  x
})


check_integrity <-  function(object, excl_unannotated_analyses) {
  #browser()
  if (nrow(object@dataset_orig) > 0 & nrow(object@annot_analyses) > 0) {
    d_xy <- length(setdiff(object@dataset_orig$raw_data_filename %>% unique(), object@annot_analyses$raw_data_filename))
    d_yx <- length(setdiff(object@annot_analyses$raw_data_filename,object@dataset_orig$raw_data_filename %>% unique()))
    if(d_xy > 0){
      if(d_xy == length(object@dataset_orig$raw_data_filename %>% unique())) stop("Error: None of the measurements/samples have matching metadata . Please check data and metadata files.")
        if(!excl_unannotated_analyses){
          stop(glue::glue("Error: {d_xy} of {object@dataset_orig$raw_data_filename %>% unique() %>% length()} measurements have no matching metadata."))
        if (d_xy < 50)
          writeLines(glue::glue("No metadata present for: {paste0(setdiff(object@dataset_orig$raw_data_filename %>% unique(), object@annot_analyses$raw_data_filename), collapse = ", ")} measurements."))
        else
          print("No metadata present for: Too many (> 50) to display")
        } else {
            writeLines(crayon::yellow(glue::glue("! Note: {d_xy} of {object@dataset_orig$raw_data_filename %>% unique() %>% length()} measurements without matching metadata were excluded.")))
        }
      } else if(d_yx > 0) {
        stop(glue::glue("{d_yx} of {object@annot_analyses$raw_data_filename %>% length()} sample metadata are not found in the measurement data."))
      } else {
      object@status_processing <- "DataMetadataLoaded"
      TRUE
      }
  }
}

##' @importFrom methods setValidity
##'
#methods::setValidity("MidarExperiment")


get_status_flag <- function(x) if_else(x, crayon::green$bold('\u2713'), crayon::red$bold('\u2717'))




#' Getter for specific slots of an MidarExperiments object
#'
#' $ syntax can be used to as a shortcut for getting specific variables and results from a MidarExperiment object
#' @return Value with a variable or a tibble
#' @param x MidarExperiment object
#' @param name MidarExperiment slot
#' @examples
#' mexp = MidarExperiment()
#' mexp$analysis_type
#' mexp$annot_analyses
#' @importFrom methods slot
#' @export
setMethod(f = "$",
  signature = c("MidarExperiment"),
  definition = function(x,name) {

  # check for other struct slots
  valid=c('analysis_type','dataset','annot_analyses', 'annot_features', 'annot_istd', 'metrics_qc', 'annot_batches', 'dataset_filtered', 'is_istd_normalized')
  if (!name %in% valid) stop('"', name, '" is not valid for this object: ', class(x)[1])
  methods::slot(x,name)
  }
)



setMethod("show", "MidarExperiment", function(object) {
  cat("\n", is(object)[[1]], "\n",
      "\n",
      "  Processing status: ",object@status_processing, "\n",
      "\n",
      "  Data: ", "\n",
      "  \u2022 Samples: ", length(unique(object@dataset$analysis_id)), "\n",
      "  \u2022 Features:  ", length(unique(object@dataset$feature_name)), "\n",
      "\n",
      "  Metadata: ", "\n",
      "  \u2022 Sample annotation: ", get_status_flag(nrow(object@annot_analyses) > 0), "\n",
      "  \u2022 Feature annotation: ", get_status_flag(nrow(object@annot_features) > 0), "\n",
      "  \u2022 Internal standard annotation: ", get_status_flag(nrow(object@annot_istd) > 0), "\n",
      "  \u2022 Response curves annotation: ", get_status_flag(nrow(object@annot_responsecurves) > 0), "\n",
      "  \u2022 Study samples annotation: ", get_status_flag(nrow(object@annot_studysamples) > 0), "\n",
      "\n",
      "  Processing: ", "\n",
      "  \u2022 Isotope corrected: ", get_status_flag(object@is_isotope_corr), "\n",
      "  \u2022 ISTD normalized: ", get_status_flag(object@is_istd_normalized), "\n",
      "  \u2022 ISTD quantitated: ", get_status_flag(object@is_quantitated), "\n",
      "  \u2022 Drift corrected: ", get_status_flag(object@is_drift_corrected), "\n",
      "  \u2022 Batch corrected: ", get_status_flag(object@is_batch_corrected), "\n",
      "\n",
      "  Outliers: ", "\n",
      "  \u2022 Technical Outliers detected ", get_status_flag(object@has_outliers_tech), "\n",
      "  \u2022 Technical Outliers excluded from filtered data ", get_status_flag(object@excl_outliers_tech), "\n"
  )
})






