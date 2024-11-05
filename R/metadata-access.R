
# Response Curves -----

#' metadata_responsecurves method
#' @description
#' Get curve IDs and sample amount information for response curves
#' @param x MidarExperiment object
#' @return Tibble
#' @export
setGeneric("metadata_responsecurves", function(x) standardGeneric("metadata_responsecurves"))

#' Get response curve metadata
#' @description
#' Get curve IDs and sample amount information for response curves
#' @param x MidarExperiment object
#' @return Tibble
#' @export
setMethod("metadata_responsecurves", "MidarExperiment", function(x) x@annot_responsecurves)


#' metadata_responsecurves method
#' @description
#' Get curve IDs and sample amount information for response curves
#' @param x MidarExperiment object
#' @param value Table with response curve metadata
#' @return Tibble
#' @export
setGeneric("metadata_responsecurves<-", function(x, value) standardGeneric("metadata_responsecurves<-"))

#' Set response curve metadata
#' @description
#' Set curve IDs and sample amount information for response curves
#' @param x MidarExperiment object
#' @param value A data.frame or tibble
#' @return MidarExperiment object
#' @export
# TODO: cleanup/remove/merge validation checks .. use assertion functions defined else in the code
setMethod("metadata_responsecurves<-", "MidarExperiment", function(x, value) {
  if(!check_data_present(x))
    cli::cli_abort(message = "No analysis data loaded. Please first import raw data.")
  if(!(is.data.frame(value)))
    cli::cli_abort(message = "`metadata` must be a data.frame or tibble.")
  if(!"analysis_id" %in% names(value))
    cli::cli_abort(message = "`metadata` must contain a column `analysis_id` defining the ")
  if(!"curve_id" %in% names(value))
    cli::cli_abort(message = "`metadata` must contain a column `curve_id` defining response curve series")
  if(!(any(c("analyzed_amount", "injection_volume") %in% names(value))))
    cli::cli_abort(message = "`metadata` must contain a column `analyzed_amount` or `injection_volume`")
  if(!all(value$analysis_id %in% x@annot_analyses$analysis_id))
    cli::cli_abort(message = "One or more analysis are not present in the analysis data. Please ensure all `analysis_id` are present in the analysis data.")
  if(any(is.na(value$curve_id)))
    cli::cli_abort(message = "One or `more curve_id` is not defined. Please check your data.")

  if(any(is.na(value$curve_id)) &
     (ifelse("analyzed_amount" %in% names(value), all(is.numeric(value$analyzed_amount)), FALSE) |
      ifelse("injection_volume" %in% names(value), all(is.numeric(value$injection_volume)), FALSE)))

    cli::cli_abort(message = "One or more `analyzed_amount` or `injection_volume` are not defined. Please ensure completness of at least one of the variables.")

  x@annot_responsecurves <- as_tibble(value)
  validObject(x)
  x
})

