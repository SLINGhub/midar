

# TODO: export functions
check_data_present <- function(data){nrow(data@dataset_orig) > 0}
check_dataset_present <- function(data){nrow(data@dataset) > 0}


#' @title Get the annotated or the originally imported analytical data
#' @param data MidarExperiment object
#' @param annotated Boolean indicating whether to return the annotated data
#' (`FALSE`) or the original imported data (`TRUE`)
#' @return A tibble with the analytical data in the long format
#' @export

get_analyticaldata <- function(data = NULL, annotated ){
  check_data(data)
  if (!annotated) {
    return(data@dataset_orig)
  } else {
    return(data@dataset)
  }
}


#' Get the number of analyses in the dataset
#'
#' Returns the number of analyses in the dataset, with an optional
#' filter based on `qc_types`.
#'
#' @param data A `MidarExperiment` object
#' @param qc_types Defines the `qc_types` to be counted. If `NULL` or `NA`,
#' all analyses will be counted.
#'
#' @return An integer with the analysis count
#'
#' @export
get_analysis_count <- function(data, qc_types = NULL) {
  if (nrow(data@dataset) == 0) {
    return(0)
  }
  if (is.null(qc_types))
    return(data@dataset|> select("analysis_id") |> distinct() |> nrow())
  else
    return(data@dataset|> filter(.data$qc_type %in% qc_types) |>
             select("analysis_id") |> distinct() |> nrow())
}


#' Get the number of features in the dataset
#'
#' Returns the number of features in the dataset, with optional
#' filters whether counted features must be internal standard and/or quantifier.
#'
#' @param data A `MidarExperiment` object
#' @param is_istd If set, then defines whether to include or exclude internal standard features. Default is `NA` means no filter for internal standards is applied.
#' @param is_quantifier If set, then defines whether to include or exclude qualifier features. Default is `NA` means no filter for qualifier features is applied.
#'
#' @return An integer with the feature count
#'
#' @export
get_feature_count <- function(data, is_istd = NA, is_quantifier = NA) {
  if (nrow(data@dataset) == 0) {
    return(0)
  }
  d <- data@dataset
  if (!is.na(is_istd)) d <- d |> filter(.data$is_istd == !!is_istd)
  if (!is.na(is_quantifier)) d <- d |> filter(.data$is_quantifier == !!is_quantifier)

  d |> select("feature_id") |> distinct() |> nrow()
}

#' Get feature IDs
#'
#' Returns a vector of annotated feature IDs (`feature_id`) present in the dataset
#'
#' @param data A `MidarExperiment` object
#' @param is_istd If set, then defines whether to include or exclude internal standard features. Default is `NA` means no filter for internal standards is applied.
#' @param is_quantifier If set, then defines whether to include or exclude qualifier features. Default is `NA` means no filter for qualifier features is applied.
#'
#' @return A character vector with `feature_id` values
#'
#' @export

get_featurelist <- function(data, is_istd = NA, is_quantifier = NA) {
  d <- data@dataset
  if (nrow(data@dataset) == 0) {
    return(NULL)
  }
  if (!is.na(is_istd)) d <- d |> filter(.data$is_istd == !!is_istd)
  if (!is.na(is_quantifier)) d <- d |> filter(.data$is_quantifier == !!is_quantifier)

  d |> select("feature_id") |> distinct() |> pull(.data$feature_id)
}

#' Get the start time of the analysis sequence
#'
#' Returns the start time of the analysis, corresponding to the earliest `acquisition_time_stamp` from the dataset.
#'
#' @param data A `MidarExperiment` object
#' @return A `POSIXct` timestamp, or `NA_POSIXct_` if the dataset is empty.
#'
#' @export
get_analyis_start <- function(data){
  if (check_data_present(data))
    return(min(data@dataset$acquisition_time_stamp))
  else
    return(lubridate::NA_POSIXct_)
}

#' Get the end time of the analysis sequence
#'
#' Returns the end time of the analysis, corresponding to the last `acquisition_time_stamp` from the dataset.
#' Note: if `estimate_analysis_end` is set to `FALSE`, the function will return
#' the timestamp of the last analysis in the dataset, corresping to the start of
#' the last analysis. Set `estimate_analysis_end` to `TRUE` to estimate the end
#' time of the analysis sequence, based on the median runtime.
#'
#' @param data A `MidarExperiment` object
#' @param estimate_sequence_end If `TRUE`, the function will estimate the end
#' time of the analysis sequence based on the median runtime. `FALSE` will
#' return to start time of last analysis in the sequence.
#' @return A `POSIXct` timestamp, or `NA_POSIXct_` if the dataset is empty.
#'
#' @export
get_analyis_end <- function(data, estimate_sequence_end){
  if (check_data_present(data)){
    if(estimate_sequence_end)
      return(max(data@dataset$acquisition_time_stamp) + get_runtime_median(data))
    else
      return(max(data@dataset$acquisition_time_stamp))
  }

  else
    return(lubridate::NA_POSIXct_)
}

#' Get the median run time
#'
#' Calculates the median run time (in seconds) based of the timestamps differences between consecutive analyses in the sequence.
#'
#' @param data A `MidarExperiment` object
#' @return A `lubridate` time period object, or `NA` if the dataset is empty.
#'
#' @export

get_runtime_median <- function(data){
  if (check_data_present(data))
    median(diff(unique(data@dataset$acquisition_time_stamp), units = "secs")) |>  lubridate::seconds_to_period()
  else
    return(NA)
}

#' Get the total duration of the analysis
#'
#' This function returns the total duration of the analysis, which is the time difference
#' between the timestamps of the first and last analyses in the sequence.
#'
#' If `estimate_sequence_end` is `TRUE`, the function will estimate the end time of the analysis sequence
#' by adding the median runtime to the timestamp of the last analysis, instead of simply using the timestamp
#' of the last analysis.
#'
#' @param data A `MidarExperiment` object
#' @param estimate_sequence_end If `TRUE`, the function will estimate the end
#' time of the analysis sequence based on the median runtime, added to the timestamp of the last analysis.
#' If `FALSE`, the function will calculate the time difference between the first and last analysis timestamps
#' without any adjustment.
#'
#' @export

get_analysis_duration <- function(data, estimate_sequence_end){
  if (check_data_present(data)){
    time_end <- max(unique(data@dataset$acquisition_time_stamp))
    if(estimate_sequence_end) time_end <- time_end + get_runtime_median(data)

    difftime(time_end,
             min(unique(data@dataset$acquisition_time_stamp)), units = "secs") |>
      lubridate::seconds_to_period()
  } else {
    return(NA)
  }
}

#' Get the number of analysis breaks in the analysis
#'
#' Counts the number of interruptions in the analysis, where an interruption is
#' defined as a time gap between consecutive acquisition timestamps that
#' exceeds a given threshold (`break_duration_minutes`).

#' @param data A `MidarExperiment` object
#' @param break_duration_minutes A numeric value specifying the minimum duration
#' (in minutes) between two consecutive analyses that qualifies as an interruption.
#'
#' @return An integer with the number of interruptions, or `NA_integer_` if the dataset is empty.
#'
#' @export
get_analysis_breaks <- function(data, break_duration_minutes){
  if (check_data_present(data)) {
    if(all(is.na(unique(data@dataset$acquisition_time_stamp)))) return(NA_integer_)
    as.integer(sum(diff(unique(data@dataset$acquisition_time_stamp), units = "secs") > break_duration_minutes * 60))
  } else {
    return(NA_integer_)
  }
}


# sets is_normalized flag, if FALSE remove normalized intensities and conc if availalble from the dataset
# is_normalized defines if the data is normalized or not, which is to be set
update_after_normalization <- function(data, is_normalized, with_message = TRUE){
  if(data@is_istd_normalized & !is_normalized) {
    data@dataset <- data@dataset |>
      select(-any_of(c("feature_norm_intensity", "feature_norm_intensity_raw")))

    if(data@is_quantitated){
      data <- update_after_quantitation(data, FALSE, FALSE)
      if(with_message)
        cli_alert_info(cli::col_yellow("The normalized intensities and concentrations have been invalidated. Please reprocess the data."))
    } else {
      if(with_message)
        cli_alert_info(cli::col_yellow("Normalized intensities have been invalidated. Please reprocess data."))
    }
  }
  data@is_istd_normalized <- is_normalized
  data
}

update_after_quantitation <- function(data, is_quantitated, with_message = TRUE){
  if(data@is_quantitated & !is_quantitated) {
    data@dataset <- data@dataset |> select(-any_of(c("feature_conc", "feature_raw_conc")))
    if(with_message) cli_alert_info(cli::col_yellow("Concentrations not valid anymore. Please reprocess data."))
  }
  data@is_quantitated <- is_quantitated
  data
}

check_var_in_dataset <- function(table, variable) {

  if(variable == "feature_conc" & !"feature_conc" %in% names(table)) cli_abort(cli::col_red("Concentrations not available, please process data or choose another variable.", show = "none", parent = NULL, call= NULL))
  if(variable == "feature_area" & !"feature_area" %in% names(table)) cli_abort(cli::col_red("Area is not available, please choose another variable.", show = "none", parent = NULL, call= NULL))
  if(variable == "feature_response" & !"feature_response" %in% names(table)) cli_abort(cli::col_red("Response is not available, please choose another variable.", show = "none", parent = NULL, call= NULL))
  if(variable == "feature_norm_intensity" & !"feature_norm_intensity" %in% names(table)) cli_abort(cli::col_red("Normalized intensities not available, please process data, or choose another variable.", show = "none", parent = NULL, call= NULL))
  if(variable == "feature_height" & !"feature_height" %in% names(table)) cli_abort(cli::col_red("Height is not available, please choose another variable.", show = "none", parent = NULL, call= NULL))
  if(variable == "feature_conc_raw" & !"feature_conc_raw" %in% names(table)) cli_abort(cli::col_red("Concentrations not available, please process data, or choose another variable.", show = "none", parent = NULL, call= NULL))
}

#' @title Get the start and end analysis numbers of specified batches
#' @description
#' Sets the analysis order (sequence) based on either (i) analysis timestamp, if available, (ii) the order in which analysis appeared in the imported raw data file, or (iii) the order in which analyses were defined in the Analysis metadata.
#' @param data MidarExperiment object
#' @param batch_ids A vector with one or two elements: the first and/or last batch ID. If NULL or invalid, the function will abort.
#' @return A vector with two elements: the lower and upper analysis number for the specified batch(es).
#' @export
get_batch_boundaries <- function(data = NULL, batch_ids = NULL) {
  check_data(data)
  if (nrow(data@annot_batches) == 0) {
    cli::cli_abort("No batches defined in the dataset.")
  }

  if (is.null(batch_ids) || length(batch_ids) == 0) {
    cli::cli_abort("No batch IDs provided.")
  }
  # Ensure batch_ids has at least one value
  if (length(batch_ids) == 1) {
    first_batch_id <- last_batch_id <- batch_ids[1]
  } else if (length(batch_ids) == 2) {
    first_batch_id <- batch_ids[1]
    last_batch_id <- batch_ids[2]
  } else {
    cli::cli_abort("Please provide a vector with one or two batch IDs.")
  }

  # Check if the specified batch IDs exist in the dataset
  existing_batches <- data@annot_batches$batch_id
  if (!all(c(first_batch_id, last_batch_id) %in% existing_batches)) {
    cli::cli_abort("One or more of the specified batch IDs do not exist in the dataset.")
  }

  # Filter the dataset for the given range of batch IDs
  d <- data@annot_batches |>
    filter(.data$batch_id >= first_batch_id & .data$batch_id <= last_batch_id)

  if (nrow(d) == 0) {
    cli::cli_abort("No batches found for the provided range of batch IDs.")
  }

  # Extract lower and upper analysis numbers
  lower_bound <- min(d$id_batch_start)
  upper_bound <- max(d$id_batch_end)

  return(c(lower_bound, upper_bound))
}




#' @title Set the analysis order
#' @description
#' Sets the analysis order (sequence), based on either (i) analysis timestamp if available, (ii) the order in which analysis appeared in the imported raw data file, or (iii) the order in which analyses were defined in the Analysis metadata
#' @param data MidarExperiment object
#' @param order_by Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
#' @return MidarExperiment object
#' @examples
#' file_path <- system.file("extdata", "MRMkit_demo.tsv", package = "midar")
#' mexp <- MidarExperiment()
#' mexp <- import_data_mrmkit(mexp, path = file_path, import_metadata = TRUE)
#' mexp <- set_analysis_order(mexp, "timestamp")

#' Internal function to set analysis order in metadata
#'
#' @param data A MidarExperiment object
#' @param order_by Character string specifying ordering method
#' @noRd
set_analysis_order_analysismetadata <- function(data = NULL,
                                                order_by = "default") {
  # Validate inputs
  check_data(data)
  order_by <- rlang::arg_match(
    arg = order_by,
    values = c("timestamp", "resultfile", "metadata", "default"),
    multiple = FALSE
  )

  # Extract relevant columns from original dataset
  d_temp <- data@dataset_orig |>
    dplyr::select("analysis_id", dplyr::any_of(c("run_seq_num", "acquisition_time_stamp"))) |>
    dplyr::distinct()

  # Handle default ordering logic
  if (order_by == "default") {
    order_by <- if (!all(is.na(d_temp$acquisition_time_stamp))) {
      "timestamp"
    } else {
      cli::cli_alert_info(
        cli::col_grey(
          "Analysis order was based on sequence of analysis results, as no timestamps were found.
           Use `set_analysis_order` to define alternative analysis orders."
        )
      )
      "resultfile"
    }
  }

  # Apply ordering based on specified method
  data@annot_analyses <- data@annot_analyses |>
    select(-"run_seq_num")

  data@annot_analyses <- switch(
    order_by,
    "timestamp" = {
      if (all(is.na(d_temp$acquisition_time_stamp))) {
        cli::cli_abort(col_red(
          "Acquisition timestamps are not present in analysis results.
           Please set argument `order_by` to either `resultfile` or `metadata`."),
          call. = FALSE
        )
      }
      data@annot_analyses |>
      dplyr::inner_join(
          d_temp |>
          dplyr::arrange(.data$acquisition_time_stamp) |>
          dplyr::mutate(run_seq_num = dplyr::row_number(), .before = 1) |>
          dplyr::select(-"acquisition_time_stamp"),
          by = "analysis_id") |>
          relocate("run_seq_num", .before = 1)
    },
    "resultfile" = {
      data@annot_analyses |>
        dplyr::inner_join(
          d_temp |>
            dplyr::mutate(run_seq_num = dplyr::row_number(), .before = 1),
        by = "analysis_id") |>
        relocate("run_seq_num", .before = 1)
    },
    "metadata" = {
      data@annot_analyses |>
        dplyr::mutate("run_seq_num" = .data$annot_order_num, .before = 1)
    }
  )

  # Clean up final dataset
  data@annot_analyses <- data@annot_analyses |>
    dplyr::select(-dplyr::any_of("acquisition_time_stamp"))

  data
}

#' Set Analysis Order
#' @description
#' Determines the sequence of analyses using either instrument timestamps,
#' the order in the imported raw data file, or the order defined in the Analysis metadata.
#' Note: After changing the analysis order, all post processing steps must be rerun.
#'
#' @param data A `MidarExperiment` object
#' @param order_by Character string specifying the ordering method.
#'   Must be one of "timestamp" (requires timestamp data in imported results),
#'   "resultfile" (uses order from imported data file), or
#'   "metadata" (uses order from analysis metadata)
#' @return An updated `MidarExperiment` object with ordered analyses
#' @export
#'
#' @examples
#' file_path <- system.file("extdata", "MRMkit_demo.tsv", package = "midar")
#' mexp <- MidarExperiment()
#' mexp <- import_data_mrmkit(mexp, path = file_path, import_metadata = TRUE)
#'
#' # Order by timestamp (if available)
#' mexp <- set_analysis_order(mexp, "timestamp")
#'
#' # Order by metadata definition
#' mexp <- set_analysis_order(mexp, "metadata")


set_analysis_order <- function(data = NULL, order_by =  c("timestamp", "resultfile", "metadata")){
  check_data(data)
  order_by <- rlang::arg_match(arg = order_by, c("timestamp", "resultfile", "metadata"), multiple = FALSE)
  data <- set_analysis_order_analysismetadata(data, order_by)
  data <- link_data_metadata(data)

  cli::cli_alert_success(cli::col_green("Analysis order set to {.val {order_by}}"))

  if (data@is_isotope_corr | data@is_filtered | data@is_istd_normalized | data@is_quantitated | any(data@var_batch_corrected) | any(data@var_drift_corrected))
  cli::cli_alert_info(col_yellow(c(
    "All data processing has been reset. ",
    "i" = "Please rerun processing steps"
  )))
  data
}






# Link DATA with METADATA and create DATASET table. =================
## - Only valid analyses and features will be added
## - Only key information will be added
link_data_metadata <- function(data = NULL, minimal_info = TRUE){
  check_data(data)
  data@dataset <- data@dataset_orig |>
    select(
      "run_seq_num",
      "analysis_id",
      "acquisition_time_stamp",
      "feature_id",
      starts_with("method_"),
      starts_with("feature_")
    )

  if (nrow(data@annot_analyses) > 0) {
  data@dataset <- data@dataset |>
    select(-any_of("run_seq_num")) |>
    inner_join(data@annot_analyses, by = "analysis_id") |>
    filter(.data$valid_analysis)
  }

  if (nrow(data@annot_features) > 0) {
    data@dataset <- data@dataset |>
      inner_join(data@annot_features, by = "feature_id") |>
      filter(.data$valid_feature)
  }

  data@dataset <- dplyr::bind_rows(pkg.env$table_templates$dataset_template, data@dataset)
  data@dataset <- data@dataset |>
    select(any_of(c(
      "run_seq_num",
      "analysis_id",
      "acquisition_time_stamp",
      "qc_type",
      "batch_id",
      "sample_id",
      "replicate_no",
      "feature_id",
      "feature_class",
      "is_istd",
      "is_quantifier",
      "specimen")),
      starts_with("method_"),
      starts_with("feature_"),
    ) |>
    relocate("feature_intensity", .after = dplyr::last_col()) |>
    relocate("feature_label", .after = "feature_class") |>
    relocate("sample_id", .after = "batch_id") |>
    relocate("specimen", .before = "feature_id") |>
    relocate("replicate_no", .after = "sample_id")

  if(minimal_info)
    data@dataset <- data@dataset |> select(-starts_with("method_"), -starts_with("feature_int_"))



  data@dataset <- data@dataset |> mutate(feature_intensity = !!(sym(data@feature_intensity_var)))


  data@is_isotope_corr <- FALSE
  #data@is_istd_normalized <- FALSE
  #data@is_quantitated <- FALSE
  data@var_drift_corrected <- c(feature_intensity = FALSE, feature_norm_intensity = FALSE, feature_conc = FALSE)
  data@var_drift_corrected <- c(feature_intensity = FALSE, feature_norm_intensity = FALSE, feature_conc = FALSE)
  data@is_filtered <- FALSE


  data@status_processing <- glue::glue("Annotated raw {toupper(str_remove(data@feature_intensity_var, 'feature_'))} values")

  data@metrics_qc <- data@metrics_qc[FALSE,]

  # Arrange run_seq_num and then by feature_id, as they appear in the metadata
  data@dataset <- data@dataset |>
    dplyr::arrange(match(.data$feature_id, data@annot_features$feature_id)) |>
    arrange(.data$run_seq_num)

  data
}

#' @title Set default variable to be used as feature raw signal value
#' @description
#' Sets the raw signal variable used for calculations starting from raw signal
#' values (i.e., normalization) Note that this set variable must be part of the
#' orginally imported data. Processed data variables (e.g., normalized intensities
#' and concentrations) can not be set as default feature intensity variable.
#' @param data MidarExperiment object
#' @param variable_name Feature variable to be used as default feature intensity for downstream processing.
#' @param auto_select If `TRUE` then the first available of these will be used as default: "intensity", "response", "area", "height".
#' @param warnings Suppress warnings
#' @param ... Feature variables to best search for one-by-one when `auto-detect = TRUE`
#' @return MidarExperiment object
#' @export

set_intensity_var <- function(data = NULL, variable_name, auto_select = FALSE, warnings = TRUE, ...){
  check_data(data)
  variable_strip <- str_remove(variable_name, "feature_")
  #rlang::arg_match(variable_strip, c("area", "height", "conc"))
  variable_name <- stringr::str_c("feature_", variable_strip)
  if (auto_select) {
    var_list <- unlist(rlang::list2(...), use.names = FALSE)
    id_all <- match(var_list,names(data@dataset_orig))
    idx <- which(!is.na(id_all))[1]
    if (!is.na(idx)) {
      data@feature_intensity_var = var_list[idx]
      cli_alert_info(text = cli::col_grey("{.var {var_list[idx]}} selected as default feature intensity. Modify with {.fn set_intensity_var}."))
      variable_name <- var_list[idx]
      } else {
      cli_alert_warning(text = cli::col_yellow("No typical feature intensity variable found in the data. Use {.fn set_intensity_var} to set it.}}."))
      return(data)
      }
  } else {
    #TODO: Double check behavior if there a feature_intensity in the raw data file
    if (! variable_name %in% names(data@dataset_orig))
      cli_abort(c("x" = "{.var {variable_name}} is not present in the raw data."))

    if (! variable_name %in% c("feature_intensity", "feature_response", "feature_area", "feature_height")){
      if(warnings) cli_alert_warning(cli::col_yellow("Note: {.var {variable_strip}} is not a typically used raw signal (i.e., area, height, intensity)."))
    }
      }
  data@feature_intensity_var <- variable_name

  if (check_dataset_present(data)) {
    calc_cols <- c("featue_norm_intensity", "feature_conc", "feature_amount", "feature_raw_conc")
    if (any(calc_cols %in% names(data@dataset))){
      data@dataset <- data@dataset |> select(-any_of(calc_cols))
      cli_alert_info(cli::col_yellow("New feature intensity variable (`{variable_name}`) defined, please reprocess data."))
    } else
    {
      cli::cli_alert_success(cli::col_green("Default feature intensity variable set to {.val {variable_name}}"))

    }
    data <- link_data_metadata(data)
  }

  if(variable_name == "feature_conc") data@is_quantitated <- TRUE
  data
}


#' @title Exclude analyses from the dataset
#'
#' @description
#' This function excludes specified analyses from a `MidarExperiment` object, either by
#' marking them as invalid for downstream processing.
#' The function also alloows to reset the exclusions.
#'
#' @param data A `MidarExperiment` object
#' @param analyses_to_exclude A character vector of analysis IDs (case-sensitive) to be excluded from the dataset.
#' If this is `NA` or an empty vector, the exclusion behavior will be handled as set via the `replace_existing` flag.
#' @param replace_existing A logical value. If `TRUE`, existing `valid_analysis` flags will be overwritten. If `FALSE`,
#' the exclusions will be appended, preserving any existing invalidated analyses.
#'
#' @return A modified `MidarExperiment` object with the specified analyses defined as excluded.
#' @export

exclude_analyses <- function(data = NULL, analyses_to_exclude, replace_existing ){
  check_data(data)
  if (all(is.na(analyses_to_exclude)) | length(analyses_to_exclude) == 0) {
    if(!replace_existing){
      cli_abort(cli::col_red("No `analysis_id` provided. To (re)include all analyses, use `analysis_ids_exlude = NA` and `replace_existing = TRUE`."))
    } else{
      cli::cli_alert_info(cli::col_green("All exclusions removed, and thus all analyses are now included for subsequent steps. Please reprocess data."))
      data@analyses_excluded <- NA
      data@annot_analyses <- data@annot_analyses |> mutate(valid_analysis = TRUE)
      data <- link_data_metadata(data)
      return(data)
      }
  }
  if (any(!c(analyses_to_exclude) %in% data@annot_analyses$analysis_id) > 0) {
    cli_abort(cli::col_red("One or more provided `analysis_id` to exclude are not present. Please check the analysis metadata."))
  }
  if(!replace_existing){
    data@annot_analyses <- data@annot_analyses |>
      mutate(valid_analysis = !(.data$analysis_id %in% analyses_to_exclude) & .data$valid_analysis)
    cli_alert_info(cli::col_green("A total of {data@annot_analyses |> filter(!.data$valid_analysis) |> nrow()} analyses are now excluded for downstream processing. Please reprocess data."))
    data@analyses_excluded <- data@annot_analyses |> filter(!.data$valid_analysis) |> pull(.data$analysis_id)
  } else {
    data@annot_analyses <- data@annot_analyses |>
      mutate(valid_analysis = !(.data$analysis_id %in% analyses_to_exclude))
    cli_alert_info(cli::col_green("{data@annot_analyses |> filter(!.data$valid_analysis) |> nrow()} analyses were excluded for downstream processing. Please reprocess data."))
    data@analyses_excluded <- data@annot_analyses |> filter(!.data$valid_analysis) |> pull(.data$analysis_id)
    }

  data <- link_data_metadata(data)

  data
}





#' @title Exclude features from the dataset
#'
#' @description
#' This function excludes specified features from a `MidarExperiment` object, either by
#' marking them as invalid for downstream processing.
#' The function also alloows to reset the exclusions.
#'
#' @param data A `MidarExperiment` object
#' @param features_to_exclude A character vector of feature IDs (case-sensitive) to be excluded from the dataset.
#' If this is `NA` or an empty vector, the exclusion behavior will be handled as set via the `replace_existing` flag.
#' @param replace_existing A logical value. If `TRUE`, existing `valid_analysis` flags will be overwritten. If `FALSE`,
#' the exclusions will be appended, preserving any existing invalidated features
#'
#' @return A modified `MidarExperiment` object with the specified analyses defined as excluded.
#' @export

exclude_features <- function(data = NULL, features_to_exclude, replace_existing ){
  check_data(data)

  if (all(is.na(features_to_exclude)) | length(features_to_exclude) == 0) {
    if(!replace_existing){
      cli_abort(cli::col_red("No `feature_id` provided. To (re)include all analyses, use `feature_ids_exlude = NA` and `replace_existing = TRUE`."))
    } else{
      cli::cli_alert_info(cli::col_green("All exlusions were removed, i.e. all features are included. Please reprocess data."))
      return(data)
    }
  }
  if (any(!c(features_to_exclude) %in% data@annot_features$feature_id) > 0) {
    cli_abort(cli::col_red("One or more provided `feature_id` are not present. Please check the feature metadata."))
  }
  if(!replace_existing){
    data@annot_features <- data@annot_features |>
      mutate(valid_feature = !(.data$feature_id %in% features_to_exclude) & .data$valid_feature)
    cli_alert_info(cli::col_green("A total of {data@annot_features |> filter(!.data$valid_feature) |> nrow()} features are now excluded for downstream processing. Please reprocess data."))
  }
  else {
    data@annot_features <- data@annot_features |>
      mutate(valid_feature = !(.data$feature_id %in% features_to_exclude))
    cli_alert_info(cli::col_green("{data@annot_features |> filter(!.data$valid_feature) |> nrow()} features were excluded for downstream processing. Please reprocess data."))
  }

  data <- link_data_metadata(data)

  data
}



