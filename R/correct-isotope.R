#' Manual correction for interference contributed by another feature
#'
#' @description
#' The interference is substracted using following formula:
#' \deqn{Value_{Corrected} = Value_{Feature} - Factor_{Contribution} * Value_{Interfering Feature}}
#'
#' @param data MidarExperiment object
#' @param feature Name of feature to be corrected
#' @param interfering_feature Name of feature that is interfering, i.e. contributing to the signal of `feature`
#' @param relative_contribution Relative portion of the interfering feature to contribute to the feature signal. Must be between 0 and 1.
#' @param new_feature_id Optional. New name of corrected feature. If empty then feature name will not change.
#' @param variable Default: `feature_intensity`. Name of Variable to be corrected.
#' @return MidarExperiment object


#  Example:  mexp <- correct_interference_manually(mexp, "feature_intensity", "PC 32:0 | SM 36:1 M+3", "SM 36:1", 0.0106924, "PC 32:0")

correct_interference_manually <- function(data, variable, feature, interfering_feature, relative_contribution, new_feature_id = NULL) {
  variable_var <- rlang::ensym(variable)

  new_feature_id <- ifelse(is.null(new_feature_id) | is.na(new_feature_id), "", new_feature_id)

  # Validate input
  if (!feature %in% data@annot_features$feature_id) cli::cli_abort("Feature is not present in the dataset")
  if (!interfering_feature %in% data@annot_features$feature_id) cli::cli_abort("Interfering feature is not present in the dataset")
  if (!variable %in% names(data@dataset)) stop(glue::glue("Variable `{variable` is not defined in the dataset"))
  if (relative_contribution < 0 | relative_contribution >= 1) cli::cli_abort("`relative_contribution` must be between 0 and 1")
  if (new_feature_id %in% data@annot_features$feature_id) cli::cli_abort("Mew fFeature name must not present already be present in the dataset")

  browser()

  # Correction
  data@dataset <- data@dataset |>
    group_by(.data$analysis_id) |>
    mutate(
      !!variable_var := if_else(.data$feature_id == feature,
        (!!variable_var)[.data$feature_id == feature] - relative_contribution * (!!variable_var)[.data$feature_id == interfering_feature],
        !!variable_var
      ),
      interference_corrected = if_else(.data$feature_id == feature, TRUE, .data$interference_corrected)
    )

  data@dataset <- data@dataset |>
    mutate(feature_id = if_else(.data$feature_id == feature & nchar(stringr::str_squish(.data$new_feature_id)) > 0, new_feature_id, .data$feature_id))

  data
}


#' Substract interferences contributed by another feature
#'
#' @description
#' The interference (e.g. isotope overlap or in-source fragments) is subtracted using following formula:
#' \deqn{Value_corrected = Value_Feature - relative_contribution * Value_InterferingFeature}
#'
#  The interfering features an their relative contributions to the interference are defined in the feature annotation table
#' @param data MidarExperiment object
#' @param variable Name of Variable to be corrected. Default: `feature_intensity`.
#' @return MidarExperiment object
#' @export

correct_interferences <- function(data, variable = "feature_intensity") {
  if (variable != "feature_intensity") cli::cli_abort("Currently only correction for raw intensities suspported, thus must be set to `feature_intensity` or not defined.")

  if (data@is_isotope_corr & (c("feature_intensity_raw") %in% names(data@dataset))) {
    cli_alert_info(glue::glue("Data is already interference-corrected. Correction will be reapplied to raw intensities."))
    data@dataset <- data@dataset |>
      mutate(
        feature_intensity = .data$feature_intensity_raw,
        interference_corrected = FALSE)
  } else {
      data@dataset <- data@dataset |>
        mutate(feature_intensity_raw = .data$feature_intensity,
               interference_corrected = FALSE,
               .before = "feature_intensity")
  }


  # ToDo: check if interering feature is in the feature list
  d_corrected_features <- data@dataset |>
    left_join(data@annot_features |> select("feature_id", "interference_feature_id", "interference_proportion"), by = "feature_id") |>
    mutate(
      is_interfering = .data$feature_id %in% data@annot_features$interference_feature_id,
      interference_group = if_else(.data$is_interfering, .data$feature_id, .data$interference_feature_id)
    ) |>
    filter(!is.na(.data$interference_group)) |>
    arrange(desc(.data$interference_feature_id)) |>
    group_by(.data$analysis_id, .data$interference_group) |>
    mutate(
      corr_intensity = if_else(!.data$is_interfering & !is.na(.data$interference_feature_id), .data$feature_intensity_raw - .data$feature_intensity_raw[.data$is_interfering] * .data$interference_proportion, NA_real_),
      interference_corrected_temp = TRUE
    ) |>
    filter(!.data$is_interfering) |>
    ungroup() |>
    select("analysis_id", "feature_id", "corr_intensity", "interference_corrected_temp")

    data@dataset <- data@dataset |>
    left_join(d_corrected_features, by = c("analysis_id", "feature_id")) |>
    mutate(interference_corrected = dplyr::coalesce(.data$interference_corrected, .data$interference_corrected_temp)) |>
    mutate(feature_intensity = if_else(.data$interference_corrected, .data$corr_intensity, .data$feature_intensity_raw)) |>
    select(-"corr_intensity", -"interference_corrected_temp")

  data@is_isotope_corr <- TRUE
  data@status_processing <- "Isotope-corrected raw data"
  data@is_istd_normalized <- FALSE
  data@is_quantitated <- FALSE
  data@is_drift_corrected <- FALSE
  data@is_batch_corrected <- FALSE
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]

  n_corr <- data@annot_features |>
    filter(.data$valid_feature) |>
    filter(!is.na(.data$interference_feature_id)) |>
    nrow()
  cli_alert_success(col_green(glue::glue("Interference-correction has been applied to {n_corr} of the {get_feature_count(data)} features.")))

  data
}
