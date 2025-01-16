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
#' @param updated_feature_id Optional. New name of corrected feature. If empty then feature name will not change.
#' @param variable Default: `feature_intensity`. Name of Variable to be corrected.
#' @return MidarExperiment object
#' @noRd

#  Example:  mexp <- correct_interference_manually(mexp, "feature_intensity", "PC 32:0 | SM 36:1 M+3", "SM 36:1", 0.0106924, "PC 32:0")

correct_interference_manually <- function(data = NULL, variable, feature, interfering_feature, relative_contribution, updated_feature_id = NULL) {

  check_data(data)

  variable_var <- rlang::ensym(variable)

  updated_feature_id <- ifelse(is.null(updated_feature_id) | is.na(updated_feature_id), "", updated_feature_id)

  # Validate input
  if (!feature %in% data@annot_features$feature_id) cli::cli_abort("Feature is not present in the dataset")
  if (!interfering_feature %in% data@annot_features$feature_id) cli::cli_abort("Interfering feature is not present in the dataset")
  if (!variable %in% names(data@dataset)) stop(glue::glue("Variable `{variable` is not defined in the dataset"))
  if (relative_contribution < 0 | relative_contribution >= 1) cli::cli_abort("`relative_contribution` must be between 0 and 1")
  if (updated_feature_id %in% data@annot_features$feature_id) cli::cli_abort("Mew fFeature name must not present already be present in the dataset")

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
    mutate(feature_id = if_else(.data$feature_id == feature & nchar(stringr::str_squish(.data$updated_feature_id)) > 0, updated_feature_id, .data$feature_id))

  data
}


#' Subtract interference contributed by another feature
#'
#' @description
#' The interference (e.g. isotope overlap or in-source fragments) is subtracted using following formula:
#' \deqn{value_corrected = value_raw - value_raw_interfering_feature * proportion_interference}
#'
#  The interfering features an their relative contributions to the interference are defined in the feature annotation table
#' @param data MidarExperiment object
#' @param variable Name of Variable to be corrected. Default: `feature_intensity`.
#' @return MidarExperiment object with feature intensities corrected for interferences
#' @export

correct_interferences <- function(data = NULL, variable = "feature_intensity") {

  check_data(data)

  if (variable != "feature_intensity") cli::cli_abort("Currently only correction for raw intensities suspported, thus must be set to `feature_intensity` or not defined.")

  if (data@is_isotope_corr & (c("feature_intensity_orig") %in% names(data@dataset))) {
    cli_alert_info(glue::glue("Data is already interference-corrected. Corrections will be reapplied to raw intensities."))
    data@dataset <- data@dataset |>
      mutate(
        feature_intensity = .data$feature_intensity_orig,
        interference_corrected = FALSE)
  } else {
      data@dataset <- data@dataset |>
        mutate(feature_intensity_orig = .data$feature_intensity,
               interference_corrected = FALSE,
               .before = "feature_intensity")
  }

  ############

  # Initial data processing
  d_correct <- data@dataset |>
    left_join(
      data@annot_features |> select("feature_id", "interference_feature_id", "interference_proportion"),
      by = "feature_id") |>
    # mutate(
    #   is_interfering = .data$feature_id %in% data@annot_features$interference_feature_id,
    #   is_interfered = !is.na(.data$interference_feature_id)
    # ) |>
    select("analysis_id", "feature_id", "feature_intensity_orig", "interference_feature_id", "interference_proportion")

  features_to_correct <- d_correct |>
    filter(!is.na(.data$interference_feature_id)) |>
    select("feature_id", "interference_feature_id") |>
    distinct() |>
    arrange(rev(.data$feature_id))


  # Function to apply correction for each feature set in features_to_correct
  correct_feature_intensity <- function(data,features_to_correct, i) {
    d <- data %>%
      group_by(.data$analysis_id) |>
      mutate(
        raw_target = if_else(.data$feature_id == features_to_correct$feature_id[i],
                    .data$feature_intensity_orig[.data$feature_id == features_to_correct$feature_id[i]], NA_real_),
        raw_source = if_else(.data$feature_id == features_to_correct$feature_id[i],
                    .data$feature_intensity_orig[.data$feature_id == features_to_correct$interference_feature_id[i]], NA_real_),
        rel_interference = if_else(.data$feature_id == features_to_correct$feature_id[i],
                    .data$interference_proportion[.data$feature_id == features_to_correct$feature_id[i]], NA_real_),
        corr_intensity = if_else(.data$feature_id == features_to_correct$feature_id[i],
                                 raw_target - raw_source * rel_interference, NA_real_)
      ) |>
      ungroup() |>
      select(-"raw_target", -"raw_source", -"rel_interference") |>
      mutate(feature_intensity_orig = if_else(!is.na(.data$corr_intensity), .data$corr_intensity, .data$feature_intensity_orig))
    d
  }

  # Initialize corrected data with original data
  d_corrected <- d_correct

  # Apply the correction function iteratively, using the result from the previous iteration
  for (i in 1:nrow(features_to_correct)) {
    d_corrected <- correct_feature_intensity(d_corrected, features_to_correct, i)
  }

  ############

  data@dataset <- data@dataset |>
    left_join(d_corrected |> select("analysis_id", "feature_id", intensity_corrected  = "feature_intensity_orig"), by = c("analysis_id", "feature_id")) |>
    mutate(interference_corrected = .data$feature_intensity != .data$intensity_corrected) |>
    mutate(feature_intensity = .data$intensity_corrected) |>
    select(-"intensity_corrected", -"feature_intensity_orig")


  data@is_isotope_corr <- TRUE
  data@status_processing <- "Isotope-corrected raw data"
  data <- update_after_normalization(data, FALSE)
  data@var_drift_corrected <- c(feature_intensity = FALSE, feature_norm_intensity = FALSE, feature_conc = FALSE)
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]

  n_corr <- nrow(features_to_correct)
  cli_alert_success(col_green(glue::glue("Interference-correction has been applied to {n_corr} of the {get_feature_count(data)} features.")))

  data
}
