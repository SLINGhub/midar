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
#' @param new_feature_name Optional. New name of corrected feature. If empty then feature name will not change.
#' @param variable Default: `feature_intensity`. Name of Variable to be corrected.
#' @return MidarExperiment object
#' @export

#  Example:  mexp <- correct_interference_manually(mexp, "feature_intensity", "PC 32:0 | SM 36:1 M+3", "SM 36:1", 0.0106924, "PC 32:0")

correct_interference_manually <- function(data, variable, feature, interfering_feature, relative_contribution, new_feature_name = NULL)  {

  variable_var <- rlang::ensym(variable)

  new_feature_name <- ifelse(is.null(new_feature_name) | is.na(new_feature_name), "" , new_feature_name)

    # Validate input
  if(!feature %in% data@annot_features$feature_name) stop("Feature is not present in the dataset")
  if(!interfering_feature %in% data@annot_features$feature_name) stop("Interfering feature is not present in the dataset")
  if(!variable %in% names(data@dataset)) stop(glue::glue("Variable `{variable` is not defined in the dataset"))
  if(relative_contribution <0 | relative_contribution >= 1) stop("`relative_contribution` must be between 0 and 1")
  if(new_feature_name %in% data@annot_features$feature_name) stop("Mew fFeature name must not present already be present in the dataset")

  browser()

  #Correction
  data@dataset <- data@dataset |>
    group_by(.data$analysis_id) |>
    mutate(
      !!variable_var := if_else(.data$feature_name == feature,
                                (!!variable_var)[.data$feature_name == feature] - relative_contribution * (!!variable_var)[.data$feature_name == interfering_feature],
                             !!variable_var),
      corrected_interference = if_else(.data$feature_name == feature, TRUE, .data$corrected_interference)
      )

  data@dataset <- data@dataset |>
    mutate(feature_name = if_else(.data$feature_name == feature & nchar(stringr::str_squish(.data$new_feature_name)) > 0, new_feature_name, .data$feature_name))

  data
}


#' Substract interferences contributed by another feature
#'
#' @description
#' The interference (e.g. isotope overlap or in-source fragments) is substracted using following formula:
#' \deqn{Value_corrected = Value_Feature - relative_contribution * Value_InterferingFeature}
#'
#  The intefering features an their relative contributions to the interferences are defined in the feature annotation table
#' @param data MidarExperiment object
#' @param variable Name of Variable to be corrected. Default: `feature_intensity`.
#' @return MidarExperiment object
#' @export

correct_interferences <- function(data, variable = "feature_intensity")  {

  if (data@is_isotope_corr) crayon::yellow(glue::glue("Note: Data is already corrected for interferences. Correction will be reapplied to raw data."))
  if(variable != "feature_intensity") stop("Currently only correction for raw intensities suspported, thus must be set to `feature_intensity` or not defined.")

  #variable_var <- rlang::ensym(variable)

  if (!"feature_intensity_raw" %in% names(data@dataset))
    data@dataset <- data@dataset |>
      mutate(feature_intensity_raw = .data$feature_intensity, .before = "feature_intensity")

  #ToDo: check if interering feature is in the feature list
  d_corrected_features <- data@dataset |>
    mutate(is_interfering = .data$feature_name %in% data@annot_features$interference_feature_name,
           interference_group = if_else(.data$is_interfering, .data$feature_name, .data$interference_feature_name)) |>
    filter(!is.na(.data$interference_group)) |>
    arrange(desc(.data$interference_feature_name)) |>
    group_by(.data$analysis_id, .data$interference_group) |>
    mutate(value_corr_interf = if_else(!.data$is_interfering & !is.na(.data$interference_feature_name), .data$feature_intensity_raw - .data$feature_intensity_raw[.data$is_interfering] * .data$interference_proportion, NA_real_),
           corrected_interference = TRUE) |>
    select("analysis_id", "feature_name", "interference_feature_name", "is_interfering", "interference_group", "feature_intensity", "interference_proportion", "value_corr_interf", "corrected_interference") |>
    filter(!.data$is_interfering) |>
    ungroup()

  data@dataset <- data@dataset |>
    full_join(d_corrected_features |> select(c("analysis_id", "feature_name","value_corr_interf", "corrected_interference")), by = c("analysis_id", "feature_name")) |>
    mutate(corrected_interference = dplyr::coalesce(.data$corrected_interference.x, .data$corrected_interference.y)) |>
    mutate(feature_intensity = if_else(.data$corrected_interference, .data$value_corr_interf, .data$feature_intensity_raw)) |>
    select(!"value_corr_interf", !"corrected_interference.x", !"corrected_interference.y")

  data@is_isotope_corr <- TRUE
  n_corr <- data@annot_features |> filter(!is.na(.data$interference_feature_name)) |> nrow()
  writeLines(crayon::green(glue::glue("\u2713 Interference/isotope corection has been applied to {n_corr} of {nrow(data@annot_features)} features.")))

  data
  }
