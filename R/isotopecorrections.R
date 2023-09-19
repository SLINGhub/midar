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
#' @param variable Default: `Intensity`. Name of Variable to be corrected.
#' @return MidarExperiment object
#' @export

#  Example:  mexp <- correct_interference_manually(mexp, "Intensity", "PC 32:0 | SM 36:1 M+3", "SM 36:1", 0.0106924, "PC 32:0")

correct_interference_manually <- function(data, variable, feature, interfering_feature, relative_contribution, new_feature_name = NULL)  {

  variable_var <- rlang::ensym(variable)

  new_feature_name <- ifelse(is.null(new_feature_name) | is.na(new_feature_name), "" , new_feature_name)

    # Validate input
  if(!feature %in% data@annot_features$FEATURE_NAME) stop("Feature is not present in the dataset")
  if(!interfering_feature %in% data@annot_features$FEATURE_NAME) stop("Interfering feature is not present in the dataset")
  if(!variable %in% names(data@dataset)) stop(glue::glue("Variable `{variable` is not defined in the dataset"))
  if(relative_contribution <0 | relative_contribution >= 1) stop("`relative_contribution` must be between 0 and 1")
  if(new_feature_name %in% data@annot_features$FEATURE_NAME) stop("Mew fFeature name must not present already be present in the dataset")

  browser()

  #Correction
  data@dataset <- data@dataset |>
    group_by(.data$ANALYSIS_ID) |>
    mutate(
      !!variable_var := if_else(.data$FEATURE_NAME == feature,
                                (!!variable_var)[.data$FEATURE_NAME == feature] - relative_contribution * (!!variable_var)[.data$FEATURE_NAME == interfering_feature],
                             !!variable_var),
      Corrected_Interference = if_else(.data$FEATURE_NAME == feature, TRUE, .data$Corrected_Interference)
      )

  data@dataset <- data@dataset |>
    mutate(FEATURE_NAME = if_else(.data$FEATURE_NAME == feature & nchar(stringr::str_squish(.data$new_feature_name)) > 0, new_feature_name, .data$FEATURE_NAME))

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
#' @param variable Name of Variable to be corrected. Default: `Intensity`.
#' @return MidarExperiment object
#' @export

correct_interferences <- function(data, variable = "Intensity")  {

  if (data@is_isotope_corr) crayon::yellow(glue::glue("Note: Data is already corrected for interferences. Correction will be reapplied to raw data."))
  if(variable != "Intensity") stop("Currently only correction for raw intensities suspported, thus must be set to `Intensity` or not defined.")

  #variable_var <- rlang::ensym(variable)

  if (!"Intensity_Raw" %in% names(data@dataset))
    data@dataset <- data@dataset |>
      mutate(Intensity_Raw = .data$Intensity, .before = "Intensity")

  #ToDo: check if interering feature is in the feature list
  d_corrected_features <- data@dataset |>
    mutate(is_interfering = .data$FEATURE_NAME %in% data@annot_features$INTERFERING_FEATURE,
           interference_group = if_else(.data$is_interfering, .data$FEATURE_NAME, .data$INTERFERING_FEATURE)) |>
    filter(!is.na(.data$interference_group)) |>
    arrange(desc(.data$INTERFERING_FEATURE)) |>
    group_by(.data$ANALYSIS_ID, .data$interference_group) |>
    mutate(value_corr_interf = if_else(!.data$is_interfering & !is.na(.data$INTERFERING_FEATURE), .data$Intensity_Raw - .data$Intensity_Raw[.data$is_interfering] * .data$INTERFERANCE_PROPORTION, NA_real_),
           Corrected_Interference = TRUE) |>
    select("ANALYSIS_ID", "FEATURE_NAME", "INTERFERING_FEATURE", "is_interfering", "interference_group", "Intensity", "INTERFERANCE_PROPORTION", "value_corr_interf", "Corrected_Interference") |>
    filter(!.data$is_interfering) |>
    ungroup()

  data@dataset <- data@dataset |>
    full_join(d_corrected_features |> select(c("ANALYSIS_ID", "FEATURE_NAME","value_corr_interf", "Corrected_Interference")), by = c("ANALYSIS_ID", "FEATURE_NAME")) |>
    mutate(Corrected_Interference = dplyr::coalesce(.data$Corrected_Interference.x, .data$Corrected_Interference.y)) |>
    mutate(Intensity = if_else(.data$Corrected_Interference, .data$value_corr_interf, .data$Intensity_Raw)) |>
    select(!"value_corr_interf", !"Corrected_Interference.x", !"Corrected_Interference.y")

  data@is_isotope_corr <- TRUE
  n_corr <- data@annot_features |> filter(!is.na(.data$INTERFERING_FEATURE)) |> nrow()
  writeLines(crayon::green(glue::glue("\u2713 Interference/isotope corection has been applied to {n_corr} of {nrow(data@annot_features)} features.")))

  data
  }
