#' Manual isotopic interference correction
#'
#' @description The interference is subtracted using following formula:
#' \deqn{Value_{Corrected} = Value_{Feature} - Factor_{Contribution} *
#' Value_{Interfering Feature}}
#'
#' @param data MidarExperiment object
#' @param feature Name of feature to be corrected
#' @param interfering_feature Name of feature that is interfering, i.e.
#'   contributing to the signal of `feature`
#' @param interference_contribution Relative portion of the interfering feature to
#'   contribute to the feature signal. Must be between 0 and 1.
#' @param updated_feature_id Optional. New name of corrected feature. If empty
#'   then feature name will not change.
#' @param variable Default: `feature_intensity`. Name of Variable to be
#'   corrected.
#' @return MidarExperiment object
#' @export


#  Example:  mexp <- correct_interference_manual(mexp, "feature_intensity", "PC 32:0 | SM 36:1 M+3", "SM 36:1", 0.0106924, "PC 32:0")

correct_interference_manual <- function(data = NULL, variable, feature, interfering_feature, interference_contribution, updated_feature_id = NA) {

  check_data(data)
  variable_var <- rlang::ensym(variable)

  updated_feature_id <- ifelse(is.null(updated_feature_id) | is.na(updated_feature_id), NA_character_, updated_feature_id)

  # Validate input
  if (is.na(feature) | !feature %in% data@annot_features$feature_id) cli::cli_abort(col_red("Selected feature is not present in the dataset. Please verify data and `feature` argument."))
  if (is.na(interfering_feature) | !interfering_feature %in% data@annot_features$feature_id) cli::cli_abort(col_red("Selected interfering feature is not present in the dataset. Please verify data and `feature` argument."))
  if (is.na(variable) |!variable %in% names(data@dataset)) cli::cli_abort(col_red("Variable `{variable}` is not defined in the dataset"))
  if (is.na(interference_contribution) | !is.numeric(interference_contribution) | interference_contribution <= 0 ) cli::cli_abort(col_red("`interference_contribution` must be a number larger than 0"))
  if (!is.na(updated_feature_id) & updated_feature_id %in% data@annot_features$feature_id) cli::cli_abort(col_red("Selected new feature id `{updated_feature_id}` is already present in the dataset. Please chose a new unique ID."))

  if (!"interference_corrected" %in% names(data@dataset)) {
    data@dataset <- data@dataset |>
      mutate(interference_corrected = FALSE)
  }

  # Correction
  data@dataset <- data@dataset |>
    group_by(.data$analysis_id) |>
    mutate(
      !!variable_var := if_else(.data$feature_id == feature,
        (!!variable_var)[.data$feature_id == feature] - interference_contribution * (!!variable_var)[.data$feature_id == interfering_feature],
        !!variable_var
      ),
      interference_corrected = if_else(.data$feature_id == feature, TRUE, .data$interference_corrected)
    )

  if(!is.na(updated_feature_id)) {
    data@dataset <- data@dataset |>
      mutate(feature_id = if_else(.data$feature_id == feature & .data$interference_corrected, updated_feature_id, .data$feature_id))
  }

  data@is_isotope_corr <- TRUE
  data@status_processing <- "Isotope-corrected raw data"
  data <- update_after_normalization(data, FALSE)
  data@var_drift_corrected <- c(feature_intensity = FALSE, feature_norm_intensity = FALSE, feature_conc = FALSE)
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]

  cli_alert_success(col_green(glue::glue("Interference-correction was manually applied to the feature `{variable}`.")))

  data
}


#' Apply interference correction
#'
#' @description This function corrects lipidomics feature intensities by
#' subtracting interference (e.g., isotope overlap or in-source fragments). The
#' correction is applied using the following formula: \deqn{value\_corrected =
#' value\_raw - value\_raw\_interfering\_feature \times
#' proportion\_interference}
#'
#' The interfering features and their relative contributions must be defined in
#' the feature metadata.
#'
#' By default, sequential series of interferences (e.g., isotopic M+2
#' interferences of PC 34:2 > PC 34:``1 > PC 34:0) will be corrected in a
#' sequential manner. This means that the correction will applied iteratively,
#' starting with the most downstream feature in the series. To disable this
#' behavior, basing each correction on the raw signal of the interfering feature
#' set `sequential_correction = FALSE`
#'
#' @details For isotopic interference correction of MRM/PRM data, the relative
#' isotope abundances needed for the calculation ('proportion_interference') can
#' be calculated using the LICAR application (Gao et al., 2021), see below..
#'
#' @param data MidarExperiment object containing lipidomics data.
#' @param variable Name of the variable to be corrected. Default:
#'   `feature_intensity`.
#' @param sequential_correction A logical indicating whether to apply
#'   corrections sequentially, starting with the most downstream feature.  If
#'   `FALSE`, corrections are based on the raw signal of the interfering
#'   features. If `FALSE`, the correction will be based on the raw signal of the
#'   interfering feature.
#' @return MidarExperiment object with feature intensities corrected for
#'   interferences.
#' @export
#' @references Gao L., Ji S, Burla B, Wenk MR, Torta F, Wenk MR, &
#' Cazenave-Gassiot A (2021). LICAR: An Application for Isotopic Correction of
#' Targeted Lipidomic Data Acquired with Class-Based Chromatographic Separations
#' Using Multiple Reaction Monitoring. *Analytical Chemistry*, 93(6), 3163-3171.
#' \url{https://doi.org/10.1021/acs.analchem.0c04565}


correct_interferences <- function(data = NULL, variable = "feature_intensity", sequential_correction = TRUE) {

  check_data(data)

  if (variable != "feature_intensity") cli::cli_abort("Currently only correction for raw intensities suspported, thus must be set to `feature_intensity` or not defined.")

  # Check if data is already interference-corrected
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

  # Join with feature metadata for interference information
  d_correct <- data@dataset |>
    left_join(
      data@annot_features |> select("feature_id", "interference_feature_id", "interference_contribution"),
      by = "feature_id")

  # Get a table with features to correct this table is ordered based on chained
  # relationships between feature_id and interference_feature_id, starting with
  # the most downstream feature in the chain. If the corrections are
  # independent, the order will be based on the order of the features in the
  # dataset. This code also checks for circular dependencies in the interference
  # correction, like LPC 18:2 > LPC 18:1 > LPC 18:0 > LPC 18:2
  features_to_correct <- d_correct |>
    filter(!is.na(.data$interference_feature_id)) |>
    select("feature_id", "interference_feature_id", "interference_contribution") |>
    distinct() |>
    order_chained_columns_tbl("feature_id", "interference_feature_id",
                              include_chain_id = FALSE, disconnected_action = "keep")

  # Check if there are incomplete interference data
  if(!all(complete.cases(features_to_correct))) {
    cli_abort("Some features have incomplete interference information (i.e., `interference_contribution` or `interference_contribution` missing. Please verify feature metadata.")
  }

  has_overlapping_interferences <- any(features_to_correct$interference_feature_id %in% features_to_correct$feature_id)

  # if sequential_correction is TRUE, reorder features_to_correct to ensure that
  # the most downstream feature is corrected first
  if (sequential_correction) {
    tryCatch({
      features_to_correct <- features_to_correct |>
        arrange(desc(row_number()))
    }, error = function(e) {
      if(e$message == "Circular dependency detected."){
        cli_abort(col_red("One or more circular correction(s) detected. Please verify the interfernece correction details defined in feature metadata."))
      }
    })
  }

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
                                   features_to_correct$interference_contribution[i], NA_real_),
        corr_intensity = if_else(.data$feature_id == features_to_correct$feature_id[i],
                                 raw_target - raw_source * rel_interference, NA_real_)
      ) |>
      ungroup() |>
      select(-"raw_target", -"raw_source", -"rel_interference") |>
      mutate(
        feature_intensity_orig = if_else(!is.na(.data$corr_intensity),
                                         .data$corr_intensity,
                                         .data$feature_intensity_orig)
        )
    d
  }

  # Apply corrections iteratively, using the result from the previous iteration
  d_corrected <- d_correct
  for (i in 1:nrow(features_to_correct)) {
    d_corrected <- correct_feature_intensity(d_corrected, features_to_correct, i)
  }

  # copy corrected intensities to original dataset

  data@dataset <- data@dataset |>
    left_join(d_corrected |>
                select("analysis_id", "feature_id",
                       intensity_corrected  = "feature_intensity_orig"),
              by = c("analysis_id", "feature_id")) |>
    mutate(
      interference_corrected = .data$feature_intensity != .data$intensity_corrected,
      feature_intensity = .data$intensity_corrected
    ) |>
    select(-"intensity_corrected", -"feature_intensity_orig")


  # Update MidarExperiment flags

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
