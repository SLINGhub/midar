#' Normalize Feature Intensities Using Internal Standards
#'
#' Normalize feature intensities by dividing them by the intensities
#' of the corresponding internal standards (ISTDs). Each feature
#' must have a defined internal standard (ISTD) in the `feature` metadata.
#'

#' @param data A `MidarExperiment` object.
#' @param fail_missing_annotation If `TRUE`, the function
#'   will raise an error when an ISTD is not defined for one or more features (excluding the ISTDs themselves).
#'   If `FALSE`, features with missing ISTD annotations will have NA values in the normalized intensities.
#'
#' @return A `MidarExperiment` object with normalized feature intensities
#'
#' @export

normalize_by_istd <- function(data = NULL, fail_missing_annotation = TRUE) {

  check_data(data)

  if (nrow(data@annot_features) < 1) cli::cli_abort("No feature metadata available. Please add matching feature metadata.")
  if (any(!is.na(data@annot_features$interference_feature_id) & !data@is_isotope_corr))
    cli::cli_alert_warning(cli::col_yellow("Interfering features defined in metadata, but no correction was applied. Use `correct_interferences()` to correct."))


  # Check if data is already ISTD normalized
  if ("feature_norm_intensity" %in% names(data@dataset)) {
    if (!all(is.na(data@dataset$feature_norm_intensity))) cli::cli_alert_warning(cli::col_yellow("Replacing previously normalized feature intensities."))
    data@dataset <- data@dataset |> select(-dplyr::any_of(c("feature_norm_intensity", "feature_pmol_total", "feature_conc", "CONC_DRIFT_ADJ", "CONC_ADJ")))  #TODO fields
  }


  # Check if all ISTDs are defined as distinct feature in the feature metadata
  d_annot <- data@annot_features |> select("feature_id", "istd_feature_id")
  all_istds <- unique(na.omit(d_annot$istd_feature_id))

  if(!all(is.na(all_istds))){
    istd_not_defined <- setdiff(all_istds, d_annot$feature_id)
    if (length(istd_not_defined) > 0) {
      cli::cli_abort(cli::col_red("{length(istd_not_defined)} ISTD(s) were not defined as individual feature(s). Please check feature metadata."))
    }
  } else {
    cli::cli_abort(cli::col_red("No ISTDs defined in feature metadata. Please define ISTDs for each feature in feature metadata."))
  }

  # check if ISTDs are defined for all features (except ISTDs that are not defined for themselves)
  features_no_istd <- data@annot_features |>
    filter(.data$valid_feature, !.data$is_istd, is.na(.data$istd_feature_id)) |>
    dplyr::semi_join(data@dataset, by = c("feature_id"))

  if (nrow(features_no_istd) > 0) {
    if(!fail_missing_annotation)
      cli::cli_alert_warning(cli::col_yellow("For {nrow(features_no_istd)} feature(s) no ISTD was defined, normalized intensities will be `NA` for these features. "))
    else
      cli::cli_abort(cli::col_red("For {nrow(features_no_istd)} feature(s) no ISTD was defined. Please ammend feature metadata or set `fail_missing_annotation = FALSE`."))
  }

  # Add ISTD intensities to temporary dataset
  d_temp <- data@dataset |>
    dplyr::left_join(d_annot, by = c("feature_id" = "feature_id"))

  # Normalize intensities
  d_temp <- d_temp |>
    dplyr::group_by(.data$istd_feature_id, .data$analysis_id) |>
    dplyr::mutate(
      feature_norm_intensity = ifelse(
        !is.na(.data$istd_feature_id),
        .data$feature_intensity / .data$feature_intensity[.data$is_istd],
        NA_real_
      )
    ) |>
    dplyr::ungroup()

  # Add normalized intensities to dataset table
  data@dataset <- data@dataset |>
    dplyr::inner_join(d_temp |> dplyr::select("analysis_id", "feature_id", "feature_norm_intensity"), by = c("analysis_id", "feature_id"))

  # Print summary
  n_features <- length(unique(d_temp$feature_id)) - nrow(features_no_istd)

  istds <- data@annot_features |> filter(.data$is_istd) |> pull(.data$feature_id)
  n_used_istds <- length(intersect(istds, data@annot_features |> filter(!.data$is_istd) |> pull(.data$quant_istd_feature_id) |> unique()))



  cli_alert_success(col_green(glue::glue("{n_features - length(istds)} features normalized with {n_used_istds} ISTDs in {get_analysis_count(data)} analyses.")))

  # Update status
  data@status_processing <- "ISTD-normalized ata"
  data <- update_after_normalization(data, TRUE)
  data <- update_after_quantitation(data, FALSE)
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]
  data
}


#' Calculate Analyte Concentrations Based on Internal Standards
#'
#' Analyte concentrations are calculated based on ISTD-normalized intensities
#' and the corresponding spiked-in ISTD amount, normalized by the sample amount
#' (single-point calibration approach). The concentration unit is determined
#' from the `sample_amount_unit` field in the analysis metadata of the `MidarExperiment` object.
#'
#' The calculated concentrations are added to the `dataset` table as a new column
#' named `feature_conc`.
#'
#' @param data A `MidarExperiment` object
#' @param fail_missing_annotation If `TRUE`, an error will be raised if any of the following information is missing: ISTD concentration, ISTD mix volume, and sample amounts for any feature.
#'   If `FALSE`, missing annotations will be ignored, and resulting feature concentration will be `NA`
#'
#' @return A `MidarExperiment` object with the calculated analyte concentrations added to the
#'   `dataset` table in the `feature_conc` column.
#'
#' @seealso [quantify_by_calibration()] for calculating concentrations based on external calibration curves.
#'
#' @export

quantify_by_istd <- function(data = NULL, fail_missing_annotation = TRUE) {

  check_data(data)

  if (nrow(data@annot_istds) < 1) cli::cli_abort("ISTD concentrations are missing...please import ISTD metadata first.")
  if (!(c("feature_norm_intensity") %in% names(data@dataset))) cli::cli_abort("Data needs to be ISTD normalized, please run 'normalize_by_istd' first.")


  # Check if sample and ISTD amounts are defined for all analyses
  samples_no_amounts <- data@annot_analyses |>
    filter(.data$valid_analysis, is.na(.data$sample_amount) | is.na(.data$istd_volume)) |>
    dplyr::semi_join(data@dataset, by = c("analysis_id"))

  if (nrow(samples_no_amounts) > 0) {
    if(!fail_missing_annotation)
      cli::cli_alert_warning(cli::col_yellow("Sample and/or ISTD solution amount(s) for {nrow(samples_no_amounts)} analyses missing, concentrations of all features for these analyses will be `NA`"))
    else
      cli::cli_abort(cli::col_red("Sample and/or ISTD amount(s) for {nrow(samples_no_amounts)} analyses missing. Please ammend analysis metadata or set `fail_missing_annotation = FALSE`."))
  }

  # Add ISTD concentrations and sample amouts to temporary dataset
  d_temp <- data@dataset |>
    select(!any_of(c("sample_amount", "sample_amount_unit", "istd_volume", "feature_pmol_total", "feature_conc", "CONC_DRIFT_ADJ", "CONC_ADJ"))) |>
    dplyr::left_join(data@annot_analyses |> dplyr::select("analysis_id", "sample_amount", "istd_volume"), by = c("analysis_id")) |>
    dplyr::left_join(data@annot_features |> dplyr::select("feature_id", "quant_istd_feature_id", "response_factor"), by = c("feature_id")) |>
    dplyr::left_join(data@annot_istds, by = c("quant_istd_feature_id"))

  istd_no_conc <- setdiff( unique(d_temp[!d_temp$is_istd,]$quant_istd_feature_id), unique(data@annot_istds$quant_istd_feature_id))
  istds <- data@annot_features |> filter(.data$is_istd) |> pull(.data$feature_id)

  # Check if ISTD concentrations in spiked-in mix  are defined for all ISTDs
  if (length(istd_no_conc) > 0){
    if(!fail_missing_annotation) {
      cli::cli_alert_warning(cli::col_yellow("Spiked-in concentrations of {length(istd_no_conc)} ISTD(s) missing, calculated concentrations of affected features will be `NA`."))
    } else {
      cli::cli_abort(cli::col_red("Concentrations of {length(istd_no_conc)} ISTD(s) missing. Please ammend ISTD metadata or set `fail_missing_annotation = FALSE`."))
    }
  }


  # Calculate concentrations
  d_temp <- d_temp |> mutate(feature_pmol_total = (.data$feature_norm_intensity) * (.data$istd_volume * (.data$istd_conc_nmolar)) * .data$response_factor / 1000)
  d_temp <- d_temp |> mutate(feature_conc = .data$feature_pmol_total / .data$sample_amount)

  if ("feature_conc" %in% names(data@dataset)) {
    data@dataset <- data@dataset |> select(-dplyr::any_of(c("feature_pmol_total", "feature_conc")))
    cli::cli_alert_warning(cli::col_yellow("Replacing previously calculated concentrations."))
  }
  # Add calculated concentrations to dataset table
  data@dataset <- data@dataset |>
    dplyr::left_join(d_temp |> dplyr::select("analysis_id", "feature_id", "feature_pmol_total", "feature_conc"), by = c("analysis_id", "feature_id"))

  n_features <- length(unique(d_temp$feature_id))
  n_istd_with_conc <- intersect(data@annot_features$quant_istd_feature_id, data@annot_istds$quant_istd_feature_id)

  n_features_with_conc <- d_temp |> filter(!.data$is_istd, .data$quant_istd_feature_id %in% n_istd_with_conc) |> select("feature_id") |> distinct() |> nrow()
  n_istd <- length(unique(d_temp$quant_istd_feature_id)) - length(istd_no_conc)

  conc_unit <- get_conc_unit(data@annot_analyses$sample_amount_unit)

  cli_alert_success(cli::col_green("{n_features_with_conc} feature concentrations calculated based on {n_istd} ISTDs and sample amounts of {get_analysis_count(data) - nrow(samples_no_amounts)} analyses."))
  cli::cli_alert_info("Concentrations are given in {conc_unit}.")

  data@status_processing <- "ISTD-quantitated data"

  data <- update_after_normalization(data, TRUE)
  data <- update_after_quantitation(data, TRUE)
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]

  data
}
