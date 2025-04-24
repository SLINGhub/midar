#' Calibrate Features Values Using a Reference Sample
#'
#' @description
#' Calibrates features abundances using a reference sample, either by re-calibrating
#' absolute concentrations using known concentrations in the reference samples or
#' by calculating relative ratios to reference samples.
#'
#' Two calibration approaches are supported:
#'
#' 1. Absolute calibration (`absolute_calibration = TRUE`):
#'    Recalibrates features to concentrations using known reference sample values.
#'
#'    Requirements:
#'    - `analyte_id` and `sample_id` must be defined in analysis metadata
#'    - Known analyte concentrations must be defined for the reference sample in QC concentration metadata
#'
#'    The function:
#'    - Overwrites original values (`intensity`, `norm_intensity`, or `conc`) with recalibrated values
#'    - Stores original values in `feature_[VARIABLE]_beforecal`, where `[VARIABLE]` is `conc`, `norm_intensity`, or `intensity`
#'    - Can process `conc`, `norm_intensity`, or `intensity` data
#'    - Always outputs concentration (`conc`) in concentration unit defined for the analytes in the reference sample in the QC concentration metadata.
#'
#'    When multiple analyses of the same reference sample are present, their concentrations are summarized using
#'    `summarize_fun` ("mean" or "median", default: "mean").
#'
#'    Options for undefined analyte concentrations for the reference sample (`undefined_conc_handling`):
#'    - "original": Keep original values (meaning the non-calibrated values will be returned)
#'    - "na": Set undefined features to NA
#'    - "error": Raise error for undefined features
#'   An error will be raised  if no concentrations are defined for any features.
#'
#'   To export:
#'    - Use `save_dataset_csv()` with `variable = "[VARIABLE]`
#'    - For MiDAR XLSX report, use `save_report_xlsx()` also with  `variable = "[VARIABLE]` to
#'    include the re-calibrated values of the filtered dataset in the report.
#'    - Set `[VARIABLE]` = `conc`, `norm_intensity`, or `intensity` to export re-calibrated values, or
#'     `conc_beforecal`, `norm_intensity_beforecal`, or `intensity_beforecal` to export the original values.
#'
#' 2. Relative calibration (`absolute_calibration = FALSE`):
#'    Expresses sample features values as ratios relative to corresponding features in the reference sample.
#'
#'    - Applies uniformly to all features
#'    - Works with `conc`, `norm_intensity`, or `intensity` data
#'    - Stores results in `feature_[VARIABLE]_normalized`, where `[VARIABLE]` is `conc`, `norm_intensity`, or `intensity`
#'
#'    To export:
#'    - Use `save_dataset_csv()` with `variable = "feature_[VARIABLE]_normalized"`
#'    - For MiDAR XLSX report, use `save_report_xlsx()` with same variable setting as for `save_dataset_csv()` to
#'    include the normalized values of the filtered dataset in the report.
#'    Unfiltered normalized feature values will be included by default in the report if present
#'
#'
#' @param data A `MidarExperiment` object containing the metabolomics data to be normalized
#' @param variable Character string indicating which data type to calibrate Must be
#'   one of: "intensity", "norm_intensity", or "conc"
#' @param reference_sample_id Character vector specifying the sample ID(s) to use as
#'   reference(s) or standards
#' @param absolute_calibration Logical indicating whether to perform absolute calibration using
#'   known concentrations of the reference sample (TRUE) or relative calibration (FALSE).
#' @param summarize_fun Either "mean" or "median". If `absolute_calibration = TRUE`,
#' this function is used to summarize the reference sample concentrations across analyses of specified `reference_sample_id`. Default is "mean".
#' @param undefined_conc_action Character string specifying how to handle features
#'   without defined concentrations in reference samples when `absolute_calibration = TRUE`.
#'   Must be one of: "original" (keep original values), "na" (set to NA), or "error". Default is "keep".
#'
#' @return A `MidarExperiment` object with calibrated data
#'
#' @examples
#' \dontrun{
#' # Relative calibration
#' mexp <- calibrate_by_reference(
#'   data = mexp,
#'   variable = "conc",
#'   reference_sample_id = "NIST-SRM1950",
#'   absolute_calibration = TRUE
#'   summarize_fun = "mean",
#'   undefined_conc_action = "original"
#' )
#'
#' # Export relative calibration results
#' save_dataset_csv(mexp, "mexp_calibrated.csv", variable = "conc")
#'
##' # Export non-calibrated results
#' save_dataset_csv(mexp, "mexp_calibrated.csv", variable = "conc_beforecal")
#'
#' # Export relative calibration results in the MiDAR XLSX report
#' save_report_xlsx(mexp, "report.xlsx", variable = "conc")
#'
#' # Absolute calibration
#' mexp <- calibrate_by_reference(
#'   data = mexp,
#'   variable = "conc",
#'   reference_sample_id = "NIST-SRM1950",
#'   absolute_calibration = FALSE
#' )
#'
#' # Export relative calibration results
#' save_dataset_csv(mexp, "mexp_calibrated.csv", variable = "conc_normalized")
#'
#' # Export relative calibration results in the MiDAR XLSX report
#' save_report_xlsx(mexp, "report.xlsx", variable = "conc_normalized")
#' }
#'
#' @seealso
#' [normalize_by_istd()], [quantify_by_istd()], [quantify_by_calibration()]
#'
#' @export

calibrate_by_reference <- function(data, variable, reference_sample_id, absolute_calibration, summarize_fun = "mean", undefined_conc_action = NULL) {

  check_data(data)

  if(absolute_calibration && is.null(undefined_conc_action))
    cli::cli_abort(col_red("When using `absolute_calibration = TRUE`, then `undefined_conc_action` must be specified, as either 'original', 'na', or 'error'"))

  # Check if the reference sample is present in the dataset
  if (!any(data@dataset$sample_id %in% reference_sample_id)) {
    cli::cli_abort(col_red("The specified `reference_sample_id` is not present in the dataset. Please verify `reference_sample_id` and analysis metadata (column `sample_id`)."))
  }

  rlang::arg_match(summarize_fun, c("mean", "median"))

  variable_strip <- str_remove(variable, "feature_")
  rlang::arg_match(variable_strip, c("intensity", "norm_intensity", "conc"))
  variable <- stringr::str_c("feature_", variable_strip)
  variable_sym <- rlang::sym(variable)
  variable_norm_sym <- rlang::sym(stringr::str_c("feature_", variable_strip, "_normalized"))
  variable_beforecal_sym <- rlang::sym(stringr::str_c("feature_", variable_strip, "_beforecal"))
  check_var_in_dataset(data@dataset, variable)



  d_temp <- data@dataset |>
    select(any_of(c("analysis_id", "sample_id", "feature_id", "analyte_id", variable)))

  d_temp <- d_temp |>
    dplyr::mutate(ref_sample = reference_sample_id) |>
    dplyr::group_by(.data$feature_id) |>
    dplyr::mutate(
      var_calibrated = !!variable_sym / match.fun(summarize_fun)(pull(filter(pick(everything()), .data$sample_id == .data$ref_sample), !!variable_sym), na.rm = TRUE)) |>
    ungroup()


  if(absolute_calibration) {
    # Check if the reference sample has concentration values
    if (absolute_calibration && !any(reference_sample_id %in% unique(data@annot_qcconcentrations$sample_id)))
      cli::cli_abort(col_red("No concentration values found for the reference sample `{reference_sample_id}`. Please verify QC concentration metadata."))

    # Check if the reference sample has the same concentration unit for all features

    ref_feature_conc_unit <- unique(data@annot_qcconcentrations[data@annot_qcconcentrations$sample_id == reference_sample_id,]$concentration_unit)

    if (absolute_calibration && all(is.na(ref_feature_conc_unit)))
      cli::cli_abort(col_red("No concentration units found for features of reference sample '{reference_sample_id}'. Please check concentration units in QC concentration metadata."))
    else if (absolute_calibration && length(ref_feature_conc_unit) > 1)
      cli::cli_abort(col_red("Different unit used for feature concentrations in reference sample '{reference_sample_id}'. Please check concentration units in QC concentration metadata."))


    if(length(setdiff(unique(data@dataset[data@dataset$sample_id == reference_sample_id,]$analyte_id),
                                      unique(data@annot_qcconcentrations[data@annot_qcconcentrations$sample_id == reference_sample_id,]$analyte_id)) > 0)){

      if(undefined_conc_action == "error")
        cli::cli_abort(col_red("One or more feature concentration are not defined in the reference sample {reference_sample_id}. Please verify QC concentration metadata or modify `undefined_conc_action` argument" ))
      else if(undefined_conc_action == "na") {
        cli::cli_alert_warning(col_yellow("One or more feature concentration are not defined in the reference sample {reference_sample_id}. `NA` will be returned for these features. Modify `undefined_conc_action` argument to change behaviour." ))
      } else if(undefined_conc_action == "original") {
        cli::cli_alert_warning(col_yellow("One or more feature concentration are not defined in the reference sample {reference_sample_id}. Original values will be returned for these. Modify `undefined_conc_action` argument to change behaviour." ))
      } else
        cli::cli_abort(col_red("Invalid value for `undefined_conc_action`. Must be one of: 'original', 'na', or 'error'"))
    }

    # Get the reference sample concentrations

    if(data@is_filtered){
      cli::cli_alert_warning(col_yellow("Previously filtered dataset is no longer valid and has been cleared. Please re-apply feature filtering."))
      data@is_filtered <- FALSE
      data@metrics_qc <- data@metrics_qc[FALSE,]
    }

    d_temp <- d_temp |>
      dplyr::left_join(data@annot_qcconcentrations |>
                         filter (.data$sample_id == reference_sample_id) |>
                         dplyr::select( "analyte_id", ref_conc= "concentration", ref_conc_unit = "concentration_unit"),
                       by = c("analyte_id")) |>
      dplyr::group_by(.data$feature_id) |>
      mutate(
        var_calibrated = if_else(!is.na(.data$ref_conc),
                                 .data$var_calibrated * .data$ref_conc,
                                 if(undefined_conc_action == "original") !!variable_sym else NA_real_)) |>
      dplyr::select(-"ref_conc")

    data@dataset <- data@dataset |>
      dplyr::left_join(d_temp |> dplyr::select("analysis_id", "feature_id", "var_calibrated"), by = c("analysis_id", "feature_id")) |>
      mutate(feature_conc_beforecal = .data$feature_conc) |>
      dplyr::mutate(feature_conc = .data$var_calibrated) |>
      dplyr::select(-"var_calibrated")

    n_features_with_conc <- data@annot_qcconcentrations |>
      filter(.data$sample_id == reference_sample_id) |>
      nrow()

    if(variable_strip == "conc")
      cli_alert_success(cli::col_green("{n_features_with_conc} feature concentrations re-calibrated using the reference sample {reference_sample_id}."))
    else
      cli_alert_success(cli::col_green("{n_features_with_conc} feature concentrations calculated using the defined reference sample concentrations."))


    cli::cli_alert_info(col_green("Concentrations are given in {ref_feature_conc_unit}."))

    data <- update_after_normalization(data, is_normalized = TRUE)
    data <- update_after_quantitation(data, is_quantitated = TRUE)


    data@is_filtered <- FALSE
    data@metrics_qc <- data@metrics_qc[FALSE,]

    data@status_processing <- paste0("Re-calibrated", " ", data@status_processing )


  } else {
    # If absolute_calibration is FALSE, create a new variable variable_noncalib to backup values before calibration and then join data@dataset with d_temp and rename var_calibrated to variable.
    data@dataset <- data@dataset |>
      dplyr::left_join(d_temp |> dplyr::select("feature_id", "analysis_id",  "var_calibrated"), by = c("analysis_id", "feature_id")) |>
      mutate(!!variable_norm_sym := .data$var_calibrated) |>
      dplyr::select(-"var_calibrated")

    cli_alert_success(cli::col_green("All features normalized with reference sample {reference_sample_id} features."))
    cli::cli_alert_info(col_green("Unit is: sample [{variable_strip}] / {reference_sample_id} [{variable_strip}]"))

    data@status_processing <- paste0("Reference sample normalized", " ", data@status_processing )
  }
  data
}
