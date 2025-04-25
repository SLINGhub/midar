#' @title Calibrate Features Values Using Reference Sample
#'
#' @description
#' Calibrates feature abundances using a specified reference sample. Two approaches are supported:
#'
#' 1. Absolute calibration (when `absolute_calibration = TRUE`)
#'
#'    Recalibrates feature abundances based on known concentrations of the corresponding features
#'    defined for a reference sample.
#'
#'    The input variable can either `conc`, `norm_intensity`, or `intensity, whereas the result will
#'    always be stored under the variable`conc` (concentration), in the unit defined
#'    for the feature concentrations in the reference sample.
#'
#'
#'    Metadata requirements:
#'    - `sample_id` and `analyte_id` must be defined for the reference sample and features in the analysis and feature metadata, respectively.
#'    - Known analyte concentrations must be defined in the `QC concentration` metadata for the for the reference sample
#'    - An error will be raised  if no concentrations are defined for any features
#'
#'    Missing analyte concentrations for the reference sample can be handed via `undefined_conc_handling` with following options:
#'    - `original`: Keep original feature values, i.e. the non-calibrated values will be returned. *Note*: this is only available when `variable = conc`. Use with caution to avoid mixing units.
#'    - `na`: Set affected features values to `NA`
#'    - `error` (default): The function stops with error in case of any undefined reference sample feature concentration.
#'    - In case all feature concentrations are undefined, the function will stop with an error.
#'
#'   The re-calibrated feature concentrations are stored as `conc`, overwriting existing `conc` values.
#'   The original `conc` values are stored as `conc_beforecal`.
#'
#'   To export the calibrated concentrations use `save_dataset_csv()` with `variable = "conc",
#'   or to export non-calibrated values with `variable = "conc_beforecal"`.
#'   When saving the MiDAR XLSX report, the calibrated concentrations will also be stored as `conc`.
#'
#' 2. Normalization (relative calibration, `absolute_calibration = FALSE`)
#'
#'    Normalizes features abundances with corresponding feature abundances in a reference sample,
#'    resulting in ratios. Any available feature abundance variable
#'    (i.e., `conc`, `norm_intensity`, or `intensity`) can be used as input. The normalization is calculate for all present features.
#'    The resulting output will be stored as `[VARIABLE]_normalized`, whereby `[VARIABLE]` is the input variable, e.g., `conc_normalized`.

#'
#'    To export the normalized abundances , use `save_dataset_csv()` with `variable = "[VARIABLE]_normalized"`
#'    For MiDAR XLSX report, use `save_report_xlsx()` with same variable setting as for `save_dataset_csv()` to
#'    When saving the MiDAR XLSX report via `save_report_xlsx()`, availble unfiltered normalized feature abundances
#'    will be included by default. To include filtered normalized feature abundances, set `filtered_variable = "[VARIABLE]_normalized"`.
#'
#' Note: When multiple measurements (analyses) of the same reference sample are present,
#'    their concentrations are summarized by `mean` (default) or `median`. This can be set
#'    via the `summarize_fun` argument.
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
#' @examples
#'
#' dat_file = system.file("extdata", "S1P_MHQuant.csv", package = "midar")
#' meta_file = system.file("extdata", "S1P_metadata_tables.xlsx", package = "midar")

#' # Load data and metadata
#' mexp <- MidarExperiment()
#' mexp <- import_data_masshunter(mexp, dat_file, import_metadata = FALSE)
#' mexp <- import_metadata_analyses(mexp, path = meta_file, sheet = "Analyses")
#' mexp <- import_metadata_features(mexp, path = meta_file, sheet = "Features")
#' mexp <- import_metadata_istds(mexp, path = meta_file, sheet = "ISTDs")
#'
#' # Load known feature concentrations in the reference sample
#' mexp <- import_metadata_qcconcentrations(mexp, path = meta_file, sheet = "QCconcentrations")
#' mexp <- normalize_by_istd(mexp)
#' mexp <- quantify_by_istd(mexp)
#'
#' # Absolute calibration
#' # --------------------
#'
#'   mexp <- calibrate_by_reference(
#'     data = mexp,
#'     variable = "conc",
#'     reference_sample_id = "SRM1950",
#'     absolute_calibration = TRUE,
#'     summarize_fun = "mean",
#'     undefined_conc_action = "original"
#'   )
#'
#'   # Export absolute calibration concentrations
#'   save_dataset_csv(mexp, "calibrated.csv", variable = "conc", filter_data = FALSE)
#'
#'   # Export non-calibrated concentrations
#'   save_dataset_csv(mexp, "noncalibrated.csv", variable = "conc_beforecal", filter_data = FALSE)
#'
#'   # Create XLSX report with calibrated concentrations as filtered dataset
#'   save_report_xlsx(mexp, "report.xlsx", filtered_variable = "conc")
#'
#' # Relative calibration
#' # --------------------
#'
#'   mexp <- calibrate_by_reference(
#'     data = mexp,
#'     variable = "conc",
#'     reference_sample_id = "SRM1950",
#'     absolute_calibration = FALSE
#'   )
#'
#'   # Export SRM1950-normalized concentrations
#'   save_dataset_csv(mexp, "normalized.csv", variable = "conc_normalized", filter_data = FALSE)
#'
#'   # Create XLSX report with SRM1950-normalized concentrations as filtered dataset
#'   save_report_xlsx(mexp, "report.xlsx", filtered_variable = "conc_normalized")
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

 if(variable != "feature_conc" && !is.null(undefined_conc_action) && undefined_conc_action == "original")
    cli::cli_abort(col_red("When using `undefined_conc_action = 'original'`, the variable must be 'conc'. See the function's documentation for more details."))

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
        cli::cli_alert_warning(col_yellow("One or more feature concentration are not defined in the reference sample {reference_sample_id}. `NA` will be returned for these features. To change this, modify `undefined_conc_action` argument." ))
      } else if(undefined_conc_action == "original") {
        cli::cli_alert_warning(col_yellow("One or more feature concentration are not defined in the reference sample {reference_sample_id}. Original values will be returned for these. To change this, modify `undefined_conc_action` argument." ))
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

    if(variable == "feature_conc") {
      data@dataset <- data@dataset |>
        mutate(feature_conc_beforecal = .data$feature_conc)
    }

    data@dataset <- data@dataset |>
      dplyr::left_join(d_temp |> dplyr::select("analysis_id", "feature_id", "var_calibrated"), by = c("analysis_id", "feature_id")) |>
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
