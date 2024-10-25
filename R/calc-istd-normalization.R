#' get_conc_unit
#'
#' @param sample_amount_unit MidarExperiment object
#' @return string with feature_conc unit
#' @noRd

get_conc_unit <- function(sample_amount_unit) {
  units <- tolower(unique(sample_amount_unit))

  if (length(units) > 1) {
    conc_unit <- "pmol/sample amount unit (multiple units)"
  } else if (units == "ul" | units == "\U003BCl") {
    conc_unit <- "\U003BCmol/L"
  } else {
    conc_unit <- glue::glue("pmol/{units}")
  }
  conc_unit
}



#' Normalize Intensities with corresponding ISTD Intensities
#'
#' @param data MidarExperiment object
#' @param error_missing_info Raise error when ISTD is not defined for 1 or more features (excluding ISTD themselves) . Default: `TRUE`.
# #' @param ignore_unused_istds Ignore ISTDs with missing concentrations that are not used in any feature quantitation. Default: `FALSE`.
#' @return MidarExperiment object
#' @export

# TODO: interference correction not done here
calc_normalize_by_istd <- function(data = NULL, error_missing_info = TRUE) {

  check_data(data)

  if (nrow(data@annot_features) < 1) cli::cli_abort("No feature metadata available. Please add matching feature metadata.")
  if (any(!is.na(data@annot_features$interference_feature_id) & !data@is_isotope_corr))
    cli::cli_alert_warning(cli::col_yellow("Interfering feature intensities defined in metadata, but no interference correction was applied. Use `correct_interferences()` to correct."))


  # Check if data is already ISTD normalized
  if ("feature_norm_intensity" %in% names(data@dataset)) {
    if (!all(is.na(data@dataset$feature_norm_intensity))) cli::cli_alert_warning(cli::col_yellow("Overwriting previously normalized feature intensities."))
    data@dataset <- data@dataset |> select(-dplyr::any_of(c("feature_norm_intensity", "pmol_total", "feature_conc", "CONC_DRIFT_ADJ", "CONC_ADJ")))  #TODO fields
  }


  # Check if all ISTDs are defined as distinct feature in the feature metadata
  d_annot <- data@annot_features |> select("feature_id", "norm_istd_feature_id")
  all_istds <- unique(d_annot$norm_istd_feature_id)

  istd_not_defined <- setdiff(all_istds, d_annot$feature_id)
  if (length(istd_not_defined) > 0) {
    cli::cli_abort(cli::col_red("All ISTDs must be defined as feature in the feature metadata, {nrow(istd_not_defined)} ISTD(s) were not. Please check metadata."))
  }

  # check if ISTDs are defined for all features (except ISTDs that are not defined for themselves)
  features_no_istd <- data@annot_features |>
    filter(.data$valid_feature, !.data$is_istd, is.na(.data$norm_istd_feature_id)) |>
    dplyr::semi_join(data@dataset, by = c("feature_id"))

  if (nrow(features_no_istd) > 0) {
    if(!error_missing_info)
      cli::cli_alert_warning(cli::col_yellow("No ISTD was defined for {nrow(features_no_istd)} features, normalized intensities will be `NA` for these. "))
    else
      cli::cli_abort(cli::col_red("No ISTD was defined for {nrow(features_no_istd)} features. Please ammend feature metadata or set `error_missing_info = FALSE`."))
  }

  # Add ISTD intensities to temporary dataset
  d_temp <- data@dataset |>
    dplyr::left_join(d_annot, by = c("feature_id" = "feature_id"))

  # Normalize intensities
  d_temp <- d_temp |>
    dplyr::group_by(.data$norm_istd_feature_id, .data$analysis_id) |>
    dplyr::mutate(feature_norm_intensity = .data$feature_intensity / .data$feature_intensity[.data$is_istd]) |>
    dplyr::ungroup()

  # Add normalized intensities to dataset table
  data@dataset <- data@dataset |>
    dplyr::inner_join(d_temp |> dplyr::select("analysis_id", "feature_id", "feature_norm_intensity"), by = c("analysis_id", "feature_id"))

  # Print summary
  n_features <- length(unique(d_temp$feature_id))

  istds <- data@annot_features |> filter(.data$is_istd) |> pull(.data$feature_id)
  n_used_istds <- length(intersect(istds, data@annot_features |> filter(!.data$is_istd) |> pull(.data$quant_istd_feature_id) |> unique()))



  cli_alert_success(col_green(glue::glue("{n_features - length(istds)} features normalized with {n_used_istds} ISTDs in {get_analysis_count(data)} analyses.")))

  # Update status
  data@status_processing <- "ISTD-normalized ata"
  data@is_istd_normalized <- TRUE
  data@is_quantitated <- FALSE
  data@is_drift_corrected <- FALSE
  data@is_batch_corrected <- FALSE
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]
  data
}

#' Calculate analyte concentrations
#'
#' Calculation is based on ISTD-normalized intensities and corresponding
#' sample and spiked-in ISTD amounts. The determined concentration unit is based
#' on the `sample_amount_unit` provided in the analysis metadata.
#'  If the spiked-in ISTD concentrations are missing, the concentrations of the
#' @param data MidarExperiment object
#' @param error_missing_info Raise error when 1 or more ISTD concentratios and sample/ISTD amounts are missing. Default: `TRUE`.
#' @param ignore_unused_istds Ignore ISTDs with missing concentrations that are not used in any feature quantitation. Default: `FALSE`.
#' @return MidarExperiment object
#' @export
calc_quant_by_istd <- function(data = NULL, error_missing_info = TRUE, ignore_unused_istds = TRUE) {

  check_data(data)

  if (nrow(data@annot_istd) < 1) cli::cli_abort("ISTD concentrations are missing...please import ISTD metadata first.")
  if (!(c("feature_norm_intensity") %in% names(data@dataset))) cli::cli_abort("Data needs first to be ISTD normalized. Please run 'calc_normalize_by_istd' first.")

  istd_no_conc <- setdiff( data@annot_features$quant_istd_feature_id, data@annot_istd$quant_istd_feature_id)
  istds <- data@annot_features |> filter(.data$is_istd) |> pull(.data$feature_id)
  unused_istds <- setdiff(istds, data@annot_features |> filter(!.data$is_istd) |> pull(.data$quant_istd_feature_id) |> unique())

  # Check if ISTD concentrations in spiked-in mix  are defined for all ISTDs
  if (length(istd_no_conc) > 0 & !ignore_unused_istds) {
    if(!error_missing_info)
      cli::cli_alert_warning(cli::col_yellow("Spiked-in concentrations of {length(istd_no_conc)} ISTD(s) missing, calculated concentrations of affected features will be `NA`."))
    else
      cli::cli_abort(cli::col_red("Concentrations of {length(istd_no_conc)} ISTD(s) missing. Please ammend ISTD metadata or set `error_missing_info = FALSE`."))
  }

  # Check if sample and ISTD amounts are defined for all analyses
  samples_no_amounts <- data@annot_analyses |>
    filter(.data$valid_analysis, is.na(.data$sample_amount), is.na(.data$istd_volume)) |>
    dplyr::semi_join(data@dataset, by = c("analysis_id"))

  if (nrow(samples_no_amounts) > 0) {
    if(!error_missing_info)
      cli::cli_alert_warning(cli::col_yellow("Sample and/or ISTD solution amount(s) for {length(samples_no_amounts)} analyses missing, concentrations of all features of these analyses will be `NA`"))
    else
      cli::cli_abort(cli::col_red("Sample and/or ISTD amount(s) for {nrow(samples_no_amounts)} analyses missing. Please ammend analysis metadata or set `error_missing_info = FALSE`."))
  }

  # Add ISTD concentrations and sample amouts to temporary dataset
  d_temp <- data@dataset |>
    select(!any_of(c("sample_amount", "sample_amount_unit", "istd_volume", "pmol_total", "feature_conc", "CONC_DRIFT_ADJ", "CONC_ADJ"))) |>
    dplyr::left_join(data@annot_analyses |> dplyr::select("analysis_id", "sample_amount", "istd_volume"), by = c("analysis_id")) |>
    dplyr::left_join(data@annot_features |> dplyr::select("feature_id", "quant_istd_feature_id", "response_factor"), by = c("feature_id")) |>
    dplyr::left_join(data@annot_istd, by = c("quant_istd_feature_id"))

  # Calculate concentrations
  d_temp <- d_temp |> mutate(pmol_total = (.data$feature_norm_intensity) * (.data$istd_volume * (.data$istd_conc_nmolar)) * .data$response_factor / 1000)
  d_temp <- d_temp |> mutate(feature_conc = .data$pmol_total / .data$sample_amount)

  if ("feature_conc" %in% names(data@dataset)) {
    data@dataset <- data@dataset |> select(-dplyr::any_of(c("pmol_total", "feature_conc")))
    cli::cli_alert_warning(cli::col_yellow("Overwriting previously calculated concentrations."))
  }
  # Add calculated concentrations to dataset table
  data@dataset <- data@dataset |>
    dplyr::left_join(d_temp |> dplyr::select("analysis_id", "feature_id", "pmol_total", "feature_conc"), by = c("analysis_id", "feature_id"))

  # Copy of the raw concentrations (after normalization and quantitation)
  data@dataset$feature_conc_raw <- data@dataset$feature_conc

  n_features <- length(unique(d_temp$feature_id))
  n_istd_with_conc <- intersect(data@annot_features$quant_istd_feature_id, data@annot_istd$quant_istd_feature_id)

  n_features_with_conc <- d_temp |> filter(!.data$is_istd, .data$quant_istd_feature_id %in% n_istd_with_conc) |> select(.data$feature_id) |> distinct() |> nrow()
  n_istd <- length(unique(d_temp$quant_istd_feature_id)) - length(istd_no_conc)

  conc_unit <- get_conc_unit(data@annot_analyses$sample_amount_unit)

  cli_alert_success(cli::col_green("{n_features_with_conc} feature concentrations calculated based on {n_istd} ISTDs and sample amounts of {get_analysis_count(data)} analyses."))
  cli::cli_alert_info("Concentrations are given in {conc_unit}.")

  data@status_processing <- "ISTD-quantitated data"
  data@is_istd_normalized <- TRUE
  data@is_quantitated <- TRUE
  data@is_drift_corrected <- FALSE
  data@is_batch_corrected <- FALSE
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]

  data
}
