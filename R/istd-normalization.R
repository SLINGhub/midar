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
#' @param interference_correction Apply interference (e.g. isotope) correction as defined in the feature metadata
#' @return MidarExperiment object

#' @importFrom glue glue
#' @export

normalize_by_istd <- function(data, interference_correction = TRUE) {
  if (nrow(data@annot_features) < 1) stop("ISTD map is missing...please import transition annotations.")

  if (any(!is.na(data@annot_features$interference_feature_id) & interference_correction)) data <- correct_interferences(data)

  if ("feature_norm_intensity" %in% names(data@dataset)) {
    if (!all(is.na(data@dataset$feature_norm_intensity))) warning("Overwriting exiting normalized Intensities")
    data@dataset <- data@dataset %>% select(-dplyr::any_of(c("feature_norm_intensity", "pmol_total", "feature_conc", "CONC_DRIFT_ADJ", "CONC_ADJ")))
  }

  d_temp <- data@dataset # %>%     dplyr::full_join(data@annot_features, by = c("feature_id" = "feature_id"),)
  d_temp <- d_temp %>%
    dplyr::group_by(.data$norm_istd_feature_id, .data$analysis_id) %>%
    dplyr::mutate(feature_norm_intensity = .data$feature_intensity / .data$feature_intensity[.data$is_istd]) %>%
    dplyr::ungroup()

  data@dataset <- data@dataset %>%
    dplyr::inner_join(d_temp %>% dplyr::select("analysis_id", "feature_id", "is_istd", "feature_norm_intensity"), by = c("analysis_id", "feature_id", "is_istd"))

  n_features <- length(data@dataset$feature_id |> unique())
  n_istd <- length(unique(data@dataset$norm_istd_feature_id))
  cli_alert_success(col_green(glue::glue("{n_features} features normalized with {n_istd} ISTDs. \n")))
  data@status_processing <- "ISTD-normalized Data"
  data@is_istd_normalized <- TRUE
  data@is_quantitated <- FALSE
  data@is_drift_corrected <- FALSE
  data@is_batch_corrected <- FALSE
  data
}

#' Quantitate using sample and spiked ISTD amounts
#'
#' @param data MidarExperiment object
#' @return MidarExperiment object
#' @importFrom glue glue
#' @export
quantitate_by_istd <- function(data) {
  if (nrow(data@annot_istd) < 1) stop("ISTD concetrations are missing...please import them first.")
  if (!(c("feature_norm_intensity") %in% names(data@dataset))) stop("Data needs first to be ISTD normalized. Please apply function 'normalize_by_istd' first.")
  d_temp <- data@dataset %>%
    select(!any_of(c("sample_amount", "sample_amount_unit", "istd_volume", "pmol_total", "feature_conc", "CONC_DRIFT_ADJ", "CONC_ADJ"))) |>
    dplyr::left_join(data@annot_analyses %>% dplyr::select("analysis_id", "sample_amount", "istd_volume"), by = c("analysis_id")) %>%
    # dplyr::left_join(data@annot_features %>% dplyr::select("feature_id", "quant_istd_feature_id"), by = c("feature_id")) %>%
    dplyr::left_join(data@annot_istd, by = c("quant_istd_feature_id"))

  d_temp <- d_temp %>% mutate(pmol_total = (.data$feature_norm_intensity) * (.data$istd_volume * (.data$istd_conc_nmolar)) * .data$feature_response_factor / 1000)
  d_temp <- d_temp %>% mutate(feature_conc = .data$pmol_total / .data$sample_amount)

  if ("feature_conc" %in% names(data@dataset)) {
    data@dataset <- data@dataset %>% select(-dplyr::any_of(c("pmol_total", "feature_conc")))
    warning("Overwriting exiting concs")
  }

  data@dataset <- data@dataset %>%
    dplyr::inner_join(d_temp %>% dplyr::select("analysis_id", "feature_id", "pmol_total", "feature_conc"), by = c("analysis_id", "feature_id"))

  data@dataset$conc_raw <- data@dataset$feature_conc

  n_features <- length(unique(data@dataset$feature_id))
  n_istd <- length(unique(data@dataset$norm_istd_feature_id))

  conc_unit <- get_conc_unit(data@annot_analyses$sample_amount_unit)

  cli_alert_success(col_green(glue::glue("{n_features} features quantitated in {nrow(data@annot_analyses)} samples using {n_istd} spiked-in ISTDs and sample amounts.
                   feature_conc unit: [{conc_unit}].")))

  data@status_processing <- "Quantitated Data"
  data@is_istd_normalized <- TRUE
  data@is_quantitated <- TRUE
  data@is_drift_corrected <- FALSE
  data@is_batch_corrected <- FALSE

  data
}
