#' get_conc_unit
#'
#' @param sample_amount_unit MidarExperiment object
#' @return string with feature_conc unit
#' @noRd

get_conc_unit <- function(sample_amount_unit){

  units <- tolower(unique(sample_amount_unit))

  if (length(units) > 1)
    conc_unit <- "pmol/sample amount unit (multiple units)"
  else if (units == "ul" | units == "\U003BCl")
    conc_unit <- "\U003BCmol/L"
  else
    conc_unit <- glue::glue("pmol/{units}")
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
  if(nrow(data@annot_features) < 1) stop("ISTD map is missing...please import transition annotations.")

  if (any(!is.na(data@annot_features$interference_feature_name) & interference_correction)) data <- correct_interferences(data)

  if("feature_norm_intensity" %in% names(data@dataset)) {
    if (!all(is.na(data@dataset$feature_norm_intensity))) warning("Overwriting exiting normalized Intensities")
    data@dataset <- data@dataset %>% select(-dplyr::any_of(c("feature_norm_intensity", "pmol_total", "feature_conc", "CONC_DRIFT_ADJ", "CONC_ADJ")))
  }

  d_temp <- data@dataset #%>%     dplyr::full_join(data@annot_features, by = c("feature_name" = "feature_name"),)
  d_temp <- d_temp  %>%
    dplyr::group_by(.data$norm_istd_feature_name, .data$analysis_id) %>%
    dplyr::mutate(feature_norm_intensity = .data$feature_intensity/.data$feature_intensity[.data$is_istd]) %>%
    dplyr::ungroup()

  data@dataset <- data@dataset %>%
    dplyr::inner_join(d_temp %>% dplyr::select("analysis_id", "feature_name", "is_istd", "feature_norm_intensity"), by = c("analysis_id", "feature_name", "is_istd"))

  n_features <- length(data@dataset$feature_name |> unique())
  n_istd <- length(unique(data@dataset$norm_istd_feature_name))
  writeLines(crayon::green(glue::glue("\u2713 {n_features} features normalized with {n_istd} ISTDs. \n")))
  data@status_processing <- "ISTD-normalized Data"
  data@is_istd_normalized <- TRUE
  data@is_quantitated <- FALSE
  data@is_drift_corrected = FALSE
  data@is_batch_corrected = FALSE
  data
}

#' Quantitate using sample and spiked ISTD amounts
#'
#' @param data MidarExperiment object
#' @return MidarExperiment object
#' @importFrom glue glue
#' @export
quantitate_by_istd <- function(data) {
  if(nrow(data@annot_istd) < 1) stop("ISTD concetrations are missing...please import them first.")
  if(!(c("feature_norm_intensity") %in% names(data@dataset))) stop("Data needs first to be ISTD normalized. Please apply function 'normalize_by_istd' first.")
  d_temp <- data@dataset %>%
    dplyr::left_join(data@annot_analyses %>% dplyr::select("analysis_id", "sample_amount", "istd_volume"), by = c("analysis_id")) %>%
    #dplyr::left_join(data@annot_features %>% dplyr::select("feature_name", "quant_istd_feature_name"), by = c("feature_name")) %>%
    dplyr::left_join(data@annot_istd, by = c("quant_istd_feature_name"))

  d_temp <- d_temp %>% mutate(pmol_total = (.data$feature_norm_intensity)*(.data$istd_volume*(.data$istd_conc_nmolar)) * .data$feature_response_factor /1000)
  d_temp <- d_temp %>% mutate(feature_conc = .data$pmol_total/.data$sample_amount)

  if("feature_conc" %in% names(data@dataset)) {
    data@dataset <- data@dataset %>% select(-dplyr::any_of(c("pmol_total", "feature_conc")))
    warning("Overwriting exiting concs")
  }

  data@dataset <- data@dataset %>%
    dplyr::inner_join(d_temp %>% dplyr::select("analysis_id", "feature_name", "sample_amount",  "istd_volume", "pmol_total", "feature_conc"), by = c("analysis_id", "feature_name"))

  data@dataset$conc_raw <- data@dataset$feature_conc

  n_features <- length(unique(data@dataset$feature_name))
  n_istd <- length(unique(data@dataset$norm_istd_feature_name))

  conc_unit <- get_conc_unit(data@annot_analyses$sample_amount_unit)

  writeLines(crayon::green(glue::glue("\u2713 {n_features} features quantitated in {nrow(data@annot_analyses)} samples using {n_istd} spiked-in ISTDs and sample amounts.
                   feature_conc unit: [{conc_unit}].")))

  data@status_processing <- "Quantitated Data"
  data@is_istd_normalized <- TRUE
  data@is_quantitated <- TRUE
  data@is_drift_corrected = FALSE
  data@is_batch_corrected = FALSE

  data
}




#' Export any parameter to a wide-format table
#'
#' @param data MidarExperiment object
#' @param variable Variable to be exported
#' @param filename File name with path of exported CSV file
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider
#' @export
exportWideCSV <- function(data, variable, filename) {

  var <- dplyr::sym(variable)

  if (!(variable %in% names(data@dataset))) stop("Variable '", variable,  "' does not (yet) exist in dataset.")

  ds <- data@dataset |>
    dplyr::select("analysis_id", "qc_type", "acquisition_time_stamp", "feature_name", !!var) %>%
    tidyr::pivot_wider(names_from = .data$feature_name, values_from = !!var)

  readr::write_csv(ds, file = filename, num_threads = 4, col_names = TRUE)
  invisible(ds)

}

#' Calculate QC metrics for each feature
#'
#' @param data MidarExperiment object
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider nest unnest

#' @importFrom purrr map
#' @importFrom broom glance
#' @importFrom dplyr summarise
#' @return MidarExperiment object
#' @export
calculate_qc_metrics <- function(data) {

  #if(!(c("feature_norm_intensity") %in% names(data@dataset))) warning("No normali is not normalized")


  #TODO: remove later when fixed
  if (tolower(data@analysis_type) == "lipidomics") data <- lipidomics_get_lipid_class_names(data)

  ds1 <- data@dataset %>%
      dplyr::filter(.data$qc_type %in% c("SPL", "NIST", "LTR", "BQC", "TQC", "PBLK", "SBLK", "UBLK")) %>%
      dplyr::group_by(.data$feature_name, .data$feature_class) %>%
      dplyr::summarise(
        #PrecursorMz = paste0(unique(.data$precursor_mz), collapse = ","),
        #ProductMz = paste0(unique(.data$product_mz), collapse = ","),
        valid_integration = unique(.data$valid_integration),
        is_quantifier = unique(.data$is_quantifier),
        is_istd = unique(.data$is_istd),
        norm_istd = unique(.data$norm_istd_feature_name),
        quant_istd = unique(.data$quant_istd_feature_name),
        feature_response_factor = unique(.data$feature_response_factor),
        Int_med_PBLK = median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE),
        Int_med_SPL = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
        Int_med_BQC = median(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
        Int_med_TQC = median(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
        Int_med_NIST = median(.data$feature_intensity[.data$qc_type == "NIST"], na.rm = TRUE),
        Int_med_LTR = median(.data$feature_intensity[.data$qc_type == "LTR"], na.rm = TRUE),

        conc_median_TQC = median(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE),
        conc_median_BQC = median(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE),
        conc_median_SPL = median(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE),
        conc_median_NIST = median(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE),
        conc_median_LTR = median(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE),

        SB_Ratio_Q10 = quantile(.data$feature_intensity[.data$qc_type == "SPL"], probs  = 0.1, na.rm = TRUE, names = FALSE)/median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE, names = FALSE),
        SB_Ratio_median = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE, names = FALSE)/median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE, names = FALSE),

        Int_CV_TQC = sd(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE)/mean(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE) * 100,
        Int_CV_BQC = sd(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE)/mean(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE) * 100,
        Int_CV_SPL = sd(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE)/mean(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE) * 100,

        normInt_CV_TQC = sd(.data$feature_norm_intensity[.data$qc_type == "TQC"], na.rm = TRUE)/mean(.data$feature_norm_intensity[.data$qc_type == "TQC"], na.rm = TRUE) * 100,
        normInt_CV_BQC = sd(.data$feature_norm_intensity[.data$qc_type == "BQC"], na.rm = TRUE)/mean(.data$feature_norm_intensity[.data$qc_type == "BQC"], na.rm = TRUE) * 100,
        normInt_CV_SPL = sd(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE)/mean(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE) * 100,

        conc_CV_TQC = sd(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE)/mean(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE) * 100,
        conc_CV_BQC = sd(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE)/mean(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE) * 100,
        conc_CV_SPL = sd(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE)/mean(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
        conc_CV_NIST = sd(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE)/mean(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE) * 100,
        conc_CV_LTR = sd(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE)/mean(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE) * 100,

        conc_dratio_cv_bqc = conc_CV_BQC/conc_CV_SPL,
        conc_dratio_cv_tqc = conc_CV_TQC/conc_CV_SPL,

        na_in_all_spl = all(is.na(.data$feature_conc[.data$qc_type == "SPL"])))


  data@metrics_qc <- ds1

   if ("RQC" %in% data@dataset$qc_type){
    model <- as.formula("feature_intensity ~ relative_sample_amount")
    ds2 <- data@dataset %>%
      dplyr::filter(.data$qc_type %in% c("RQC")) %>%
      dplyr::full_join(data@annot_responsecurves, by = "analysis_id") %>%
      dplyr::group_by(.data$feature_name, .data$rqc_series_id) %>%
      dplyr::filter(!all(is.na(.data$feature_intensity))) %>%
      tidyr::nest() %>%
        mutate(
          models = purrr::map(data, function(x) lm(model, data = x, na.action = na.exclude)),
          #mandel = map(data, \(x) DCVtestkit::calculate_mandel(x, "relative_sample_amount", "feature_intensity")),
          #ppa = map(data, \(x) DCVtestkit::calculate_pra_linear(x, "relative_sample_amount", "feature_intensity")),
          tidy = purrr::map(.data$models, function(x) broom::glance(x))) %>%
      tidyr::unnest(c("tidy")) %>%
      dplyr::select("feature_name", "rqc_series_id", R2 = "r.squared", Y0 = "sigma") %>%
      tidyr::pivot_wider(names_from = "rqc_series_id", values_from = c("R2", "Y0"), names_prefix = "RQC_") |>
      ungroup()

    data@metrics_qc <- data@metrics_qc  %>%
      dplyr::left_join(ds2, by = "feature_name") |>
      ungroup()
  }



  data

}


#' Save the QC table to a CSV file
#'
#' @param data MidarExperiment object
#' @param filename File name with path of exported CSV file
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider
#' @export

saveQCinfo <- function(data, filename) {

  if (nrow(data@metrics_qc)==0) stop("QC info has not yet been calculated. Please apply 'calculate_qc_metrics' first.")

  readr::write_csv(data@metrics_qc, file = filename, num_threads = 4, col_names = TRUE)
  invisible(data@metrics_qc)

}

#' Filter dataset according to QC and other criteria
#' @description
#' Filter dataset according to QC parameter criteria, remove features that are internal standards (ISTDs) or not annotated as quantifier (optional).
#' Exclude features and analyses that were annotated as not valid in the metadata (valid_integration, valid_analysis).
#'
#' Note: When `exclude_istds` is FALSE, then `min_signal_blank_ratio` is ignored, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#'
#' @param data MidarExperiment object
#' @param exclude_technical_outliers Remove samples classified as outliers
#' @param min_intensity_bqc Minimum median signal intensity of BQC
#' @param min_cv_conc_bqc = Maximum %CV of BQC
#' @param min_intensity_tqc Minimum median signal intensity of TQC
#' @param min_intensity_spl Minimum median signal intensity of study samples (SPL)
#' @param min_cv_conc_tqc Maximum %CV of TQC
#' @param max_dratio_conc_bqc D-ratio defined as CV_BQC/CV_SPL
#' @param max_dratio_conc_tqc D-ratio defined as CV_TQC/CV_SPL
#' @param min_signal_blank_ratio = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Process Blank (PBLK)
#' @param min_response_rsquare = Minimum r squared of RQC curve defined under `rqc_curve_used_for_filt`
#' @param rqc_curve_used_for_filt Name of RQC curve as string, or index number of curve to use for filtering (first curve is 1)
#' @param keep_quantifier_only Remove features where Quantifier is set to FALSE.
#' @param exclude_istds Remove Internal Standards (ISTD). If set to FALSE, meaning ISTDs will be included, then `min_signal_blank_ratio` is ignored, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#' @param features_to_keep Features that must be kept, even if they did not meet the given QC criteria
#' @return MidarExperiment object
#' @export

apply_qc_filter <-  function(data,
                             exclude_technical_outliers,
                             min_intensity_bqc = NA,
                             min_intensity_tqc = NA,
                             min_intensity_spl = NA,
                             min_signal_blank_ratio = NA,
                             min_cv_conc_bqc = NA,
                             min_cv_conc_tqc = NA,
                             max_dratio_conc_bqc = NA,
                             max_dratio_conc_tqc = NA,
                             min_response_rsquare = NA,
                             rqc_curve_used_for_filt = NA,
                             keep_quantifier_only = TRUE,
                             exclude_istds = TRUE,
                             features_to_keep = NULL) {


  if ((!is.na(min_response_rsquare)) & is.na(rqc_curve_used_for_filt)  & nrow(data@annot_responsecurves) > 0) stop("RQC Curve ID not defined! Please set rqc_curve_used_for_filt parameter or set min_response_rsquare to NA if you which not to filter based on RQC r2 values.")
  if (((!is.na(min_response_rsquare)) | !is.na(rqc_curve_used_for_filt))  & nrow(data@annot_responsecurves) == 0) stop("No RQC curves were defined in the metadata. Please reprocess with updated metadata, or to ignore linearity filtering, remove or set rqc_curve_used_for_filt and min_response_rsquare to NA")


  if (nrow(data@metrics_qc)== 0){
    stop("QC info has not yet been calculated. Please apply 'calculate_qc_metrics' first.")
    }
  if (!is.null(features_to_keep)) {
    keepers_not_defined <- setdiff( features_to_keep , unique(data@dataset$feature_name))
    txt <- glue::glue_collapse(keepers_not_defined, sep = ', ', last = ', and ')
    if (length(keepers_not_defined) > 0) stop(glue::glue("Following defined in features_to_keep are not present in this dataset: {txt}"))
    }
  data <- calculate_qc_metrics(data)  #ToDo: Run when needed

  if(is.na(min_intensity_bqc)) min_intensity_bqc <- -Inf
  if(is.na(min_intensity_tqc)) min_intensity_tqc <- -Inf
  if(is.na(min_intensity_spl)) min_intensity_spl <- -Inf
  if(is.na(min_cv_conc_bqc)) min_cv_conc_bqc <- Inf
  if(is.na(min_cv_conc_tqc)) min_cv_conc_tqc <- Inf
  if(is.na(max_dratio_conc_bqc)) max_dratio_conc_bqc <- Inf
  if(is.na(max_dratio_conc_tqc)) max_dratio_conc_tqc <- Inf
  if(is.na(min_signal_blank_ratio)) min_signal_blank_ratio <- 0
  if(is.na(min_response_rsquare)) min_response_rsquare <- 0

  # TODO: fix some of the param below ie. features_to_keep
  data@parameters_processing <- data@parameters_processing |>
    mutate(
      exclude_technical_outliers = exclude_technical_outliers,
      min_intensity_bqc = min_intensity_bqc,
      min_intensity_tqc = min_intensity_tqc,
      min_intensity_spl = min_intensity_spl,
      min_cv_conc_bqc = min_cv_conc_bqc,
      min_cv_conc_tqc = min_cv_conc_tqc,
      max_dratio_conc_bqc = max_dratio_conc_bqc,
      max_dratio_conc_tqc = max_dratio_conc_tqc,
      min_signal_blank_ratio = min_signal_blank_ratio,
      rqc_curve_used_for_filt_used_for_filt = rqc_curve_used_for_filt,
      min_response_rsquare = min_response_rsquare,
      min_response_yintersect = 0.4,
      keep_quantifier_only = keep_quantifier_only,
      exclude_istds = exclude_istds,
      features_to_keep = NA )


  #TODO: Support missing values  filtering

  data@metrics_qc <- data@metrics_qc |>
    mutate(
      pass_lod = (is.na(.data$Int_med_BQC)|.data$Int_med_BQC > min_intensity_bqc) &
        (is.na(.data$Int_med_TQC)|.data$Int_med_TQC > min_intensity_tqc) &
        (is.na(.data$Int_med_SPL)|.data$Int_med_SPL > min_intensity_spl),
      pass_sb =  (is.na(.data$SB_Ratio_median)|(.data$SB_Ratio_median > min_signal_blank_ratio|(.data$is_istd & !exclude_istds))),
      pass_cva = (is.na(.data$conc_CV_BQC)|.data$conc_CV_BQC < min_cv_conc_bqc) & (is.na(.data$conc_CV_TQC)|.data$conc_CV_TQC < min_cv_conc_tqc),
      pass_dratio = (is.na(.data$conc_dratio_cv_bqc)|.data$conc_dratio_cv_bqc < max_dratio_conc_bqc) & (is.na(.data$conc_dratio_cv_tqc)|.data$conc_dratio_cv_tqc < max_dratio_conc_tqc),
      pass_linearity = TRUE,
      pass_no_na = !(.data$na_in_all_spl)
    )


  d_filt <-  data@metrics_qc %>% filter((pass_lod & pass_sb & pass_cva) | (.data$feature_name %in% features_to_keep))

  if(is.numeric(rqc_curve_used_for_filt)) {
    rqc_r2_col_names <- names(data@metrics_qc)[which(stringr::str_detect(names(data@metrics_qc), "R2_RQC"))]
    rqc_r2_col <- rqc_r2_col_names[rqc_curve_used_for_filt]
    if(is.na(rqc_r2_col)) stop(glue::glue("RQC curve index exceeds the {length(rqc_r2_col_names)} present RQC curves in the dataset. Please resivit `rqc_curve_used_for_filt` value."))
  } else {
    rqc_r2_col <- paste0("R2_RQC_", rqc_curve_used_for_filt)
  }
  if(rqc_r2_col %in% names(data@metrics_qc)) {
    data@metrics_qc <- data@metrics_qc |>
      mutate(
        pass_linearity  = if_else(is.na(!!ensym(rqc_r2_col)), NA,  !!ensym(rqc_r2_col) > min_response_rsquare)
      )

    d_filt <- d_filt %>% filter(is.na(!!ensym(rqc_r2_col)) | pass_linearity | (.data$feature_name %in% features_to_keep))
  }


  if(exclude_istds) d_filt <- d_filt |> filter(!.data$is_istd)
  if(keep_quantifier_only) d_filt <- d_filt |> filter(.data$is_quantifier)

  d_filt <- d_filt |> filter(.data$valid_integration)


  data@metrics_qc <- data@metrics_qc |>
    mutate(
      qc_pass  = pass_no_na & pass_lod & pass_sb & pass_cva & pass_linearity
    )

  n_valid <- data@metrics_qc |> filter(.data$valid_integration)
  if(keep_quantifier_only) n_valid <- n_valid |> filter(.data$is_quantifier)
  if(exclude_istds) n_valid <- n_valid |> filter(!.data$is_istd)

  writeLines(crayon::green(glue::glue("\u2713 QC filtering applied: {nrow(d_filt)} of {nrow(n_valid)} valid features passed QC criteria")))
  data@dataset_filtered <- data@dataset %>%
    dplyr::right_join(d_filt|> dplyr::select("feature_name"), by = "feature_name") |>
    filter(.data$valid_analysis, !(.data$outlier_technical & exclude_technical_outliers))
  data
}



