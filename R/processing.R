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
        conc_CV_LTR = sd(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE)/mean(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE) * 100)

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
      tidyr::pivot_wider(names_from = "rqc_series_id", values_from = c("R2", "Y0"), names_prefix = "RQC_")

    data@metrics_qc <- data@metrics_qc  %>% dplyr::left_join(ds2, by = "feature_name")
  }

  if (tolower(data@analysis_type) == "lipidomics") data@metrics_qc <- data@metrics_qc |> add_lipid_class_transition()
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
#' Note: When `exclude_istds` is FALSE, then `SB_RATIO_min` is ignored, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#'
#' @param data MidarExperiment object
#' @param exclude_technical_outliers Remove samples classified as outliers
#' @param Intensity_BQC_min Minimum median signal intensity of BQC
#' @param CV_BQC_max = Maximum %CV of BQC
#' @param Intensity_TQC_min Minimum median signal intensity of TQC
#' @param CV_TQC_max Maximum %CV of TQC
#' @param SB_RATIO_min = Signal-to-Blank ratio. Calculated from the 10% percentile of all samples and the median of the Process Blank (PBLK)
#' @param R2_min = Minimum r squared of RQC curve defined under `RQC_CURVE`
#' @param RQC_CURVE Name of RQC curve as string, or index number of curve to use for filtering (first curve is 1)
#' @param quantifier_only Remove features where Quantifier is set to FALSE.
#' @param exclude_istds Remove Internal Standards (ISTD). If set to FALSE, meaning ISTDs will be included, then `SB_RATIO_min` is ignored, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#' @param features_to_keep Features that must be kept, even if they did not meet the given QC criteria
#' @return MidarExperiment object
#' @export

apply_qc_filter <-  function(data,
                             exclude_technical_outliers,
                             Intensity_BQC_min = NA,
                             CV_BQC_max = NA,
                             Intensity_TQC_min = NA,
                             CV_TQC_max = NA,
                             SB_RATIO_min = NA,
                             R2_min = NA,
                             RQC_CURVE = NA,
                             quantifier_only = TRUE,
                             exclude_istds = TRUE,
                             features_to_keep = NULL) {


  if ((!is.na(R2_min)) & is.na(RQC_CURVE)  & nrow(data@annot_responsecurves) > 0) stop("RQC Curve ID not defined! Please set RQC_CURVE parameter or set R2_min to NA if you which not to filter based on RQC r2 values.")
  if (((!is.na(R2_min)) | !is.na(RQC_CURVE))  & nrow(data@annot_responsecurves) == 0) stop("No RQC curves were defined in the metadata. Please reprocess with updated metadata, or to ignore linearity filtering, remove or set RQC_CURVE and R2_min to NA")


  if (nrow(data@metrics_qc)== 0){
    stop("QC info has not yet been calculated. Please apply 'calculate_qc_metrics' first.")
    }
  if (!is.null(features_to_keep)) {
    keepers_not_defined <- setdiff( features_to_keep , unique(data@dataset$feature_name))
    txt <- glue::glue_collapse(keepers_not_defined, sep = ', ', last = ', and ')
    if (length(keepers_not_defined) > 0) stop(glue::glue("Following defined in features_to_keep are not present in this dataset: {txt}"))
    }
  data <- calculate_qc_metrics(data)  #ToDo: Run when needed

  if(is.na(Intensity_BQC_min)) Intensity_BQC_min <- -Inf
  if(is.na(CV_BQC_max)) CV_BQC_max <- Inf
  if(is.na(Intensity_TQC_min)) Intensity_TQC_min <- -Inf
  if(is.na(CV_TQC_max)) CV_TQC_max <- Inf
  if(is.na(SB_RATIO_min)) SB_RATIO_min <- 0
  if(is.na(R2_min)) R2_min <- 0

  d_filt <-  data@metrics_qc %>% filter((
                                  (is.na(.data$Int_med_BQC)|.data$Int_med_BQC > Intensity_BQC_min) &
                                  (is.na(.data$Int_med_TQC)|.data$Int_med_TQC > Intensity_TQC_min) &
                                  (is.na(.data$conc_CV_BQC)|.data$conc_CV_BQC < CV_BQC_max) &
                                  (is.na(.data$conc_CV_TQC)|.data$conc_CV_TQC < CV_TQC_max) &
                                  (is.na(.data$SB_Ratio_median)|(.data$SB_Ratio_median > SB_RATIO_min|(.data$is_istd & !exclude_istds))))|
                                  (.data$feature_name %in% features_to_keep))

  if(is.numeric(RQC_CURVE)) {
    rqc_r2_col_names <- names(data@metrics_qc)[which(stringr::str_detect(names(data@metrics_qc), "R2_RQC"))]
    rqc_r2_col <- rqc_r2_col_names[RQC_CURVE]
    if(is.na(rqc_r2_col)) stop(glue::glue("RQC curve index exceeds the {length(rqc_r2_col_names)} present RQC curves in the dataset. Please resivit `RQC_CURVE` value."))
  } else {
    rqc_r2_col <- paste0("R2_RQC_", RQC_CURVE)
  }
  if(rqc_r2_col %in% names(data@metrics_qc)) d_filt <- d_filt %>% filter(is.na(!!ensym(rqc_r2_col)) | !!ensym(rqc_r2_col) > R2_min | (.data$feature_name %in% features_to_keep))
  if(exclude_istds) d_filt <- d_filt |> filter(!.data$is_istd)
  if(quantifier_only) d_filt <- d_filt |> filter(.data$is_quantifier)

  d_filt <- d_filt |> filter(.data$valid_integration)

  n_valid <- data@metrics_qc |> filter(.data$valid_integration)
  if(quantifier_only) n_valid <- n_valid |> filter(.data$is_quantifier)
  if(exclude_istds) n_valid <- n_valid |> filter(!.data$is_istd)

  writeLines(crayon::green(glue::glue("\u2713 QC filtering applied: {nrow(d_filt)} of {nrow(n_valid)} valid features passed QC criteria")))
  data@dataset_QC_filtered <- data@dataset %>%
    dplyr::right_join(d_filt|> dplyr::select("feature_name"), by = "feature_name") |>
    filter(.data$valid_analysis, !(.data$outlier_technical & exclude_technical_outliers))
  data
}



