#' get_conc_unit
#'
#' @param sample_amount_unit MidarExperiment object
#' @return string with concentration unit
#' @noRd

get_conc_unit <- function(sample_amount_unit){
  if (length(unique(sample_amount_unit)) > 1)
    conc_unit <- "pmol/sample amount unit (multiple units)"
  else if (sample_amount_unit[1] == "uL" | sample_amount_unit[1] == "\U003BCL")
    conc_unit <- "\U003BCmol/L"
  else
    conc_unit <- glue::glue("pmol/{sample_amount_unit}")
  conc_unit
}



#' Normalize Intensities with corresponding ISTD Intensities
#'
#' @param data MidarExperiment object
#' @return MidarExperiment object
#' @importFrom glue glue
#' @export

normalize_by_istd <- function(data) {
  if(nrow(data@annot_features) < 1) stop("ISTD map is missing...please import transition annotations.")
  if("normIntensity" %in% names(data@dataset)) {
    if (!all(is.na(data@dataset$normIntensity))) warning("Overwriting exiting normalized Intensities")
    data@dataset <- data@dataset %>% select(-dplyr::any_of(c("normIntensity", "pmol_total", "Concentration", "CONC_DRIFT_ADJ", "CONC_FINAL")))
  }

  d_temp <- data@dataset #%>%     dplyr::full_join(data@annot_features, by = c("FEATURE_NAME" = "FEATURE_NAME"),)
  d_temp <- d_temp  %>%
    dplyr::group_by(.data$NORM_ISTD_FEATURE_NAME, .data$ANALYSIS_ID) %>%
    dplyr::mutate(normIntensity = .data$Intensity/.data$Intensity[.data$isISTD]) %>%
    dplyr::ungroup()

  data@dataset <- data@dataset %>%
    dplyr::inner_join(d_temp %>% dplyr::select("ANALYSIS_ID", "FEATURE_NAME", "isISTD", "normIntensity"), by = c("ANALYSIS_ID", "FEATURE_NAME"))

  n_features <- length(unique(data@annot_features$FEATURE_NAME))
  n_ISTDs <- length(unique(data@annot_features$NORM_ISTD_FEATURE_NAME))
  print(glue::glue("{n_features} features normalized with {n_ISTDs} ISTDs. \n"))
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

  if(!(c("normIntensity") %in% names(data@dataset))) stop("Data needs first to be ISTD normalized. Please apply function 'normalize_by_istd' first.")
  d_temp <- data@dataset %>%
    dplyr::left_join(data@annot_analyses %>% dplyr::select("ANALYSIS_ID", "SAMPLE_AMOUNT", "ISTD_VOL"), by = c("ANALYSIS_ID")) %>%
    dplyr::left_join(data@annot_features %>% dplyr::select("FEATURE_NAME", "QUANT_ISTD_FEATURE_NAME"), by = c("FEATURE_NAME")) %>%
    dplyr::left_join(data@annot_istd, by = c("QUANT_ISTD_FEATURE_NAME"))

  d_temp <- d_temp %>% mutate(pmol_total = (.data$normIntensity)*(.data$ISTD_VOL*(.data$ISTD_CONC_nM))/1000)
  d_temp <- d_temp %>% mutate(Concentration = .data$pmol_total/.data$SAMPLE_AMOUNT)

  if("Concentration" %in% names(data@dataset)) {
    data@dataset <- data@dataset %>% select(-dplyr::any_of(c("pmol_total", "Concentration")))
    warning("Overwriting exiting Concentrations")
  }

  data@dataset <- data@dataset %>%
    dplyr::inner_join(d_temp %>% dplyr::select("ANALYSIS_ID", "FEATURE_NAME", "pmol_total", "Concentration"), by = c("ANALYSIS_ID", "FEATURE_NAME"))

  n_features <- length(unique(data@annot_features$FEATURE_NAME))
  n_ISTDs <- length(unique(data@annot_features$NORM_ISTD_FEATURE_NAME))

  conc_unit <- get_conc_unit(data@annot_analyses$SAMPLE_AMOUNT_UNIT)

  print(glue::glue("{n_features} compounds quantitated in {nrow(data@annot_analyses)} samples using {n_ISTDs} spiked ISTDs.
                   Concentration unit: [{conc_unit}]"))

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

  if (!(variable %in% names(data@dataset))) stop("Variable '", variable,  "' does not (yet) exist in dataset")

  ds <- data@dataset |>
    dplyr::select("ANALYSIS_ID", "QC_TYPE", "AcqTimeStamp", "FEATURE_NAME", !!var) %>%
    tidyr::pivot_wider(names_from = .data$FEATURE_NAME, values_from = !!var)

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

  #if(!(c("normIntensity") %in% names(data@dataset))) warning("No normali is not normalized")

  ds1 <- data@dataset %>%
    dplyr::filter(.data$QC_TYPE %in% c("SPL", "NIST", "LTR", "BQC", "TQC", "PBLK", "SBLK", "UBLK")) %>%
    dplyr::group_by(.data$FEATURE_NAME) %>%
    dplyr::summarise(
      #PrecursorMz = paste0(unique(.data$PRECURSOR_MZ), collapse = ","),
      #ProductMz = paste0(unique(.data$PRODUCT_MZ), collapse = ","),
      Quantifier = unique(.data$isQUANTIFIER),
      ISTD = unique(.data$NORM_ISTD_FEATURE_NAME),
      Int_med_PBLK = median(.data$Intensity[.data$QC_TYPE == "PBLK"], na.rm = TRUE),
      Int_med_SPL = median(.data$Intensity[.data$QC_TYPE == "SPL"], na.rm = TRUE),
      Int_med_BQC = median(.data$Intensity[.data$QC_TYPE == "BQC"], na.rm = TRUE),
      Int_med_TQC = median(.data$Intensity[.data$QC_TYPE == "TQC"], na.rm = TRUE),
      Int_med_NIST = median(.data$Intensity[.data$QC_TYPE == "NIST"], na.rm = TRUE),
      Int_med_LTR = median(.data$Intensity[.data$QC_TYPE == "LTR"], na.rm = TRUE),

      conc_median_TQC = median(.data$Concentration[.data$QC_TYPE == "TQC"], na.rm = TRUE),
      conc_median_BQC = median(.data$Concentration[.data$QC_TYPE == "BQC"], na.rm = TRUE),
      conc_median_SPL = median(.data$Concentration[.data$QC_TYPE == "SPL"], na.rm = TRUE),
      conc_median_NIST = median(.data$Concentration[.data$QC_TYPE == "NIST"], na.rm = TRUE),
      conc_median_LTR = median(.data$Concentration[.data$QC_TYPE == "LTR"], na.rm = TRUE),

      SB_Ratio_Q10 = quantile(.data$Intensity[.data$QC_TYPE == "SPL"], probs  = 0.1, na.rm = TRUE, names = FALSE)/median(.data$Intensity[.data$QC_TYPE == "PBLK"], na.rm = TRUE, names = FALSE),

      Int_CV_TQC = sd(.data$Intensity[.data$QC_TYPE == "TQC"], na.rm = TRUE)/mean(.data$Intensity[.data$QC_TYPE == "TQC"], na.rm = TRUE) * 100,
      Int_CV_BQC = sd(.data$Intensity[.data$QC_TYPE == "BQC"], na.rm = TRUE)/mean(.data$Intensity[.data$QC_TYPE == "BQC"], na.rm = TRUE) * 100,
      Int_CV_SPL = sd(.data$Intensity[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$Intensity[.data$QC_TYPE == "SPL"], na.rm = TRUE) * 100,

      normInt_CV_TQC = sd(.data$normIntensity[.data$QC_TYPE == "TQC"], na.rm = TRUE)/mean(.data$normIntensity[.data$QC_TYPE == "TQC"], na.rm = TRUE) * 100,
      normInt_CV_BQC = sd(.data$normIntensity[.data$QC_TYPE == "BQC"], na.rm = TRUE)/mean(.data$normIntensity[.data$QC_TYPE == "BQC"], na.rm = TRUE) * 100,
      normInt_CV_SPL = sd(.data$normIntensity[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$normIntensity[.data$QC_TYPE == "SPL"], na.rm = TRUE) * 100,

      conc_CV_TQC = sd(.data$Concentration[.data$QC_TYPE == "TQC"], na.rm = TRUE)/mean(.data$Concentration[.data$QC_TYPE == "TQC"], na.rm = TRUE) * 100,
      conc_CV_BQC = sd(.data$Concentration[.data$QC_TYPE == "BQC"], na.rm = TRUE)/mean(.data$Concentration[.data$QC_TYPE == "BQC"], na.rm = TRUE) * 100,
      conc_CV_SPL = sd(.data$Concentration[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$Concentration[.data$QC_TYPE == "SPL"], na.rm = TRUE) * 100,
      conc_CV_NIST = sd(.data$Concentration[.data$QC_TYPE == "NIST"], na.rm = TRUE)/mean(.data$Concentration[.data$QC_TYPE == "NIST"], na.rm = TRUE) * 100,
      conc_CV_LTR = sd(.data$Concentration[.data$QC_TYPE == "LTR"], na.rm = TRUE)/mean(.data$Concentration[.data$QC_TYPE == "LTR"], na.rm = TRUE) * 100)

  data@metrics_qc <- ds1

  if ("RQC" %in% data@dataset$QC_TYPE){
    model <- as.formula("Intensity ~ RELATIVE_SAMPLE_AMOUNT")
    ds2 <- data@dataset %>%
      dplyr::filter(.data$QC_TYPE %in% c("RQC")) %>%
      dplyr::full_join(data@annot_responsecurves, by = "ANALYSIS_ID") %>%
      dplyr::group_by(.data$FEATURE_NAME, .data$RQC_SERIES_ID) %>%
      dplyr::filter(!all(is.na(.data$Intensity))) %>%
      tidyr::nest() %>%
        mutate(
          models = purrr::map(data, function(x) lm(model, data = x, na.action = na.exclude)),
          #mandel = map(data, \(x) DCVtestkit::calculate_mandel(x, "RELATIVE_SAMPLE_AMOUNT", "Intensity")),
          #ppa = map(data, \(x) DCVtestkit::calculate_pra_linear(x, "RELATIVE_SAMPLE_AMOUNT", "Intensity")),
          tidy = purrr::map(.data$models, function(x) broom::glance(x))) %>%
      tidyr::unnest(c("tidy")) %>%
      dplyr::select("FEATURE_NAME", "RQC_SERIES_ID", R2 = "r.squared", Y0 = "sigma") %>%
      tidyr::pivot_wider(names_from = "RQC_SERIES_ID", values_from = c("R2", "Y0"), names_prefix = "RQC_")

    data@metrics_qc <- data@metrics_qc  %>% dplyr::left_join(ds2, by = "FEATURE_NAME")
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

  if (nrow(data@metrics_qc)== 0) stop("QC info has not yet been calculated. Please apply 'calculate_qc_metrics' first.")

  readr::write_csv(data@metrics_qc, file = filename, num_threads = 4, col_names = TRUE)
  invisible(data@metrics_qc)

}

#' QC-filter dataset
#'
#' @param data MidarExperiment object
#' @param Intensity_BQC_min Intensity_BQC_min
#' @param CV_BQC_max = CV_BQC_max
#' @param Intensity_TQC_min = Intensity_TQC_min
#' @param CV_TQC_max = CV_TQC_max
#' @param SB_RATIO_min = SB_RATIO_min
#' @param R2_min = R2_min
#' @param RQC_CURVE = RQC_CURVE
#' @return MidarExperiment object
#' @export

apply_qc_filter <-  function(data,
                   Intensity_BQC_min = NA,
                   CV_BQC_max = NA,
                   Intensity_TQC_min = NA,
                   CV_TQC_max = NA,
                   SB_RATIO_min = NA,
                   R2_min = NA,
                   RQC_CURVE = NA
                   ){


  if (nrow(data@metrics_qc)== 0) stop("QC info has not yet been calculated. Please apply 'calculate_qc_metrics' first.")

  if(is.na(Intensity_BQC_min)) Intensity_BQC_min <- -Inf
  if(is.na(CV_BQC_max)) CV_BQC_max <- Inf
  if(is.na(Intensity_TQC_min)) Intensity_TQC_min <- -Inf
  if(is.na(CV_TQC_max)) CV_TQC_max <- Inf
  if(is.na(SB_RATIO_min)) SB_RATIO_min <- 0
  if(is.na(R2_min)) R2_min <- 0

  d_filt <-  data@metrics_qc %>% filter(is.na(.data$Int_med_BQC)|.data$Int_med_BQC > Intensity_BQC_min,
                                  is.na(.data$Int_med_TQC)|.data$Int_med_TQC > Intensity_TQC_min,
                                  is.na(.data$conc_CV_BQC)|.data$conc_CV_BQC < CV_BQC_max,
                                  is.na(.data$conc_CV_TQC)|.data$conc_CV_TQC < CV_TQC_max,
                                  is.na(.data$SB_Ratio_Q10)|.data$SB_Ratio_Q10 > SB_RATIO_min)


  if ((!is.na(R2_min))&is.na(RQC_CURVE)) stop("RQC Curve ID not defined! Please set RQC_CURVE parameter or set R2_min to NA if you which not to filter based on RQC r2 values")
  if(is.numeric(RQC_CURVE)) {
    rqc_r2_col <- names(data@metrics_qc)[which(stringr::str_detect(names(data@metrics_qc), "R2_RQC"))[RQC_CURVE]]
  } else {
    rqc_r2_col <- paste0("R2_RQC_", RQC_CURVE)
  }

  if(rqc_r2_col %in% names(data@metrics_qc)) d_filt <- d_filt %>% filter(is.na(!!ensym(rqc_r2_col))|!!ensym(rqc_r2_col) > R2_min)

  print(glue::glue("{nrow(d_filt)} of {nrow(data@metrics_qc)} features passed QC filtering."))
  data@dataset_QC_filtered <- data@dataset %>% dplyr::right_join(d_filt|> dplyr::select("FEATURE_NAME"), by = "FEATURE_NAME")
  data
}



