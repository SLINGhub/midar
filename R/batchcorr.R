##########################################################################################
# Gaussian Kernel Smoother
# based on approach and scripts by Hyung Won Choi, National University of Singapore
# Teo et al, Analytical Chemistry, 2020

fun_corr_gaussiankernel = function(data, qc_types, span_width, ...) {

  arguments <- list(...)

  d_subset <- data[data$QC_TYPE %in% qc_types, ] |> tidyr::drop_na(.data$y)
  res <- tryCatch({
    fit <- KernSmooth::locpoly(d_subset$x, d_subset$y, bandwidth = span_width, gridsize = nrow(data), range.x = c(min(data$x), max(data$x)))
    fit$y
  },
  error = function(e) {
    #print(e$message)
    return(rep(NA_real_, length(data$x)))})
  list(res = res, fit_error = all(is.na(res)))
}

fun_corr_loess <- function(d, qc_types, span_width, ...) {
  arguments <- list(...)
  surface <- ifelse(arguments$extrapolate, "direct", "interpolate")
  res <- tryCatch({
    stats::loess(y ~ x,
                 span = span_width, family = "gaussian", degree = 2, normalize=FALSE, iterations=4, surface = surface,
                 data = d[d$QC_TYPE %in% qc_types, ]) %>%
      stats::predict(tibble::tibble(x = seq(min(d$x), max(d$x), 1))) %>% as.numeric()},
    error = function(e) {
      #print(e$message) # will be shown for each feature/batch...
      return(rep(NA_real_, length(d$x)))})
  list(res = res, fit_error = all(is.na(res)))
}

#' Drift Correction by LOESS Smoothing
#' @description
#' Function to correct for run-order drifts within or across batches using loess smoothing
#' @details
#' Note that using extrapolation for loess smoothing is generally not recommended. Use this only if you must include samples or QCs that are outside of the range spanned by the QCs used for smoothing.
#' Cases where this may be necessary are when specific drifts occur in analysis sequence segments that are not spanned by QC, e.g. when the instrument broke down or suddenly changed its sensitivity.
#' Only use for samples  adjacent to the first or last QC and consult runscatter plots.
#'
#'
#' @param data MidarExperiment object
#' @param qc_types QC types used for drift correction
#' @param span Loess span width (default is 0.75)
#' @param log2_transform log2 transform data during correction (Default is TRUE). Log transformation is only used for the fit, results will not be transformed.
#' @param within_batch Correct each batch separately (Default is TRUE)
#' @param apply_conditionally Apply drift correction to all species or conditionally based on 'min_sample_cv_ratio_before_after'
#' @param apply_conditionally_per_batch Apply correction conditionally using min_sample_cv_ratio_before_after criteriaper batch or across batches
#' @param min_sample_cv_ratio_before_after Maximum sample CV change for correction to be applied
#' @param feature_list Apply correction only to species matching (RegEx).
#' @param extrapolate Extrapolate loess smoothing. WARNING: Use with casre, it is generally not recommended to extrapolate outside of the range spanned by the QCs used for smoothing. See details below.
#' @return MidarExperiment object
#' @export
corr_drift_loess <- function(data, qc_types, within_batch, apply_conditionally, apply_conditionally_per_batch = TRUE,
                           log2_transform = TRUE, span = 0.75, feature_list = NULL, min_sample_cv_ratio_before_after = 1, extrapolate = FALSE){

  corr_drift_fun(data=data, smooth_fun = "fun_corr_loess", qc_types=qc_types, within_batch=within_batch, apply_conditionally=apply_conditionally, apply_conditionally_per_batch=apply_conditionally_per_batch,
                             log2_transform=log2_transform, span_width = span, feature_list = feature_list, min_sample_cv_ratio_before_after = min_sample_cv_ratio_before_after, extrapolate = extrapolate)

}

#' Drift Correction by Gaussian Kernel Smoothing
#' @description
#' Function to correct for run-order drifts within or across batches using gaussian kernel smoothing. This is typically used to smooth based on the study samples.
#' #' @details
#'
#' @param data MidarExperiment object
#' @param qc_types QC types used for drift correction. Typically includes the study samples (`SPL`).
#' @param bandwidth Kernel bandwidth
#' @param log2_transform log2 transform data during correction (Default is TRUE). Log transformation is only used for the fit, results will not be transformed.
#' @param within_batch Correct each batch separately (Default is TRUE)
#' @param apply_conditionally Apply drift correction to all species or conditionally based on 'min_sample_cv_ratio_before_after'
#' @param apply_conditionally_per_batch Apply correction conditionally using min_sample_cv_ratio_before_after criteriaper batch or across batches
#' @param min_sample_cv_ratio_before_after Maximum sample CV change for correction to be applied
#' @param feature_list Apply correction only to species matching (RegEx)
#' @param use_uncorrected_if_fit_failed If fit function fails, use orginal (uncorrected) data or return NA for all analyses of this feature
#' @return MidarExperiment object
#' @export
corr_drift_gaussiankernel <- function(data, qc_types, bandwidth, log2_transform = TRUE, within_batch, apply_conditionally, apply_conditionally_per_batch = TRUE,
                              feature_list = NULL, min_sample_cv_ratio_before_after = 1, use_uncorrected_if_fit_failed=FALSE){

  corr_drift_fun(data=data, smooth_fun = "fun_corr_gaussiankernel", qc_types=qc_types, within_batch=within_batch, apply_conditionally=apply_conditionally, apply_conditionally_per_batch=apply_conditionally_per_batch,
                 log2_transform=log2_transform, span_width = bandwidth, feature_list = feature_list, min_sample_cv_ratio_before_after = min_sample_cv_ratio_before_after, use_uncorrected_if_fit_failed = use_uncorrected_if_fit_failed, options)

}

#' Drift Correction by Custom Function
#' @description
#' Function to correct for run-order drifts within or across batches via a provided custom function
#' #' @details
#' The drift correction function needs to be provided by the user. See `smooth_fun` for details.
#' @param data MidarExperiment object
#' @param smooth_fun Function that performs drift correction. Function need to have following parameter `data` (`MidarExperiment`), `QC_TYPES` (one or more strings), and `span_width` (numerical).
#' Function needs to return a numerical vector with the length of number of rows in `data`. In case functions fails a vector with NA_real_ needs be returned
#' @param qc_types QC types used for drift correction
#' @param span_width Width used by smoothing function
#' @param log2_transform log2 transform data during correction (Default is TRUE). Log transformation is only used for the fit, results will not be transformed.
#' @param within_batch Correct each batch separately (Default is TRUE)
#' @param apply_conditionally Apply drift correction to all species or conditionally based on 'min_sample_cv_ratio_before_after'
#' @param apply_conditionally_per_batch Apply correction conditionally using min_sample_cv_ratio_before_after criteriaper batch or across batches
#' @param min_sample_cv_ratio_before_after Maximum sample CV change for correction to be applied
#' @param use_uncorrected_if_fit_failed If TRUE, then uncorrected values are returned when fit failed, if FALSE then NA is returned.
#' @param feature_list Apply correction only to species matching (RegEx).
#' @param ... Parameters specific for the smoothing function
#' @return MidarExperiment object
#' @export
corr_drift_fun <- function(data, smooth_fun, qc_types, log2_transform = TRUE, span_width, within_batch, apply_conditionally, apply_conditionally_per_batch = TRUE,
                             min_sample_cv_ratio_before_after = 1, use_uncorrected_if_fit_failed = TRUE, feature_list = NULL, ...){


  if(is.null(feature_list))
    ds <- data@dataset
  else
    ds <- data@dataset %>% dplyr::filter(stringr::str_detect(.data$FEATURE_NAME, feature_list))
  ds$x <- ds$RUN_ID
  ds$y <- ds$CONC_RAW
  if(log2_transform) suppressWarnings(ds$y <- log2(ds$y))
  if(within_batch) adj_groups <- c("FEATURE_NAME", "BATCH_ID") else adj_groups <- c("FEATURE_NAME")
  suppressWarnings(
  d <- ds %>%
    group_by(group_by(across(all_of(adj_groups)))) %>%
    nest() %>%
    mutate(
      RES = purrr::map(data, \(x) do.call(smooth_fun,list(x, qc_types, span_width, ...))),
      Y_PREDICTED = purrr::map(.data$RES, \(x) x$res),
      fit_error = purrr::map(.data$RES, \(x) x$fit_error)) |>
    unnest(cols = c(data, .data$Y_PREDICTED, .data$fit_error))
  )

  # Get interpolated values
  if (log2_transform)
  {
    d <- d %>%
      group_by(group_by(across(all_of(adj_groups)))) %>%
      mutate(Y_PREDICTED = .data$Y_PREDICTED - median(.data$Y_PREDICTED, na.rm = TRUE),
             Y_ADJ = 2^(.data$y - .data$Y_PREDICTED))
  } else {
    d <- d %>%
      group_by(group_by(across(all_of(adj_groups)))) %>%
      mutate(Y_PREDICTED = .data$Y_PREDICTED/median(.data$Y_PREDICTED, na.rm = TRUE),
             Y_ADJ = .data$y/.data$Y_PREDICTED)
  }

  # Calculate CVs and apply to all or conditionally
  if(apply_conditionally_per_batch & within_batch) filter_groups <- c("FEATURE_NAME", "BATCH_ID") else filter_groups <- c("FEATURE_NAME")
  d <- d %>%
    #filter(!is.na(.data$Y_PREDICTED)) |>
    group_by(across(all_of(filter_groups))) %>%
    mutate(CV_RAW_SPL = sd(.data$CONC_RAW[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$CONC_RAW[.data$QC_TYPE == "SPL"], na.rm = TRUE) *100,
           CV_ADJ_SPL = sd(.data$Y_ADJ[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$Y_ADJ[.data$QC_TYPE == "SPL"], na.rm = TRUE) *100) %>%
    mutate(
      DRIFT_CORRECTED = ((( .data$CV_ADJ_SPL/.data$CV_RAW_SPL ) <  min_sample_cv_ratio_before_after) | !apply_conditionally) & !is.na(.data$Y_ADJ) & !.data$fit_error,
      Y_FINAL = dplyr::if_else(.data$DRIFT_CORRECTED, .data$Y_ADJ, dplyr::if_else(is.na(.data$Y_ADJ) & !use_uncorrected_if_fit_failed , NA_real_, .data$CONC_RAW))) |>
    ungroup()

  data@dataset <- data@dataset %>% dplyr::left_join(
    d %>% dplyr::select("ANALYSIS_ID", "FEATURE_NAME", CURVE_Y_PREDICTED = "Y_PREDICTED", CONC_DRIFT_ADJ = "Y_ADJ",
                        "CV_RAW_SPL", "CV_ADJ_SPL", "DRIFT_CORRECTED", FIT_ERROR = "fit_error", CONC_ADJ = "Y_FINAL"), by = c("ANALYSIS_ID", "FEATURE_NAME"))


  # ToDo: below if some species have invalid numbers (e.g. negative) they will be ignore, but how about if eg 90% are invalide, Still correct result?
  d_sum_adj <- data@dataset |>
    filter(!is.na(.data$CURVE_Y_PREDICTED)) |>
    group_by(across(all_of(filter_groups))) |>
    summarise(DRIFT_CORRECTED = any(.data$DRIFT_CORRECTED & !is.na(.data$CONC_RAW))|all(.data$FIT_ERROR))

  d_fit <- data@dataset |>
    group_by(.data$FEATURE_NAME) |>
    summarise(FIT_ERROR = any(.data$FIT_ERROR),
              CV_RAW_SPL = sd(.data$CONC_RAW[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$CONC_RAW[.data$QC_TYPE == "SPL"], na.rm = TRUE) *100,
              CV_ADJ_SPL = sd(.data$CONC_ADJ[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$CONC_ADJ[.data$QC_TYPE == "SPL"], na.rm = TRUE) *100)
  fit_errors <- sum(d_fit$FIT_ERROR)
  cv_median_raw <- round(median(d_fit$CV_RAW_SPL, na.rm = TRUE),2)
  cv_median_adj <- round(median(d_fit$CV_ADJ_SPL, na.rm = TRUE),2)
  features_with_fiterror <- d_fit |> filter(.data$FIT_ERROR) |> pull(.data$FEATURE_NAME)

  if(length(features_with_fiterror) > 9)
    features_with_fiterror_text <- glue::glue("{glue::glue_collapse(features_with_fiterror[1:9], sep = ", "}, and {fit_errors-9} more features.")
  else
    features_with_fiterror_text <- glue::glue_collapse(features_with_fiterror, ", ", last = " and ")

    if (apply_conditionally_per_batch & within_batch){
    d_sum_adj <- d_sum_adj |>
      summarise(DRIFT_CORRECTED = any(.data$DRIFT_CORRECTED))
  }
  if(!apply_conditionally)
    count_feature_text <- glue::glue("to {sum(d_sum_adj$DRIFT_CORRECTED, na.rm = TRUE)} of {nrow(d_sum_adj)} features.")
  else
    if(apply_conditionally_per_batch & within_batch)
      count_feature_text <- glue::glue("of at least one batch for {sum(d_sum_adj$DRIFT_CORRECTED, na.rm = TRUE)} of total {nrow(d_sum_adj)} features.")
  else
    count_feature_text <- glue::glue("to {sum(d_sum_adj$DRIFT_CORRECTED, na.rm = TRUE)} of {nrow(d_sum_adj)} features.")
  if(within_batch) mode_text <- "was applied batch-wise" else mode_text <- "across all batches was applied"

  writeLines(crayon::green(glue::glue("\u2713 Drift correction {mode_text} to raw concentrations {count_feature_text} Median study-sample CV of all features before and after correction: {cv_median_raw}% and {cv_median_adj}%.")))



  d_sum_span <- d |>
    group_by(across(all_of(c("ANALYSIS_ID", "BATCH_ID", "QC_TYPE")))) |>
    summarise(WITHIN_QC_SPAN = any(!is.na(.data$Y_PREDICTED)))

  n_spl_excl = sum(d_sum_span$QC_TYPE == "SPL" & !d_sum_span$WITHIN_QC_SPAN)
  n_nist_excl = sum(d_sum_span$QC_TYPE == "NIST" & !d_sum_span$WITHIN_QC_SPAN)
  n_bqc_excl = sum(d_sum_span$QC_TYPE == "BQC" & !d_sum_span$WITHIN_QC_SPAN)
  n_tqc_excl = sum(d_sum_span$QC_TYPE == "TQC" & !d_sum_span$WITHIN_QC_SPAN)
  n_ltr_excl = sum(d_sum_span$QC_TYPE == "LTR" & !d_sum_span$WITHIN_QC_SPAN)

  txt_1 <- txt_2 <- txt_3 <- txt_4 <- txt_5 <- character()
  if(n_spl_excl>0) txt_1  <- paste0(n_spl_excl, " of ", sum(d_sum_span$QC_TYPE == "SPL"), " study samples (SPL)")
  if(n_nist_excl>0) txt_2   <- paste0(n_nist_excl, " of ", sum(d_sum_span$QC_TYPE == "NIST"), " NISTs")
  if(n_bqc_excl>0) txt_3   <- paste0(n_bqc_excl,  " of ", sum(d_sum_span$QC_TYPE == "BQC"), " BQCs")
  if(n_tqc_excl>0) txt_4   <- paste0(n_tqc_excl, " of ", sum(d_sum_span$QC_TYPE == "TQC"), " TQCs")
  if(n_ltr_excl>0) txt_5   <- paste0(n_ltr_excl, " of ", sum(d_sum_span$QC_TYPE == "LTR"), " LTRs")

  txt_final <- paste(c(txt_1, txt_2, txt_3, txt_4, txt_5), collapse = ", ")
  if(txt_final != "") writeLines(crayon::yellow(glue::glue("Warning: {txt_final} excluded from correction (beyond regions spanned by QCs).")))

  data@status_processing <- "Adjusted Quantitated Data"
  data@is_drift_corrected <- TRUE
  data@is_batch_corrected <- FALSE

  if(data@is_batch_corrected) writeLines(crayon::yellow(glue::glue("Note: previous batch correction has been removed.")))
  if(fit_errors > 0) writeLines(crayon::yellow(glue::glue("Warning: No smoothing applied for {fit_errors} features because the fit algorithm failed (insufficient or invalid data points): {features_with_fiterror_text}")))
  data@dataset$Concentration <- data@dataset$CONC_ADJ
  data
}



#' Batch centering
#' @description
#' A short description...
#' #' @details
#' Additional details...
#' @param data MidarExperiment object
#' @param qc_types QC types used for batch correction
#' @param use_raw_concentrations Apply to unadjusted (raw) concentration. Default is FALSE, which means previously drift-corrected concentrations will be used if available, otherwise unadjusted concentrations will be used
#' @param center_fun Function used to center. Default is "median".
#' @return MidarExperiment object
#' @importFrom glue glue
#' @export
corr_batch_centering <- function(data, qc_types, use_raw_concentrations = FALSE, center_fun = "median"){

  ds <- data@dataset
  if (!data@is_drift_corrected | use_raw_concentrations) var <- rlang::sym("CONC_RAW") else var <- rlang::sym("CONC_ADJ")
  # Normalize by the median (or user-defined function)
  ds <- ds %>%
    dplyr::group_by(.data$FEATURE_NAME,  .data$BATCH_ID) %>%
    dplyr::mutate(CONC_ADJ_NEW = {{var}}/do.call(center_fun,list(({{var}}[.data$QC_TYPE %in% qc_types]), na.rm = TRUE))) |>
    dplyr::ungroup()

  # Re-level data to the median of all batches
  ds <- ds %>%
    dplyr::group_by(.data$FEATURE_NAME) %>%
    dplyr::mutate(CONC_ADJ =  .data$CONC_ADJ_NEW * do.call(center_fun,list({{var}}[.data$QC_TYPE %in% qc_types], na.rm = TRUE))) |>
    dplyr::select(-"CONC_ADJ_NEW") |>
    dplyr::ungroup()

  data@dataset <- ds
  if(data@is_drift_corrected)
    writeLines(crayon::green(glue::glue("\u2713 Batch correction was applied to drift-corrected concentrations of all {nrow(data@annot_features)} features.")))
  else
    writeLines(crayon::green(glue::glue("\u2713 Batch correction was applied to raw concentrations of all {nrow(data@annot_features)} features.")))

  if(data@is_drift_corrected & use_raw_concentrations) writeLines(crayon::yellow(glue::glue("Note: previous drift correction has been removed.\n")))
  data@status_processing <- "Adjusted Quantitated Data"
  data@is_batch_corrected = TRUE
  data@dataset$Concentration <- data@dataset$CONC_ADJ
  data
}
