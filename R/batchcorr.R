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
    writeLines(crayon::green(glue::glue("\u2713 Batch correction has been applied to drift-corrected concentrations of all {nrow(data@annot_features)} features.")))
  else
    writeLines(crayon::green(glue::glue("\u2713 Batch correction has been applied to raw concentrations of all {nrow(data@annot_features)} features.")))

  if(data@is_drift_corrected & use_raw_concentrations) message("Note: previous drift correction has been removed.\n")
  data@status_processing <- "Adjusted Quantitated Data"
  data@is_batch_corrected = TRUE
  data
}


#' Run-rder Drift Correction
#' @description
#' A short description...
#' #' @details
#' Additional details...
#' @param data MidarExperiment object
#' @param qc_types QC types used for drift correction
#' @param within_batch Correct each batch separately (Default is TRUE)
#' @param apply_conditionally Apply drift correction to all species or conditionally based on 'min_sample_cv_ratio_before_after'
#' @param log2_transform log2 transform data during correction (Default is TRUE)
#' @param loess_span Loess span width (default is 0.75)
#' @param feature_list Apply correction only to species matching (RegEx)
#' @param apply_conditionally_per_batch Apply correction conditionally using min_sample_cv_ratio_before_after criteriaper batch or across batches
#' @param min_sample_cv_ratio_before_after Maximum sample CV change for correction to be applied
#' @return MidarExperiment object
#' @export
corr_drift_loess <- function(data, qc_types, within_batch, apply_conditionally, apply_conditionally_per_batch = TRUE,
                             log2_transform = TRUE, loess_span = 0.75, feature_list = NULL, min_sample_cv_ratio_before_after = 1){

  get_loess <- function(d, qc_types,loess_span) {
    res <- tryCatch({
      stats::loess(y ~ x,
                          span = loess_span, family = "gaussian", degree = 2, normalize=FALSE, iterations=4,
                          data = d[d$QC_TYPE %in% qc_types, ]) %>%
        stats::predict(tibble::tibble(x = seq(min(d$x), max(d$x), 1))) %>% as.numeric()},
    error = function(e) {
      #message("Not enough data for the Loess fit"). # will be shown for each feature/batch...
      return(rep(NA_real_, length(d$x)))})
    list(res = res, fit_error = all(is.na(res)))
  }

  if(is.null(feature_list))
    ds <- data@dataset
  else
    ds <- data@dataset %>% dplyr::filter(stringr::str_detect(.data$FEATURE_NAME, feature_list))

  ds$x <- ds$RUN_ID
  ds$y <- ds$CONC_RAW
  if(log2_transform) ds$y <- log2(ds$y)
  if(within_batch) adj_groups <- c("FEATURE_NAME", "BATCH_ID") else adj_groups <- c("FEATURE_NAME")

  suppressWarnings(
    d <- ds %>%
      group_by(group_by(across(all_of(adj_groups)))) %>%
      nest() %>%
      mutate(
        RES = purrr::map(data, \(x) get_loess(x, qc_types, loess_span)),
        Y_PREDICTED = purrr::map(RES, \(x) x$res),
        fit_error = purrr::map(RES, \(x) x$fit_error)) |>
        # Y_PREDICTED =  MODEL_RES$res,
        # MODEL_ERR = MODEL_RES$err) |>
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

  # Calculate CV and apply to all or conditionally
  if(apply_conditionally_per_batch & within_batch) filter_groups <- c("FEATURE_NAME", "BATCH_ID") else filter_groups <- c("FEATURE_NAME")
  d <- d %>%
    #filter(!is.na(.data$Y_PREDICTED)) |>
    group_by(across(all_of(filter_groups))) %>%
    mutate(CV_RAW_SPL = sd(.data$CONC_RAW[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$CONC_RAW[.data$QC_TYPE == "SPL"], na.rm = TRUE) *100,
           CV_ADJ_SPL = sd(.data$Y_ADJ[.data$QC_TYPE == "SPL"], na.rm = TRUE)/mean(.data$Y_ADJ[.data$QC_TYPE == "SPL"], na.rm = TRUE) *100) %>%
    mutate(
      DRIFT_CORRECTED = ((( .data$CV_ADJ_SPL/.data$CV_RAW_SPL ) <  min_sample_cv_ratio_before_after) | !apply_conditionally) & !is.na(.data$Y_ADJ),
      Y_FINAL = dplyr::if_else(.data$DRIFT_CORRECTED, .data$Y_ADJ, dplyr::if_else(is.na(.data$Y_ADJ), NA_real_, .data$CONC_RAW))) |>
    ungroup()


  data@dataset <- data@dataset %>% dplyr::left_join(
    d %>% dplyr::select("ANALYSIS_ID", "FEATURE_NAME", CURVE_Y_PREDICTED = "Y_PREDICTED", CONC_DRIFT_ADJ = "Y_ADJ",
                        "CV_RAW_SPL", "CV_ADJ_SPL", "DRIFT_CORRECTED", FIT_ERROR = "fit_error", CONC_ADJ = "Y_FINAL"), by = c("ANALYSIS_ID", "FEATURE_NAME"))

  d_fit <- d |>
    group_by(FEATURE_NAME) |>
    summarise(FIT_ERROR = any(fit_error))
  fit_errors <- sum(d_fit$FIT_ERROR)

  d_sum_adj <- data@dataset |>
    filter(!is.na(.data$CURVE_Y_PREDICTED)) |>
    group_by(across(all_of(filter_groups))) |>
    summarise(DRIFT_CORRECTED = all(.data$DRIFT_CORRECTED))

  if (apply_conditionally_per_batch & within_batch){
    d_sum_adj <- d_sum_adj |>
      summarise(DRIFT_CORRECTED = any(DRIFT_CORRECTED))
  }
  if(!apply_conditionally)
    count_feature_text <- glue::glue("for {sum(d_sum_adj$DRIFT_CORRECTED, na.rm = TRUE)-fit_errors} of total {nrow(d_sum_adj)} features.")
  else
    if(apply_conditionally_per_batch & within_batch)
      count_feature_text <- glue::glue("of at least one batch for {sum(d_sum_adj$DRIFT_CORRECTED, na.rm = TRUE)-fit_errors} of total {nrow(d_sum_adj)} features.")
  else
    count_feature_text <- glue::glue("of all batches for {sum(d_sum_adj$DRIFT_CORRECTED, na.rm = TRUE)-fit_errors} of total {nrow(d_sum_adj)} features.")
  if(within_batch) mode_text <- "batch-wise" else mode_text <- "spanning all batches"
  writeLines(crayon::green(glue::glue("\u2713 Drift correction has been applied {mode_text} to raw concentrations {count_feature_text}.")))



  d_sum_span <- d |>
    group_by(across(all_of(c("ANALYSIS_ID", "BATCH_ID", "QC_TYPE")))) |>
    summarise(WITHIN_QC_SPAN = any(!is.na(.data$Y_PREDICTED)))

  n_spl_excl = sum(d_sum_span$QC_TYPE == "SPL" & !d_sum_span$WITHIN_QC_SPAN)
  n_nist_excl = sum(d_sum_span$QC_TYPE == "NIST" & !d_sum_span$WITHIN_QC_SPAN)
  n_bqc_excl = sum(d_sum_span$QC_TYPE == "BQC" & !d_sum_span$WITHIN_QC_SPAN)

  txt_1 <- txt_2 <- txt_3 <- ""
  if(n_spl_excl>0) txt_1  <- paste0(n_spl_excl, " study samples (SPL)")
  if(n_nist_excl>0) txt_2   <- paste0(n_nist_excl, " NISTs")
  if(n_bqc_excl>0) txt_3   <- paste0(n_bqc_excl, " BQCs")

  txt_final <- stringr::str_sub(paste(list(txt_1, txt_2, txt_3), collapse = ", "),1, -3)
  if(length(txt_final)>0) writeLines(crayon::red(glue::glue("Warning: {txt_final} were excluded from the correction (outside of regions spanned by QC used by LOESS).")))

  data@status_processing <- "Adjusted Quantitated Data"
  data@is_drift_corrected <- TRUE
  data@is_batch_corrected <- FALSE

  if(data@is_batch_corrected) writeLines(crayon::blue(glue::glue("Note: previous batch correction has been removed.")))
  if(fit_errors > 0) writeLines(crayon::red(glue::glue("Warning: Fit failed for {fit_errors} features (span too small and not enough data).")))
  data
}
