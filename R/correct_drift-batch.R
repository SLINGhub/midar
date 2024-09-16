#' Gaussian Kernel smoothing helper function
#' @description
#' Function for Gaussian kernel-based smoothing, for use by `corr_drift_fun`.
#' @author Hyung Won Choi
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param qc_types QC types used for the smoothing (fit) by loess
#' @param kernel_size Bandwidth of the gaussian kernel function
#' @param ... Additional parameters
#' @return List with a `data.frame` containing original `x` and the smoothed `y` values, and a `boolean` value indicting whether the fit failed or not not.

### Function to clean up serial trend within a batch
fun_gauss.kernel.smooth = function(tbl, qc_types, ...) {
  arg <- list(...)

  xx <- seq(1, length(tbl$x))
  yy <- tbl$y
  s.train <- tbl$qc_type %in% qc_types # Use sample for training

  n = length(yy) ## number of data points


  res <- tryCatch(
    {
      ## If outlier filter is turned on, mark an outlier as NA
      ## ksd = K times standard deviation of data distribution
      yy.train <- yy - median(yy, na.rm=TRUE)

      if(arg$outlier_filter) {
        mean.yy <- mean(yy.train, na.rm=TRUE)
        sd.yy <- sd(yy.train, na.rm=TRUE)
        oid <- (abs(yy.train - mean.yy) / sd.yy) > arg$outlier_ksd
        oid[is.na(oid)] <- FALSE
        yy.train[oid] <- NA
      }

      yy.train[!s.train] <- NA

      yy.est = yy
      ## Location parameter smoothing
      if(arg$location_smooth) {
        for(i in 1:n) {
          wt <- (xx - xx[i]) / arg$kernel_size
          wt <- dnorm(wt, 0, 1)
          wt[is.na(yy.train)] <- NA
          yy.est[i] <- yy[i] - sum(wt * yy.train, na.rm=TRUE) / sum(wt, na.rm=TRUE)
        }
      }


      ## If scale parameter smoothing is on, mean-center data,
      ## estimate point-wise weighted stdev's and scale them
      ## and add back the overall mean
      yy.final <- yy.est

      if(arg$scale_smooth) {
        v <- rep(NA, n) ## point-wise weighted variances
        yy.mean <- mean(yy.est, na.rm=TRUE)
        yy.est <- yy.est - yy.mean
        for(i in 1:n) {
          if(!is.na(yy.est[i])) {
            wt <- (xx - xx[i]) / arg$kernel_size
            wt <- dnorm(wt, 0, 1)
            wt[is.na(yy.est)] <- NA
            v[i] <- sum(wt * yy.est^2, na.rm=TRUE) / sum(wt, na.rm=TRUE)
          }
        }
        v.mean <- mean(v, na.rm=TRUE) ## average weighted variances across the data points
        s <- sqrt(v)
        s.mean <- mean(s, na.rm=TRUE)
        yy.final <- yy.mean + yy.est * s.mean / s
      }
      yy.final
    } ,
    error = function(e) {
      print(e$message)
      return(rep(NA_real_, length(n)))
    }
  )
  ## Report the final values
  list(res = res, fit_error = all(is.na(res)), data_adjusted = TRUE)
}

#' Gaussian Kernel smoothing helper function
#' @description
#' Function for Gaussian kernel-based smoothing, for use by `corr_drift_fun`
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param qc_types QC types used for the smoothing (fit) by loess
#' @param span_width Bandwidth of the gaussian kernel function
#' @param ... Additional parameters forwarded to KernSmooth::locpoly
#' @return List with a `data.frame` containing original `x` and the smoothed `y` values, and a `boolean` value indicting whether the fit failed or not not.

fun_gaussiankernel_old <- function(tbl, qc_types, ...) {
  arg <- list(...)

  # browser()
  d_subset <- tbl[tbl$qc_type %in% qc_types, ] |> tidyr::drop_na(.data$y)
  res <- tryCatch(
    {
      fit <- KernSmooth::locpoly(d_subset$x, d_subset$y, bandwidth = arg$span_width, gridsize = nrow(tbl), range.x = c(min(tbl$x), max(tbl$x)))
      fit$y
    },
    error = function(e) {
      print(e$message)
      return(rep(NA_real_, length(tbl$x)))
    }
  )
  list(res = res, fit_error = all(is.na(res)))
}


#' Loess smoothing helper function
#' @description
#' Function for loess-based smoothing, for use by `corr_drift_fun`
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param qc_types QC types used for the smoothing (fit) by loess
#' @param span_width Loess span width
#' @param ... Additional parameters forwarded to Loess
#' @return List with a `data.frame` containing original x and the smoothed y values, and a `boolean` value indicting whether the fit failed or not not.
fun_loess <- function(tbl, qc_types,...) {
  #tbl <- tbl |> filter(!.data$outlier_technical)
  arg <- list(...)
  surface <- ifelse(arguments$extrapolate, "direct", "interpolate")
  res <- tryCatch(
    {
      stats::loess(y ~ x,
        span = arg$span_width, family = "symmetric", degree = 2, normalize = FALSE, iterations = 4, surface = surface, ...,
        data = tbl[tbl$qc_type %in% qc_types, ]
      ) %>%
        stats::predict(tibble::tibble(x = seq(min(tbl$x), max(tbl$x), 1))) %>%
        as.numeric()
    },
    error = function(e) {
      # print(e$message) # will be shown for each feature/batch...
      return(rep(NA_real_, length(tbl$x)))
    }
  )
  list(res = res, fit_error = all(is.na(res)))
}


#' Drift Correction by Custom Function
#' @description
#' Function to correct for run-order drifts within or across batches via a provided custom function
#' #' @details
#' The drift correction function needs to be provided by the user. See `smooth_fun` for details.
#' @param data MidarExperiment object
#' @param smooth_fun Function that performs drift correction. Function need to have following parameter `data` (`MidarExperiment`), `qc_typeS` (one or more strings), and `span_width` (numerical).
#' Function needs to return a numerical vector with the length of number of rows in `data`. In case functions fails a vector with NA_real_ needs be returned
#' @param qc_types QC types used for drift correction
#' @param within_batch Apply to each batch separately if `TRUE` (the default)
#' @param log2_transform Log transform the data for correction when `TRUE` (the default). Note: log transformation is solely applied internally for smoothing, results will not be be log-transformed. Log transformation may result in more robust smoothing that is less sensitive to outlier.
#' @param apply_conditionally Apply drift correction to all species if `TRUE`, or only when sample CV after smoothing changes below a threshold defined via `max_cv_ratio_before_after`
#' @param apply_conditionally_per_batch When `apply_conditionally = TRUE`, correction is conditionally applied per batch when `TRUE` and across all batches when `FALSE`
#' @param max_cv_ratio_before_after Only used when `apply_conditionally = TRUE`. Maximum allowed ratio of sample CV change before and after smoothing for the correction to be applied.
#' A value of 1 (the default) indicates the CV needs to improve or remain unchanged after smoothing so that the conditional smoothing is applied. A value of < 1 means that CV needs to improve, a value of e.g. 1.20 that the CV need to improve or get worse by max 1.20-fold after smoothing.
#' @param feature_list Subset the features for correction whose names matches the specified text using regular expression. Default is `NULL` which means all features are selected.
#' @param use_uncorrected_if_fail In case the smoothing function fails for a species, then use original (uncorrected) data when `TRUE` (the default) or return `NA` for all analyses of the feature where the fit failed.
#' @param ... Parameters specific for the smoothing function
#' @return MidarExperiment object
#' @export
corr_drift_fun <- function(data, smooth_fun, qc_types, log2_transform = TRUE, within_batch, apply_conditionally, apply_conditionally_per_batch = TRUE,
                           max_cv_ratio_before_after = 1, use_uncorrected_if_fail = TRUE, feature_list = NULL, ...) {
  if (is.null(feature_list)) {
    ds <- data@dataset
  } else {
    ds <- data@dataset %>% dplyr::filter(stringr::str_detect(.data$feature_id, feature_list))
  }

  ds$y <- ds$conc_raw
  ds <- ds |> mutate(x = dplyr::row_number(), .by = "feature_id")

  if (log2_transform) suppressWarnings(ds$y <- log2(ds$y))
  if (within_batch) adj_groups <- c("feature_id", "batch_id") else adj_groups <- c("feature_id")


  d_smooth <- ds %>%
      select("analysis_id", "qc_type", "feature_id", "batch_id", "conc_raw", "x", "y") |>
      filter(.data$qc_type %in% c("SPL", "BQC", "TQC", "NIST", "LTR", "PBLK", "SBPK", "RQC")) |>
      group_by(across(all_of(adj_groups))) %>%
      nest() %>%
      mutate(
        RES = purrr::map(data, \(x) do.call(smooth_fun, list(x, qc_types, ...))),
        Y_PREDICTED = purrr::map(.data$RES, \(x) x$res),
        fit_error = purrr::map(.data$RES, \(x) x$fit_error),
        data_adjusted = purrr::map(.data$RES, \(x) x$data_adjusted)
      ) |>
      unnest(cols = c(data, .data$Y_PREDICTED, .data$fit_error,.data$data_adjusted))


  # Get data and backtransform if log-scaled.
  # TODODOD: split calc of Y_MEDIAN, dont add same table

  adjusted <- all(d_smooth$data_adjusted) #TODOTODO

  if (log2_transform) {
    if(!adjusted){
      d_smooth <- d_smooth %>%
        group_by(across(all_of(adj_groups))) %>%
        mutate(
          Y_MEDIAN = median(2^(.data$Y_PREDICTED), na.rm = TRUE),
          Y_ADJ = 2^(.data$y - .data$Y_PREDICTED + median(.data$Y_PREDICTED, na.rm = TRUE))
        )
    } else {
      d_smooth <- d_smooth %>%
        group_by(across(all_of(adj_groups))) %>%
        mutate(
          Y_MEDIAN = median(2^(.data$Y_PREDICTED), na.rm = TRUE),
          Y_ADJ = 2^(.data$Y_PREDICTED)
        )
    }
  } else {
    if(!adjusted){
      d_smooth <- d_smooth %>%
        group_by(across(all_of(adj_groups))) %>%
        mutate(
          Y_MEDIAN = median(.data$Y_PREDICTED),
          Y_ADJ = .data$y / .data$Y_PREDICTED * median(.data$Y_PREDICTED, na.rm = TRUE)
        )
      } else {
        d_smooth <- d_smooth %>%
        group_by(across(all_of(adj_groups))) %>%
        mutate(
          Y_MEDIAN = median(.data$Y_PREDICTED, na.rm = TRUE),
          Y_ADJ = .data$Y_PREDICTED
        )
        }
    }


  # Calculate CVs and apply to all or conditionally
  if (within_batch) filter_groups <- c("feature_id", "batch_id") else filter_groups <- c("feature_id")
  d_smooth_final <- d_smooth %>%
    # filter(!is.na(.data$Y_PREDICTED)) |>
    group_by(across(all_of(filter_groups))) %>%
    mutate(
      CV_RAW_SPL = sd(.data$conc_raw[.data$qc_type == "SPL"], na.rm = TRUE) / mean(.data$conc_raw[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
      CV_ADJ_SPL = sd(.data$Y_ADJ[.data$qc_type == "SPL"], na.rm = TRUE) / mean(.data$Y_ADJ[.data$qc_type == "SPL"], na.rm = TRUE) * 100
    ) %>%
    mutate(
      DRIFT_CORRECTED = (((.data$CV_ADJ_SPL / .data$CV_RAW_SPL) < max_cv_ratio_before_after) | !apply_conditionally) & !is.na(.data$Y_ADJ) & !.data$fit_error,
      Y_FINAL = dplyr::if_else(.data$DRIFT_CORRECTED, .data$Y_ADJ, dplyr::if_else(is.na(.data$Y_ADJ) & !use_uncorrected_if_fail, NA_real_, .data$conc_raw))
    ) |>
    ungroup()

  data@dataset <- data@dataset %>%
    # filter(!outlier_technical) |>
    dplyr::left_join(
      d_smooth_final %>% dplyr::select(
        "analysis_id", "feature_id",
        CURVE_Y_PREDICTED = "Y_PREDICTED",
        Y_MEDIAN = "Y_MEDIAN",
        CONC_DRIFT_ADJ = "Y_ADJ",
        "CV_RAW_SPL", "CV_ADJ_SPL", "DRIFT_CORRECTED",
        FIT_ERROR = "fit_error",
        CONC_ADJ = "Y_FINAL"
      ),
      by = c("analysis_id", "feature_id")
    )


  # ToDo: below if some species have invalid numbers (e.g. negative) they will be ignore, but how about if eg 90% are invalide, Still correct result?
  d_sum_adj <- data@dataset |>
    filter(!is.na(.data$CURVE_Y_PREDICTED)) |>
    group_by(across(all_of(filter_groups))) |>
    summarise(DRIFT_CORRECTED = any(.data$DRIFT_CORRECTED & !is.na(.data$conc_raw)) | all(.data$FIT_ERROR)) |>
    group_by(feature_id) |>
    summarise(DRIFT_CORRECTED = any(.data$DRIFT_CORRECTED)) |>
    ungroup()

  d_fit <- data@dataset |>
    # filter(!outlier_technical) |>
    group_by(across(all_of(filter_groups))) %>%
    summarise(
      FIT_ERROR = any(.data$FIT_ERROR, na.rm = TRUE),
      CV_RAW_SPL = sd(.data$conc_raw[.data$qc_type == "SPL"], na.rm = TRUE) / mean(.data$conc_raw[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
      CV_ADJ_SPL = sd(.data$CONC_ADJ[.data$qc_type == "SPL"], na.rm = TRUE) / mean(.data$CONC_ADJ[.data$qc_type == "SPL"], na.rm = TRUE) * 100
    ) |>
    group_by(feature_id) |>
    summarise(
      FIT_ERROR = any(.data$FIT_ERROR, na.rm = TRUE),
      CV_RAW_SPL = median(.data$CV_RAW_SPL, na.rm = TRUE),
      CV_ADJ_SPL = median(.data$CV_ADJ_SPL, na.rm = TRUE)
    )


  fit_errors <- sum(d_fit$FIT_ERROR)
  cv_median_raw <- round(mean(d_fit$CV_RAW_SPL, na.rm = TRUE), 3)
  cv_median_adj <- round(mean(d_fit$CV_ADJ_SPL, na.rm = TRUE), 3)
  features_with_fiterror <- d_fit |>
    filter(.data$FIT_ERROR) |>
    pull(.data$feature_id)
    features_with_fiterror_text <- glue::glue_collapse(features_with_fiterror, ", ", last = " and ", width = 160)

  if (apply_conditionally_per_batch & within_batch) {
    d_sum_adj <- d_sum_adj |>
      summarise(DRIFT_CORRECTED = any(.data$DRIFT_CORRECTED))
  }
  if (!apply_conditionally) {
    count_feature_text <- glue::glue("{sum(d_sum_adj$DRIFT_CORRECTED, na.rm = TRUE)} of {nrow(d_sum_adj)} features")
  } else if (within_batch) {
    count_feature_text <- glue::glue("of at least one batch for {sum(d_sum_adj$DRIFT_CORRECTED, na.rm = TRUE)} of {nrow(d_sum_adj)} features")
  } else {
    count_feature_text <- glue::glue("{sum(d_sum_adj$DRIFT_CORRECTED, na.rm = TRUE)} of {nrow(d_sum_adj)} features.=")
  }
  if (within_batch) mode_text <- "(batch-wise)" else mode_text <- "(across all batches) "

  cli_alert_success(col_green(c("Drift correction {mode_text} of raw concentrations was applied to {count_feature_text}.
                                         Average CV of all features across samples before/after correction: {.strong {cv_median_raw}%} and {.strong {cv_median_adj}%}.")))



  d_sum_span <- d_smooth |>
    group_by(across(all_of(c("analysis_id", "batch_id", "qc_type")))) |>
    summarise(WITHIN_QC_SPAN = any(!is.na(.data$Y_PREDICTED)))

  n_spl_excl <- sum(d_sum_span$qc_type == "SPL" & !d_sum_span$WITHIN_QC_SPAN)
  n_nist_excl <- sum(d_sum_span$qc_type == "NIST" & !d_sum_span$WITHIN_QC_SPAN)
  n_bqc_excl <- sum(d_sum_span$qc_type == "BQC" & !d_sum_span$WITHIN_QC_SPAN)
  n_tqc_excl <- sum(d_sum_span$qc_type == "TQC" & !d_sum_span$WITHIN_QC_SPAN)
  n_ltr_excl <- sum(d_sum_span$qc_type == "LTR" & !d_sum_span$WITHIN_QC_SPAN)

  txt_1 <- txt_2 <- txt_3 <- txt_4 <- txt_5 <- character()
  if (n_spl_excl > 0) txt_1 <- paste0(n_spl_excl, " of ", sum(d_sum_span$qc_type == "SPL"), " study samples (SPL)")
  if (n_nist_excl > 0) txt_2 <- paste0(n_nist_excl, " of ", sum(d_sum_span$qc_type == "NIST"), " NISTs")
  if (n_bqc_excl > 0) txt_3 <- paste0(n_bqc_excl, " of ", sum(d_sum_span$qc_type == "BQC"), " BQCs")
  if (n_tqc_excl > 0) txt_4 <- paste0(n_tqc_excl, " of ", sum(d_sum_span$qc_type == "TQC"), " TQCs")
  if (n_ltr_excl > 0) txt_5 <- paste0(n_ltr_excl, " of ", sum(d_sum_span$qc_type == "LTR"), " LTRs")

  txt_final <- paste(c(txt_1, txt_2, txt_3, txt_4, txt_5), collapse = ", ")
  if (txt_final != "") cli_alert_warning(col_yellow(glue::glue("{txt_final} excluded from correction (beyond regions spanned by QCs).")))

  data@status_processing <- "Adjusted Quantitated Data"
  data@is_drift_corrected <- TRUE
  data@is_batch_corrected <- FALSE



  if (data@is_batch_corrected) cli_alert_info(col_blue(glue::glue("Previous batch correction has been removed.")))
  if (fit_errors > 0) cli_alert_warning(col_yellow(glue::glue("No smoothing applied for {fit_errors} features due failure of the fit (insufficient/invalid data): {features_with_fiterror_text}")))
  data@dataset$feature_conc <- data@dataset$CONC_ADJ
  data
}

#' Drift Correction by Gaussian Kernel Smoothing
#' @description
#' Function to correct for run-order drifts within or across batches using gaussian kernel smoothing (see *Tan et al. (2020)*).
#' This is typically used to smooth based on the study samples. To avoid local biases and artefacts, this function should only be applied to analyses wit sufficient number of samples that were well randomized.
#' @param data MidarExperiment object
#' @param qc_types QC types used for drift correction. Typically includes the study samples (`SPL`).
#' @param kernel_size Kernel bandwidth
#' @param outlier_filter Kernel Outlier filter
#' @param outlier_ksd Kernel K times standard deviation of data distribution
#' @param location_smooth Location parameter smoothing
#' @param scale_smooth Scale parameter smoothing
#' @param batch_wise Apply to each batch separately if `TRUE` (the default)
#' @param log2_transform Log transform the data for correction when `TRUE` (the default). Note: log transformation is solely applied internally for smoothing, results will not be be log-transformed. Log transformation may result in more robust smoothing that is less sensitive to outlier.
#' @param apply_conditionally Apply drift correction to all species if `TRUE`, or only when sample CV after smoothing changes below a threshold defined via `max_cv_ratio_before_after`
#' @param apply_conditionally_per_batch When `apply_conditionally = TRUE`, correction is conditionally applied per batch when `TRUE` and across all batches when `FALSE`
#' @param max_cv_ratio_before_after Only used when `apply_conditionally = TRUE`. Maximum allowed ratio of sample CV change before and after smoothing for the correction to be applied.
#' A value of 1 (the default) indicates the CV needs to improve or remain unchanged after smoothing so that the conditional smoothing is applied. A value of < 1 means that CV needs to improve, a value of e.g. 1.20 that the CV need to improve or get worse by max 1.20-fold after smoothing.
#' @param feature_list Subset the features for correction whose names matches the specified text using regular expression. Default is `NULL` which means all features are selected.
#' @param use_uncorrected_if_fail In case the smoothing function fails for a species, then use original (uncorrected) data when `TRUE` (the default) or return `NA` for all analyses of the feature where the fit failed.
#' @return MidarExperiment object
#' @references
#' Teo G., Chew WS, Burla B, Herr D, Tai ES, Wenk MR, Torta F, & Choi H (2020). MRMkit: Automated Data Processing for Large-Scale Targeted Metabolomics Analysis. *Analytical Chemistry*, 92(20), 13677–13682. \url{https://doi.org/10.1021/acs.analchem.0c03060}

#' @export
corr_drift_gaussiankernel <- function(data,
                                      qc_types,
                                      batch_wise = TRUE,
                                      kernel_size,
                                      outlier_filter = FALSE,
                                      outlier_ksd = 5,
                                      location_smooth = TRUE,
                                      scale_smooth = TRUE,
                                      log2_transform = TRUE,
                                      apply_conditionally = FALSE,
                                      apply_conditionally_per_batch = TRUE,
                                      feature_list = NULL,
                                      max_cv_ratio_before_after = 1,
                                      use_uncorrected_if_fail = FALSE
) {

  corr_drift_fun(
    data = data,
    smooth_fun = "fun_gauss.kernel.smooth",
    qc_types = qc_types,
    within_batch = batch_wise ,
    apply_conditionally = apply_conditionally,
    apply_conditionally_per_batch = apply_conditionally_per_batch,
    log2_transform = log2_transform,
    max_cv_ratio_before_after = max_cv_ratio_before_after,
    use_uncorrected_if_fail = use_uncorrected_if_fail,
    feature_list = feature_list,
    outlier_filter = outlier_filter,
    outlier_ksd = outlier_ksd,
    location_smooth = location_smooth,
    scale_smooth = location_smooth,
    kernel_size = kernel_size
  )
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
#' @param within_batch Apply to each batch separately if `TRUE` (the default)
#' @param log2_transform Log transform the data for correction when `TRUE` (the default). Note: log transformation is solely applied internally for smoothing, results will not be be log-transformed. Log transformation may result in more robust smoothing that is less sensitive to outlier.
#' @param apply_conditionally Apply drift correction to all species if `FALSE`, or only when sample CV after smoothing changes below a threshold defined via `max_cv_ratio_before_after`
#' @param apply_conditionally_per_batch When `apply_conditionally = TRUE`, correction is conditionally applied per batch when `TRUE` and across all batches when `FALSE`
#' @param max_cv_ratio_before_after Only used when `apply_conditionally = TRUE`. Maximum allowed ratio of sample CV change before and after smoothing for the correction to be applied.
#' A value of 1 (the default) indicates the CV needs to improve or remain unchanged after smoothing so that the conditional smoothing is applied. A value of < 1 means that CV needs to improve, a value of e.g. 1.20 that the CV need to improve or get worse by max 1.20-fold after smoothing.
#' @param feature_list Subset the features for correction whose names matches the specified text using regular expression. Default is `NULL` which means all features are selected.
#' @param use_uncorrected_if_fail In case the smoothing function fails for a species, then use original (uncorrected) data when `TRUE` (the default) or return `NA` for all analyses of the feature where the fit failed.
#' @param extrapolate Extrapolate loess smoothing. WARNING: It is generally not recommended to extrapolate outside of the range spanned by the QCs used for smoothing. See details below.
#' @return MidarExperiment object
#' @export
corr_drift_loess <- function(data, qc_types, within_batch = TRUE, apply_conditionally = FALSE, apply_conditionally_per_batch = TRUE,
                             log2_transform = TRUE, span = 0.75, feature_list = NULL, max_cv_ratio_before_after = 1, use_uncorrected_if_fail = TRUE, extrapolate = FALSE) {
  corr_drift_fun(
    data = data, smooth_fun = "fun_loess", qc_types = qc_types, within_batch = within_batch, apply_conditionally = apply_conditionally, apply_conditionally_per_batch = apply_conditionally_per_batch,
    log2_transform = log2_transform, feature_list = feature_list, max_cv_ratio_before_after = max_cv_ratio_before_after, use_uncorrected_if_fail = use_uncorrected_if_fail, extrapolate = extrapolate,  span_width = span
  )
}

#' Batch centering
#' @description
#' A short description...
#' #' @details
#' Additional details...
#' @param data MidarExperiment object
#' @param qc_types QC types used for batch correction
#' @param use_raw_concs Apply to unadjusted (raw) feature_conc. Default is FALSE, which means previously drift-corrected concs will be used if available, otherwise unadjusted concs will be used
#' @param center_fun Function used to center. Default is "median".
#' @return MidarExperiment object
#' @importFrom glue glue
#' @export
corr_batch_centering <- function(data, qc_types, use_raw_concs = FALSE, center_fun = "median") {
  ds <- data@dataset
  if (!data@is_drift_corrected | use_raw_concs) var <- rlang::sym("conc_raw") else var <- rlang::sym("CONC_ADJ")
  # Normalize by the median (or user-defined function)
  ds <- ds %>%
    dplyr::group_by(.data$feature_id, .data$batch_id) %>%
    dplyr::mutate(CONC_ADJ_NEW = {{ var }} / do.call(center_fun, list(({{ var }}[.data$qc_type %in% qc_types]), na.rm = TRUE))) |>
    dplyr::ungroup()

  # Re-level data to the median of all batches
  ds <- ds %>%
    dplyr::group_by(.data$feature_id) %>%
    dplyr::mutate(CONC_ADJ = .data$CONC_ADJ_NEW * do.call(center_fun, list({{ var }}[.data$qc_type %in% qc_types], na.rm = TRUE))) |>
    dplyr::select(-"CONC_ADJ_NEW") |>
    dplyr::ungroup()

  data@dataset <- ds
  if (data@is_drift_corrected) {
    cli_alert_success(col_green(glue::glue("Batch correction was applied to drift-corrected concs of all {nrow(data@annot_features)} features.")))
  } else {
    cli_alert_success(col_green(glue::glue("Batch correction was applied to raw concs of all {nrow(data@annot_features)} features.")))
  }

  if (data@is_drift_corrected & use_raw_concs) cli_alert_warning(col_yellow(glue::glue("Note: previous drift correction has been removed.")))
  data@status_processing <- "Adjusted Quantitated Data"
  data@is_batch_corrected <- TRUE
  data@dataset$feature_conc <- data@dataset$CONC_ADJ
  data
}
