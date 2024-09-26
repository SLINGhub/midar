#' Gaussian Kernel smoothing helper function
#' @description
#' Function for Gaussian kernel-based smoothing, for use by `corr_drift_fun`.
#' @author Hyung Won Choi
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param qc_types QC types used for the smoothing (fit) by loess
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
            wt[is.na(yy.train)] <- NA
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
#' @param ... Additional parameters forwarded to KernSmooth::locpoly
#' @return List with a `data.frame` containing original `x` and the smoothed `y` values, and a `boolean` value indicting whether the fit failed or not not.

fun_gaussiankernel_old <- function(tbl, qc_types, ...) {
  arg <- list(...)

  if (!requireNamespace("KernSmooth", quietly = TRUE)) {
    stop(
      "Package {KernSmooth} must be installed when using this function. It is available from CRAN via `install.packages()`",
      call. = FALSE
    )
  }

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
#' @param ... Additional parameters forwarded to Loess
#' @return List with a `data.frame` containing original x and the smoothed y values, and a `boolean` value indicting whether the fit failed or not not.
fun_loess <- function(tbl, qc_types,...) {
  #tbl <- tbl |> filter(!.data$outlier_technical)
  arg <- list(...)
  surface <- ifelse(arg$extrapolate, "direct", "interpolate")

  #sample_weights <- as.numeric(tbl$qc_type %in% qc_types) # only train on selected sample types

  res <- tryCatch(
    {
      stats::loess(y ~ x,
        span = arg$span_width, family = "symmetric", degree = 2, normalize = FALSE, iterations = 4, surface = surface, ...,
        data = tbl
      ) |>
        stats::predict(tibble::tibble(x = seq(1, nrow(tbl)))) |>
        as.numeric()
    },
    error = function(e) {
      # print(e$message) # will be shown for each feature/batch...
      return(rep(NA_real_, length(tbl$x)))
    }
  )
  list(res = res, fit_error = all(is.na(res)), data_adjusted = FALSE)
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
#' @param ignore_istd Do not apply corrections to ISTDs
#' @param feature_list Subset the features for correction whose names matches the specified text using regular expression. Default is `NULL` which means all features are selected.
#' @param use_uncorrected_if_fail In case the smoothing function fails for a species, then use original (uncorrected) data when `TRUE` (the default) or return `NA` for all analyses of the feature where the fit failed.
#' @param ... Parameters specific for the smoothing function
#' @return MidarExperiment object
#' @export
corr_drift_fun <- function(data, smooth_fun, qc_types, log2_transform = TRUE, within_batch, apply_conditionally, apply_conditionally_per_batch = TRUE,
                           max_cv_ratio_before_after = 1, use_uncorrected_if_fail = TRUE, ignore_istd = TRUE, feature_list = NULL, ...) {

  # Clear all previous calculations
  data@dataset <- data@dataset |> select(-any_of(c("CURVE_y_predicted", "y_predicted_median",
                                                   "CONC_DRIFT_ADJ", "cv_raw_spl", "cv_adj_spl", "drift_correct",
                                                   "fit_error", "feature_conc_adj")))

  if(!"conc_raw" %in% names(data@dataset))
    cli_abort(col_red(glue::glue("Drift correction currently only implemented for concentration values. Please processes first with `quantitate_by_istd`")))

    # Subset features
  ds <- data@dataset |> select("analysis_id", "qc_type", "batch_id", "feature_id", "is_istd", y_original = "conc_raw")


  if (!is.null(feature_list))
    ds <- ds |> dplyr::filter(stringr::str_detect(.data$feature_id, feature_list))

  if(ignore_istd) ds <- ds |> filter(!.data$is_istd)
  ds <- ds |> mutate(x = dplyr::row_number(), .by = "feature_id")
  ds$y <- ds$y_original
  if (log2_transform) (ds$y <- log2(ds$y))   #TODO: suppressWarnings?

  if (within_batch) adj_groups <- c("feature_id", "batch_id") else adj_groups <- c("feature_id")

  d_smooth_res <- ds |>
      select("analysis_id", "qc_type", "feature_id", "batch_id", "y_original", "x", "y") |>
      #filter(.data$qc_type %in% c("SPL", "BQC", "TQC", "NIST", "LTR", "PBLK", "SBPK", "RQC")) |>
      group_by(across(all_of(adj_groups))) |>
      nest() |>
      mutate(
        res = purrr::map(data, \(x) do.call(smooth_fun, list(x, qc_types, ...))),
        y_predicted = purrr::map(.data$res, \(x) x$res),
        fit_error = purrr::map(.data$res, \(x) x$fit_error),
        data_adjusted = purrr::map(.data$res, \(x) x$data_adjusted)
      ) |>
      unnest(cols = c(data, .data$y_predicted, .data$fit_error,.data$data_adjusted))

  # Get flag if data needs to be back-transformed (from smoothing function)
  is_adjusted <- all(d_smooth_res$data_adjusted) #TODOTODO

  # Backtransform data
  if (log2_transform)
    d_smooth_res <- d_smooth_res |> mutate(y_predicted = 2^(.data$y_predicted))

  if(is_adjusted){
    d_smooth_res <- d_smooth_res |>
      left_join(ds, by = c("analysis_id", "feature_id", "batch_id", "qc_type", "y_original")) |>
      group_by(across(all_of(adj_groups))) |>
      mutate(
        y_predicted_median = median(.data$y_predicted, na.rm = TRUE),
        y_adj = .data$y_predicted)
    } else {
      d_smooth_res <- d_smooth_res |>
        group_by(across(all_of(adj_groups))) |>
        mutate(
          y_predicted_median = median(.data$y_predicted, na.rm = TRUE),
          y_adj = .data$y_original / .data$y_predicted * .data$y_predicted_median)
    }
  # Summarize which species to apply drift correction to and assign final concentrations

  d_smooth_summary <- d_smooth_res |>
    group_by(across(all_of(adj_groups))) |>
    summarise(
      any_fit_error = any(.data$fit_error, na.rm = TRUE),
      cv_raw_spl = cv(.data$y_original[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
      cv_adj_spl = cv(.data$y_adj[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
      cv_change = .data$cv_adj_spl - .data$cv_raw_spl,
      drift_correct =
        !apply_conditionally | (.data$cv_adj_spl / .data$cv_raw_spl) < max_cv_ratio_before_after)|>
    group_by(.data$feature_id) |>
    summarise(
      any_fit_error = any(.data$any_fit_error, na.rm = TRUE),
      drift_correct = any(.data$drift_correct, na.rm = TRUE),
      cv_raw_spl = median(.data$cv_raw_spl, na.rm = TRUE),
      cv_adj_spl = median(.data$cv_adj_spl, na.rm = TRUE),
      cv_change = median(.data$cv_change, na.rm = TRUE)) |>
      ungroup()



  # Replace values with drift-corrected values
  ## TODO check what is needed in terms of check to assign value
  d_smooth_final <- d_smooth_res |>
    dplyr::left_join(d_smooth_summary, by = c("feature_id"))  |>
    group_by(.data$feature_id) |>
    mutate(
      y_final = case_when(
        is.na(.data$y_adj) ~ NA_real_,
        .data$fit_error & use_uncorrected_if_fail ~ .data$y_original,
        .data$drift_correct & !.data$fit_error  ~ .data$y_adj,
        TRUE ~ NA_real_)
    ) |>
    ungroup() |>
    select("analysis_id", "feature_id", "fit_error", "drift_correct", feature_conc_adj = "y_final")

  # Add drift-corrected data to the datast
  data@dataset <- data@dataset |>
    left_join(d_smooth_final, by = c("analysis_id", "feature_id")) |>
    mutate(feature_conc = .data$feature_conc_adj)


  # Prepare info/texts for command line output
  features_with_fit_errors <- sum(d_smooth_summary$any_fit_error)
  #features_corrected <- sum(d_smooth_summary$drift_correct) - features_with_fit_errors

  features_corrected <- data@dataset |>
    filter(!.data$is_istd, .data$drift_correct) |>
    pull(.data$feature_id) |>
    unique() |>
    length()

  # Calculate median CVs for print summary
  cv_median_raw <- round(mean(d_smooth_summary$cv_raw_spl, na.rm = TRUE), 1)
  cv_median_adj <- round(mean(d_smooth_summary$cv_adj_spl, na.rm = TRUE), 1)
  cv_difference_median <- round(median(d_smooth_summary$cv_change, na.rm = TRUE), 1)
  cv_difference_q1 <- format(round(cv_difference_median - IQR(d_smooth_summary$cv_change, na.rm = TRUE), 1), nsmall = 1)
  cv_difference_q3 <- format(round(cv_difference_median + IQR(d_smooth_summary$cv_change, na.rm = TRUE), 1), nsmall = 1)

  features_with_fit_errors_text <- glue::glue_collapse(
    d_smooth_summary$feature_id[d_smooth_summary$any_fit_error], ", ",
    last = " and ", width = 160)


  nfeat <- get_feature_count(data, istd = FALSE)

  if(apply_conditionally & within_batch)
    count_feature_text <- glue::glue("of at least one batch for {features_corrected}
                                     of {nfeat} features")
  else
    count_feature_text <- glue::glue("{features_corrected} of {nfeat} features")

  mode_text <- ifelse(within_batch, "(batch-wise)", "(across all batches)")
  mode_text2 <- ifelse(within_batch, "in study samples (median of batches)", "in study samples (across batches)")

  text_change <- ifelse(cv_difference_median > 0, "increased", "decreased")

  cli_alert_success(
    col_green(
      c("Drift correction {mode_text} was applied to raw concentrations of {count_feature_text}.")))

  cli_alert_info(cli::col_grey(
           "The median CV of all features {mode_text2} {.strong {text_change}} by {cv_difference_median}% ({cv_difference_q1} to {cv_difference_q3}%) to {format(cv_median_adj, nsmall = 1)}%."))



  d_sum_span <- d_smooth_res |>
    group_by(across(all_of(c("analysis_id", "batch_id", "qc_type")))) |>
    summarise(WITHIN_QC_SPAN = any(!is.na(.data$y_predicted)))

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



  if (data@is_batch_corrected) cli_alert_info(col_yellow(glue::glue("Previous batch correction has been removed.")))
  if (features_with_fit_errors > 0) cli_alert_warning(col_yellow(glue::glue("No smoothing applied for {features_with_fit_errors} feature(s) due failure(s) of the fitting (insufficient/invalid data): {features_with_fit_errors_text}")))

  data@status_processing <- "Drift-corrected concentrations"
  data@is_drift_corrected <- TRUE
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
                                      scale_smooth = FALSE,
                                      log2_transform = TRUE,
                                      apply_conditionally = FALSE,
                                      apply_conditionally_per_batch = FALSE,
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
    scale_smooth = scale_smooth,
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
corr_drift_loess <- function(data, qc_types, within_batch = TRUE, span = 0.75, apply_conditionally = FALSE, apply_conditionally_per_batch = TRUE,
                             log2_transform = TRUE,  feature_list = NULL, max_cv_ratio_before_after = 1, use_uncorrected_if_fail = TRUE, extrapolate = FALSE) {
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
corr_batch_centering <- function(data, qc_types, use_raw_concs = FALSE, center_fun = "median") {
  ds <- data@dataset
  if (!data@is_drift_corrected | use_raw_concs) var <- rlang::sym("conc_raw") else var <- rlang::sym("CONC_ADJ")
  # Normalize by the median (or user-defined function)
  ds <- ds |>
    dplyr::group_by(.data$feature_id, .data$batch_id) |>
    dplyr::mutate(CONC_ADJ_NEW = {{ var }} / do.call(center_fun, list(({{ var }}[.data$qc_type %in% qc_types]), na.rm = TRUE))) |>
    dplyr::ungroup()

  # Re-level data to the median of all batches
  ds <- ds |>
    dplyr::group_by(.data$feature_id) |>
    dplyr::mutate(CONC_ADJ = .data$CONC_ADJ_NEW * do.call(center_fun, list({{ var }}[.data$qc_type %in% qc_types], na.rm = TRUE))) |>
    dplyr::select(-"CONC_ADJ_NEW") |>
    dplyr::ungroup()

  data@dataset <- ds
  if (data@is_drift_corrected) {
    cli_alert_success(col_green(glue::glue("Batch correction was applied to drift-corrected concs of all {get_feature_count(data)} features.")))
  } else {
    cli_alert_success(col_green(glue::glue("Batch correction was applied to raw concs of all {get_feature_count(data)} features.")))
  }

  if (data@is_drift_corrected & use_raw_concs) cli_alert_warning(col_yellow(glue::glue("Note: previous drift correction has been removed.")))
  data@status_processing <- "Batch-adjusted concentrations"
  data@is_batch_corrected <- TRUE
  data@dataset$feature_conc <- data@dataset$CONC_ADJ
  data
}


#' Batch centering
#' @description
#' A short description...
#' #' @details
#' Additional details...
#' @param data MidarExperiment object
#' @param qc_types QC types used for batch correction
#' @param correct_location Align locations (median) of batches
#' @param correct_scale Scale batches to the same level
#' @param overwrite_batch_corr Overwrite previous batch correction or apply on top existing batch correction
#' @param log_transform Log transform the data internally for correction. Will not transform the final results.
#' @param ... Other parameters for batch correction function. Currently not in use.
#'
#' @return MidarExperiment object
#' @export
corr_batcheffects <- function(data, qc_types, correct_location = TRUE, correct_scale = FALSE, overwrite_batch_corr = TRUE, log_transform = TRUE, ...) {
  ds <- data@dataset |> select("analysis_id", "feature_id", "qc_type", "batch_id", "feature_conc")
  nbatches <- length(unique(ds$batch_id))

  if(nbatches < 2) {
    cli_abort(col_yellow(glue::glue("Batch correction was not applied as there is only one batch.")))
    return(data)
  }
  if (data@is_batch_corrected){
      if(overwrite_batch_corr){
        cli_alert_warning(col_yellow(glue::glue("Previous batch-correction was overwritten.")))
        data@dataset$feature_conc <- data@dataset$feature_conc_adj_raw
  } else {
        cli_alert_warning(col_yellow(glue::glue("Batch correction was applied onto previous batch-correction!")))
      }
  } else {
      data@dataset$feature_conc_adj_raw <- data@dataset$feature_conc
  }
  # TODO var <- rlang::sym("conc_raw") else var <- rlang::sym("CONC_ADJ")
  #if (!data@is_drift_corrected)

  d_res <- ds |>
  group_by(.data$feature_id) |>
    nest() |>
    mutate(
      res = purrr::map(data, \(x) do.call("batch.correction", list(x, qc_types, correct_location, correct_scale, ...))),
    ) |>
    unnest(cols = c(.data$res)) |>
    select(-"data")

  d_res_sum <- d_res |>
    group_by(.data$feature_id) |>
    summarise(
      cv_before = cv(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
      cv_after = cv(.data$y_adjusted[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
      cv_diff = .data$cv_after - .data$cv_before,
    ) |> ungroup() |>
    summarise(
      cv_before = mean(.data$cv_before, na.rm = TRUE),
      cv_after = mean(.data$cv_after, na.rm = TRUE),
      cv_diff_median = median(.data$cv_diff, na.rm = TRUE),
      cv_diff_q1 = format(round(.data$cv_diff_median - IQR(.data$cv_diff, na.rm = TRUE), 1), nsmall = 1),
      cv_diff_q3 = format(round(.data$cv_diff_median + IQR(.data$cv_diff, na.rm = TRUE), 1), nsmall = 1),
      cv_diff_text = format(round(.data$cv_diff_median, 1), nsmall = 1)
    )

  nfeat <- get_feature_count(data, istd = FALSE)

  # Print summary
  if (data@is_drift_corrected) {
    cli_alert_success(col_green(glue::glue("Batch median-centering of {nbatches} batches was applied to drift-corrected concentrations of all {nfeat} features.")))
    data@status_processing <- "Batch- and drift-corrected concentrations"
  } else {
    cli_alert_success(col_green(glue::glue("Batch median-centering of {nbatches} batches was applied to raw concentrations of all {nfeat} features.")))
    data@status_processing <- "Batch-corrected concentrations"
  }
  # Print stats
  text_change <- ifelse(d_res_sum$cv_diff_median > 0, "increased", "decreased")

  cli_alert_info(cli::col_grey(
      c("The median CV of all features in study samples across batches {.strong {text_change}} by {d_res_sum$cv_diff_text}% ({d_res_sum$cv_diff_q1} to {d_res_sum$cv_diff_q3}%) to {format(round(d_res_sum$cv_after,1), nsmall = 1)}%.")))

  # Return data

  data@dataset <- data@dataset |>
    left_join(d_res |> select(-"feature_conc"),
              by = c("analysis_id", "feature_id", "qc_type", "batch_id")) |>
    mutate(feature_conc = .data$y_adjusted) |>
    select(-"y_adjusted")

  data@is_batch_corrected <- TRUE
  data
}



batch.correction = function(tab,
                            qc_types,
                            correct_location,
                            correct_scale,
                            log_transform = TRUE, ...) {

  batch <- tab$batch_id
  batch.order <- seq(1, nrow(tab))


  val <- tab$feature_conc
  if(log_transform) val <- log2(val)
  sample.for.loc <- tab$qc_type %in% qc_types # Use sample for location median

  ubatch <- unique(batch)
  nbatch <- length(ubatch)

  val.clean <- val ## placeholder

  ### Cross-batch scale normalization
  if(correct_location) {
    tmp <- val.clean
    tmp_for_loc <-val.clean
    tmp_for_loc[!sample.for.loc] <- NA_real_
    loc.batch <- rep(NA, nbatch)
    for(b in 1:nbatch) {
      id <- which(batch == ubatch[b])
      loc.batch[b] <- median(tmp_for_loc[id], na.rm=TRUE)
    }
    loc.batch.mean = mean(loc.batch)
    for(b in 1:nbatch) {
      id <- which(batch == ubatch[b])
      xloc <- loc.batch[b]
      if(log_transform)
        val.clean[id] <- (tmp[id] - xloc) + loc.batch.mean
      else
        val.clean[id] <- (tmp[id] / xloc) * loc.batch.mean
    }
  }

  if(correct_scale) {
    tmp <- val.clean
    tmp_for_loc <- val.clean
    tmp_for_loc[!sample.for.loc] <- NA_real_
    loc.batch <- rep(NA, nbatch)
    sca.batch <- rep(NA, nbatch)
    for(b in 1:nbatch) {
      id = which(batch == ubatch[b])
      loc.batch[b] <- median(tmp_for_loc[id], na.rm=TRUE)
      sca.batch[b] <- mad(tmp_for_loc[id], na.rm=TRUE)
    }
    loc.batch.mean <- mean(loc.batch)
    sca.batch.mean <- mean(sca.batch)
    for(b in 1:nbatch) {
      id <- which(batch == ubatch[b])
      xloc <- loc.batch[b]
      if(log_transform)
        val.clean[id] <- (tmp[id] - xloc) / sca.batch[b] * sca.batch.mean + loc.batch.mean
      else
        stop("Non-log batch scaling not implemented yet")
        val.clean[id] <- (tmp[id] - xloc) / sca.batch[b] * sca.batch.mean + loc.batch.mean
    }
  }

    tab$y_adjusted <- if(log_transform) 2^val.clean else val.clean
  tab

}
