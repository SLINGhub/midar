#' Gaussian Kernel smoothing helper function
#' @description
#' Function for Gaussian kernel-based smoothing, for use by `corr_drift_fun`.
#' @author Hyung Won Choi
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param reference_qc_types QC types used for the smoothing (fit) by loess
#' @param ... Additional parameters
#' @return List with a `data.frame` containing original `x` and the smoothed `y` values, and a `boolean` value indicting whether the fit failed or not not.

### Function to clean up serial trend within a batch
fun_gauss.kernel.smooth = function(tbl, reference_qc_types, ...) {
  arg <- list(...)

  xx <- seq(1, length(tbl$x))
  yy <- tbl$y
  s.train <- tbl$qc_type %in% reference_qc_types # Use sample for training

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

      # TODO: Also check HW. When to remove NA or warn there are too many NAs?
      yy.train[is.infinite(yy.train)] <- NA
      yy.train[is.nan(yy.train)] <- NA

      yy.est <- yy
      yy.corr <- yy #Added BB
      ## Location parameter smoothing

      if(arg$location_smooth) {
        for(i in 1:n) {
          wt <- (xx - xx[i]) / arg$kernel_size
          wt <- dnorm(wt, 0, 1)
          wt[is.na(yy.train)] <- NA
          yy.corr[i] <- median(yy, na.rm=TRUE) + sum(wt * yy.train, na.rm=TRUE) / sum(wt, na.rm=TRUE)  #Added by BB to get the fit
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
      list(y_fit = yy.corr, y_predicted = yy.final)
    } ,
    error = function(e) {
      print(e$message)
      return(rep(NA_real_, length(n)))
    }
  )

  ## Report the final values
  list(analysis_id = tbl$analysis_id, feature_id = tbl$feature_id, batch_id = tbl$batch_id, qc_type = tbl$qc_type, x = tbl$x, y_fit = res$y_fit , y_predicted = res$y_predicted, fit_error = all(is.na(res)), data_adjusted = TRUE)
}

#' Gaussian Kernel smoothing helper function
#' @description
#' Function for Gaussian kernel-based smoothing, for use by `corr_drift_fun`
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param reference_qc_types QC types used for the smoothing (fit) by loess
#' @param ... Additional parameters forwarded to KernSmooth::locpoly
#' @return List with a `data.frame` containing original `x` and the smoothed `y` values, and a `boolean` value indicting whether the fit failed or not not.
#' @noRd
fun_gaussiankernel_old <- function(tbl, reference_qc_types, ...) {
  arg <- list(...)

  if (!requireNamespace("KernSmooth", quietly = TRUE)) {
    stop(
      "Package {KernSmooth} must be installed when using this function. It is available from CRAN via `install.packages()`",
      call. = FALSE
    )
  }

  # browser()
  d_subset <- tbl[tbl$qc_type %in% reference_qc_types, ] |> tidyr::drop_na(.data$y)
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
#' @param reference_qc_types QC types used for the smoothing (fit) by loess
#' @param ... Additional parameters forwarded to Loess
#' @return List with a `data.frame` containing original x and the smoothed y values, and a `boolean` value indicting whether the fit failed or not not.
fun_loess <- function(tbl, reference_qc_types,...) {
  arg <- list(...)
  surface <- ifelse(arg$extrapolate, "direct", "interpolate")

  #sample_weights <- as.numeric(tbl$qc_type %in% reference_qc_types) # only train on selected sample types

  res <- tryCatch(
    {
      res_fit <- stats::loess(y ~ x,
        span = arg$span_width, family = "symmetric", degree = 2, normalize = FALSE, iterations = 4, surface = surface, ...,
        data = tbl
      )
      res <- stats::predict(res_fit, tibble::tibble(x = seq(1, nrow(tbl)))) |>
        as.numeric()
      list(y_fit = res_fit$fit, y_predicted = res)
    },
    error = function(e) {
      # print(e$message) # will be shown for each feature/batch...
      return(rep(NA_real_, length(tbl$x)))
    }
  )
  list(analysis_id = tbl$analysis_id, feature_id = tbl$feature_id, batch_id = tbl$batch_id, qc_type = tbl$qc_type, y_fit = res$y_fit, y_predicted = res$y_predicted, fit_error = all(is.na(res)), data_adjusted = TRUE)

}


#' Drift Correction by Custom Function
#' @description
#' Function to correct for run-order drifts within or across batches via a provided custom function
#' #' @details
#' The drift correction function needs to be provided by the user. See `smooth_fun` for details.
#' @param data MidarExperiment object
#' @param smooth_fun Function that performs drift correction. Function need to have following parameter `data` (`MidarExperiment`), `reference_qc_types` (one or more strings), and `span_width` (numerical).
#' Function needs to return a numerical vector with the length of number of rows in `data`. In case functions fails a vector with NA_real_ needs be returned
#' @param variable  The variable to be corrected for drift effects. Must be one of "intensity", "norm_intensity", or "conc"
#' @param reference_qc_types QC types used for drift correction
#' @param within_batch Apply to each batch separately if `TRUE` (the default)
#' @param replace_previous Logical. Replace previous correction (`TRUE`), or adds on top of previous correction (`FALSE`). Default is `TRUE`.
#' @param log_transform_internal Apply log transformation internally for smoothing if `TRUE` (default). This enhances robustness against outliers but does not affect the final data, which remains untransformed.
#' @param conditional_correction Determines whether drift correction should be applied to all species unconditionally (`TRUE`) or only when the sample CV after smoothing is below the ratio threshold specified by `cv_ratio_threshold`.
#' @param cv_ratio_threshold This parameter defines the maximum allowable change (ratio) in the coefficient of variation (CV) ratio of samples before and after smoothing for the correction to be applied.
#' A value of 1 (the default) requires the CV to improve, while a value above 1 allows the CV to improve, or to degrade by a maximum of the indicated fold-difference times post-smoothing.
#' @param ignore_istd Do not apply corrections to ISTDs
#' @param feature_list Subset the features for correction whose names matches the specified text using regular expression. Default is `NULL` which means all features are selected.
#' @param use_original_if_fail Determines the action when smoothing fails for a feature. If TRUE (default), the original data is used; if FALSE, the result for each analysis is NA.
#' @param recalc_trend_after Recalculate trend post-drift correction for `plot_qc_runscatter()`. This will double calculation time.
#' @param show_progress Show progress bar. Default = `TRUE.
#' @param ... Parameters specific for the smoothing function
#' @return MidarExperiment object
#' @export
corr_drift_fun <- function(data = NULL,
                           smooth_fun,
                           variable,
                           reference_qc_types,
                           within_batch,
                           replace_previous = TRUE,
                           log_transform_internal = TRUE,
                           conditional_correction = FALSE,
                           cv_ratio_threshold = 1,
                           use_original_if_fail = TRUE,
                           ignore_istd = TRUE,
                           feature_list = NULL,
                           recalc_trend_after = FALSE,
                           show_progress = TRUE,
                           ...) {


  check_data(data)

  variable_strip <- str_remove(variable, "feature_")
  rlang::arg_match(variable_strip, c("intensity", "norm_intensity", "conc"))
  variable <- stringr::str_c("feature_", variable_strip)
  variable_sym <- rlang::sym(variable)

  check_var_in_dataset(data@dataset, variable)

  variable_raw <- paste0(variable, "_raw")

   # Clear all previous calculations
  data@dataset <- data@dataset |> select(-any_of(c("curve_y_predicted", "y_fit", "y_fit_after", "y_predicted_median",
                                                   "cv_raw_spl", "cv_adj_spl", "drift_correct",
                                                   "fit_error", "var_adj")))

  # start from raw data, previous drift correction will be overwritten
  # if no drift correction has been applied yet, make a copy of the original (raw) data of the specified variable

  txt1 <- ifelse(data@var_drift_corrected[[variable]], "drift", NA)
  txt2 <- ifelse(data@var_batch_corrected[[variable]], "batch", NA)
  txt <- stringr::str_flatten(c(txt1, txt2), collapse = " and ", na.rm = TRUE)

  if(data@var_drift_corrected[[variable]]){
    if(!replace_previous){
      message(col_yellow(glue::glue("Adding correction on top of previous `{variable_strip}` {txt} corrections...")))
    } else { # make a copy of the original data
      message(col_yellow(glue::glue("Replacing previous `{variable_strip}` {txt} corrections...")))
      data@dataset[[variable]] <- data@dataset[[variable_raw]]
    }
  } else {
    data@dataset[[variable_raw]] <- data@dataset[[variable]]
    message(cli::col_green(glue::glue("Applying `{variable_strip}` drift correction...")))
  }

  # if variable has been batch_corrected, then it will also be overwritten
  if (data@var_drift_corrected[[variable]]){
    cli_alert_info(col_yellow(glue::glue("Previous drift and batch corrections are being overwritten!")))
    data@var_drift_corrected <- c(feature_intensity = FALSE, feature_norm_intensity = FALSE, feature_conc = FALSE)
    data@dataset[[variable]] <- data@dataset[[variable_raw]]
  }

    # Subset features
  ds <- data@dataset |> select("analysis_id", "qc_type", "batch_id", "feature_id", "is_istd", y_original = variable)


  if (!is.null(feature_list))
    ds <- ds |> dplyr::filter(stringr::str_detect(.data$feature_id, feature_list))

  if(ignore_istd) ds <- ds |> filter(!.data$is_istd)
  ds <- ds |> mutate(x = dplyr::row_number(), .by = "feature_id")
  ds$y <- ds$y_original
  if (log_transform_internal) (ds$y <- log2(ds$y))   #TODO: suppressWarnings?

  if (within_batch) adj_groups <- c("feature_id", "batch_id") else adj_groups <- c("feature_id")


  # call smooth function for each feature and each batch
  d_smooth_res <- ds |>
    select("analysis_id", "qc_type", "feature_id", "batch_id", "y_original", "x", "y") |>
    group_split(pick(adj_groups))

  # Calculate the total number of groups and set up the progress bar
  total_groups <- length(d_smooth_res)
  update_frequency <- ceiling(total_groups * 0.02)

  # Setup progress bar, as the build-in  based on cli does not work in quarto
  update_progress <- function(i) {
    if (i %% update_frequency == 0) {setTxtProgressBar(pb, i)}
  }

  if(show_progress) pb <- txtProgressBar(min = 0, max = total_groups, style = 3, width = 44)

  d_smooth_res <- d_smooth_res |>
    purrr::map2(seq_along(d_smooth_res), \(x, i) {
      if(show_progress) update_progress(i)
      do.call(smooth_fun, list(x, reference_qc_types, ...))
    }, .progress = FALSE) |>
    bind_rows()

  if(show_progress) {
    setTxtProgressBar(pb, total_groups)
    cat(cli::col_green(" - trend smoothing done!"))
    close(pb)
    }

  if(recalc_trend_after){
    if(show_progress) pb <- txtProgressBar(min = 0, max = total_groups, style = 3, width = 44)
     d_smooth_recalc <- d_smooth_res |>
      rename(y = .data$y_predicted) |>
      select("analysis_id", "qc_type", "feature_id", "batch_id", "x", "y") |>
      group_split(pick(adj_groups))
     d_smooth_recalc <- d_smooth_recalc |>
     purrr::map2(seq_along(d_smooth_recalc), \(x, i) {
       if(show_progress) update_progress(i)
       do.call(smooth_fun, list(x, reference_qc_types, ...))
     }, .progress = FALSE) |>
     bind_rows() |>
      select("analysis_id", "qc_type", "feature_id", "batch_id", y_fit_after = "y_fit")

      if(show_progress) {
        setTxtProgressBar(pb, total_groups)
        cat(cli::col_green(" - trend recalc done!"))
        close(pb)
      }

     d_smooth_res <- d_smooth_res |>
       left_join(d_smooth_recalc, by = c("analysis_id", "feature_id", "batch_id", "qc_type"))
  } else {
    d_smooth_res <- d_smooth_res |> mutate(y_fit_after = NA_real_)
  }
  # Get flag if data needs to be back-transformed (from smoothing function)
  is_adjusted <- all(d_smooth_res$data_adjusted) #TODOTODO

  # Backtransform data
  if (log_transform_internal)
    d_smooth_res <- d_smooth_res |> mutate(y_predicted = 2^(.data$y_predicted),
                                           y_fit = 2^(.data$y_fit),
                                           y_fit_after = 2^(.data$y_fit_after))

  if(is_adjusted){
    d_smooth_res <- d_smooth_res |>
      left_join(ds, by = c("analysis_id", "feature_id", "batch_id", "qc_type")) |>
      group_by(dplyr::pick(adj_groups)) |>
      mutate(
        y_predicted_median = median(.data$y_predicted, na.rm = TRUE),
        y_fit = .data$y_fit,
        y_fit_after = .data$y_fit_after,
        y_adj = .data$y_predicted)
    } else {
      d_smooth_res <- d_smooth_res |>
        group_by(dplyr::pick(adj_groups)) |>
        mutate(
          y_predicted_median = median(.data$y_predicted, na.rm = TRUE),
          y_fit = .data$y_fit, # TODODO: CHECK if correct,
          y_fit_after = .data$y_fit_after, # TODODO: CHECK if correct,
          y_adj = .data$y_original / .data$y_predicted * .data$y_predicted_median) # TODODO: CHECK if correct, log vs lin division
    }

  # Summarize which species to apply drift correction to and assign final concentrations
  d_smooth_summary <- d_smooth_res |>
      group_by(dplyr::pick(adj_groups)) |>
    summarise(
      any_fit_error = any(.data$fit_error, na.rm = TRUE),
      cv_raw_spl = cv(.data$y_original[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
      cv_adj_spl = cv(.data$y_adj[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
      cv_change = .data$cv_adj_spl - .data$cv_raw_spl,
      drift_correct =
        !conditional_correction | (.data$cv_adj_spl / .data$cv_raw_spl) < cv_ratio_threshold)|>
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
        .data$fit_error & use_original_if_fail ~ .data$y_original,
        .data$drift_correct & !.data$fit_error  ~ .data$y_adj,
        TRUE ~ NA_real_)
    ) |>
    ungroup() |>
    select("analysis_id", "feature_id", "fit_error", "y_fit",  "y_fit_after", "drift_correct", var_adj = "y_final")

  # Add drift-corrected data to the dataset

  variable_fit_sym <- rlang::sym(paste0(variable_raw, "_fit"))
  variable_fit_after_sym <- rlang::sym(paste0(variable, "_fit_after"))
  variable_corrected <- rlang::sym(paste0(variable, "_drift_correct"))
  variable_fit_error <- rlang::sym(paste0(variable, "_fit_error"))

  data@dataset <- data@dataset |>
    left_join(d_smooth_final, by = c("analysis_id", "feature_id")) |>
    mutate(!!variable_sym := .data$var_adj,
           !!variable_fit_sym := .data$y_fit,
           !!variable_fit_after_sym := .data$y_fit_after,
           !!variable_corrected := .data$drift_correct,
           !!variable_fit_error := .data$fit_error)

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


  nfeat <- get_feature_count(data, isistd = FALSE)

  if(conditional_correction & within_batch)
    count_feature_text <- glue::glue("of at least one batch for {features_corrected}
                                     of {nfeat} features")
  else
    count_feature_text <- glue::glue("{features_corrected} of {nfeat} features")

  mode_text <- ifelse(within_batch, "(batch-wise)", "(across all batches)")
  mode_text2 <- ifelse(within_batch, "in study samples (batch medians)", "in study samples (across batches)")

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


  if (features_with_fit_errors > 0) cli_alert_warning(col_yellow(glue::glue("No smoothing applied for {features_with_fit_errors} feature(s) due failure(s) of the fitting (insufficient/invalid data): {features_with_fit_errors_text}")))

  # Invalidate downstream processed data
  if(variable == "feature_intensity"){
    data <- change_is_normalized(data, FALSE)
  } else if(variable == "feature_norm_intensity") {
    data <- change_is_quantitated(data, FALSE)
  }

  data@status_processing <- "Drift-corrected concentrations"
  data@var_drift_corrected[[variable]] <- TRUE
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]
  data
}

#' Drift Correction by Gaussian Kernel Smoothing
#' @description
#' Function to correct for run-order drifts within or across batches using gaussian kernel smoothing (see *Tan et al. (2020)*).
#' This is typically used to smooth based on the study samples. To avoid local biases and artefacts, this function should only be applied to analyses wit sufficient number of samples that were well randomized.
#' @param data MidarExperiment object
#' @param variable  The variable to be corrected for drift effects. Must be one of "intensity", "norm_intensity", or "conc"
#' @param reference_qc_types QC types used for drift correction. Typically includes the study samples (`SPL`).
#' @param within_batch Apply to each batch separately if `TRUE` (the default), or across all batches if `FALSE`.
#' @param replace_previous Logical. Replace previous correction (`TRUE`), or adds on top of previous correction (`FALSE`). Default is `TRUE`.
#' @param kernel_size Kernel bandwidth
#' @param outlier_filter Kernel Outlier filter
#' @param outlier_ksd Kernel K times standard deviation of data distribution
#' @param location_smooth Location parameter smoothing
#' @param scale_smooth Scale parameter smoothing
#' @param log_transform_internal Apply log transformation internally for smoothing if `TRUE` (default). This enhances robustness against outliers but does not affect the final data, which remains untransformed.
#' @param recalc_trend_after Recalculate trends after smoothing, used for plotting (e.g., in `plot_runscatter()`)
#' @param ignore_istd Do not apply corrections to ISTDs
#' @param conditional_correction Apply drift correction to all species if `TRUE`, or only when sample CV after smoothing changes below a threshold defined via `cv_ratio_threshold`
#' @param cv_ratio_threshold Only used when `conditional_correction = TRUE`. Maximum allowed ratio of sample CV change before and after smoothing for the correction to be applied.
#' A value of 1 (the default) indicates the CV needs to improve or remain unchanged after smoothing so that the conditional smoothing is applied. A value of < 1 means that CV needs to improve, a value of e.g. 1.20 that the CV need to improve or get worse by max 1.20-fold after smoothing.
#' @param feature_list Subset the features for correction whose names matches the specified text using regular expression. Default is `NULL` which means all features are selected.
#' @param use_original_if_fail Determines the action when smoothing fails for a feature. If TRUE (default), the original data is used; if FALSE, the result for each analysis is NA.
#' @param show_progress Show progress bars. Set this to `FALSE` when rendering the notebook.
#' @return MidarExperiment object
#' @references
#' Teo G., Chew WS, Burla B, Herr D, Tai ES, Wenk MR, Torta F, & Choi H (2020). MRMkit: Automated Data Processing for Large-Scale Targeted Metabolomics Analysis. *Analytical Chemistry*, 92(20), 13677–13682. \url{https://doi.org/10.1021/acs.analchem.0c03060}

#' @export
correct_drift_gaussiankernel <- function(data = NULL,
                                      variable,
                                      reference_qc_types,
                                      within_batch = TRUE,
                                      replace_previous = TRUE,
                                      kernel_size,
                                      outlier_filter = FALSE,
                                      outlier_ksd = 5,
                                      location_smooth = TRUE,
                                      scale_smooth = FALSE,
                                      log_transform_internal = TRUE,
                                      recalc_trend_after = TRUE,
                                      conditional_correction = FALSE,
                                      feature_list = NULL,
                                      ignore_istd = TRUE,
                                      cv_ratio_threshold = 1,
                                      use_original_if_fail = FALSE,
                                      show_progress = TRUE
) {

  check_data(data)

  corr_drift_fun(
    data = data,
    smooth_fun = "fun_gauss.kernel.smooth",
    variable = variable,
    reference_qc_types = reference_qc_types,
    within_batch = within_batch ,
    replace_previous = replace_previous,
    conditional_correction = conditional_correction,
    log_transform_internal = log_transform_internal,
    cv_ratio_threshold = cv_ratio_threshold,
    use_original_if_fail = use_original_if_fail,
    feature_list = feature_list,
    recalc_trend_after = recalc_trend_after,
    ignore_istd = ignore_istd,
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
#' @param reference_qc_types QC types used for drift correction
#' @param variable  The variable to be corrected for drift effects. Must be one of "intensity", "norm_intensity", or "conc"
#' @param span Loess span width (default is 0.75)
#' @param within_batch Apply to each batch separately if `TRUE` (the default), or across all batches if `FALSE`.
#' @param replace_previous Logical. Replace previous correction (`TRUE`), or adds on top of previous correction (`FALSE`). Default is `TRUE`.
#' @param log_transform_internal Log transform the data for correction when `TRUE` (the default). Note: log transformation is solely applied internally for smoothing, results will not be be log-transformed. Log transformation may result in more robust smoothing that is less sensitive to outlier.
#' @param conditional_correction Determines whether drift correction should be applied to all species unconditionally (`TRUE`) or only when the sample CV after smoothing is below the ratio threshold specified by `cv_ratio_threshold`.
#' @param cv_ratio_threshold This parameter defines the maximum allowable change (ratio) in the coefficient of variation (CV) ratio of samples before and after smoothing for the correction to be applied.
#' A value of 1 (the default) requires the CV to improve, while a value above 1 allows the CV to improve, or to degrade by a maximum of the indicated fold-difference times post-smoothing.
#' @param feature_list Subset the features for correction whose names matches the specified text using regular expression. Default is `NULL` which means all features are selected.
#' @param ignore_istd Do not apply corrections to ISTDs
#' @param use_original_if_fail Determines the action when smoothing fails for a feature. If TRUE (default), the original data is used; if FALSE, the result for each analysis is NA.
#' @param extrapolate Extrapolate loess smoothing. WARNING: It is generally not recommended to extrapolate outside of the range spanned by the QCs used for smoothing. See details below.
#' @return MidarExperiment object
#' @export
correct_drift_loess <- function(data = NULL,
                                variable,
                                reference_qc_types,
                                within_batch = TRUE,
                                replace_previous = TRUE,
                                span = 0.75,
                                conditional_correction = FALSE,
                                log_transform_internal = TRUE,
                                feature_list = NULL,
                                ignore_istd = TRUE,
                                cv_ratio_threshold = 1,
                                use_original_if_fail = TRUE,
                                extrapolate = FALSE) {
  check_data(data)
  corr_drift_fun(
    data = data,
    smooth_fun = "fun_loess",
    reference_qc_types = reference_qc_types,
    variable = variable,
    within_batch = within_batch,
    replace_previous = TRUE,
    conditional_correction = conditional_correction,
    log_transform_internal = log_transform_internal,
    feature_list = feature_list,
    ignore_istd = ignore_istd,
    cv_ratio_threshold = cv_ratio_threshold,
    use_original_if_fail = use_original_if_fail,
    extrapolate = extrapolate,
    span_width = span
  )
}




#' Batch Centering Correction
#'
#' This function performs batch centering correction of each feature. Optionally,
#' the scale of batches can be equalized across batches.The
#' selected QC types (reference_qc_types) are used to calculate the
#' medians and align all other qc types based on them. The correction can be
#' applied to one of three variables: "intensity", "norm_intensity", or "conc".
#' The correction can be applied on top of previous
#' corrections or  replace all previous batch corrections.
#'
#' @param data A MidarExperiment object containing the data to be corrected.
#'   This object must include information about QC types and measurements.
#' @param variable The variable to be corrected. Must be one of "intensity",
#'   "norm_intensity", or "conc".
#' @param reference_qc_types A character vector specifying the QC types to be
#'   used as reference for batch centering.
#' @param correct_location A logical value indicating whether to align the
#'   median locations (centers) of the batches. Defaults to TRUE.
#' @param correct_scale A logical value indicating whether to scale the batches
#'   to the same level (equalize scale). Defaults to FALSE.
#' @param replace_previous A logical value indicating whether to replace any
#'   previous batch correction or apply the new correction on top. Defaults to
#'   TRUE (replace).
#' @param log_transform_internal A logical value indicating whether to
#'   log-transform the data internally during correction. Defaults to TRUE.
#' @param ... Additional parameters that can be passed to the batch correction
#'   function.
#' @return A MidarExperiment object containing the corrected data.
#' @seealso plot_runscatter for visualizing the correction before and after.
#' @export
correct_batch_centering <- function(data = NULL,
                                    variable,
                                    reference_qc_types,
                                    correct_location = TRUE,
                                    correct_scale = FALSE,
                                    replace_previous = TRUE,
                                    log_transform_internal = TRUE,
                                    ...) {
  check_data(data)

  variable_strip <- str_remove(variable, "feature_")
  rlang::arg_match(variable_strip, c("intensity", "norm_intensity", "conc"))
  variable <- stringr::str_c("feature_", variable_strip)
  variable_sym <- rlang::sym(variable)
  variable_raw <- paste0(variable, "_raw")
  variable_fit <- paste0(variable_raw, "_fit")
  variable_fit_sym <- rlang::sym(variable_fit)
  variable_fit_after <- paste0(variable, "_fit_after")
  variable_fit_after_sym <- rlang::sym(variable_fit_after)
  variable_smoothed_fit_after <- paste0(variable, "_smoothed_fit_after")
  variable_smoothed<- paste0(variable, "_smoothed")

  check_var_in_dataset(data@dataset, variable)

  # start from raw data, previous drift correction will be overwritten
  # if no drift correction has been applied yet, make a copy of the original (raw) data of the specified variable

  txt1 <- ifelse(data@var_drift_corrected[[variable]], "drift", NA)
  txt2 <- ifelse(data@var_drift_corrected[[variable]], "batch", NA)
  txt <- stringr::str_flatten(c(txt1, txt2), collapse = " and ", na.rm = TRUE)

  if(data@var_drift_corrected[[variable]] ){
    # Drift-corrected data
    if (data@var_batch_corrected[[variable]] ){
      # Drift but AND batch-corrected data
      if(replace_previous){
        message(col_yellow(glue::glue("Replacing previous `{variable_strip}` batch correction of drift-corrected data...")))
        # use drift corrected data and apply correction on top, replacing any previous batch correction
        data@dataset[[variable]] <- data@dataset[[variable_smoothed]]
        data@dataset[[variable_fit_after]] <- data@dataset[[variable_smoothed_fit_after]]
      } else {
        message(col_yellow(glue::glue("Adding batch correction on top of previous `{variable_strip}` drift and batch corrections...")))
        # use main variable (e.g. conc) and add correction on top
        #data@dataset[[variable]] <- data@dataset[[variable]]
      }
    } else {
      # Drift but NOT batch-corrected data
      # replace_previous will have no effect as no previous batch correction has beem applied


      data@dataset[[variable_smoothed]] <- data@dataset[[variable]]
      data@dataset[[variable_smoothed_fit_after]] <- data@dataset[[variable_fit_after]]

      message(col_yellow(glue::glue("Adding batch correction on top of `{variable_strip}` drift-correction...")))


    }
  } else {
    # Data is NOT drift corrected
    if (data@var_batch_corrected[[variable]]){
        # Only batch-corrected data (no previous drift correction)
        if(replace_previous){
            message(col_yellow(glue::glue("Replacing previous `{variable_strip}` batch correction...")))
            # use drift corrected data and apply correction on top, replacing any previous batch correction
            data@dataset[[variable]] <- data@dataset[[variable_raw]]
            data@dataset[[variable_fit_after]] <- data@dataset[[variable_fit]]

        } else {
            message(col_yellow(glue::glue("Adding batch correction on top of previous `{variable_strip}` batch correction...")))
            # use main variable (e.g. conc) and add correction on top
            #data@dataset[[variable]] <- data@dataset[[variable]]
        }
    } else {
        # Uncorrected data (raw)
        message(col_yellow(glue::glue("Adding batch correction to `{variable_strip}` data...")))
        #Store raw data
        data@dataset[[variable_raw]] <- data@dataset[[variable]]
        # Use median horizontal lines to indicate fits if data was not drift corrected before
        data@dataset <- data@dataset |>
            mutate(!!variable_fit_sym := median(.data[[variable_sym]][.data$qc_type == "SPL"], na.rm = TRUE), .by = c("feature_id", "batch_id"))
        data@dataset[[variable_fit_after]] <- data@dataset[[variable_fit]]
      }

    }

  ds <- data@dataset |>
    select(any_of(c("analysis_id", "feature_id", "qc_type", "batch_id", y_fit_after = variable_fit_after, y = variable)))
  nbatches <- length(unique(ds$batch_id))
  if(nbatches < 2) {
    cli_abort(col_yellow(glue::glue("Batch correction was not applied as there is only one batch.")))
    return(data)
  }

  d_res <- ds |>
    group_by(.data$feature_id) |>
      nest() |>
      mutate(
        res = purrr::map(data, \(x) do.call("batch.correction", list(x, reference_qc_types, correct_location, correct_scale, ...))),
      ) |>
      unnest(cols = c(.data$res)) |>
      select(-"data")

  d_res_sum <- d_res |>
    group_by(.data$feature_id) |>
    summarise(
      cv_before = cv(.data$y_fit_after[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
      cv_after = cv(.data$y[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
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

  nfeat <- get_feature_count(data, isistd = FALSE)

  # Print summary
  if (data@var_drift_corrected[[variable]]) {
    cli_alert_success(col_green(glue::glue("Batch median-centering of {nbatches} batches was applied to drift-corrected concentrations of all {nfeat} features.")))
    data@status_processing <- "Batch- and drift-corrected concentrations"
  } else {
    cli_alert_success(col_green(glue::glue("Batch median-centering of {nbatches} batches was applied to raw concentrations of all {nfeat} features.")))
    data@status_processing <- "Batch-corrected concentrations"
  }
  # Print stats
  text_change <- ifelse(d_res_sum$cv_diff_median > 0, "increased", "decreased")

  cli_alert_info(cli::col_grey(
      c("The median CV of features in the study samples across batches {.strong {text_change}} by {d_res_sum$cv_diff_text}% ({d_res_sum$cv_diff_q1} to {d_res_sum$cv_diff_q3}%) to {format(round(d_res_sum$cv_after,1), nsmall = 1)}%.")))

  # Return data

  data@dataset <- data@dataset |>
    left_join(d_res |> select(-y, -"y_fit_after"),
              by = c("analysis_id", "feature_id", "qc_type", "batch_id")) |>
    mutate(!!variable_sym := .data$y_adjusted,
           !!variable_fit_sym := .data[[variable_fit]],
           !!variable_fit_after_sym := .data$y_fit_after_adjusted) |>
    select(-"y_adjusted", -"y_fit_after_adjusted")


  data@var_batch_corrected[[variable]] <- TRUE
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]
  data
}



batch.correction = function(tab,
                            reference_qc_types,
                            correct_location,
                            correct_scale,
                            log_transform_internal = TRUE, ...) {

  batch <- tab$batch_id
  batch.order <- seq(1, nrow(tab))
  val <- tab$y
  y_fit_after <- tab$y_fit_after  # BB
  if(log_transform_internal){
    val <- log2(val)
    y_fit_after <- log2(y_fit_after)
  }
  sample.for.loc <- tab$qc_type %in% reference_qc_types # Use sample for location median

  ubatch <- unique(batch)
  nbatch <- length(ubatch)

  val.clean <- val ## placeholder
  y_fit_after.clean <- y_fit_after # BB



  ### Cross-batch scale normalization
  if(correct_location) {
    tmp <- val.clean
    tmp_fit_after <- y_fit_after.clean # BB
    tmp_for_loc <-val.clean
    tmp_for_loc[!sample.for.loc] <- NA_real_
    loc.batch <- rep(NA, nbatch)
    for(b in 1:nbatch) {
      id <- which(batch == ubatch[b])
      loc.batch[b] <- median(tmp_for_loc[id], na.rm=TRUE)
    }
    loc.batch.mean = mean(loc.batch, na.rm=TRUE)
    for(b in 1:nbatch) {
      id <- which(batch == ubatch[b])
      xloc <- loc.batch[b]
      if(log_transform_internal){
        val.clean[id] <- (tmp[id] - xloc) + loc.batch.mean
        y_fit_after.clean[id] <- (tmp_fit_after[id] - xloc) + loc.batch.mean
      } else {
        val.clean[id] <- (tmp[id] / xloc) * loc.batch.mean
        y_fit_after.clean[id] <- (tmp_fit_after[id] / xloc) * loc.batch.mean
      }
    }
  }

  if(correct_scale) {
    tmp <- val.clean
    tmp_fit_after <- y_fit_after.clean # BB
    tmp_for_loc <- val.clean
    tmp_for_loc[!sample.for.loc] <- NA_real_

    # TODO: confirm if this is ok to do and its consequences. Once batch is NA, then shoud we report all?
    tmp_for_loc[is.infinite(tmp_for_loc)] <- NA_real_
    tmp_for_loc[is.nan(tmp_for_loc)] <- NA_real_


    loc.batch <- rep(NA, nbatch)
    sca.batch <- rep(NA, nbatch)
    for(b in 1:nbatch) {
      id = which(batch == ubatch[b])
      loc.batch[b] <- median(tmp_for_loc[id], na.rm=TRUE)
      sca.batch[b] <- mad(tmp_for_loc[id], na.rm=TRUE)
    }
    loc.batch.mean <- mean(loc.batch, na.rm=TRUE)
    sca.batch.mean <- mean(sca.batch, na.rm=TRUE)
    for(b in 1:nbatch) {
      id <- which(batch == ubatch[b])
      xloc <- loc.batch[b]
      if(log_transform_internal) {
        val.clean[id] <- (tmp[id] - xloc) / sca.batch[b] * sca.batch.mean + loc.batch.mean
        y_fit_after.clean[id] <- (tmp_fit_after[id] - xloc) / sca.batch[b] * sca.batch.mean + loc.batch.mean
      } else {
        stop("Non-log batch scaling not implemented yet")
        val.clean[id] <- (tmp[id] - xloc) / sca.batch[b] * sca.batch.mean + loc.batch.mean
      }
    }
  }

    tab$y_adjusted <- if(log_transform_internal) 2^val.clean else val.clean
    tab$y_fit_after_adjusted <- if(log_transform_internal) 2^y_fit_after.clean else  y_fit_after.clean
  tab

}

