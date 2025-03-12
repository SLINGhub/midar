#' Gaussian Kernel smoothing helper function
#' @description
#' Function for Gaussian kernel-based smoothing, for use by `fun_correct_drift`.
#' @author Hyung Won Choi
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param reference_qc_types QC types used for the smoothing (fit) by loess
#' @param log_transform_internal Apply log transformation internally for smoothing if `TRUE` (default). This does not affect the final data, which remains untransformed.
#' @param ... Additional arguments
#' @return List with a `data.frame` containing original `x` and the smoothed `y` values, and a `boolean` value indicting whether the fit failed or not not.

### Function to clean up serial trend within a batch
fun_gauss.kernel.smooth = function(tbl, reference_qc_types, log_transform_internal, ...) {
  arg <- list(...)
  xx <- seq(1, length(tbl$x))

  if (log_transform_internal) tbl$y <- log10(tbl$y)
  yy <- tbl$y
  s.train <- tbl$qc_type %in% reference_qc_types # Use sample for training

  n = length(yy) ## number of data points

  warnings_list <- NULL

  res <- tryCatch(
    {
      withCallingHandlers({

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

        if(log_transform_internal) {
          y_predicted <- 10^yy.final
          y_fit <- 10^yy.corr
        } else {
          y_predicted <-  yy.final
          y_fit = yy.corr
        }
      },
      warning = function(w) {
        warnings_list <<- c(warnings_list, paste(conditionMessage(w), collapse = "|"))
        invokeRestart("muffleWarning")
      }
    )
    list(y_fit = y_fit, y_predicted = y_predicted, has_error = FALSE)
  },
    error = function(e) {
      print(e$message)
      return(list(y_fit = NA_real_, y_adj = NA_real_, has_error = TRUE))
    }
  )

  ## Report the final values
  res$fit_warning <- !is.null(warnings_list)
  list(analysis_id = tbl$analysis_id, feature_id = tbl$feature_id, batch_id = tbl$batch_id, qc_type = tbl$qc_type, x = tbl$x, y_fit = res$y_fit , y_adj = res$y_predicted, fit_error = res$has_error, fit_warning = res$fit_warning)
}


#' Loess smoothing helper function
#' @description
#' Function for loess-based smoothing, for use by `fun_correct_drift`
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param reference_qc_types QC types used for the smoothing (fit) by loess
#' @param log_transform_internal Apply log transformation internally for smoothing if `TRUE` (default). This does not affect the final data, which remains untransformed.
#' @param ... Additional arguments forwarded to Loess
#' @return List with a `data.frame` containing original x and the smoothed y values, and a `boolean` value indicting whether the fit failed or not not.
fun_loess <- function(tbl, reference_qc_types, log_transform_internal, ...) {
  arg <- list(...)
  surface <- ifelse(arg$extrapolate, "direct", "interpolate")

  if (log_transform_internal) tbl$y <- log10(tbl$y)   #TODO: suppressWarnings?

  warnings_list <- NULL
  res <- tryCatch(
    {
      withCallingHandlers({
        tbl_train <- tbl[tbl$qc_type %in% reference_qc_types,]
        res_fit <- stats::loess(y ~ x,
          span = arg$span, family = "symmetric", degree = arg$degree, normalize = FALSE, iterations = 4, surface = surface,
          data = tbl_train
        )
        y_fit <- stats::predict(res_fit,
                                dplyr::tibble(x = seq(min(tbl$x), max(tbl$x),1))) |>
          as.numeric()

        if(log_transform_internal) {
          y_predicted <- tbl$y - y_fit + median(y_fit, na.rm = TRUE)
          y_predicted <- 10^y_predicted
          y_fit <- 10^y_fit
        } else {
          y_predicted <- tbl$y / y_fit * median(y_fit, na.rm = TRUE)
        }
      },
      warning = function(w) {

        warnings_list <<- c(warnings_list, conditionMessage(w))

        invokeRestart("muffleWarning")
      }

    )

    list(y_fit = y_fit, y_predicted = y_predicted, has_error = FALSE)
    },
    error = function(e) {
      # print(e$message) # will be shown for each feature/batch...

      return(list(y_fit = NA_real_, y_predicted = NA_real_, has_error = TRUE))
    }
  )

  res$fit_warning <- !is.null(warnings_list) && !res$has_error

  list(analysis_id = tbl$analysis_id, feature_id = tbl$feature_id, batch_id = tbl$batch_id, qc_type = tbl$qc_type, x = tbl$x, y_fit = res$y_fit, y_adj = res$y_predicted, fit_error = res$has_error, fit_warning = res$fit_warning)
}

#' Cubic spline smoothing helper function
#' @description
#' Function for cubic spline-based smoothing with optional cross-validation, for use by `fun_correct_drift`
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param reference_qc_types QC types used for the smoothing (fit) by cubic spline
#' @param log_transform_internal Apply log transformation internally for smoothing if `TRUE` (default). This does not affect the final data, which remains untransformed.
#' @param ... Additional arguments forwarded to `smooth.spline`
#' @return List with a `data.frame` containing original x and the smoothed y values, and a `boolean` value indicating whether the fit failed or not.
fun_cspline <- function(tbl, reference_qc_types, log_transform_internal, ...) {
  arg <- list(...)

  # Apply log transformation if required (used internally for better smoothing performance)
  if (log_transform_internal) tbl$y <- log10(tbl$y)

  warnings_list <- NULL
  res <- tryCatch(
    {
      withCallingHandlers({
        # Subset the data to include only reference QC types for training
        tbl_train <- tbl[tbl$qc_type %in% reference_qc_types,]
        # Fit cubic spline model, using cross-validation if specified

        if (is.null(arg$lambda))
          res_fit <- stats::smooth.spline(x = tbl_train$x, y = tbl_train$y, cv = arg$cv, spar = arg$spar, penalty = arg$penalty,  all.knots = FALSE)
        else
          res_fit <- stats::smooth.spline(x = tbl_train$x, y = tbl_train$y, cv = arg$cv, lambda = arg$lambda, penalty = arg$penalty,  all.knots = FALSE)

        # Predict smoothed values across the full range of x values
        y_fit <- predict(res_fit, x = seq(min(tbl$x), max(tbl$x), 1))$y

        # Adjust predictions based on log transformation if applied
        if (log_transform_internal) {
          y_predicted <- tbl$y - y_fit + median(y_fit, na.rm = TRUE)  # Center residuals
          y_predicted <- 10^y_predicted  # Convert back from log scale
          y_fit <- 10^y_fit  # Convert fitted values back from log scale
        } else {
          y_predicted <- tbl$y / y_fit * median(y_fit, na.rm = TRUE)  # Scale back predictions
        }
      },
      warning = function(w) {
        # Capture warnings without stopping execution
        warnings_list <<- c(warnings_list, conditionMessage(w))
        invokeRestart("muffleWarning")
      })

      list(y_fit = y_fit, y_predicted = y_predicted, has_error = FALSE)
    },
    error = function(e) {

      # Handle errors gracefully by returning NA values
      return(list(y_fit = NA_real_, y_predicted = NA_real_, has_error = TRUE))
    }
  )

  # Indicate whether warnings occurred during fitting
  res$fit_warning <- !is.null(warnings_list) && !res$has_error

  # Return final results as a structured list
  list(analysis_id = tbl$analysis_id, feature_id = tbl$feature_id, batch_id = tbl$batch_id, qc_type = tbl$qc_type, x = tbl$x, y_fit = res$y_fit, y_adj = res$y_predicted, fit_error = res$has_error, fit_warning = res$fit_warning)
}


#' Generalized Additive Model (GAM) smoothing helper function
#' @description
#' Function for penalized spline-based smoothing using GAM, for use by `fun_correct_drift`
#' @param tbl Table (`tibble` or `data.frame`) containing the fields `qc_type`, `x` (run order number), and `y` (variable)
#' @param reference_qc_types QC types used for the smoothing (fit) by GAM
#' @param log_transform_internal Apply log transformation internally for smoothing if `TRUE` (default). This does not affect the final data, which remains untransformed.
#' @param ... Additional arguments forwarded to `mgcv::gam`
#' @return List with a `data.frame` containing original x and the smoothed y values, and a `boolean` value indicating whether the fit failed or not.
fun_gam_smooth <- function(tbl, reference_qc_types, log_transform_internal = TRUE, ...) {

  arg <- list(...)

  # Apply log transformation if required (used internally for better smoothing performance)
  if (log_transform_internal) tbl$y <- log10(tbl$y)

  warnings_list <- NULL
  res <- tryCatch(
    {
      withCallingHandlers({
        # Subset the data to include only reference QC types for training
        tbl_train <- tbl[tbl$qc_type %in% reference_qc_types,]

        # Fit GAM with user-defined spline type, basis functions (k), and smoothing parameter (sp)
        res_fit <- mgcv::gam(y ~ s(x, bs = arg$bs, k = arg$k), data = tbl_train, method = "REML", sp = arg$sp)

        # Predict smoothed values across the full range of x values
        y_fit <- predict(res_fit, newdata = data.frame(x = tbl$x), type = "response")

        # Adjust predictions based on log transformation if applied
        if (log_transform_internal) {
          y_predicted <- tbl$y - y_fit + median(y_fit, na.rm = TRUE)  # Center residuals
          y_predicted <- 10^y_predicted  # Convert back from log scale
          y_fit <- 10^y_fit  # Convert fitted values back from log scale
        } else {
          y_predicted <- tbl$y / y_fit * median(y_fit, na.rm = TRUE)  # Scale back predictions
        }
      },
      warning = function(w) {
        # Capture warnings without stopping execution
        warnings_list <<- c(warnings_list, conditionMessage(w))
        invokeRestart("muffleWarning")
      })

      list(y_fit = y_fit, y_predicted = y_predicted, has_error = FALSE)
    },
    error = function(e) {
      # Handle error by returning NA values
      return(list(y_fit = NA_real_, y_predicted = NA_real_, has_error = TRUE))
    }
  )

  # Indicate whether warnings occurred during fitting
  res$fit_warning <- !is.null(warnings_list) && !res$has_error

  # Return final results as a structured list
  list(analysis_id = tbl$analysis_id, feature_id = tbl$feature_id, batch_id = tbl$batch_id,
       qc_type = tbl$qc_type, x = tbl$x, y_fit = res$y_fit, y_adj = res$y_predicted,
       fit_error = res$has_error, fit_warning = res$fit_warning)
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
#' @param batch_wise Apply to each batch separately if `TRUE` (the default)
#' @param replace_previous Logical. Replace previous correction (`TRUE`), or adds on top of previous correction (`FALSE`). Default is `TRUE`.
#' @param log_transform_internal Apply log transformation internally for smoothing if `TRUE` (default). This enhances robustness against outliers but does not affect the final data, which remains untransformed.
#' @param conditional_correction Determines whether drift correction should be applied to all features unconditionally (`TRUE`) or
#' only when the difference of sample CV before vs after smoothing is below the threshold specified by `cv_diff_threshold`.
#' @param cv_diff_threshold This parameter defines the maximum allowable change (difference) in the coefficient of variation (CV) ratio of samples before and after smoothing for the correction to be applied.
#' A value of 0 (the default) requires the CV to improve, while a value above 0 allows the CV to also become worse by a maximum of the defined difference.
#' @param ignore_istd Do not apply corrections to ISTDs
#' @param feature_list Sets specific features for correction only. Can be character vector or regular expression. Default is `NULL` which means all features are selected.
#' @param use_original_if_fail Determines the action when smoothing fails or results in invalid values for a feature. If TRUE (default), the original data is used; if FALSE, the result for each analysis is NA.
#' @param recalc_trend_after Recalculate trend post-drift correction for `plot_qc_runscatter()`. This will double calculation time.
#' @param show_progress Show progress bar. Default = `TRUE.
#' @param ... Arguments specific for the smoothing function
#' @return MidarExperiment object
#' @export
fun_correct_drift <- function(data = NULL,
                           smooth_fun,
                           variable,
                           reference_qc_types,
                           batch_wise,
                           replace_previous = TRUE,
                           log_transform_internal = TRUE,
                           conditional_correction = FALSE,
                           cv_diff_threshold = 0,
                           use_original_if_fail = TRUE,
                           ignore_istd = TRUE,
                           feature_list = NULL,
                           recalc_trend_after = FALSE,
                           show_progress = TRUE,
                           ...) {


  check_data(data)

  if (!all(reference_qc_types %in% unique(data@dataset$qc_type))) {
    cli::cli_abort(col_red("One or more specified `qc_types` are not present in the dataset. Please verify data or analysis metadata."))
  }

  variable_strip <- str_remove(variable, "feature_")
  rlang::arg_match(variable_strip, c("intensity", "norm_intensity", "conc"))
  variable <- stringr::str_c("feature_", variable_strip)
  variable_before <- stringr::str_c("feature_", variable_strip, "_before")
  variable_sym <- rlang::sym(variable)
  variable_before_sym <- rlang::sym(variable_before)

  check_var_in_dataset(data@dataset, variable)

  variable_raw <- paste0(variable, "_raw")
  variable_before <- paste0(variable, "_before")


   # Clear all previous calculations
  data@dataset <- data@dataset |> select(-any_of(c("curve_y_predicted", "y_fit", "y_fit_after", "y_predicted_median",
                                                   "cv_raw_spl", "cv_adj_spl", "drift_correct",
                                                   "fit_error", "fit_warning", "var_adj")))

  # start from raw data, previous drift correction will be overwritten
  # if no drift correction has been applied yet, make a copy of the original (raw) data of the specified variable

  txt1 <- ifelse(data@var_drift_corrected[[variable]], "drift", NA)
  txt2 <- ifelse(data@var_batch_corrected[[variable]], "batch", NA)
  txt <- stringr::str_flatten(c(txt1, txt2), collapse = " and ", na.rm = TRUE)

  is_first_correction <- FALSE

  if(data@var_drift_corrected[[variable]] | data@var_batch_corrected[[variable]]){
    if(!replace_previous){
      cli::cli_alert_info(col_yellow("Adding correction on top of previous '{.emph {variable_strip}}' {txt} corrections."))

    } else { # make a copy of the original data
      cli::cli_alert_info(col_yellow("Replacing previous `{variable_strip}` {txt} corrections..."))
      data@var_drift_corrected <- c(feature_intensity = FALSE, feature_norm_intensity = FALSE, feature_conc = FALSE)
      data@var_batch_corrected <- c(feature_intensity = FALSE, feature_norm_intensity = FALSE, feature_conc = FALSE)
      data@dataset[[variable]] <- data@dataset[[variable_raw]]
    }
  } else {
    # data was not corrected before, make a copy of the original data as "_raw"
    is_first_correction <- TRUE
    data@dataset[[variable_raw]] <- data@dataset[[variable]]
    cli::cli_alert_info(cli::col_green(glue::glue("Applying `{variable_strip}` drift correction...")))
  }


# Subset features
  ds <- data@dataset |> select("analysis_id", "qc_type", "batch_id", "feature_id", "is_istd", "y_original" = all_of(variable))


  if (!is.null(feature_list)){
    if (length(feature_list) == 1) {
      ds <- ds  |> dplyr::filter(stringr::str_detect(.data$feature_id, feature_list))
      if(nrow(ds) == 0){
        cli::cli_abort(col_red("The feature filter set via `feature_list` does not match any feature in the dataset."))
      }
    } else {

      if(!all(feature_list %in% unique(ds$feature_id))){
        cli::cli_abort(col_red("One or more feature(s) specified with `feature_list` are not present in the dataset."))
      }

      ds <- ds  |> dplyr::filter(.data$feature_id %in% feature_list)
    }
  }

  if(ignore_istd) ds <- ds |> filter(!.data$is_istd)

  ds <- ds |> mutate(x = dplyr::row_number(), .by = "feature_id")
  ds$y <- ds$y_original

  if(log_transform_internal) {
    count_negative_or_zero <- ds |>
      group_by(.data$feature_id) |>
      summarise(count = sum(.data$y <= 0)) |>
      filter(.data$count > 0)
    if(nrow(count_negative_or_zero) > 0) {
      ds$y[ds$y <= 0] <- NA_real_
      cli::cli_alert_info(
        col_yellow(
          "{nrow(count_negative_or_zero)} feature(s) contain one or more zero or negative `{variable_strip}` values. Verify your data or use `log_transform_internal = FALSE`."
        )
      )
    }
  }



  if (batch_wise) adj_groups <- c("feature_id", "batch_id") else adj_groups <- c("feature_id")


  # call smooth function for each feature and each batch
  d_smooth_res <- ds |>
    select("analysis_id", "qc_type", "feature_id", "batch_id", "y_original", "x", "y") |>
    group_split(!!!syms(adj_groups))

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
      do.call(smooth_fun, list(x, reference_qc_types, log_transform_internal, ...))
    }, .progress = FALSE) |>
    bind_rows()

  if(show_progress) {
    setTxtProgressBar(pb, total_groups)
    cat(cli::col_green(" - trend smoothing done!"))
    close(pb)
    }

  if(recalc_trend_after){
    if(show_progress) pb <- txtProgressBar(min = 0, max = total_groups, style = 3, width = 44)

    d_smooth_recalc <- d_smooth_res |> rename(y = "y_adj")

    d_smooth_recalc <- d_smooth_recalc |>
      select("analysis_id", "qc_type", "feature_id", "batch_id", "x", "y") |>
       #group_by(across(all_of(adj_groups))) |>
       group_split(!!!syms(adj_groups))

    d_smooth_recalc <- d_smooth_recalc |>
      purrr::map2(seq_along(d_smooth_recalc), \(x, i) {
        if(show_progress) update_progress(i)
        do.call(smooth_fun, list(x, reference_qc_types, log_transform_internal, ...))
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


  # add back the original data
  d_smooth_res <- d_smooth_res |>
    left_join(ds, by = c("analysis_id", "feature_id", "batch_id", "qc_type", "x")) |>
    ungroup()

  # Summarize which features to apply drift correction to and assign final concentrations
  # 1. First the change per feature per BATCH is calculated
  # 2. Then the median change per feature is calculated for all or only those that meet the criteria
  d_smooth_summary_bybatch <- d_smooth_res |>
    group_by(!!!syms(adj_groups)) |>
    summarise(
      #.by = !!!syms(adj_groups),
      any_fit_error = any(.data$fit_error, na.rm = TRUE),
      any_fit_warning = any(.data$fit_warning, na.rm = TRUE),
      cv_raw_spl = cv(.data$y_original[.data$qc_type == "SPL"], na.rm = TRUE),
      cv_adj_spl = cv(.data$y_adj[.data$qc_type == "SPL"], na.rm = TRUE),
      cv_change = if_else(.data$cv_raw_spl > 0 | .data$cv_adj_spl > 0, .data$cv_adj_spl - .data$cv_raw_spl, NA_real_),
      cv_change_valid =  !is.na(.data$cv_change) & (is.nan(.data$cv_adj_spl) | is.na(.data$cv_adj_spl) | is.infinite(.data$cv_adj_spl)  | .data$cv_adj_spl > 0),
      drift_correct =
        (!conditional_correction) | (.data$cv_change  < cv_diff_threshold)
      ) |> ungroup()

  d_smooth_summary <- d_smooth_summary_bybatch |>
    group_by(.data$feature_id) |>
    summarise(
      any_fit_error_summary = any(.data$any_fit_error, na.rm = TRUE),
      all_fit_error_summary = all(.data$any_fit_error, na.rm = TRUE),
      any_fit_warning_summary = any(.data$any_fit_warning, na.rm = TRUE),
      drift_correct_summary = any(.data$drift_correct, na.rm = TRUE),
      cv_raw_spl_mean = mean(.data$cv_raw_spl[.data$cv_raw_spl > 0], na.rm = TRUE),
      cv_adj_spl_mean = mean(.data$cv_adj_spl[.data$cv_adj_spl > 0], na.rm = TRUE),
      cv_change_mean = if_else(!all(is.na(.data$cv_change)), mean(.data$cv_change[.data$drift_correct], na.rm = TRUE), NA_real_),
      cv_change_valid = any(.data$cv_change_valid),
    )|>
      ungroup()


  invariant_features <- sum((d_smooth_summary$cv_raw_spl_mean == 0 | is.nan(d_smooth_summary$cv_raw_spl_mean)))


  if(invariant_features > 0){
    txt <- if(!ignore_istd && any(ds$is_istd)) "To ignore ISTDs set argument `ignore_istd = TRUE`." else ""
    cli_alert_warning(cli::col_yellow(
      "{invariant_features} features showed no variation in the study sample's original values across analyses. {txt}"))
  }

  if(sum(!d_smooth_summary$cv_change_valid > 0 )){

    txt <- if(use_original_if_fail) "The original values were kept for these features" else "NA will be be returned for all values of these faetures. Set `use_original_if_fail = FALSE to return orginal values."

    cli_alert_warning(cli::col_yellow(
      "{sum(!d_smooth_summary$cv_change_valid)} features have invalid values after smoothing. {txt}."))

  }

  # Prepare info/texts for command line output
  features_with_fit_errors <- sum(d_smooth_summary$any_fit_error_summary)
  features_with_fit_errors_allbatches <- sum(d_smooth_summary$all_fit_error_summary)
  features_with_fit_warnings <- sum(d_smooth_summary$any_fit_warning_summary)


  features_corrected <- d_smooth_summary|>
    filter(.data$drift_correct_summary, !is.na(.data$cv_change_mean)) |>
    nrow()

  if (features_with_fit_errors_allbatches > 0)
  {

    cli_abort(col_red("Smoothing failed for {features_with_fit_errors_allbatches} feature(s) in all batches. Please check data, metadata, and fit parameters."))
  }
  features_with_fit_errors_text <- glue::glue_collapse(
    d_smooth_summary$feature_id[d_smooth_summary$any_fit_error_summary], ", ",
    last = " and ", width = 80)

  features_with_fit_warnings_text <- glue::glue_collapse(
    d_smooth_summary$feature_id[d_smooth_summary$any_fit_warning_summary], ", ",
    last = " and ", width = 80)

  if (features_with_fit_errors > 0) cli_alert_warning(col_yellow("Smoothing failed for {features_with_fit_errors} feature(s) in at least one batch: {features_with_fit_errors_text}. Please check data, metadata and fit parameters."))
  if (features_with_fit_warnings > 0) {
    if(features_with_fit_warnings == features_corrected)
      cli_alert_warning(col_yellow("Issues (warnings) occured during smoothing of all features in at least one batch. Please inspect arguments and data, and consider plotting with `plot_runscatter()`."))
    else
      cli_alert_warning(col_yellow("Issues (warnings) occured during the smoothing of {features_with_fit_warnings} feature(s) in at least one batch: {features_with_fit_warnings_text}. Please inspect argument and data, and consider plotting with `plot_runscatter()`."))

  }


  # Replace values with drift-corrected values
  ## TODO check what is needed in terms of check to assign value
  d_smooth_final <- d_smooth_res |>
    dplyr::left_join(d_smooth_summary_bybatch |> select(-"cv_change_valid"), by = adj_groups)  |>
    dplyr::left_join(d_smooth_summary |> select("feature_id", "any_fit_error_summary", "cv_change_valid"), by = "feature_id")  |>
    group_by(!!!syms(adj_groups)) |>
    mutate(
      #.by = !!!syms(adj_groups),
      y_final = case_when(
        is.na(.data$y_adj) ~ NA_real_,
        !.data$cv_change_valid & !use_original_if_fail  ~ NA_real_,
        !.data$drift_correct | ((.data$any_fit_error_summary | !.data$cv_change_valid) & use_original_if_fail)  ~ .data$y_original,
        .data$drift_correct & !.data$any_fit_error_summary  ~ .data$y_adj,
        TRUE ~ NA_real_)
    ) |>
    ungroup() |>
    select("analysis_id", "qc_type", "batch_id", "feature_id", "fit_error", "fit_warning", "y_fit",  "y_fit_after", "drift_correct", "y_original", var_adj = "y_final")

  d_stats <- d_smooth_final |>
    group_by(!!!syms(adj_groups)) |>
    summarise(
    #  .by = !!!syms(adj_groups),
      cv_raw_spl = cv(.data$y_original[.data$qc_type == "SPL"], na.rm = TRUE),
      cv_adj_spl = cv(.data$var_adj[.data$qc_type == "SPL"], na.rm = TRUE)
    ) |>
    group_by(.data$feature_id) |>
    summarise(
      cv_raw_spl_mean = mean(.data$cv_raw_spl[.data$cv_raw_spl > 0], na.rm = TRUE),
      cv_adj_spl_mean = mean(.data$cv_adj_spl[.data$cv_adj_spl > 0], na.rm = TRUE)
    ) |> ungroup()


  # Add drift-corrected data to the dataset

  variable_fit_sym <- rlang::sym(paste0(variable_before, "_fit"))
  variable_fit_after_sym <- rlang::sym(paste0(variable, "_fit_after"))
  variable_raw_fit_sym <- rlang::sym(paste0(variable_raw, "_fit"))
  variable_corrected <- rlang::sym(paste0(variable, "_drift_correct"))
  variable_fit_error <- rlang::sym(paste0(variable, "_fit_error"))
  variable_fit_warning <- rlang::sym(paste0(variable, "_fit_warning"))

  data@dataset <- data@dataset |>
    select(-any_of(c("y_original", variable_before))) |>
    left_join(d_smooth_final, by = c("analysis_id", "qc_type", "feature_id", "batch_id")) |>
    mutate(!!variable_sym := .data$var_adj,
           !!variable_before_sym := .data$y_original,
           !!variable_raw_fit_sym := if(is_first_correction) .data$y_fit else !!variable_raw_fit_sym,
           !!variable_fit_sym := .data$y_fit,
           !!variable_fit_after_sym := .data$y_fit_after,
           !!variable_corrected := .data$drift_correct,
           !!variable_fit_error := .data$fit_error,
           !!variable_fit_warning := .data$fit_warning)


  # Calculate median CVs across features for print summary.
  # In case of conditional correction, only consider features that were corrected
  # for the median change
  cv_median_raw <- median(d_stats$cv_raw_spl_mean, na.rm = TRUE)
  cv_median_adj <- median(d_stats$cv_adj_spl_mean, na.rm = TRUE)
  cv_difference_median <- median(d_smooth_summary$cv_change_mean, na.rm = TRUE)
  cv_difference_low <- safe_min(d_smooth_summary$cv_change_mean, na.rm = TRUE)
  cv_difference_high <- safe_max(d_smooth_summary$cv_change_mean, na.rm = TRUE)

  d_sum_span <- d_smooth_res |>
    group_by(across(all_of(c("analysis_id", "batch_id", "qc_type")))) |>
    summarise(WITHIN_QC_SPAN = any(!is.na(.data$y_fit))) |>
    ungroup()

  if(any(d_sum_span$WITHIN_QC_SPAN)){

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
    if (txt_final != "") {
      txt <- glue::glue_collapse(reference_qc_types, sep = ", ", last = ", and ")
      cli_alert_warning(col_yellow("{txt_final} were excluded from correction as they fall outside the regions spanned by the QCs/samples used for smoothing ({.emph {txt}})."))
    }
  }

  nfeat <- get_feature_count(data, is_istd = ifelse(ignore_istd, FALSE, NA))

  if(conditional_correction & batch_wise)
    count_feature_text <- glue::glue("{features_corrected}
                                     of {nfeat} features in at least one batch per feature")
  else
    count_feature_text <- glue::glue("{features_corrected} of {nfeat} features")

  mode_text <- ifelse(batch_wise, "(batch-wise)", "(across all batches)")
  mode_text2 <- ifelse(batch_wise, "(batch medians)", "(across batches)")

  text_change <- case_when(
     round(cv_median_adj,2) - round(cv_median_raw ,2) >= 0.01 ~ "increased from",
     round(cv_median_adj,2) - round(cv_median_raw ,2) <= -0.01 ~ "decreased from",
     round(cv_median_adj,2) - round(cv_median_raw ,2) == 0 ~ "remained the same at",
     TRUE ~ "remained similar at")

  text_start <- if(conditional_correction) "Conditional drift" else "Drift"

  text_batchwise <- if(batch_wise) " across batches " else " "

  cli_alert_success(
    col_green(
      c("{text_start} correction was applied to {count_feature_text} {mode_text}.")))

  cli_alert_info(cli::col_grey(
           "The median CV change of all features in study samples was {.strong {formatC(cv_difference_median, format = 'f', digits = 2)}%} (range: {formatC(cv_difference_low, format = 'f', digits = 2)}% to {formatC(cv_difference_high, format = 'f', digits = 2)}%). The median absolute CV of all features{text_batchwise}{.strong {text_change}} {.strong {formatC(cv_median_raw, format = 'f', digits = 2)}% {ifelse(str_detect(text_change, 'remained'), '',  paste0('to ', formatC(cv_median_adj, format = 'f', digits = 2),'%'))}}."))


  # Invalidate downstream processed data
  if(variable == "feature_intensity"){
    data <- update_after_normalization(data, FALSE)
  } else if(variable == "feature_norm_intensity") {
    data <- update_after_quantitation(data, FALSE)
  }

  data@status_processing <- "Drift-corrected concentrations"
  data@var_drift_corrected[[variable]] <- TRUE
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]
  data
}

#' Drift Correction by Gaussian Kernel Smoothing
#'
#' @description
#' Performs drift correction for run-order effects within or across batches
#' using Gaussian kernel smoothing, as detailed in Teo et al. (2020). The
#' Gaussian kernel estimates the local data trend, with bandwidth defined by
#' the `kernel_size` parameter. This smoothing approach is mostly used with
#' study samples and should only be applied to datasets with sufficiently randomized or
#' stratified samples to avoid local biases and artifacts. The smoothing can be
#' applied  to `concentration`, `norm_intensity`, and `intensity` data.
#'
#' Corrections can be applied on a batch-by-batch basis (`batch_wise = TRUE`,
#' default) or across all batches (`batch_wise = FALSE`). The correction can
#' either replace existing drift or batch corrections (`replace_previous =
#' `TRUE`, default) or applied on top of existing corrections (`replace_previous = FALSE`).
#'
#' Drift correction can be applied to all features (`conditional_correction = FALSE`)
#' or conditionally, based on whether the sample CV difference before and
#' after correction is below a defined threshold (`cv_diff_threshold`). The
#' conditional correction is applied separately for each batch if
#' `batch_wise = TRUE`, .
#'
#' It is recommended to visually inspect the correction using the
#' [plot_runscatter()] function. Set the argument
#' `recalc_trend_after = TRUE` so that the trends after correction are also
#' available for plotting. For further details, refer to the description
#' of [plot_runscatter()]. This will double the processing time.
#'
#' **Note**: The function outputs a message indicating the median CV change
#' and the mean absolute CV before and after correction for all samples.
#' However, these metrics are experimental and should not be used as
#' definitive criteria for correction (see Details below).
#'
#' @details
#' In the output message, the median CV change is computed as the median of CV changes for all
#' features in global correction or for
#' features where the correction passed the defined CV difference treshold in
#' case of conditional correction  (`conditional_correction = FALSE`).
#' For batch-wise correction, the change is calculated per batch, with the final median CV
#' change being the median of these batch medians across features.
#'
#' @param data A MidarExperiment object.
#' @param variable The target variable for drift correction; options include
#' "intensity", "norm_intensity", or "conc".
#' @param reference_qc_types QC types used for drift correction, typically
#' including study samples (`SPL`).
#' @param batch_wise Logical. Apply the correction to each batch separately
#' (`TRUE`, default) or across all batches (`FALSE`).
#' @param ignore_istd Logical. Exclude internal standards (ISTDs) from
#' correction if `TRUE`.
#' @param replace_previous Logical. Replace existing correction (`TRUE`,
#' default) or layer on top of it (`FALSE`).
#' @param kernel_size Numeric. Defines the Gaussian kernel's bandwidth.
#' @param outlier_filter Logical. Enable kernel outlier filtering if `TRUE`.
#' @param outlier_ksd Numeric. Set the kernel's k times standard deviation for
#' outlier detection.
#' @param location_smooth Logical. Apply smoothing to the location parameter if
#' `TRUE`.
#' @param scale_smooth Logical. Apply smoothing to the scale parameter if
#' `TRUE`.
#' @param log_transform_internal Logical. Conduct log transformation internally
#' for enhanced outlier robustness if `TRUE` (default); does not alter the
#' output data.
#' @param recalc_trend_after Logical. Recalculate trends post-smoothing for
#' visualization (e.g., `plot_runscatter()`).

#' @param conditional_correction Determines whether drift correction should be applied to all features unconditionally (`TRUE`) or
#' only when the difference of sample CV before vs after smoothing is below the threshold specified by `cv_diff_threshold`.
#' @param cv_diff_threshold This parameter defines the maximum allowable change (difference) in the coefficient of variation (CV) ratio of samples before and after smoothing for the correction to be applied.
#' A value of 0 (the default) requires the CV to improve, while a value above 0 allows the CV to also become worse by a maximum of the defined difference.
#' @param feature_list Character vector. Regular expression pattern to select
#' specific features for correction. Default is `NULL` for all features.
#' @param use_original_if_fail Determines the action when smoothing fails or results in invalid values for a feature. If `FALSE` (default), the result for each feature will `NA` for all batches, if `TRUE`, the original data is kept.
#' @param show_progress Logical. Display progress bars if `TRUE`; disable for
#' notebook rendering by setting to `FALSE`.
#' @return Returns a MidarExperiment object.
#' @references
#' Teo G., Chew WS, Burla B, Herr D, Tai ES, Wenk MR, Torta F, & Choi H
#' (2020). MRMkit: Automated Data Processing for Large-Scale Targeted
#' Metabolomics Analysis. *Analytical Chemistry*, 92(20), 13677–13682.
#' \url{https://doi.org/10.1021/acs.analchem.0c03060}
#' @export
correct_drift_gaussiankernel <- function(data = NULL,
                                      variable,
                                      reference_qc_types,
                                      batch_wise = TRUE,
                                      ignore_istd = TRUE,
                                      replace_previous = TRUE,
                                      kernel_size = 10,
                                      outlier_filter = FALSE,
                                      outlier_ksd = 5,
                                      location_smooth = TRUE,
                                      scale_smooth = FALSE,
                                      log_transform_internal = TRUE,
                                      conditional_correction = FALSE,
                                      cv_diff_threshold = 0,
                                      recalc_trend_after = FALSE,
                                      feature_list = NULL,
                                      use_original_if_fail = FALSE,
                                      show_progress = TRUE
) {

  check_data(data)

  if(is.na(kernel_size) || kernel_size <= 0) {
    cli_abort(col_red("Argument `kernel_size` must larger than 0."))
  }

  if(is.na(outlier_ksd) || outlier_ksd <= 0) {
    cli_abort(col_red("Argument `outlier_ksd` must larger than 0."))
  }

  if(!log_transform_internal) {
    cli_abort(col_red("Currently only `log_transform_internal = TRUE` is supported."))
  }


  fun_correct_drift(
    data = data,
    smooth_fun = "fun_gauss.kernel.smooth",
    variable = variable,
    reference_qc_types = reference_qc_types,
    batch_wise = batch_wise ,
    replace_previous = replace_previous,
    conditional_correction = conditional_correction,
    log_transform_internal = log_transform_internal,
    cv_diff_threshold = cv_diff_threshold,
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
#' This function corrects for run-order drifts within or across batches using
#' LOESS (Locally Estimated Scatterplot Smoothing). The correction is typically
#' based on QC (Quality Control) samples that were measured at specific intervals
#' throughout the run sequence. The smoothed curve derived from the QC samples is
#' then used to adjust all other samples in the dataset. The correction can be applied
#' to "intensity", "norm_intensity", or "conc" data.
#'
#' The degree of smoothing is controlled by the span parameter `span` (default is 0.75).
#' Additionally,
#' the `degree` parameter can be specified to
#' control the degree of the polynomial used in the local regression (default is 2)
#'
#' It is recommended to visually inspect the correction using the
#' [plot_runscatter()] function. Set the argument
#' `recalc_trend_after = TRUE` so that the trends after correction are also
#' available for plotting. For further details, refer to the description
#' of [plot_runscatter()].
#'
#' The LOESS correction only applies to samples that lie within the span of the
#' QC samples used for smoothing. Extrapolation outside this range is
#' not recommended, as it can lead to unreliable corrections or artefacts in the extrapolated regions. However,
#' extrapolation can be activated by setting `extrapolate = TRUE`. This may be
#' useful in cases where specific drifts occur in segments of the analysis
#' sequence that are not spanned by the QC samples, such as when the
#' analysis was interrupted or the instrument had rapid changes in performance.
#'
#' The corrections can be applied on a batch-by-batch basis (`batch_wise =
#' TRUE`, default) or across all batches (`batch_wise = FALSE`). Existing
#' corrections are either replaced (`replace_previous = TRUE`) or added on top
#' of them (`replace_previous = FALSE`).
#'
#' Furthermore, drift correction can be applied unconditionally
#' (`conditional_correction = FALSE`) or conditionally, based on whether the
#' sample CV change before and after correction is below a defined
#' threshold (`cv_diff_threshold`). This conditional correction is assessed
#' independently for each batch if `batch_wise = TRUE`, where the median of
#' the CV changes across the batch is compared with the threshold.
#'
#' **Note**: The function outputs a message indicating the median CV change
#' and the mean absolute CV before and after correction for all samples.
#' However, these metrics are experimental and should not be used as
#' definitive criteria for correction (see Details below).
#'
#' This LOESS method is implemented using the base R function
#' [stats::loess()].
#'
#' @details
#' In the output message, the median CV change is computed as the median of CV changes for all
#' features in global correction or for
#' features where the correction passed the defined CV difference treshold in
#' case of conditional correction  (`conditional_correction = FALSE`).
#' For batch-wise correction, the change is calculated per batch, with the final median CV
#' change being the median of these batch medians across features.
#'
#' @param data MidarExperiment object
#' @param reference_qc_types QC types used for drift correction
#' @param variable  The variable to be corrected for drift effects. Must be one of "intensity", "norm_intensity", or "conc"
#' @param batch_wise Logical. Apply the correction to each batch separately
#' (`TRUE`, default) or across all batches (`FALSE`).
#' @param ignore_istd Logical. Exclude internal standards (ISTDs) from
#' correction if `TRUE`.
#' @param replace_previous Logical. Replace existing correction (`TRUE`,
#' default) or layer on top of it (`FALSE`).
#' @param span Loess span width (default is 0.75)
#' @param degree Degree of the polynomial to be used in the loess smoothing, normally 1 (default) or 2
#' @param extrapolate Extrapolate loess smoothing. WARNING: It is generally not recommended to extrapolate outside of the range spanned by the QCs used for smoothing. See details below.
#' @param log_transform_internal Log transform the data for correction when `TRUE` (the default). Note: log transformation is solely applied internally for smoothing, results will not be be log-transformed. Log transformation may result in more robust smoothing that is less sensitive to outlier.
#' @param conditional_correction Determines whether drift correction should be applied to all features unconditionally (`TRUE`) or
#' only when the difference of sample CV before vs after smoothing is below the threshold specified by `cv_diff_threshold`.
#' @param recalc_trend_after Recalculate trend post-drift correction for `plot_qc_runscatter()`. This will double calculation time.
#' @param cv_diff_threshold This parameter defines the maximum allowable change (difference) in the coefficient of variation (CV) ratio of samples before and after smoothing for the correction to be applied.
#' A value of 0 (the default) requires the CV to improve, while a value above 0 allows the CV to also become worse by a maximum of the defined difference.
#' @param feature_list Subset the features for correction whose names matches the specified text using regular expression. Default is `NULL` which means all features are selected.
#' @param ignore_istd Do not apply corrections to ISTDs
#' @param use_original_if_fail Determines the action when smoothing fails or results in invalid values for a feature. If `FALSE` (default), the result for each feature will `NA` for all batches, if `TRUE`, the original data is kept.
#' @return MidarExperiment object
#' @export
correct_drift_loess <- function(data = NULL,
                                variable,
                                reference_qc_types,
                                batch_wise = TRUE,
                                ignore_istd = TRUE,
                                replace_previous = TRUE,
                                conditional_correction = FALSE,
                                recalc_trend_after = FALSE,
                                log_transform_internal = TRUE,
                                feature_list = NULL,
                                cv_diff_threshold = 0,
                                use_original_if_fail = FALSE,
                                extrapolate = FALSE,
                                span = 0.75,
                                degree = 2) {
  check_data(data)

    if(is.na(span) || span <= 0) {
      cli_abort(col_red("Argument `span` must larger than 0, typically between 0.25 and 1.0."))
    }

  if(is.na(degree) || degree < 1 || degree > 2) {
    cli_abort(col_red("Argument `degree` must be 1 or 2."))
  }

  fun_correct_drift(
    data = data,
    smooth_fun = "fun_loess",
    reference_qc_types = reference_qc_types,
    variable = variable,
    batch_wise = batch_wise,
    replace_previous = replace_previous,
    extrapolate = extrapolate,
    span = span,
    degree = degree,
    conditional_correction = conditional_correction,
    log_transform_internal = log_transform_internal,
    feature_list = feature_list,
    ignore_istd = ignore_istd,
    recalc_trend_after = recalc_trend_after,
    cv_diff_threshold = cv_diff_threshold,
    use_original_if_fail = use_original_if_fail,

  )
}

#' Drift Correction by Cubic Spline Smoothing
#' @description
#' This function corrects for run-order drifts within or across batches using
#' cubic spline smoothing. The correction is typically based on QC (Quality
#' Control) samples that are measured at specific intervals throughout the run
#' sequence. The smoothed curve derived from the QC samples is then used to
#' adjust all other samples in the dataset. The correction can be applied to "intensity",
#' "norm_intensity", or "conc" data.
#'
#' The cubic spline smoothing approach, particularly when used with the
#' regularization parameter `lambda`, is similar but not identical to previously described
#' QC-based drift correction methods, such as **QC-RSC (Quality Control Regularized
#' Spline Correction)**, described in Dunn et al. (Nat Protoc, 2011) and Kirwan et al.
#' (Anal Bioanal Chem, 2014).
#'
#' By default, the smoothing parameter is determined using cross-validation,
#' which can lead to overfitting. To reduce overfitting
#' the regularization parameter `lambda` may be defined, with a good starting point
#' being `lambda = 0.01`. Additionally, the global
#' smoothing parameter can be specified via `spar`.
#'
#' It is recommended to visually inspect the correction using the
#' [plot_runscatter()] function. Set the argument
#' `recalc_trend_after = TRUE` so that the trends after correction are also
#' available for plotting. For further details, refer to the description
#' of [plot_runscatter()].
#'
#' The corrections can be applied on a batch-by-batch basis (`batch_wise =
#' TRUE`, default) or across all batches (`batch_wise = FALSE`). Existing
#' corrections are either replaced (`replace_previous = TRUE`) or added on top
#' of them (`replace_previous = FALSE`).
#'
#' Furthermore, drift correction can be applied unconditionally
#' (`conditional_correction = FALSE`) or conditionally, based on whether the
#' sample CV change before and after correction is below a defined
#' threshold (`cv_diff_threshold`). This conditional correction is assessed
#' independently for each batch if `batch_wise = TRUE`, where the median of
#' the CV changes across the batch is compared with the threshold.
#'
#' **Note**: The function outputs a message indicating the median CV change
#' and the mean absolute CV before and after correction for all samples.
#' However, these metrics are experimental and should not be used as
#' definitive criteria for correction (see Details below).
#'
#' This cubic spline method is implemented using the base R function
#  \code{\link[stats]{smooth.spline}}
#'
#' @details
#' In the output message, the median CV change is computed as the median of CV changes for all
#' features in global correction or for
#' features where the correction passed the defined CV difference treshold in
#' case of conditional correction  (`conditional_correction = FALSE`).
#' For batch-wise correction, the change is calculated per batch, with the final median CV
#' change being the median of these batch medians across features.
#'

#'
#' @param data MidarExperiment object
#' @param reference_qc_types QC types used for drift correction
#' @param variable  The variable to be corrected for drift effects. Must be one of "intensity", "norm_intensity", or "conc"
#' @param batch_wise Logical. Apply the correction to each batch separately (`TRUE`, default) or across all batches (`FALSE`).
#' @param ignore_istd Logical. Exclude internal standards (ISTDs) from correction if `TRUE`.
#' @param replace_previous Logical. Replace existing correction (`TRUE`, default) or layer on top of it (`FALSE`).
#' @param cv Ordinary leave-one-out (TRUE) or ‘generalized’ cross-validation (GCV) when FALSE; is used for smoothing parameter computation only when spar is not specified
#' @param spar Smoothing parameter for cubic spline smoothing. If not specified or `NULL`, the smoothing parameter is computed using the specified cv method. Typically (but not necessarily) in in (0,1].
#' @param lambda Regularization parameter for cubic spline smoothing. Default is 0, which means no regularization.
#' @param penalty The coefficient of the penalty for degrees of freedom in the GCV criterion.
#' @param log_transform_internal Log transform the data for correction when `TRUE` (the default). Note: log transformation is solely applied internally for smoothing, results will not be log-transformed.
#' @param conditional_correction Determines whether drift correction should be applied to all features unconditionally (`TRUE`) or conditionally, based on sample CV change.
#' @param recalc_trend_after Recalculate trend post-drift correction for `plot_qc_runscatter()`. This will double calculation time.
#' @param cv_diff_threshold Maximum allowable change in CV ratio before and after smoothing for correction to be applied.
#' @param feature_list Subset the features for correction whose names match the specified text using regular expression. Default is `NULL`.
#' @param use_original_if_fail Determines the action when smoothing fails or results in invalid values for a feature. If `FALSE` (default), the result for each feature will `NA` for all batches, if `TRUE`, the original data is kept.
#' @return MidarExperiment object
#' @export
#' @references  Dunn, W., Broadhurst, D., Begley, P. et al. Procedures for
#' large-scale metabolic profiling of serum and plasma using gas chromatography
#' and liquid chromatography coupled to mass spectrometry.
#' Nat Protoc 6, 1060–1083 (2011).
#' \url{https://doi.org/10.1038/nprot.2011.335}
#'
#' Kirwan, J.A., Broadhurst, D.I., Davidson, R.L. et al. Characterising and correcting
#' batch variation in an automated direct infusion mass spectrometry (DIMS) metabolomics workflow.
#' Anal Bioanal Chem 405, 5147–5157 (2013). https://doi-org.libproxy1.nus.edu.sg/10.1007/s00216-013-6856-7
#'
#' @seealso \code{\link[stats]{smooth.spline}}
correct_drift_cubicspline <- function(data = NULL,
                                 variable,
                                 reference_qc_types,
                                 batch_wise = TRUE,
                                 ignore_istd = TRUE,
                                 replace_previous = TRUE,
                                 cv = TRUE,
                                 spar = NULL,
                                 lambda = NULL,
                                 penalty = 1,
                                 conditional_correction = FALSE,
                                 recalc_trend_after = FALSE,
                                 log_transform_internal = TRUE,
                                 feature_list = NULL,
                                 cv_diff_threshold = 0,
                                 use_original_if_fail = FALSE
                                 ) {
  check_data(data)

  if(!is.null(spar) && !(is.numeric(spar))) {
    cli_abort(col_red("Argument `spar` must be NULL or numeric, typically between 0 and 1."))
  }

  if(!is.null(spar) && !(is.null(lambda))) {
    cli_abort(col_red("Either `spar` or `lambda` can be specified, not both."))
  }

  fun_correct_drift(
    data = data,
    reference_qc_types = reference_qc_types,
    variable = variable,
    batch_wise = batch_wise,
    replace_previous = replace_previous,
    conditional_correction = conditional_correction,
    log_transform_internal = log_transform_internal,
    feature_list = feature_list,
    ignore_istd = ignore_istd,
    recalc_trend_after = recalc_trend_after,
    cv_diff_threshold = cv_diff_threshold,
    use_original_if_fail = use_original_if_fail,
    smooth_fun = "fun_cspline",
    cv = cv,
    spar = spar,
    lambda = lambda,
    penalty = penalty
  )
}

#' Drift Correction by Generalized Additive Model (GAM) Smoothing
#' @description
#' This function corrects for run-order drifts within or across batches using
#' Generalized Additive Models (GAMs). The correction uses penalized splines,
#' with automatic selection of smoothing parameters based on cross-validation or
#' penalized likelihood. It is typically based on QC (Quality Control) samples
#' measured at specific intervals throughout the run sequence.
#' The correction can be applied to "intensity",
#' "norm_intensity", or "conc" data.
#'
#' It is recommended to visually inspect the correction using the
#' [plot_runscatter()] function. Set the argument
#' `recalc_trend_after = TRUE` so that the trends after correction are also
#' available for plotting. For further details, refer to the description
#' of [plot_runscatter()].
#'
#' The corrections can be applied on a batch-by-batch basis (`batch_wise =
#' TRUE`, default) or across all batches (`batch_wise = FALSE`). Existing
#' corrections are either replaced (`replace_previous = TRUE`) or added on top
#' of them (`replace_previous = FALSE`).
#'
#' Furthermore, drift correction can be applied unconditionally
#' (`conditional_correction = FALSE`) or conditionally, based on whether the
#' sample CV change before and after correction is below a defined
#' threshold (`cv_diff_threshold`). This conditional correction is assessed
#' independently for each batch if `batch_wise = TRUE`, where the median of
#' the CV changes across the batch is compared with the threshold.
#'
#' **Note**: The function outputs a message indicating the median CV change
#' and the mean absolute CV before and after correction for all samples.
#' However, these metrics are experimental and should not be used as
#' definitive criteria for correction (see Details below).
#'
#' This cubic spline method is implemented using the base R function
#' [stats::spline()].
#'
#' @details
#' In the output message, the median CV change is computed as the median of CV changes for all
#' features in global correction or for
#' features where the correction passed the defined CV difference treshold in
#' case of conditional correction  (`conditional_correction = FALSE`).
#' For batch-wise correction, the change is calculated per batch, with the final median CV
#' change being the median of these batch medians across features.
#'
#' This smoothing is based on Generalized Additive Models (GAM) using penalized splines, implemented via `mgcv::gam()`.
#'
#' @param data MidarExperiment object
#' @param reference_qc_types QC types used for drift correction
#' @param variable  The variable to be corrected for drift effects. Must be one of "intensity", "norm_intensity", or "conc"
#' @param batch_wise Logical. Apply the correction to each batch separately (`TRUE`, default) or across all batches (`FALSE`).
#' @param ignore_istd Logical. Exclude internal standards (ISTDs) from correction if `TRUE`.
#' @param replace_previous Logical. Replace existing correction (`TRUE`, default) or layer on top of it (`FALSE`).
#' @param bs Basis type for the spline: `"ps"` (penalized spline, default) or `"tp"` (thin plate spline).
#' @param k Number of basis functions (default: `-1`, automatically chosen by GAM).
#' @param sp Smoothing parameter (`NULL` by default, estimated automatically).
#' @param log_transform_internal Log transform the data for correction when `TRUE` (the default). Note: log transformation is solely applied internally for smoothing, results will not be log-transformed.
#' @param conditional_correction Determines whether drift correction should be applied to all features unconditionally (`TRUE`) or conditionally, based on sample CV change.
#' @param recalc_trend_after Recalculate trend post-drift correction for `plot_qc_runscatter()`. This will double calculation time.
#' @param cv_diff_threshold Maximum allowable change in CV ratio before and after smoothing for correction to be applied.
#' @param feature_list Subset the features for correction whose names match the specified text using regular expression. Default is `NULL`.
#' @param use_original_if_fail Determines the action when smoothing fails or results in invalid values for a feature. If `FALSE` (default), the result for each feature will `NA` for all batches, if `TRUE`, the original data is kept.
#' @return MidarExperiment object
#' @export
#'
#' @seealso [mgcv::gam()]
correct_drift_gam <- function(data = NULL,
                              variable,
                              reference_qc_types,
                              batch_wise = TRUE,
                              ignore_istd = TRUE,
                              replace_previous = TRUE,
                              bs = "ps",
                              k = -1,
                              sp = NULL,
                              log_transform_internal = TRUE,
                              conditional_correction = FALSE,
                              recalc_trend_after = FALSE,
                              feature_list = NULL,
                              cv_diff_threshold = 0,
                              use_original_if_fail = FALSE) {

  # {ggpmisc} neeeded for plots
  check_installed("mgcv")

  check_data(data)

  if (!is.null(sp) && !(is.numeric(sp))) {
    cli_abort(col_red("Argument `sp` must be NULL or numeric."))
  }

  fun_correct_drift(
    data = data,
    reference_qc_types = reference_qc_types,
    variable = variable,
    batch_wise = batch_wise,
    replace_previous = replace_previous,
    conditional_correction = conditional_correction,
    log_transform_internal = log_transform_internal,
    feature_list = feature_list,
    ignore_istd = ignore_istd,
    recalc_trend_after = recalc_trend_after,
    cv_diff_threshold = cv_diff_threshold,
    use_original_if_fail = use_original_if_fail,
    smooth_fun = "fun_gam_smooth",
    bs = bs,
    k = k,
    sp = sp
  )
}


#' Batch Centering Correction
#' @description
#' This function performs batch centering correction on each feature.
#' Optionally, the scale of the batches can be equalized across batches.
#' The selected QC types (`reference_qc_types`) are used to calculate
#' the medians, which are then used to align all other samples. The
#' correction can be applied to one of three variables: "intensity",
#' "norm_intensity", or "conc". The correction can either be applied
#' on top of previous corrections or replace all prior batch corrections.
#'
#' @param data A `MidarExperiment` object containing the data to be corrected.
#'   This object must include information about QC types and measurements.
#' @param variable The variable to be corrected. Must be one of "intensity",
#'   "norm_intensity", or "conc".
#' @param reference_qc_types A character vector specifying the QC types to be
#'   used as references for batch centering.
#' @param correct_scale A logical value indicating whether to equalize the scale
#' of the batches in addition to center them. Defaults to `FALSE`.
#' @param replace_previous A logical value indicating whether to replace any
#'   previous batch corrections or apply the new correction on top. Defaults to
#'   `TRUE` (replace).
#' @param log_transform_internal A logical value indicating whether to log-transform
#'   the data internally during correction. Defaults to `TRUE`.
#' @param replace_exisiting_trendcurves A logical value indicating whether to replace
#' trend curves from previous corrections. This is only use for plotting using `plot_runscatter()`. Default is `FALSE`.
#' @param ... Additional arguments that can be passed to the batch correction
#'   function.
#' @return A `MidarExperiment` object containing the corrected data.
#' @seealso `plot_runscatter` for visualizing the correction before and after.
#' @export

correct_batch_centering <- function(data = NULL,
                                    variable,
                                    reference_qc_types,
                                   #correct_location = TRUE,
                                    correct_scale = FALSE,
                                    replace_previous = TRUE,
                                    log_transform_internal = TRUE,
                                    replace_exisiting_trendcurves = FALSE,
                                    ...) {
  check_data(data)

  if (!all(reference_qc_types %in% unique(data@dataset$qc_type))) {
    cli::cli_abort(col_red("One or more specified `qc_types` are not present in the dataset. Please verify data or analysis metadata."))
  }

  if(!log_transform_internal && correct_scale){
    cli_abort(col_red("Currently data must be log-transformed for batch scaling. Please set `log_transform_internal = TRUE`"))
  }

  variable_strip <- str_remove(variable, "feature_")
  rlang::arg_match(variable_strip, c("intensity", "norm_intensity", "conc"))
  variable <- stringr::str_c("feature_", variable_strip)
  variable_sym <- rlang::sym(variable)
  variable_raw <- paste0(variable, "_raw")
  variable_before <- paste0(variable, "_before")
  variable_before_sym <- rlang::sym(variable_before)
  variable_fit <- paste0(variable_before, "_fit")
  variable_raw_fit <- paste0(variable_raw, "_fit")
  variable_fit_sym <- rlang::sym(variable_fit)
  variable_raw_fit_sym <- rlang::sym(variable_raw_fit)
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

  is_first_correction <- FALSE

  if(data@var_drift_corrected[[variable]] ){
    # Drift-corrected data
    if (data@var_batch_corrected[[variable]] ){
      # Drift  AND batch-corrected data
      if(replace_previous){
        cli_alert_info(col_yellow(glue::glue("Replacing previous `{variable_strip}` batch correction of drift-corrected data.")))
        # use drift corrected data and apply correction on top, replacing any previous batch correction
        data@dataset[[variable]] <- data@dataset[[variable_smoothed]]
        data@dataset[[variable_fit_after]] <- data@dataset[[variable_smoothed_fit_after]]
      } else {
        cli_alert_info(col_yellow(glue::glue("Adding batch correction on top of previous `{variable_strip}` drift and batch corrections.")))
        # use main variable (e.g. conc) and add correction on top
        #data@dataset[[variable]] <- data@dataset[[variable]]
      }
    } else {
      # Drift but NOT batch-corrected data
      # replace_previous will have no effect as no previous batch correction has beem applied


      data@dataset[[variable_smoothed]] <- data@dataset[[variable]]
      data@dataset[[variable_smoothed_fit_after]] <- data@dataset[[variable_fit_after]]

      cli_alert_info(col_yellow(glue::glue("Adding batch correction on top of `{variable_strip}` drift-correction.")))

    }
  } else {
    # Data is NOT drift corrected
    if (data@var_batch_corrected[[variable]]){
        # Only batch-corrected data (no previous drift correction)
        if(replace_previous){
          cli_alert_info(col_yellow(glue::glue("Replacing previous `{variable_strip}` batch correction.")))
            # use drift corrected data and apply correction on top, replacing any previous batch correction
            data@dataset[[variable]] <- data@dataset[[variable_before]]
            data@dataset[[variable_fit_after]] <- data@dataset[[variable_fit]]

        } else {
          cli_alert_info(col_yellow(glue::glue("Adding batch correction on top of previous `{variable_strip}` batch correction.")))
            # use main variable (e.g. conc) and add correction on top
            #data@dataset[[variable]] <- data@dataset[[variable]]
        }

    } else {
        # Uncorrected data (raw)
      cli_alert_info(col_yellow(glue::glue("Adding batch correction to `{variable_strip}` data.")))
        #Store raw data
        is_first_correction <- TRUE
        data@dataset[[variable_raw]] <- data@dataset[[variable]]


        # Use median horizontal lines to indicate fits if data was not drift corrected before
        data@dataset <- data@dataset |>
            mutate(!!variable_fit_sym := median(.data[[variable_sym]][.data$qc_type %in% reference_qc_types ], na.rm = TRUE), .by = c("feature_id", "batch_id"))
        data@dataset[[variable_fit_after]] <- data@dataset[[variable_fit]]
      }

  }

  if(replace_exisiting_trendcurves){
    data@dataset <- data@dataset |>
      mutate(!!variable_fit_sym := median(.data[[variable_sym]][.data$qc_type %in% reference_qc_types ], na.rm = TRUE), .by = c("feature_id", "batch_id"))
    data@dataset[[variable_fit_after]] <- data@dataset[[variable_fit]]
  }

  ds <- data@dataset |>
    select(any_of(c("analysis_id", "feature_id", "qc_type", "batch_id", y_fit_after = variable_fit_after, y = variable)))
  nbatches <- length(unique(ds$batch_id))
  if(nbatches < 2) {
    cli_abort(col_red(glue::glue("Batch correction was not applied as there is only one batch.")))
  }

  d_res <- ds |>
    group_by(.data$feature_id) |>
      nest() |>
      mutate(
        res = purrr::map(data, \(x) do.call("fun_batch.correction", list(x, reference_qc_types, correct_scale, log_transform_internal,...))),
      ) |>
      unnest(cols = c("res")) |>
      select(-"data")

  d_res_sum <- d_res |>
    group_by(.data$feature_id) |>
    summarise(
      cv_before = cv(.data$y[.data$qc_type == "SPL"], na.rm = TRUE),
      cv_after = cv(.data$y_adj[.data$qc_type == "SPL"], na.rm = TRUE),
      cv_diff = .data$cv_after - .data$cv_before,
    ) |>
    ungroup() |>
    summarise(
      cv_before = median(.data$cv_before, na.rm = TRUE),
      cv_after = median(.data$cv_after, na.rm = TRUE),
      cv_diff_median = median(.data$cv_diff, na.rm = TRUE),
      cv_diff_min = format(round(min(.data$cv_diff, na.rm = TRUE), 1), nsmall = 2),
      cv_diff_max = format(round(max(.data$cv_diff, na.rm = TRUE), 1), nsmall = 2),
      cv_diff_text = format(round(.data$cv_diff_median, 1), nsmall = 1)
    ) |> ungroup()

  nfeat <- get_feature_count(data, is_istd = FALSE)

  # Print summary
  if (data@var_drift_corrected[[variable]]) {
    cli_alert_success(col_green(glue::glue("Batch median-centering of {nbatches} batches was applied to drift-corrected concentrations of all {nfeat} features.")))
    data@status_processing <- "Batch- and drift-corrected concentrations"
  } else {
    cli_alert_success(col_green(glue::glue("Batch median-centering of {nbatches} batches was applied to raw concentrations of all {nfeat} features.")))
    data@status_processing <- "Batch-corrected concentrations"
  }
  # Print stats


  text_change <- case_when(
    round(d_res_sum$cv_after,2) - round(d_res_sum$cv_before ,2) >= 0.01 ~ "increased from",
    round(d_res_sum$cv_after,2) - round(d_res_sum$cv_before ,2) <= -0.01 ~ "decreased from",
    round(d_res_sum$cv_after,2) - round(d_res_sum$cv_before ,2) == 0 ~ "remained the same at",
    TRUE ~ "remained similar at")

  # cli_alert_info(cli::col_grey(
  #     c("The median CV of features in the study samples across batches {.strong {text_change}} by {d_res_sum$cv_diff_text}% (range: {d_res_sum$cv_diff_min}% to {d_res_sum$cv_diff_max}%), from {format(round(d_res_sum$cv_before,1), nsmall = 1)}%  to {format(round(d_res_sum$cv_after,1), nsmall = 1)}%.")))

  cli_alert_info(cli::col_grey(
    "The median CV change of all features in study samples was {.strong {formatC(d_res_sum$cv_diff_median, format = 'f', digits = 2)}%} (range: {formatC(d_res_sum$cv_diff_min, format = 'f', digits = 2)}% to {formatC(d_res_sum$cv_diff_max, format = 'f', digits = 2)}%).  The median absolute CV of all features {.strong {text_change}} {.strong {formatC(d_res_sum$cv_before, format = 'f', digits = 2)}% {ifelse(str_detect(text_change, 'remained'), '',  paste0('to ', formatC(d_res_sum$cv_after, format = 'f', digits = 2),'%'))}}."))


  # Return data


  data@dataset <- data@dataset |>
    left_join(d_res |> select(  -"y_fit_after"),
              by = c("analysis_id", "feature_id", "qc_type", "batch_id")) |>
    mutate(!!variable_sym := .data$y_adj,
           !!variable_before_sym := .data$y,
           !!variable_fit_sym := .data[[variable_fit]],
           !!variable_raw_fit_sym := if(is_first_correction) .data[[variable_fit_after]] else !!variable_raw_fit_sym,
           !!variable_fit_after_sym := .data$y_fit_after_adj) |>
    select(-"y_adj", -"y_fit_after_adj", -"y")




  # Invalidate downstream processed data
  if(variable == "feature_intensity"){
    data <- update_after_normalization(data, FALSE)
  } else if(variable == "feature_norm_intensity") {
    data <- update_after_quantitation(data, FALSE)
  }

  data@status_processing <- if(any(data@var_drift_corrected)) "Drift-Batch-corrected concentrations" else "Batch-corrected concentrations"
  data@var_batch_corrected[[variable]] <- TRUE
  data@is_filtered <- FALSE
  data@metrics_qc <- data@metrics_qc[FALSE,]
  data
}



fun_batch.correction = function(tab,
                            reference_qc_types,
                            correct_scale,
                            log_transform_internal, ...) {

  batch <- tab$batch_id
  batch.order <- seq(1, nrow(tab))
  val <- tab$y
  y_fit_after <- tab$y_fit_after  # BB
  if(log_transform_internal){
    val <- log10(val)
    y_fit_after <- log10(y_fit_after)
  }
  sample.for.loc <- tab$qc_type %in% reference_qc_types # Use sample for location median

  ubatch <- unique(batch)
  nbatch <- length(ubatch)

  val.clean <- val ## placeholder
  y_fit_after.clean <- y_fit_after # BB



  ### Cross-batch scale normalization
  if(!correct_scale) {
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
  } else {
    tmp <- val.clean
    tmp_fit_after <- y_fit_after.clean # BB
    tmp_for_loc <- val.clean
    tmp_for_loc[!sample.for.loc] <- NA_real_

    # TODO: confirm if this is ok to do and its consequences. Once batch is NA, then shoudnt we report all?
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
        cli_abort(col_red("Currently data must be log-transformed for batch scaling. Please set `log_transform_internal = TRUE`"))
      }
    }
  }

    tab$y_adj <- if(log_transform_internal) 10^val.clean else val.clean
    tab$y_fit_after_adj <- if(log_transform_internal) 10^y_fit_after.clean else  y_fit_after.clean
  tab

}

