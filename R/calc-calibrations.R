#' Calculate concentrations based on external calibration
#'
#' Concentrations of all features in all analyses are determined using ISTD-normalized intensities and corresponding external calibration curves.
#' Calibration curves are calculated for each feature based on calibration sample concentrations defined in the `qc_concentrations` metadata.
#' The regression fit model (linear or quadratic) and the weighting method (either "none", "1/x", or "1/x^2") can be defined globally via
#' the arguments `fit_method` and `fit_weighting` for all features, if `overwrite_metadata` is `TRUE`.
#' Alternatively, the model and weighting can be defined individually for each feature in the `feature` metadata (columns `curve_fit_method` and `fit_weighting`).
#' If these details are missing in the metadata, the default values provided via `fit_method` and `fit_weighting` will be used.
#'
#' The concentrations are added to the `dataset` table as `feature_conc` column. The results of the regression and the calculated LoD and LoQ values are stored in the `metrics_calibration` table of the returned `MidarExperiment` object.

#' @param data A `MidarExperiment` object
#' @param overwrite_metadata If `TRUE`,  any fit method and weighting settings defined in the metadata will be ignored and instead the `fit_method` and `fit_weighting` are used for all features
#' @param fit_method A character string indicating the default regression fit method to use for the calibration curve. Must be one of `"linear"` or `"quadratic"`. This method will be used if no specific fit method is defined for a feature in the metadata.
#' @param fit_weighting A character string indicating the default weighting method for the regression points in the calibration curve. Must be one of `"none"`, `"1/x"`, or `"1/x^2"`. If no specific weighting method is defined for a feature in the metadata, this method will be used.
#' @param error_failed_calibration If `TRUE`, an error will be raised if the calibration curve fitting failed for any feature. If `FALSE`, failed calibration curve fitting will be ignored, and resulting feature concentration will be `NA`.
#' @param fail_missing_annotation Raise error if any of the following information is missing: calibration curve data, ISTD mix volume and sample amounts for any feature.
#'   If `FALSE`, missing annotations will be ignored, and resulting feature concentration will be `NA`
#' @return A modified `MidarExperiment` object with updated concentration values.
#'
#' @seealso [calc_calibration_results()] for calculating the calibration curve results including LoD and LoQ.
#' @seealso [quantify_by_istd()] for calculation of concentrations based on spiked-in internal standard concentration.
#' @export
quantify_by_calibration <- function(data = NULL,
                                     overwrite_metadata = FALSE,
                                     fit_method = c("linear", "quadratic"),
                                     fit_weighting = c("none", "1/x", "1/x^2"),
                                     error_failed_calibration = TRUE,
                                     fail_missing_annotation = TRUE
                                     ) {

  check_data(data)

  rlang::arg_match(fit_method, c("linear", "quadratic"))
  rlang::arg_match(fit_weighting, c("none", "1/x", "1/x^2"))


   data <- calc_calibration_results(data = data,
                                   overwrite_metadata = overwrite_metadata,
                                   fit_method = fit_method,
                                   fit_weighting = fit_weighting,
                                   fail_missing_annotation = fail_missing_annotation)
  d_calib <- data@metrics_calibration

  features_no_calib <- setdiff(d_calib$feature_id, get_featurelist(data, is_istd = FALSE, is_quantifier = TRUE))
  # Check if calibration curve data is missing for any feature
  if (length(features_no_calib) > 0){
    if (!fail_missing_annotation) {
      cli::cli_alert_warning(cli::col_yellow("Calibration curve results for {length(feartures_no_calib)} features missing. Calculated concentrations of affected features will be `NA`."))
    } else {
      cli::cli_abort(cli::col_red("Calibration curve results for {length(feartures_no_calib)} features missing. Please update metadata or set `fail_missing_annotation = FALSE`."))
    }
  }

  features_failed_calib <- sum(d_calib$reg_failed)
  # Check if calibration curve data is missing for any feature
  if (features_failed_calib > 0){
    if(!fail_missing_annotation) {
      cli::cli_alert_warning(cli::col_yellow("Calibration curve fitting failed for {length(features_failed_calib)} features. Calculated concentrations of affected features will be `NA`."))
    } else {
      cli::cli_abort(cli::col_red("Calibration curve fitting failed for {length(features_failed_calib)} features. Please inspect calibration curve details, e.g. by plotting using `plot_calibrationcurves()`."))
    }
  }

  d_stats_calc <- d_calib |>
    dplyr::select("feature_id", "fit_model", "coef_a_cal_1", "coef_b_cal_1", "coef_c_cal_1")

  data@dataset <- data@dataset |>
    left_join(d_stats_calc |> select("feature_id", "fit_model", "coef_a_cal_1", "coef_b_cal_1", "coef_c_cal_1"), by = c("feature_id" = "feature_id")) |>
    mutate(feature_conc = case_when(
      fit_model != "quadratic" ~ (.data$feature_norm_intensity - .data$coef_b_cal_1) / .data$coef_a_cal_1,
      TRUE ~ purrr::pmap_dbl(
        list(coef_c_cal_1, coef_b_cal_1, coef_a_cal_1, feature_norm_intensity),
        function(c, b, a, x) {
          # Check for invalid coefficients that would lead to a degenerate polynomial
          if (is.na(c) || is.na(b) || is.na(a) || is.na(x) || a == 0) {
            return(NA_real_)  # Return NA if coefficients are invalid
          }

          # extract the real root
          root <- tryCatch({
            polyroot(c(c - x, b, a))
          }, error = function(e) {
            # Catch errors in polyroot and return NA in case of error
            return(NA_complex_)
          })

          # If the root is a complex number, check if it is valid (not NA)
          if (length(root) > 0 && !is.na(Re(root[1]))) {
            return(Re(root[1]))  # Return the real part of the first root
          } else {
            return(NA_real_)  # Return NA if no valid root
          }
        }
      )
    )) |> select(-c("coef_a_cal_1", "coef_b_cal_1", "coef_c_cal_1"))
  data
}





#' Calculate external calibration curve results
#'
#' Calibration curves are calculated for each feature using ISTD-normalized intensities and the corresponding concentrations of calibration samples, as defined in the `qc_concentrations` metadata.
#' The regression fit model (linear or quadratic) and the weighting method (either "none", "1/x", or "1/x^2") can be defined globally via
#' the arguments `fit_method` and `fit_weighting` for all features, if `overwrite_metadata` is `TRUE`.
#' Alternatively, the model and weighting can be defined individually for each feature in the `feature` metadata (columns `curve_fit_method` and `fit_weighting`).
#' If these details are missing in the metadata, the default values provided via `fit_method` and `fit_weighting` will be used.
#'
#' Additionally, the limit of detection (LoD) and limit of quantification (LoQ) are calculated for each feature based on the calibration curve.
#' LoD is calculated as 3 times the sample standard error of the regression residuals divided by the regression slope,
#' and LoQ is 10 times the same ratio. In the case of a quadratic fit, LoD and LoQ are calculated using the slope at 2x the concentration of lowest calibration point.
#'
#' The results of the regression and the calculated LoD and LoQ values are stored in the `metrics_calibration` table of the returned `MidarExperiment` object.
#'
#' @param data A `MidarExperiment` object containing the data to be used for calibration.
#' @param overwrite_metadata A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will ignore any fit method and weighting settings defined in the metadata and use the provided `fit_method` and `fit_weighting` values for all analytes.
#' @param fit_method A character string specifying the default regression fit method to use for the calibration curve. Must be one of `"linear"` or `"quadratic"`. This method will be applied if no specific fit method is defined for a feature in the metadata.
#' @param fit_weighting A character string specifying the default weighting method for the regression points in the calibration curve. Must be one of `"none"`, `"1/x"`, or `"1/x^2"`. This method will be applied if no specific weighting method is defined for a feature in the metadata.
#' @param fail_missing_annotation If `TRUE`, an error will be raised if any of the following information is missing: calibration curve data, ISTD mix volume, and sample amounts for any feature.
#' @param include_fit_object If `TRUE`, the function will return the full regression fit objects for each feature in the `metrics_calibration` table.
#'
#' @return A modified `MidarExperiment` object with an updated `metrics_calibration` table containing the calibration curve results, including concentrations, LoD, and LoQ values for each feature.
#'
#' @seealso [quantify_by_calibration()] for calculating concentrations based on external calibration curves.
#'
#' @export

#TODO: implement error handling  fail_missing_annotation
calc_calibration_results <- function(data = NULL,
                                    overwrite_metadata = FALSE,
                                    fit_method = c("linear", "quadratic"),
                                    fit_weighting = c("none", "1/x", "1/x^2"),
                                    fail_missing_annotation = TRUE,
                                    include_fit_object = FALSE
                                    ) {

  check_data(data)

  rlang::arg_match(fit_method, c("linear", "quadratic"))
  rlang::arg_match(fit_weighting, c("none", "1/x", "1/x^2"))

  if (!(c("feature_norm_intensity") %in% names(data@dataset))) cli::cli_abort("Data needs to be ISTD-normalized, please run 'normalize_by_istd' before.")
  if ( !any("CAL" %in% data@dataset$qc_type) & nrow(data@annot_qcconcentrations) == 0) cli::cli_abort("Calibration curve data missing...please check data and correct annotation om `analyis` and `qc_concentration` metadata. See this function's documentation.")


  calc_lm <- function(dt){
    tryCatch(
      {
        dt <- dt |>
          mutate(weight = case_when(
            fit_weighting[1] == "none" ~ 1,
            fit_weighting[1] == "1/x" ~ 1 / concentration,
            fit_weighting[1] == "1/x^2" ~ 1 / concentration^2,
            fit_weighting[1] == "1/sqrt(x)" ~ 1 / sqrt(concentration),
            TRUE ~ NA_real_  # Default case to handle any unexpected values
          ))

        formula <- ifelse(dt$fit_method[1] == "linear",
                                  "feature_norm_intensity ~ concentration",
                                  "feature_norm_intensity ~ poly(concentration, 2, raw = TRUE)")

        res <- lm(formula = formula, weights = dt$weight, data = dt, na.action = na.exclude)

        r.squared <- summary(res)$r.squared
        sigma <- summary(res)$sigma

        if(dt$fit_method[1] == "quadratic"){
            reg_failed <- is.na(res$coefficients[[3]]) | is.na(res$coefficients[[2]]) | is.na(res$coefficients[[1]])
            result <- list(
              feature_id = dt$feature_id[1],
              curve_id = dt$curve_id[1],
              fit_model = dt$fit_method[1],
              weighting = dt$fit_weighting[1],
              lowest_cal = sort(dt$concentration[dt$concentration != 0])[1],
              r.squared = r.squared,
              coef_a = res$coefficients[[3]],
              coef_b = res$coefficients[[2]],
              coef_c = res$coefficients[[1]],
              sigma = sigma,
              reg_failed = reg_failed,
              fit = if (include_fit_object) list(res) else list(NULL)
            )
            return(result)
          } else {
              reg_failed <- is.na(res$coefficients[[2]]) | is.na(res$coefficients[[1]])
              result <- list(
                feature_id = dt$feature_id[1],
                curve_id = dt$curve_id[1],
                fit_model = dt$fit_method[1],
                weighting = dt$fit_weighting[1],
                lowest_cal = sort(dt$concentration[dt$concentration != 0])[1],
                r.squared = r.squared,
                coef_a = res$coefficients[[2]],
                coef_b = res$coefficients[[1]],
                coef_c = NA_real_,
                sigma = sigma,
                reg_failed = reg_failed,
                fit = if (include_fit_object) list(res) else list(NULL)
              )
              return(result)
            }
      },
      error = function(e) {
        if(dt$fit_method[1] == "quadratic"){
          return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], fit_model =  dt$fit_method[1], weighting = dt$fit_weighting[1], lowest_cal = sort(dt$concentration[dt$concentration != 0])[1], r.squared = NA_real_ , coef_a = NA_real_, coef_b = NA_real_, coef_c = res$coefficients[[1]], sigma = NA_real_, reg_failed = TRUE, fit = list(NULL)))
        } else {
          return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], fit_model =  dt$fit_method[1], weighting = dt$fit_weighting[1], lowest_cal = sort(dt$concentration[dt$concentration != 0])[1], r.squared = NA_real_ , coef_a = NA_real_, coef_b = NA_real_, coef_c = NA_real_, sigma = NA_real_, reg_failed = TRUE, fit = list(NULL)))
        }
      }
    )
  }

  # mult_lowest_calib refert to multiplication factor of the lowest calibration
  # point used when calculate LoD and LoQ with a quadratic model
  # TODO: add reference
  add_quantlimits <- function(data, mult_lowest_calib = 2) {
    data <- data |>
      mutate(slope_at_conc = if_else(.data$fit_model == "quadratic",
                                     .data$coef_a + 2 * .data$coef_b * mult_lowest_calib,
                                     .data$coef_a))

    data <- data |>
      mutate(
        lod = 3 * .data$sigma / .data$slope_at_conc,
        loq = 10 * .data$sigma / .data$slope_at_conc
      )

    data
  }

  d_calib <- data@dataset |>
    dplyr::ungroup() |>
    dplyr::select(tidyselect::any_of(
      c("analysis_id", "sample_id", "qc_type", "feature_id", "is_istd", "is_quantifier", "feature_norm_intensity")
    )) |>
    filter(str_detect(.data$qc_type, "^CAL$"), !.data$is_istd, .data$is_quantifier) |>
    dplyr::left_join(data@annot_qcconcentrations, by = c("sample_id" = "sample_id", "feature_id" = "feature_id")) |>
    mutate(curve_id = "1")

  if (!overwrite_metadata) {
    d_calib <- d_calib |>
      dplyr::left_join(data@annot_features |> select("feature_id", "curve_fit_method", "curve_fit_weighting"), by = c("feature_id" = "feature_id")) |>
      mutate(fit_method = if_else(is.na(.data$curve_fit_method), .data$fit_method, .data$curve_fit_method),
             fit_weighting = if_else(is.na(.data$curve_fit_weighting), .data$fit_weighting, .data$curve_fit_weighting)) |>
      select(-"curve_fit_method", -"curve_fit_weighting")
  } else {
    d_calib <- d_calib |>
      mutate(fit_method = .data$fit_method,
             fit_weighting = .data$fit_weighting)

  }


    d_calib <- d_calib |>
      dplyr::group_split(.data$feature_id, .data$curve_id)

  d_stats <- map(d_calib, function(x) calc_lm(x)) |> bind_rows()
  d_stats <- add_quantlimits(d_stats, mult_lowest_calib = 2)

 d_stats <- d_stats |>
    dplyr::select("feature_id", "fit_model", "weighting", "lowest_cal", "curve_id", "reg_failed", r2 = "r.squared", "coef_a", "coef_b", "coef_c", "sigma", "lowest_cal", "lod", "loq", "fit") |>
    tidyr::pivot_wider(names_from = "curve_id", values_from = c("reg_failed", "r2", "coef_a", "coef_b", "coef_c", "sigma", "lowest_cal", "lod", "loq", "fit"), names_prefix = "cal_")
  data@metrics_calibration <- d_stats
  data
  }
