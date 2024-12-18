#' Calculate concentrations based on external calibration
#'

#' Concentrations are determined using ISTD-normalized intensities and external calibration curves. The determined concentration unit is based on values defined in the column `concentration_unit` of the qc_concentration metadata.
#'
#' The fit model (e.g., linear or quadratic) and the weighting method (e.g., none, 1/x, 1/x^2) used for the calibration curve can either be defined for each analyte in the metadata, or, if missing, the default values provided as arguments will be used.
#'
#' @param data A `MidarExperiment` object
#' @param overwrite_metadata_fit A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will ignore any fit method and weighting settings defined in the metadata and use the `fit_method` and `fit_weighting` values for all analytes.
#' @param fit_method A character string indicating the default regression fit method to use for the calibration curve. Must be one of `"linear"` or `"quadratic"`. This method will be used if no specific fit method is defined for a feature in the metadata.
#' @param fit_weighting A character string indicating the default weighting method for the regression points in the calibration curve. Must be one of `"none"`, `"1/x"`, or `"1/x^2"`. If no specific weighting method is defined for a feature in the metadata, this method will be used.
#' @param missing_error A logical value (`TRUE` or `FALSE`). If `TRUE`, an error will be raised if one or more ISTD concentrations and sample/ISTD amounts are missing. Default is `TRUE`.
#' @param ignore_istd A logical value (`TRUE` or `FALSE`). If `TRUE`, ISTD values with missing concentrations that are not used in any feature quantification will be ignored. Default is `FALSE`.
#'
#' @return A modified `MidarExperiment` object with updated concentration values for the analytes, based on the calibration curve calculations.
#'
#' @export
quantify_by_calibration <- function(data = NULL,
                                     overwrite_metadata_fit = FALSE,
                                     fit_method = c("linear", "quadratic"),
                                     fit_weighting = c(NA, "none", "1/x", "1/x^2"),
                                     missing_error = TRUE,
                                     ignore_istd = FALSE) {

  check_data(data)

  rlang::arg_match(fit_method, c(NA, "linear", "quadratic"))
  rlang::arg_match(fit_weighting, c(NA, "none", "1/x", "1/x^2"))

  data <- calc_calibration_results(data = data,
                                   overwrite_metadata_fit = overwrite_metadata_fit,
                                   fit_method = fit_method,
                                   fit_weighting = fit_weighting,
                                   missing_error = missing_error,
                                   ignore_istd = ignore_istd)

  d_stats_calc <- data@metrics_calibration |>
    dplyr::select("feature_id", "fit_model", "coef_a_cal_1", "coef_b_cal_1", "coef_c_cal_1")

  data@dataset <- data@dataset |>
    left_join(d_stats_calc, by = c("feature_id" = "feature_id")) |>
    mutate(feature_conc = case_when(
      fit_model != "quadratic" ~ (.data$feature_norm_intensity - .data$coef_b_cal_1) / .data$coef_a_cal_1,
      TRUE ~ purrr::pmap_dbl(
        list(coef_c_cal_1, coef_b_cal_1, coef_a_cal_1, feature_norm_intensity),
        function(c, b, a, x) {
          # Check for invalid coefficients that would lead to a degenerate polynomial
          if (is.na(c) || is.na(b) || is.na(a) || is.na(x) || a == 0) {
            return(NA_real_)  # Return NA if coefficients are invalid
          }

          # Apply polyroot and extract the real root
          root <- tryCatch({
            polyroot(c(c - x, b, a))
          }, error = function(e) {
            # Catch any errors in polyroot and return NA in case of error
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

#' Concentrations are determined using ISTD-normalized intensities and external calibration curves. The determined concentration unit is based on values defined in the column `concentration_unit` of the qc_concentration metadata.
#'
#' The fit model (e.g., linear or quadratic) and the weighting method (e.g., none, 1/x, 1/x^2) used for the calibration curve can either be defined for each analyte in the metadata, or, if missing, the default values provided as arguments will be used.
#'
#' @param data A `MidarExperiment` object
#' @param overwrite_metadata_fit A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will ignore any fit method and weighting settings defined in the metadata and use the `fit_method` and `fit_weighting` values for all analytes.
#' @param fit_method A character string indicating the default regression fit method to use for the calibration curve. Must be one of `"linear"` or `"quadratic"`. This method will be used if no specific fit method is defined for a feature in the metadata.
#' @param fit_weighting A character string indicating the default weighting method for the regression points in the calibration curve. Must be one of `"none"`, `"1/x"`, or `"1/x^2"`. If no specific weighting method is defined for a feature in the metadata, this method will be used.
#' @param missing_error A logical value (`TRUE` or `FALSE`). If `TRUE`, an error will be raised if one or more ISTD concentrations and sample/ISTD amounts are missing. Default is `TRUE`.
#' @param ignore_istd A logical value (`TRUE` or `FALSE`). If `TRUE`, ISTD values with missing concentrations that are not used in any feature quantification will be ignored. Default is `FALSE`.
#'
#' @return A modified `MidarExperiment` object with updated concentration values for the analytes, based on the calibration curve calculations.
#'
#' @export
calc_calibration_results <- function(data = NULL,
                                    overwrite_metadata_fit = FALSE,
                                    fit_method = c("linear", "quadratic"),
                                    fit_weighting = c(NA, "none", "1/x", "1/x^2"),
                                    missing_error = TRUE,
                                    ignore_istd = FALSE) {

  check_data(data)

  rlang::arg_match(fit_method, c(NA, "linear", "quadratic"))
  rlang::arg_match(fit_weighting, c(NA, "none", "1/x", "1/x^2"))



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

        res <- lm(formula = formula, weights = weight, data = dt, na.action = na.exclude)

        r.squared <- summary(res)$r.squared
        sigma <- summary(res)$sigma
        if(dt$fit_method[1] == "quadratic"){
          return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], fit_model =  dt$fit_method[1], weighting = dt$fit_weighting[1], lowest_cal = sort(dt$concentration[dt$concentration != 0])[1], r.squared = r.squared , coef_a = res$coefficients[[3]], coef_b = res$coefficients[[2]], coef_c = res$coefficients[[1]], sigma = sigma))
        } else {
          return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], fit_model =  dt$fit_method[1], weighting = dt$fit_weighting[1], lowest_cal = sort(dt$concentration[dt$concentration != 0])[1], r.squared = r.squared , coef_a = res$coefficients[[2]], coef_b = res$coefficients[1], coef_c = NA_real_, sigma = sigma))
        }
      },
      error = function(e) {
        if(dt$fit_method[1] == "quadratic"){
          return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], fit_model =  dt$fit_method[1], weighting = dt$fit_weighting[1], lowest_cal = sort(dt$concentration[dt$concentration != 0])[1], r.squared = NA_real_ , coef_a = NA_real_, coef_b = NA_real_, coef_c = res$coefficients[[1]], sigma = NA_real_))
        } else {
          return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], fit_model =  dt$fit_method[1], weighting = dt$fit_weighting[1], lowest_cal = sort(dt$concentration[dt$concentration != 0])[1], r.squared = NA_real_ , coef_a = NA_real_, coef_b = NA_real_, coef_c = NA_real_, sigma = NA_real_))
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
        lod = 3 * .data$sigma / slope_at_conc,
        loq = 10 * .data$sigma / slope_at_conc
      )

    data
  }


  d_calib <- data@dataset |>
    dplyr::ungroup() |>
    dplyr::select(tidyselect::any_of(
      c("analysis_id", "sample_id", "qc_type", "feature_id", "is_istd", "is_quantifier", "feature_norm_intensity")
    )) |>
    filter(str_detect(qc_type, "^CAL$"), !is_istd, is_quantifier) |>
    dplyr::left_join(data@annot_qcconcentrations, by = c("sample_id" = "sample_id", "feature_id" = "feature_id")) |>
    mutate(curve_id = "1") |>
    mutate(fit_method = fit_method,
           fit_weighting = fit_weighting)



    d_calib <- d_calib |>
      dplyr::group_split(.data$feature_id, .data$curve_id)

  d_stats <- map(d_calib, function(x) calc_lm(x)) |> bind_rows()


  d_stats <- add_quantlimits(d_stats, mult_lowest_calib = 3)

  d_stats <- d_stats |>
    dplyr::select("feature_id", "fit_model", "weighting", "lowest_cal", "curve_id", r2 = "r.squared", "coef_a", "coef_b", "coef_c", "sigma", "lowest_cal", "lod", "loq") |>
    tidyr::pivot_wider(names_from = "curve_id", values_from = c("r2", "coef_a", "coef_b", "coef_c", "sigma", "lowest_cal", "lod", "loq"), names_prefix = "cal_")
  data@metrics_calibration <- d_stats

  data

  }
