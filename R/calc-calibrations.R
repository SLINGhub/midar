#' Calculate concentrations based on external calibration
#'
#' This function calculates the concentrations of analytes in samples based on a calibration curve derived from external calibration samples.
#' Concentrations are determined using ISTD-normalized intensities, and the determined concentration unit is based on the `sample_amount_unit` provided in the analysis metadata.
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



  calc_lm <- function(dt){
    tryCatch(
      {
        res <- lm(formula = dt$formula[1], weights = weight, data = dt, na.action = na.exclude)

        r.squared <- summary(res)$r.squared
        slope <- res$coefficients[[2]]
        intercept <- res$coefficients[1]
        sigma <- summary(res)$sigma
        return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], r.squared = r.squared , slope = slope, intercept = intercept, sigma = sigma))

      },
      error = function(e) {
        return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], r.squared = NA_real_ , slope = NA_real_, intercept = NA_real_))
      }
    )
  }

  d_calib <- data@dataset |>
    dplyr::ungroup() |>
    dplyr::select(tidyselect::any_of(
      c("analysis_id", "sample_id", "qc_type", "feature_id", "is_istd", "is_quantifier", "feature_norm_intensity")
    )) |>
    filter(str_detect(qc_type, "^CAL$"), !is_istd, is_quantifier) |>
    dplyr::left_join(data@annot_qcconcentrations, by = c("sample_id" = "sample_id", "feature_id" = "feature_id")) |>
    mutate(curve_id = "1")

  d_calib <- d_calib |>
    mutate(weight = case_when(
      fit_weighting == "none" ~ 1,
      fit_weighting == "1/x" ~ 1 / concentration,
      fit_weighting == "1/x^2" ~ 1 / concentration^2,
      fit_weighting == "1/sqrt(x)" ~ 1 / sqrt(concentration),
      TRUE ~ NA_real_  # Default case to handle any unexpected values
    ))

  d_calib$formula <- ifelse(fit_method == "linear", "feature_norm_intensity ~ concentration", "feature_norm_intensity ~ concentration + I(concentration^2)")


    d_calib <- d_calib |>
      dplyr::group_split(.data$feature_id, .data$curve_id)


  d_stats <- map(d_calib, function(x) calc_lm(x))

  d_stats <- d_stats |> bind_rows() |>
    dplyr::mutate(slopenorm = .data$slope,
                  y0 = .data$intercept,
                  lod =  3 * .data$sigma / .data$slope,
                  loq =  10 * .data$sigma / .data$slope) |>
    dplyr::select("feature_id", "curve_id", r2 = "r.squared", "slope", "y0", "sigma", "lod", "loq") |>
    tidyr::pivot_wider(names_from = "curve_id", values_from = c("r2", "slope", "y0", "sigma", "lod", "loq"), names_prefix = "cal_")



  }
