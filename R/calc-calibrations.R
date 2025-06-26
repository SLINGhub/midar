#' Calculate concentrations based on external calibration
#'
#' Concentrations of all features in all analyses are determined using ISTD-normalized intensities and corresponding external calibration curves.
#' Calibration curves are calculated for each feature based on calibration sample concentrations defined in the `qc_concentrations` metadata.
#' The regression fit model (linear or quadratic) and the weighting method (either "none", "1/x", or "1/x^2") can be defined globally via
#' the arguments `fit_model` and `fit_weighting` for all features, if `fit_overwrite` is `TRUE`.
#' Alternatively, the model and weighting can be defined individually for each feature in the `feature` metadata (columns `curve_fit_model` and `fit_weighting`).
#' If these details are missing in the metadata, the default values provided via `fit_model` and `fit_weighting` will be used.
#'
#' The concentrations are added to the `dataset` table as `feature_conc` column. The results of the regression and the calculated LoD and LoQ values are stored in the `metrics_calibration` table of the returned `MidarExperiment` object.

#' @param data A `MidarExperiment` object
#' @param include_qualifier A logical value. If `TRUE`, the function will include quantifier features in the calibration curve calculations.
#' @param fit_overwrite If `TRUE`,
#'   the function will use the provided `fit_model` and `fit_weighting` values
#'   for all analytes and ignore any fit method and weighting settings defined in
#'   the metadata.
#' @param fit_model A character string specifying the default regression fit
#'   method to use for the calibration curve. Must be one of `"linear"` or
#'   `"quadratic"`. This method will be applied if no specific fit method is
#'   defined for a feature in the metadata, or
#'   when `fit_overwrite = TRUE`.
#' @param fit_weighting A character string specifying the default weighting
#'   method for the regression points in the calibration curve. Must be one of
#'   `"none"`, `"1/x"`, or `"1/x^2"`. This method will be applied if no
#'   specific weighting method is defined for a feature in the metadata, or
#'   when `fit_overwrite = TRUE`.
#'@param ignore_missing_annotation If `FALSE`, raises error if any of the following information is missing: calibration curve data, ISTD mix volume and sample amounts for any feature.
#'   If `TRUE`, missing annotations will be ignored, and resulting feature concentration will be `NA`
#'@param ignore_failed_calibration If `FALSE`, raises error if calibration curve fit fails for any feature. If `TRUE`, failed fits will be ignored, and resulting feature concentration will be `NA`.
#' @return A modified `MidarExperiment` object with updated concentration values.
#'
#' @seealso [calc_calibration_results()] for calculating the calibration curve results including LoD and LoQ.
#' @seealso [quantify_by_istd()] for calculation of concentrations based on spiked-in internal standard concentration.
#' @export
quantify_by_calibration <- function(data = NULL,
                                    include_qualifier = TRUE,
                                     fit_overwrite,
                                     fit_model = c("linear", "quadratic"),
                                     fit_weighting = c("none", "1/x", "1/x^2"),
                                     ignore_failed_calibration = FALSE,
                                     ignore_missing_annotation = FALSE
                                     ) {

  check_data(data)

  rlang::arg_match(fit_model, c("linear", "quadratic"))
  rlang::arg_match(fit_weighting, c("none", "1/x", "1/x^2"))


   data <- calc_calibration_results(data = data,
                                    variable = "feature_norm_intensity",
                                    include_qualifier = include_qualifier,
                                   fit_overwrite = fit_overwrite,
                                   fit_model = fit_model,
                                   fit_weighting = fit_weighting,
                                   ignore_missing_annotation = ignore_missing_annotation)
  d_calib <- data@metrics_calibration

  features_no_calib <- setdiff( get_featurelist(data, is_istd = FALSE, is_quantifier = ifelse(include_qualifier, NA,TRUE)), d_calib$feature_id)
  # Check if calibration curve data is missing for any feature
  if (length(features_no_calib) > 0){
    if (ignore_missing_annotation) {
      cli::cli_alert_warning(cli::col_yellow("Calibration curve annotations for {length(features_no_calib)} features are missing. Calculated concentrations for these features will be `NA`."))
    } else {
      cli::cli_abort(cli::col_red("Calibration curve annotations for {length(features_no_calib)} features are missing. Please verify data and QC-concentration metadata, or ignore by setting `ignore_missing_annotation = TRUE`."))
    }
  }

  features_failed_calib <- sum(d_calib$reg_failed_cal_1)
  # Check if calibration curve data is missing for any feature
  if (features_failed_calib > 0){
    if(ignore_failed_calibration) {
      cli::cli_alert_warning(cli::col_yellow("Calibration curve fit failed for {features_failed_calib} features. Calculated concentrations for these features will be `NA`."))
    } else {
      cli::cli_abort(cli::col_red("Calibration curve fit failed for {features_failed_calib} features. Please verify data and QC-concentration metadata, or ignore by setting `ignore_failed_calibration = TRUE`."))
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


  n_features <- length(unique(d_stats_calc$feature_id))

  n_features_with_conc <- data@dataset  |> filter(!.data$is_istd, !is.na(.data$feature_conc)) |> select("feature_id") |> distinct() |> nrow()

  conc_unit <- unique(data@annot_qcconcentrations$concentration_unit)
  if (length(conc_unit) > 1)
    cli::cli_abort(col_red("Multiple concentration units found in `qc_concentrations` metadata. Please verify and correct QC-concentration metadata."))

  conc_unit <- get_conc_unit(data@annot_analyses$sample_amount_unit, conc_unit)

  samples_no_amounts <- data@annot_analyses |>
    filter(.data$valid_analysis, is.na(.data$sample_amount) | is.na(.data$istd_volume)) |>
    dplyr::semi_join(data@dataset, by = c("analysis_id"))

  count_quant_pass <- sum(!d_calib$reg_failed_cal_1[d_calib$is_quantifier])
  count_qual_pass <- sum(!d_calib$reg_failed_cal_1[d_calib$is_quantifier])
  count_quant_all <- sum(d_calib$is_quantifier)
  count_qual_all <- sum(!d_calib$is_quantifier)

  text_failed <- ifelse(features_failed_calib > 0, "the other", "these")

  if (include_qualifier && any(!d_calib$is_quantifier)) {
        cli_alert_success(cli::col_green("Concentrations of {text_failed} features were calculated for {get_analysis_count(data) - nrow(samples_no_amounts)} analyses."))
  } else {
        cli_alert_success(cli::col_green("Concentrations of {text_failed} features were calculated for {get_analysis_count(data) - nrow(samples_no_amounts)} analyses."))
  }


  cli::cli_alert_info(col_green("Concentrations are given in {conc_unit}."))

  data@status_processing <- "Calibration-quantitated data"

  data <- update_after_normalization(data, TRUE)
  data <- update_after_quantitation(data, TRUE)
  data@is_filtered <- FALSE

  data@metrics_qc <- data@metrics_qc[FALSE,]


  data
}





#' Calculate external calibration curve results
#'
#' Calibration curves are calculated for each feature using ISTD-normalized
#' intensities and the corresponding concentrations of calibration samples, as
#' defined in the `qc_concentrations` metadata. The regression fit model (linear
#' or quadratic) and the weighting method (either "none", "1/x", or "1/x^2")
#' can be defined globally via the arguments `fit_model` and `fit_weighting`
#' for all features, if `fit_overwrite` is `TRUE`. Alternatively, the
#' model and weighting can be defined individually for each feature in the
#' `feature` metadata (columns `curve_fit_model` and `fit_weighting`). If
#' these details are missing in the metadata, the default values provided via
#' `fit_model` and `fit_weighting` will be used.
#'
#' Additionally, the limit of detection (LoD) and limit of quantification (LoQ)
#' are calculated for each feature based on the calibration curve. LoD is
#' calculated as 3 times the sample standard error of the regression residuals
#' divided by the regression slope, and LoQ is 10 times the same ratio. In the
#' case of a quadratic fit, LoD and LoQ are calculated using the slope at
#' the concentration of the lowest calibration point.
#'
#' The results of the regression and the calculated LoD and LoQ values are
#' stored in the `metrics_calibration` table of the returned `MidarExperiment`
#' object.
#'
#' @param data A `MidarExperiment` object containing the data to be used for
#'   calibration.
#' @param variable A character string specifying the variable for calibration.
#'   Use `"feature_norm_intensity"` for typical scenarios involving internal
#'   standardization. When performing only external standardization, without
#'   internal standardization, use `"feature_intensity"`.
#' @param include_qualifier A logical value. If `TRUE`, the function will
#'   include quantifier features in the calibration curve calculations.
#' @param fit_overwrite If `TRUE`,
#'   the function will se the provided `fit_model` and `fit_weighting` values
#'   for all analytes and ignore any fit method and weighting settings defined in
#'   the metadata.
#' @param fit_model A character string specifying the default regression fit
#'   method to use for the calibration curve. Must be one of `"linear"` or
#'   `"quadratic"`. This method will be applied if no specific fit method is
#'   defined for a feature in the metadata, or
#'   when `fit_overwrite = TRUE`.
#' @param fit_weighting A character string specifying the default weighting
#'   method for the regression points in the calibration curve. Must be one of
#'   `"none"`, `"1/x"`, or `"1/x^2"`. This method will be applied if no
#'   specific weighting method is defined for a feature in the metadata, or
#'   when `fit_overwrite = TRUE`.
#' @param ignore_missing_annotation If `FALSE`, an error will be raised if any of
#'   the following information is missing: calibration curve data, ISTD mix
#'   volume, and sample amounts for any feature.
#' @param include_fit_object If `TRUE`, the function will return the full
#'   regression fit objects for each feature in the `metrics_calibration` table.
#'
#' @return A modified `MidarExperiment` object with an updated
#'   `metrics_calibration` table containing the calibration curve results,
#'   including concentrations, LoD, and LoQ values for each feature.
#'
#' @seealso [quantify_by_calibration()] for calculating concentrations based on
#'   external calibration curves.
#'
#' @export


calc_calibration_results <- function(data = NULL,
                                    variable = "feature_norm_intensity",
                                    include_qualifier = TRUE,
                                    fit_overwrite,
                                    fit_model,
                                    fit_weighting,
                                    ignore_missing_annotation = FALSE,
                                    include_fit_object = FALSE
                                    ) {

  check_data(data)

  rlang::arg_match(fit_model, c("linear", "quadratic"))
  rlang::arg_match(fit_weighting, c("none", "1/x", "1/x^2"))
  rlang::arg_match(variable, c("feature_intensity", "feature_norm_intensity"))


  if (variable == "feature_norm_intensity" && !any("feature_norm_intensity" %in% names(data@dataset))) cli::cli_abort("Data needs to be ISTD-normalized, please run 'normalize_by_istd' before.")

  if ( !any("CAL" %in% data@dataset$qc_type) & nrow(data@annot_qcconcentrations) == 0) cli::cli_abort("Calibration curve data missing...Please verify data and correct annotation om `analyis` and `qc_concentration` metadata. See this function's documentation.")



  calc_lm <- function(dt){

    tryCatch(
      {
        dt <- dt |>
          mutate(weight = case_when(
            fit_weighting[1] == "none" ~ 1,
            fit_weighting[1] == "1/x" ~ 1 / .data$concentration,
            fit_weighting[1] == "1/x^2" ~ 1 / .data$concentration^2,
            fit_weighting[1] == "1/sqrt(x)" ~ 1 / sqrt(.data$concentration),
            TRUE ~ NA_real_  # Default case to handle any unexpected values
          ))

        formula <- ifelse(dt$fit_model[1] == "linear",
                                  paste0(variable, " ~ concentration"),
                                  paste0(variable, " ~ poly(concentration, 2, raw = TRUE)"))

        res <- lm(formula = formula, weights = weight, data = dt, na.action = na.exclude)

        r.squared <- summary(res)$r.squared
        sigma <- summary(res)$sigma

        if(dt$fit_model[1] == "quadratic"){
            reg_failed <- is.na(res$coefficients[[3]]) | is.na(res$coefficients[[2]]) | is.na(res$coefficients[[1]])
            result <- list(
              feature_id = dt$feature_id[1],
              is_quantifier = dt$is_quantifier [1],
              curve_id = dt$curve_id[1],
              fit_model = dt$fit_model[1],
              fit_weighting = dt$fit_weighting[1],
              lowest_cal = sort(dt$concentration[dt$concentration != 0])[1],
              highest_cal = sort(dt$concentration[dt$concentration != 0],decreasing = TRUE)[1],
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
                is_quantifier = dt$is_quantifier [1],
                curve_id = dt$curve_id[1],
                fit_model = dt$fit_model[1],
                fit_weighting = dt$fit_weighting[1],
                lowest_cal = sort(dt$concentration[dt$concentration != 0])[1],
                highest_cal = sort(dt$concentration[dt$concentration != 0],decreasing = TRUE)[1],
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
        if(dt$fit_model[1] == "quadratic"){
          return(list(feature_id = dt$feature_id[1], is_quantifier = dt$is_quantifier[1], curve_id = dt$curve_id[1], fit_model =  dt$fit_model[1], fit_weighting = dt$fit_weighting[1], lowest_cal = sort(dt$concentration[dt$concentration != 0])[1],highest_cal = sort(dt$concentration[dt$concentration != 0], decreasing=TRUE)[1], r.squared = NA_real_ , coef_a = NA_real_, coef_b = NA_real_, coef_c = NA_real_, sigma = NA_real_, reg_failed = TRUE, fit = list(NULL)))
        } else {
          return(list(feature_id = dt$feature_id[1], is_quantifier = dt$is_quantifier[1], curve_id = dt$curve_id[1], fit_model =  dt$fit_model[1], fit_weighting = dt$fit_weighting[1], lowest_cal = sort(dt$concentration[dt$concentration != 0])[1],highest_cal = sort(dt$concentration[dt$concentration != 0], decreasing=TRUE)[1], r.squared = NA_real_ , coef_a = NA_real_, coef_b = NA_real_, coef_c = NA_real_, sigma = NA_real_, reg_failed = TRUE, fit = list(NULL)))
        }
      }
    )
  }


  # mult_lowest_calib refert to multiplication factor of the lowest calibration
  # point used when calculate LoD and LoQ with a quadratic model
  # TODO: add reference
  add_quantlimits <- function(data, mult_lowest_calib = 1) {
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
    dplyr::select(any_of(
      c("analysis_id", "sample_id", "qc_type", "feature_id", "analyte_id", "is_istd", "is_quantifier", variable)
    )) |>
    filter(.data$qc_type == "CAL", !.data$is_istd)

  if(!include_qualifier)
      d_calib <- d_calib |> dplyr::filter(.data$is_quantifier)

  d_calib <- d_calib |>
    dplyr::inner_join(data@annot_qcconcentrations, by = c("sample_id" = "sample_id", "analyte_id" = "analyte_id")) |>
    mutate(curve_id = "1")

  if (!fit_overwrite) {
    d_calib <- d_calib |>
      dplyr::left_join(data@annot_features |> select("feature_id", "curve_fit_model", "curve_fit_weighting"), by = c("feature_id" = "feature_id")) |>
      mutate(fit_model = if_else(is.na(.data$curve_fit_model), fit_model, .data$curve_fit_model),
             fit_weighting = if_else(is.na(.data$curve_fit_weighting), fit_weighting, .data$curve_fit_weighting)) |>
      select(-"curve_fit_model", -"curve_fit_weighting")
  } else {
    d_calib <- d_calib |>
      mutate(fit_model = fit_model,
             fit_weighting = fit_weighting)

  }

    d_calib <- d_calib |>
      dplyr::group_split(.data$feature_id, .data$curve_id)

  d_stats <- map(d_calib, function(x) calc_lm(x)) |> bind_rows()
  d_stats <- add_quantlimits(d_stats, mult_lowest_calib = 2)

 d_stats <- d_stats |>
    dplyr::select("feature_id", "is_quantifier", "fit_model", "fit_weighting", "lowest_cal", "curve_id", "reg_failed", r2 = "r.squared", "coef_a", "coef_b", "coef_c", "sigma", "lowest_cal", "highest_cal", "lod", "loq", "fit") |>
    tidyr::pivot_wider(names_from = "curve_id", values_from = c("reg_failed", "r2", "coef_a", "coef_b", "coef_c", "sigma", "lowest_cal", "highest_cal", "lod", "loq", "fit"), names_prefix = "cal_")
  data@metrics_calibration <- d_stats



  features_no_calib <- setdiff( get_featurelist(data, is_istd = FALSE, is_quantifier = ifelse(include_qualifier, NA,TRUE)), d_stats$feature_id)
  # Check if calibration curve data is missing for any feature
  if (length(features_no_calib) > 0){
      cli::cli_alert_warning(cli::col_yellow("Calibration curve annotations for {length(features_no_calib)} features are missing."))
  }

  text_missing <- ifelse(length(features_no_calib) > 0, "all annotated ", "all ")

  count_quant_pass <- sum(!d_stats$reg_failed_cal_1[d_stats$is_quantifier])
  count_qual_pass <- sum(!d_stats$reg_failed_cal_1[d_stats$is_quantifier])
  count_quant_all <- sum(d_stats$is_quantifier)
  count_qual_all <- sum(!d_stats$is_quantifier)

  text_total_quant <- ifelse(count_quant_pass == count_quant_all, paste0(text_missing, count_quant_pass), paste0(count_quant_pass, " (of ", count_quant_all, ")"))
  text_total_qual <- ifelse(count_qual_pass == count_qual_all, paste0("", count_qual_pass), paste0(count_qual_pass, " (of ", count_qual_all, ")"))





  if (include_qualifier && any(!d_stats$is_quantifier)) {
    if(count_quant_pass == 0){
      cli::cli_abort(cli::col_red("All calibration curve fits for quantifier features failed. Please check data, and feature/qc-concentration metadata."))
  }
    cli::cli_alert_success(
      cli::col_green(
        "Calibration curve fits calculated for {text_total_quant} quantifier and {text_total_qual} qualifier features. Average r-squared: {sprintf('%.4f', mean(d_stats$r2_cal_1[d_stats$is_quantifier], na.rm = TRUE))} and {sprintf('%.4f', mean(d_stats$r2_cal_1[!d_stats$is_quantifier], na.rm = TRUE))}."
      )
    )
  } else {
    if(count_quant_pass == 0) cli::cli_abort(cli::col_red("All calibration curve fits failed. Please check data, and feature/qc-concentration metadata."))
    cli::cli_alert_success(
      cli::col_green(
        "Calibration curve fits calculated for {text_total_quant} quantifier features. Average r-squared: {sprintf('%.4f', mean(d_stats$r2_cal_1[d_stats$is_quantifier], na.rm = TRUE))}."
      )
    )
  }
  data
}




#' Retrieve Calibration Regression Results
#'
#' This function retrieves calibration curve regression results from a `MidarExperiment` object.
#' It returns a summary of quality control (QC) metrics for specified QC samples.
#' including bias, absolute bias, and intra-assay coefficient of variation (CV).
#' The standard deviation of bias and percentage bias are also included unless the
#' it is `NA` for all analytes, i.e. when no replicates were measured.
#'
#' The standard deviation of concentration is also included unless the number of replicates was 1.
#'
#' @param data A `MidarExperiment` object containing the dataset and necessary annotations for calibration analysis.
#' @param qc_types A character vector specifying the QC types to include in the results, in addition to `CAL`. If not specified, all applicable QC types are included by default.
#' @param sample_ids A character vector specifying the sample IDs to include in the results. If not specified, all analyses regardless of their sample IDs are included by default.
#' @param wide_format Format of the output table. Must be one of `"none"`, `"features"`, or `"samples"`.
#' If `"none"`, the output is in long format. If `"features"`, the output is in wide format with features as columns.
#' If `"samples"`, the output is in wide format with samples as columns.
#' @param include_qualifier Logical. If `TRUE`, includes qualifier features in the results. Defaults to `FALSE`.
#' @param with_conc Logical. If `TRUE`, includes target and measured mean concentrations in the results. Defaults to `TRUE`.
#' @param with_conc_target Logical. If `TRUE`, includes target (know) concentration of the QC sample in the results. Defaults to `TRUE`.
#' @param with_bias Logical. If `TRUE`, includes percentage bias in the results. Defaults to `TRUE`.
#' @param with_bias_abs Logical. If `TRUE`, includes absolute bias in concentration units in the results. Defaults to `FALSE`.
#' @param with_conc_ratio Logical. If `TRUE`, includes the ratio of measured to target concentration in the results. Defaults to `TRUE`.
#' @param with_cv_intra Logical. If `TRUE`, includes intra-assay coefficient of variation (CV) for the in the results. Defaults to `TRUE`.
#'
#' @return A data frame containing the calibration results, including metrics such as bias, percentage bias, and intra-assay CV based on specified parameters.
#'
#' @details
#' The function uses data from the `MidarExperiment` object and filters it according to the specified QC types and other parameters. It then calculates summary statistics for each feature, such as bias and CV, and organizes the data into a user-specified format.
#'
#' @export
get_qc_bias_variability <- function(data,
                          qc_types = NA,
                          sample_ids = NA,
                          wide_format = "none",
                          include_qualifier = FALSE,
                          with_conc = TRUE,
                          with_conc_target = TRUE,
                          with_bias = TRUE,
                          with_bias_abs = FALSE,
                          with_conc_ratio = FALSE,
                          with_cv_intra = TRUE
                          ) {
  check_data(data)


  if (!is.character(wide_format) || length(wide_format) != 1) {
    cli::cli_abort("`wide_format` must be one of {.val none}, {.val features}, or {.val samples}.")
  }
  rlang::arg_match(wide_format, c("none", "features", "samples"))

  d_qc_summary <- data@dataset |>
    filter(!.data$is_istd) |>
    select("analysis_id", "qc_type", "sample_id", "feature_id", "is_quantifier", "analyte_id", "feature_conc") |>
    inner_join(data@annot_qcconcentrations |> select("sample_id", "analyte_id", target_concentration = "concentration"), by = c("sample_id", "analyte_id"))

  if(all(is.na(qc_types))){
    qc_types <- unique(d_qc_summary$qc_type)
  } else {
    if(length(setdiff(qc_types,
                      unique(d_qc_summary$qc_type))) > 0){
      cli::cli_abort(cli::col_red(paste("One or more selected `qc_types` are not present in the data or have no defined analyte concentrations. Please verify the analyses, feature and QC-concentration metadata, or select other `qc_types`.")))
    }
  }

  d_qc_summary <- d_qc_summary|> filter(.data$qc_type %in% qc_types)

  if(!all(is.na(sample_ids))){
    if(length(setdiff(sample_ids,
                      unique(d_qc_summary$qc_type))) > 0){
      cli::cli_abort(cli::col_red(paste("One or more selected `sample_id` are not present in the data or have no defined analyte concentrations. Please verify the analyses, feature and QC-concentration metadata, or select other `qc_types`.")))
    }
  }


  # Check if qc type and sample id are not paired resulting in no selected analyses
  if(!all(is.na(sample_ids))){
    d_qc_summary <- d_qc_summary |> filter(.data$sample_id %in% sample_ids)
    if(nrow(d_qc_summary) == 0){
      cli::cli_abort(cli::col_red(paste("No analyses with the selected `sample_id` and `qc_types` were found. Please verify the argument values, and corresponding feature metadata.")))
    }
  }

  if(!include_qualifier)
    d_qc_summary <- d_qc_summary |> filter(.data$is_quantifier)

  d_qc_summary <- d_qc_summary |>
     relocate("target_concentration", .after = "analyte_id") |>
    mutate(
      bias_abs_val = .data$feature_conc - .data$target_concentration,
      bias_val = (.data$feature_conc - .data$target_concentration) / .data$target_concentration * 100,
      conc_ratio = .data$feature_conc / .data$target_concentration,
    ) |>
    summarise(
      n = dplyr::n(),
      conc_target = mean(.data$target_concentration, na.rm = FALSE),
      conc_mean = mean(.data$feature_conc, na.rm = TRUE),
      conc_sd = sd(.data$feature_conc, na.rm = TRUE),
      cv_intra = .data$conc_sd / .data$conc_mean * 100,
      bias = mean(.data$bias_val, na.rm = TRUE),
      bias_abs = mean(.data$bias_abs_val, na.rm = TRUE),
      conc_ratio = mean(.data$conc_ratio, na.rm = TRUE),
      conc_ratio_sd = sd(.data$conc_ratio, na.rm = TRUE),
      .by = c("sample_id", "qc_type", "feature_id"))
  d_qc_summary <- d_qc_summary |>
  select(
      "feature_id",
      "qc_type",
      "sample_id",
      "n",
      if (with_conc_target) "conc_target",
      if (with_conc) "conc_mean",
      if (with_conc && !all(is.na(d_qc_summary$conc_sd))) "conc_sd",
      if (with_cv_intra) "cv_intra",
      if (with_bias) "bias",
      if (with_bias_abs) "bias_abs",
      if (with_conc_ratio) "conc_ratio",
      if (with_conc_ratio  && !all(is.na(d_qc_summary$conc_ratio_sd))) "conc_ratio_sd")

  if(wide_format != "none"){
    optinal_columns <- c("n", "conc_target", "conc_mean", "conc_sd","cv_intra","bias", "bias_abs", "conc_ratio", "conc_ratio_sd")
    available_columns <- intersect(optinal_columns, names(d_qc_summary))

    if(wide_format == "features"){
        d_qc_summary <- d_qc_summary |>
          tidyr::pivot_wider(names_from = "feature_id", values_from = any_of(available_columns),
                             names_sort = TRUE, names_glue = "{feature_id}_{.value}")
        d_qc_summary <- d_qc_summary |>
          select(order(colnames(d_qc_summary))) |>
          relocate("sample_id", "qc_type", .before  = 1)
      } else {
        d_qc_summary <- d_qc_summary |>
          select(-"qc_type") |>
          tidyr::pivot_wider(names_from = "sample_id", values_from = any_of(available_columns),
                             names_sort = TRUE, names_glue = "{sample_id}_{.value}")
        d_qc_summary <- d_qc_summary |>
          select(order(colnames(d_qc_summary))) |>
          relocate("feature_id", .before  = 1) |>
          arrange(.data$feature_id)
      }


  } else
  {
    d_qc_summary <- d_qc_summary |>
      relocate("feature_id", "sample_id", "qc_type", .before  = 1) |>
      arrange(.data$feature_id)
  }




}

#' Get Calibration Metrics
#'
#' Extracts calibration fit metrics from a `MidarExperiment` object.
#'
#' Requires prior computation of regression results using [`calc_calibration_results()`].
#' See its documentation for details.
#'
#' ## Returned Details and Metrics
#' - `feature_id`: Feature identifier.
#' - `is_quantifier`: Logical, indicates if the feature is a quantifier.
#' - `fit_model`: Regression model used for fitting.
#' - `weighting`: Weighting method used in fitting.
#' - `lowest_cal`: Lowest nonzero calibration concentration.
#' - `highest_cal`: Highest  calibration concentration.
#' - `r.squared`: R-squared value, indicating goodness of fit.
#' - `coef_a`: Slope of the regression line (**linear**) or coefficient of the quadratic term (`x^2`) (**quadratic**).
#' - `coef_b`: Intercept of the regression line (**linear**) or coefficient of the linear term (`x`) (**quadratic**).
#' - `coef_c`: Intercept of the regression equation (**quadratic**). Set to `NA` for **linear** models.
#' - `sigma`: Residual standard error of the model.
#' - `reg_failed`: `TRUE` if regression fitting failed.
#' - `LoD` = 3× the sample standard error of residuals / slope of the regression.
#' - `LoQ` = 10× the sample standard error of residuals / slope of the regression.
#'
#' **Note:** For LoD/LoQ calculations, the slope used in the formula is calculated
#' at the lowest nonzero calibration point for **quadratic** fits.

#'
#' @param data A `MidarExperiment` object with QC metrics.
#' @param with_lod Whether to include LoD in output. Default is `TRUE`.
#' @param with_loq Whether to include LoQ in output. Default is `TRUE`.
#' @param with_bias Whether to include bias in output. Default is `TRUE`.
#' @param with_coefficients Whether to include regression coefficients. Default is `TRUE`.
#' @param with_sigma Whether to include sigma in output. Default is `TRUE`.
#' @return A tibble with exported calibration metrics.
#' @export

get_calibration_metrics <- function(
    data = NULL,
    with_lod = TRUE,
    with_loq = TRUE,
    with_bias = TRUE,
    with_coefficients = TRUE,
    with_sigma = TRUE) {

  check_data(data)

  # Verify that the QC metrics have been calculated
  if (nrow(data@metrics_calibration) == 0) {
    cli::cli_abort(col_red("Calibration metrics has not yet been calculated. Please run `calc_calibration_results()` first."))
  }

  cal <- data@metrics_calibration |>
    dplyr::rename_with(~ str_replace(., "_cal_1", "")) |>  # Remove _cal_1 from all column names
    select(
      "feature_id",
      "is_quantifier",
      "fit_model",
      "fit_weighting",
      "reg_failed",
      "r2",
      "lowest_cal",
      "highest_cal",
      if (with_coefficients) "coef_a",
      if (with_coefficients) "coef_b",
      if (with_coefficients) "coef_c",
      if (with_lod) "lod",
      if (with_loq) "loq",
      if (with_sigma) "sigma"
    )

  # Return the metrics
  cal
}



