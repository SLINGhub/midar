library(testthat)
library(dplyr)

mexp <- quant_lcms_dataset
mexp_norm <- normalize_by_istd(mexp)

test_that("calc_calibration_results works", {
  expect_message(
    mexp_res <- calc_calibration_results(mexp_norm,
                                         overwrite_fit_param = TRUE,
                                         fit_model = "linear",
                                         fit_weighting = "1/x"),
    "Calibration curve fits calculated for all 4 quantifier and 4 qualifier features")

  res <- mexp_res@metrics_calibration
  expect_equal(dim(res), c(8, 15))
  expect_equal(unique(res$fit_model), "linear")
  expect_equal(unique(res$fit_weighting), "1/x")


  expect_message(
    mexp_res <- calc_calibration_results(mexp_norm,
                                         include_qualifier = FALSE,
                                         overwrite_fit_param = TRUE,
                                         fit_model = "linear",
                                         fit_weighting = "1/x"),
    "Calibration curve fits calculated for all 4 quantifier features")

  res <- mexp_res@metrics_calibration
  expect_equal(dim(res), c(4, 15))
  expect_equal(unique(res$fit_model), "linear")
  expect_equal(unique(res$fit_weighting), "1/x")


  mexp_res <- calc_calibration_results(mexp_norm,
                                       overwrite_fit_param = FALSE,
                                       fit_model = "linear",
                                       fit_weighting = "1/x")
  res <- mexp_res@metrics_calibration
  expect_equal(unique(res$fit_model), c("quadratic", "linear"))
  expect_equal(unique(res$fit_weighting), "1/x")
  expect_equal(mean(res$r2_cal_1), 0.97875279)
  expect_equal(mean(res$lowest_cal_cal_1),  2.01350)
  expect_equal(mean(res$loq_cal_1),  4.690940)

  # Missing fit parameter replaced with defauls provided with fit_ args.
  mexp_temp <- mexp_norm
  mexp_temp@annot_features$curve_fit_model[c(1,3, 5, 7)] <- NA
  mexp_temp@annot_features$curve_fit_weighting[c(1,3, 5, 7)] <- NA

  mexp_res <- calc_calibration_results(mexp_temp,
                                       overwrite_fit_param = FALSE,
                                       fit_model = "linear",
                                       fit_weighting = "none")
  res <- mexp_res@metrics_calibration
  expect_equal(unique(res$fit_model[c(1,2, 3, 4)]), c("linear"))
  expect_equal(unique(res$fit_weighting[c(1,2, 3, 4)]), c("none", "1/x"))
})

test_that("calc_calibration_results error handling works", {

  mexp_temp <- mexp_norm

  mexp_temp@annot_qcconcentrations$concentration <- NA


  expect_error(
    mexp_res <- calc_calibration_results(mexp_temp,
                                         overwrite_fit_param = TRUE,
                                         fit_model = "linear",
                                         fit_weighting = "1/x"),
    "All calibration curve fits for quantifier features")

  expect_error(
    mexp_res <- calc_calibration_results(mexp_temp,
                                         overwrite_fit_param = TRUE,
                                         include_qualifier = FALSE,
                                         fit_model = "linear",
                                         fit_weighting = "1/x"),
    "All calibration curve fits failed")

})


test_that("quantify_by_calibration works", {
  expect_message(
    mexp_res <- quantify_by_calibration(mexp_norm,
                                         overwrite_fit_param = FALSE,
                                         include_qualifier = TRUE,
                                         fit_model = "quadratic",
                                         fit_weighting = "1/x"),
    "Concentrations of these features were calculated for 25 analyses")

  res <- mexp_res@dataset |> filter (analysis_id == "CalE", !is_istd)
  expect_equal(mean(res$feature_conc, na.rm = TRUE), 100.8283531)

  res <- mexp_res@metrics_calibration
  expect_equal(unique(res$fit_model), c("quadratic", "linear"))

  expect_message(
    mexp_res <- quantify_by_calibration(mexp_norm,
                                        overwrite_fit_param = FALSE,
                                        include_qualifier = FALSE,
                                        fit_model = "quadratic",
                                        fit_weighting = "1/x"),
    "Concentrations of these features were calculated for 25 analyses")

  res <- mexp_res@dataset |> filter (analysis_id == "CalE", !is_istd)
  expect_equal(mean(res$feature_conc, na.rm = TRUE),  101.651349)

})

test_that("quantify_by_calibration handles errors", {

  mexp_temp <- mexp_norm
  mexp_temp@annot_qcconcentrations <- mexp_temp@annot_qcconcentrations |>
    mutate(concentration = if_else(str_detect(analyte_id, "Cortiso"), NA_real_, concentration))

  expect_error(
    mexp_res <- quantify_by_calibration(mexp_temp,
                                        overwrite_fit_param = FALSE,
                                        include_qualifier = FALSE,
                                        ignore_failed_calibration = FALSE,
                                        fit_model = "quadratic",
                                        fit_weighting = "1/x"),
  "Calibration curve fit failed for 2 features")

  expect_message(
    mexp_res <- quantify_by_calibration(mexp_temp,
                                        overwrite_fit_param = FALSE,
                                        include_qualifier = FALSE,
                                        ignore_failed_calibration = TRUE,
                                        fit_model = "quadratic",
                                        fit_weighting = "1/x"),
    "Calibration curve fit failed for 2 features")

  expect_message(
    mexp_res <- quantify_by_calibration(mexp_temp,
                                        overwrite_fit_param = FALSE,
                                        include_qualifier = FALSE,
                                        ignore_failed_calibration = TRUE,
                                        fit_model = "quadratic",
                                        fit_weighting = "1/x"),
  "Concentrations of the other features were calculated")

  mexp_temp <- mexp_norm
  mexp_temp@annot_qcconcentrations <- mexp_temp@annot_qcconcentrations |> filter(!str_detect(analyte_id, "Cortiso"))

  expect_error(
    mexp_res <- quantify_by_calibration(mexp_temp,
                                        overwrite_fit_param = FALSE,
                                        include_qualifier = FALSE,
                                        ignore_failed_calibration = FALSE,
                                        fit_model = "quadratic",
                                        fit_weighting = "1/x"),
    "Calibration curve annotations for 2 features are missing.")

  expect_message(
    mexp_res <- quantify_by_calibration(mexp_temp,
                                        overwrite_fit_param = FALSE,
                                        include_qualifier = FALSE,
                                        ignore_failed_calibration = FALSE,
                                        ignore_missing_annotation = TRUE,
                                        fit_model = "quadratic",
                                        fit_weighting = "1/x"),
    "Calibration curve annotations for 2 features are missing.")
})







mexp_quant <- quant_lcms_dataset
mexp_quant_norm <- normalize_by_istd(mexp_quant)
mexp_quant_norm <- calc_calibration_results(mexp_quant_norm,fit_model = "quadratic",fit_weighting = "1/x")
mexp_quant_norm <- quantify_by_calibration(mexp_quant_norm,
                                           fit_model = "quadratic",fit_weighting = "1/x")


 test_that("get_qc_bias_variability returns correct data", {
   result <- get_qc_bias_variability(mexp_quant_norm, qc_types = c("CAL", "LQC", "HQC"))
   expect_s3_class(result, "data.frame")
   expect_equal(names(result), c("feature_id", "sample_id", "qc_type", "conc_target", "conc_mean", "conc_sd", "bias","bias_sd","bias_perc","bias_perc_sd","cv_intra"))
   expect_equal(nrow(result), 32)

   result <- get_qc_bias_variability(mexp_quant_norm, qc_types = NA,
                                          with_conc = FALSE,
                                          with_bias = FALSE,
                                          with_bias_perc = FALSE,
                                          with_cv_intra = FALSE)
   expect_equal(names(result), c("feature_id","sample_id","qc_type"))
   expect_equal(unique(result$qc_type), c("CAL","HQC","LQC"))

   result <- get_qc_bias_variability(mexp_quant_norm, include_qualifier = TRUE)
   expect_equal(nrow(result), 64)

   result <- get_qc_bias_variability(mexp_quant_norm, wide_format = "features")
   expect_equal(names(result)[1:3], c("sample_id","qc_type","Aldosterone_bias"))
   expect_equal(nrow(result), 8)

   result <- get_qc_bias_variability(mexp_quant_norm, wide_format = "samples")
   expect_equal(names(result)[1:3], c("feature_id","CAL-A_bias","CAL-A_bias_perc"))
   expect_equal(nrow(result), 4)
   })

 test_that("get_qc_bias_variability handles errors", {
   expect_error(
     get_qc_bias_variability(mexp_quant_norm, qc_types = c("CAL", "LQC", "HQC", "EQA")),
                "One or more selected \\`qc_types\\`")

   expect_error(
     get_qc_bias_variability(mexp_quant_norm, wide_format = FALSE),
     "\\`wide_format\\` must be one of")
 })

#
#
 test_that("get_calibration_metrics returns correct data", {
  result <- get_calibration_metrics(mexp_quant_norm)

  expect_s3_class(result, "data.frame")
  expect_equal(names(result), c("feature_id","is_quantifier","fit_model","fit_weighting","reg_failed","r2","lowest_cal","highest_cal","coef_a","coef_b","coef_c","lod","loq", "sigma"))

  result <- get_calibration_metrics(mexp_quant_norm,
                                    with_lod = FALSE,
                                    with_loq = FALSE,
                                    with_bias = FALSE,
                                    with_coefficients = FALSE,
                                    with_sigma = FALSE)

  expect_equal(names(result), c("feature_id","is_quantifier","fit_model","fit_weighting","reg_failed","r2","lowest_cal","highest_cal"))

})


 test_that("get_calibration_metrics handles errors", {
   expect_error(
     get_calibration_metrics(mexp_quant),
     "Calibration metrics has not yet been calculated")
 })
