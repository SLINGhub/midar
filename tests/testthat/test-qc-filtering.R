library(testthat)
library(dplyr)

mexp_orig <- lipidomics_dataset
mexp <- exclude_analyses(mexp_orig, analyses = "Longit_batch6_51", clear_existing = TRUE )
mexp <- normalize_by_istd(mexp_orig)
mexp <- quantify_by_istd(mexp)
mexp_proc <- calc_qc_metrics(mexp,  use_batch_medians = FALSE)




test_that("calc_qc_metrics works for all qc groups", {
    mexp_res <- calc_qc_metrics(mexp,
                                use_batch_medians = FALSE)

    expect_s4_class(mexp_res, "MidarExperiment")
    expect_equal(dim(mexp_res@metrics_qc), c(29,79))

    expect_equal(max(mexp_res@metrics_qc$product_mz), 829.4)
    expect_equal(min(mexp_res@metrics_qc$missing_intensity_prop_spl),0)
    expect_equal(sum(mexp_res@metrics_qc$in_data),29)
    expect_equal(sum(mexp_res@metrics_qc$is_quantifier),28)
    expect_equal(sum(mexp_res@metrics_qc$is_istd),9)
    expect_equal(sum(mexp_res@metrics_qc$valid_feature),29)
    expect_equal(max(mexp_res@metrics_qc$collision_energy),30)
    expect_equal(sum(mexp_res@metrics_qc$na_in_all),0)
    expect_equal(max(mexp_res@metrics_qc$rt_median_SPL),7.311)
    expect_equal(max(mexp_res@metrics_qc$intensity_median_SPL),38240524)
    expect_equal(max(mexp_res@metrics_qc$intensity_min_SPL),3678267.3)
    expect_equal(max(mexp_res@metrics_qc$intensity_max_SPL),50621380.0)
    expect_equal(max(mexp_res@metrics_qc$intensity_cv_SPL),76.06928681)
    expect_equal(max(mexp_res@metrics_qc$intensity_cv_TQC),30.1009462)
    expect_equal(max(mexp_res@metrics_qc$intensity_cv_BQC),30.29650738)
    expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_SPL),103.7784475)
    expect_equal(max(mexp_res@metrics_qc$conc_cv_SPL),103.7784475)
    expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_TQC, na.rm = T),32.9608359974)
    expect_equal(max(mexp_res@metrics_qc$conc_cv_TQC, na.rm = T),32.9608359974)
    expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_BQC, na.rm = T), 31.93429867)
    expect_equal(max(mexp_res@metrics_qc$conc_cv_BQC, na.rm = T),31.93429867)
    expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_bqc_conc, na.rm = T), 0.5136919558)
    expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_bqc_normint, na.rm = T),0.5136919558)
    expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_tqc_conc, na.rm = T),0.5077020488)
    expect_equal(median(mexp_res@metrics_qc$conc_dratio_mad_bqc_conc, na.rm = T),0.6261527021)
    expect_equal(min(mexp_res@metrics_qc$r2_rqc_A),0.91931047)
    expect_equal(min(mexp_res@metrics_qc$r2_rqc_B),0.85693787)
    expect_equal(min(mexp_res@metrics_qc$slopenorm_rqc_A),0.69281368)
    expect_equal(min(mexp_res@metrics_qc$slopenorm_rqc_B),0.65001359)
})

test_that("calc_qc_metrics batch-wise works for all qc groups", {
  mexp_res <- calc_qc_metrics(mexp,
                              use_batch_medians = TRUE)

  expect_s4_class(mexp_res, "MidarExperiment")
  expect_equal(dim(mexp_res@metrics_qc), c(29,79))

  expect_equal(max(mexp_res@metrics_qc$product_mz), 829.4)
  expect_equal(min(mexp_res@metrics_qc$missing_intensity_prop_spl),0)
  expect_equal(sum(mexp_res@metrics_qc$in_data),29)
  expect_equal(sum(mexp_res@metrics_qc$is_quantifier),28)
  expect_equal(sum(mexp_res@metrics_qc$is_istd),9)
  expect_equal(sum(mexp_res@metrics_qc$valid_feature),29)
  expect_equal(max(mexp_res@metrics_qc$collision_energy),30)
  expect_equal(sum(mexp_res@metrics_qc$na_in_all),0)
  expect_equal(max(mexp_res@metrics_qc$rt_median_SPL),7.311)
  expect_equal(max(mexp_res@metrics_qc$intensity_median_SPL),37547322)
  expect_equal(max(mexp_res@metrics_qc$intensity_min_SPL),25756142)
  expect_equal(max(mexp_res@metrics_qc$intensity_max_SPL),45840582)
  expect_equal(max(mexp_res@metrics_qc$intensity_cv_SPL), 71.0697513)
  expect_equal(max(mexp_res@metrics_qc$intensity_cv_TQC),20.4880956525)
  expect_equal(max(mexp_res@metrics_qc$intensity_cv_BQC),16.3277340356)
  expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_SPL),96.9355336955)
  expect_equal(max(mexp_res@metrics_qc$conc_cv_SPL),96.9355336955)
  expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_TQC),22.8694786877)
  expect_equal(max(mexp_res@metrics_qc$conc_cv_TQC),22.8694786877)
  expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_BQC), 16.64318769)
  expect_equal(max(mexp_res@metrics_qc$conc_cv_BQC),16.64318769)
  expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_bqc_conc, na.rm = T), 0.3984480619)
  expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_bqc_normint, na.rm = T),0.3984480619)
  expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_tqc_conc, na.rm = T),0.4245940972)
  expect_equal(median(mexp_res@metrics_qc$conc_dratio_mad_bqc_conc, na.rm = T),0.4060582407)
  expect_equal(min(mexp_res@metrics_qc$r2_rqc_A),0.91931047)
  expect_equal(min(mexp_res@metrics_qc$r2_rqc_B),0.85693787)
  expect_equal(min(mexp_res@metrics_qc$slopenorm_rqc_A),0.69281368)
  expect_equal(min(mexp_res@metrics_qc$slopenorm_rqc_B),0.65001359)
})

test_that("calc_qc_metrics batch-wise works for with calibration metrics ", {
  mexp_quant <- quant_lcms_dataset
  mexp_quant_norm <- normalize_by_istd(mexp_quant)
  mexp_quant_norm <- calc_calibration_results(mexp_quant_norm,fit_model = "quadratic",fit_weighting = "1/x")

  mexp_res <- calc_qc_metrics(mexp_quant_norm,
                              use_batch_medians = TRUE,
                              include_norm_intensity_stats = FALSE,
                              include_conc_stats = FALSE,
                              include_response_stats = FALSE,
                              include_calibration_results = TRUE)

  expect_true(all(c("fit_model", "fit_weighting", "reg_failed_cal", "r2_cal") %in% names(mexp_res@metrics_qc)))

  })


test_that("calc_qc_metrics batch-wise works for all with all incl FALSE ", {
  mexp_res <- calc_qc_metrics(mexp,
                              use_batch_medians = TRUE,
                              include_norm_intensity_stats = FALSE,
                              include_conc_stats = FALSE,
                              include_response_stats = FALSE,
                              include_calibration_results = FALSE)

  expect_s4_class(mexp_res, "MidarExperiment")
  expect_equal(dim(mexp_res@metrics_qc), c(29,50))

  expect_false("norm_intensity_cv_SPL" %in% colnames(mexp_res@metrics_qc))
  expect_false("conc_cv_SPL" %in% colnames(mexp_res@metrics_qc))
  expect_false("r2_rqc_A" %in% colnames(mexp_res@metrics_qc))
  expect_false("fit_model" %in% colnames(mexp_res@metrics_qc))
})

test_that("calc_qc_metrics batch-wise works for all with all incl FALSE across batches ", {
  mexp_res <- calc_qc_metrics(mexp,
                              use_batch_medians = FALSE,
                              include_norm_intensity_stats = FALSE,
                              include_conc_stats = FALSE,
                              include_response_stats = FALSE,
                              include_calibration_results = FALSE)

  expect_equal(dim(mexp_res@metrics_qc), c(29,50))

  expect_false("norm_intensity_cv_SPL" %in% colnames(mexp_res@metrics_qc))
  expect_false("conc_cv_SPL" %in% colnames(mexp_res@metrics_qc))
  expect_false("r2_rqc_A" %in% colnames(mexp_res@metrics_qc))
  expect_false("fit_model" %in% colnames(mexp_res@metrics_qc))
})

test_that("calc_qc_metrics batch-wise works for some incl FALSE ", {
  mexp_res <- calc_qc_metrics(mexp,
                              use_batch_medians = TRUE,
                              include_norm_intensity_stats = TRUE,
                              include_conc_stats = FALSE,
                              include_response_stats = TRUE,
                              include_calibration_results = FALSE)
  expect_equal(dim(mexp_res@metrics_qc), c(29,61))

  expect_true("norm_intensity_cv_SPL" %in% colnames(mexp_res@metrics_qc))
  expect_false("conc_cv_SPL" %in% colnames(mexp_res@metrics_qc))
  expect_true("r2_rqc_A" %in% colnames(mexp_res@metrics_qc))
  expect_false("fit_model" %in% colnames(mexp_res@metrics_qc))
})

test_that("calc_qc_metrics batch-wise works for some other incl FALSE ", {
  mexp_res <- calc_qc_metrics(mexp,
                              use_batch_medians = TRUE,
                              include_norm_intensity_stats = FALSE,
                              include_conc_stats = TRUE,
                              include_response_stats = TRUE,
                              include_calibration_results = FALSE)
  expect_equal(dim(mexp_res@metrics_qc), c(29,74))

  expect_false("norm_intensity_cv_SPL" %in% colnames(mexp_res@metrics_qc))
  expect_true("conc_cv_SPL" %in% colnames(mexp_res@metrics_qc))
  expect_true("r2_rqc_A" %in% colnames(mexp_res@metrics_qc))
  expect_false("fit_model" %in% colnames(mexp_res@metrics_qc))
})


test_that("calc_qc_metrics batch-wise works at different processing status ", {
  mexp_temp <- mexp
  mexp_temp@is_istd_normalized <- FALSE
  mexp_temp@is_quantitated <- FALSE
  # delete all rows of tibble below

  #mexp_temp@annot_responsecurves <- mexp_temp@annot_responsecurves[0,]
  mexp_res <- calc_qc_metrics(mexp_temp,
                              use_batch_medians = TRUE)
  expect_equal(dim(mexp_res@metrics_qc), c(29,56))
  expect_false("norm_intensity_cv_SPL" %in% colnames(mexp_res@metrics_qc))
  expect_false("conc_cv_SPL" %in% colnames(mexp_res@metrics_qc))
})

test_that("calc_qc_metrics no method data works", {
  mexp_temp <- mexp
  mexp_temp@dataset_orig$method_precursor_mz <- NULL
  mexp_temp@dataset_orig$method_product_mz <- NULL
  mexp_temp@dataset_orig$method_collision_energy <- NULL

  #mexp_temp@annot_responsecurves <- mexp_temp@annot_responsecurves[0,]
  mexp_res <- calc_qc_metrics(mexp_temp,
                              use_batch_medians = TRUE)
  expect_true(all(is.na(mexp_res@metrics_qc$precursor_mz)))
})

test_that("calc_qc_metrics batch-wise raise error correctly when data missing ", {
  mexp_temp <- mexp
  mexp_temp@is_istd_normalized <- FALSE
  mexp_temp@is_quantitated <- FALSE
  # delete all rows of tibble below

  #mexp_temp@annot_responsecurves <- mexp_temp@annot_responsecurves[0,]
  expect_error(calc_qc_metrics(mexp_temp,
                              use_batch_medians = TRUE,
                              include_norm_intensity_stats = TRUE,
                              include_conc_stats = TRUE,
                              include_response_stats = TRUE,
                              include_calibration_results = TRUE),
               "Normalized intensity data is missing")

  expect_error(calc_qc_metrics(mexp_temp,
                               use_batch_medians = TRUE,
                               include_norm_intensity_stats = NA,
                               include_conc_stats = TRUE,
                               include_response_stats = TRUE,
                               include_calibration_results = TRUE),
               "Concentration data is missing")


  mexp_temp@is_istd_normalized <- TRUE
  expect_error(calc_qc_metrics(mexp_temp,
                               use_batch_medians = TRUE,
                               include_norm_intensity_stats = FALSE,
                               include_conc_stats = TRUE,
                               include_response_stats = TRUE,
                               include_calibration_results = TRUE),
               "Concentration data is missing")

  expect_error(calc_qc_metrics(mexp_temp,
                               use_batch_medians = TRUE,
                               include_norm_intensity_stats = FALSE,
                               include_conc_stats = TRUE,
                               include_response_stats = TRUE,
                               include_calibration_results = TRUE),
               "Concentration data is missing")

  mexp_temp@is_istd_normalized <- FALSE
  mexp_temp@is_quantitated <- TRUE

  expect_error(calc_qc_metrics(mexp_temp,
                               use_batch_medians = TRUE,
                               include_norm_intensity_stats = TRUE,
                               include_conc_stats = TRUE,
                               include_response_stats = TRUE,
                               include_calibration_results = TRUE),
               "Normalized intensity data is missing")
})

test_that("calc_qc_metrics handles missing/missmatching info for response curve stats ", {
  mexp_temp <- mexp
  mexp_temp@annot_responsecurves <- mexp_temp@annot_responsecurves[0,]

  expect_s4_class(calc_qc_metrics(mexp_temp,
                               use_batch_medians = TRUE,
                               include_norm_intensity_stats = TRUE,
                               include_conc_stats = TRUE,
                               include_response_stats = NA,
                               include_calibration_results = NA),
               "MidarExperiment")

  expect_error(calc_qc_metrics(mexp_temp,
                              use_batch_medians = TRUE,
                              include_norm_intensity_stats = TRUE,
                              include_conc_stats = TRUE,
                              include_response_stats = TRUE,
                              include_calibration_results = FALSE),
               "No response curve metadata found")

  expect_error(calc_qc_metrics(mexp_temp,
                               use_batch_medians = TRUE,
                               include_norm_intensity_stats = TRUE,
                               include_conc_stats = TRUE,
                               include_response_stats = FALSE,
                               include_calibration_results = TRUE),
               "Calibration metrics are missing")

  mexp_temp <- mexp
  mexp_temp@annot_responsecurves$analysis_id[1] <- "unknown1"
  mexp_temp@annot_responsecurves$analysis_id[3] <- "unknown2"

  expect_s4_class(calc_qc_metrics(mexp_temp,
                                  use_batch_medians = TRUE,
                                  include_norm_intensity_stats = TRUE,
                                  include_conc_stats = TRUE,
                                  include_response_stats = NA,
                                  include_calibration_results = NA),
                  "MidarExperiment")

  expect_error(calc_qc_metrics(mexp_temp,
                               use_batch_medians = TRUE,
                               include_norm_intensity_stats = TRUE,
                               include_conc_stats = TRUE,
                               include_response_stats = TRUE,
                               include_calibration_results = NA),
               "One or more analysis IDs")


})

test_that("filter_features_qc works with istd and qualifier subsetting", {
  expect_message(
    mexp_res <- filter_features_qc(
      mexp_proc,
      clear_existing = TRUE,
      include_qualifier = FALSE,
      include_istd = FALSE,
      min.intensity.median.bqc = 0),
    "19 of 19 quantifier features meet QC criteria \\(not including the 9 quantifier ISTD features\\)"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 19)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_proc,
      clear_existing = TRUE,
      include_qualifier = TRUE,
      include_istd = FALSE,
      min.intensity.median.bqc = 0),
    "19 of 19 quantifier and 1 of 1 qualifier features meet QC criteria \\(not including the 9 quantifier and 0 qualifier ISTD features\\)"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 20)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_proc,
      clear_existing = TRUE,
      include_qualifier = FALSE,
      include_istd = TRUE,
      min.intensity.median.bqc = 0),
    "28 of 28 quantifier features meet QC criteria \\(including the 9 quantifier ISTD features\\)"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 28)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_proc,
      clear_existing = TRUE,
      include_qualifier = TRUE,
      include_istd = TRUE,
      min.intensity.median.bqc = 0),
    "28 of 28 quantifier and 1 of 1 qualifier features meet QC criteria \\(including the 9 quantifier and 0 qualifier ISTD features\\)"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 29)
})


test_that("filter_features_qc works on selected criteria", {
  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.intensity.median.bqc = 1E5)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 16)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.intensity.median.tqc = 1E5)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 16)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.intensity.median.bqc = 1E5,
    min.intensity.median.tqc = 1E5)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 16)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.signalblank.median.spl.pblk = 100)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 12)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.intensity.median.bqc = 1E5,
    min.signalblank.median.spl.pblk = 100)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 10)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    max.cv.conc.bqc = 20)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 18)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.signalblank.median.spl.pblk = 100,
    max.cv.conc.bqc = 20)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.signalblank.median.spl.pblk = 100,
    max.cv.normintensity.bqc = 20)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    max.dratio.sd.conc.bqc = 0.5)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 10)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    max.prop.missing.conc.spl = 0)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 19)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = 1,
    response.curves.summary = "mean",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 8)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = 2,
    response.curves.summary = "mean",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 9)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = c(1,2),
    response.curves.summary = "mean",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 6)


  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = c("A","B"),
    response.curves.summary = "mean",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 6)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = c(1,2),
    response.curves.summary = "median",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 6)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = c(1,2),
    response.curves.summary = "best",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 12)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = c(1,2),
    response.curves.summary = "worst",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = c(1,2),
    response.curves.summary = "mean",
    min.slope.response = 0.9)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = c(1,2),
    response.curves.summary = "worst",
    min.slope.response = 0.9)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 9)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = c(1,2),
    response.curves.summary = "best",
    max.slope.response = 1.005)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 17)

  mexp_res <- filter_features_qc(
    mexp_proc, clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    response.curves.selection = c(1,2),
    response.curves.summary = "worst",
    max.slope.response = 1.005)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 12)
})


# Confirm overwriting of QC criteria works


test_that("Confirm overwriting of QC criteria works", {
  mexp_res <- filter_features_qc(
    mexp_proc,
    clear_existing = FALSE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.signalblank.median.spl.pblk = 300)

  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 9)


  expect_message(
   mexp_res <- filter_features_qc(
    mexp_res,
    clear_existing = FALSE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.signalblank.median.spl.pblk = 100),
   "Replaced following previously defined QC filters: Signal-to-Blank"
   )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 12)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_res,
      include_qualifier = FALSE,
      include_istd = FALSE,
      clear_existing = FALSE,
      max.cv.conc.bqc = 20),
    "Feature QC filters were updated\\: 11 \\(before 12\\)"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_res,
      include_qualifier = FALSE,
      include_istd = FALSE,
      clear_existing = FALSE,
      max.dratio.sd.conc.bqc = 0.5),
    "Feature QC filters were updated\\: 4 \\(before 11\\)"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 4)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_res,
      clear_existing = FALSE,
      include_qualifier = FALSE,
      include_istd = FALSE,
      max.cv.conc.bqc = 25),
    "Replaced following previously defined QC filters\\: \\%CV"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)),4)

    mexp_res <- filter_features_qc(
      mexp_res,
      clear_existing = TRUE,
      include_qualifier = FALSE,
      include_istd = FALSE,
      min.signalblank.median.spl.pblk = 100,
      max.dratio.sd.conc.bqc = 0.8,
      max.cv.conc.bqc = 20)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

    mexp_res <- filter_features_qc(
      mexp_res,
      clear_existing = TRUE,
      include_qualifier = FALSE,
      include_istd = FALSE,
      min.signalblank.median.spl.pblk = 100,
      max.dratio.sd.conc.bqc = 0.8,
      max.cv.conc.bqc = 25)

    expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

      mexp_res <- filter_features_qc(
        mexp_res,
        clear_existing = FALSE,
        include_qualifier = FALSE,
        include_istd = FALSE,
        max.cv.conc.bqc = 15)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

      mexp_res <- filter_features_qc(
        mexp_proc,
        clear_existing = FALSE,
        include_qualifier = FALSE,
        include_istd = FALSE,
        min.signalblank.median.spl.pblk = 100,
        max.cv.conc.bqc = 23)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)


      mexp_res <- filter_features_qc(
        mexp_res,
        clear_existing = FALSE,
        include_qualifier = FALSE,
        include_istd = FALSE,
        response.curves.selection = c(1,2),
        response.curves.summary = "worst",
        max.slope.response = 1.005)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

      mexp_res <- filter_features_qc(
        mexp_res,
        clear_existing = FALSE,
        include_qualifier = FALSE,
        include_istd = FALSE,
        response.curves.selection = c(1,2),
        response.curves.summary = "best",
        max.slope.response = 1.005)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 9)

      mexp_res <- filter_features_qc(
        mexp_res,
        clear_existing = FALSE,
        include_qualifier = FALSE,
        include_istd = FALSE,
        response.curves.selection = c(1,2),
        response.curves.summary = "best",
        max.slope.response = 1.002)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 9)


})

# Confirm clearing QC filter works
test_that("Confirm clearing QC filter works", {

  mexp_res2 <- filter_features_qc(
    mexp,
    include_qualifier = FALSE,
    include_istd = FALSE,
    clear_existing = TRUE,
    min.signalblank.median.spl.pblk = 40,
    min.intensity.median.bqc = 1E2,
    max.cv.conc.bqc = 26,
    max.dratio.sd.conc.bqc = 0.8,
    max.prop.missing.conc.spl = 0.2,
    response.curves.selection = c(1),
    response.curves.summary = "mean",
    max.slope.response = 1.05)

  expect_message(
    mexp_res3 <- filter_features_qc(
      mexp_res2,
      clear_existing = TRUE,
      include_qualifier = FALSE,
      include_istd = FALSE),
    "all feature QC filters\\! All 19 quantifier")
  expect_equal(length(unique(mexp_res3@dataset_filtered$feature_id)), 19)

})

test_that("Confirm replace QC criteria category works", {

  expect_message(
    mexp_res2 <- filter_features_qc(
      mexp_proc,
      include_qualifier = FALSE,
      include_istd = FALSE,
      clear_existing = FALSE,
      min.signalblank.median.spl.pblk = 100,
      min.intensity.median.bqc = 1E3,
      max.cv.conc.bqc = 23,
      max.dratio.sd.conc.bqc = 0.7,
      max.prop.missing.conc.spl = 0.1,
      response.curves.selection = c(1,2),
      response.curves.summary = "best",
      max.slope.response = 1.05),
    "New feature QC filters were defined: 10 of 19 quantifier"
  )
  expect_equal(length(unique(mexp_res2@dataset_filtered$feature_id)), 10)

  expect_message(
    mexp_res3 <- filter_features_qc(
      mexp_res2,
      include_qualifier = FALSE,
      include_istd = FALSE,
      clear_existing = FALSE,
      min.signalblank.median.spl.pblk = 40,
      min.intensity.median.bqc = 1E2,
      max.cv.conc.bqc = 26,
      max.dratio.sd.conc.bqc = 0.8,
      max.prop.missing.conc.spl = 0.2,
      response.curves.selection = c(1),
      response.curves.summary = "mean",
      max.slope.response = 1.05),
    "Replaced following previously defined QC filters\\: Missing Values, Min-Intensity, Signal-to-Blank, \\%CV, D-ratio, and Linearity"
  )

  expect_equal(length(unique(mexp_res3@dataset_filtered$feature_id)), 10)

  expect_message(
    mexp_res4 <- filter_features_qc(
      mexp_res3,
      include_qualifier = TRUE,
      include_istd = TRUE,
      clear_existing = FALSE,
      min.signalblank.median.spl.pblk = 40,
      min.intensity.median.bqc = 1E2,
      max.cv.conc.bqc = 26,
      max.dratio.sd.conc.bqc = 0.8,
      max.prop.missing.conc.spl = 0.2,
      response.curves.selection = c(1),
      response.curves.summary = "mean",
      max.slope.response = 1.05),
    "Replaced following previously defined QC filters\\: Missing Values, Min-Intensity, Signal-to-Blank, \\%CV, D-ratio, Linearity, ISTD, and Qualifier"
  )

  expect_message(
    mexp_res4 <- filter_features_qc(
      mexp_res3,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      min.signalblank.median.spl.pblk = 40,
      min.intensity.median.bqc = 1E2,
      max.cv.conc.bqc = 26,
      max.dratio.sd.conc.bqc = 0.8,
      max.prop.missing.conc.spl = 0.2,
      response.curves.selection = c(1),
      response.curves.summary = "mean",
      max.slope.response = 1.05),
    "10 \\(before 10\\) of 19 quantifier and 1 of 1 qualifier features meet QC criteria \\(not including the 9 quantifier and 0 qualifier ISTD features"
  )

  expect_message(
    mexp_res2 <- filter_features_qc(
      mexp_proc,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = TRUE,
      min.signalblank.median.spl.pblk = 100,
      min.intensity.median.bqc = 1E3,
      max.cv.conc.bqc = 23,
      max.dratio.sd.conc.bqc = 0.7,
      max.prop.missing.conc.spl = 0.1,
      response.curves.selection = c(1,2),
      response.curves.summary = "best",
      max.slope.response = 1.05),
    "New feature QC filters were defined\\: 10 of 19 quantifier and 1 of 1 qualifier"
  )
  expect_equal(length(unique(mexp_res2@dataset_filtered$feature_id)), 11)

  expect_message(
    mexp_res3 <- filter_features_qc(
      mexp_res2,
      clear_existing = FALSE,
      include_qualifier = FALSE,
      include_istd = TRUE,
      min.signalblank.median.spl.pblk = 40,
      min.intensity.median.bqc = 1E2,
      max.cv.conc.bqc = 26,
      max.dratio.sd.conc.bqc = 0.8,
      max.prop.missing.conc.spl = 0.2,
      response.curves.selection = c(1),
      response.curves.summary = "mean",
      max.slope.response = 1.05),
    "Replaced following previously defined QC filters\\: Missing Values, Min-Intensity, Signal-to-Blank, \\%CV, D-ratio, Linearity, ISTD, and Qualifier"
  )

  expect_equal(length(unique(mexp_res3@dataset_filtered$feature_id)), 10)

  expect_message(
    mexp_res3 <- filter_features_qc(
      mexp_res2,
      include_istd = FALSE,
      include_qualifier = TRUE,
      clear_existing = FALSE,
      min.signalblank.median.spl.pblk = 40,
      min.intensity.median.bqc = 1E2,
      max.cv.conc.bqc = 14,
      max.dratio.sd.conc.bqc = 0.8,
      max.prop.missing.conc.spl = 0.2,
      response.curves.selection = c(1),
      response.curves.summary = "mean",
      max.slope.response = 1.05),
    "4 (before 10) of 19 quantifier and 0 of 1 qualifier", fixed = TRUE
    )


  expect_equal(length(unique(mexp_res3@dataset_filtered$feature_id)), 4)

  #clear all filters
  expect_message(
    mexp_res_cleared <- filter_features_qc(
      mexp_res2,
      clear_existing = TRUE,
      include_istd = FALSE,
      include_qualifier = TRUE),
    "all feature QC filters\\! All 19 quantifier and all 1 qualifier"
  )

  expect_equal(length(unique(mexp_res_cleared@dataset_filtered$feature_id)), 20)

  })


test_that("using filters without underlying data",{

  mexp_orig <- lipidomics_dataset
  mexpp <- calc_qc_metrics(mexp_orig,  use_batch_medians = FALSE)
  expect_error(
    mexp_res <- filter_features_qc(
      mexpp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      max.cv.conc.bqc = 20),
    "Cannot filter by `max.cv.conc.bqc` because concentration data")


})

test_that("filter_features_qc requires handles inconistent reponse filter setttings",{

  expect_error(
    mexp_res <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      min.rsquare.response = 0.8),
    "No response curves selected")

  expect_error(
    mexp_res <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      response.curves.selection = 1),
    "No response filters were defined")

  expect_error(
    mexp_res <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      response.curves.selection = 111,
    min.rsquare.response = 0.8),
    "The specified response curve index exceeds")

  expect_error(
    mexp_res <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      response.curves.selection = c(1,2),
      min.rsquare.response = 0.8),
    "Please set `response.curves.summary` to define curve")


  expect_error(
    mexp_res <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      response.curves.selection = "NotAValidSelection",
      min.rsquare.response = 0.8),
    "The following response curves are not defined in the metadata\\: NotAValidSelection")


  mexp_temp <- mexp
  mexp_temp@annot_responsecurves <- mexp_temp@annot_responsecurves[0,]
  expect_error(
    mexp_res <- filter_features_qc(
      mexp_temp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      response.curves.selection = 1,
      min.rsquare.response = 0.8
     ),
    "No response curve metadata found.")

  expect_no_error(
    mexp_temp2 <- filter_features_qc(
      mexp_temp,
      include_qualifier = TRUE,
      include_istd = TRUE,
      clear_existing = FALSE))

})

test_that("filter_features_qc requires specific args set",{

  expect_error(
    mexp_res <- filter_features_qc(
      mexp,
      clear_existing = FALSE,
      max.cv.conc.bqc = 20),
    "Argument `include_qualifier` is missing")

  expect_error(
    mexp_res <- filter_features_qc(
      mexp, include_qualifier = TRUE,
      clear_existing = FALSE,
      max.cv.conc.bqc = 20),
    "Argument `include_istd` is missing")

  expect_error(
    mexp_res <- filter_features_qc(
      mexp,
      include_istd = TRUE,
      clear_existing = FALSE,
      max.cv.conc.bqc = 20),
    "Argument `include_qualifier` is missing")
})


test_that("filter_features_qc handles user_defined_keepers",{
  expect_error(
    mexp_res <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = TRUE,
      max.cv.conc.bqc = 10,
      features.to.keep = c("CE 18:3")),
    "Following features defined via `features.to.keep` are not present in this dataset\\: CE 18\\:3")

  expect_error(
    mexp_res <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = TRUE,
      max.cv.conc.bqc = 10,
      features.to.keep = c("CE 18:1", "CE 18:3", "PC 40:6", "PC 40:9")),
    "Following features defined via `features\\.to\\.keep` are not present in this dataset\\: CE 18\\:3\\, and PC 40\\:9")

  expect_no_error(
    mexp_res <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = TRUE,
      max.cv.conc.bqc = 18,
      features.to.keep = NA))

  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 14)
  expect_false("CE 18:3" %in% mexp_res@dataset_filtered$feature_id)


    mexp_res <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = TRUE,
      max.cv.conc.bqc = 18,
      features.to.keep = c("CE 18:1"))

  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 14)
  expect_true("CE 18:1" %in% mexp_res@dataset_filtered$feature_id)

  expect_message(
    mexp_res2 <- filter_features_qc(
      mexp_res,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      max.cv.conc.bqc = 18,
      features.to.keep = c("PC 40:6")),
    "Replaced following previously defined QC filters\\: \\%CV\\, and Keepers")

  expect_message(
    mexp_res2 <- filter_features_qc(
      mexp_res2,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      max.cv.conc.bqc = 19,
      features.to.keep = c("PC 40:6")),
    "Replaced following previously defined QC filters\\: \\%CV")


  expect_message(
    mexp_res2 <- filter_features_qc(
      mexp,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = TRUE,
      max.cv.conc.bqc = 12,
      features.to.keep = c("CE 18:1", "PC 40:6")),
  "The following features were forced to be retained despite not meeting filtering criteria\\: CE 18\\:1\\, and PC 40\\:6")

  expect_equal(length(unique(mexp_res2@dataset_filtered$feature_id)), 5)
  expect_true(all(c("CE 18:1", "PC 40:6") %in% mexp_res2@dataset_filtered$feature_id))


    messages <- capture.output(
      mexp_res2 <- filter_features_qc(
        mexp,
        include_qualifier = TRUE,
        include_istd = FALSE,
        clear_existing = TRUE,
        max.cv.conc.bqc = 19999,
        features.to.keep = c("CE 18:1", "PC 40:6")),
      type = "message")
    expect_false(any(grepl("The following features were forced to be retained despite n", messages)))



  expect_equal(length(unique(mexp_res2@dataset_filtered$feature_id)), 20)
  expect_true(all(c("CE 18:1", "PC 40:6") %in% mexp_res2@dataset_filtered$feature_id))

  expect_message(
    mexp_res3 <- filter_features_qc(
      mexp_res2,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      max.cv.conc.bqc = 18),
    "Replaced following previously defined QC filters\\: \\%CV\\, and Keepers")

  expect_message(
    mexp_res3 <- filter_features_qc(
      mexp_res2,
      include_qualifier = TRUE,
      include_istd = FALSE,
      clear_existing = FALSE,
      max.cv.conc.bqc = 11),
    "1 \\(before 19\\) of 19 quantifier and 0 of 1 qualifier features meet QC criteria \\(not including the 9 quantifier and 0 qualifier ISTD features\\)")

  expect_equal(length(unique(mexp_res3@dataset_filtered$feature_id)), 1)
  expect_false(any(c("CE 18:1", "PC 40:6") %in% mexp_res3@dataset_filtered$feature_id))

})


