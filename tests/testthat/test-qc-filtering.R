library(testthat)
library(dplyr)

mexp_orig <- lipidomics_dataset
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
    expect_equal(max(mexp_res@metrics_qc$intensity_cv_SPL),76.0692868)
    expect_equal(max(mexp_res@metrics_qc$intensity_cv_TQC),30.1009462)
    expect_equal(max(mexp_res@metrics_qc$intensity_cv_BQC),33.296624)
    expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_SPL),103.7784475)
    expect_equal(max(mexp_res@metrics_qc$conc_cv_SPL),103.7784475)
    expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_TQC, na.rm = T),32.9608359974)
    expect_equal(max(mexp_res@metrics_qc$conc_cv_TQC, na.rm = T),32.9608359974)
    expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_BQC, na.rm = T), 34.8191142734)
    expect_equal(max(mexp_res@metrics_qc$conc_cv_BQC, na.rm = T),34.8191142734)
    expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_bqc_conc, na.rm = T), 0.69213784)
    expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_bqc_normint, na.rm = T),0.69213784)
    expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_tqc_conc, na.rm = T),0.507702049)
    expect_equal(median(mexp_res@metrics_qc$conc_dratio_mad_bqc_conc, na.rm = T),0.619781631)
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
  expect_equal(max(mexp_res@metrics_qc$norm_intensity_cv_BQC), 18.1254294266)
  expect_equal(max(mexp_res@metrics_qc$conc_cv_BQC),18.1254294266)
  expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_bqc_conc, na.rm = T), 0.438675798469)
  expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_bqc_normint, na.rm = T),0.438675798469)
  expect_equal(median(mexp_res@metrics_qc$conc_dratio_sd_tqc_conc, na.rm = T),0.42459409722)
  expect_equal(median(mexp_res@metrics_qc$conc_dratio_mad_bqc_conc, na.rm = T),0.409148315867)
  expect_equal(min(mexp_res@metrics_qc$r2_rqc_A),0.91931047)
  expect_equal(min(mexp_res@metrics_qc$r2_rqc_B),0.85693787)
  expect_equal(min(mexp_res@metrics_qc$slopenorm_rqc_A),0.69281368)
  expect_equal(min(mexp_res@metrics_qc$slopenorm_rqc_B),0.65001359)
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
                               include_calibration_results = TRUE),
               "MidarExperiment")

  expect_error(calc_qc_metrics(mexp_temp,
                              use_batch_medians = TRUE,
                              include_norm_intensity_stats = TRUE,
                              include_conc_stats = TRUE,
                              include_response_stats = TRUE,
                              include_calibration_results = TRUE),
               "No response curve metadata found")

  mexp_temp <- mexp
  mexp_temp@annot_responsecurves$analysis_id[1] <- "unknown1"
  mexp_temp@annot_responsecurves$analysis_id[3] <- "unknown2"

  expect_s4_class(calc_qc_metrics(mexp_temp,
                                  use_batch_medians = TRUE,
                                  include_norm_intensity_stats = TRUE,
                                  include_conc_stats = TRUE,
                                  include_response_stats = NA,
                                  include_calibration_results = TRUE),
                  "MidarExperiment")

  expect_error(calc_qc_metrics(mexp_temp,
                               use_batch_medians = TRUE,
                               include_norm_intensity_stats = TRUE,
                               include_conc_stats = TRUE,
                               include_response_stats = TRUE,
                               include_calibration_results = TRUE),
               "One or more analysis IDs")

})

test_that("filter_features_qc works on selected criteria", {
  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    min.intensity.median.bqc = 1E5)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 15)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    min.intensity.median.tqc = 1E5)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 16)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    min.intensity.median.bqc = 1E5,
    min.intensity.median.tqc = 1E5)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 15)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    min.signalblank.median.spl.pblk = 100)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    min.intensity.median.bqc = 1E5,
    min.signalblank.median.spl.pblk = 100)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 9)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    max.cv.conc.bqc = 20)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    min.signalblank.median.spl.pblk = 100,
    max.cv.conc.bqc = 20)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    min.signalblank.median.spl.pblk = 100,
    max.cv.normintensity.bqc = 20)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    max.dratio.sd.conc.bqc = 0.5)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 8)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    max.prop.missing.conc.spl = 0)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 19)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = 1,
    response.curves.summary = "mean",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 8)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = 2,
    response.curves.summary = "mean",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 9)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = c(1,2),
    response.curves.summary = "mean",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 6)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = c(1,2),
    response.curves.summary = "median",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 6)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = c(1,2),
    response.curves.summary = "best",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 12)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = c(1,2),
    response.curves.summary = "worst",
    min.rsquare.response = 0.98)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = c(1,2),
    response.curves.summary = "mean",
    min.slope.response = 0.9)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = c(1,2),
    response.curves.summary = "worst",
    min.slope.response = 0.9)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 9)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = c(1,2),
    response.curves.summary = "best",
    max.slope.response = 1.005)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 17)

  mexp_res <- filter_features_qc(
    mexp_proc, replace_existing = TRUE,
    response.curves.selection = c(1,2),
    response.curves.summary = "worst",
    max.slope.response = 1.005)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 12)
})


# Confirm overwriting of QC criteria works


test_that("Confirm overwriting of QC criteria works", {
  mexp_res <- filter_features_qc(
    mexp_proc,
    replace_existing = FALSE,
    min.signalblank.median.spl.pblk = 300)

  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 8)


  expect_message(
   mexp_res <- filter_features_qc(
    mexp_res,
    replace_existing = FALSE,
    min.signalblank.median.spl.pblk = 100),
   "Following previously defined QC filters were replaced: Signal-to-Blank"
   )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 11)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_res,
      replace_existing = FALSE,
      max.cv.conc.bqc = 20),
    "QC filter criteria were updated\\: 5 \\(before 11\\)"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_res,
      replace_existing = FALSE,
      max.dratio.sd.conc.bqc = 0.8),
    "QC filter criteria were updated\\: 3 \\(before 5\\)"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 3)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_res,
      replace_existing = FALSE,
      max.cv.conc.bqc = 25),
    "Following previously defined QC filters were replaced"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

    mexp_res <- filter_features_qc(
      mexp_res,
      replace_existing = TRUE,
      min.signalblank.median.spl.pblk = 100,
      max.dratio.sd.conc.bqc = 0.8,
      max.cv.conc.bqc = 20)
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 3)

    mexp_res <- filter_features_qc(
      mexp_res,
      replace_existing = TRUE,
      min.signalblank.median.spl.pblk = 100,
      max.dratio.sd.conc.bqc = 0.8,
      max.cv.conc.bqc = 25)

    expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

      mexp_res <- filter_features_qc(
        mexp_res,
        replace_existing = FALSE,
        max.cv.conc.bqc = 20)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 3)

      mexp_res <- filter_features_qc(
        mexp_proc,
        replace_existing = FALSE,
        min.signalblank.median.spl.pblk = 100,
        max.cv.conc.bqc = 23)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 9)


      mexp_res <- filter_features_qc(
        mexp_res,
        replace_existing = FALSE,
        response.curves.selection = c(1,2),
        response.curves.summary = "worst",
        max.slope.response = 1.005)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 5)

      mexp_res <- filter_features_qc(
        mexp_res,
        replace_existing = FALSE,
        response.curves.selection = c(1,2),
        response.curves.summary = "best",
        max.slope.response = 1.005)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 8)

      mexp_res <- filter_features_qc(
        mexp_res,
        replace_existing = FALSE,
        response.curves.selection = c(1,2),
        response.curves.summary = "best",
        max.slope.response = 1.002)
      expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 8)

})

test_that("Confirm replace QC criteria category works", {

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_proc,
      replace_existing = FALSE,
      min.signalblank.median.spl.pblk = 100,
      min.intensity.median.bqc = 1E3,
      max.cv.conc.bqc = 23,
      max.dratio.sd.conc.bqc = 0.7,
      max.prop.missing.conc.spl = 0.1,
      response.curves.selection = c(1,2),
      response.curves.summary = "best",
      max.slope.response = 1.05),
    "New QC filter criteria were defined: 3 of 19 quantifier"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 3)

  expect_message(
    mexp_res2 <- filter_features_qc(
      mexp_res,
      replace_existing = FALSE,
      min.signalblank.median.spl.pblk = 40,
      min.intensity.median.bqc = 1E2,
      max.cv.conc.bqc = 26,
      max.dratio.sd.conc.bqc = 0.8,
      max.prop.missing.conc.spl = 0.2,
      response.curves.selection = c(1),
      response.curves.summary = "mean",
      max.slope.response = 1.05),
    "Following previously defined QC filters were replaced\\: Missing Values, Min-Intensity, Signal-to-Blank, \\%CV, D-ratio, and Linearity"
  )

  expect_equal(length(unique(mexp_res2@dataset_filtered$feature_id)), 6)

  expect_message(
    mexp_res <- filter_features_qc(
      mexp_proc,
      qualifier.include = TRUE,
      replace_existing = TRUE,
      min.signalblank.median.spl.pblk = 100,
      min.intensity.median.bqc = 1E3,
      max.cv.conc.bqc = 23,
      max.dratio.sd.conc.bqc = 0.7,
      max.prop.missing.conc.spl = 0.1,
      response.curves.selection = c(1,2),
      response.curves.summary = "best",
      max.slope.response = 1.05),
    "New QC filter criteria were defined: 3 of 19 quantifier and 0 of 1 qualifier"
  )
  expect_equal(length(unique(mexp_res@dataset_filtered$feature_id)), 3)

  expect_message(
    mexp_res2 <- filter_features_qc(
      mexp_res,
      qualifier.include = TRUE,
      replace_existing = FALSE,
      min.signalblank.median.spl.pblk = 40,
      min.intensity.median.bqc = 1E2,
      max.cv.conc.bqc = 26,
      max.dratio.sd.conc.bqc = 0.8,
      max.prop.missing.conc.spl = 0.2,
      response.curves.selection = c(1),
      response.curves.summary = "mean",
      max.slope.response = 1.05),
    "Following previously defined QC filters were replaced\\: Missing Values, Min-Intensity, Signal-to-Blank, \\%CV, D-ratio, and Linearity"
  )

  expect_equal(length(unique(mexp_res2@dataset_filtered$feature_id)), 7)

  expect_message(
    mexp_res2 <- filter_features_qc(
      mexp_res,
      qualifier.include = TRUE,
      replace_existing = FALSE,
      min.signalblank.median.spl.pblk = 40,
      min.intensity.median.bqc = 1E2,
      max.cv.conc.bqc = 26,
      max.dratio.sd.conc.bqc = 0.8,
      max.prop.missing.conc.spl = 0.2,
      response.curves.selection = c(1),
      response.curves.summary = "mean",
      max.slope.response = 1.05),
    "updated: 6 \\(before 3\\) of 19 quantifier and 1 of 1 qualifier"
    )


  expect_equal(length(unique(mexp_res2@dataset_filtered$feature_id)), 7)


  })



