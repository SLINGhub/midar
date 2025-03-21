# library(fs)
# library(vdiffr)
# library(ggplot2)
# library(testthat)
set.seed(123)


# TODO: add tests with a (published) reference dataset and results to compare specific functions (such as LOESS)

mexp_orig <- lipidomics_dataset

mexp_raw <- exclude_analyses(mexp_orig, analyses = c("Longit_batch6_51"), clear_existing  = TRUE)
mexp_norm <- normalize_by_istd(mexp_raw)
mexp <- quantify_by_istd(mexp_norm)

mexp_err <- mexp_raw
mexp_err@dataset$qc_type <- if_else(str_detect(mexp_err@dataset$analysis_id, "ISTD"), "BQC", mexp_err@dataset$qc_type)
mexp_err <- normalize_by_istd(mexp_err)
mexp_err <- quantify_by_istd(mexp_err)

mexp_err2 <- mexp_raw
mexp_err2 <- normalize_by_istd(mexp_err2)
mexp_err2 <- quantify_by_istd(mexp_err2)

mexp_err2@dataset$feature_conc[581] <- 0
mexp_err2@dataset$feature_conc[611] <- 0
mexp_err2@dataset$feature_conc[1123] <- 0




test_that("correct_drift_gaussiankernel works", {

  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "conc",
      conditional_correction = FALSE,
      kernel_size = 10,
      batch_wise = TRUE,
      ref_qc_types = "SPL",
      ignore_istd = FALSE),
    "Drift correction was applied to 20 of 29 features (batch-wise)", fixed = TRUE
  )

    # drift is not calculated after correction (recalc_trend_after = FALSE)
  expect_true(
    all(is.na(mexp_drift@dataset$feature_conc_fit_after)))


  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "conc",
      conditional_correction = FALSE,
      kernel_size = 10,
      batch_wise = TRUE,
      ref_qc_types = "SPL",
      ignore_istd = FALSE),
    "Drift correction was applied to 20 of 29 features (batch-wise)", fixed = TRUE
  )

    expect_message(
      mexp_drift <- correct_drift_gaussiankernel(
        mexp,
        variable = "conc",
        conditional_correction = FALSE,
        outlier_ksd = 1,
        outlier_filter = TRUE,
        kernel_size = 10,
        batch_wise = TRUE,
        ref_qc_types = "SPL",
        ignore_istd = TRUE),
    "-1.70% to 2.39%",
    fixed = TRUE
  )

  # Check conc and intensity
  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "conc",
      conditional_correction = FALSE,
      kernel_size = 10,
      batch_wise = TRUE,
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "(range: -2.47% to 0.04%)",
    fixed = TRUE
  )

  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "intensity",
      conditional_correction = F,
      kernel_size = 10,
      batch_wise = TRUE,
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "(range: -2.10% to -0.27%)",
    fixed = TRUE
  )

  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "intensity",
      conditional_correction = F,
      kernel_size = 10,
      batch_wise = TRUE,
      ref_qc_types = "SPL",
      ignore_istd = FALSE),
    "Drift correction was applied to 29 of 29 features",
    fixed = TRUE
  )

  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "intensity",
      conditional_correction = TRUE,
      cv_diff_threshold = 0,
      kernel_size = 10,
      batch_wise = TRUE,
      ref_qc_types = "SPL",
      ignore_istd = FALSE),
    "was applied to 29 of 29 features",
    fixed = TRUE
  )


  # check within batch FALSE
  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "conc",
      conditional_correction = FALSE,
      kernel_size = 10,
      batch_wise = FALSE,
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "(range: -11.46% to -0.74%)",
    fixed = TRUE
  )

  # including istd: no difference, as ISTD conc is constant, and this exclude from calc
  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "conc",
      conditional_correction = TRUE,
      kernel_size = 10,
      batch_wise = FALSE,
      ref_qc_types = "SPL",
      ignore_istd = FALSE),
    "(range: -11.46% to -0.74%)",
    fixed = TRUE
  )

  # however with intensity, there must be a difference when including/excluding ISTD
  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "intensity",
      conditional_correction = T,
      kernel_size = 10,
      batch_wise = FALSE,
      ref_qc_types = "SPL",
      ignore_istd = FALSE),
    "(range: -19.06% to -0.51%)",
    fixed = TRUE
  )

  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      variable = "intensity",
      conditional_correction = T,
      kernel_size = 10,
      batch_wise = FALSE,
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "(range: -12.78% to -0.99%)",
    fixed = TRUE
  )
})

test_that("using sample types other than SPL", {
  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = TRUE,
      variable = "conc",
      ref_qc_types = c("BQC", "TQC"),
      ignore_istd = TRUE),
    "(range: -1.39% to 0.92%)",
    fixed = TRUE
  )

  # using sample types other than SPL, CV of SPL increases
  expect_message(
    mexp_drift <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = TRUE,
      variable = "conc",
      ref_qc_types = c("BQC", "TQC"),
      ignore_istd = TRUE),
    "decreased from",
    fixed = TRUE
  )
})

test_that("replace_previous FALSE works", {
  mexp_drift2 <- correct_drift_gaussiankernel(
    mexp,
    batch_wise = FALSE,
    replace_previous = TRUE,
    variable = "conc",
    conditional_correction = F,
    kernel_size = 10,
    ref_qc_types = "BQC",
    ignore_istd = TRUE)

  expect_message(
    mexp_drift3 <- correct_drift_gaussiankernel(
      mexp_drift2,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      conditional_correction = F,
      kernel_size = 10,
      ref_qc_types = "SPL"),
    "Replacing previous `conc` drift corrections",
    fixed = TRUE
  )

  expect_message(
    mexp_drift3 <- correct_drift_gaussiankernel(
      mexp_drift2,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      conditional_correction = F,
      kernel_size = 10,
      ignore_istd = TRUE,
      ref_qc_types = "SPL"),
    "-11.46% to -0.74%)",
    fixed = TRUE
  )

  # Check replace_previous TRUE again
  expect_message(
    mexp_drift4 <- correct_drift_gaussiankernel(
      mexp_drift3,
      batch_wise = FALSE,
      replace_previous = FALSE,
      variable = "conc",
      conditional_correction = F,
      kernel_size = 5,
      ignore_istd = TRUE,
      ref_qc_types = "SPL"),
    "Adding correction on top of previous",
    fixed = TRUE
  )

  expect_message(
    mexp_drift4 <- correct_drift_gaussiankernel(
      mexp_drift3,
      batch_wise = FALSE,
      replace_previous = FALSE,
      variable = "conc",
      conditional_correction = F,
      kernel_size = 5,
      ignore_istd = TRUE,
      ref_qc_types = "SPL"),
    "-3.66% to -0.57%)",
    fixed = TRUE
  )
})

  # applying corrections to a variable of 'lower' processing order, will invalidate all processing steps that are based on this variable
test_that("applying corrections to a variable of 'lower' processing order is working as it should", {
  expect_message(
    mexp_drift4 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = FALSE,
      replace_previous = FALSE,
      variable = "intensity",
      conditional_correction = F,
      kernel_size = 10,
      ignore_istd = TRUE,
      ref_qc_types = "SPL"),
    "normalized intensities and concentrations are no longer valid. Please reprocess",
    fixed = TRUE
  )

  expect_message(
    mexp_drift4 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = FALSE,
      replace_previous = FALSE,
      variable = "norm_intensity",
      conditional_correction = F,
      kernel_size = 10,
      ignore_istd = TRUE,
      ref_qc_types = "SPL"),
    "Concentrations are no longer valid. Please reprocess",
    fixed = TRUE
  )

  expect_message(
    mexp_drift4 <- correct_drift_gaussiankernel(
      mexp_norm,
      batch_wise = FALSE,
      replace_previous = FALSE,
      variable = "intensity",
      conditional_correction = F,
      kernel_size = 10,
      ignore_istd = TRUE,
      ref_qc_types = "SPL"),
    "Normalized intensities are no longer valid. Please reprocess",
    fixed = TRUE
  )
  })

test_that("recalc_trend_after works", {
  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = TRUE,
      replace_previous = TRUE,
      recalc_trend_after = TRUE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = FALSE),
    "-2.47% to 0.04%)",
    fixed = TRUE
  )

  expect_equal(
    max(mexp_drift2@dataset$feature_conc_fit_after), 1.359779940)


  p <- plot_runscatter(mexp_drift2, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("gaussiankernel runscatter plot before 1 ", p[[1]])

  p <- plot_runscatter(mexp_drift2, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("gaussiankernel runscatter plot after 1 ", p[[1]])

})

test_that("Scale smooth works", {

  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      scale_smooth = TRUE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "decreased",
    fixed = TRUE
  )

  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      recalc_trend_after = TRUE,
      scale_smooth = TRUE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = FALSE),
    "range: -14.86% to -1.31%",
    fixed = TRUE
  )
  p <- plot_runscatter(mexp_drift2, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("gaussiankernel_runscatter_scalesmooth_after_1 ", p[[1]])

})

# conditional correction
# result when correcting all
#The median CV change of all features in study samples was -1.91% (range: -11.46% to -0.74%). The median absolute CV of all features decreased from 33.81% to 32.22%.

test_that("conditional correction works", {
  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      conditional_correction = TRUE,
      cv_diff_threshold = 0,
      scale_smooth = FALSE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "-11.46% to -0.74%",
    fixed = TRUE
  )

  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      conditional_correction = TRUE,
      cv_diff_threshold = 2,
      scale_smooth = FALSE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "-11.46% to -0.74%",
    fixed = TRUE
  )

  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      conditional_correction = TRUE,
      cv_diff_threshold = -1,
      scale_smooth = FALSE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "-11.46% to -1.10",
    fixed = TRUE
  )

  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      conditional_correction = TRUE,
      cv_diff_threshold = -6,
      scale_smooth = FALSE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "-11.46% to -6.18%)",
    fixed = TRUE
  )

  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      conditional_correction = TRUE,
      cv_diff_threshold = -6,
      scale_smooth = FALSE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "applied to 3 of 20 features (across all batches)",
    fixed = TRUE
  )

  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = TRUE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      conditional_correction = TRUE,
      cv_diff_threshold = 0,
      scale_smooth = FALSE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "-4.48% to -0.43%)",
    fixed = TRUE
  )

  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = TRUE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      conditional_correction = TRUE,
      cv_diff_threshold = 11,
      scale_smooth = FALSE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "-2.47% to 0.04%)",
    fixed = TRUE
  )
  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = TRUE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      conditional_correction = TRUE,
      cv_diff_threshold = -3,
      scale_smooth = FALSE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = TRUE),
    "-7.13% to -3.28%",
    fixed = TRUE
  )


  expect_message(
    mexp_drift2 <- correct_drift_gaussiankernel(
      mexp,
      batch_wise = TRUE,
      replace_previous = TRUE,
      recalc_trend_after = FALSE,
      conditional_correction = TRUE,
      cv_diff_threshold = -3,
      scale_smooth = FALSE,
      variable = "conc",
      ref_qc_types = "SPL",
      ignore_istd = FALSE),
    "-7.13% to -3.28%",
    fixed = TRUE
  )
})

# test_that("correct_drift_loess works", {
#   expect_message(
#     mexp_drift <- correct_drift_loess(
#       mexp,
#       batch_wise = FALSE,
#       variable = "conc",
#       ref_qc_types = "SPL",
#       ignore_istd = FALSE),
#     "(-60.7 to 15.9%)", fixed = T)
#
#
#   expect_message(
#     mexp_drift_batch <- correct_batch_centering(
#       mexp,
#       correct_location = FALSE,
#       correct_scale = FALSE,
#       ref_qc_types = "BQC",
#       variable = "conc"),
#   "(-14.1 to 59.1%)")
#
#   mexp_drift_batch <- correct_batch_centering(
#     mexp,
#     correct_location = TRUE,
#     correct_scale = TRUE,
#     ref_qc_types = "BQC",
#     variable = "conc")
#
# })


test_that("correct_drift_loess works", {
  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp,
      span = 0.75,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "-0.29% to 1.40%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_drift_loess_before ", p[[2]])

  p <- plot_runscatter(mexp_drift1, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_drift_loess_after", p[[2]])

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp,
      span = 0.75,
      batch_wise = TRUE,
      replace_previous = TRUE,
      log_transform_internal = FALSE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "-0.56% to 1.42%",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp,
      span = 0.75,
      batch_wise = TRUE,
      replace_previous = TRUE,
      log_transform_internal = FALSE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      feature_list = c("CE 18:1", "PC 40:8"),
      ignore_istd = TRUE),
    "Drift correction was applied to 2 of 20 features (batch-wise)",
    fixed = TRUE)

  expect_error(
    mexp_drift1 <- correct_drift_loess(
      mexp,
      span = 0.75,
      batch_wise = TRUE,
      replace_previous = TRUE,
      log_transform_internal = FALSE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      feature_list = c("CE 18:1", "NO PC 40:9","NOPE"),
      ignore_istd = TRUE),
    "One or more feature(s) specified with `feature_list` are not present in the dataset",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp,
      span = 0.75,
      batch_wise = TRUE,
      replace_previous = TRUE,
      log_transform_internal = FALSE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      feature_list = "^PC",
      ignore_istd = TRUE),
    "Drift correction was applied to 4 of 20 features (batch-wise)",
    fixed = TRUE)

  expect_error(
    mexp_drift1 <- correct_drift_loess(
      mexp,
      span = 0.75,
      batch_wise = TRUE,
      replace_previous = TRUE,
      log_transform_internal = FALSE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      feature_list = "^NOPE",
      ignore_istd = TRUE),
    "The feature filter set via `feature_list` does not match any feature in the dataset",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp_err2,
      span = 0.75,
      batch_wise = TRUE,
      replace_previous = TRUE,
      log_transform_internal = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "2 feature(s) contain one or more zero or negative `conc` values",
    fixed = TRUE)
})

test_that("correct_drift_loess warnigs report work", {
  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0.5,
      degree = 2,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Issues (warnings) occured during smoothing of all features",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0.75,
      degree = 2,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
  "Issues (warnings) occured during the smoothing of 11 feature(s)",
  fixed = TRUE)


  expect_error(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0.75,
      degree = 2,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = c("BQC","NOQC"),
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "One or more specified `qc_types` are not present",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0.75,
      degree = 2,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = c("SPL"),
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "27 of 62 BQCs, 6 of 41 TQCs, 3 of 3 LTRs",
    fixed = TRUE)
})


test_that("correct_drift_gaussiankernel fit error are handeled", {

  expect_error(
    mexp_drift1 <- correct_drift_gaussiankernel(
      mexp_err,
      log_transform_internal = FALSE, # not supported at moment
      kernel_size = 10,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Currently only `log_transform_internal = TRUE` is supported.",
    fixed = TRUE)


  expect_error(
    mexp_drift1 <- correct_drift_gaussiankernel(
      mexp_err,
      kernel_size = 0,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Argument `kernel_size` must larger than 0",
    fixed = TRUE)

  expect_error(
    mexp_drift1 <- correct_drift_gaussiankernel(
      mexp_err,
      kernel_size = 10,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      outlier_ksd = 0,
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Argument `outlier_ksd` must larger than 0",
    fixed = TRUE)

  # Could not find data or arguments to make this function fail or throw a warning
})

test_that("correct_drift_loess fit error are handeled", {

  expect_error(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0,
      degree = 1,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Argument `span` must larger than 0",
    fixed = TRUE)

  expect_error(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 1,
      degree = 3,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Argument `degree` must be 1 or 2",
    fixed = TRUE)

  expect_error(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0.3,
      degree = 2,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Smoothing failed for 4 feature(s) in all batches",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0.5,
      degree = 1,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Issues (warnings) occured during the smoothing of 17 feature(s)",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0.5,
      degree = 1,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "2 of 41 TQCs, 2 of 3 LTRs were excluded from correction as",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0.5,
      degree = 2,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Issues (warnings) occured during smoothing of all features ",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_loess(
      mexp_err,
      span = 0.75,
      degree = 2,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      ignore_istd = TRUE),
    "Issues (warnings) occured during the smoothing of 11 feature(s)",
    fixed = TRUE)

})

test_that("fits resulting in invalid values are handeled", {
  expect_message(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      cv = FALSE,
      penalty = 1.23,   # penalty too high for some features resulting in extreme values for these
      spar = NULL,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "4 features have invalid values after smoothing. NA will be be returned ",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_drift_cubicspline_withinvalid_smooths_1 ", p[[2]])


  expect_message(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      cv = FALSE,
      penalty = 1.23,   # penalty too high for some features resulting in extreme values for these
      spar = NULL,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = TRUE,
      ignore_istd = TRUE),
    "4 features have invalid values after smoothing. The original values were kept for these features",
    fixed = TRUE)

  expect_message(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      cv = FALSE,
      penalty = 1.23,   # penalty too high for some features resulting in extreme values for these
      spar = NULL,
      batch_wise = TRUE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = TRUE,
      ignore_istd = TRUE),
  "-0.31% to 6.99%",
  fixed = TRUE)

})

test_that("correct_drift_cubicspline works", {
  expect_message(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.21% to 3.56%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_drift_cubicspline_basic_after", p[[2]])

  expect_message(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      cv = FALSE,  # use ‘generalized’ cross-validation (GCV)
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.31% to 2.90%)",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_drift_cubicspline_basic_cvfalse_after", p[[2]])


  expect_message(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      log_transform_internal = FALSE,
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.42% to 3.63%",
    fixed = TRUE)


  expect_message(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      cv = TRUE,  # use ‘generalized’ cross-validation (GCV)
      lambda = 0.01,   # define a fixed lambda
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-7.87% to 1.93%)",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("drift_cubicspline_withlambda_bef", p[[2]])

  expect_message(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      cv = TRUE,  # use ‘generalized’ cross-validation (GCV)
      penalty = 1.1,   # define a fixed lambda
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.21% to 3.56%)",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("cubicspline_withpanalty_bef", p[[2]])


  expect_message(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      cv = TRUE,
      spar = 0.8,
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.46% to 1.79%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("cubicspline_withspar_bef", p[[2]])

  expect_error(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      cv = TRUE,
      spar = NA,# use ‘generalized’ cross-validation (GCV)
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "must be NULL or numeric",
    fixed = TRUE)

  expect_error(
    mexp_drift1 <- correct_drift_cubicspline(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      cv = TRUE,
      spar = 0.7,# use ‘generalized’ cross-validation (GCV)
      lambda = 0.1,   # define a fixed lambda
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "Either `spar` or `lambda` can be specified, not both",
    fixed = TRUE)

})


test_that("correct_drift_gam works", {
  expect_message(
    mexp_drift1 <- correct_drift_gam(
      mexp,bs = "ps",
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.53% to 0.83%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_drift_gam_ps_before1", p[[2]])

  expect_message(
    mexp_drift1 <- correct_drift_gam(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      bs = "tp",  # thin plate
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.52% to 0.73%)",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_drift_gam_tp_before1", p[[2]])

  expect_message(
    mexp_drift1 <- correct_drift_gam(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      sp = 0.01,
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.19% to 2.41%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_drift_gam_sp001_before1", p[[2]])

  expect_message(
    mexp_drift1 <- correct_drift_gam(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      k = 10,
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.53% to 0.83%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_drift1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_drift_gam_k10_before1", p[[2]])


  expect_message(
    mexp_drift1 <- correct_drift_gam(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      log_transform_internal = FALSE,
      k = 10,
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "-8.84% to 0.79%",
    fixed = TRUE)

  expect_error(
    mexp_drift1 <- correct_drift_gam(
      mexp,
      batch_wise = FALSE,
      replace_previous = TRUE,
      variable = "conc",
      log_transform_internal = FALSE,
      sp = FALSE,
      ref_qc_types = "BQC",
      recalc_trend_after = TRUE,
      use_original_if_fail = FALSE,
      ignore_istd = TRUE),
    "Argument `sp` must be NULL or numeric",
    fixed = TRUE)
})


mexp_dcorr <- correct_drift_gaussiankernel(
  mexp,
  variable = "conc",
  conditional_correction = FALSE,
  kernel_size = 10,
  batch_wise = TRUE,recalc_trend_after = TRUE,
  ref_qc_types = "SPL",
  ignore_istd = FALSE)


test_that("correct_batch_centering works", {

  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      ref_qc_types = "SPL",
      variable = "conc"),
    "Adding batch correction on top of `conc` drift-correction",
    fixed = TRUE)

  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      ref_qc_types = "SPL",
      variable = "conc"),
    "Batch median-centering of 6 batches was applied to drift-corrected concentrations of all 20 features",
    fixed = TRUE)

  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      ref_qc_types = "SPL",
      variable = "conc"),
    "range: -8.40% to 2.20%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_batch1, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("batch_centering_batchcenter1", p[[2]])

 #replacing the trend curves from gaussiankernel with the new batch centering trend lines
  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      ref_qc_types = "SPL",
      variable = "conc",
      replace_exisiting_trendcurves = TRUE),
    "-8.40% to 2.20%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_batch1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("correct_batch_centering_replacetrends ", p[[2]])

  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp,
      correct_scale = FALSE,
      ref_qc_types = "SPL",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "-9.40% to 2.60%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_batch1, variable = "conc_before", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("batch_centering_nodriftbefore ", p[[2]])

  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp,
      correct_scale = TRUE,
      ref_qc_types = "SPL",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "-11.00% to 2.80%",
    fixed = TRUE)

  p <- plot_runscatter(mexp_batch1, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("batch_centering_correctscalelocation ", p[[2]])


  #TODO: log transform has no impcat on the scale correction?
  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      replace_previous = TRUE,
      ref_qc_types = c("BQC", "TQC"),
      variable = "conc",
      log_transform_internal = FALSE,
      replace_exisiting_trendcurves = FALSE),
    "-5.90% to 4.90%)",
    fixed = TRUE)

  expect_error(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = TRUE,
      replace_previous = TRUE,
      ref_qc_types = c("BQC", "TQC"),
      variable = "conc",
      log_transform_internal = FALSE,
      replace_exisiting_trendcurves = FALSE),
    "Currently data must be log-transformed for batch scaling",
    fixed = TRUE)

  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      replace_previous = TRUE,
      ref_qc_types = c("BQC", "TQC"),
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "-5.90% to 4.90%)",
    fixed = TRUE)

  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      replace_previous = TRUE,
      ref_qc_types = c("BQC", "TQC"),
      variable = "intensity",
      replace_exisiting_trendcurves = FALSE),
    "-16.40% to 1.40%)",
    fixed = TRUE)

  expect_error(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      replace_previous = TRUE,
      ref_qc_types = c("BQC", "EQC"),
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "One or more specified `qc_types` are not present ",
    fixed = TRUE)

})

test_that("correct_batch_centering works with replace_previous", {

  # add on top of previous with only drift correction before is same result as replace previous
  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      replace_previous = FALSE,
      ref_qc_types = "BQC",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "-4.30% to 3.10%)",
    fixed = TRUE)

  p <- plot_runscatter(mexp_batch1, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  expect_message(
    mexp_batch2 <- correct_batch_centering(
      mexp_batch1,
      correct_scale = FALSE,
      replace_previous = FALSE,
      ref_qc_types = "SPL",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "-4.10% to -0.10%)",
    fixed = TRUE)

  p <- plot_runscatter(mexp_batch2, variable = "conc", qc_types = c("SPL", "BQC"),
                       show_trend = T, include_istd = FALSE, return_plots = TRUE)

  vdiffr::expect_doppelganger("batch_centering_replaceprevious", p[[2]])

  expect_message(
    mexp_batch2 <- correct_batch_centering(
      mexp_batch1,
      correct_scale = FALSE,
      replace_previous = TRUE,
      ref_qc_types = "BQC",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "-4.30% to 3.10%)",
    fixed = TRUE)

  expect_message(
    mexp_batch2 <- correct_batch_centering(
      mexp_batch1,
      correct_scale = FALSE,
      replace_previous = TRUE,
      ref_qc_types = "BQC",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "Replacing previous `conc` batch correction",
    fixed = TRUE)

  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp,
      correct_scale = FALSE,
      replace_previous = FALSE,
      ref_qc_types = "BQC",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "-7.00% to 2.70%)",
    fixed = TRUE)

  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp,
      correct_scale = FALSE,
      replace_previous = FALSE,
      ref_qc_types = "BQC",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "Adding batch correction to `conc` data",
    fixed = TRUE)

  expect_message(
    mexp_batch2 <- correct_batch_centering(
      mexp_batch1,
      correct_scale = FALSE,
      replace_previous = FALSE,
      ref_qc_types = "TQC",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "Adding batch correction on top of previous `conc` batch correction.",
    fixed = TRUE)

  expect_message(
    mexp_batch2 <- correct_batch_centering(
      mexp_batch1,
      correct_scale = FALSE,
      replace_previous = TRUE,
      ref_qc_types = "TQC",
      variable = "conc",
      replace_exisiting_trendcurves = FALSE),
    "Replacing previous `conc` batch correction.",
    fixed = TRUE)

})


test_that("correct_batch_centering invalidates downstream states when correcting upstream variable", {

  # add on top of previous with only drift correction before is same result as replace previous
  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      replace_previous = FALSE,
      ref_qc_types = "BQC",
      variable = "intensity",
      replace_exisiting_trendcurves = FALSE),
    "The normalized intensities and concentrations are no longer valid. Please reprocess the data",
    fixed = TRUE)


  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      replace_previous = FALSE,
      ref_qc_types = "BQC",
      variable = "norm_intensity",
      replace_exisiting_trendcurves = FALSE),
    "Concentrations are no longer valid. Please reprocess the data.",
    fixed = TRUE)

})

test_that("correct_batch_centering handels other errors", {

  # add on top of previous with only drift correction before is same result as replace previous

  mexp_dcorr_tmp <- mexp_dcorr
  mexp_dcorr_tmp@dataset$batch_id = "1"


  expect_error(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr_tmp,
      correct_scale = FALSE,
      replace_previous = FALSE,
      ref_qc_types = "BQC",
      variable = "intensity",
      replace_exisiting_trendcurves = FALSE),
    "Batch correction was not applied as there is only one batch.",
    fixed = TRUE)


  expect_message(
    mexp_batch1 <- correct_batch_centering(
      mexp_dcorr,
      correct_scale = FALSE,
      replace_previous = FALSE,
      ref_qc_types = "BQC",
      variable = "norm_intensity",
      replace_exisiting_trendcurves = FALSE),
    "Concentrations are no longer valid. Please reprocess the data.",
    fixed = TRUE)

})

test_that("fun_batch.correction handles non log setting when batch scaling", {
  expect_error(
    fun_batch.correction(tibble(x = 1:10, y = 1:10, batch_id = 1, y_fit_after = 1:10, qc_type = "BQC"),
                     log_transform_internal = FALSE, ref_qc_types = "BQC", correct_scale = TRUE),
    "Currently data must be log-transformed for batch scaling")
})
