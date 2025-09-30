library(testthat)
library(ggplot2)
library(dplyr)
library(cli)
library(stringr) # Added for tests


# --- 1. Setup  ---

mexp_orig <- lipidomics_dataset
mexp <- exclude_analyses(
  mexp_orig,
  c("Longit_batch6_51", "Longit_batch6_B-ISTD 09"),
  clear_existing = TRUE
)
mexp <- normalize_by_istd(mexp)
mexp <- quantify_by_istd(mexp)
mexp <- calc_qc_metrics(mexp) # Ensure calc_qc_metrics is executed before
mexp <- filter_features_qc(
  mexp,
  include_qualifier = FALSE,
  include_istd = FALSE,
  max.cv.conc.bqc = 13
)

# --- 2. Start Testing ---

test_that("Core functionality works correctly", {
  p_default <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = TRUE
  )
  expect_s3_class(p_default, "ggplot")
  expect_silent(ggplot_build(p_default))
  vdiffr::expect_doppelganger("plot_abundanceprofile-default", p_default)

  p_qc <- plot_abundanceprofile(
    data = mexp,
    variable = "conc_median_spl",
    qc_types = "SPL",
    log_scale = TRUE,
    use_qc_metrics = TRUE
  )
  expect_s3_class(p_qc, "ggplot")
  expect_silent(ggplot_build(p_qc))
})

test_that("Argument validation and errors are handled", {
  expect_error(
    plot_abundanceprofile(
      data = mexp,
      variable = "invalid_var",
      qc_types = "SPL",
      log_scale = TRUE
    ),
    "must be one of"
  )
  expect_error(
    plot_abundanceprofile(
      data = mexp,
      variable = "invalid_var",
      qc_types = "SPL",
      log_scale = TRUE,
      use_qc_metrics = TRUE
    ),
    "Could not find a summary column"
  )
  expect_warning(
    plot_abundanceprofile(
      data = mexp,
      variable = "conc_median_spl",
      qc_types = c("SPL", "BQC"),
      log_scale = TRUE,
      use_qc_metrics = TRUE
    ),
    "will be ignored"
  )
})

test_that("Filtering flags work correctly", {
  # Test `include_istd`
  p_no_istd <- plot_abundanceprofile(
    data = mexp,
    variable = "intensity",
    qc_types = "SPL",
    log_scale = TRUE,
    include_istd = FALSE
  )
  d_no_istd <- ggplot_build(p_no_istd)$data[[2]] # geom_segment data
  expect_equal(nrow(d_no_istd), 19) # F7 (ISTD) should be removed

  p_with_istd <- plot_abundanceprofile(
    data = mexp,
    variable = "intensity",
    qc_types = "SPL",
    log_scale = TRUE,
    include_istd = TRUE
  )
  d_with_istd <- ggplot_build(p_with_istd)$data[[2]] # geom_segment data
  expect_equal(nrow(d_with_istd), 28) # F7 (ISTD) should be removed

  p_no_istd_qc <- plot_abundanceprofile(
    data = mexp,
    variable = "intensity_median_spl",
    qc_types = "SPL",
    log_scale = TRUE,
    use_qc_metrics = TRUE,
    include_istd = FALSE
  )
  d_no_istd_qc <- ggplot_build(p_no_istd_qc)$data[[2]]
  expect_equal(nrow(d_no_istd_qc), 19) # F7 (ISTD) should be removed

  p_with_istd <- plot_abundanceprofile(
    data = mexp,
    variable = "intensity_median_spl",
    qc_types = "SPL",
    log_scale = TRUE,
    use_qc_metrics = TRUE,
    include_istd = TRUE
  )
  d_with_istd <- ggplot_build(p_with_istd)$data[[2]]
  expect_equal(nrow(d_with_istd), 28) # F7 (ISTD) should be removed

  # Test `include_qualifier`
  p_with_qual <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = TRUE,
    include_qualifier = TRUE
  )
  d_with_qual <- ggplot_build(p_with_qual)$data[[2]]
  expect_equal(nrow(d_with_qual), 20) # F4 (qualifier) should be removed

  p_with_qual <- plot_abundanceprofile(
    data = mexp,
    variable = "conc_median_spl",
    qc_types = "SPL",
    log_scale = TRUE,
    use_qc_metrics = TRUE,
    include_qualifier = TRUE
  )
  d_with_qual <- ggplot_build(p_with_qual)$data[[2]]
  expect_equal(nrow(d_with_qual), 20) # F4 (qualifier) should be removed

  # Test `exclude qualifier`
  p_no_qual <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = TRUE,
    include_qualifier = FALSE
  )
  d_no_qual <- ggplot_build(p_no_qual)$data[[2]]
  expect_equal(nrow(d_no_qual), 19) # F4 (qualifier) should be removed

  p_no_qual_qc <- plot_abundanceprofile(
    data = mexp,
    variable = "conc_median_spl",
    qc_types = "SPL",
    log_scale = TRUE,
    use_qc_metrics = TRUE,
    include_qualifier = FALSE
  )
  d_no_qual_qc <- ggplot_build(p_no_qual_qc)$data[[2]]
  expect_equal(nrow(d_no_qual_qc), 19) # F4 (qualifier) should be removed

  # Test `exclude_classes`
  p_exclude <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = TRUE,
    exclude_classes = "TG"
  )
  labels_exclude <- ggplot_build(p_exclude)$layout$panel_params[[
    1
  ]]$y$get_labels()
  expect_false("TG" %in% labels_exclude)

  # Test `filter_data` with use_qc_metrics = TRUE
  p_filter_qc <- plot_abundanceprofile(
    data = mexp,
    variable = "conc_median_spl",
    qc_types = "SPL",
    log_scale = TRUE,
    use_qc_metrics = TRUE,
    filter_data = TRUE
  )
  d_filter_qc <- ggplot_build(p_filter_qc)$data[[2]]
  expect_equal(nrow(d_filter_qc), 6) # F3 (all_filter_pass=FALSE) should be removed
})

test_that("Feature map and y-axis logic is correct", {
  # Test string shortcut
  p_shortcut <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = TRUE,
    feature_map = "lipidomics"
  )
  expect_s3_class(p_shortcut, "ggplot")

  # Test `drop_empty_classes` = FALSE
  p_keep <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = TRUE,
    drop_empty_classes = FALSE
  )
  labels_keep <- ggplot_build(p_keep)$layout$panel_params[[1]]$y$get_labels()
  expect_true("LPE" %in% labels_keep)

  # Test `drop_empty_classes` = TRUE (default)
  p_drop <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = TRUE,
    drop_empty_classes = TRUE
  )
  labels_drop <- ggplot_build(p_drop)$layout$panel_params[[1]]$y$get_labels()
  expect_false("LPE" %in% labels_drop)
})

test_that("X-axis limit handling and warnings work correctly", {
  # Set limits that will exclude one feature (the LPC at conc=100)
  p_warn <- expect_message(
    plot_abundanceprofile(
      data = mexp,
      variable = "conc",
      qc_types = "SPL",
      log_scale = FALSE,
      x_lim = c(200, 3000)
    ),
    "Some data points fall outside the `x_lim` range and were removed: 19 features and 9 class sums",
  )

  # Check that ggplot does NOT issue its own warning due to pre-filtering
  p_warn <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = FALSE,
    x_lim = c(200, 3000)
  )
  expect_no_warning(ggplot_build(p_warn))
})

test_that("Dynamic show_sum and x_label logic works", {
  # `show_sum` should be FALSE for "rt" variable
  p_rt <- plot_abundanceprofile(
    data = mexp,
    variable = "rt_median_spl",
    qc_types = "SPL",
    log_scale = FALSE,
    use_qc_metrics = TRUE
  )
  build_rt <- ggplot_build(p_rt)
  # Should only have geom_rect and geom_segment layers
  expect_length(build_rt$data, 2)
  expect_equal(p_rt$labels$x, "Retention Time")

  # `x_label` should be specific for "cv"
  p_cv <- plot_abundanceprofile(
    data = mexp,
    variable = "conc_cv_spl",
    qc_types = "SPL",
    log_scale = FALSE,
    use_qc_metrics = TRUE
  )
  expect_equal(p_cv$labels$x, "Coefficient of Variation (%)")

  # `x_label` should be specific for "cv"
  p_sb <- plot_abundanceprofile(
    data = mexp,
    variable = "sb_ratio_pblk",
    qc_types = "SPL",
    log_scale = FALSE,
    use_qc_metrics = TRUE
  )
  expect_equal(p_sb$labels$x, "Signal/Blank Ratio")

  p_rqc <- plot_abundanceprofile(
    data = mexp,
    variable = "r2_rqc_A",
    qc_types = "SPL",
    log_scale = FALSE,
    use_qc_metrics = TRUE
  )
  expect_equal(p_rqc$labels$x, "R2 of Response Curve")

  p_dratio <- plot_abundanceprofile(
    data = mexp,
    variable = "conc_dratio_sd_bqc",
    qc_types = "SPL",
    log_scale = FALSE,
    use_qc_metrics = TRUE
  )
  expect_equal(p_dratio$labels$x, "D-ratio")

  expect_equal(p_rqc$labels$x, "R2 of Response Curve")

  p_na <- plot_abundanceprofile(
    data = mexp,
    variable = "missing_conc_prop_spl",
    qc_types = "SPL",
    log_scale = FALSE,
    use_qc_metrics = TRUE
  )
  expect_equal(p_na$labels$x, "Missingness (Proportion of Samples)")
})

test_that("plot_abundanceprofiledensity strip", {
  p_density <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = FALSE,
    density_strip = TRUE
  )
  expect_s3_class(p_density, "patchwork")
  vdiffr::expect_doppelganger("plot_abundanceprofile p_density", p_density$plot)
})

test_that("plot_abundanceprofiledensity strip log", {
  p_density <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = TRUE,
    density_strip = TRUE
  )
  expect_s3_class(p_density, "patchwork")
  vdiffr::expect_doppelganger("plot_abundanceprofile p_density log", p_density$plot)
})

test_that("plot_abundanceprofiledensity strip log range", {
  p_density <- plot_abundanceprofile(
    data = mexp,
    variable = "conc",
    qc_types = "SPL",
    log_scale = TRUE,
    analysis_range = c(100, 120),
    density_strip = TRUE
  )
  expect_s3_class(p_density, "patchwork")
  vdiffr::expect_doppelganger(
    "plot_abundanceprofile p_density rangelog",
    p_density$plot
  )
})
