library(vdiffr)
library(ggplot2)

mexp <- lipidomics_dataset
mexp_nonorm <- mexp
mexp <- normalize_by_istd(mexp)
mexp <- calc_qc_metrics(mexp) # Ensure calc_qc_metrics is executed before
mexp <- filter_features_qc(
  mexp,
  include_qualifier = TRUE,
  include_istd = FALSE,
  max.cv.intensity.bqc = 10
)

test_that("plot_normalization_qc() generates a plot", {
  # Test with valid arguments
  p <- plot_normalization_qc(
    data = mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    plot_type = "diff",
    qc_type = NA,
    facet_by_class = TRUE
  )
  # Check if a ggplot object is returned
  expect_s3_class(p, "gg")

  # Extract the plot's data (data frame used for the first layer)
  plot_data <- ggplot2::ggplot_build(p)$data[[1]]
  # Test if the number of points in the plot matches the expected value
  expect_equal(nrow(plot_data), 9)

  expect_error(
    plot_normalization_qc(
      data = mexp,
      before_norm_var = "conc",
      after_norm_var = "norm_intensity",
      plot_type = "diff",
      qc_type = "BQC",
      facet_by_class = TRUE
    ),
    "`before_norm_var` must be one of"
  )
  expect_error(
    plot_normalization_qc(
      data = mexp,
      before_norm_var = "norm_intensity",
      after_norm_var = "intensity",
      plot_type = "diff",
      qc_type = "BQC",
      facet_by_class = TRUE
    ),
    "`after_norm_var` must be one of"
  )

  expect_error(
    plot_normalization_qc(
      data = mexp,
      before_norm_var = "conc_raw",
      after_norm_var = "conc_raw",
      plot_type = "diff",
      qc_type = "BQC",
      facet_by_class = TRUE
    ),
    "`before_norm_var` and `after_norm_var` cannot be the same",
    fixed = TRUE
  )

  expect_error(
    plot_normalization_qc(
      data = mexp,
      before_norm_var = "conc_raw",
      after_norm_var = "intensity",
      qc_type = "BQC",
      plot_type = "diff",
      facet_by_class = TRUE
    ),
    "`after_norm_var` must be one of",
    fixed = TRUE
  )

  expect_error(
    plot_normalization_qc(
      data = mexp,
      before_norm_var = "intensity",
      after_norm_var = "norm_intensity",
      qc_type = "UNdefind",
      facet_by_class = TRUE,
      plot_type = "diff"
    ),
    "One or more specified `qc_types` are not present in the dataset",
    fixed = TRUE
  )
})

# tests/testthat/test-plot-qc-comparisons.R

# --- Visual Regression Tests for plot_normalization_qc ---

test_that("plot_normalization_qc generates correct plot types", {
  p_scatter <- plot_normalization_qc(
    data = mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    plot_type = "scatter",
    qc_types = "BQC"
  )
  vdiffr::expect_doppelganger("norm-qc-scatter", p_scatter)

  p_diff <- plot_normalization_qc(
    data = mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    plot_type = "diff",
    qc_types = "BQC"
  )
  vdiffr::expect_doppelganger("norm-qc-diff", p_diff)

  p_ratio <- plot_normalization_qc(
    data = mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    plot_type = "ratio",
    qc_types = "BQC"
  )
  vdiffr::expect_doppelganger("norm-qc-ratio", p_ratio)
})

test_that("plot_normalization_qc faceting works", {
  p_faceted <- plot_normalization_qc(
    data = mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    plot_type = "scatter",
    qc_types = "SPL",
    facet_by_class = TRUE
  )
  vdiffr::expect_doppelganger("norm-qc-faceted", p_faceted)

  # Object check for faceting
  built_plot <- ggplot_build(p_faceted)
  expect_gt(length(built_plot$layout$panel_params), 1)
})


test_that("plot_normalization_qc pre-flight checks for data state work", {
  # Test on data that has not been normalized
  expect_error(
    plot_normalization_qc(
      mexp_nonorm,
      plot_type = "diff",
      before_norm_var = "intensity",
      after_norm_var = "norm_intensity",
      qc_types = "BQC"
    ),
    "Data has not yet been normalized",
    fixed = TRUE
  )

  # Test on data without QC metrics
  mexp_no_metrics <- mexp
  mexp_no_metrics@metrics_qc <- mexp_no_metrics@metrics_qc[0, ] # Empty the table
  expect_error(
    plot_normalization_qc(
      mexp_no_metrics,
      before_norm_var = "intensity",
      after_norm_var = "norm_intensity",
      plot_type = "diff"
    ),
    "No QC metrics available yet",
    fixed = TRUE
  )

  expect_error(
    plot_normalization_qc(
      mexp_no_metrics,
      after_norm_var = "norm_intensity",
      plot_type = "diff"
    ),
    "`before_norm_var` and `after_norm_var` must be supplied",
    fixed = TRUE
  )
  expect_error(
    plot_normalization_qc(
      mexp_no_metrics,
      before_norm_var = "norm_intensity",
      plot_type = "diff"
    ),
    "`before_norm_var` and `after_norm_var` must be supplied",
    fixed = TRUE
  )
  expect_error(
    plot_normalization_qc(
      mexp_no_metrics,
      before_norm_var = "intensity",
      after_norm_var = "norm_intensity"
    ),
    "`plot_type` must be supplied ('scatter', 'diff', or 'ratio')",
    fixed = TRUE
  )
})

test_that("plot_normalization_qc constructs correct axis labels", {
  # This object check verifies the internal logic without a visual diff
  p_diff <- plot_normalization_qc(
    mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    plot_type = "diff",
    qc_types = "BQC"
  )

  built_plot <- ggplot_build(p_diff)

  expect_equal(
    built_plot$plot$scales$get_scales("y")$name,
    "norm_intensity_cv - intensity_cv"
  )
  expect_equal(
    built_plot$plot$scales$get_scales("x")$name,
    "Mean of intensity_cv and norm_intensity_cv"
  )

  p_ratio <- plot_normalization_qc(
    mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    plot_type = "ratio",
    qc_types = "BQC"
  )
  built_plot <- ggplot_build(p_ratio)
  expect_equal(
    built_plot$plot$scales$get_scales("y")$name,
    "log2( norm_intensity_cv / intensity_cv )"
  )

  p_ratio <- plot_normalization_qc(
    mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    plot_type = "scatter",
    qc_types = "BQC"
  )
  built_plot <- ggplot_build(p_ratio)
  expect_equal(
    built_plot$plot$scales$get_scales("y")$name,
    "QC metric: norm_intensity_cv"
  )
})


# --- Tests for plot_qcmetrics_comparison ---

test_that("plot_qcmetrics_comparison plot looks as expected", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "scatter",
    y_shared = TRUE,
    equality_line = FALSE,
    facet_by_class = TRUE
  )
  vdiffr::expect_doppelganger("default plot_qcmetrics_comparison plot", p)
})

test_that("plot_qcmetrics_comparison plot no facets", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "scatter",
    y_shared = TRUE,
    equality_line = FALSE,
    facet_by_class = FALSE
  )
  vdiffr::expect_doppelganger("nofacet plot_qcmetrics_comparison plot", p)
})

# this comparison doesnt make sense, but it tests the plotting function
test_that("plot_qcmetrics_comparison plot looks as expected", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "diff",
    y_shared = TRUE,
    equality_line = FALSE,
    facet_by_class = TRUE
  )
  vdiffr::expect_doppelganger("diff plot_qcmetrics_comparison plot", p)
})

# this comparison doesnt make sense, but it tests the plotting function
test_that("plot_qcmetrics_comparison plot looks as expected", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "ratio",
    y_shared = TRUE,
    equality_line = FALSE,
    facet_by_class = TRUE
  )
  vdiffr::expect_doppelganger("plot_qcmetrics_comparison ratio ", p)
})

test_that("plot_qcmetrics_comparison plot looks as expected", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "scatter",
    y_shared = FALSE,
    equality_line = FALSE,
    filter_data = FALSE,
    include_qualifier = FALSE,
    facet_by_class = TRUE
  )
  # Extract the plot's data (data frame used for the first layer)
  plot_data <- ggplot2::ggplot_build(p)$data[[1]]
  # Test if the number of points in the plot matches the expected value
  expect_equal(nrow(plot_data), 19)
  vdiffr::expect_doppelganger("plot_qcmetrics_comparison 1panel", p)
})

test_that("plot_qcmetrics_comparison plot looks as expected", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "scatter",
    y_shared = FALSE,
    equality_line = FALSE,
    filter_data = FALSE,
    include_qualifier = TRUE,
    facet_by_class = TRUE
  )
  # Extract the plot's data (data frame used for the first layer)
  plot_data <- ggplot2::ggplot_build(p)$data[[1]]
  # Test if the number of points in the plot matches the expected value
  expect_equal(nrow(plot_data), 20)
})

test_that("plot_qcmetrics_comparison data filter works", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "scatter",
    filter_data = TRUE,
    y_shared = FALSE,
    equality_line = FALSE,
    facet_by_class = TRUE
  )
  # Extract the plot's data (data frame used for the first layer)
  plot_data <- ggplot2::ggplot_build(p)$data[[1]]
  # Test if the number of points in the plot matches the expected value
  expect_equal(nrow(plot_data), 6)
  vdiffr::expect_doppelganger("plot_qcmetrics_comparison filter", p)
})

test_that("plot_qcmetrics_comparison plot threshold", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "scatter",
    y_shared = FALSE,
    threshold_value = 10,
    facet_by_class = TRUE
  )
  vdiffr::expect_doppelganger("plot_qcmetrics_comparison threshold1", p)
})

test_that("plot_qcmetrics_comparison plot log", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "scatter",
    y_shared = FALSE,
    threshold_value = c(NA, 10),
    facet_by_class = TRUE
  )
  vdiffr::expect_doppelganger("plot_qcmetrics_comparison threshold2", p)
})


test_that("plot_qcmetrics_comparison plot looks as expected", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "rt_median_bqc",
    y_variable = "norm_intensity_cv_bqc",
    plot_type = "scatter",
    y_shared = FALSE,
    threshold_value = c(NA, 10),
    facet_by_class = TRUE,
    log_scale = TRUE
  )
  vdiffr::expect_doppelganger("plot_qcmetrics_comparison thresholdlog", p)
})


test_that("plot_qcmetrics_comparison error handling works", {
  expect_error(
    p <- plot_qcmetrics_comparison(
      data = mexp,
      x_variable = "rt_median_bqc",
      y_variable = "norm_intensity_cv_bqc",
      plot_type = "scatter",
      y_shared = FALSE,
      y_lim = c(0, 10),
      facet_by_class = TRUE,
      log_scale = TRUE
    ),
    "Log scale cannot be used with zero, negative, infinite, or NA axis limits",
    fixed = TRUE
  )

  # Test faceting when feature_class column is missing
  mexp_no_class <- mexp
  mexp_no_class@metrics_qc$feature_class <- NULL
  expect_error(
    p <- plot_qcmetrics_comparison(
      data = mexp_no_class,
      x_variable = "rt_median_bqc",
      y_variable = "norm_intensity_cv_bqc",
      plot_type = "scatter",
      y_shared = FALSE,
      facet_by_class = TRUE,
      log_scale = TRUE
    ),
    "`feature_class` to be defined in the metadata",
    fixed = TRUE
  )
})


test_that("get_feature_correlations returns correct filtered long format", {
  # Create example data
  set.seed(123)
  df <- tibble(
    analysis_id = paste0("sample", 1:5),
    qc_type = rep("BQC", 5),
    feat1 = rnorm(5),
    feat2 = rnorm(5),
    feat3 = rnorm(5)
  )

  # Run the function with thresholds
  res <- get_feature_correlations(df, cor_min_neg = -0.2, cor_min = 0.2)

  # Check that result is a tibble/data.frame
  expect_s3_class(res, "data.frame")

  # Check required columns exist
  expect_true(all(c("var1", "var2", "value") %in% colnames(res)))

  # Check only upper triangle correlations are included
  expect_true(all(res$var1 < res$var2))

  # Check thresholds
  expect_true(all(res$value <= -0.2 | res$value >= 0.2))

  # Check number of rows is <= ncol(numeric choose 2)
  numeric_cols <- df |> select(where(is.numeric)) |> colnames()
  n_pairs <- choose(length(numeric_cols), 2)
  expect_lte(nrow(res), n_pairs)
})

test_that("function works when no correlations pass thresholds", {
  set.seed(42)
  df <- tibble(
    analysis_id = paste0("sample", 1:3),
    qc_type = "BQC",
    x = c(1, 2, 3),
    y = c(4, 5, 6)
  )

  # Use high thresholds to filter everything out
  res <- get_feature_correlations(df, cor_min_neg = -1, cor_min = 2)

  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 0) # no correlations pass thresholds
})
