library(vdiffr)
library(ggplot2)

mexp <- lipidomics_dataset
mexp_nonorm <- mexp
mexp <- normalize_by_istd(mexp)
mexp <- calc_qc_metrics(mexp)  # Ensure calc_qc_metrics is executed before

test_that("plot_qcmetrics_comparison plot looks as expected", {
  p <- plot_qcmetrics_comparison(
    data = mexp,
    x_variable = "precursor_mz",
    y_variable = "rt_median_SPL",
    equality_line = TRUE,
    facet_by_class = TRUE
  )
  expect_doppelganger("default plot_qcmetrics_comparison plot", p)
})


test_that("plot_qcmetrics_comparison plot looks as expected", {
  p <- plot_normalization_qc(
    data = mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    qc_type = "BQC",
    facet_by_class = TRUE,
  )
  expect_doppelganger("default plot_normalization_qc plot", p)
})


test_that("plot_qcmetrics_comparison() generates a plot", {

  # Test with valid arguments
  p <- plot_qcmetrics_comparison(
        data = mexp,
        x_variable = "precursor_mz",
        y_variable = "rt_median_SPL",
        facet_by_class = TRUE
    )
  # Check if a ggplot object is returned
  expect_s3_class(p, "gg")

  # Extract the plot's data (data frame used for the first layer)
  plot_data <- ggplot_build(p)$data[[1]]
    # Test if the number of points in the plot matches the expected value
  expect_equal(nrow(plot_data), 29)

  # Test with non-existing variables
  expect_error(
    p <- plot_qcmetrics_comparison(
      data = mexp,
      x_variable = "unknownA",
      y_variable = "unknownB",
      facet_by_class = FALSE
    ),
    "One or both of the specified "
  )

  # Check if a ggplot object is returned
  expect_s3_class(p, "gg")


  mexp_no_class <- mexp
  mexp_no_class@metrics_qc$feature_class <- NA  # Simulate missing feature class
  expect_no_error(
    plot_qcmetrics_comparison(
      data = mexp_no_class,
      x_variable = "precursor_mz",
      y_variable = "rt_median_SPL",
      facet_by_class = FALSE
    )
  )

  expect_error(
    plot_qcmetrics_comparison(
      data = mexp_no_class,
      x_variable = "precursor_mz",
      y_variable = "rt_median_SPL",
      facet_by_class = TRUE
    ),
    "`feature_class` to be defined in the metadata"
  )

  expect_error(
    plot_qcmetrics_comparison(
      data = mexp,
      x_variable = "precursor_mz",
      y_variable = "rt_median_SPL",
      facet_by_class = FALSE,
      filter_data = TRUE

    ),
    "Data has not yet been QC-filtered"
  )

  mexp_filt <- mexp
  mexp_filt <- filter_features_qc(mexp_filt,
                                  include_qualifier = TRUE, include_istd = FALSE,
                                  max.cv.intensity.bqc = 25)


  # Check if qc filter data works
  expect_no_error(
    p <- plot_qcmetrics_comparison(
      data = mexp_filt,
      x_variable = "precursor_mz",
      y_variable = "rt_median_SPL",
      facet_by_class = FALSE,
      filter_data = TRUE

    )
  )
  # Extract the plot's data (data frame used for the first layer)
  plot_data <- ggplot_build(p)$data[[1]]
  # Test if the number of points in the plot matches the expected value
  expect_equal(nrow(plot_data), 19)


  # Check if include_qualifier = FALSE works
  p <- plot_qcmetrics_comparison(
    data = mexp_filt,
    x_variable = "precursor_mz",
    y_variable = "rt_median_SPL",
    facet_by_class = FALSE,
    filter_data = TRUE, include_qualifier = FALSE
  )

  # Extract the plot's data (data frame used for the first layer)
  plot_data <- ggplot_build(p)$data[[1]]
  # Test if the number of points in the plot matches the expected value
  expect_equal(nrow(plot_data), 18)

})


test_that("plot_normalization_qc() generates a plot", {

  # Test with valid arguments
  p <- plot_normalization_qc(
    data = mexp,
    before_norm_var = "intensity",
    after_norm_var = "norm_intensity",
    qc_type = "BQC",
    facet_by_class = TRUE,
  )
  # Check if a ggplot object is returned
  expect_s3_class(p, "gg")

  # Extract the plot's data (data frame used for the first layer)
  plot_data <- ggplot_build(p)$data[[1]]
  # Test if the number of points in the plot matches the expected value
  expect_equal(nrow(plot_data), 28)

  expect_error(
    plot_normalization_qc(
      data = mexp,
      before_norm_var = "conc",
      after_norm_var = "norm_intensity",
      qc_type = "BQC",
      facet_by_class = TRUE,
    ),
    "`before_norm_var` must be one of"
  )
  expect_error(
    plot_normalization_qc(
      data = mexp,
      before_norm_var = "norm_intensity",
      after_norm_var = "intensity",
      qc_type = "BQC",
      facet_by_class = TRUE,
    ),
    "`after_norm_var` must be one of"
  )

  expect_error(
    plot_normalization_qc(
      data = mexp,
      before_norm_var = "norm_intensity",
      after_norm_var = "norm_intensity",
      qc_type = "BQC",
      facet_by_class = TRUE,
    ),
    "`before_norm_var` and `after_norm_var` cannot be the same"
  )

  expect_error(
    plot_normalization_qc(
      data = mexp,
      before_norm_var = "intensity",
      after_norm_var = "norm_intensity",
      qc_type = "UNdefind",
      facet_by_class = TRUE,
    ),
    "`qc_type` must be one of"
  )

})

