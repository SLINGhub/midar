# library(fs)
# library(vdiffr)
# library(ggplot2)
# library(testthat)
# library(scales)

mexp <- lipidomics_dataset

mexp <- normalize_by_istd(mexp)
mexp <- calc_qc_metrics(mexp)

test_that("plot_responsecurves generates a plot", {

  # Test with valid arguments
  p <- plot_responsecurves(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_s3_class(p[[1]], "gg")
  # Check how many pages
  expect_equal(length(p), 3)
  vdiffr::expect_doppelganger("default plot_responsecurves plot", p[[1]])

  # Test if the number of points in the plot matches the expected value
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[1]]
  expect_equal(nrow(plot_data), 1920)

  temp_pdf_path <- file.path(tempdir(), "midar_test_responsecurve.pdf")
  p <- plot_responsecurves(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    output_pdf = TRUE,
    path = temp_pdf_path,
    return_plots = FALSE
  )
  expect_null(p)
  expect_true(file_exists(temp_pdf_path), info = "PDF file was not created.")
  expect_equal(as.character(fs::file_size(temp_pdf_path)), "66K")
  fs::file_delete(temp_pdf_path)

  temp_pdf_path <- file.path(tempdir(), "midar_test_responsecurve.pdf")
  p <- plot_responsecurves(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    output_pdf = TRUE,
    specific_page = 3,
    path = temp_pdf_path,
    return_plots = FALSE
  )

  expect_silent(p)
  expect_true(file_exists(temp_pdf_path), info = "PDF file was not created.")
  expect_equal(as.character(fs::file_size(temp_pdf_path)), "15.4K")
  fs::file_delete(temp_pdf_path)

  expect_error(
    p <- plot_responsecurves(
      data = mexp,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      output_pdf = FALSE,
      specific_page = 4,
      return_plots = TRUE
    ),
    "Selected page exceeds "
  )
  p <- plot_responsecurves(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    output_pdf = FALSE,
    specific_page = 3,
    return_plots = TRUE
  )
  expect_equal(length(p), 1)

  expect_error(
    p <- plot_responsecurves(
      data = mexp,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      output_pdf = TRUE,
      specific_page = 4,
      return_plots = TRUE
    ),
    "The argument "
  )

})


test_that("plot_responsecurves handles missing data", {

  expect_error(
    p <- plot_responsecurves(
      data = mexp,
      variable = "intensity",
      output_pdf = TRUE,
      return_plots = TRUE
    ),
    "The argument "
  )

  mexp_defect <- mexp
  mexp_defect@dataset <- mexp_defect@dataset |> dplyr::slice_head(n = 0)

  expect_error(
    p <- plot_responsecurves(
      data = mexp_defect,
      variable = "intensity",
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "No data available in "
  )

  mexp_defect <- mexp
  mexp_defect@dataset$qc_type <- "SPL"

  expect_error(
    p <- plot_responsecurves(
      data = mexp_defect,
      variable = "intensity",
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "No QC type "
  )

  mexp_defect <- mexp
  mexp_defect@annot_responsecurves <- mexp_defect@annot_responsecurves  |> dplyr::slice_head(n = 0)

  expect_error(
    p <- plot_responsecurves(
      data = mexp_defect,
      variable = "intensity",
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "No response curve metadata is available"
  )


})

test_that("curve color definition works", {
  expect_error(
    p <- plot_responsecurves(
      data = mexp,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      output_pdf = FALSE,
      color_curves = "red",
      specific_page = 3,
      return_plots = TRUE
    ),
  "Insufficient colors in `color_curves`. Provide at least 2 unique"
  )

    p <- plot_responsecurves(
      data = mexp,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      output_pdf = FALSE,
      color_curves = c("red", "blue"),
      specific_page = 3,
      return_plots = TRUE
    )

    smooth_data <- ggplot2::ggplot_build(p[[1]])$data
    expect_equal(unique(smooth_data[[1]]$colour), c("red", "blue"))

    p <- plot_responsecurves(
      data = mexp,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      output_pdf = FALSE,
      color_curves = NA,
      specific_page = 3,
      return_plots = TRUE
    )
    smooth_data <- ggplot2::ggplot_build(p[[1]])$data
    expect_equal(unique(smooth_data[[1]]$colour), c("#34629e","#91bfdb"))

    mexp_temp <- mexp
    mexp_temp@annot_responsecurves$curve_id <- rep(1:6, each = 2)
    p <- plot_responsecurves(
      data = mexp_temp,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      output_pdf = FALSE,
      color_curves = NA,
      specific_page = 3,
      return_plots = TRUE)

    smooth_data <- ggplot2::ggplot_build(p[[1]])$data
    expect_equal(unique(smooth_data[[1]]$colour), c("#F8766D" ,"#B79F00", "#00BA38","#00BFC4" ,"#619CFF" ,"#F564E3"))
})


test_that("`max_regression_value` works", {

  mex_reg_val <- 80

   p <- plot_responsecurves(
    data = mexp,
    variable = "intensity",
    max_regression_value = mex_reg_val,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )

  # Extract the regression data from the ggplot object
  # We need to extract the smooth line data for the regression (method = "lm")
  smooth_data <- ggplot2::ggplot_build(p[[1]])$data[[1]]

  # Check that the 'analyzed_amount' values used for regression are <= max_reg_value
  expect_true(all(smooth_data$x <= 80), info = "Regression data exceeds max_regression_value")
})

test_that("plot_responsecurves feature filters work", {

  # Test with valid arguments
  p <- plot_responsecurves(
    data = mexp,
    variable = "intensity",
    include_feature_filter = "PC",
    exclude_feature_filter = "ISTD",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )

  # Check how many pages
  expect_equal(length(p), 1)

  # Test if the number of points in the plot matches the expected value
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[1]]
  expect_equal(nrow(plot_data), 960)
  expect_equal(mean(plot_data$y), 977995.276)

  expect_error(
    p <- plot_responsecurves(
      data = mexp,
      filter_data = TRUE,
      variable = "norm_intensity",
      include_feature_filter = "PC",
      exclude_feature_filter = "ISTD",
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE
    ),
    "Data has not been QC-filtered"
  )
  mexp_filt <-  mexp
  mexp_filt <- filter_features_qc(mexp_filt, include_qualifier = FALSE, include_istd = FALSE, min.intensity.median.spl = 100000)


  p <- plot_responsecurves(
    data = mexp_filt,
    filter_data = FALSE,
    variable = "norm_intensity",
    include_feature_filter = "PC",
    exclude_feature_filter = "ISTD",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE)

  # Test if the number of points in the plot matches the expected value
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[1]]
  expect_equal(nrow(plot_data), 960)
  expect_equal(mean(plot_data$y), 0.35804775)

  p <- plot_responsecurves(
    data = mexp_filt,
    filter_data = TRUE,
    variable = "norm_intensity",
    include_feature_filter = c("PC 40:6", "PC 40:8"),
    exclude_feature_filter = c("PC 40:1", "PC 40:2"),
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE)

  # Test if the number of points in the plot matches the expected value
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[1]]
  expect_equal(nrow(plot_data), 320)
  expect_equal(mean(plot_data$y), 0.31212162)

  expect_error(
    p <- plot_responsecurves(
      data = mexp_filt,
      filter_data = TRUE,
      variable = "norm_intensity",
      include_feature_filter = c("PC 40:6", "PC 40:8"),
      exclude_feature_filter = "PC",
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE),
    "defined feature filter criteria resulted in no"
  )

  mexp_def <- mexp_filt
  mexp_def@dataset$qc_type <- ifelse(mexp_def@dataset$qc_type == "RQC", "SPL", mexp_def@dataset$qc_type )
  expect_error(
    p <- plot_responsecurves(
      data = mexp_def,
      filter_data = FALSE,
      variable = "norm_intensity",
      include_feature_filter = c("PC 40:6", "PC 40:8"),
      exclude_feature_filter = "PC",
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE),
    "No QC type 'RQC'"
  )

  mexp_def <- mexp_filt
  mexp_def@annot_responsecurves$analysis_id <- paste0(mexp_def@annot_responsecurves$analysis_id, "_no")
  expect_error(
    p <- plot_responsecurves(
      data = mexp_def,
      filter_data = FALSE,
      variable = "norm_intensity",
      include_feature_filter = c("PC 40:6", "PC 40:8"),
      exclude_feature_filter = NA,
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE),
    "Missmatch between data and response curve metadata"
  )
})



