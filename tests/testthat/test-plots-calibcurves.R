library(fs)
library(vdiffr)
library(ggplot2)
library(testthat)
library(scales)

mexp <- quant_lcms_dataset

mexp <- normalize_by_istd(mexp)
mexp <- calc_calibration_results(mexp, fit_model = "quadratic", fit_weighting = "1/x")

test_that("plot_responsecurves generates a plot", {

  # Test with valid arguments
  p <- plot_calibrationcurves(
    data = mexp,
    fit_model = "quadratic",
    fit_weighting = "1/x",
    rows_page = 2,
    cols_page = 2,
    return_plots = TRUE
  )
  expect_s3_class(p[[1]], "gg")
  # Check how many pages
  expect_equal(length(p), 2)
  expect_doppelganger("default plot_calibration plot 1", p[[1]])

  p <- plot_calibrationcurves(
    data = mexp,
    fit_model = "quadratic",
    show_confidence_interval = FALSE,
    fit_weighting = "1/x",
    rows_page = 2,
    cols_page = 2,
    return_plots = TRUE
  )
  expect_s3_class(p[[1]], "gg")
  # Check how many pages
  expect_equal(length(p), 2)
  expect_doppelganger("no ci plot_calibration plot ", p[[1]])

  expect_no_error(
    p <- plot_calibrationcurves(
      data = mexp,
      log_axes = TRUE,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      return_plots = TRUE
    )
  )
  expect_doppelganger("log-log plot_calibration plot default ", p[[1]])

  expect_message(
    p <- plot_calibrationcurves(
      data = mexp,
      log_axes = TRUE,
      fit_model = "quadratic",
      show_confidence_interval = TRUE,
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      return_plots = TRUE
  ), "Regions of the regression confidence intervals are partially")
    expect_doppelganger("log-log plot_calibration plot with ci ", p[[1]])


    mexp_temp <- mexp
    mexp_temp@dataset<- mexp_temp@dataset|>
      mutate(feature_norm_intensity = if_else(str_detect(analyte_id, "Cortiso") & analysis_id =="CalA", 0.000001, feature_norm_intensity))
    mexp_temp@dataset<- mexp_temp@dataset|>
      mutate(feature_norm_intensity = if_else(str_detect(analyte_id, "Cortiso") & analysis_id =="CalB", 0.000001, feature_norm_intensity))


    # Test with valid arguments
    expect_message(
      p <- plot_calibrationcurves(
        data = mexp_temp,
        overwrite_fit_param = TRUE,
        fit_model = "linear",
        fit_weighting = "1/x",
        log_axes = TRUE,
        rows_page = 2,
        cols_page = 2,
        return_plots = TRUE
      ),
      "Regions of the regression curve are partially")

    # Test with valid arguments
    expect_message(
      p <- plot_calibrationcurves(
        data = mexp_temp,
        overwrite_fit_param = TRUE,
        show_confidence_interval = TRUE,
        fit_model = "linear",
        fit_weighting = "1/x",
        log_axes = TRUE,
        rows_page = 2,
        cols_page = 2,
        return_plots = TRUE
      ),
      "Regions of the regression curve and confidence intervals are partially")


  # Test if the number of points in the plot matches the expected value
  plot_data <- ggplot_build(p[[1]])$data[[1]]
  expect_equal(nrow(plot_data), 400)

  temp_pdf_path <- file.path(tempdir(), "midar_test_calibcurve.pdf")
  p <- plot_calibrationcurves(
    data = mexp,
    fit_model = "quadratic",
    fit_weighting = "1/x",
    rows_page = 2,
    cols_page = 2,
    output_pdf = TRUE,
    path = temp_pdf_path,
    return_plots = FALSE
  )
  expect_null(p)
  expect_true(file_exists(temp_pdf_path), info = "PDF file was not created.")
  expect_equal(as.character(fs::file_size(temp_pdf_path)), "29K")
  fs::file_delete(temp_pdf_path)

  temp_pdf_path <- file.path(tempdir(), "midar_test_responsecurve.pdf")
  expect_no_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      output_pdf = TRUE,
      path = temp_pdf_path,
      return_plots = FALSE
    )
  )

  expect_silent(p)
  expect_true(file_exists(temp_pdf_path), info = "PDF file was not created.")
  expect_equal(as.character(fs::file_size(temp_pdf_path)), "29K")
  fs::file_delete(temp_pdf_path)

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      output_pdf = FALSE,
      specific_page = 3,
      path = temp_pdf_path,
      return_plots = FALSE
    ),
    "Selected page exceeds "
  )
  p <- plot_calibrationcurves(
    data = mexp,
    fit_model = "quadratic",
    fit_weighting = "1/x",
    rows_page = 2,
    cols_page = 2,
    output_pdf = FALSE,
    specific_page = 2,
    path = temp_pdf_path,
    return_plots = TRUE
  )
  expect_equal(length(p), 1)

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      output_pdf = TRUE,
      specific_page = 2,
      return_plots = TRUE
    ),
    "The argument "
  )

})


test_that("plot_responsecurves handles missing data", {

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      output_pdf = TRUE,
      return_plots = TRUE
    ),
    "The argument "
  )

  mexp_defect <- mexp
  mexp_defect@dataset <- mexp_defect@dataset |> dplyr::slice_head(n = 0)

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp_defect,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "No data available in "
  )

  mexp_defect <- mexp
  mexp_defect@dataset$qc_type <- "SPL"

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp_defect,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "No QC type "
  )

  mexp_defect <- mexp
  mexp_defect@annot_qcconcentrations <- mexp_defect@annot_qcconcentrations  |> dplyr::slice_head(n = 0)

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp_defect,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "No QC-concentration metadata"
  )
})


test_that("curve color definition works", {
  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      point_color = c("red","green"),
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "Insufficient colors in \\`point_colors\\`"
  )

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      point_fill = c("red","green"),
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "Insufficient fill colors in \\`point_fill\\`"
  )

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      point_shape = c(1,2),
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "Insufficient shape codes in \\`point_shape\\`"
  )

  expect_no_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      point_fill = c("red","green", "blue"),
      output_pdf = FALSE,
      return_plots = TRUE
    ))
  p_data <- ggplot_build(p[[1]])$data
  expect_equal(unique(p_data[[4]]$fill), c("red", "green", "blue"))



  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      point_color = c("CAL" = "red","green", "blue"),
      point_fill = c("CAL" = "red","green", "blue"),
      point_shape = c("CAL" = 22,21, 23),

            output_pdf = FALSE,
      return_plots = TRUE
    ),
    "The names in \\`point_color\\`")

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      point_fill = c("CAL" = "red","green", "blue"),
      point_shape = c("CAL" = 22,21, 23),

      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "The names in \\`point_fill\\`")

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      point_shape = c("CAL" = 22,21, 23),

      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "The names in \\`point_shape\\`")

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      point_color = c("CAL" = "red","QC"="green", "LQC"="blue"),
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "The names in \\`point_color\`")


  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      qc_types = c("CAL", "QC", "LQC"),
      rows_page = 2,
      cols_page = 2,
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "One or more specified \\`qc_types\\`")

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      qc_types = c("CAL", "SPL", "LQC"),
      rows_page = 2,
      cols_page = 2,
      output_pdf = FALSE,
      return_plots = TRUE
    ),
    "One or more selected \\`qc_types\\`")

  expect_error(
    p <- plot_calibrationcurves(
      data = mexp,
      filter_data = TRUE,
      fit_model = "quadratic",
      fit_weighting = "1/x",
      rows_page = 2,
      cols_page = 2,
      return_plots = TRUE
    ),
    "Data has not been QC-filtered")

})


test_that("plot_responsecurves generates a plot with calib failes", {

  mexp_temp <- mexp
  mexp_temp@annot_qcconcentrations <- mexp_temp@annot_qcconcentrations |>
    mutate(concentration = if_else(str_detect(analyte_id, "Cortiso") & str_detect(sample_id, "CAL"), NA_real_, concentration))

  # Test with valid arguments
  expect_message(
  p <- plot_calibrationcurves(
    data = mexp_temp,
    fit_model = "quadratic",
    fit_weighting = "1/x",
    rows_page = 2,
    cols_page = 2,
    return_plots = TRUE
  ),
  "Regression failed for 4 features")
  expect_s3_class(p[[1]], "gg")
  expect_equal(length(p), 2)
  expect_doppelganger("default plot_calibration plot log_axes 1", p[[1]])


  p <- plot_calibrationcurves(
    data = mexp_temp,
    fit_model = "quadratic",
    fit_weighting = "1/x",
    log_axes = TRUE,
    rows_page = 2,
    cols_page = 2,
    return_plots = TRUE
  )
  expect_s3_class(p[[1]], "gg")
  expect_equal(length(p), 2)
  expect_doppelganger("default plot_calibration plot log_axes 2", p[[1]])

})


