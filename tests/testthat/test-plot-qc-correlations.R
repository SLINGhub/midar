library(vdiffr)
library(ggplot2)

mexp <- lipidomics_dataset


test_that("get_feature_correlations works correctly", {
  # Create test data
  test_data <- data.frame(
    analysis_id = 1:100,
    qc_type = "QC",
    feature1 = rnorm(100),
    feature2 = NA,
    feature3 = rnorm(100)
  )
  test_data$feature2 <- test_data$feature1  # Perfect correlation

  # Test basic functionality
  cors <- get_feature_correlations(test_data, cor_min_neg = -0.9, cor_min = 0.9)
  expect_s3_class(cors, "data.frame")
  expect_equal(names(cors), c("var1", "var2", "value"))
  expect_true(any(cors$value > 0.9))  # Should find high correlation
})

test_that("plot_feature_correlations handles invalid inputs", {
  # Test invalid variable
  expect_error(
    plot_feature_correlations(mexp, variable = "invalid_var"),
    "`variable` must be one of", fixed = TRUE
  )

  # Test invalid correlation thresholds
  expect_error(
    plot_feature_correlations(
      mexp,
      variable = "area",
      cor_min_neg = 0.9,
      cor_min = 0.8
    ),
    "Lower correlation threshold must be less than upper threshold"
  )
})

test_that("plot_feature_correlations handles empty results", {


  # Should return NULL with message
  expect_null(
    plot_feature_correlations(
      mexp,
      variable = "intensity",
      cor_min = 0.99,
      cor_min_neg = -0.99
    )
  )

  expect_message(
    plot_feature_correlations(
      mexp,
      variable = "intensity",
      cor_min = 0.99,
      cor_min_neg = -0.99
    ),
    "No correlations found exceeding thresholds", fixed = TRUE
  )
})

test_that("plot_feature_correlations respects QC types", {

  # Test QC type filtering
  p <- plot_feature_correlations(
    mexp,
    variable = "area",
    qc_types = c("BQC", "SPL", "RQC"),
    cor_min = 0.8,
    cor_min_neg = -0.9, return_plot = TRUE
  )

  # Check that only QC samples are included
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[3]]
  expect_equal(plot_data[1,"label"], "r = 0.971")

  # Test QC type filtering
  p <- plot_feature_correlations(
    mexp,
    variable = "area",
    cor_min = 0.8,
    cor_min_neg = -0.9, return_plot = TRUE
  )

  # Check that only QC samples are included
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[3]]
  expect_equal(plot_data[1,"label"], "r = 0.969")

  vdiffr::expect_doppelganger("default plot_feature_correlations plot", p)

  # Test QC type filtering
  p <- plot_feature_correlations(
    mexp,
    variable = "area",
    cor_min = 0.8,
    cor_min_neg = -0.9, log_scale = TRUE, return_plot = TRUE
  )

  # Check that only QC samples are included
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[3]]
  expect_equal(plot_data[1,"label"], "r = 0.969")

  vdiffr::expect_doppelganger("log plot_feature_correlations plot", p)

  p <- plot_feature_correlations(
    mexp,
    variable = "area",
    cor_min = 0.8, sort_by_corr = FALSE,
    cor_min_neg = -0.9, log_scale = FALSE, return_plot = TRUE
  )

  # Check that only QC samples are included
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[3]]
  expect_equal(plot_data[1,"label"], "r = 0.860")

  p <- plot_feature_correlations(
    mexp,
    variable = "area",
    cor_min = 0.8, sort_by_corr = FALSE, cols_page = 2, rows_page = 2,
    cor_min_neg = -0.9, log_scale = FALSE, return_plot = TRUE
  )

  # Check that only QC samples are included
  plot_data <- ggplot2::ggplot_build(p[[3]])$data[[3]]
  expect_equal(plot_data[1,"label"], "r = 0.949")

  p <- plot_feature_correlations(
    mexp,
    variable = "area",
    cor_min = 0.8, sort_by_corr = FALSE, cols_page = 2, rows_page = 2,specific_page = 3,
    cor_min_neg = -0.9, log_scale = FALSE, return_plot = TRUE
  )

  # Check that only QC samples are included
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[3]]
  expect_equal(plot_data[1,"label"], "r = 0.949")


  expect_error(
    p <- plot_feature_correlations(
      mexp,
      variable = "area",
      cor_min = 0.8, sort_by_corr = FALSE, cols_page = 2, rows_page = 2,specific_page = 4,
      cor_min_neg = -0.9, log_scale = FALSE, return_plot = TRUE
    ),
      "Selected page exceeds the total number of pages"
  )

})


test_that("plot aesthetics are correctly set", {

  # Test custom aesthetics
  p <- plot_feature_correlations(
    mexp,
    variable = "intensity",
    cor_min = 0.85,
    point_size = 2,
    point_alpha = 0.5,
    line_color = "blue",return_plot = TRUE,
    font_base_size = 10
  )

  # Check plot elements
  plot_build <- ggplot2::ggplot_build(p[[1]])

  # Check point size
  expect_equal(plot_build$data[[1]]$size[1], 2)

  # Check point alpha
  expect_equal(plot_build$data[[1]]$alpha[1], 0.5)

  # Check line color
  expect_equal(plot_build$data[[2]]$colour[1], "blue")
})

test_that("scientific notation formatting works", {
  x <- c(0, 10, 100, 10000)
  formatted <- scientific_format_end(x)

  expect_equal(formatted[1], "0")
  expect_equal(formatted[2:3], c("", ""))
  expect_match(formatted[4], "1e\\+4")
})

test_that("save plots", {

  temp_pdf_path <- file.path(tempdir(), "midar_test_responsecurve.pdf")


  p <- plot_feature_correlations(
    mexp,
    variable = "intensity",
    cor_min = 0.85,
    point_size = 2,
    point_alpha = 0.5,
    line_color = "blue",
    output_pdf = TRUE,
    path = temp_pdf_path,

    return_plot = FALSE,
    font_base_size = 10
  )

  expect_null(p)
  expect_true(file_exists(temp_pdf_path), info = "PDF file was not created.")
  expect_equal(as.character(fs::file_size(temp_pdf_path)), "239K")
  fs::file_delete(temp_pdf_path)

})
