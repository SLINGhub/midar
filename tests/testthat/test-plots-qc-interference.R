library(ggplot2)

mexp_orig <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))
mexp <- mexp_orig

mexp2 <- lipidomics_dataset

mexp2@annot_features$interference_contribution[9] <- 0.5


mexp_corrected <-
  correct_interferences(
    mexp,
    variable = "feature_intensity",
    sequential_correction = FALSE
  )


target_feature <- "S1P d18\\:0"


# Baseline test
test_that("Basic plot_qc_interferences looks as expected", {
  p <- plot_qc_interferences(mexp_corrected)
  # The name "qc-interferences-default" is more descriptive
  vdiffr::expect_doppelganger("qc-interferences-default", p)
})


# --- Tests for Data Filtering Parameters ---

test_that("Filtering by qc_types works", {
  # The default includes "SPL". Let's plot only two.
  p <- plot_qc_interferences(mexp_corrected, qc_types = c("SPL"))
  vdiffr::expect_doppelganger("qc-interferences-filter-qcs", p)
})

test_that("Excluding ISTDs works", {
  # The default is include_istd = TRUE. This should remove the ISTD features.
  p <- plot_qc_interferences(mexp_corrected, include_istd = FALSE)
  vdiffr::expect_doppelganger("qc-interferences-no-istd", p)
})

test_that("Including a specific feature works", {
  # This should result in a plot with only one feature on the x-axis.
  # We pick a known corrected ISTD from the data.
  p <- plot_qc_interferences(
    mexp_corrected,
    include_feature_filter = target_feature
  )
  vdiffr::expect_doppelganger("qc-interferences-include-filter", p)
})

test_that("Excluding a specific feature works", {
  # The plot should be missing the "LPC(17:0) (IS)" feature compared to the default.
  p <- plot_qc_interferences(
    mexp_corrected,
    exclude_feature_filter = target_feature
  )
  vdiffr::expect_doppelganger("qc-interferences-exclude-filter", p)
})

test_that("min_median_value filter works visually", {
  # From inspecting the data, a value of 9000 will filter out some but not all features.
  # This tests that the filtering logic is applied correctly.
  p <- plot_qc_interferences(mexp_corrected, min_median_value = 49000)
  vdiffr::expect_doppelganger("qc-interferences-min-median", p)
})

test_that("min_median_value throws error when no features remain", {
  # A very high value should filter out all features and trigger the error.
  expect_error(
    plot_qc_interferences(mexp_corrected, min_median_value = 99999999),
    "No features passed the `min_median_value` filter"
  )
})


# --- Tests for Aesthetic Parameters ---

test_that("Aesthetic parameters are applied correctly", {
  # Change several visual parameters at once to create a distinct plot.
  p <- plot_qc_interferences(
    mexp_corrected,
    y_lim = c(80, 120), # Zoom in on the y-axis
    point_size = 2, # Larger points
    point_alpha = 0.8, # More opaque points
    angle_x = 0, # Horizontal x-axis labels
    font_base_size = 12 # Larger font
  )
  vdiffr::expect_doppelganger("qc-interferences-aesthetics", p)
})

test_that("Plot works with NA qc_types to auto-detect", {
  # This tests the initial `if (all(is.na(qc_types)))` block.
  # The result should be identical to the default plot in this case.
  p <- plot_qc_interferences(mexp_corrected, qc_types = NA)
  vdiffr::expect_doppelganger("qc-interferences-na-qcs", p)
})


test_that("Object check: include_feature_filter correctly filters the data layer", {
  # Define the feature we want to isolate

  # Generate the plot
  p <- plot_qc_interferences(
    mexp_corrected,
    include_feature_filter = target_feature
  )

  # 1. Build the plot to access the data used for rendering
  built_plot <- ggplot_build(p)
  plot_data <- built_plot$data[[1]]
  x_axis_labels <- built_plot$layout$panel_params[[1]]$x$get_labels()

  expect_length(x_axis_labels, 1)
  expect_equal(x_axis_labels, "S1P d18:0 [M>60]")
})
