library(fs)
library(vdiffr)
library(ggplot2)
# library(testthat)

mexp_orig <- lipidomics_dataset

mexp <- exclude_analyses(
  mexp_orig,
  c("Longit_batch6_51", "Longit_batch6_B-ISTD 09"),
  clear_existing = TRUE
)
mexp <- normalize_by_istd(mexp)
mexp <- quantify_by_istd(mexp)
mexp <- calc_qc_metrics(mexp) # Ensure calc_qc_metrics is executed before

test_that("plot_pca works", {
  set.seed(123)
  p <- expect_message(
    class(plot_pca(
      mexp,
      variable = "intensity",
      ellipse_variable = "qc_type",
      ellipse_alpha = 0.2,
      filter_data = FALSE
    )),
    "The PCA was calculated based on \\`feature_intensity\\` values of 19 features"
  )

  vdiffr::expect_doppelganger("default plot_pca plot 1", p)

  p <- plot_pca(
    mexp,
    variable = "intensity",
    ellipse_variable = "none",
    labels_threshold_mad = 1,
    ellipse_alpha = 0.2,
    filter_data = FALSE
  )

  expect_warning(print(p), "unlabeled data points")

  p <- plot_pca(
    mexp,
    variable = "intensity",
    ellipse_variable = "batch_id",
    ellipse_alpha = 0.2,
    filter_data = FALSE
  )

  vdiffr::expect_doppelganger("default plot_pca plot 2", p)

  p <- plot_pca(
    mexp,
    variable = "intensity",
    ellipse_variable = "batch_id",
    ellipse_alpha = 0.02,
    show_labels = FALSE,
    filter_data = FALSE
  )

  vdiffr::expect_doppelganger("default plot_pca plot 3 no labels", p)

  p <- plot_pca(
    mexp,
    variable = "intensity",
    ellipse_variable = "batch_id",
    ellipse_alpha = 0.05,
    log_transform = TRUE,
    show_labels = TRUE,
    shared_labeltext_hide = "Longit_",
    filter_data = FALSE
  )

  vdiffr::expect_doppelganger("default plot_pca plot 4 label-shared-rm", p)
  plot_data <- ggplot2::ggplot_build(p)$data
  expect_equal(max(plot_data[[3]]$x, na.rm = T), 6.838079)

  p <- plot_pca(
    mexp,
    variable = "intensity",
    ellipse_variable = "batch_id",
    ellipse_alpha = 0.05,
    log_transform = FALSE,
    filter_data = FALSE
  )

  plot_data <- ggplot2::ggplot_build(p)$data
  expect_equal(max(plot_data[[3]]$x, na.rm = T), 3.8357779)

  plot_pca(
    mexp,
    variable = "intensity",
    ellipse_variable = "qc_type",
    ellipse_alpha = 0.2,
    filter_data = FALSE,
    ellipse_fill = TRUE,
    ellipse_levels = c("BQC", "TQC")
  )
  plot_pca(
    mexp,
    variable = "intensity",
    ellipse_variable = "qc_type",
    ellipse_alpha = 0.2,
    filter_data = FALSE,
    ellipse_fill = TRUE,
    ellipse_levels = c("BQC", "TQC"),
    ellipse_fillcolor = c("cyan", "blue")
  )

  vdiffr::expect_doppelganger("default plot_pca plot 5 user fill colors", p)

  p <- plot_pca(
    mexp,
    variable = "intensity",
    filter_data = FALSE,
    label_font_size = 5,
    ellipse_variable = "qc_type",
    ellipse_alpha = 0.2,

    ellipse_fill = TRUE,
    ellipse_levels = c("BQC", "TQC"),
    ellipse_fillcolor = c("BQC" = "cyan", "TQC" = "blue")
  )

  vdiffr::expect_doppelganger(
    "default plot_pca plot 5 user mapped fill colors",
    p
  )
})

test_that("plot_pca filter work", {
  expect_message(
    p <- plot_pca(
      mexp,
      variable = "intensity",
      ellipse_variable = "batch_id",
      ellipse_alpha = 0.05,
      min_median_value = 1000000,
      log_transform = FALSE,
      filter_data = FALSE
    ),
    "values of 13 features."
  )

  plot_data <- ggplot2::ggplot_build(p)$data
  expect_equal(max(plot_data[[3]]$x, na.rm = T), 3.50274432)

  # check missing values
  mexp_temp <- mexp
  mexp_temp@dataset$feature_intensity[c(411, 2212, 5133)] <- NA
  expect_message(
    p <- plot_pca(
      mexp_temp,
      variable = "intensity",
      ellipse_variable = "batch_id",
      ellipse_alpha = 0.05,
      min_median_value = 1000000,
      log_transform = FALSE,
      filter_data = FALSE
    ),
    "2 features contained missing or non-numeric values and were exluded"
  )

  expect_message(
    p <- plot_pca(
      mexp_temp,
      variable = "intensity",
      ellipse_variable = "batch_id",
      ellipse_alpha = 0.05,
      min_median_value = 1000000,
      log_transform = FALSE,
      filter_data = FALSE
    ),
    "values of 13 features"
  )
})

test_that("plot_pca error handler work", {
  set.seed(123)
  expect_error(
    plot_pca(
      mexp,
      variable = "intensity",
      ellipse_variable = "batch_id",
      ellipse_alpha = 0.2,
      filter_data = FALSE,
      ellipse_fill = TRUE,
      ellipse_levels = c("BQC", "TQC", "SPL")
    ),
    "One or more levels in "
  )

  expect_error(
    plot_pca(
      mexp,
      variable = "intensity",
      ellipse_variable = "qc_type",
      ellipse_alpha = 0.2,
      filter_data = FALSE,
      ellipse_fill = TRUE,
      ellipse_levels = c("BQC", "TQC", "SPL"),
      ellipse_fillcolor = c("cyan", "blue")
    ),
    "Insufficient colors "
  )
  expect_error(
    p <- plot_pca(
      mexp,
      variable = "intensity",
      ellipse_variable = "batch_id",
      ellipse_alpha = 0.05,
      min_median_value = 32645924,
      log_transform = FALSE,
      filter_data = FALSE
    ),
    "Only 1 feature passed the"
  )

  expect_error(
    p <- plot_pca(
      mexp,
      variable = "intensity",
      ellipse_variable = "batch_id",
      ellipse_alpha = 0.05,
      min_median_value = 1132645924,
      log_transform = FALSE,
      filter_data = FALSE
    ),
    "No features passed the"
  )

  expect_error(
    p <- plot_pca(
      mexp,
      variable = "intensity",
      ellipse_variable = "batch_id",
      ellipse_alpha = 0.05,
      log_transform = TRUE,
      show_labels = TRUE,
      shared_labeltext_hide = "batch[0-9]",
      filter_data = FALSE
    ),
    "duplicate labels"
  )
})


# --- Baseline visual snapshot ---

test_that("Default plot_pca_loading looks as expected", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    # keep defaults: log_transform = TRUE, abs_loading = TRUE, vertical_bars = FALSE
    pca_dims = c(1, 2, 3, 4),
    top_n = 15 # smaller to keep snapshot readable and stable
  )
  vdiffr::expect_doppelganger("pca-loading-default", p)
})

# --- Core logic and parameterization snapshots ---

test_that("abs_loading = FALSE shows signed loadings", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = c(1, 2),
    abs_loading = FALSE,
    top_n = 20
  )
  vdiffr::expect_doppelganger("pca-loading-signed", p)
})

test_that("vertical_bars = TRUE changes orientation and layout", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = c(1, 2, 3),
    vertical_bars = TRUE,
    top_n = 10
  )
  vdiffr::expect_doppelganger("pca-loading-vertical-bars", p)
})

test_that("log_transform = FALSE produces a different loading plot", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = c(1, 2),
    log_transform = FALSE,
    top_n = 15
  )
  vdiffr::expect_doppelganger("pca-loading-no-log", p)
})

test_that("Custom pca_dims subset plots only requested PCs", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = c(1, 3),
    top_n = 10
  )
  vdiffr::expect_doppelganger("pca-loading-pcs-1-and-3", p)
})

test_that("Filtering qc_types changes the PCs input data", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    qc_types = c("SPL", "TQC"),
    pca_dims = c(1, 2),
    top_n = 10
  )
  vdiffr::expect_doppelganger("pca-loading-qctypes-spl-tqc", p)
})

test_that("Auto-detected qc_types works (qc_types = NA)", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    qc_types = NA,
    pca_dims = c(1, 2),
    top_n = 10
  )
  vdiffr::expect_doppelganger("pca-loading-qctypes-auto", p)
})

# --- Data filtering with min_median_value ---

test_that("min_median_value filters features and still produces a plot", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = c(1, 2),
    min_median_value = 5e5, # moderate threshold: some features filtered out
    top_n = 10
  )
  vdiffr::expect_doppelganger("pca-loading-min-median", p)
})

test_that("min_median_value throws error when no features pass", {
  expect_error(
    plot_pca_loading(
      data = mexp,
      variable = "intensity",
      pca_dims = c(1, 2),
      min_median_value = 1e12, # extremely high to filter all
      top_n = 10
    ),
    "No features passed the `min_median_value` filter"
  )
})

test_that("min_median_value throws error when only one feature passes", {
  # Choose a threshold high enough to leave 1 feature; adjust if needed for your dataset.
  # If your dataset doesn't hit exactly 1, you can select a narrower QC set + threshold.
  expect_error(
    plot_pca_loading(
      data = mexp,
      variable = "intensity",
      qc_types = c("SPL", "TQC"),
      pca_dims = c(1, 2),
      min_median_value = 1.5e7,
      top_n = 10
    ),
    "Only 1 feature passed the `min_median_value` filter"
  )
})

# --- Object-level checks (unit tests on ggplot object) ---

test_that("Number of bars equals top_n * number of PCs", {
  top_n <- 7
  pcs <- c(1, 2, 3)
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = pcs,
    top_n = top_n
  )
  bp <- ggplot_build(p)
  # First layer is geom_col
  n_rows <- nrow(bp$data[[1]])
  expect_equal(n_rows, top_n * length(pcs))
})

test_that("Facet panel count equals number of requested PCs", {
  pcs <- c(1, 4)
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = pcs,
    top_n = 5
  )
  bp <- ggplot_build(p)
  # panel_params length equals number of panels
  expect_equal(length(bp$layout$panel_params), length(pcs))
})

test_that("vertical_bars affects coordinate system (no flip) and columns in facet wrap", {
  pcs <- c(1, 2, 3)
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = pcs,
    vertical_bars = TRUE,
    top_n = 6
  )
  # No coord_flip when vertical_bars = TRUE
  expect_false(inherits(p$coordinates, "CoordFlip"))
  bp <- ggplot_build(p)
  # facet_wrap uses ncol = 1 in this mode; so COL should be 1 for all panels
  lay <- bp$layout$layout
  expect_equal(length(unique(lay$COL)), 1L)
})

test_that("When abs_loading = FALSE, both positive and negative directions can appear", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = c(1, 2),
    abs_loading = FALSE,
    top_n = 20
  )
  bp <- ggplot_build(p)
  # Extract data used for the bars
  d1 <- bp$data[[1]]
  # direction is mapped to color/fill; in the data, negative bars should have negative 'x' or 'y' values
  # depending on coord_flip. If not flipped (vertical_bars = FALSE -> coord_flip), negatives will be in 'x'.
  if (inherits(p$coordinates, "CoordFlip")) {
    expect_true(any(d1$ymin < 0))
  } else {
    expect_true(any(d1$xmin < 0))
  }
})

test_that("When abs_loading = FALSE, bars are non-negative in plotted axis", {
  p <- plot_pca_loading(
    data = mexp,
    variable = "intensity",
    pca_dims = c(1, 2),
    abs_loading = TRUE,
    top_n = 20
  )
  bp <- ggplot_build(p)
  d1 <- bp$data[[1]]
  if (inherits(p$coordinates, "CoordFlip")) {
    expect_true(all(d1$x >= 0))
  } else {
    expect_true(all(d1$y >= 0))
  }
})
