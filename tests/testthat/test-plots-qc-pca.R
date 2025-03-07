library(fs)
library(vdiffr)
library(ggplot2)
library(testthat)



mexp_orig <- lipidomics_dataset

mexp <- exclude_analyses(mexp_orig, c("Longit_batch6_51","Longit_batch6_B-ISTD 09"), clear_existing = TRUE)
mexp <- normalize_by_istd(mexp)
mexp <- quantify_by_istd(mexp)
mexp <- calc_qc_metrics(mexp)  # Ensure calc_qc_metrics is executed before

test_that("plot_pca works", {
  set.seed(123)
  p <- expect_message(class(plot_pca(mexp,
           variable = "intensity",
           ellipse_variable = "qc_type",
           ellipse_alpha = 0.2,
           filter_data = FALSE)),
  "The PCA was calculated based on \\`feature_intensity\\` values of 19 features")

  expect_doppelganger("default plot_pca plot 1", p)


    p <- plot_pca(mexp,
                  variable = "intensity",
                  ellipse_variable = "none",
                  labels_threshold_mad = 1,
                  ellipse_alpha = 0.2,
                  filter_data = FALSE)

    expect_warning(print(p), "unlabeled data points")


  p <- plot_pca(mexp,
                variable = "intensity",
                ellipse_variable = "batch_id",
                ellipse_alpha = 0.2,
                filter_data = FALSE)

  expect_doppelganger("default plot_pca plot 2", p)

  p <- plot_pca(mexp,
                variable = "intensity",
                ellipse_variable = "batch_id",
                ellipse_alpha = 0.02,
                show_labels = FALSE,
                filter_data = FALSE)

  expect_doppelganger("default plot_pca plot 3 no labels", p)

  p <- plot_pca(mexp,
                variable = "intensity",
                ellipse_variable = "batch_id",
                ellipse_alpha = 0.05,
                log_transform = TRUE,
                show_labels = TRUE,
                shared_labeltext_hide = "Longit_",
                filter_data = FALSE)

  expect_doppelganger("default plot_pca plot 4 label-shared-rm", p)
  plot_data <- ggplot_build(p)$data
  expect_equal(max(plot_data[[3]]$x, na.rm = T),6.838079)

  p <- plot_pca(mexp,
                variable = "intensity",
                ellipse_variable = "batch_id",
                ellipse_alpha = 0.05,
                log_transform = FALSE,
                filter_data = FALSE)

  plot_data <- ggplot_build(p)$data
  expect_equal(max(plot_data[[3]]$x, na.rm = T),3.8357779)

  plot_pca(mexp,
           variable = "intensity",
           ellipse_variable = "qc_type",
           ellipse_alpha = 0.2,
           filter_data = FALSE,
           ellipse_fill = TRUE,
           ellipse_levels = c("BQC", "TQC"))
  plot_pca(mexp,
           variable = "intensity",
           ellipse_variable = "qc_type",
           ellipse_alpha = 0.2,
           filter_data = FALSE,
           ellipse_fill = TRUE,
           ellipse_levels = c("BQC", "TQC"),
           ellipse_fillcolor = c("cyan", "blue"))

  expect_doppelganger("default plot_pca plot 5 user fill colors", p)

  p <- plot_pca(mexp,
           variable = "intensity",
           filter_data = FALSE,
           label_font_size = 5,
           ellipse_variable = "qc_type",
           ellipse_alpha = 0.2,

           ellipse_fill = TRUE,
           ellipse_levels = c("BQC", "TQC"),
           ellipse_fillcolor = c("BQC" = "cyan", "TQC" = "blue"))

  expect_doppelganger("default plot_pca plot 5 user mapped fill colors", p)

})

test_that("plot_pca filter work", {

  expect_message(
    p <- plot_pca(mexp,
                  variable = "intensity",
                  ellipse_variable = "batch_id",
                  ellipse_alpha = 0.05,
                  min_median_value =1000000,
                  log_transform = FALSE,
                  filter_data = FALSE),
    "values of 13 features.")

  plot_data <- ggplot_build(p)$data
  expect_equal(max(plot_data[[3]]$x, na.rm = T),3.50274432)

  # check missing values
  mexp_temp <- mexp
  mexp_temp@dataset$feature_intensity[c(411,2212,5133)] <- NA
  expect_message(
    p <- plot_pca(mexp_temp,
                  variable = "intensity",
                  ellipse_variable = "batch_id",
                  ellipse_alpha = 0.05,
                  min_median_value =1000000,
                  log_transform = FALSE,
                  filter_data = FALSE),
    "2 features contained missing or non-numeric values and were exluded")

  expect_message(
    p <- plot_pca(mexp_temp,
                  variable = "intensity",
                  ellipse_variable = "batch_id",
                  ellipse_alpha = 0.05,
                  min_median_value =1000000,
                  log_transform = FALSE,
                  filter_data = FALSE),
    "values of 13 features")

})

test_that("plot_pca error handler work", {
  set.seed(123)
  expect_error(plot_pca(mexp,
             variable = "intensity",
             ellipse_variable = "batch_id",
             ellipse_alpha = 0.2,
             filter_data = FALSE,
             ellipse_fill = TRUE,
             ellipse_levels = c("BQC", "TQC", "SPL")),
             "One or more levels in ")

  expect_error(  plot_pca(mexp,
                          variable = "intensity",
                          ellipse_variable = "qc_type",
                          ellipse_alpha = 0.2,
                          filter_data = FALSE,
                          ellipse_fill = TRUE,
                          ellipse_levels = c("BQC", "TQC", "SPL"),
                          ellipse_fillcolor = c("cyan", "blue")),
                          "Insufficient colors "
                          )
  expect_error(
    p <- plot_pca(mexp,
                 variable = "intensity",
                 ellipse_variable = "batch_id",
                 ellipse_alpha = 0.05,
                 min_median_value =32645924,
                 log_transform = FALSE,
                 filter_data = FALSE),
    "Only 1 feature passed the")

  expect_error(
    p <- plot_pca(mexp,
                  variable = "intensity",
                  ellipse_variable = "batch_id",
                  ellipse_alpha = 0.05,
                  min_median_value =1132645924,
                  log_transform = FALSE,
                  filter_data = FALSE),
    "No features passed the")

  expect_error(
    p <- plot_pca(mexp,
                  variable = "intensity",
                  ellipse_variable = "batch_id",
                  ellipse_alpha = 0.05,
                  log_transform = TRUE,
                  show_labels = TRUE,
                  shared_labeltext_hide = "batch[0-9]",
                  filter_data = FALSE),
  "duplicate labels")

})
