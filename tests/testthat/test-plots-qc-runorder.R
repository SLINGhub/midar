# library(fs)
# library(vdiffr)
# library(ggplot2)
# library(testthat)
set.seed(123)
testthat::local_edition(3)

mexp_orig <- lipidomics_dataset

mexp <- normalize_by_istd(mexp_orig)
mexp <- quantify_by_istd(mexp)
mexp <- calc_qc_metrics(mexp)  # Ensure calc_qc_metrics is executed before

test_that("plot_runsequence works with basic parameters", {
  p <- plot_runsequence(mexp, show_batches = TRUE)
  expect_s3_class(p, "gg")
  plot_data <- ggplot2::ggplot_build(p)$data
  expect_equal(plot_data[[1]]$xintercept[1],93.5) # check no date as x axis
  expect_equal(dim(plot_data[[2]]),c(499,10))

  p <- plot_runsequence(data = mexp, show_timestamp = TRUE)
  expect_true(any(grepl("Acquisition Time", p$labels$x))) # Check x-axis label

  p <- plot_runsequence(data = mexp, single_row = TRUE)
  expect_true(any(grepl("y", p$mapping))) # Check if y axis is not present

  p <- plot_runsequence(data = mexp, show_batches = FALSE)
  plot_data <- ggplot2::ggplot_build(p)$data
  expect_equal(length(plot_data),1) # additional layer for batches geoms


  p <- plot_runsequence(data = mexp, show_batches = TRUE,
                        batch_zebra_stripe = FALSE)
  plot_data <- ggplot2::ggplot_build(p)$data
  expect_equal(length(plot_data),2) # additional layer for batches geoms
  expect_equal(dim(plot_data[[1]]),c(6,7)) # rows for each line

  p <- plot_runsequence(data = mexp, show_batches = TRUE, batch_zebra_stripe = TRUE)
  plot_data <- ggplot2::ggplot_build(p)$data
  expect_equal(length(plot_data),2) # additional layer for batches geoms
  expect_equal(dim(plot_data[[1]]),c(3,11)) # rows for each stripe

  p <- plot_runsequence(mexp, show_batches = TRUE, show_timestamp = TRUE)
  plot_data <- ggplot2::ggplot_build(p)$data
  ts <- as.POSIXct(plot_data[[1]]$xintercept[1],
                   origin = "1970-01-01",
                  tz = "Asia/Singapore") 
  expect_equal(format(ts, "%Y-%m-%d %H:%M:%S %Z"),
               "2017-10-20 22:15:36 +08") # check x axis uses data

  p <- plot_runsequence(mexp,qc_types = c("SPL", "BQC", "RQC"), show_batches = TRUE)
  plot_data <- ggplot2::ggplot_build(p)$data
  expect_equal(dim(plot_data[[2]]),c(444,10))

  p <- plot_runsequence(mexp,qc_types = c("SPL|BQC|RQC"), show_batches = TRUE)
  plot_data <- ggplot2::ggplot_build(p)$data
  # expect_equal(dim(plot_data[[2]]),c(445,10))

})


test_that("plot_rla_boxplot works with basic parameters", {
  p <- plot_rla_boxplot(mexp,
                        rla_type_batch = "within",
                        variable = "intensity",
                        show_batches = TRUE)
  expect_s3_class(p$plot, "gg")
  plot_data <- ggplot2::ggplot_build(p$plot)$data
  expect_equal(plot_data[[1]]$xintercept[1],93.5) # check batches shown
  expect_equal(dim(plot_data[[2]]),c(499,29)) # more columns for box plots
  expect_equal(mean(plot_data[[2]]$middle),-0.193117043184098)
  vdiffr::expect_doppelganger("default plot_rla_boxplot plot", p$plot)


p <- plot_rla_boxplot(mexp,
                      rla_type_batch = "across",
                      variable = "intensity",
                      show_batches = TRUE)

plot_data <- ggplot2::ggplot_build(p$plot)$data
expect_equal(mean(plot_data[[2]]$middle),-0.185713556)
expect_equal(max(plot_data[[2]]$x),499)

expect_no_error(
  p <- plot_rla_boxplot(mexp,
                        rla_type_batch = "across",
                        variable = "intensity",
                        show_batches = TRUE,
                        batch_zebra_stripe = FALSE)
)
expect_no_error(
  p <- plot_rla_boxplot(mexp,
                        rla_type_batch = "across",
                        variable = "intensity",
                        show_batches = TRUE,
                        batch_zebra_stripe = TRUE)
)

plot_data <- ggplot2::ggplot_build(p$plot)$data
expect_equal(dim(plot_data[[1]]),c(3,11))
expect_equal(max(plot_data[[2]]$x),499)

p <- plot_rla_boxplot(mexp,
                      rla_type_batch = "across",
                      variable = "intensity",
                      show_batches = TRUE,
                      batch_zebra_stripe = TRUE)
plot_data <- ggplot2::ggplot_build(p$plot)$data
expect_equal(dim(plot_data[[1]]),c(3,11)) # rows for each stripe

p <- plot_rla_boxplot(mexp,
                      rla_type_batch = "across",
                      variable = "intensity",
                      show_batches = TRUE,
                      show_timestamp = TRUE)
plot_data <- ggplot2::ggplot_build(p$plot)$data
vdiffr::expect_doppelganger("dt plot_rla_boxplot", p$plot)

p <- plot_rla_boxplot(mexp,
                      rla_type_batch = "across",
                      variable = "intensity",
                      show_batches = TRUE,
                      show_timestamp = FALSE,
                      x_gridlines = TRUE)
vdiffr::expect_doppelganger("x_gridlines_rla", p$plot)


p <- plot_rla_boxplot(mexp,
                      rla_type_batch = "within",
                      variable = "intensity",
                      show_batches = FALSE)


  # no batches shown
vdiffr::expect_doppelganger("nobatches_rla", p$plot)
  
 # Ensure plot range is applied correctly, not changing calculations 
p <- plot_rla_boxplot(mexp,
                      rla_type_batch = "across",
                      variable = "intensity",
                      plot_range = c(100,400),
                      show_batches = TRUE)

plot_data <- ggplot2::ggplot_build(p$plot)$data
expect_equal(mean(plot_data[[2]]$middle),-0.185713556)
expect_equal(min(plot_data[[2]]$x),1)
expect_equal(max(plot_data[[2]]$x),499)
vdiffr::expect_doppelganger("plotrange_rla", p$plot)

p <- plot_rla_boxplot(mexp,
                      rla_type_batch = "across",
                      variable = "intensity",
                      outlier_hide = TRUE,
                      show_batches = TRUE)

vdiffr::expect_doppelganger("outlierrem_rla", p$plot)

p <- plot_rla_boxplot(mexp,
                      rla_type_batch = "across",
                      variable = "intensity",
                      y_lim = c(-2,3),
                      show_batches = TRUE)

vdiffr::expect_doppelganger("ylimrla", p$plot)

})
