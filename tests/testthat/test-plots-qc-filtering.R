library(fs)
library(vdiffr)
library(ggplot2)
library(testthat)
library(scales)

mexp_orig <- lipidomics_dataset
mexp <- normalize_by_istd(mexp_orig)
mexp <- quantify_by_istd(mexp)
mexp_proc <- calc_qc_metrics(mexp,  use_batch_medians = FALSE)

mexp_filt_all <- filter_features_qc(
  mexp_proc,
  clear_existing = TRUE,
  include_qualifier = FALSE,
  include_istd = FALSE,
  min.intensity.lowest.tqc = 0)

get_feature_n <- function(plt) {
  df <- ggplot_build(plt)$data |> as.data.frame() |> filter(group == 1)
  sum(df$y)
}


test_that("plot_qc_summary_byclass plots correctly", {
  mexp_res <- filter_features_qc(
    mexp_proc,
    clear_existing = TRUE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    min.intensity.lowest.tqc = 0)


  p <- plot_qc_summary_byclass(mexp_res )
  expect_equal( get_feature_n(p), 19)

  p <- plot_qc_summary_overall(mexp_res )
  expect_doppelganger("plot_qc_summary_summ with no fails ", p)

  mexp_res <- filter_features_qc(
    mexp_proc,
    clear_existing = TRUE,
    include_qualifier = TRUE,
    include_istd = TRUE,
    min.intensity.lowest.tqc = 0)

  p <- plot_qc_summary_byclass(mexp_res)
  expect_equal( get_feature_n(p), 29)
  p <- plot_qc_summary_overall(mexp_res )
  expect_doppelganger("plot_qc_summary_summ with no fails with qual ", p)

  mexp_res2 <- filter_features_qc(
    mexp_proc,
    include_qualifier = FALSE,
    include_istd = FALSE,
    clear_existing = TRUE,
    min.signalblank.median.spl.pblk = 100,
    min.intensity.median.bqc = 1E3,
    max.cv.conc.bqc = 23,
    max.dratio.sd.conc.bqc = 0.7,
    max.prop.missing.conc.spl = 0.1,
    response.curves.selection = c(1,2),
    response.curves.summary = "best",
    max.slope.response = 1.05)

  p <- plot_qc_summary_byclass(mexp_res2)
  expect_equal( get_feature_n(p), 19)
  expect_doppelganger("plot_qc_summary_byclass with fails 1", p)

  p <- plot_qc_summary_overall(mexp_res2 )
  expect_doppelganger("plot_qc_summary_summ with fails ", p)

  p <- plot_qc_summary_overall(mexp_res2,with_venn = FALSE )
  expect_doppelganger("plot_qc_summary_summ with fails no venn", p)

})



test_that("plot_qc_summary_x handle errors", {

  expect_error(
    plot_qc_summary_byclass(mexp_proc),
    "Feature QC filter has not yet been applied"
  )

  expect_error(
    plot_qc_summary_overall(mexp_proc),
    "Feature QC filter has not yet been applied"
  )



  # No feature class defined
  mexp_temp <- mexp_filt_all
  mexp_temp@metrics_qc$feature_class <- NA

  expect_error(
    plot_qc_summary_byclass(mexp_temp),
    "This plot requires the `feature_class` to be defined in the data"
  )

  expect_no_error(
    plot_qc_summary_overall(mexp_temp)
  )

})
