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


test_that("plot_runscatter generates a plot", {

  # Test with valid arguments
  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  # Check if a ggplot object is returned
  expect_s3_class(p[[1]], "gg")

  # Check how many pages
  expect_equal(length(p), 3)

  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[2]]
  expect_equal(nrow(plot_data), 5988)
  vdiffr::expect_doppelganger("def_runscatter", p)

  # log y axis
  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    show_batches = TRUE,
    batch_zebra_stripe = TRUE,
    log_scale = TRUE,
    return_plots = TRUE
  )
  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(max(plot_data[[2]]$y),7.09307642)
  vdiffr::expect_doppelganger("log_runscat", p)

  mexp_withzero <- mexp
  mexp_withzero@dataset$feature_intensity[sample(1:499, 10)] <- 0
  expect_message(
    p <- plot_runscatter(
      data = mexp_withzero,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      show_batches = TRUE,
      batch_zebra_stripe = TRUE,
      log_scale = TRUE,
      return_plots = TRUE
    ),
    "Zero or negative values were replaced")


  # Test with valid arguments
  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    specific_page = 3,
    return_plots = TRUE
  )
  # Check how many pages
  expect_equal(length(p), 1)

  # Test with valid arguments
  p <- plot_runscatter(
    data = mexp,
    variable = "conc",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_false(is.null(p))
  expect_equal(max(p[[1]]$data$value_mod),18.2841845)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[2]]
  expect_equal(nrow(plot_data), 5988)

  temp_pdf_path <- file.path(tempdir(), "midar_test_runscay.pdf")
  p <- plot_runscatter(
    data = mexp,
    variable = "conc",
    rows_page = 3,
    cols_page = 4,
    return_plots = FALSE,
    output_pdf = TRUE,
    path = temp_pdf_path
  )
  expect_null(p)
  expect_true(file_exists(temp_pdf_path), info = "PDF file was not created.")
  expect_equal(as.character(fs::file_size(temp_pdf_path)), "118K")
  fs::file_delete(temp_pdf_path)

})

# check diverse feature filters
test_that("plot_runscatter filter work", {
  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    qc_types = c("BQC", "SPL"),
    include_qualifier = FALSE,
    include_feature_filter = "PC|PE",
    exclude_feature_filter = "ISTD",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect(length(p), 1)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(nrow(plot_data[[2]]), 3456)
  vdiffr::expect_doppelganger("filteredrunscat", p)

  # Test range filter and also regex for qc type
  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    qc_types = c("BQC|SPL"),
    include_qualifier = FALSE,
    include_feature_filter = "PC|PE",
    exclude_feature_filter = "ISTD|SIM",
    plot_range = c(100,400),
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect(length(p), 1)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(nrow(plot_data[[2]]), 3456)
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  vdiffr::expect_doppelganger("rangerunscat", p)

  # Test include e
  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    qc_types = c("BQC|SPL"),
    include_qualifier = FALSE,
    include_feature_filter = c("PC 40:6", "PC 40:8"),
    plot_range = c(100,400),
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect(length(p), 1)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(nrow(plot_data[[2]]), 864)
  expect_equal(max(p[[1]]$data$value_mod),9364398.0)
  vdiffr::expect_doppelganger("filtrunscat2", p)

  expect_error(
    p <- plot_runscatter(
      data = mexp,
      variable = "intensity",
      qc_types = c("BQC|SPL"),
      include_qualifier = FALSE,
      include_feature_filter = "40",
      exclude_feature_filter = "PC",
      plot_range = c(100,400),
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE
    ),
    "filter criteria resulted in no selected features "
  )

 # TODO: Remove after confirm that purrr progress bar works

  # captured_output <- capture_output(
  #   plot_runscatter(
  #     data = mexp,
  #     variable = "conc",
  #     rows_page = 3,
  #     cols_page = 4,
  #     return_plots = FALSE,
  #     show_progress = TRUE
  # ))

  # expect_false(any(grepl("======================",
  #                       captured_output)))

  # captured_output <- capture_output(
  #   plot_runscatter(
  #     data = mexp,
  #     variable = "conc",
  #     rows_page = 3,
  #     cols_page = 4,
  #     return_plots = FALSE,
  #     show_progress = FALSE
  #   ))

  # expect_false(any(grepl("==",
  #                       captured_output)))

})

# check diverse feature filters
test_that("plot_runscatter with unknown qc_types", {

  mexp_newqc <- import_data_csv(data = MidarExperiment(),
                          path = test_path("testdata/plain_wide_dataset2_22rows_unknownQC.csv"),
                          variable_name = "conc",
                          analysis_id_col = "analysis_id",
                          import_metadata = TRUE)


 expect_message(
   p <- plot_runscatter(
    data = mexp_newqc,
    variable = "conc",
    include_qualifier = FALSE,
    rows_page = 3,
    cols_page = 4,
    point_size = 10,
    return_plots = TRUE
  ),
  "QC types 'XYX, MyQC, BLK' not predefined in Midar and will be displayed in black with auto-assigned shapes", fixed = TRUE)


  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[2]]
  expect_equal(length(unique(plot_data$shape)), 6)
  vdiffr::expect_doppelganger("runscatunknownqc", p)
})





test_that("plot_runscatter show batches works", {

  # Test with valid arguments
  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    show_batches = FALSE,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  # No batches lines/shapes, so data is in first item
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[1]]
  expect_equal(nrow(plot_data), 5988)

  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    show_batches = TRUE,
    batch_zebra_stripe = TRUE,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  # No batches lines/shapes, so data is in first item
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[1]]
  expect_equal(nrow(plot_data), 36)

})


test_that("plot_runscatter outlier cap works", {

  # Test with valid arguments
  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE,
    cap_outliers = TRUE,
    cap_sample_k_mad = 2,
    cap_qc_k_mad = 2
  )
  vdiffr::expect_doppelganger("runscattercapoutlier", p)

  expect_equal(max(p[[1]]$data$value_mod),6481466.4)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[3]] # not fully understand this test
  expect_equal(nrow(plot_data), 5988)

  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE,
    cap_outliers = TRUE,
    cap_sample_k_mad =3,
    cap_qc_k_mad = 3
  )
  expect_equal(max(p[[1]]$data$value_mod),7694153.4)


  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE,
    cap_outliers = TRUE,
    cap_sample_k_mad = NA,
    cap_qc_k_mad = NA,
    cap_top_n_outliers = 0,
  )
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)

  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE,
    cap_outliers = TRUE,
    cap_sample_k_mad = NA,
    cap_qc_k_mad = NA,
    cap_top_n_outliers = 1,
  )
  expect_equal(max(p[[1]]$data$value_mod),10147223.2)

  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE,
    cap_outliers = TRUE,
    cap_sample_k_mad = NA,
    cap_qc_k_mad = NA,
    cap_top_n_outliers = 30,
  )
  expect_equal(max(p[[1]]$data$value_mod),6068838.3)

  expect_error(
    p <- plot_runscatter(
      data = mexp,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE,
      cap_outliers = TRUE,
      cap_sample_k_mad = NA,
      cap_qc_k_mad = NA,
      cap_top_n_outliers = NA,
    ),
    "One or more of `cap_sample_k_mad`"
  )

  })


test_that("plot_runscatter show reference lines works", {

  # No batch-wise reference lines
  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    show_batches = TRUE,
    show_reference_lines = TRUE,
    ref_qc_types = "BQC",
    reference_batchwise = FALSE,
    reference_sd_shade = FALSE,
    reference_k_sd = 2,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[4]]
  expect_equal(unique(plot_data$linetype),"dashed")



  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    show_batches = TRUE,
    show_reference_lines = TRUE,
    ref_qc_types = "BQC",
    reference_batchwise = FALSE,
    reference_sd_shade = TRUE,
    reference_k_sd = 2,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[2]]
  expect_equal(unique(plot_data$alpha),0.15)

  # batch-wise reference lines

  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    show_batches = TRUE,
    show_reference_lines = TRUE,
    ref_qc_types = "BQC",
    reference_batchwise = TRUE,
    reference_sd_shade = FALSE,
    reference_k_sd = 2,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(length(plot_data),6)
  expect_equal(mean(plot_data[[3]]$yend),2103609.1)


  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    show_batches = TRUE,
    show_reference_lines = TRUE,
    ref_qc_types = "BQC",
    reference_batchwise = TRUE,
    reference_sd_shade = TRUE,
    reference_k_sd = 2,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(length(plot_data),5)
  expect_equal(unique(plot_data[[2]]$alpha),0.15)
  expect_equal(mean(plot_data[[2]]$ymax),2408485.4)
  vdiffr::expect_doppelganger("extrunscatterref", p)

  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    show_batches = TRUE,
    show_reference_lines = TRUE,
    reference_batchwise = FALSE,
    ref_qc_types = "BQC",
    reference_k_sd = NA,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(length(plot_data),4)
  expect_equal(mean(plot_data[[3]]$yintercept), 2103633.9)

  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    show_batches = TRUE,
    show_reference_lines = TRUE,
    reference_batchwise = TRUE,
    ref_qc_types = "SPL",
    reference_k_sd = 2,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_equal(max(p[[1]]$data$value_mod),12390146.0)
  plot_data <- ggplot2::ggplot_build(p[[1]])$data[[5]]
  expect_equal(dim(plot_data),c(72,10))


  expect_error(
    p <- plot_runscatter(
      data = mexp,
      variable = "intensity",
      show_reference_lines = TRUE,
      rows_page = 3,
      cols_page = 4,
      show_trend = TRUE,
      return_plots = TRUE
    ),
    "Please define a QC to show reference lines"
  )

  expect_error(
    p <- plot_runscatter(
      data = mexp,
      variable = "intensity",
      show_reference_lines = TRUE,
      ref_qc_types = "CAL",
      rows_page = 3,
      cols_page = 4,
      show_trend = TRUE,
      return_plots = TRUE
    ),
    "Selected `ref_qc_types` not"
  )
})


test_that("plot_runscatter show trend works", {

  expect_error(
    p <- plot_runscatter(
      data = mexp,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      show_trend = TRUE,
      return_plots = TRUE
    ),
    "Drift or batch correction is currently required to show"
  )

  mexp_drift <- correct_drift_gaussiankernel(
    mexp_orig,
    variable = "intensity",
    recalc_trend_after = TRUE,

    ref_qc_types = "SPL",
    ignore_istd = FALSE)

  p <- plot_runscatter(
    data = mexp_drift,
    variable = "intensity",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = TRUE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184,10))
  expect_equal(mean(plot_data[[2]]$y),2054300.962) #  data points
  expect_equal(dim(plot_data[[3]]),c(5184,9)) # smoothed ref data points
  expect_equal(mean(plot_data[[3]]$y),1985330.456) # smoothed ref data points

  p <- plot_runscatter(
    data = mexp_drift,
    variable = "intensity_raw",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = TRUE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184, 10))
  expect_equal(dim(plot_data[[3]]),c(5184, 9)) # smoothed data points
  expect_equal(mean(plot_data[[2]]$y),2044597.1) # data points
  expect_equal(mean(plot_data[[3]]$y),1982570.1) # ref data points (batch-wise)

  p <- plot_runscatter(
    data = mexp_drift,
    variable = "intensity_before",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = TRUE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184, 10))
  expect_equal(dim(plot_data[[3]]),c(5184, 9)) # smoothed data points
  expect_equal(mean(plot_data[[2]]$y),2044597.1) # data points
  expect_equal(mean(plot_data[[3]]$y),1982570.1) # ref data points (batch-wise)

  mexp_drift <- correct_batch_centering(
    mexp_orig,
    ref_qc_types = "BQC",
    variable = "intensity")

  p <- plot_runscatter(
    data = mexp_drift,
    variable = "intensity",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = TRUE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184, 10))
  expect_equal(dim(plot_data[[3]]),c(5184, 9)) # ref data points
  expect_equal(mean(plot_data[[2]]$y),2039591.3) #  data points
  expect_equal(mean(plot_data[[3]]$y),2084560.8) #  ref data points (batch-wise)

  p <- plot_runscatter(
    data = mexp_drift,
    variable = "intensity_raw",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = TRUE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184, 10))
  expect_equal(dim(plot_data[[3]]),c(5184, 9)) # smoothed data points
  expect_equal(mean(plot_data[[2]]$y),2044597.1) # data points
  expect_equal(mean(plot_data[[3]]$y),2093909.456) # ref data points (batch-wise)

  p <- plot_runscatter(
    data = mexp_drift,
    variable = "intensity_raw",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = FALSE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184, 10))
  expect_equal(mean(plot_data[[2]]$y),2044597.1) # data points

  p <- plot_runscatter(
    data = mexp_drift,
    variable = "intensity_before",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = TRUE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184, 10))
  expect_equal(dim(plot_data[[3]]),c(5184, 9)) # smoothed data points
  expect_equal(mean(plot_data[[2]]$y),2044597.1) # data points
  expect_equal(mean(plot_data[[3]]$y),2093909.456) # ref data points (batch-wise)

  # compare again with orginal data from mexp orginal
  p <- plot_runscatter(
    data = mexp_orig,
    variable = "intensity",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = FALSE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184,10))
  expect_equal(mean(plot_data[[2]]$y),2044597.1) # data points

  mexp_drift <- correct_drift_gaussiankernel(
    mexp_orig,
    variable = "intensity",
    recalc_trend_after = TRUE,

    ref_qc_types = "SPL",
    ignore_istd = FALSE)

  mexp_drift2 <- correct_batch_centering(
    mexp_drift,
    ref_qc_types = "SPL",
    variable = "intensity")


  p <- plot_runscatter(
    data = mexp_drift2,
    variable = "intensity_raw",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = TRUE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184, 10))
  expect_equal(dim(plot_data[[3]]),c(5184, 9)) # smoothed data points
  expect_equal(mean(plot_data[[2]]$y),2044597.1) # data points
  expect_equal(mean(plot_data[[3]]$y),1982570.1) # ref data points (batch-wise)

  p <- plot_runscatter(
    data = mexp_drift2,
    variable = "intensity_before",
    qc_types = c("BQC", "SPL"),
    rows_page = 3,
    cols_page = 4,
    show_trend = TRUE,
    return_plots = TRUE
  )

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(5184, 10))
  expect_equal(dim(plot_data[[3]]),c(5184, 9)) # smoothed data points
  expect_equal(mean(plot_data[[2]]$y),2054300.962) # data points
  expect_equal(mean(plot_data[[3]]$y),1982570.1) # ref data points (batch-wise)


  expect_error(
    p <- plot_runscatter(
      data = mexp_orig,
      variable = "intensity_before",
      qc_types = c("BQC", "SPL"),
      rows_page = 3,
      cols_page = 4,
      show_trend = FALSE,
      return_plots = TRUE
    ),
  "Variables `_before` and `_raw` after only available after drift/batch corrections",
  fixed = TRUE)

  expect_error(
    p <- plot_runscatter(
      data = mexp_orig,
      variable = "intensity",
      qc_types = c("BQC", "SPL"),
      rows_page = 3,
      cols_page = 4,
      show_trend = TRUE,
      return_plots = TRUE
    ),
    "Drift or batch correction is currently required to show trend lines",
    fixed = TRUE)

})

test_that("plot_runscatter handles filter and missing data", {

expect_error(
  p <- plot_runscatter(
    data = MidarExperiment(),
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE,
    cap_outliers = TRUE,
    cap_sample_k_mad = NA,
    cap_qc_k_mad = NA,
    cap_top_n_outliers = NA,
  ),
  "No data available"
)



  expect_error(
    p <- plot_runscatter(
      data = mexp,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE,
      y_min = "NO"
    ),
    "`y_min` must be a numeric value or", fixed = TRUE
  )

  expect_error(
    p <- plot_runscatter(
      data = mexp,
      filter_data = TRUE,
      variable = "intensity",
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE,
    ),
    "Data has not been QC-filtered"
  )
  mexp_filt <- midar::filter_features_qc(mexp, include_qualifier = FALSE, include_istd = FALSE,

                                         min.intensity.median.spl = 1000000)

  # uses filtered dataset
  p <- plot_runscatter(
    data = mexp_filt,
    variable = "intensity",
    filter_data = TRUE,
    rows_page = 3,
    cols_page = 4,
    return_plots = TRUE
  )
  expect_equal(length(p),2)

  # returns error if variable not present in dataset
  expect_error(
    p <- plot_runscatter(
      data = mexp_orig,
      variable = "conc",
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE,
    ),
    "Concentration data are not available"
  )

  # returns error if variable not present in dataset
  expect_error(
    p <- plot_runscatter(
      data = mexp,
      variable = "conc_raw",
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE,
    ),
    "Concentration data are not available"
  )

  # returns error if variable not present in dataset
  expect_error(
    p <- plot_runscatter(
      data = mexp,
      variable = "conc_raw",
      output_pdf = TRUE,
      rows_page = 3,
      cols_page = 4,
      return_plots = TRUE,
    ),
    "No valid path defined"
  )


  p <- plot_runscatter(
    data = mexp,
    variable = "intensity",
    rows_page = 3,
    cols_page = 4,
    plot_range = c(100,400),
    return_plots = TRUE,
  )
 vdiffr::expect_doppelganger("runscatterrange_filter",p)
})