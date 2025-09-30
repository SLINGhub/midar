
mexp <- midar::MidarExperiment(title = "sPerfect")

data_path <- test_path("testdata/FullPanelFewSamples_MRMkit.csv")
mexp <- import_data_mrmkit(data = mexp, path = data_path, import_metadata = TRUE)


test_that("plot_rt_vs_chain works", {
  p <- plot_rt_vs_chain(mexp, qc_types = "SPL", x_var = "total_c")
   vdiffr::expect_doppelganger("default plot_rt_vs_chain", p)
})

# Regre
test_that("plot_rt_vs_chain works with x axis db" , {
  p <- plot_rt_vs_chain(mexp, qc_types = "SPL", x_var = "total_db")
   vdiffr::expect_doppelganger("plot_rt_vs_chain xaxis db", p)
})

test_that("plot_rt_vs_chain works with x axis db" , {
  p <- plot_rt_vs_chain(mexp, qc_types = "SPL", x_var = "ecn")
   vdiffr::expect_doppelganger("plot_rt_vs_chain xaxis ecn", p)
})

test_that("plot_rt_vs_chain no robust regress", {
  p <- plot_rt_vs_chain(
    mexp, 
    qc_types = "SPL", 
    x_var = "total_c",
    robust_regression = FALSE)
   vdiffr::expect_doppelganger(" plot_rtr_vs_chain norobustreg", p)
})

test_that("plot_rt_vs_chain include qualifier", {

  mexp_temp <- mexp
  mexp_temp@dataset[str_detect(mexp_temp@dataset$feature_id, "PC 3"), ]$is_quantifier  <- FALSE
  p <- plot_rt_vs_chain(
    mexp_temp, 
    qc_types = "SPL", 
    x_var = "total_c",
  include_qualifier = TRUE)
   vdiffr::expect_doppelganger("plot_rt_vs_chain with qual", p)
  
    p <- plot_rt_vs_chain(
    mexp_temp, 
    qc_types = "SPL", 
    x_var = "total_c",
  include_qualifier = FALSE)
   vdiffr::expect_doppelganger("plot_rt_vs_chain no qual", p)
})

test_that("plot_rt_vs_chain no robust regress", {
  p <- plot_rt_vs_chain(
    mexp, 
    qc_types = "SPL", 
    x_var = "total_c",
    outliers_highlight = FALSE)
   vdiffr::expect_doppelganger(" plot_rtr_vs_chain hide outlierpoint", p)
})


test_that("plot_rt_vs_chain no outlier report", {
  expect_message(
    p <- plot_rt_vs_chain(
      mexp, 
      qc_types = "SPL", 
      x_var = "total_c", outlier_print = TRUE
    ),
    "potential annotation outliers: Hex2Cer d18:1/16:0 d3", fixed = TRUE)
  
  expect_no_message(
    p <- plot_rt_vs_chain(
      mexp, 
      qc_types = "SPL", 
      x_var = "total_c", outlier_print = FALSE
    ))
  
})


test_that("plot_rt_vs_chain no outlier report", {
  df <- tibble::tibble(
    feature = c("LPC 18:1", "LPC 20:1", "LPC 22:1",
                "LPE 18:1", "LPE 20:1", "LPE 22:1"),
    value   = c(1.1, 2.2, 2.8,
                2.1, 3.2, 4.8)
  )
  expect_message(
    p <- plot_rt_vs_chain(
      df, 
      qc_types = "SPL", 
      x_var = "total_c"),
      "No potential annotation outliers were detected", fixed = TRUE)
  
  vdiffr::expect_doppelganger(" plot_rtr_vs_chain tibble data", p)
})