testthat::test_that("import_msorganizer_xlm: Template with all information", {
  # mexp <- midar::MidarExperiment()
  # mexp <- midar::data_import_agilent(mexp, "21_Test_MH.csv")
  # mexp <- midar::metadata_import_midarxlm(mexp, "20_MSTemplate_Creator_forTest.xlsm", excl_unannotated_analyses = TRUE)
  #
  # mexp <- midar::calc_normalize_by_istd(mexp)
  # mexp <- midar::calc_quant_by_istd(mexp)
  # mexp <- midar::qc_calc_metrics(mexp)
  #
  # dd <- readRDS("21_MidarExperiment_1.rds")
  # #expect_equal(all.equal(mexp@metrics_qc, dd@metrics_qc), TRUE) &
  # #expect_equal(all.equal(mexp@dataset, dd@dataset), TRUE)
  #
  # announce_snapshot_file(name = "mexp")
  #
  # r_to_rds <- function(x) {
  #   path <- tempfile(fileext = "red")
  #   readr::write_rds(x, path)
  #   path
  # }
  #
  # expect_snapshot(mexp,)
  # expect_snapshot_file(r_to_rds(mexp), "mexp.rds")
  # testthat::snapshot_review(path = "mexp.rds")
})
