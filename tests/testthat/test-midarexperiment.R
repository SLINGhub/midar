testthat::test_that("import_MSOrganizerXLM: Template with all information", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::loadMasshunterCSV(mexp, "21_Test_MH.csv")
  mexp <- midar::loadMSOrganizerXLM(mexp, "20_MSTemplate_Creator_forTest.xlsm")

  mexp <- midar::normalize_by_istd(mexp)
  mexp <- midar::quantitate_by_istd(mexp)
  mexp <- midar::calculate_qc_metrics(mexp)

  dd <- readRDS("21_MidarExperiment_1.rds")
expect_equal(all.equal(mexp@d_QC, dd@d_QC), TRUE) &
  expect_equal(all.equal(mexp@dataset, dd@dataset), TRUE)
})

