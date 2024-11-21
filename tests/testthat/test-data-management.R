test_that("Orders analyses according to set criteria and absence/presence of timestamp", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("23_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq_notimestamp.csv"),
                                        include_metadata = FALSE)
  mexp <- midar::import_metadata_midarxlm(mexp,
                                          path = testthat::test_path("MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1_reorder.xlsm"),
                                          excl_unmatched_analyses = FALSE)

  #Orders according to the order in the input csv file (no timestamp available)
  expect_equal(mexp@dataset[[1, "analysis_id"]], "020_SPL_S001")

  # Error as no timestamp available
  expect_error(midar::set_analysis_order(mexp, order_by = "timestamp"), "Acquisition timestamps are not present")

  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("22_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq.csv"),
                                        include_metadata = FALSE)
  mexp <- midar::import_metadata_midarxlm(mexp,
                                          path = testthat::test_path("MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1_reorder.xlsm"),
                                          excl_unmatched_analyses = FALSE)
  #Orders by default according to the timestamp, if available
  expect_equal(mexp@dataset[[1, "analysis_id"]], "006_EBLK_Extracted Blank+ISTD01")
  mexp <- midar::set_analysis_order(mexp, order_by = "metadata")
  expect_equal(mexp@dataset[[1, "analysis_id"]], "207_SOLV_Blank02")
  mexp <- midar::set_analysis_order(mexp, order_by = "resultfile")
  expect_equal(mexp@dataset[[1, "analysis_id"]], "020_SPL_S001")
  mexp <- midar::set_analysis_order(mexp, order_by = "timestamp")
  expect_equal(mexp@dataset[[1, "analysis_id"]], "006_EBLK_Extracted Blank+ISTD01")
})
