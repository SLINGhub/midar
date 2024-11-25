test_that("Add metadata table by table, the normalize and quantify based on ISTD", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(
    mexp,
    path = testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_DefaultSampleInfo_RT-Areas-FWHM.csv"),
    import_metadata = FALSE)
  path <- testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_metadata_tables.xlsx")
  mexp <- midar:::import_metadata_analyses(mexp, path = path, sheet = "Analyses", ignore_warnings = FALSE, excl_unmatched_analyses = TRUE)
  mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features", ignore_warnings = TRUE)
  mexp <- midar:::import_metadata_istds(mexp, path = path, sheet = "ISTDs", ignore_warnings = FALSE)
  mexp <- midar::normalize_by_istd(mexp)
  expect_equal(mexp@dataset[[5, "feature_intensity" ]], 43545)
  expect_equal(mexp@dataset[[434, "feature_norm_intensity"]], 0.64371386)
  expect_equal(mexp@dataset[[583, "feature_norm_intensity"]], 1.0)  # ISTD norm by itself
  expect_equal(mexp@dataset[[583, "is_istd"]], TRUE) # ISTD norm by itself
  expect_equal(mexp@dataset[[583, "is_istd"]], TRUE) # ISTD norm by itself
  mexp <- midar::quantify_by_istd(mexp)
  expect_equal(mexp@dataset[[583, "feature_pmol_total"]], 4.0)
  expect_equal(mexp@dataset[[583, "feature_conc"]], 0.2)  # ISTD norm by itself
  expect_equal(mexp@dataset[[582, "feature_pmol_total"]], 15.6093933)
  expect_equal(mexp@dataset[[582, "feature_conc"]], 0.78046967)  # ISTD norm by itself
})
