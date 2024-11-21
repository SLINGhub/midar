test_that("Add metadata table by table, the normalize and quantify based on ISTD", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(
    mexp,
    path = testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_DefaultSampleInfo_RT-Areas-FWHM.csv"),
    include_metadata = FALSE)
  path <- testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_metadata_tables.xlsx")
  mexp <- midar:::import_metadata_analyses(mexp, path = path, sheet = "Analyses", ignore_warnings = FALSE, excl_unmatched_analyses = TRUE)
  mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features", ignore_warnings = TRUE)
  mexp <- midar:::import_metadata_istds(mexp, path = path, sheet = "ISTDs", ignore_warnings = FALSE)

})
