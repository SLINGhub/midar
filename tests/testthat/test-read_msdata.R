testthat::test_that("read_MassHunterCSV: 1_DefaultSampleInfo_AreaOnly", {

    d <- midar::read_MassHunterCSV("1_Testdata_MHQuant_DefaultSampleInfo_AreaOnly.csv")
    dd <- readRDS("1_Testdata_MHQuant_DefaultSampleInfo_AreaOnly.rds")
    expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 2_DefaultSampleInfo_RT-Areas-FWHM.csv", {

  d <- midar::read_MassHunterCSV("2_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM.csv")
  dd <- readRDS("2_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM.rds")
  expect_identical(d, dd)
})


testthat::test_that("read_MassHunterCSV: 3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults.csv", {

  d <- midar::read_MassHunterCSV("3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults.csv")
  dd <- readRDS("3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults.rds")
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.rds", {

  d <- midar::read_MassHunterCSV("4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.csv")
  dd <- readRDS("4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.rds")
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 5_Testdata_MHQuant_DetailedSampleInfo-RT-Areas-FWHM.csv", {

  d <- midar::read_MassHunterCSV("5_Testdata_MHQuant_DetailedSampleInfo-RT-Areas-FWHM.csv")
  dd <- readRDS("5_Testdata_MHQuant_DetailedSampleInfo-RT-Areas-FWHM.rds")
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 6_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoOutlierSum.csv", {

  d <- midar::read_MassHunterCSV("6_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoOutlierSum.csv")
  dd <- readRDS("6_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoOutlierSum.rds")
  expect_identical(d, dd)
})


testthat::test_that("read_MassHunterCSV: 7_Testdata_MHQuant_NoOutlierSum-noQuantMsgSum", {

  d <- midar::read_MassHunterCSV("7_Testdata_MHQuant_NoOutlierSum-noQuantMsgSum.csv")
  dd <- readRDS("7_Testdata_MHQuant_NoOutlierSum-noQuantMsgSum.rds")
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 8_Testdata_MHQuant_Corrupt_OutlierQuantMsgSumDeleted", {

  d <- midar::read_MassHunterCSV("8_Testdata_MHQuant_Corrupt_OutlierQuantMsgSumDeleted.csv")
  dd <- readRDS("8_Testdata_MHQuant_Corrupt_OutlierQuantMsgSumDeleted.rds")
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 9_Testdata_MHQuant_withQuantMethods_withQualifierMethResults", {

  d <- midar::read_MassHunterCSV("9_Testdata_MHQuant_withQuantMethods_withQualifierMethResults.csv")
  dd <- readRDS("9_Testdata_MHQuant_withQuantMethods_withQualifierMethResults.rds")
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 10_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoQuantMsgSum", {

  d <- midar::read_MassHunterCSV("10_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoQuantMsgSum.csv")
  dd <- readRDS("10_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoQuantMsgSum.rds")
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 11_Testdata_MHQuant_DefaultSampleInfo-noAcqDataTime_RT-Areas-FWHM", {

  d <- midar::read_MassHunterCSV("11_Testdata_MHQuant_DefaultSampleInfo-noAcqDataTime_RT-Areas-FWHM.csv")
  dd <- readRDS("11_Testdata_MHQuant_DefaultSampleInfo-noAcqDataTime_RT-Areas-FWHM.rds")
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 12_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-withQuantMsg.", {
  expect_error(
    midar::read_MassHunterCSV("12_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-withQuantMsg.csv"), "Quantitation Message")

})

testthat::test_that("read_MassHunterCSV: 13_Testdata_MHQuant_CompoundTable_DefaultSampleInfo_RT-Areas-FWHM.", {
  expect_error(
    midar::read_MassHunterCSV("13_Testdata_MHQuant_CompoundTable_DefaultSampleInfo_RT-Areas-FWHM.csv"),
    regexp = "Compound table")

})

testthat::test_that("read_MassHunterCSV: 14_Testdata_MHQuant_Corrupt_RowAreaDeleted", {
  expect_error(
    midar::read_MassHunterCSV("14_Testdata_MHQuant_Corrupt_RowAreaDeleted.csv"),
    regexp = "Unknown format")

})

testthat::test_that("read_MassHunterCSV: 15_Testdata_MHQuant_Corrupt_ExtraTopLine", {
  d <- midar::read_MassHunterCSV("15_Testdata_MHQuant_Corrupt_ExtraTopLine.csv")
  dd <- readRDS("15_Testdata_MHQuant_Corrupt_ExtraTopLine.rds")
  expect_identical(d, dd)

})

testthat::test_that("read_MassHunterCSV_wide: 3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults-Area-flat.rds", {
  d <- midar::read_MassHunterCSV_wide("3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults.csv", field = "Area")
  dd <- readRDS("3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults-Area-flat.rds")
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 17_Testdata_Lipidomics_GermanSystem", {
  d <- midar::read_MassHunterCSV("17_Testdata_Lipidomics_GermanSystem.csv")
  dd <- readRDS("17_Testdata_Lipidomics_GermanSystem.rds")
  print(dd)
  expect_identical(d, dd)
})

testthat::test_that("read_MassHunterCSV: 18_AMPM", {
  d <- midar::read_MassHunterCSV("18_Testdata_MHQuant_DefaultSampleInfo_AreaOnly_AMPM.csv")
  dd <- readRDS("18_Testdata_MHQuant_DefaultSampleInfo_AreaOnly_AMPM.rds")
  print(dd)
  expect_identical(d, dd)
})
