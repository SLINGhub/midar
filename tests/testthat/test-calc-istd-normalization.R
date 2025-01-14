test_that("istd-based normalization is correct and overwites previous if present", {
  mexp <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))

  testthat::expect_message(mexp <- normalize_by_istd(mexp,fail_missing_annotation = FALSE), "14 features normalized with 2 ISTDs in 65 analyses")

  # istd normalized by itself should be 1
  testthat::expect_equal(mexp@dataset$feature_norm_intensity[5], 1)

  testthat::expect_equal(mexp@dataset$feature_norm_intensity[100], 3.31972989)

  # Repeated normalization should give a notification
  testthat::expect_message(mexp <- normalize_by_istd(mexp,fail_missing_annotation = FALSE), "Replacing previously normalized feature intensities")
  testthat::expect_message(mexp <- normalize_by_istd(mexp,fail_missing_annotation = FALSE), "14 features normalized with 2 ISTDs in 65 analyses")
})

test_that("istd-based normalization handles missing info correctly", {
  mexp <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))

  mexp_mod <- mexp
  mexp_mod@annot_features <- mexp_mod@annot_features[-13,]
  testthat::expect_error(mexp_mod <- normalize_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                         "1 ISTD\\(s\\) were not defined as individual feature")

  mexp_mod <- mexp
  mexp_mod@annot_features$istd_feature_id <- NA
  testthat::expect_error(mexp_mod <- normalize_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                         "No ISTDs defined in feature metadata")

  mexp_mod <- mexp
  mexp_mod@annot_features$istd_feature_id[1] <- NA
  testthat::expect_message(mexp_mod <- normalize_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                         "For 1 feature\\(s\\) no ISTD was defined, normalized intensities will be \\`NA")
  testthat::expect_message(mexp_mod <- normalize_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                           "13 features normalized with 2 ISTDs in 65")

  testthat::expect_error(mexp_mod <- normalize_by_istd(mexp_mod,fail_missing_annotation = TRUE),
                         "For 1 feature\\(s\\) no ISTD was defined. Please ammend feature")

})


test_that("Add metadata table by table, the normalize and quantify based on ISTD", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(
    mexp,
    path = testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_DefaultSampleInfo_RT-Areas-FWHM.csv"),
    import_metadata = FALSE)
  path <- testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_metadata_tables.xlsx")
  expect_message(mexp <- midar:::import_metadata_analyses(mexp, path = path, sheet = "Analyses", ignore_warnings = FALSE, excl_unmatched_analyses = TRUE),
                 "Analysis metadata associated with 64 analyses")
  expect_message(mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features", ignore_warnings = TRUE),
                 "Feature metadata associated with 15 features")
  expect_message(mexp <- midar:::import_metadata_istds(mexp, path = path, sheet = "ISTDs", ignore_warnings = FALSE),
                 "Internal Standard metadata associated with 2 ISTDs")
  expect_message(mexp <- midar:::import_metadata_responsecurves(mexp, path = path, sheet = "RQCs", ignore_warnings = FALSE),
                 "Response curve metadata associated with 12 analyses")
  expect_message(mexp <- midar:::import_metadata_qcconcentrations(mexp, path = path, sheet = "QCconc", ignore_warnings = FALSE),
                 "QC concentration metadata associated with 2 QC samples and 3 features")

  mexp <- midar::normalize_by_istd(mexp)
  expect_equal(mexp@dataset[[5, "feature_intensity" ]], 43545)
  expect_equal(mexp@dataset[[434, "feature_norm_intensity"]], 0.64371386)
  expect_equal(mexp@dataset[[583, "feature_norm_intensity"]], 1.0)  # ISTD norm by itself
  expect_equal(mexp@dataset[[583, "is_istd"]], TRUE) # ISTD norm by itself
  expect_equal(mexp@dataset[[583, "is_istd"]], TRUE) # ISTD norm by itself
  mexp <- midar::quantify_by_istd(mexp)
  # expect_equal(mexp@dataset[[583, "feature_pmol_total"]], 4.0)
  expect_equal(mexp@dataset[[583, "feature_conc"]], 0.2)  # ISTD norm by itself
  expect_equal(mexp@dataset[[582, "feature_pmol_total"]], 15.6093933)
  expect_equal(mexp@dataset[[582, "feature_conc"]], 0.78046967)  # ISTD norm by itself

  expect_equal(dim(mexp@annot_responsecurves), c(12,5))
  expect_equal(dim(mexp@annot_qcconcentrations), c(6,5))

})



test_that("istd-based quantification is correct and overwites previous if present", {
  mexp <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))
  mexp <- normalize_by_istd(mexp,fail_missing_annotation = FALSE)

  testthat::expect_message(mexp <- quantify_by_istd(mexp,fail_missing_annotation = FALSE),
                           "14 feature concentrations calculated based on 2 ISTDs and sample amounts of 65 analyses")

  testthat::expect_message(mexp <- quantify_by_istd(mexp,fail_missing_annotation = FALSE),
                           "Concentrations are given in μmol/L")

  # istd normalized by itself and quantified should be istd conc corrected for the dilution factor
  testthat::expect_equal(mexp@dataset$feature_conc[5], 0.2)

  testthat::expect_equal(mexp@dataset$feature_conc[100], 0.663945978 )

  # Repeated normalization should give a notification
  testthat::expect_message(mexp <- quantify_by_istd(mexp,fail_missing_annotation = FALSE),
                           "Replacing previously calculated concentrations")
  testthat::expect_message(mexp <- quantify_by_istd(mexp,fail_missing_annotation = FALSE),
                           "14 feature concentrations calculated based on 2 ISTDs and sample amounts of 65 analyses")
})



test_that("istd-based quantification handles missing info correctly", {
  mexp <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))
  mexp <- normalize_by_istd(mexp,fail_missing_annotation = FALSE)

  mexp_mod <- mexp
  mexp_mod@annot_analyses$sample_amount[11] <- NA
  testthat::expect_message(mexp_mod <- quantify_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                         "Sample and/or ISTD solution amount\\(s\\) for 1 analyses missing, concentrations of all features for these analyses will be \\`NA")
  testthat::expect_message(mexp_mod <- quantify_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                           "14 feature concentrations calculated based on 2 ISTDs and sample amounts of 64 analyses")


   testthat::expect_error(mexp_mod <- quantify_by_istd(mexp_mod,fail_missing_annotation = TRUE),
                           "Sample and\\/or ISTD amount\\(s\\) for 1 analyses missing. Please ammend")


  mexp_mod <- mexp
  mexp_mod@annot_istds <- mexp_mod@annot_istds[-2,]
  testthat::expect_error(mexp_mod <- quantify_by_istd(mexp_mod,fail_missing_annotation = TRUE),
                         "Concentrations of 1 ISTD\\(s\\) missing. Please ammend ISTD")

  testthat::expect_message(mexp_mod <- quantify_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                         "Spiked-in concentrations of 1 ISTD\\(s\\) missing, calculated concentrations of affected features will be \\`NA")

  mexp_mod <- mexp
  mexp_mod@annot_features$istd_feature_id <- NA
  testthat::expect_error(mexp_mod <- normalize_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                         "No ISTDs defined in feature metadata")

  mexp_mod <- mexp
  mexp_mod@annot_features$istd_feature_id[1] <- NA
  testthat::expect_message(mexp_mod <- normalize_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                           "For 1 feature\\(s\\) no ISTD was defined, normalized intensities will be \\`NA")
  testthat::expect_message(mexp_mod <- normalize_by_istd(mexp_mod,fail_missing_annotation = FALSE),
                           "13 features normalized with 2 ISTDs in 65")

  testthat::expect_error(mexp_mod <- normalize_by_istd(mexp_mod,fail_missing_annotation = TRUE),
                         "For 1 feature\\(s\\) no ISTD was defined. Please ammend feature")

})
