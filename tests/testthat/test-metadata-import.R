# library(testthat)
# library(tibble)
# library(cli)
# library(dplyr)

test_that("Imports/associates data and metadata, orders analyses by dataset (timestamp missing)", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("23_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq_notimestamp.csv"),
                                     import_metadata = FALSE)
  mexp <- midar::import_metadata_msorganiser(mexp,
    path = testthat::test_path("testdata/MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1.xlsm"),
    excl_unmatched_analyses = FALSE
  )
  expect_equal(mexp@dataset[[1, "analysis_id"]], "020_SPL_S001")
  expect_equal(dim(mexp@annot_analyses),c(65, 13))
  expect_equal(dim(mexp@annot_features),c(16, 17))
  expect_equal(dim(mexp@annot_istds),c(2, 5))
  expect_equal(dim(mexp@annot_responsecurves),c(12, 5))
  expect_equal(dim(mexp@annot_batches),c(3, 4))
  expect_equal(dim(mexp@annot_qcconcentrations),c(6, 6))
  expect_in(c("analysis_order", "batch_id", "is_quantifier", "qc_type", "is_istd", "is_quantifier"), names(mexp@dataset))
})

test_that("import_metadata_msorganiser handles missing / invalid files", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("22_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq-noalphafeat.csv"), import_metadata = FALSE)
  expect_error(mexp <- midar::import_metadata_msorganiser(mexp,
                                                          path = testthat::test_path("Does not exist.xlsm"),
                                                          excl_unmatched_analyses = FALSE),
               "File not found. Please verify path"
  )

  expect_error(mexp <- midar::import_metadata_msorganiser(mexp,
                                                          path = testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_DefaultSampleInfo_RT-Areas-FWHM.csv"),
                                                          excl_unmatched_analyses = FALSE),
               "Invalid file type not. A MSOrganiser template", fixed = TRUE
  )

  expect_error(mexp <- midar::import_metadata_msorganiser(mexp,
                                                          path = testthat::test_path("testdata/MiDAR_Metadata_Template_wrongabout.xlsm"),
                                                          excl_unmatched_analyses = FALSE),
               "This appears to be an invalid or unsupported MSOrganiser template file.", fixed = TRUE
  )

  expect_error(mexp <- midar::import_metadata_msorganiser(mexp,
                                                          path = testthat::test_path("testdata/MiDAR_Metadata_Template_missingabout.xlsm"),
                                                          excl_unmatched_analyses = FALSE),
               "This appears to be an invalid or unsupported MSOrganiser template file, without `About`", fixed = TRUE
  )

   expect_error(mexp <- midar::import_metadata_msorganiser(mexp,
                                                           path = testthat::test_path("testdata/MiDAR_Metadata_Template_lowerversion.xlsm"),
                                                           excl_unmatched_analyses = FALSE),
                "Unsupported MSOrganiser template version. Please use an MSOrganiser template v0.2 or higher.", fixed = TRUE
  )
   expect_error(mexp <- midar::import_metadata_msorganiser(mexp,
                                                           path = testthat::test_path("testdata/MiDAR_Metadata_Template_higherversion.xlsm"),
                                                           excl_unmatched_analyses = FALSE),
                "Unsupported MSOrganiser template version. Please use an MSOrganiser template v0.2 or higher.", fixed = TRUE
   )

   expect_error(mexp <- midar::import_metadata_msorganiser(mexp,
                                                           path = testthat::test_path("testdata/MiDAR_Metadata_Template_invalidversion.xlsm"),
                                                           excl_unmatched_analyses = FALSE),
                "Invalid version number found in the template. Please", fixed = TRUE
   )
})


test_that("Imports/associates data and metadata, orders features by default according to order in metadata", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("22_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq-noalphafeat.csv"), import_metadata = FALSE)
  mexp <- midar::import_metadata_msorganiser(mexp,
    path = testthat::test_path("testdata/MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1.xlsm"),
    excl_unmatched_analyses = FALSE
  )
  expect_equal(mexp@dataset_orig[1, ] |> pull("feature_id"), "S1P d18:0 [M>113]")
  expect_equal(mexp@dataset[1, ] |> pull("feature_id"), "S1P d16:1 [M>60]")
})


test_that("Raise data assertion warning and stops with not all analyses defined in analysis metadata ", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("22_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq-noalphafeat.csv"), import_metadata = FALSE)
  expect_error(mexp <- midar::import_metadata_msorganiser(mexp,
                                          path = testthat::test_path("MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_analysissubset.xlsm"),
                                          excl_unmatched_analyses = FALSE),
               "Not all analyses listed"
  )
})

test_that("Shows analyses defined in metadata but missing in data as Note in assertion table, instead of Warning and proceeds", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("22_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq-noalphafeat.csv"), import_metadata = FALSE)
  mexp <- midar::import_metadata_msorganiser(mexp,
                                          path = testthat::test_path("MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_analysissubset.xlsm"),
                                          excl_unmatched_analyses = TRUE)

  expect_equal(mexp@dataset_orig[1, ] |> pull("feature_id"), "S1P d18:0 [M>113]")
  expect_equal(mexp@dataset[1, ] |> pull("feature_id"), "S1P d16:1 [M>60]")
})

test_that("Ignores warnings after metadata import and proceeds", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("22_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq-noalphafeat.csv"), import_metadata = FALSE)
  mexp <- midar::import_metadata_msorganiser(mexp,
                                          path = testthat::test_path("MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_analysissubset.xlsm"),
                                          ignore_warnings = TRUE)

  expect_equal(mexp@dataset_orig[1, ] |> pull("feature_id"), "S1P d18:0 [M>113]")
  expect_equal(mexp@dataset[1, ] |> pull("feature_id"), "S1P d16:1 [M>60]")
})


test_that("Reads metadata from csv file", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.csv")
  tbl <- get_metadata_table(path = path)
  expect_equal(tbl[[1,"analysis_id"]], "Longit_BLANK-01 (Eluent A)")
  expect_equal(tbl[[1,"sample_amount"]], 10)
})

test_that("Reads metadata from XLSX sheet", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.xlsx")
  tbl <- get_metadata_table(path = path, sheet = "Analyses")
  expect_equal(tbl[[1,"analysis_id"]], "Longit_BLANK-01 (Eluent A)")
  expect_equal(tbl[[1,"sample_amount"]], 10)
})

test_that("Reads metadata from XLSX sheet", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.xlsx")
  tbl <- get_metadata_table(path = path, sheet = "Features")
  expect_equal(tbl[[1,"feature_id"]], "CE 14:0")
  expect_equal(tbl[[224,"interference_contribution"]], 0.00774513)
})

test_that("Reads metadata from given data frame", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.csv")
  df <- readr::read_csv(file = path, show_col_types = FALSE)
  tbl <- get_metadata_table(dataset = df)
  expect_equal(tbl[[1,"analysis_id"]], "Longit_BLANK-01 (Eluent A)")
  expect_equal(tbl[[1,"sample_amount"]], 10)
})


test_that("Reads metadata from given data frame", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.csv")
  df <- readr::read_csv(file = path, show_col_types = FALSE)
  expect_error(get_metadata_table(path = path, dataset = df), regexp = "cannot be specified at the same time")
})


test_that("Prepare analysis metadata from given data file", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.csv")
  tbl <- get_metadata_table(path = path) |> dplyr::select(-"batch_id")
  metadata <- clean_analysis_metadata(tbl)
  expect_in(c("batch_id", "replicate_no", "valid_analysis"), names(metadata))
  expect_equal(metadata[[1,"valid_analysis"]], TRUE)
})

test_that("Prepare analysis metadata from given data file", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.csv")
  tbl <- get_metadata_table(path = path) |> select(-"batch_id")
  metadata <- clean_analysis_metadata(tbl)
  expect_in(c("batch_id", "replicate_no", "valid_analysis"), names(metadata))
  expect_equal(metadata[[1,"valid_analysis"]], TRUE)
  expect_error(clean_analysis_metadata(tbl |> select(-"analysis_id")),
               regexp = "Analysis (Sample) metadata must have the `analysis_id` columnn",
  fixed = TRUE)
})

test_that("Prepare feature metadata from given table imported from an XLSX sheet", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.xlsx")
  tbl <- get_metadata_table(path = path, sheet = "Features")
  metadata <- clean_feature_metadata(tbl)
  expect_type(metadata$response_factor, "double")
  expect_type(metadata$interference_contribution, "double")
  expect_type(metadata$is_quantifier, "logical")
  expect_type(metadata$valid_feature, "logical")
  expect_equal(metadata[[1,"feature_id"]], "CE 14:0")
  expect_error(clean_istd_metadata(tbl |> select(-"feature_id")), regexp = "must have following columns")
})

test_that("Prepare istd metadata from given table imported from an XLSX sheet", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.xlsx")
  tbl <- get_metadata_table(path = path, sheet = "ISTDs")
  metadata <- clean_istd_metadata(tbl)
  expect_type(metadata$quant_istd_feature_id, "character")
  expect_type(metadata$istd_conc_nmolar, "double")
  expect_type(metadata$remarks, "character")
  expect_equal(metadata[[1,"quant_istd_feature_id"]], "CE 18:1 d7 (ISTD)")
  expect_equal(nrow(metadata), 19)
  expect_error(clean_istd_metadata(tbl |> select(-"istd_feature_id")), regexp = "must have following columns")
})

test_that("Prepare rqc metadata from given table imported from an XLSX sheet", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.xlsx")
  tbl <- get_metadata_table(path = path, sheet = "RQCs")
  metadata <- clean_response_metadata(tbl)
  expect_type(metadata$analysis_id, "character")
  expect_type(metadata$analyzed_amount, "double")
  expect_type(metadata$analyzed_amount_unit, "character")
  expect_type(metadata$remarks, "character")
  expect_equal(metadata[[1,"analysis_id"]], "Longit_TQC-10%")
  expect_equal(nrow(metadata), 12)
  expect_error(clean_response_metadata(tbl |> select(-"curve_id")), regexp = "must have following columns")
})

test_that("Prepare qc concentration metadata from given table imported from an XLSX sheet", {
  path <- testthat::test_path("testdata/sperfect_metadata_tables.xlsx")
  tbl <- get_metadata_table(path = path, sheet = "QCconc")
  metadata <- clean_qcconc_metadata(tbl)
  expect_type(metadata$sample_id, "character")
  expect_type(metadata$analyte_id, "character")
  expect_type(metadata$concentration, "double")
  expect_type(metadata$concentration_unit, "character")
  expect_type(metadata$remarks, "character")
  expect_equal(metadata[[1,"sample_id"]], "199_NIST_NIST04")
  expect_equal(metadata[[2,"sample_id"]], "198_LTR_LTR04")
  expect_equal(metadata[[3,"concentration"]], 0.1)
  expect_equal(nrow(metadata), 6)
  expect_equal(ncol(metadata), 6)
  expect_error(clean_qcconc_metadata(tbl |> select(-"analyte_id")), regexp = "must have following columns")
})

test_that("Add indidual metadata types to data, first analyses then features", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(
    mexp,
    path = testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_DefaultSampleInfo_RT-Areas-FWHM.csv"),
    import_metadata = FALSE)
  path <- testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_metadata_tables.xlsx")
  mexp <- midar:::import_metadata_analyses(mexp, path = path, sheet = "Analyses", excl_unmatched_analyses = TRUE)
  expect_equal(mexp@dataset[[10, "feature_intensity" ]], 43545)
  expect_in(c("sample_id", "feature_class", "is_istd", "is_quantifier"), names(mexp@dataset))
  mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features",ignore_warnings = TRUE)
  #expect_equal(mexp@annot_features[[245, "interference_contribution" ]], 43545)
  expect_equal(mexp@dataset[[5, "feature_intensity" ]], 43545)
  expect_equal(mexp@dataset[[5, "feature_class" ]], "SPBP")
  expect_equal(mexp@dataset[[5, "is_istd" ]], TRUE)
  expect_in(c("sample_id", "feature_class", "is_istd", "is_quantifier"), names(mexp@dataset))

})


test_that("Add indidual metadata types to data, first features then analyses", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(
    mexp,
    path = testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_DefaultSampleInfo_RT-Areas-FWHM.csv"),
    import_metadata = FALSE)
  path <- testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_metadata_tables.xlsx")
  mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features",ignore_warnings = TRUE)
  expect_equal(mexp@dataset[[5, "feature_intensity" ]], 43545)
  expect_equal(mexp@dataset[[5, "feature_class" ]], "SPBP")
  expect_true(is.na(mexp@dataset[[1, "batch_id" ]]))
  expect_true(is.na(mexp@dataset[[1, "qc_type" ]]))
  mexp <- midar:::import_metadata_analyses(mexp, path = path, sheet = "Analyses", excl_unmatched_analyses = TRUE,ignore_warnings = TRUE)
  expect_equal(mexp@dataset[[5, "feature_intensity" ]], 43545)
  expect_in(c("sample_id", "feature_class", "is_istd", "is_quantifier"), names(mexp@dataset))
  expect_equal(mexp@dataset[[5, "qc_type" ]], "PBLK")
  expect_equal(mexp@dataset[[5, "batch_id" ]], "1")

})

test_that("Check import of inconsitent metadata", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(
    mexp,
    path = testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_DefaultSampleInfo_RT-Areas-FWHM.csv"),
    import_metadata = FALSE)
  path <- testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_metadata_tables.xlsx")
  expect_error(mexp <- midar:::import_metadata_analyses(mexp, path = path, sheet = "Analyses_Inconsistent", excl_unmatched_analyses = TRUE),
               "`valid_analysis` is inconsistently defined")
  expect_error(mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features_Inconsist_Both"),
               "`valid_feature` is inconsistently defined, i.e., not for one or more features. Please")
  expect_error(mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features_Inconsist_qual"),
               "`is_quantifier` is inconsistently defined, i.e., not for one or more features. Please")
  expect_error(mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features_Inconsist_val"),
               "`valid_feature` is inconsistently defined, i.e., not for one or more features. Please")
  expect_error(mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features_inval_interf"))
  expect_error(mexp <- midar:::import_metadata_features(mexp, path = path, sheet = "Features_miss_interf"))
})

test_that("Replacing specific undefined metadata", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(
    mexp,
    path = testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_DefaultSampleInfo_RT-Areas-FWHM.csv"),
    import_metadata = FALSE)
  path <- testthat::test_path("testdata/MIDAR_TestData_MHQuant_S1P_metadata_tables.xlsx")
  mexp <- midar:::import_metadata_analyses(mexp, path = path, sheet = "Analyses_missing_val", excl_unmatched_analyses = TRUE)
  expect_true(all(mexp@annot_analyses$valid_analysis))

  mexp2 <- midar:::import_metadata_features(mexp, path = path, sheet = "Features_missing_val",ignore_warnings = TRUE)
  expect_true(all(mexp@annot_features$valid_feature))

  mexp2 <- midar:::import_metadata_features(mexp, path = path, sheet = "Features_missing_quan",ignore_warnings = TRUE)
  expect_true(all(mexp@annot_features$is_quantifier))

})
