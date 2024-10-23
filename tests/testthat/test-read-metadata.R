test_that("Orders samples in daraset without timestamp field according to parameter analysis_sequence", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::data_import_agilent(mexp, path = testthat::test_path("23_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq_notimestamp.csv"),
                                     use_metadata = FALSE)
  mexp <- midar::metadata_import_midarxlm(mexp,
    path = testthat::test_path("MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1.xlsm"),
    excl_unannotated_analyses = FALSE
  )
  expect_equal(mexp@dataset[1, ] |> pull("analysis_id"), "020_SPL_S001")
})

test_that("Orders features by default according to order in metadata", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::data_import_agilent(mexp, path = testthat::test_path("22_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq-noalphafeat.csv"), use_metadata = FALSE)
  mexp <- midar::metadata_import_midarxlm(mexp,
    path = testthat::test_path("MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1.xlsm"),
    excl_unannotated_analyses = FALSE
  )
  expect_equal(mexp@dataset_orig[1, ] |> pull("feature_id"), "S1P d18:0 [M>113]")
  expect_equal(mexp@dataset[1, ] |> pull("feature_id"), "S1P d16:1 [M>60]")
})
