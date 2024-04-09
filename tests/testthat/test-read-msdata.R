testthat::test_that("Parses basic Agilent MH-Quant .csv file with only peak areas", {
 d <- read_masshunter_csv(testthat::test_path("1_Testdata_MHQuant_DefaultSampleInfo_AreaOnly.csv"))

 expect_contains(names(d), c("file_run_id", "raw_data_filename", "sample_name", "sample_type", "acquisition_time_stamp", "feature_name", "feature_area"))
 expect_equal(nrow(d), 1040)
 expect_equal(mean(d$feature_area, na.rm = TRUE), 17237.244)
 expect_contains(unname(unlist(lapply(d, \(x) class(x)[[1]]))), c("integer", "character", "character", "character", "POSIXct", "character", "numeric"))
})



testthat::test_that("Parses nested Agilent MH-Quant .csv file with diverse peak variables", {
 d <- read_masshunter_csv(testthat::test_path("3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults.csv"))

 expect_identical(names(d), c("file_run_id", "raw_data_filename", "sample_name", "sample_type", "acquisition_time_stamp", "feature_name", "integration_qualifier", "feature_rt", "feature_area",
               "feature_fwhm", "feature_height", "feature_int_start", "feature_int_end", "feature_sn_ratio", "feature_symetry", "feature_width", "feature_manual_integration"))
 expect_equal(nrow(d), 1040)
 expect_equal(mean(d$feature_area, na.rm = TRUE), 17237.244)
 expect_identical(unname(unlist(lapply(d, \(x) class(x)[[1]]))), c(c("integer", "character", "character", "character", "POSIXct", "character","logical"), rep("numeric",9), "logical"))
})

testthat::test_that("Parses nested Agilent MH-Quant .csv file with detailed sample info and different peak parameters", {
 d <- read_masshunter_csv(testthat::test_path("5_Testdata_MHQuant_DetailedSampleInfo-RT-Areas-FWHM.csv"))

 expect_identical(names(d), c("file_run_id", "raw_data_filename", "sample_name", "sample_group", "sample_type", "acquisition_time_stamp", "inj_volume", "comment", "completed",
  "dilution_factor", "instrument_name", "instrument_type", "acq_method_file", "acq_method_path", "data_file_path", "feature_name", "integration_qualifier", "feature_rt", "feature_area", "feature_fwhm"))
 expect_equal(nrow(d), 1040)
 expect_equal(mean(d$feature_area, na.rm = TRUE), 17237.244)
 expect_identical(unname(unlist(lapply(d, \(x) class(x)[[1]]))), c(
  "integer", "character", "character", "character", "character", "POSIXct", "character", "character", "character",
  "character", "character", "character", "character", "character", "character", "character", "logical", "numeric", "numeric", "numeric"))
})

testthat::test_that("Parses nested MH Quant .csv file with detailed method info and different peak parameters", {
 d <- read_masshunter_csv(testthat::test_path("4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.csv"))

 expect_identical(names(d), c("file_run_id","raw_data_filename","sample_name","sample_type","acquisition_time_stamp","feature_name","integration_qualifier", "method_compound_group",
                "method_collision_energy","method_fragmentor","method_compound_id","method_integration_method","method_integration_parameters","method_polarity","method_ion_source",
                "method_multiplier","method_noise_algorithm","method_noise_raw_signal","method_precursor_mz","method_product_mz","method_peak_smoothing","method_peak_smoothing_gauss_width",
                "method_peak_smoothing_function_width","method_transition","method_time_segment","method_type","feature_rt","feature_area","feature_fwhm" ))
 expect_equal(nrow(d), 1040)
 expect_equal(mean(d$feature_area, na.rm = TRUE), 17237.244)
 expect_identical(unname(unlist(lapply(d, \(x) class(x)[[1]]))), c(
  "integer","character","character","character","POSIXct","character", "logical","character","numeric","numeric","character","character","character","factor","character","numeric","character","numeric","numeric","numeric",
  "character","character","character","character","integer","character","numeric","numeric","numeric" ))
})

testthat::test_that(desc = "Console output with correct number of samples and features", code = {
 expect_output(read_masshunter_csv(testthat::test_path("4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.csv")),
               regexp = "Imported 65 samples with 16 features")
})


testthat::test_that("Parses nested MH Quant .csv without the 'outlier' column", {
 d <- read_masshunter_csv(testthat::test_path("6_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoOutlierSum.csv"))
 expect_equal(ncol(d), 10)
 expect_equal(nrow(d), 1040)
 expect_equal(names(d)[1], "file_run_id")
 expect_equal(d[1,"feature_rt", drop = TRUE], 3.422)
 expect_equal(d[1,"feature_area", drop = TRUE], 71)
})

testthat::test_that("Parses nested MH Quant .csv without the 'outlier' and 'quant message' column", {
 d <- read_masshunter_csv(testthat::test_path("7_Testdata_MHQuant_NoOutlierSum-noQuantMsgSum.csv"))
 expect_equal(ncol(d), 10)
 expect_equal(nrow(d), 1040)
 expect_equal(names(d)[1], "file_run_id")
 expect_equal(d[1,"feature_rt", drop = TRUE], 3.422)
 expect_equal(d[1,"feature_area", drop = TRUE], 71)
})

testthat::test_that("Returns error parsing corrupted nested MH Quant .csv with 'outlier'/'quant message' columns removed including header 'Samples'", {
  expect_error (read_masshunter_csv(testthat::test_path("8_Testdata_MHQuant_Corrupt_OutlierQuantMsgSumDeleted.csv")),
                regexp = "Error parsing this file\\. It may in unsupported format")
})


testthat::test_that("Parses nested MH Quant .csv file containing QUALIFIER peak info", {
  expect_output(d <- read_masshunter_csv(
    testthat::test_path("9_Testdata_MHQuant_withQuantMethods_withQualifierMethResults.csv"),
    expand_qualifier_names = TRUE),
    "Imported 65 samples with 16 features \\(8 quantifiers\\, 8 qualifiers\\)")
  expect_equal(ncol(d), 16)
  expect_equal(nrow(d), 1040)
  expect_equal(d[1,"feature_rt", drop = TRUE], 3.422)
  expect_equal(d[1,"feature_area", drop = TRUE], 51)
  expect_equal(d |> filter(integration_qualifier) |> pull(feature_name) |> dplyr::first(),"S1P d16:1 [M>60] [QUAL 408.3 -> 113.0]" )
  expect_equal(sum(d$integration_qualifier[d$file_run_id == 1]), 8)
})

testthat::test_that("Parses nested MH Quant .csv without Quant Message Summary", {
  d <- read_masshunter_csv(testthat::test_path("10_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoQuantMsgSum.csv"))
  expect_equal(ncol(d), 10)
  expect_equal(nrow(d), 1040)
  expect_equal(names(d)[1], "file_run_id")
  expect_equal(d[1,"feature_rt", drop = TRUE], 3.422)
  expect_equal(d[1,"feature_area", drop = TRUE], 71)
})

testthat::test_that("Parses nested MH Quant .csv without acquistion time stamp", {
  d <- read_masshunter_csv(testthat::test_path("11_Testdata_MHQuant_DefaultSampleInfo-noAcqDataTime_RT-Areas-FWHM.csv"))
  expect_false(c("acquisition_time_stamp") %in% names(d))
  expect_equal(ncol(d), 9)
  expect_equal(nrow(d), 1040)
  expect_equal(d[1,"feature_rt", drop = TRUE], 3.422)
  expect_equal(d[1,"feature_area", drop = TRUE], 71)
})

testthat::test_that("Returns a defined error when reading nested MH Quant .csv containg a 'Quantitation Message' is imported", {
 expect_error(
  read_masshunter_csv(testthat::test_path("12_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-withQuantMsg.csv")),
  regexp = "Field \\'Quantitation Message\\' currently not supported")
})

testthat::test_that("Returns a defined error when reading MH Quant .csv with analytes/features as rows (Compound Table) is imported", {
 expect_error(
  read_masshunter_csv(testthat::test_path("13_Testdata_MHQuant_CompoundTable_DefaultSampleInfo_RT-Areas-FWHM.csv")),
  regexp = "unsupported format\\, e\\.g\\. with features\\/analytes in rows")
})

testthat::test_that("Returns a defined error when reading a corrupted MH Quant .csv", {
 expect_error(
  read_masshunter_csv(testthat::test_path("14_Testdata_MHQuant_Corrupt_RowAreaDeleted.csv")),
  regexp = "Error parsing this Masshunter \\.csv file")
})

testthat::test_that("Parses nested MH Quant .csv file that has am empty first row", {
 d <- read_masshunter_csv(testthat::test_path("15_Testdata_MHQuant_Corrupt_ExtraTopLine.csv"))
 expect_equal(ncol(d), 16)
 expect_equal(nrow(d), 1040)
 expect_equal(names(d)[1], "file_run_id")
 expect_equal(d[1,"feature_rt", drop = TRUE], 3.422)
 expect_equal(d[1,"feature_area", drop = TRUE], 51)
})


testthat::test_that("Parses nested MH Quant .csv file exported from German Windows system with comma as decimal point", {
 d <- read_masshunter_csv(testthat::test_path("17_Testdata_Lipidomics_GermanSystem.csv"))
 expect_equal(d[1,"feature_rt", drop = TRUE], 9.754)
 expect_equal(d[1,"feature_fwhm", drop = TRUE], 0.056)
})

testthat::test_that("Parses nested MH Quant .csv file in UTF-8 format with different languages/characters", {
  d <- read_masshunter_csv(testthat::test_path("18_Testdata_Lipidomics_MultiLanguageCharactersSamplenamesFeatures.csv"))
  expect_equal(d[1,"feature_rt", drop = TRUE], 9.754)
  expect_equal(d[1,"feature_name", drop = TRUE], "谷氨酰胺")
  expect_equal(d[2,"feature_name", drop = TRUE], "글루타민")
  expect_equal(d[3,"feature_name", drop = TRUE], "Glutaminsäure")
  expect_equal(d[4,"feature_name", drop = TRUE], "குளுட்டமின்")
  expect_equal(d[1,"raw_data_filename", drop = TRUE], "Über_Schöner_Blank")
  expect_equal(d[300,"raw_data_filename", drop = TRUE], "空白的")
  expect_equal(d[600,"raw_data_filename", drop = TRUE], "공백2")
  expect_equal(d[900,"raw_data_filename", drop = TRUE], "空白")
  expect_equal(d[1200,"raw_data_filename", drop = TRUE], "வெற்று")
})

testthat::test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
  d <- read_masshunter_csv(testthat::test_path("19_Testdata_MHQuant_MultipleQUAL_with_expectedRT.csv"), expand_qualifier_names = TRUE)
  expect_equal(d[1,"feature_rt", drop = TRUE], 6.649)
  expect_equal(d[1,"method_target_rt", drop = TRUE], 7.200)
  })

testthat::test_that("Parses nested MH Quant .csv file with special characters (e.g. !@#$%^) in feature names", {
  d <- read_masshunter_csv(testthat::test_path("20_Testdata_MHQuant_withSpecialCharsInFeatures.csv"), expand_qualifier_names = TRUE)
  expect_equal(d[1,"feature_rt", drop = TRUE], 6.649)
  expect_equal(d[1,"method_target_rt", drop = TRUE], 7.200)
  expect_equal(d[3,"feature_name", drop = TRUE], "Analyte 2* 14:0")
  expect_equal(d[5,"feature_name", drop = TRUE], "Analyte 2~%^$# 14:0 (d5)")
  expect_match(d[10,"feature_name", drop = TRUE], "Analyte 3 16\\:0-_=\\+\\\\/~!@ \\[QUAL 317\\.3 -> 299\\.3\\]")
  expect_match(d[13,"feature_name", drop = TRUE], "Analyte 4 \\?><,\\.\`\\:\"\\}\\{\\[\\]18:0 \\[QUAL 331\\.3 -> 331\\.3\\]")
})


testthat::test_that("Parses nested MH Quant .csv with . \ | in feature names", {
  d <- read_masshunter_csv(testthat::test_path("21_Testdata_MHQuant_with_dots_InFeatures.csv"), expand_qualifier_names = TRUE)
  expect_equal(d[1,"feature_rt", drop = TRUE], 3.422)
  expect_match(d[3,"feature_name", drop = TRUE], "S1P d17\\.1\\|S1P d17\\.2 \\[M>113\\]")
  expect_match(d[4,"feature_name", drop = TRUE], "S1P d17\\.1\\\\S1P d17:2 \\[M>60\\]")
  expect_match(d[5,"feature_name", drop = TRUE], "S1P d18\\.0/S1P 18:1 \\[M>113\\]")
})

# testthat::test_that("Parses nested MH Quant .csv file and matches saved copy", {
#   d <- read_masshunter_csv(testthat::test_path("1_Testdata_MHQuant_DefaultSampleInfo_AreaOnly.csv"))
#   announce_snapshot_file(name = "newtest")
#   expect_snapshot(d)
# })

testthat::test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
  d <- read_analysisresult_table(testthat::test_path("001_Generic_Results_1.csv"),value_type = "area")
  expect_equal(d[1,"feature_area", drop = TRUE], 71)
  expect_equal(d[1,"feature_name", drop = TRUE], "S1P d16:1 [M>113]")
  expect_equal(d[1,"analysis_id", drop = TRUE], "006_EBLK_Extracted Blank+ISTD01")
})

testthat::test_that("Test mult", {
  d <- .test_mult(a,b)
  expect_equal(d,   rep("numeric", 3))
})


testthat::test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
  d <- read_analysisresult_table(testthat::test_path("001_Generic_Results_1.xlsx"), value_type = "area", sheet = "Sheet1")
  expect_equal(d[1,"feature_area", drop = TRUE], 71)
  expect_equal(d[1,"feature_name", drop = TRUE], "S1P d16:1 [M>113]")
  expect_equal(d[1,"analysis_id", drop = TRUE], "006_EBLK_Extracted Blank+ISTD01")
})

 testthat::test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
   expect_error(read_analysisresult_table(testthat::test_path("001_Generic_Results_1.xlsx"), value_type = "area"),regexp = "Please define sheet name")
 })

 testthat::test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
   expect_error(read_analysisresult_table(testthat::test_path("001_Generic_Results_1.txt"), value_type = "area"),regexp = "Invalid file format")
 })

