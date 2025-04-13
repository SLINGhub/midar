# library(testthat)

test_that("Parses basic Agilent MH-Quant .csv file with only peak areas", {
  d <- parse_masshunter_csv(test_path("1_Testdata_MHQuant_DefaultSampleInfo_AreaOnly.csv"))

  expect_contains(names(d), c("file_analysis_order", "raw_data_filename", "sample_name", "sample_type", "acquisition_time_stamp", "feature_id", "feature_area"))
  expect_equal(nrow(d), 1040)
  expect_equal(mean(d$feature_area, na.rm = TRUE), 17237.244)
  expect_contains(unname(unlist(lapply(d, \(x) class(x)[[1]]))), c("integer", "character", "character", "character", "POSIXct", "character", "numeric"))
})

test_that("Parses basic Agilent MH-Quant .csv file without sample_name", {
  d <- parse_masshunter_csv(test_path("1_Testdata_MHQuant_DefaultSampleInfo_AreaOnly.csv"))

  expect_contains(names(d), c("file_analysis_order", "raw_data_filename", "sample_name", "sample_type", "acquisition_time_stamp", "feature_id", "feature_area"))
  expect_equal(nrow(d), 1040)
  expect_equal(mean(d$feature_area, na.rm = TRUE), 17237.244)
  expect_contains(unname(unlist(lapply(d, \(x) class(x)[[1]]))), c("integer", "character", "character", "character", "POSIXct", "character", "numeric"))
})

test_that("Parses nested Agilent MH-Quant .csv file with diverse peak variables", {
  d <- parse_masshunter_csv(test_path("3_Testdata_MHQuant_DefaultSampleInfo_DetailedResults.csv"))

  expect_identical(names(d), c(
    "analysis_id", "file_analysis_order", "raw_data_filename", "sample_name", "sample_type", "sample_level", "acquisition_time_stamp", "feature_id", "integration_qualifier", "feature_rt", "feature_area",
    "feature_fwhm", "feature_height", "feature_int_start", "feature_int_end", "feature_sn_ratio", "feature_symetry", "feature_width", "feature_manual_integration"
  ))
  expect_equal(nrow(d), 1040)
  expect_equal(mean(d$feature_area, na.rm = TRUE), 17237.244)
  expect_identical(unname(unlist(lapply(d, \(x) class(x)[[1]]))), c(c("character", "integer", "character", "character", "character", "character", "POSIXct", "character", "logical"), rep("numeric", 9), "logical"))
})

test_that("Parses nested Agilent MH-Quant .csv file with detailed sample info and different peak parameters", {
  d <- parse_masshunter_csv(test_path("5_Testdata_MHQuant_DetailedSampleInfo-RT-Areas-FWHM.csv"))

  expect_identical(names(d), c(
    "analysis_id", "file_analysis_order", "raw_data_filename", "sample_name", "sample_group", "sample_type", "sample_level", "acquisition_time_stamp", "inj_volume", "comment", "completed",
    "dilution_factor", "instrument_name", "instrument_type", "acq_method_file", "acq_method_path", "data_file_path", "feature_id", "integration_qualifier", "feature_rt", "feature_area", "feature_fwhm"
  ))
  expect_equal(nrow(d), 1040)
  expect_equal(mean(d$feature_area, na.rm = TRUE), 17237.244)
  expect_identical(unname(unlist(lapply(d, \(x) class(x)[[1]]))), c(
    "character", "integer", "character", "character", "character", "character", "character", "POSIXct", "numeric", "character", "character",
    "character", "character", "character", "character", "character", "character", "character", "logical", "numeric", "numeric", "numeric"
  ))
})

test_that("Parses nested MH Quant .csv file with detailed method info and different peak parameters", {
  d <- parse_masshunter_csv(test_path("4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.csv"))

  expect_identical(names(d), c(
    "analysis_id", "file_analysis_order", "raw_data_filename", "sample_name", "sample_type","sample_level", "acquisition_time_stamp", "feature_id", "integration_qualifier", "method_compound_group",
    "method_collision_energy", "method_fragmentor", "method_compound_id", "method_integration_method", "method_integration_parameters", "method_polarity", "method_ion_source",
    "method_multiplier", "method_noise_algorithm", "method_noise_raw_signal", "method_precursor_mz", "method_product_mz", "method_peak_smoothing", "method_peak_smoothing_gauss_width",
    "method_peak_smoothing_function_width", "method_transition", "method_time_segment", "method_type", "feature_rt", "feature_area", "feature_fwhm"
  ))
  expect_equal(nrow(d), 1040)
  expect_equal(mean(d$feature_area, na.rm = TRUE), 17237.244)
  expect_identical(unname(unlist(lapply(d, \(x) class(x)[[1]]))), c(
    "character","integer","character","character","character", "character", "POSIXct", "character", "logical", "character", "numeric", "numeric",
    "character", "character", "character", "factor", "character", "numeric", "character", "numeric", "numeric", "numeric", "character", "character", "character", "character",
    "integer", "character", "numeric", "numeric", "numeric"
  ))
})


test_that("Parses nested MH Quant .csv without the 'outlier' column", {
  d <- parse_masshunter_csv(test_path("6_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoOutlierSum.csv"))
  expect_equal(ncol(d), 12)
  expect_equal(nrow(d), 1040)
  expect_equal(names(d)[1], "analysis_id")
  expect_equal(d[[1, "feature_rt"]], 3.422)
  expect_equal(d[[1, "feature_area"]], 71)
})

test_that("Parses nested MH Quant .csv without the 'outlier' and 'quant message' column", {
  d <- parse_masshunter_csv(test_path("7_Testdata_MHQuant_NoOutlierSum-noQuantMsgSum.csv"))
  expect_equal(ncol(d), 12)
  expect_equal(nrow(d), 1040)
  expect_equal(names(d)[1], "analysis_id")
  expect_equal(d[[1, "feature_rt"]], 3.422)
  expect_equal(d[[1, "feature_area"]], 71)
})


test_that("Parsing nested MH Quant .csv without 'outlier'/'quant message' columns and header 'Samples' in first row/col", {

  d <- parse_masshunter_csv(test_path("8_Testdata_MHQuant_Corrupt_OutlierQuantMsgSumDeleted.csv"))
  expect_equal(ncol(d), 12)
  expect_equal(nrow(d), 1040)
  expect_equal(names(d)[1], "analysis_id")
  expect_equal(d[[1, "feature_rt"]], 3.422)
  expect_equal(d[[1, "feature_area"]], 71)
})


test_that("Parses nested MH Quant .csv file containing QUALIFIER peak info", {
  d <- parse_masshunter_csv(
    test_path("9_Testdata_MHQuant_withQuantMethods_withQualifierMethResults.csv"),
    expand_qualifier_names = TRUE
  )
  expect_equal(ncol(d), 18)
  expect_equal(nrow(d), 1040)
  expect_equal(d[[1, "feature_rt"]], 3.422)
  expect_equal(d[[1, "feature_area"]], 51)
  expect_equal(d |> filter(integration_qualifier) |> pull(feature_id) |> dplyr::first(), "S1P d16:1 [M>60] [QUAL 408.3 -> 113.0]")
  expect_equal(sum(d$integration_qualifier[d$file_analysis_order == 1]), 8)
})

test_that("Parses nested MH Quant .csv without Quant Message Summary", {
  d <- parse_masshunter_csv(test_path("10_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-NoQuantMsgSum.csv"))
  expect_equal(ncol(d), 12)
  expect_equal(nrow(d), 1040)
  expect_equal(names(d)[1], "analysis_id")
  expect_equal(d[[1, "feature_rt"]], 3.422)
  expect_equal(d[[1, "feature_area"]], 71)
})

test_that("Parses nested MH Quant .csv without acquistion time stamp", {
  d <- parse_masshunter_csv(test_path("11_Testdata_MHQuant_DefaultSampleInfo-noAcqDataTime_RT-Areas-FWHM.csv"))
  expect_false(c("acquisition_time_stamp") %in% names(d))
  expect_equal(ncol(d), 11)
  expect_equal(nrow(d), 1040)
  expect_equal(d[[1, "feature_rt"]], 3.422)
  expect_equal(d[[1, "feature_area"]], 71)
})

test_that("Returns a defined error when reading nested MH Quant .csv containg a 'Quantitation Message' is imported", {
  expect_error(
    parse_masshunter_csv(test_path("12_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM-withQuantMsg.csv")),
    regexp = "Field \\'Quantitation Message\\' currently not supported"
  )
})

test_that("Returns a defined error when reading MH Quant .csv with analytes/features as rows (Compound Table) is imported", {
  expect_error(
    parse_masshunter_csv(test_path("13_Testdata_MHQuant_CompoundTable_DefaultSampleInfo_RT-Areas-FWHM.csv")),
    regexp = "Compound table format is currently not supported", fixed = TRUE
  )
})

test_that("Returns a defined error when reading MH Quant .csv with analytes/features as rows (Compound Table) is imported", {
  expect_error(
    parse_masshunter_csv(test_path("24_Testdata_MHQuant_DefaultSampleInfo_noRawdatafilename.csv")),
    regexp = "'Data File' column is required and used as a unique identifier, but is missing or the file", fixed = TRUE
  )
})



test_that("Returns a defined error when reading a corrupted MH Quant .csv", {
  expect_error(
    parse_masshunter_csv(test_path("14_Testdata_MHQuant_Corrupt_RowAreaDeleted.csv")),
    regexp = "Data file is in an unsupported or corrupted format. Please try re-export your data in MH with compounds as columns",
    fixed = TRUE
  )
})

test_that("Parses nested MH Quant .csv file that has am empty first row", {
  d <- parse_masshunter_csv(test_path("15_Testdata_MHQuant_Corrupt_ExtraTopLine.csv"))
  expect_equal(ncol(d), 18)
  expect_equal(nrow(d), 1040)
  expect_equal(names(d)[1], "analysis_id")
  expect_equal(d[[1, "feature_rt"]], 3.422)
  expect_equal(d[[1, "feature_area"]], 51)
})


test_that("Parses nested MH Quant .csv file exported from German Windows system with comma as decimal point", {
  d <- parse_masshunter_csv(test_path("17_Testdata_Lipidomics_GermanSystem.csv"))
  expect_equal(d[[1, "feature_rt"]], 9.754)
  expect_equal(d[[1, "feature_fwhm"]], 0.056)
})

test_that("Parses nested MH Quant .csv file in UTF-8 format with different languages/characters", {
  d <- parse_masshunter_csv(test_path("18_Testdata_Lipidomics_MultiLanguageCharactersSamplenamesFeatures.csv"))
  expect_equal(d[[1, "feature_rt"]], 9.754)
  expect_equal(d[[1, "feature_id"]], "谷氨酰胺")
  expect_equal(d[[2, "feature_id"]], "글루타민")
  expect_equal(d[[3, "feature_id"]], "Glutaminsäure")
  expect_equal(d[[4, "feature_id"]], "குளுட்டமின்")
  expect_equal(d[[1, "raw_data_filename"]], "Über_Schöner_Blank")
  expect_equal(d[[300, "raw_data_filename"]], "空白的")
  expect_equal(d[[600, "raw_data_filename"]], "공백2")
  expect_equal(d[[900, "raw_data_filename"]], "空白")
  expect_equal(d[[1200, "raw_data_filename"]], "வெற்று")
})

test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
  d <- parse_masshunter_csv(test_path("19_Testdata_MHQuant_MultipleQUAL_with_expectedRT.csv"), expand_qualifier_names = TRUE)
  expect_equal(d[[1, "feature_rt"]], 6.649)
  expect_equal(d[[1, "method_target_rt"]], 7.200)
})

test_that("Parses nested MH Quant .csv file with special characters (e.g. !@#$%^) in feature names", {
  d <- parse_masshunter_csv(test_path("20_Testdata_MHQuant_withSpecialCharsInFeatures.csv"), expand_qualifier_names = TRUE)
  expect_equal(d[[1, "feature_rt"]], 6.649)
  expect_equal(d[[1, "method_target_rt"]], 7.200)
  expect_equal(d[[3, "feature_id"]], "Analyte 2* 14:0")
  expect_equal(d[[5, "feature_id"]], "Analyte 2~%^$# 14:0 (d5)")
  expect_match(d[[10, "feature_id"]], "Analyte 3 16\\:0-_=\\+\\\\/~!@ \\[QUAL 317\\.3 -> 299\\.3\\]")
  expect_match(d[[13, "feature_id"]], "Analyte 4 \\?><,\\.\`\\:\"\\}\\{\\[\\]18:0 \\[QUAL 331\\.3 -> 331\\.3\\]")
})


test_that("Parses nested MH Quant .csv with . \ | in feature names", {
  d <- parse_masshunter_csv(test_path("21_Testdata_MHQuant_with_dots_InFeatures.csv"), expand_qualifier_names = TRUE)
  expect_equal(d[[1, "feature_rt"]], 3.422)
  expect_match(d[[3, "feature_id"]], "S1P d17\\.1\\|S1P d17\\.2 \\[M>113\\]")
  expect_match(d[[4, "feature_id"]], "S1P d17\\.1\\\\S1P d17:2 \\[M>60\\]")
  expect_match(d[[5, "feature_id"]], "S1P d18\\.0/S1P 18:1 \\[M>113\\]")
})



test_that("Imports nested MH Quant .csv file containing QUALIFIER peak info into a MidarExperiment", {
  mexp <- MidarExperiment()
  mexp <- import_data_masshunter(
    mexp,
    test_path("9_Testdata_MHQuant_withQuantMethods_withQualifierMethResults.csv"),
    import_metadata = TRUE,
    expand_qualifier_names = TRUE
  )
  d <- mexp@dataset
  expect_equal(ncol(d), 18)
  expect_equal(nrow(d), 1040)
  expect_equal(d[[1, "feature_rt"]], 3.422)
  expect_equal(d[[1, "feature_area"]], 51)
  expect_equal(d[[2, "feature_id"]], "S1P d16:1 [M>60] [QUAL 408.3 -> 113.0]")
  expect_equal(nrow(mexp@annot_analyses), 65)
  expect_equal(nrow(mexp@annot_features), 16)
  expect_equal(nrow(mexp@annot_istds), 0)

})


test_that("Imports another MH Quant .csv files into one MidarExperiment", {
  mexp <- MidarExperiment()
  mexp <- import_data_masshunter(
    mexp,
    test_path("testdata/MHQuant_demo.csv"),
    import_metadata = TRUE,
    expand_qualifier_names = TRUE
  )
  d <- mexp@dataset
  expect_equal(ncol(d), 19)
  expect_equal(nrow(d), 1178)
  expect_equal(d[[1, "feature_rt"]], 7.160)
  expect_equal(d[[1, "feature_area"]], 5152996.0)
  expect_equal(d[[2, "feature_id"]], "CE 18:1 d7 (ISTD)")
  expect_equal(nrow(mexp@annot_analyses), 38)
  expect_equal(nrow(mexp@annot_features), 31)
  expect_equal(nrow(mexp@annot_istds), 0)

})

# Splitted above file into 2 files and import as folder
test_that("Imports multiple MH Quant .csv files into one MidarExperiment 1", {
  mexp <- MidarExperiment()
  mexp <- import_data_masshunter(
    mexp,
    test_path("testdata/MQquant_multiple/"),
    import_metadata = TRUE,
    expand_qualifier_names = TRUE
  )
  d <- mexp@dataset
  expect_equal(ncol(d), 19)
  expect_equal(nrow(d), 1178)
  expect_equal(d[[1, "feature_rt"]], 7.160)
  expect_equal(d[[1, "feature_area"]], 5152996.0)
  expect_equal(d[[2, "feature_id"]], "CE 18:1 d7 (ISTD)")
  expect_equal(nrow(mexp@annot_analyses), 38)
  expect_equal(nrow(mexp@annot_features), 31)
  expect_equal(nrow(mexp@annot_istds), 0)

})

# Splitted above file into 2 files with 1 overlapping (duplicated) feature and import as folder
test_that("Error duplicated reporting when import multiple MH Quant .csv files into one MidarExperiment", {
  mexp <- MidarExperiment()

  expect_error(
    mexp <- import_data_masshunter(
      mexp,
      test_path("testdata/MQquant_multiple_duplicates/"),
      import_metadata = TRUE,
      expand_qualifier_names = TRUE
    ),
    regexp = "duplicated reportings \\(analysis and feature pairs\\) with identical")
})

# Splitted above file into 2 files with 1 overlapping (duplicated) feature and import as folder
test_that("Error file not exist", {
  mexp <- MidarExperiment()

  expect_error(
    mexp <- import_data_masshunter(
      mexp,
      test_path(c("testdata/MQquant_multiple_duplicates/MHQuant_demo_Part1.csv",
                  "testdata/MQquant_multiple_duplicates/MHQuant_demo_Part3.csv")),
      import_metadata = TRUE,
      expand_qualifier_names = TRUE
    ),
    regexp = "One or more given files do not exist", fixed = TRUE)

  expect_error(
    mexp <- import_data_masshunter(
      mexp,
      test_path(c("testdata/MQquant_multiple_duplicates/MHQuant_demo_Part1.csv",
                  "testdata/MQquant_multiple_duplicates/MHQuant_demo_Part1.csv")),
      import_metadata = TRUE,
      expand_qualifier_names = TRUE
    ),
    regexp = "One or more given files are duplicated", fixed = TRUE)
})


# Splitted above file into 2 files with 1 overlapping (duplicated) feature with different values and import as folder
test_that("Imports multiple MH Quant .csv files into one MidarExperiment", {
  mexp <- MidarExperiment()

  expect_error(
    mexp <- import_data_masshunter(
      mexp,
      test_path("testdata/MQquant_multiple_duplicates2/"),
      import_metadata = TRUE,
      expand_qualifier_names = TRUE
    ),
    regexp = "duplicated reportings \\(analysis and feature pairs\\) with different"
  )

})


test_that("Imports MH with Calc. Conc or Final Conc. and Exp. Conc missing Name (sample name) and Sample header", {
  mexp <- MidarExperiment()

  expect_message(
    mexp <- import_data_masshunter(
      mexp,
      test_path("testdata/QuantLCMS_Example_MassHunter.csv"), expand_qualifier_names = TRUE),
    "Imported 25 analyses with 16 features (8 quantifiers, 8 qualifiers)", fixed = TRUE
  )

  expect_true("feature_conc_calc" %in% names(mexp@dataset_orig))
  expect_true(identical(mexp$dataset$feature_conc_final, mexp$dataset$feature_conc))

  expect_message(
    mexp <- import_data_masshunter(
      mexp,
      test_path("testdata/QuantLCMS_Example_MassHunter.csv"), expand_qualifier_names = TRUE, conc_column = "conc_calc"),
    "Imported 25 analyses with 16 features (8 quantifiers, 8 qualifiers)", fixed = TRUE
  )

  expect_true("feature_conc_calc" %in% names(mexp@dataset_orig))
  expect_true(identical(mexp$dataset$feature_conc_calc, mexp$dataset$feature_conc))

  expect_message(
    mexp <- import_data_masshunter(
      mexp,
      test_path("testdata/QuantLCMS_Example_MassHunter_FinalConc.csv"), expand_qualifier_names = TRUE, conc_column = "conc_calc"),
    "Imported 25 analyses with 16 features (8 quantifiers, 8 qualifiers)", fixed = TRUE
  )

  expect_false("feature_conc_calc" %in% names(mexp@dataset_orig))
  expect_true(identical(mexp$dataset$feature_conc_final, mexp$dataset$feature_conc))

  expect_message(
    mexp <- import_data_masshunter(
      mexp,
      test_path("testdata/QuantLCMS_Example_MassHunter_CalcConc.csv"), expand_qualifier_names = TRUE, conc_column = "conc_calc"),
    "Imported 25 analyses with 16 features (8 quantifiers, 8 qualifiers)", fixed = TRUE
  )

  expect_false("feature_conc_final" %in% names(mexp@dataset_orig))
  expect_true(identical(mexp$dataset$feature_conc_calc, mexp$dataset$feature_conc))


  expect_message(
  mexp <- import_data_masshunter(
      mexp,
      test_path("testdata/QuantLCMS_Example_MassHunter-NoHdrSampleName.csv"), expand_qualifier_names = TRUE),
  "Imported 25 analyses with 16 features (8 quantifiers, 8 qualifiers)", fixed = TRUE
  )

  expect_message(
    mexp <- import_data_masshunter(
      mexp,
      test_path("testdata/QuantLCMS_Example_MassHunter.csv"), expand_qualifier_names = TRUE),
    "Imported 25 analyses with 16 features (8 quantifiers, 8 qualifiers)", fixed = TRUE
  )
})

test_that("Imports MRMkit result file (long format) into a MidarExperiment", {
  mexp <- MidarExperiment()
  mexp <- import_data_mrmkit(
    mexp,
    test_path("testdata/MRMkit_demo.tsv"),
    import_metadata = TRUE,
  )
  d <- mexp@dataset
  expect_equal(ncol(d), 19)
  expect_equal(nrow(d), 13972.0)
  expect_equal(d[[1, "feature_rt"]], 7.2950)
  expect_equal(d[[1, "feature_area"]], 3134.16360)
  expect_equal(d[[2, "feature_id"]], "CE 18:1 d7 (ISTD)")
  expect_equal(nrow(mexp@annot_analyses), 499)
  expect_equal(nrow(mexp@annot_features), 28)
  expect_equal(nrow(mexp@annot_istds), 0)

})

test_that("Handles import_data_mrmkit errors", {
  mexp <- MidarExperiment()
  expect_error(
    mexp <- import_data_mrmkit(
      mexp,
      test_path("testdata/MRMkit_demo.txt"),
      import_metadata = TRUE,
    ),
    "Data file type/extension not supported", fixed = TRUE
    )

})

#' file_path = system.file("extdata", "MHQuant_demo.csv", package = "midar")
#'
#' mexp <- import_data_masshunter(
#'   data = mexp,
#'   path = file_path,
#'   import_metadata = TRUE,
#'   expand_qualifier_names = TRUE)


# test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
#   d <- read_data_table(test_path("001_Generic_Results_1.csv"), value_type = "area")
#   expect_equal(d[[1, "feature_area"]], 71)
#   expect_equal(d[[1, "feature_id"]], "S1P d16:1 [M>113]")
#   expect_equal(d[[1, "analysis_id"]], "006_EBLK_Extracted Blank+ISTD01")
# })
#
#
# test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
#   d <- read_data_table(test_path("001_Generic_Results_1.xlsx"), value_type = "area", sheet = "Sheet1")
#   expect_equal(d[[1, "feature_area"]], 71)
#   expect_equal(d[[1, "feature_id"]], "S1P d16:1 [M>113]")
#   expect_equal(d[[1, "analysis_id"]], "006_EBLK_Extracted Blank+ISTD01")
# })
#
# test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
#   expect_error(read_data_table(test_path("001_Generic_Results_1.xlsx"), value_type = "area"), regexp = "Please define sheet name")
# })
#
# test_that("Parses nested MH Quant .csv file with target (expected) RT and peak RT and multiple Qualifier per analyte", {
#   expect_error(read_data_table(test_path("001_Generic_Results_1.txt"), value_type = "area"), regexp = "Invalid file format")
# })


# Test parse_plain_csv

test_that("Parses plain csv file with metadata with correct column names and autodetecting analysis_id", {
  d <- parse_plain_csv(test_path("batch_effect-simdata-u1000-sd100_7batches.csv"), variable_name = "conc", import_metadata = TRUE)
  expect_identical(names(d), c("analysis_id","qc_type","batch_id","feature_id","feature_conc", "integration_qualifier"))
  expect_equal(mean(d$feature_conc[d$batch_id == 1 & d$feature_id == "Analyte-1"]), 1004.61572)

  d <- parse_plain_csv(test_path("testdata/plain_wide_nometadata.csv"), variable_name = "conc", import_metadata = TRUE)
  expect_identical(names(d), c("analysis_id","feature_id","feature_conc", "integration_qualifier"))
})




test_that("Returns error when parse_plain_csv imports other than csv", {
  expect_error(parse_plain_csv(test_path("testdata/MRMkit_demo.tsv"), variable_name = "conc", import_metadata = FALSE),
                 regexp = "Only csv files are currently supported", fixed = TRUE)
})

test_that("Returns error when plain csv file with columns containing text is read, when import_metadata = FALSE", {
  expect_message(parse_plain_csv(test_path("batch_effect-simdata-u1000-sd100_7batches.csv"), variable_name = "conc", import_metadata = FALSE),
                 regexp = "Metadata column(s) 'qc_type, batch_id' found and ignored", fixed = TRUE)
})

test_that("Returns error when plain csv file with analysis_id_col set that does not exist", {
  expect_error(parse_plain_csv(test_path("batch_effect-simdata-u1000-sd100_7batches.csv"), analysis_id_col = "sample_id", variable_name = "conc", import_metadata = TRUE),
               regexp = "No column with the name `sample_id` found in the data file.", fixed = TRUE)
})

test_that("Returns error when plain csv file with analysis_id_col  = NA and no analysis_id col present", {
  expect_error(parse_plain_csv(test_path("testdata/plain_wide_noanalysisid.csv"), variable_name = "conc", import_metadata = TRUE),
               regexp = "Column `analysis_id` not found in imported data.", fixed = TRUE)
})




test_that("Parses plain csv file with metadata and defined analysis_id_col, with correct data types", {
  d <- parse_plain_csv(test_path("batch_effect-simdata-diff_firstcol.csv"), analysis_id_col = "sample_id", variable_name = "conc", import_metadata = TRUE)
  expect_identical(names(d), c("analysis_id","qc_type","batch_id","feature_id","feature_conc", "integration_qualifier"))
  expect_identical(typeof(d$analysis_id), "character")
  expect_identical(typeof(d$batch_id), "character")
  expect_identical(typeof(d$feature_conc), "double")
})

test_that("Imports plain csv file with metadata parsing the numbers to 'analysis_id_col', with correct data types and metadata", {
  path <- test_path("testdata/plain_wide_dataset.csv")

  mexp <- MidarExperiment()
  expect_message(
    mexp <- import_data_csv(data = mexp,
                            path = path,
                            variable_name = "conc",
                            import_metadata = TRUE),
    "Imported 87 analyses with 5 features",
    fixed = TRUE)

  mexp <- MidarExperiment()
  expect_message(
    mexp <- import_data_csv(data = mexp,
                            path = path,
                            variable_name = "conc",
                            import_metadata = TRUE),
    "Metadata column(s) 'analysis_order, qc_type, batch_id' imported",
    fixed = TRUE)

  expect_in(c("analysis_id", "batch_id", "replicate_no", "is_istd", "feature_conc"), names(mexp@dataset))
  expect_equal(mexp@dataset[[111,"feature_conc"]], 892.82088)
  expect_equal(mexp@dataset[[50,"qc_type"]], "SPL")
  expect_equal(mexp@dataset[[255,"analysis_order"]], 51L)
  expect_equal(mexp@annot_analyses[[11,"analysis_id"]], "Spl11")
  expect_equal(mexp@annot_analyses[[11,"analysis_order"]], 77)
  expect_equal(mexp@annot_analyses[[10,"qc_type"]], "BQC")
  expect_equal(mexp@annot_features[[2,"feature_id"]], "S1P 18:2;O2")
  expect_type(mexp@dataset$is_istd, "logical")
  expect_type(mexp@dataset$batch_id, "character")

  mexp <- MidarExperiment()

  expect_message(
    mexp <- import_data_csv(data = mexp,
                            path = path,
                            variable_name = "conc",
                            import_metadata = FALSE),
    "Metadata column(s) 'analysis_order, qc_type, batch_id' found and ignored", fixed = TRUE
  )

  expect_in(c("analysis_id", "batch_id", "replicate_no", "is_istd", "feature_conc"), names(mexp@dataset))
  expect_equal(mexp@dataset[[111,"feature_conc"]], 897.39956)
  expect_equal(as.character(mexp@dataset[[50,"qc_type"]]), NA_character_)
  expect_equal(mexp@dataset[[255,"analysis_order"]], 51L)
  expect_equal(nrow(mexp@annot_analyses), 0)
  expect_equal(nrow(mexp@annot_features), 0)

  mexp <- MidarExperiment()
  expect_message(
    mexp <- import_data_csv(data = mexp,
                            path = path,
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = FALSE
    ),
    "Metadata column(s) 'analysis_order, qc_type, batch_id' found and ignored", fixed = TRUE
  )

  mexp <- MidarExperiment()
  expect_message(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset_no_order.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE
    ),
    "Analysis order was based on `analysis_order` column of imported data", fixed = TRUE
  )

  expect_equal(mexp@dataset[[111,"feature_conc"]], 897.39956)

  mexp <- MidarExperiment()
  expect_message(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE,
                            first_feature_column = "S1P 18:2;O2"
    ),
    "Analysis order was based on `analysis_order` column of imported data", fixed = TRUE
  )

  expect_in(c("S1P 18:1;O2"), names(mexp@dataset_orig))
  expect_false(c("S1P 18:1;O2") %in%names(mexp@dataset))
  expect_equal(mexp@dataset[[111,"feature_conc"]],  16.0568983)
  expect_equal(as.character(mexp@dataset[[50,"qc_type"]]), "SPL")
  expect_equal(mexp@dataset[[255,"analysis_order"]], 64L)

  mexp <- MidarExperiment()
  expect_error(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE,
                            first_feature_column = "NotThere"
    ),
    "Column `NotThere` not found in the data file", fixed = TRUE
  )

  mexp <- MidarExperiment()
  expect_error(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE,
                            first_feature_column = 10
    ),
    "Column index set via `first_feature_column` out of range", fixed = TRUE
  )

  mexp <- MidarExperiment()
  expect_message(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE,
                            first_feature_column = 2
    ),
    "Imported 87 analyses with 5 features", fixed = TRUE
  )
  mexp <- MidarExperiment()
  expect_error(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset_duplicate_analysisid.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE,
                            first_feature_column = 2
    ),
    "3 duplicated `analysis_id` present in the data file.", fixed = TRUE
  )
  mexp <- MidarExperiment()
  expect_error(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset_dup_featid.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE
    ),
    "1 duplicated column name(s) present in the data file", fixed = TRUE
  )
  mexp <- MidarExperiment()
  expect_error(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset_morecol.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE
    ),
    "ll columns with feature values must be numeric", fixed = TRUE
  )
  mexp <- MidarExperiment()
  expect_error(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset_morecol.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = FALSE
    ),
    "ll columns with feature values must be numeric", fixed = TRUE
  )
  mexp <- MidarExperiment()
  expect_message(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset_morecol.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = FALSE,
                            first_feature_column = "S1P 18:1;O2"),
    "Imported 87 analyses with 5 features", fixed = TRUE
  )
  mexp <- MidarExperiment()
  expect_error(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset_duplicate_orderid.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE
    ),
    "`analysis_order` contains duplicated values", fixed = TRUE
  )

  expect_error(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset_textorderid.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE
    ),
    "`analysis_order` contains duplicated values", fixed = TRUE
  )

  # Other dataset

  mexp <- MidarExperiment()
  expect_message(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset2_22rows.csv"),
                            variable_name = "area",
                            import_metadata = TRUE),
    "Metadata column(s) 'qc_type, batch_id' imported.",
    fixed = TRUE)

  expect_in(c("analysis_id", "batch_id", "replicate_no", "is_istd", "feature_area"), names(mexp@dataset))
  expect_equal(mexp@dataset[[11,"feature_area"]], 2276.88770)
  expect_equal(mexp@dataset[[50,"qc_type"]], "SPL")
  expect_equal(mexp@dataset[[13,"analysis_order"]], 2L)
  expect_equal(mexp@dataset[[13,"analysis_id"]], "P001-A02")
  expect_equal(mexp@annot_analyses[[9,"analysis_id"]], "P001-A09")
  expect_equal(mexp@annot_features[[2,"feature_id"]], "Cer d18:1/16:0 d7")
  expect_type(mexp@dataset$is_istd, "logical")
  expect_type(mexp@dataset$batch_id, "character")

  # CHeck order iD imported and batch id is string
  mexp2 <- MidarExperiment()
  expect_message(
    mexp2 <- import_data_csv(data = mexp,
                             path = test_path("testdata/plain_wide_dataset2_10rows_orderid.csv"),
                             variable_name = "area",
                             import_metadata = TRUE),
    "Metadata column(s) 'analysis_order, qc_type, batch_id' imported.",
    fixed = TRUE)

  expect_equal(mexp2@annot_analyses[[9,"analysis_id"]], "P001-A09")
  expect_equal(mexp2@dataset[[13,"analysis_order"]], 2L)
  expect_equal(mexp2@dataset[[13,"analysis_id"]], "P001-A09")
  expect_equal(mexp2@dataset[[13,"batch_id"]], "2") # must be text

  expect_message(
    mexp <-correct_drift_gaussiankernel(mexp, variable = "intensity", ref_qc_types = "SPL"),
    "-0.88% to -0.10%)", fixed = TRUE)

  expect_message(
    mexp <-correct_batch_centering(mexp, variable = "intensity", ref_qc_types = "SPL"),
    "-7.30% to 0.10%)", fixed = TRUE)

  p <- plot_runscatter(mexp, variable = "intensity", return_plot = TRUE)

 plot_runsequence(data = mexp2)

  plot_data <- ggplot2::ggplot_build(p[[1]])$data
  expect_equal(dim(plot_data[[2]]),c(176, 10))

  expect_error(
    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset2_10rows_orderidtext.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE,
                            first_feature_column = 10
    ),
    "Column `analysis_order` must contain unique numbers", fixed = TRUE
  )

    mexp <- import_data_csv(data = mexp,
                            path = test_path("testdata/plain_wide_dataset2_10rows_analysisidnumber.csv"),
                            variable_name = "conc",
                            analysis_id_col = "analysis_id",
                            import_metadata = TRUE,
                            first_feature_column = 10
    )

    expect_equal(mexp@dataset[[13,"analysis_id"]], "6") # must be text even if was number

})



