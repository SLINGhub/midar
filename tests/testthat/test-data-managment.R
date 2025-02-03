
# get RDS for tests

# mexp <- midar::MidarExperiment()
# mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.csv"), import_metadata = FALSE)
# mexp <- midar::import_metadata_midarxlm(mexp,
#                                         path = testthat::test_path("testdata/MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1.xlsm"),
#                                         excl_unmatched_analyses = FALSE)
# readr::write_rds(mexp, testthat::test_path("testdata/MHQuant_demo.rds"))

mexp_empty <- MidarExperiment()
mexp <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))
mexp_proc <- mexp
mexp_proc <- normalize_by_istd(mexp_proc)
mexp_proc <- quantify_by_istd(mexp_proc)
mexp_filt <- filter_features_qc(mexp_proc, max.cv.conc.bqc = 20, min.intensity.median.bqc = 1000)

mexp2 <- lipidomics_dataset
mexp2_filt <- filter_features_qc(mexp2, min.intensity.median.spl = 1E6)

test_that("check_data_present works", {
  testthat::expect_true(check_data_present(mexp))
  testthat::expect_false(check_data_present(mexp_empty))
})

test_that("check_dataset_present works", {
  testthat::expect_true(check_dataset_present(mexp))
  testthat::expect_false(check_dataset_present(mexp_empty))
})

test_that("get_dataset_subset returns filtered dataset for unfiltered data", {
  expect_equal(nrow(mexp2@dataset), 14471)
  result <- get_dataset_subset(
    data = mexp2,
    filter_data = FALSE,
    include_qualifier = TRUE,
    qc_types = NA,
    include_feature_filter = NA,
    exclude_feature_filter = NA
  )
  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 14471)
  expect_false(all(result$is_quantifier))
  expect_true(any(result$is_istd))

  result <- get_dataset_subset(
    data = mexp2_filt,
    filter_data = TRUE,
    include_qualifier = TRUE,
    qc_types = NA,
    include_feature_filter = NA,
    exclude_feature_filter = NA
  )
  expect_equal(nrow(result), 6487)

  result <- get_dataset_subset(
    data = mexp2,
    filter_data = FALSE,
    include_qualifier = FALSE,
    qc_types = NA,
    include_feature_filter = NA,
    exclude_feature_filter = NA
  )
  expect_equal(nrow(result), 13972)
  expect_true(all(result$is_quantifier))

  result <- get_dataset_subset(
    data = mexp2,
    filter_data = FALSE,
    include_qualifier = FALSE,
    include_istd = FALSE,
    qc_types = NA,
    include_feature_filter = NA,
    exclude_feature_filter = NA
  )
  expect_equal(nrow(result), 9481)
  expect_true(!any(result$is_istd))

  result <- get_dataset_subset(
    data = mexp2,
    filter_data = FALSE,
    include_qualifier = TRUE,
    include_istd = TRUE,
    qc_types = c("BQC", "SPL"),
    include_feature_filter = NA,
    exclude_feature_filter = NA
  )
  expect_equal(unique(result$qc_type),c("BQC","SPL"))

  result <- get_dataset_subset(
    data = mexp2,
    filter_data = FALSE,
    include_qualifier = TRUE,
    include_istd = TRUE,
    qc_types = "QC|SPL",
    include_feature_filter = NA,
    exclude_feature_filter = NA
  )
  expect_equal(unique(result$qc_type),c("RQC","TQC","BQC","SPL"))

  result <- get_dataset_subset(
    data = mexp2,
    filter_data = FALSE,
    include_qualifier = TRUE,
    include_istd = TRUE,
    qc_types = NA,
    include_feature_filter = "PC|PE|TG",
    exclude_feature_filter =  NA
  )
  expect_equal(nrow(result), 9481)

  result <- get_dataset_subset(
    data = mexp2,
    filter_data = FALSE,
    include_qualifier = TRUE,
    include_istd = TRUE,
    qc_types = NA,
    include_feature_filter = "PC|PE|TG",
    exclude_feature_filter =  "ISTD|SIM"
  )
  expect_equal(nrow(result),  5988)

  result <- get_dataset_subset(
    data = mexp2,
    filter_data = FALSE,
    include_qualifier = TRUE,
    include_istd = TRUE,
    qc_types = NA,
    include_feature_filter = c("PC 40:6", "PC 40:8"),
    exclude_feature_filter =  NA
  )
  expect_equal(nrow(result),  998)

  result <- get_dataset_subset(
    data = mexp2,
    filter_data = FALSE,
    include_qualifier = TRUE,
    include_istd = TRUE,
    qc_types = NA,
    include_feature_filter =  NA,
    exclude_feature_filter = c("PC 40:6", "PC 40:8")
  )
  expect_equal(nrow(result),  13473)
})

test_that("get_dataset_subset handles errors in data or filters", {
  expect_error(
    result <- get_dataset_subset(
      data = mexp2,
      filter_data = TRUE),
    "Data has not been QC-filtered"
  )

  expect_error(
    result <- get_dataset_subset(
      data = mexp2,
      filter_data = FALSE,
      qc_types = c("CAL", "SPL")),
    "One or more specified `qc_types`"
  )

  expect_error(
    result <- get_dataset_subset(
      data = mexp2,
      filter_data = FALSE,
      qc_types = c("CAL|QA")),
    "`qc_type` filter criteria resulted in no "
  )

  expect_error(
    result <- get_dataset_subset(
      data = mexp2,
      filter_data = FALSE,
      include_istd = FALSE,
      qc_types = NA,
      include_feature_filter = "ISTD"),
    "The defined feature filter criteria resulted in no"
  )

  expect_error(
    result <- get_dataset_subset(
      data = mexp2,
      filter_data = FALSE,
      include_istd = TRUE,
      qc_types = NA,
      include_feature_filter = "PC",
      exclude_feature_filter = "PC"),
    "contain overlapping features"
  )
})



test_that("get_analyticaldata returns correct table", {
  expect_equal(dim(get_analyticaldata(mexp, annotated = FALSE)), c(1040, 33))
  expect_equal(dim(get_analyticaldata(mexp, annotated = TRUE)), c(1040, 17))
})

test_that("get_analysis_count works", {
  testthat::expect_equal(get_analysis_count(mexp), 65)
  testthat::expect_equal(get_analysis_count(mexp_empty), 0)
  testthat::expect_equal(get_analysis_count(mexp, qc_types = "SPL"), 31)
  testthat::expect_equal(get_analysis_count(mexp, qc_types = "CAL"), 0)
  testthat::expect_equal(get_analysis_count(mexp_empty, qc_types = "SPL"), 0)
})

test_that("get_feature_count works", {
  testthat::expect_equal(get_feature_count(mexp), 16)
  testthat::expect_equal(get_feature_count(mexp_empty), 0)
  testthat::expect_equal(get_feature_count(mexp, is_istd = TRUE), 2)
  testthat::expect_equal(get_feature_count(mexp, is_istd = FALSE), 14)
  testthat::expect_equal(get_feature_count(mexp, is_quantifier = TRUE), 8)
  testthat::expect_equal(get_feature_count(mexp, is_quantifier = FALSE), 8)
})

test_that("get_featurelist works", {
  testthat::expect_equal(length(get_featurelist(mexp)), 16)
  testthat::expect_equal(get_featurelist(mexp)[2], "S1P d17:1 [M>60]")
  testthat::expect_equal(get_featurelist(mexp_empty), NULL)
})


test_that("get analysis timings works", {
  mexp_notimestamp <- mexp
  mexp_notimestamp@dataset$acquisition_time_stamp <- lubridate::NA_POSIXct_

  testthat::expect_equal(as.character(get_analyis_start(mexp)), "2018-04-19 18:37:00")
  testthat::expect_equal(get_analyis_start(mexp_notimestamp), lubridate::NA_POSIXct_)
  testthat::expect_equal(as.character(get_analyis_start(mexp_empty)), NA_character_)

  testthat::expect_equal(as.character(get_analyis_end(mexp, estimate_sequence_end = FALSE)), "2018-04-21 05:04:00")
  testthat::expect_equal(as.character(get_analyis_end(mexp, estimate_sequence_end = TRUE)), "2018-04-21 05:14:00")
  testthat::expect_equal(get_analyis_end(mexp_notimestamp, estimate_sequence_end = FALSE), lubridate::NA_POSIXct_)
  testthat::expect_equal(get_analyis_end(mexp_notimestamp, estimate_sequence_end = TRUE), lubridate::NA_POSIXct_)
  testthat::expect_equal(as.character(get_analyis_end(mexp_empty, estimate_sequence_end = FALSE)), NA_character_)
  testthat::expect_equal(as.character(get_analyis_end(mexp_empty, estimate_sequence_end = TRUE)), NA_character_)

  testthat::expect_equal(as.character(get_runtime_median(mexp)), "10M 0S")
  testthat::expect_equal(as.character(get_runtime_median(mexp_notimestamp)), NA_character_)
  testthat::expect_equal(as.character(get_runtime_median(mexp_empty)), NA_character_)

  testthat::expect_equal(as.character(get_analysis_duration(mexp, estimate_sequence_end = FALSE)), "1d 10H 27M 0S")
  testthat::expect_equal(as.character(get_analysis_duration(mexp, estimate_sequence_end = TRUE)), "1d 10H 37M 0S")
  testthat::expect_equal(as.character(get_analysis_duration(mexp_notimestamp, estimate_sequence_end = FALSE)), NA_character_)
  testthat::expect_equal(as.character(get_analysis_duration(mexp_notimestamp, estimate_sequence_end = TRUE)), NA_character_)
  testthat::expect_equal(as.character(get_analysis_duration(mexp_empty, estimate_sequence_end = FALSE)), NA_character_)
  testthat::expect_equal(as.character(get_analysis_duration(mexp_empty, estimate_sequence_end = TRUE)), NA_character_)

  testthat::expect_equal(get_analysis_breaks(mexp,break_duration_min = 30), 2)
  testthat::expect_equal(get_analysis_breaks(mexp_notimestamp,break_duration_min = 30), NA_integer_)
  testthat::expect_equal(get_analysis_breaks(mexp_empty,break_duration_min = 30), NA_integer_)

})

test_that("update_after_normalization clears norm_intensity and conc if not normalized", {

  mexp_temp_1 <- mexp_proc
  testthat::expect_true(any(c("feature_norm_intensity", "feature_conc") %in% names(mexp_temp_1@dataset)))
  expect_message(mexp_temp_1 <- update_after_normalization(mexp_temp_1, is_normalized = FALSE, with_message = TRUE),
                 "normalized intensities and concentrations have been invalidated")
  testthat::expect_false(any(c("feature_norm_intensity", "feature_conc") %in% names(mexp_temp_1@dataset)))
  testthat::expect_false(mexp_temp_1@is_istd_normalized)
  testthat::expect_false(mexp_temp_1@is_quantitated)

  mexp_temp_1 <- mexp_proc
  testthat::expect_true(any(c("feature_norm_intensity", "feature_conc") %in% names(mexp_temp_1@dataset)))
  expect_message(mexp_temp_1 <- update_after_normalization(mexp_temp_1, is_normalized = TRUE, with_message = TRUE),
                 NA)
  testthat::expect_true(any(c("feature_norm_intensity", "feature_conc") %in% names(mexp_temp_1@dataset)))
  testthat::expect_true(mexp_temp_1@is_istd_normalized)
  testthat::expect_true(mexp_temp_1@is_quantitated)
})

test_that("update_after_quantitation clears norm_intensity and conc if not quantitated", {

  mexp_temp_1 <- mexp_proc
  testthat::expect_true(any(c("feature_norm_intensity", "feature_conc") %in% names(mexp_temp_1@dataset)))
  expect_message(mexp_temp_1 <- update_after_quantitation(mexp_temp_1, is_quantitated = FALSE, with_message = TRUE),
                 "Concentrations not valid anymore. Please reprocess data.")
  testthat::expect_false(any(c("feature_conc") %in% names(mexp_temp_1@dataset)))
  testthat::expect_true(any(c("feature_norm_intensity") %in% names(mexp_temp_1@dataset)))
  testthat::expect_true(mexp_temp_1@is_istd_normalized)
  testthat::expect_false(mexp_temp_1@is_quantitated)

  mexp_temp_1 <- mexp_proc
  testthat::expect_true(any(c("feature_norm_intensity", "feature_conc") %in% names(mexp_temp_1@dataset)))
  expect_message(mexp_temp_1 <- update_after_quantitation(mexp_temp_1, is_quantitated = TRUE, with_message = TRUE),
                 NA)
  testthat::expect_true(any(c("feature_conc") %in% names(mexp_temp_1@dataset)))
  testthat::expect_true(any(c("feature_norm_intensity") %in% names(mexp_temp_1@dataset)))
  testthat::expect_true(mexp_temp_1@is_istd_normalized)
  testthat::expect_true(mexp_temp_1@is_quantitated)
})

test_that("check_var_in_dataset returns correct error if column not present", {

  expect_error(check_var_in_dataset(mexp$dataset, "feature_conc"),
                 "Concentration data are not available, please process data")

  expect_error(check_var_in_dataset(mexp$dataset, "feature_norm_intensity"),
               "Normalized intensities not available, please process")

  expect_error(check_var_in_dataset(mexp$dataset, "feature_response"),
               "Response is not available, please choose")

  tbl <- mexp$dataset |> select(-"feature_area")
  expect_error(check_var_in_dataset(tbl, "feature_area"),
               "Peak area data are not available, please choose")
  expect_error(check_var_in_dataset(tbl, "feature_height"),
               "Peak height data are not available, please choose")

})

test_that("get_batch_boundaries returns correct values", {
  expect_equal(get_batch_boundaries(mexp, 1), c(1,33))
  expect_equal(get_batch_boundaries(mexp, 2), c(34,51))
  expect_equal(get_batch_boundaries(mexp, 3), c(52,65))
  expect_equal(get_batch_boundaries(mexp, c(1,3)), c(1,65))
  expect_error(get_batch_boundaries(mexp, c(1,2,3)),
               "Please provide a vector with one or two batch IDs")
})


test_that("set_analysis_order orders according to set criteria and absence/presence of timestamp", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("23_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq_notimestamp.csv"),
                                        import_metadata = FALSE)
  mexp <- midar::import_metadata_midarxlm(mexp,
                                          path = testthat::test_path("MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1_reorder.xlsm"),
                                          excl_unmatched_analyses = FALSE)

  #Orders according to the order in the input csv file (no timestamp available)
  expect_equal(mexp@dataset[[1, "analysis_id"]], "020_SPL_S001")

  # Error as no timestamp available
  expect_error(midar::set_analysis_order(mexp, order_by = "timestamp"), "Acquisition timestamps are not present")

  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("22_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq.csv"),
                                        import_metadata = FALSE)
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

# link_data_metadata
#    Tested in other parts

test_that("set_intensity_var returns correct messages", {
  mexp_proc_temp <- mexp_proc

  expect_message(mexp_proc_temp <- set_intensity_var(mexp_proc_temp, "feature_area"),
  "New feature intensity variable \\(`feature_area`\\)")


  expect_error(mexp_proc_temp <- set_intensity_var(mexp_proc_temp, "conc"),
                 "\\`feature_conc\\` is not present in the raw data")

  expect_error(mexp_proc_temp <- set_intensity_var(mexp_proc_temp, "feature_conc"),
               "\\`feature_conc\\` is not present in the raw data")

  mexp_proc_temp_withconc <- mexp_proc
  mexp_proc_temp_withconc@dataset_orig$feature_conc <- mexp_proc_temp_withconc@dataset$feature_conc
  expect_message(res <- set_intensity_var(mexp_proc_temp_withconc, "feature_conc"),
               "`conc` is not a typically used raw signal")
  expect_message(res <- set_intensity_var(mexp_proc_temp_withconc, "feature_conc"),
                 "New feature intensity variable \\(`feature_conc`\\) defined")

  expect_message(mexp_proc_temp <- set_intensity_var(mexp_proc_temp, "feature_area", auto_select = TRUE, ... =  c("feature_area", "feature_height")),
                 "Default feature intensity variable.*feature_area")

  mexp_proc_temp@dataset_orig <- mexp_proc_temp@dataset_orig |> select(-"feature_area")

  expect_message(mexp_proc_temp <- set_intensity_var(mexp_proc_temp, "feature_area", auto_select = TRUE, ... =  c("feature_area", "feature_height")),
                 "No typical feature intensity variable found in the data")

})

test_that("exclude_analyses excludes analyses", {

  mexp_temp <- mexp
  mexp_temp@annot_analyses[mexp_temp@annot_analyses$analysis_id == "030_SPL_S010",]$valid_analysis <- FALSE
  expect_message(mexp_temp_excl <-
                   exclude_analyses(mexp_temp,
                   analyses_to_exclude = c("026_SPL_S007", "020_SPL_S001"),
                   replace_existing = FALSE),
                 "A total of 3 analyses are now excluded for downstream processing")

  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "030_SPL_S010",]$valid_analysis)
  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "026_SPL_S007",]$valid_analysis)
  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "020_SPL_S001",]$valid_analysis)

  expect_message(mexp_temp_excl <-
                   exclude_analyses(mexp_temp,
                     analyses_to_exclude = c("026_SPL_S007", "020_SPL_S001"),
                     replace_existing = TRUE),
                 "2 analyses were excluded for downstream processing")
  expect_true(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "030_SPL_S010",]$valid_analysis)
  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "026_SPL_S007",]$valid_analysis)
  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "020_SPL_S001",]$valid_analysis)

  expect_error(mexp_temp_excl <-
                   exclude_analyses(
                     mexp_temp,
                     analyses_to_exclude = c("026_SPL_S007", "020_SPL_S010"),
                     replace_existing = FALSE),
                 "One or more provided `analysis_id` to exclude are not present")

  expect_error(mexp_temp_excl <-
                   exclude_analyses(
                     mexp_temp,
                     analyses_to_exclude = NA,
                     replace_existing = FALSE),
                 "No `analysis_id` provided. To \\(re\\)include")

  expect_message(mexp_temp_excl <-
                   exclude_analyses(
                     mexp_temp,
                     analyses_to_exclude = NA,
                     replace_existing = TRUE),
                 "All exclusions removed")

})



test_that("exclude_features excludes feautures", {

  mexp_temp <- mexp
  mexp_temp@annot_features[mexp_temp@annot_features$feature_id == "S1P d20:1 [M>60]",]$valid_feature <- FALSE
  expect_message(mexp_temp_excl <-
                   exclude_features(mexp_temp,
                                    features_to_exclude = c("S1P d19:1 [M>60]", "S1P d17:1 [M>60]"),
                                    replace_existing = FALSE),
                 "A total of 3 features are now excluded for downstream")

  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "S1P d20:1 [M>60]",]$valid_feature)
  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "S1P d19:1 [M>60]",]$valid_feature)
  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "S1P d17:1 [M>60]",]$valid_feature)

  expect_message(mexp_temp_excl <-
                   exclude_features(mexp_temp,
                                    features_to_exclude = c("S1P d19:1 [M>60]", "S1P d17:1 [M>60]"),
                                    replace_existing = TRUE),
                 "2 features were excluded for downstream processing")
  expect_true(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "S1P d20:1 [M>60]",]$valid_feature)
  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "S1P d19:1 [M>60]",]$valid_feature)
  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "S1P d17:1 [M>60]",]$valid_feature)

  expect_error(mexp_temp_excl <-
                   exclude_features(
                     mexp_temp,
                     features_to_exclude = c("1P d19:1 [M>60]", "1P d17:2 [M>60]"),
                     replace_existing = FALSE),
                 "One or more provided `feature_id` are not present")

  expect_error(mexp_temp_excl <-
                   exclude_features(
                     mexp_temp,
                     features_to_exclude = NA,
                     replace_existing = FALSE),
                 "No `feature_id` provided. To \\(re\\)include")

  expect_message(mexp_temp_excl <-
                   exclude_analyses(
                     mexp_temp,
                     analyses_to_exclude = NA,
                     replace_existing = TRUE),
                 "All exclusions removed")

})

