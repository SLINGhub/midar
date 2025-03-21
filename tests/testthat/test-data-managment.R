
# get RDS for tests

# mexp <- midar::MidarExperiment()
# mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("4_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_DetailedMethods.csv"), import_metadata = FALSE)
# mexp <- midar::import_metadata_msorganiser(mexp,
#                                         path = testthat::test_path("testdata/MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1.xlsm"),
#                                         excl_unmatched_analyses = FALSE)
# readr::write_rds(mexp, testthat::test_path("testdata/MHQuant_demo.rds"))

mexp_empty <- MidarExperiment()
mexp <- lipidomics_dataset
mexp_proc <- mexp
mexp_proc <- normalize_by_istd(mexp_proc)
mexp_proc <- quantify_by_istd(mexp_proc)
mexp_filt <- filter_features_qc(mexp_proc, include_qualifier = FALSE, include_istd = FALSE, max.cv.conc.bqc = 20, min.intensity.median.bqc = 1000)

mexp2 <- lipidomics_dataset
mexp2_filt <- filter_features_qc(mexp2, include_qualifier = FALSE, include_istd = FALSE, min.intensity.median.spl = 1E6)

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
  expect_equal(dim(get_analyticaldata(mexp, annotated = FALSE)), c(14471, 20))
  expect_equal(dim(get_analyticaldata(mexp, annotated = TRUE)), c(14471, 19))
})

test_that("get_analysis_count works", {
  testthat::expect_equal(get_analysis_count(mexp), 499)
  testthat::expect_equal(get_analysis_count(mexp_empty), 0)
  testthat::expect_equal(get_analysis_count(mexp, qc_types = "SPL"), 374)
  testthat::expect_equal(get_analysis_count(mexp, qc_types = "CAL"), 0)
  testthat::expect_equal(get_analysis_count(mexp_empty, qc_types = "SPL"), 0)
})

test_that("get_feature_count works", {
  testthat::expect_equal(get_feature_count(mexp), 29)
  testthat::expect_equal(get_feature_count(mexp_empty), 0)
  testthat::expect_equal(get_feature_count(mexp, is_istd = TRUE), 9)
  testthat::expect_equal(get_feature_count(mexp, is_istd = FALSE), 20)
  testthat::expect_equal(get_feature_count(mexp, is_quantifier = TRUE), 28)
  testthat::expect_equal(get_feature_count(mexp, is_quantifier = FALSE), 1)
})

test_that("get_featurelist works", {
  testthat::expect_equal(length(get_featurelist(mexp)), 29)
  testthat::expect_equal(get_featurelist(mexp)[2], "CE 18:1 d7 (ISTD)")
  testthat::expect_equal(get_featurelist(mexp_empty), NULL)
})


test_that("get analysis timings works", {
  mexp_notimestamp <- mexp
  mexp_notimestamp@dataset$acquisition_time_stamp <- lubridate::NA_POSIXct_

  testthat::expect_equal(as.character(get_analyis_start(mexp)), "2017-10-20 14:15:36")
  testthat::expect_equal(get_analyis_start(mexp_notimestamp), lubridate::NA_POSIXct_)
  testthat::expect_equal(as.character(get_analyis_start(mexp_empty)), NA_character_)

  testthat::expect_equal(as.character(get_analyis_end(mexp, estimate_sequence_end = FALSE)), "2017-10-24 17:33:03")
  testthat::expect_equal(as.character(get_analyis_end(mexp, estimate_sequence_end = TRUE)), "2017-10-24 17:44:23")
  testthat::expect_equal(get_analyis_end(mexp_notimestamp, estimate_sequence_end = FALSE), lubridate::NA_POSIXct_)
  testthat::expect_equal(get_analyis_end(mexp_notimestamp, estimate_sequence_end = TRUE), lubridate::NA_POSIXct_)
  testthat::expect_equal(as.character(get_analyis_end(mexp_empty, estimate_sequence_end = FALSE)), NA_character_)
  testthat::expect_equal(as.character(get_analyis_end(mexp_empty, estimate_sequence_end = TRUE)), NA_character_)

  testthat::expect_equal(as.character(get_runtime_median(mexp)), "11M 20S")
  testthat::expect_equal(as.character(get_runtime_median(mexp_notimestamp)), NA_character_)
  testthat::expect_equal(as.character(get_runtime_median(mexp_empty)), NA_character_)

  testthat::expect_equal(as.character(get_analysis_duration(mexp, estimate_sequence_end = FALSE)), "4d 3H 17M 27S")
  testthat::expect_equal(as.character(get_analysis_duration(mexp, estimate_sequence_end = TRUE)), "4d 3H 28M 47S")
  testthat::expect_equal(as.character(get_analysis_duration(mexp_notimestamp, estimate_sequence_end = FALSE)), NA_character_)
  testthat::expect_equal(as.character(get_analysis_duration(mexp_notimestamp, estimate_sequence_end = TRUE)), NA_character_)
  testthat::expect_equal(as.character(get_analysis_duration(mexp_empty, estimate_sequence_end = FALSE)), NA_character_)
  testthat::expect_equal(as.character(get_analysis_duration(mexp_empty, estimate_sequence_end = TRUE)), NA_character_)

  testthat::expect_equal(get_analysis_breaks(mexp,break_duration_min = 30), 7)
  testthat::expect_equal(get_analysis_breaks(mexp_notimestamp,break_duration_min = 30), NA_integer_)
  testthat::expect_equal(get_analysis_breaks(mexp_empty,break_duration_min = 30), NA_integer_)

})

test_that("update_after_normalization clears norm_intensity and conc if not normalized", {

  mexp_temp_1 <- mexp_proc
  testthat::expect_true(any(c("feature_norm_intensity", "feature_conc") %in% names(mexp_temp_1@dataset)))
  expect_message(mexp_temp_1 <- update_after_normalization(mexp_temp_1, is_normalized = FALSE, with_message = TRUE),
                 "The normalized intensities and concentrations are no longer valid")
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
                 "Concentrations are no longer valid. Please reprocess the data.")
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


})

test_that("get_batch_boundaries returns correct values", {
  expect_equal(get_batch_boundaries(mexp, 1), c(1,93))
  expect_equal(get_batch_boundaries(mexp, 2), c(94,175))
  expect_equal(get_batch_boundaries(mexp, 3), c(176,258))
  expect_equal(get_batch_boundaries(mexp, c(1,3)), c(1,258))
  expect_error(get_batch_boundaries(mexp, c(1,2,3)),
               "Please provide a vector with one or two batch IDs")
})


test_that("set_analysis_order orders according to set criteria and absence/presence of timestamp", {
  mexp <- midar::MidarExperiment()
  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("23_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq_notimestamp.csv"),
                                        import_metadata = FALSE)
  mexp <- midar::import_metadata_msorganiser(mexp,
                                          path = testthat::test_path("MiDAR_Metadata_Template_191_20240226_MHQuant_S1P_V1_reorder.xlsm"),
                                          excl_unmatched_analyses = FALSE)

  #Orders according to the order in the input csv file (no timestamp available)
  expect_equal(mexp@dataset[[1, "analysis_id"]], "020_SPL_S001")

  # Error as no timestamp available
  expect_error(midar::set_analysis_order(mexp, order_by = "timestamp"), "Acquisition timestamps are not present")

  mexp <- midar::import_data_masshunter(mexp, path = testthat::test_path("22_Testdata_MHQuant_DefaultSampleInfo_RT-Areas-FWHM_notInSeq.csv"),
                                        import_metadata = FALSE)
  mexp <- midar::import_metadata_msorganiser(mexp,
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

  expect_message(mexp_proc_temp <- set_intensity_var(mexp_proc_temp, "feature_height", auto_select = TRUE, ... =  c("feature_area", "feature_height")),
                 "Default feature intensity variable.*feature_area")

  mexp_proc_temp@dataset_orig <- mexp_proc_temp@dataset_orig |> select(-"feature_area") |> select(-"feature_height")

  expect_message(mexp_proc_temp <- set_intensity_var(mexp_proc_temp, "feature_area", auto_select = TRUE, ... =  c("feature_area", "feature_height")),
                 "No typical feature intensity variable found in the data")

})

test_that("exclude_analyses excludes analyses", {

  mexp_temp <- mexp
  mexp_temp@annot_analyses[mexp_temp@annot_analyses$analysis_id == "Longit_TQC-80%",]$valid_analysis <- FALSE
  expect_message(mexp_temp_excl <-
                   exclude_analyses(mexp_temp,
                                    analyses = c("Longit_LTR 01", "Longit_TQC-100%"),
                   clear_existing = FALSE),
                 "A total of 3 analyses are now excluded for downstream processing")

  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "Longit_TQC-80%",]$valid_analysis)
  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "Longit_TQC-100%",]$valid_analysis)
  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "Longit_LTR 01",]$valid_analysis)

  expect_message(mexp_temp_excl <-
                   exclude_analyses(mexp_temp,
                                    analyses = c("Longit_batch6_16", "Longit_batch5_41"),
                     clear_existing = TRUE),
                 "2 analyses were excluded for downstream processing")
  expect_true(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "Longit_TQC-80%",]$valid_analysis)
  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "Longit_batch5_41",]$valid_analysis)
  expect_false(mexp_temp_excl@annot_analyses[mexp_temp_excl@annot_analyses$analysis_id == "Longit_batch5_41",]$valid_analysis)

  expect_error(mexp_temp_excl <-
                   exclude_analyses(
                     mexp_temp,
                     analyses = c("Longit_batch5_41", "020_SPL_S010"),
                     clear_existing = FALSE),
                 "One or more provided `analysis_id` to exclude are not present")

  expect_error(mexp_temp_excl <-
                   exclude_analyses(
                     mexp_temp,
                     analyses = NA,
                     clear_existing = FALSE),
                 "No `analysis_id` provided. To \\(re\\)include")

  expect_message(mexp_temp_excl <-
                   exclude_analyses(
                     mexp_temp,
                     analyses = NA,
                     clear_existing = TRUE),
                 "All exclusions removed")

})



test_that("exclude_features excludes features", {

  mexp_temp <- mexp
  mexp_temp@annot_features[mexp_temp@annot_features$feature_id == "PC 40:8",]$valid_feature <- FALSE
  expect_message(mexp_temp_excl <-
                   exclude_features(mexp_temp,
                                    features = c("PC 40:6", "PC 32:1"),
                                    clear_existing = FALSE),
                 "A total of 3 features are now excluded for downstream")

  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "PC 40:8",]$valid_feature)
  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "PC 32:1",]$valid_feature)
  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "PC 40:6",]$valid_feature)

  expect_message(mexp_temp_excl <-
                   exclude_features(mexp_temp,
                                    features = c("PC 40:8", "PC 40:6"),
                                    clear_existing = TRUE),
                 "2 features were excluded for downstream processing")
  expect_true(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "PC 32:1",]$valid_feature)
  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "PC 40:6",]$valid_feature)
  expect_false(mexp_temp_excl@annot_features[mexp_temp_excl@annot_features$feature_id == "PC 40:8",]$valid_feature)

  expect_error(mexp_temp_excl <-
                   exclude_features(
                     mexp_temp,
                     features = c("1P d19:1 [M>60]", "1P d17:2 [M>60]"),
                     clear_existing = FALSE),
                 "One or more provided `feature_id` are not present")

  expect_error(mexp_temp_excl <-
                   exclude_features(
                     mexp_temp,
                     features = NA,
                     clear_existing = FALSE),
                 "No `feature_id` provided. To \\(re\\)include")

  expect_message(mexp_temp_excl <-
                   exclude_analyses(
                     mexp_temp,
                     analyses = NA,
                     clear_existing = TRUE),
                 "All exclusions removed")

})

