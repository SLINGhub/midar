mexp_empty <- MidarExperiment(title = "Test Experiment", analysis_type = "lipidomics")
mexp <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))
mexp_proc <- mexp
mexp_proc <- normalize_by_istd(mexp_proc)
mexp_proc <- quantify_by_istd(mexp_proc)

test_that("Construct MidarExperiments", {
  expect_type(MidarExperiment(), "S4")
  expect_equal(as.character(class(MidarExperiment())), "MidarExperiment")
  mexp <- MidarExperiment(title = "Test", analysis_type = "lipidomics")
  expect_equal(mexp@title, "Test")
  expect_equal(mexp@analysis_type, "lipidomics")
  expect_error(mexp <- MidarExperiment(title = "Test", analysis_type = "undefined"),
               "Invalid analysis type")
})

test_that(" MidarExperiments setter/getter work", {
 mexp <- MidarExperiment(title = "Test", analysis_type = "lipidomics")
 expect_equal(mexp@analysis_type, "lipidomics")

 expect_error(mexp <- MidarExperiment(title = "Test", analysis_type = "undefined"),
               "Invalid analysis type")

})

test_that("MidarExperiment $ accessor works correctly", {
  mexp <- MidarExperiment(title = "Test Experiment", analysis_type = "lipidomics")

  #  accessing valid slots
  expect_equal(mexp$title, "Test Experiment")
  expect_equal(mexp$analysis_type, "lipidomics")

  # invalid slot
  expect_error(mexp$invalid_slot, "is not valid for this object")

 })



test_that("show method displays correct title", {
  # Create a MidarExperiment object with a specific title
  mexp <- MidarExperiment(title = "Test Experiment", analysis_type = "lipidomics")
  print(mexp)
  # Capture the output of the show method
  expect_message(show(mexp), "Title: Test Experiment")
})

test_that("`show` method displays processing status", {

  # Assuming status_processing is a slot you can set - set it for testing

  text_output <- toString(cli::cli_fmt(print(mexp)))

  # Capture the output of the show method
  expect_match(text_output, "feature_area")
  expect_match(text_output, "Analyses manually excluded")

  mexp <- exclude_analyses(mexp, c("030_SPL_S010", "107_SPL_S075"), replace_existing = TRUE)
  mexp <- exclude_features(mexp, c("S1P d20:1 [M>60]", "S1P d20:1 [M>113]"), replace_existing = TRUE)
  mexp
})


test_that("check_integrity_analyses works", {

  expect_no_message(check_integrity_analyses(mexp_empty))
  expect_no_message(check_integrity_analyses(mexp))
  expect_no_message(check_integrity_analyses(mexp_proc))

  expect_true(check_integrity_analyses(mexp_empty))
  expect_true(check_integrity_analyses(mexp))
  expect_true(check_integrity_analyses(mexp_proc))

  mexp_temp <- mexp
  mexp_temp@annot_analyses <- mexp_temp@annot_analyses[c(-10, -2),]
  expect_error(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = FALSE, silent = FALSE),
               "No metadata present for 2 of 65 analyses\\: 007_SOLV_Blank01\\, 015_TQCd\\-100_TQC\\-100percent")

  expect_no_error(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = FALSE, silent = TRUE))
  expect_no_error(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = TRUE, silent = FALSE))

  mexp_temp@annot_analyses <- mexp_temp@annot_analyses[c(-10, -1, -2, -4),]
  expect_error(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = FALSE, silent = FALSE, max_num_print= 2),
               "6 of 65 analyses have no matching metadata.")

  expect_true(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = TRUE, silent = TRUE))
  expect_false(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = FALSE, silent = TRUE))

  mexp_temp <- mexp
  mexp_temp@dataset_orig <- mexp_temp@dataset_orig[mexp_temp@dataset_orig$analysis_id != "007_SOLV_Blank01" &
                                                     mexp_temp@dataset_orig$analysis_id != "015_TQCd-100_TQC-100percent",]

  expect_error(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = FALSE, silent = FALSE),
               "ollowing 2 analyses defined in the metadata are not present in the measurement data\\: 007_SOLV_Blank01\\, 015_TQCd\\-100_TQC-100percent")

  expect_error(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = FALSE, silent = FALSE, max_num_print= 1),
               "2 of 65 analyses defined in the metadata are not present in the measurement data")

  expect_false(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = FALSE, silent = TRUE))

  mexp_temp <- mexp
  mexp_temp@annot_analyses$analysis_id <- "none"
  expect_error(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = FALSE, silent = FALSE, max_num_print= 2),
               "None of the measurements/samples have matching metadata")

  expect_false(check_integrity_analyses(mexp_temp, excl_unmatched_analyses = FALSE, silent = TRUE))

})

test_that("check_integrity_analyses works", {
  expect_equal(get_status_flag(TRUE), "v")
  expect_equal(get_status_flag(FALSE), "x")
})

test_that("check_data works", {
  expect_no_error(check_data(mexp))
  expect_no_error(check_data(mexp_empty))
  expect_error(check_data(NULL), "`data` cannot be NULL")
  expect_error(check_data(tibble(a=1)), "`data` must be a MidarExperiment")
})
