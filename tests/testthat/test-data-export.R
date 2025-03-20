# library(fs)
# library(ggplot2)
# library(openxlsx2)
# library(readr)
# library(testthat)
set.seed(123)

mexp_orig <- lipidomics_dataset

mexp <- normalize_by_istd(mexp_orig)
mexp <- quantify_by_istd(mexp)
mexp_empty <- MidarExperiment()
mexp_filt <- filter_features_qc(mexp, include_qualifier = FALSE, include_istd = FALSE, max.cv.conc.bqc = 20)
mexp@title <-"Test Experiment"
mexp <- calc_qc_metrics(mexp)  # Ensure calc_qc_metrics is executed before
mexp_drift <- correct_drift_gaussiankernel(
  mexp_orig,
  variable = "intensity",
  ref_qc_types = "SPL",
  ignore_istd = FALSE)




test_that("save_report_xlsx creates an Excel file", {
  temp_file <- tempfile(fileext = ".xlsx")

  expect_message(
    save_report_xlsx(data = mexp, path = temp_file),
    "The data processing report of experiment 'Test Experiment' has been saved to"
  )

  expect_true(file.exists(temp_file))
  on.exit(unlink(temp_file)) # Clean up
})

test_that("save_report_xlsx adds .xlsx extension if not provided", {
  temp_file <- tempfile()
  mexp@title <- ""
  expect_message(
    save_report_xlsx(data = mexp, path = temp_file),
    "The data processing report has been saved to"
  )

  expect_true(file.exists(paste0(temp_file, ".xlsx")))
  unlink(paste0(temp_file, ".xlsx")) # Clean up
})

test_that("save_report_xlsx creates the correct sheets", {
  temp_file <- tempfile(fileext = ".xlsx")

  save_report_xlsx(data = mexp, path = temp_file)

  # Load the workbook and check for sheets
  w_xlm <- openxlsx2::wb_load(temp_file)
  expected_sheets <- c(
    "Info", "Feature_QC_metrics", "Calibration_metrics", "QCfilt_StudySamples",
    "QCfilt_AllSamples", "Conc_FullDataset",
    "Raw_Intensity_FullDataset", "Norm_Intensity_FullDataset",
    "SampleMetadata", "FeatureMetadata", "InternalStandards", "BatchInfo"
  )

  expect_setequal(w_xlm$sheet_names, expected_sheets)
  on.exit(unlink(temp_file)) # Clean up
})

test_that("save_report_xlsx inluded feature-filtered data", {
  temp_file <- tempfile(fileext = ".xlsx")

  save_report_xlsx(data = mexp_filt, path = temp_file, filtered_variable = "norm_intensity")


  # Load the workbook and check for sheets
  w_xlm <- openxlsx2::wb_load(temp_file)
  expected_sheets <- c(
    "Info", "Feature_QC_metrics", "Calibration_metrics", "QCfilt_NormInt_StudySamples",
    "QCfilt_NormInt_AllSamples", "Conc_FullDataset",
    "Raw_Intensity_FullDataset", "Norm_Intensity_FullDataset",
    "SampleMetadata", "FeatureMetadata", "InternalStandards", "BatchInfo"
  )
  expect_setequal(w_xlm$sheet_names, expected_sheets)
  tbl <- openxlsx2::wb_to_df(temp_file, sheet = "QCfilt_NormInt_StudySamples")
  expect_equal(dim(tbl), c(374, 19))
  tbl <- openxlsx2::wb_to_df(temp_file, sheet = "QCfilt_NormInt_AllSamples")
  expect_equal(dim(tbl), c(476, 20))
  on.exit(unlink(temp_file)) # Clean up
})

test_that("save_report_xlsx handles missing data", {
  # Create a dataset with missing fields
  mexp_empty <- MidarExperiment()

  temp_file <- tempfile(fileext = ".xlsx")

  save_report_xlsx(data = mexp_empty, path = temp_file)

  # Verify the file creation
  expect_true(file.exists(temp_file))

  # Optionally, verify the content
  wb <- openxlsx2::wb_to_df(temp_file, sheet = "Raw_Intensity_FullDataset")
  expect_true(all(is.na(wb[["No annotated raw data available."]])))

  on.exit(unlink(temp_file)) # Clean up
})







test_that("Function exports correct variables", {
  temp_file <- tempfile(fileext = ".csv")
  expect_message(
    save_dataset_csv(data = mexp, path = temp_file,
                     variable = "intensity", filter_data = FALSE),
    "Intensity values for 499 analyses and 29 features"
  )

  exported_data <- readr::read_csv(temp_file)
  expect_true("analysis_id" %in% colnames(exported_data))
  expect_true("SM 36:2 d9 (ISTD)" %in% colnames(exported_data))
  expect_false("feature_intensity" %in% colnames(exported_data))
  expect_equal(mean(exported_data$`PC 40:6`), 3293741.4)

  temp_file <- tempfile(fileext = ".csv")
  expect_message(
    save_dataset_csv(data = mexp, path = temp_file,
                     variable = "conc", filter_data = FALSE),
    "Concentration values for 499 analyses and 19 features"
  )

  exported_data <- readr::read_csv(temp_file)
  expect_true("analysis_id" %in% colnames(exported_data))
  expect_false("SM 36:2 d9 (ISTD)" %in% colnames(exported_data))
  expect_false("feature_conc" %in% colnames(exported_data))
  expect_equal(mean(exported_data$`PC 40:6`), 0.082982104)
})

test_that("QC-filtered data is used when filter_data is TRUE", {
  temp_file <- tempfile(fileext = ".csv")
  expect_message(
    save_dataset_csv(data = mexp_filt, path = temp_file, variable = "area", filter_data = TRUE),
    "Area values for 499 analyses and 18 features")

  exported_data <- readr::read_csv(temp_file)
  expect_equal(dim(exported_data), c(499, 19))
})

test_that("QC-filtered data is used when filter_data is TRUE", {
  temp_file <- tempfile(fileext = ".csv")
  expect_message(
    save_dataset_csv(data = mexp_filt,
                     path = temp_file,
                     variable = "area",
                     filter_data = TRUE,
                     add_qctype = TRUE),
    "Area values for 499 analyses and 18 features")

  exported_data <- readr::read_csv(temp_file)
  expect_equal(dim(exported_data), c(499, 20))
  expect_true("qc_type" %in% colnames(exported_data))
})

test_that("Function handles non-existent variable gracefully", {
  expect_error(save_dataset_csv(data = mexp, path = tempfile(), variable = "non_existent_var", filter_data = FALSE),
               regexp = "`variable` must be one of")
})



test_that("QC types filtering works correctly", {
  temp_file <- tempfile(fileext = ".csv")
  save_dataset_csv(data = mexp,
                   path = temp_file,
                   variable = "intensity",
                   filter_data = FALSE,
                   qc_types = c("SPL", "BQC"),
                   add_qctype = TRUE)

  exported_data <- readr::read_csv(temp_file)
  expect_equal(unique(exported_data$qc_type), c("BQC", "SPL"))
})


# Test when the data is NULL
test_that("save_feature_qc_metrics handles NULL data input", {
  expect_error(save_feature_qc_metrics(mexp_empty, "output.csv"),
               "Feature QC metrics has not yet been")
})

# Test when the QC metrics are present
test_that("save_feature_qc_metrics exports QC metrics to CSV", {

  # Use a temporary file path to ensure tests do not interfere with actual files
  temp_file <- tempfile(fileext = ".csv")
  on.exit(unlink(temp_file))  # Ensure the temporary file is deleted after the test

  # Run the function and check it doesn't return errors
  expect_message(
    dat <- save_feature_qc_metrics(mexp_filt, temp_file),
                 "Feature QC metrics table was saved")

  expect_equal(dim(dat),c(29, 92))
  # Check if the file was created
  expect_true(file.exists(temp_file))

  # Read the file and compare with original data
  written_data <- readr::read_csv(temp_file)
  expect_equal(dim(written_data),c(29, 92))
})


test_that("save_metadata_templates() copies the file correctly and sends correct errors if required", {
  temp_file <- tempfile(fileext = "test.xlsx")

  expect_message(save_metadata_templates(temp_file), "Metadata table templates were saved to")
  expect_true(file.exists(temp_file))
  expect_error(save_metadata_templates(temp_file), "A file with this name already exists at the specified location.")
  unlink(temp_file)

  default_file <- "metadata_template.xlsx"
  if (file.exists(default_file)) unlink(default_file)

  expect_message(save_metadata_templates(), "Metadata table templates were saved to 'metadata_template.xlsx'")
  expect_true(file.exists(default_file))
  unlink(default_file)
})

# test_that("save_metadata_templates() returns an error if template is missing", {
#   # Temporarily change the system.file() return to an empty string
#   mock_template_path <- function(...) { "" }
#
#   with_mocked_bindings(
#     `system.file` = mock_template_path,
#     .package = "midar",
#     expect_error(save_metadata_templates(tempfile()), "Template file not found in package")
#   )
# })

test_that("save_metadata_msorganizer_template() copies the file correctly and sends correct errors if required", {
  temp_file <- tempfile(fileext = ".xlsx")

  expect_message(save_metadata_msorganizer_template(temp_file), "A MiDAR Metadata Organizer template was saved")
  expect_true(file.exists(temp_file))
  expect_error(save_metadata_msorganizer_template(temp_file), "A file with this name already exists at the specified location.")
  unlink(temp_file)

  default_file <- "metadata_msorganizer_template.xlsm"
  if (file.exists(default_file)) unlink(default_file)

  expect_message(save_metadata_msorganizer_template(), "A MiDAR Metadata Organizer template was saved to 'metadata_msorganizer_template.xlsm'")
  expect_true(file.exists(default_file))
  unlink(default_file)
})
#
# test_that("save_metadata_templates() returns an error if template is missing", {
#   # Temporarily change the system.file() return to an empty string
#   mock_template_path <- function(...) { "" }
#
#   with_mocked_bindings(
#     `system.file` = mock_template_path,
#     .package = "midar",
#     expect_error(save_metadata_msorganizer_template(tempfile()), "Template file not found in package")
#   )
# })
