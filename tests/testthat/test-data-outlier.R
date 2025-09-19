mexp_empty <- MidarExperiment()
mexp <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))
mexp_proc <- mexp
mexp_proc <- normalize_by_istd(mexp_proc)
mexp_proc <- quantify_by_istd(mexp_proc)


test_that("detect_outlier_pca works", {
  outliers <- detect_outlier_pca(mexp_proc, variable = "intensity", filter_data = FALSE,qc_types = "BQC",
                             outlier_detection = "mad", pca_component = 1,
                             fence_multiplicator = 2, log_transform = TRUE)
  expect_equal(outliers, c("194_BQC_PQC_B 06", "195_BQC_PQC_B 07"))

  outliers <- detect_outlier_pca(mexp_proc, variable = "intensity", filter_data = FALSE,
                             outlier_detection = "mad", pca_component = 1,
                             fence_multiplicator = 2, log_transform = TRUE)
  expect_contains(outliers, c("194_BQC_PQC_B 06", "195_BQC_PQC_B 07", "024_SPL_S005", "033_SPL_S013", "198_LTR_LTR04"))

  outliers <- detect_outlier_pca(mexp_proc, variable = "intensity", filter_data = FALSE,
                             outlier_detection = "mad", pca_component = 1,qc_types = "BQC",
                             log_transform = TRUE,
                             fence_multiplicator = 1.3)
  expect_equal(outliers, c("018_BQC_PQC01", "194_BQC_PQC_B 06", "195_BQC_PQC_B 07"))

  outliers <- detect_outlier_pca(mexp_proc, variable = "intensity", filter_data = FALSE,
                             outlier_detection = "mad", pca_component = 1,qc_types = "BQC",
                             log_transform = FALSE,
                             fence_multiplicator = 1.3)
  expect_equal(outliers, c("018_BQC_PQC01", "194_BQC_PQC_B 06", "195_BQC_PQC_B 07"))


  outliers <- detect_outlier_pca(mexp_proc, variable = "intensity", filter_data = FALSE,
                             outlier_detection = "sd", pca_component = 1,qc_types = "BQC",
                             log_transform = FALSE,
                             fence_multiplicator = 1.3)
  expect_equal(outliers, c("194_BQC_PQC_B 06", "195_BQC_PQC_B 07"))

  outliers <- detect_outlier_pca(mexp_proc, variable = "conc", filter_data = FALSE,
                             outlier_detection = "mad", pca_component = 1,qc_types = "BQC",
                             log_transform = TRUE,
                             fence_multiplicator = 1.3)
  expect_equal(outliers, c("194_BQC_PQC_B 06", "195_BQC_PQC_B 07"))

  outliers <- detect_outlier_pca(mexp_proc, variable = "intensity", filter_data = FALSE,
                             outlier_detection = "mad", pca_component = 2,qc_types = "BQC",
                             log_transform = TRUE,
                             fence_multiplicator = 1.3)
  expect_equal(outliers, NULL)
  expect_error(
    outliers <- detect_outlier_pca(mexp_proc, variable = "intensity", filter_data = TRUE,
                               outlier_detection = "mad", pca_component = 2,
                               log_transform = TRUE,
                               fence_multiplicator = 1.3),
  "Data has not been qc filtered")

  expect_error(
    outliers <- detect_outlier_pca(mexp_proc, variable = "intensity", filter_data = FALSE,
                               outlier_detection = "mad", pca_component = 2, summarize_fun = "rma",
                               log_transform = TRUE,
                               fence_multiplicator = 1.3),
    "Relative Mean Abundance")

    expect_error(
      outliers <- detect_outlier_pca(mexp_proc, variable = "intensity", filter_data = FALSE,
                             outlier_detection = "mad", pca_component = 1, qc_types =  c("BQC", "XYZ"),
                             fence_multiplicator = 2, log_transform = TRUE),
    "Following specified QC types are missing in the dataset: XYZ")
})
