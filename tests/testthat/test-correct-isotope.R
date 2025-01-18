mexp_orig <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))

mexp <- mexp_orig

test_that("correct_interferences corrects overlapping interferences", {

  # d18:2 is interfering with d18:1, which in turn is interfering with d18:0
  # the code does not correct for M+4 isotope interference

  expect_equal(mexp@dataset |>
                 filter(feature_id == "S1P d18:2 [M>60]",
                        analysis_id == "008_LTR_LTR01") |>
                 pull(feature_intensity),
               31526)

  expect_equal(mexp@dataset |>
                 filter(feature_id == "S1P d18:1 [M>60]",
                        analysis_id == "008_LTR_LTR01") |>
                    pull(feature_intensity),
               85299)

  expect_equal(mexp@dataset |>
                 filter(feature_id == "S1P d18:0 [M>60]",
                        analysis_id == "008_LTR_LTR01") |>
                 pull(feature_intensity),
               9919)

  expect_message(mexp <- correct_interferences(mexp, variable = "feature_intensity"),
                 "Interference-correction has been applied to 4 of the 16 features")

  expect_equal(mexp@dataset |>
                 filter(feature_id == "S1P d18:2 [M>60]",
                        analysis_id == "008_LTR_LTR01") |>
                 pull(feature_intensity),
               31526)

  expect_equal(mexp@dataset |>
                 filter(feature_id == "S1P d18:1 [M>60]",
                        analysis_id == "008_LTR_LTR01") |>
                 pull(feature_intensity),
               84304.7172490)

  expect_equal(mexp@dataset |>
                 filter(feature_id == "S1P d18:0 [M>60]",
                        analysis_id == "008_LTR_LTR01") |>
                 pull(feature_intensity),
               7256.0500353124)

})





test_that("correct_interferences corrects overlapping interferences", {

  # d18:2 is interfering with d18:1, which in turn is interfering with d18:0
  # the code does not correct for M+4 isotope interference


  mexp2 <- mexp_orig


  expect_equal(mexp2@dataset |>
                 filter(feature_id == "S1P d18:1 [M>60]",
                        analysis_id == "008_LTR_LTR01") |>
                 pull(feature_intensity),
               85299)


  mexp2 <- correct_interference_manually(mexp2,
                                         variable = "feature_intensity",
                                         feature = "S1P d18:1 [M>60]",
                                         interfering_feature = "S1P d18:2 [M>60]",
                                         relative_contribution = 0.0315385)



  expect_equal(mexp2@dataset |>
                 filter(feature_id == "S1P d18:1 [M>60]",
                        analysis_id == "008_LTR_LTR01") |>
                 pull(feature_intensity),
               84304.7172490)

  # test renaming of feature and additional correction
  mexp2 <- correct_interference_manually(mexp2,
                                         variable = "feature_intensity",
                                         feature = "S1P d18:0 [M>60]",
                                         interfering_feature = "S1P d18:1 [M>60]",
                                         relative_contribution = 0.0315872, updated_feature_id = "S1P d18:0 [M>60] corrected")

  expect_equal(mexp2@dataset |>
                 filter(feature_id == "S1P d18:0 [M>60] corrected",
                        analysis_id == "008_LTR_LTR01") |>
                 pull(feature_intensity),
               7256.0500353124)


  expect_error(  mexp2 <- correct_interference_manually(mexp2,
                                                        variable = "feature_intensity",
                                                        feature = "S1P d18:0 [M>601]",
                                                        interfering_feature = "S1P d18:1 [M>60]",
                                                        relative_contribution = 0.0315872, updated_feature_id = "S1P d18:0 [M>60] corrected"),
                 "Feature is not present in the dataset, please adjust `feature` argument"
  )
})
