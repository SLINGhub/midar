mexp_orig <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))
mexp <- mexp_orig

mexp2 <- lipidomics_dataset

mexp2@annot_features$interference_contribution[9] <- 0.5

test_that("correct_interferences corrects overlapping interferences", {
  # d18:2 is interfering with d18:1, which in turn is interfering with d18:0
  # the code does not correct for M+4 isotope interference


  # Check initial uncorrected values

  expect_equal(
    mexp@dataset |>
      filter(
        feature_id == "S1P d18:2 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    31526
  )

  expect_equal(
    mexp@dataset |>
      filter(
        feature_id == "S1P d18:1 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    85299
  )

  expect_equal(
    mexp@dataset |>
      filter(
        feature_id == "S1P d18:0 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    9919
  )

  # Apply correction and check corrected values, with sequential correction enabled

  expect_message(
    mexp_res <-
      correct_interferences(
        mexp,
        variable = "feature_intensity",
        sequential_correction = TRUE
      ),
    "Interference-correction has been applied to 4 of the 16 features"
  )



  expect_equal(
    mexp_res@dataset |>
      filter(
        feature_id == "S1P d18:2 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    31526
  )

  expect_equal(
    mexp_res@dataset |>
      filter(
        feature_id == "S1P d18:1 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    84304.7172490
  )

  # corrrected with corrected S1P d18:1
  expect_equal(
    mexp_res@dataset |>
      filter(
        feature_id == "S1P d18:0 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    7256.0500353124
  )

  # reapply
  expect_message(
    mexp_res <-
      correct_interferences(
        mexp_res,
        variable = "feature_intensity",
        sequential_correction = TRUE
      ),
    "Interference-correction has been applied to 4 of the 16 features"
  )

  expect_equal(
    mexp_res@dataset |>
      filter(
        feature_id == "S1P d18:1 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    84304.7172490
  )

  # corrrected with corrected S1P d18:1
  expect_equal(
    mexp_res@dataset |>
      filter(
        feature_id == "S1P d18:0 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    7256.0500353124
  )

  # Apply correction and check corrected values, with sequential correction disabled

  expect_message(
    mexp_res <-
      correct_interferences(
        mexp,
        variable = "feature_intensity",
        sequential_correction = FALSE
      ),
    "Interference-correction has been applied to 4 of the 16 features"
  )

  expect_equal(
    mexp_res@dataset |>
      filter(
        feature_id == "S1P d18:2 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    31526
  )

  expect_equal(
    mexp_res@dataset |>
      filter(
        feature_id == "S1P d18:1 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    84304.7172490
  )

  # corrrected based on raw S1P d18:1 intensity
  expect_equal(
    mexp_res@dataset |>
      filter(
        feature_id == "S1P d18:0 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    7224.6434272
  )

})





test_that("correct_interferences corrects overlapping interferences", {
  # d18:2 is interfering with d18:1, which in turn is interfering with d18:0
  # the code does not correct for M+4 isotope interference


  mexp2 <- mexp_orig


  expect_equal(
    mexp2@dataset |>
      filter(
        feature_id == "S1P d18:1 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    85299
  )


  expect_message(
    mexp2 <- correct_interference_manual(
      mexp2,
      variable = "feature_intensity",
      feature = "S1P d18:1 [M>60]",
      interfering_feature = "S1P d18:2 [M>60]",
      interference_contribution = 0.0315385
    ),
    "Interference-correction was manually applied to "
  )



  expect_equal(
    mexp2@dataset |>
      filter(
        feature_id == "S1P d18:1 [M>60]",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    84304.7172490
  )

  # test renaming of feature and additional correction
  mexp2 <- correct_interference_manual(
    mexp2,
    variable = "feature_intensity",
    feature = "S1P d18:0 [M>60]",
    interfering_feature = "S1P d18:1 [M>60]",
    interference_contribution = 0.0315872,
    updated_feature_id = "S1P d18:0 [M>60] corrected"
  )

  expect_equal(
    mexp2@dataset |>
      filter(
        feature_id == "S1P d18:0 [M>60] corrected",
        analysis_id == "008_LTR_LTR01"
      ) |>
      pull(feature_intensity),
    7256.0500353124
  )


  expect_error(
    mexp2 <- correct_interference_manual(
      mexp2,
      variable = "var_undefined",
      feature = "S1P d18:0 [M>60]",
      interfering_feature = "S1P d18:1 [M>60]",
      interference_contribution = 0.0315872,
      updated_feature_id = "S1P d18:0 [M>60] corrected"
    ),
    "Variable `var_undefined` is not"
  )

  expect_error(
    mexp2 <- correct_interference_manual(
      mexp2,
      variable = "feature_intensity",
      feature = "S1P d18:0 [M>601]",
      interfering_feature = "S1P d18:1 [M>60]",
      interference_contribution = 0.0315872,
      updated_feature_id = "S1P d18:0 [M>60] corrected"
    ),
    "Selected feature is not present in the dataset"
  )

  expect_error(
    mexp2 <- correct_interference_manual(
      mexp2,
      variable = "feature_intensity",
      feature = "S1P d18:0 [M>60]",
      interfering_feature = "S1P d18:3 [M>60]",
      interference_contribution = 0.0315872,
      updated_feature_id = "S1P d18:0 [M>60] corrected"
    ),
    "Selected interfering feature is not present"
  )

  expect_error(
    mexp2 <- correct_interference_manual(
      mexp2,
      variable = "feature_intensity",
      feature = "S1P d18:0 [M>60]",
      interfering_feature = "S1P d18:1 [M>60]",
      interference_contribution = NA,
      updated_feature_id = "S1P d18:0 [M>60] corrected"
    ),
    "must be a number larger than 0"
  )

  expect_error(
    mexp2 <- correct_interference_manual(
      mexp2,
      variable = "feature_intensity",
      feature = "S1P d18:0 [M>60]",
      interfering_feature = "S1P d18:1 [M>60]",
      interference_contribution = 0.1,
      updated_feature_id = "S1P d18:2 [M>60]"
    ),
    "is already present in the dataset"
  )
})

test_that("Handles corrections that lead to negative values",{
  expect_message(
    mexp3 <- correct_interference_manual(
      mexp2,
      variable = "feature_intensity",
      feature = "PC 28:0|SM 32:1 M+3",
      interfering_feature = "SM 32:1",
      interference_contribution = 0.1,
      updated_feature_id = "PC 28:0"
    ),
    "Interference correction led to 478 negative or zero values in samples/QCs. Please verify",
    fixed = TRUE
  )

  a <- mexp3@dataset |>
    group_by(feature_id) |>
    summarise(nas = sum(is.na(feature_intensity)))

  expect_equal(max(a$nas), 0)

  expect_message(
    mexp3 <- correct_interference_manual(
      mexp2,
      variable = "feature_intensity",
      feature = "PC 28:0|SM 32:1 M+3",
      interfering_feature = "SM 32:1",
      interference_contribution = 0.1,
      neg_to_na = TRUE,
      updated_feature_id = "PC 28:0"
    ),
    "Interference correction led to 478 negative or zero values in samples/QCs. All negative/zero values",
    fixed = TRUE
  )

  a <- mexp3@dataset |>
    group_by(feature_id) |>
    summarise(nas = sum(is.na(feature_intensity)))

    expect_equal(max(a$nas), 479)

  expect_message(
    mexp_res <-
      correct_interferences(
        mexp2,
        variable = "feature_intensity",
        sequential_correction = TRUE
      ),
    "Interference correction led to negative or zero values in 1 feature(s) in samples/QCs. Please verify ",
    fixed = TRUE
  )

  expect_message(
    mexp_res <-
      correct_interferences(
        mexp2,
        variable = "feature_intensity",
        sequential_correction = TRUE,
        neg_to_na = TRUE
      ),
    "Interference correction led to negative or zero values in 1 feature(s) in samples/QCs. All negative/zero values ",
    fixed = TRUE
  )


})
