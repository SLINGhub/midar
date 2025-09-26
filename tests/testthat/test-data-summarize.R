library(ggplot2)

mexp_orig <- lipidomics_dataset
mexp_orig <- normalize_by_istd(mexp_orig)
mexp_orig <- calc_qc_metrics(mexp_orig)

mexp <- mexp_orig 
mexp@annot_features[str_detect(mexp@annot_features$feature_id, "LPC 18:1 \\((a|b)\\)"), ]$analyte_id <- "LPC 18:1"

mexp <- midar:::link_data_metadata(mexp)

test_that("Default plot_qc_matrixeffects looks as expected", {
  mexp_dedup <- data_sum_features(mexp)
  expect_true("LPC 18:1" %in% mexp_dedup@annot_features$feature_id)
  expect_false("LPC 18:1 (a)" %in% mexp_dedup@annot_features$feature_id)
  expect_false("LPC 18:1 (b)" %in% mexp_dedup@annot_features$feature_id)
  expect_true("LPC 18:1 (ab) d7 (ISTD)" %in% mexp_dedup@annot_features$feature_id)
  expect_true("LPC 18:1" %in% unique(mexp_dedup@dataset$feature_id))
  expect_false("LPC 18:1 (a)" %in% unique(mexp_dedup@dataset$feature_id))
  expect_false("LPC 18:1 (ab)" %in% unique(mexp_dedup@dataset$feature_id))

})

mexp2 <- mexp_orig
mexp2@annot_features[str_detect(mexp2@annot_features$feature_id, "^PC"), ]$analyte_id <- "PC"
mexp2@annot_features[str_detect(mexp2@annot_features$feature_id, "^PC 4"), ]$is_quantifier <- FALSE
mexp2 <- midar:::link_data_metadata(mexp2)

test_that("Default plot_qc_matrixeffects looks as expected", {
  mexp2_dedup <- data_sum_features(mexp2, qualifier_action = "separate")
  expect_true("PC" %in% unique(mexp2_dedup@dataset$feature_id))
  expect_false("PC 40:6" %in% unique(mexp2_dedup@dataset$feature_id))
  expect_false("PC 32:1" %in% unique(mexp2_dedup@dataset$feature_id))

  sum_pc <-  sum(mexp2_dedup@dataset[mexp2_dedup@dataset$feature_id == "PC", ]$feature_intensity)
  expect_equal(sum_pc, 2937988066.1)
  sum_pc <-  sum(mexp2_dedup@dataset[mexp2_dedup@dataset$feature_id == "PC_qual", ]$feature_intensity)
  expect_equal(sum_pc, 1715685212.3)

  mexp2_dedup <- data_sum_features(mexp2, qualifier_action = "include")
  expect_true("PC" %in% unique(mexp2_dedup@dataset$feature_id))
  expect_false("PC 40:6" %in% unique(mexp2_dedup@dataset$feature_id))
  expect_false("PC 32:1" %in% unique(mexp2_dedup@dataset$feature_id))

  sum_pc <-  sum(mexp2_dedup@dataset[mexp2_dedup@dataset$feature_id == "PC", ]$feature_intensity)
  expect_equal(sum_pc, 4653673278.4)

  mexp2_dedup <- data_sum_features(mexp2, qualifier_action = "exclude")
  expect_true("PC" %in% unique(mexp2_dedup@dataset$feature_id))
  expect_false("PC 40:6" %in% unique(mexp2_dedup@dataset$feature_id))
  expect_false("PC 32:1" %in% unique(mexp2_dedup@dataset$feature_id))

  sum_pc <-  sum(mexp2_dedup@dataset[mexp2_dedup@dataset$feature_id == "PC", ]$feature_intensity)
  expect_equal(sum_pc, 2937988066.1)

})