mexp <- readRDS(file = testthat::test_path("testdata/MHQuant_demo.rds"))

test_that("get_lipid_class_names works", {
  mexp_temp <- mexp
  mexp_temp <- get_lipid_class_names(mexp_temp, use_as_feature_class = "lipid_class", add_transition_names = TRUE)
  expect_true(all(c("analyte_name", "lipid_class", "lipid_class_by_lcb", "lipid_class_base", "transition_name", "transition_group") %in% colnames(mexp_temp@dataset)))

  expect_equal(mexp_temp@dataset$analyte_name[1], "S1P d16:1")
  expect_equal(mexp_temp@dataset$lipid_class[1], "S1P d")
  expect_equal(mexp_temp@dataset$lipid_class_by_lcb[1], "S1P d16:1")
  expect_equal(mexp_temp@dataset$lipid_class_base[1], "S1P")
  expect_equal(mexp_temp@dataset$transition_name[1], "M>60")
  expect_equal(mexp_temp@dataset$transition_group[1], 2)
  expect_equal(mexp_temp@dataset$feature_class[1], "S1P d")

  mexp_temp <- get_lipid_class_names(mexp_temp, use_as_feature_class = "lipid_class", add_transition_names = FALSE)
  expect_true(all(c("analyte_name", "lipid_class", "lipid_class_by_lcb", "lipid_class_base") %in% colnames(mexp_temp@dataset)))
})
