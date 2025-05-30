mexp <- midar::MidarExperiment(title = "")

data_path <- test_path("testdata/FullPanelFewSamples_MRMkit.csv")
mexp <- import_data_mrmkit(data = mexp, path = data_path, import_metadata = TRUE)

test_that("parse_lipid_feature_names works", {
  mexp_temp <- mexp

  mexp_temp@dataset <- parse_lipid_feature_names(mexp_temp@dataset, use_as_feature_class = "lipid_class_lcb", add_transition_names = TRUE)
  expect_true(all(c("analyte_name", "lipid_class", "lipid_class_lcb", "lipid_class_base", "transition_name", "transition_group_id") %in% colnames(mexp_temp@dataset)))

  lipids <- mexp_temp@dataset |> filter(analysis_id == "Longit_batch5_TQC38")


  expect_equal(lipids$analyte_name[44], "Cer 18:1;O2/22:0")
  expect_equal(lipids$analyte_name[70], "DG 16:0_18:1")
  expect_equal(lipids$analyte_name[155], "LPC 18:1/0:0")
  expect_equal(lipids$analyte_name[323], "PE P-16:0/20:4")
  expect_equal(lipids$analyte_name[440], "TG 18:0_32:2")
  expect_equal(lipids$analyte_name[436], "TG 50:0")
  expect_equal(lipids$lipid_class[44], "Cer")
  expect_equal(lipids$lipid_class_lcb[44], "Cer;O2")
  expect_equal(lipids$feature_class[44], "Cer;O2")
  expect_equal(lipids$lipid_class_base[44], "SP")
  expect_equal(lipids$transition_name[323], "-FA-HG")
  expect_equal(lipids$transition_group_id[323], 1)
  expect_equal(lipids$analyte_name[322], "PE P-16:0/18:2")
  expect_equal(lipids$transition_name[322], "-FA")
  expect_equal(lipids$transition_group_id[322],2)

  mexp_temp@dataset <- parse_lipid_feature_names(mexp_temp@dataset,
                                    use_as_feature_class = "lipid_class",
                                    add_transition_names = FALSE,
                                    add_chain_composition = FALSE )
  expect_true(all(c("analyte_name", "lipid_class", "lipid_class_lcb", "lipid_class_base") %in% colnames(mexp_temp@dataset)))
  expect_equal(mexp_temp@dataset$feature_class[44], "Cer")
  expect_false(any(c("transition_name", "transition_group_id") %in% colnames(mexp_temp@dataset)))

})
