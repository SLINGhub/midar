test_that("multiplication works", {


  myexp <- midar::MidarExperiment(title = "sPerfect")

  data_path <- test_path("testdata/FullPanelFewSamples_MRMkit.csv")
  myexp <- import_data_mrmkit(data = myexp, path = data_path, import_metadata = TRUE)


  plot_rt_vs_chain(myexp, qc_types = "SPL", x_var = "total_c")

})



