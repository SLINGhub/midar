test_that("data_load_example works", {
  mexp <- data_load_example()
  expect_true(is.data.frame(mexp@dataset))
  expect_equal(class(mexp) |> as.character(), "MidarExperiment")
  expect_equal(dim(mexp@dataset), c(14471, 18))

  expect_error(mexp <- data_load_example(dataset = 2), "Only dataset 1 is currently available")
})
