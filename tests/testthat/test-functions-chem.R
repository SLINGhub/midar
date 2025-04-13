test_that("calc_average_molweight works", {
  expect_equal(calc_average_molweight(c("C6H12O6", "C16H34NO5P", "C16[13]C2H36D2NO5P")), c(180.15612566, 351.419229593, 383.470024240))
  expect_error(calc_average_molweight(c("16H12O6", "C16H34NO5P")), "Following invalid chemical formula defined: 16H12O6")
})
test_that("calc_average_molweight handles empty input", {
  expect_error(calc_average_molweight(character(0)),
               "No chemical formula provided. Please provide on or more valid chemical formula.", fixed = TRUE)
})
