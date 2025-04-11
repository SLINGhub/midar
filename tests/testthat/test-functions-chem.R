test_that("calc_averagemolmass works", {
  expect_equal(calc_averagemolmass(c("C6H12O6", "C16H34NO5P", "C16[13C]2H36[3H]2NO5P")), c(180.15612566, 351.419229593, 379.47248))
  expect_error(calc_averagemolmass(c("16H12O6", "C16H34NO5P")), "Following invalid chemical formula defined: 16H12O6")
})
test_that("calc_averagemolmass handles empty input", {
  expect_equal(calc_averagemolmass(character(0)), numeric(0))
})
