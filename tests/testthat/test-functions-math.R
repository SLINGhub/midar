test_that("cv works", {
  expect_equal(cv(c(5, 6, 3, 4, 5, NA), na.rm = TRUE), 24.7864223)
  expect_equal(cv(c(5, 6, 3, 4, 5, NA), na.rm = FALSE), NA_real_)
  expect_equal(cv(NA, na.rm = FALSE), NA_real_)
  expect_equal(cv(1, na.rm = TRUE), NA_real_)
})

test_that("cv_log works", {
  expect_equal(cv_log(c(5, 6, 3, 4, 5, NA), na.rm = TRUE), 27.0819963)
  expect_equal(cv_log(c(5, 6, 3, 4, 5, NA), na.rm = FALSE), NA_real_)
  expect_equal(cv_log(NA, na.rm = FALSE), NA_real_)
})



test_that("get_mad_tails returns correct values", {
  k <- 2
  x <- c(-100, 1, 1.2, 2, 3, 4, 5, 200)
  expect_equal(get_mad_tails(x, k), c(1, 5))

  expect_equal(get_mad_tails(c(4), k, TRUE), c(NA_real_, NA_real_))
  expect_equal(get_mad_tails(c(4), k, FALSE), c(NA_real_, NA_real_))

  x <- c(-100, NA, 1, 1.2, 2, 3, 4, 5, 200)
  expect_equal(get_mad_tails(x, k, na.rm = TRUE), c(1, 5))

  expect_equal(get_mad_tails(x, k, FALSE), c(NA_real_, NA_real_))

})
