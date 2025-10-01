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


test_that("handles short or empty vectors", {
  expect_equal(get_outlier_bounds(1), c(NA_real_, NA_real_))
  expect_equal(get_outlier_bounds(numeric(0)), c(NA_real_, NA_real_))
})

test_that("handles NA values correctly with na.rm", {
  x_na <- c(1, 2, 3, 4, 100, NA)

  # With na.rm = TRUE, NA should be removed and calculation should proceed
  expect_equal(get_outlier_bounds(x_na, "iqr", na.rm = TRUE), c(1, 4))

  # With na.rm = FALSE (default), calculations involving NAs should result in NA bounds
  expect_equal(
    get_outlier_bounds(x_na, "iqr", na.rm = FALSE),
    c(NA_real_, NA_real_)
  )
})

test_that("throws an error for an invalid method string", {
  expect_error(
    get_outlier_bounds(1:10, method = "some_invalid_method"),
    "should be one of"
  )
})

test_that("handles vectors with all identical values", {
  x <- rep(5, 10)
  # For identical values, IQR, SD, and MAD are 0, so bounds should be the value itself.
  expect_equal(get_outlier_bounds(x, "iqr"), c(5, 5))
  expect_equal(get_outlier_bounds(x, "sd"), c(5, 5))
  expect_equal(get_outlier_bounds(x, "mad"), c(5, 5))
  expect_equal(get_outlier_bounds(x, "quantile"), c(5, 5))
  expect_equal(get_outlier_bounds(x, "z_robust"), c(5, 5))
})


test_that("method 'iqr' calculates bounds correctly", {
  x <- c(1, 2, 3, 4, 100)
  # Q1=2, Q3=4, IQR=2. Lower fence = 2 - 1.5*2 = -1. Upper fence = 4 + 1.5*2 = 7.
  # Smallest value >= -1 is 1. Largest value <= 7 is 4.
  expect_equal(get_outlier_bounds(x, "iqr"), c(1, 4))

  # Test with a custom k that includes the outlier
  # Upper fence = 4 + 50*2 = 104.
  expect_equal(get_outlier_bounds(x, "iqr", k = 50), c(1, 100))
})

test_that("methods 'sd' and 'z_normal' calculate bounds correctly", {
  x <- c(1, 2, 3, 4, 100)
  # mean=22, sd=43.618. With k=3, fences are approx. -108.85 and 152.85.
  # All values are within these fences.
  expect_equal(get_outlier_bounds(x, "sd"), c(1, 100), tolerance = 1e-6)
  expect_equal(get_outlier_bounds(x, "z_normal"), c(1, 100), tolerance = 1e-6)

  # Test with a custom k that excludes the outlier
  # With k=1, fences are approx. -21.6 and 65.6. Largest value <= 65.6 is 4.
  expect_equal(get_outlier_bounds(x, "sd", k = 1), c(1, 4), tolerance = 1e-6)
})

test_that("method 'mad' calculates bounds correctly", {
  x <- c(1, 2, 3, 4, 100)
  # median=3, mad=1. With k=3, fences are 0 and 6.
  # Smallest value >= 0 is 1. Largest value <= 6 is 4.
  expect_equal(get_outlier_bounds(x, "mad"), c(1, 4))

  # Test with a custom k that includes the outlier
  # With k=100, fences are -97 and 103.
  expect_equal(get_outlier_bounds(x, "mad", k = 100), c(1, 100))
})

test_that("method 'quantile' calculates bounds correctly", {
  x <- c(1, 2, 3, 4, 100)
  # Default k=0.01. Quantiles are approx. 1.04 and 96.16.
  # Smallest value >= 1.04 is 2. Largest value <= 96.16 is 4.
  expect_equal(get_outlier_bounds(x, "quantile"), c(2, 4))

  # Test with a custom k
  x_long <- 1:100
  # k=0.05. Quantiles are approx. 5.95 and 95.05.
  # Smallest value >= 5.95 is 6. Largest value <= 95.05 is 95.
  expect_equal(get_outlier_bounds(x_long, "quantile", k = 0.05), c(6, 95))
})

test_that("method 'z_robust' returns min/max of inliers", {
  x <- c(1, 2, 3, 4, 100)
  # This method returns min/max of values whose modified z-score is <= k.
  # With k=3.5, only 1, 2, 3, 4 are included.
  expect_equal(get_outlier_bounds(x, "z_robust"), c(1, 4))

  # With a large k, all values are included.
  expect_equal(get_outlier_bounds(x, "z_robust", k = 70), c(1, 100))
})

test_that("method 'fold_change' calculates bounds correctly", {
  x <- log10(c(1, 2, 4, 8)) # Example from documentation
  # median=0.4515. k=2 -> delta=log2(2)=1. Fences are -0.5485 and 1.4515.
  # All values are within the fences.
  expect_equal(get_outlier_bounds(x, "fold_change"), range(x), tolerance = 1e-6)

  # Test with a custom k
  x2 <- log10(c(1, 2, 4, 8, 100))
  # median=log10(4)=0.602. k=3 -> delta=log2(3)=1.585.
  # Fences are approx. -0.98 and 2.187.
  # log10(100) = 2. All values are within these fences.
  expect_equal(
    get_outlier_bounds(x2, "fold_change", k = 3),
    range(x2),
    tolerance = 1e-6
  )

  expect_equal(
    get_outlier_bounds(x2, "fold_change", k = c(-3,3)),
    range(x2),
    tolerance = 1e-6
  )
})


test_that("outlier_log transforms data before calculation", {
  x <- c(10, 100, 1000, 10000, 10000000)
  x_log <- c(1, 2, 3, 4, 7)
  # Calculation should be on log10 data: c(1, 2, 3, 4, 7)
  # Using 'iqr': Q1=2, Q3=4, IQR=2. Fences are -1 and 7.
  # The function returns bounds on the transformed scale.
  # Smallest value >= -1 is 1. Largest value <= 7 is 7.
  expect_equal(get_outlier_bounds(x, "iqr", outlier_log = TRUE), c(1, 7))
})

test_that("outlier_log throws an error for non-positive values", {
  msg <- "All values must be positive for log transformation."
  expect_error(get_outlier_bounds(c(-1, 1, 2), "iqr", outlier_log = TRUE), msg)
  expect_error(get_outlier_bounds(c(0, 1, 2), "iqr", outlier_log = TRUE), msg)
})


# Assuming the corrected get_mad_tails function is loaded
# source("R/get_mad_tails.R")

test_that("handles short or empty vectors", {
  expect_equal(get_mad_tails(numeric(0), k = 1.5), c(NA_real_, NA_real_))
  expect_equal(get_mad_tails(5, k = 1.5), c(NA_real_, NA_real_))
})

test_that("handles NA values correctly", {
  x_na <- c(-100, 1, 2, 3, 4, 100, NA)
  x <- c(-100, 1, 2, 3, 4, 100)
  k <- 1.5

  # With na.rm = TRUE, should be identical to calculation without NA
  expect_equal(get_mad_tails(x_na, k, na.rm = TRUE), get_mad_tails(x, k))

  # With na.rm = FALSE (default), should propagate NA
  expect_equal(get_mad_tails(x_na, k, na.rm = FALSE), c(NA_real_, NA_real_))
})

test_that("handles vectors with all identical values", {
  x <- rep(5, 10)
  # mad is 0, fences are both 5. No values are > 5 or < 5.
  # Corrected function should return NA for both bounds.
  expect_equal(get_mad_tails(x, k = 1.5), c(NA_real_, NA_real_))
})

test_that("handles vectors with only two distinct values", {
  x <- c(1, 1, 1, 10, 10, 10)
  k <- 1.5
  # median=5.5, mad=4.5. With default constant, mad is ~6.67.
  # Fences will be outside the range of data.
  # lo should be min(x) and up should be max(x).
  expect_equal(get_mad_tails(x, k), c(1, 10))
})


test_that("calculates tails correctly for a standard case", {
  x <- c(-100, 1, 2, 3, 4, 100)
  k <- 1.5

  # Manual calculation:
  # median = 2.5
  # mad = 2.2239 (using R's default constant)
  # lower_fence = 2.5 - 1.5 * 2.2239 = -0.83585
  # upper_fence = 2.5 + 1.5 * 2.2239 = 5.83585
  # lo = min(x[x > -0.83585]) = 1
  # up = max(x[x < 5.83585]) = 4
  expect_equal(get_mad_tails(x, k), c(1, 4), tolerance = 1e-6)
})

test_that("behaves correctly with a large k value", {
  x <- c(-100, 1, 2, 3, 4, 100)
  k <- 100 # A very large multiplier

  # With a large k, fences will be far outside the data range.
  # All points will be > lower_fence and < upper_fence.
  # Therefore, lo should be min(x) and up should be max(x).
  expect_equal(get_mad_tails(x, k), c(-100, 100))
})

test_that("behaves correctly with a small k value", {
  x <- c(-100, 1, 2, 3, 4, 100)
  k <- 0.1 # A very small multiplier

  # Manual calculation:
  # median = 2.5
  # mad = 2.2239
  # lower_fence = 2.5 - 0.1 * 2.2239 = 2.27761
  # upper_fence = 2.5 + 0.1 * 2.2239 = 2.72239
  # lo = min(x[x > 2.27761]) = 3
  # up = max(x[x < 2.72239]) = 2
  # Note: The lower bound is > upper bound, which is valid.
  expect_equal(get_mad_tails(x, k), c(3, 2), tolerance = 1e-6)
})

test_that("works with a simple symmetric vector", {
  x <- 1:11
  k <- 1.5

  # median = 6
  # mad = 3 (unscaled), ~4.4478 (scaled)
  # lower_fence = 6 - 1.5 * 4.4478 = -0.6717
  # upper_fence = 6 + 1.5 * 4.4478 = 12.6717
  # lo = min(x[x > -0.6717]) = 1
  # up = max(x[x < 12.6717]) = 11
  expect_equal(get_mad_tails(x, k), c(1, 11), tolerance = 1e-6)
})


test_that("handles short or empty vectors", {
  expect_equal(get_iqr_tails(numeric(0)), c(NA_real_, NA_real_))
  expect_equal(get_iqr_tails(10), c(NA_real_, NA_real_))
})

test_that("handles NA values correctly", {
  x_na <- c(-100, -1, -2, -3, 1, 2, 3, 4, 100, NA)
  x <- c(-100, -1, -2, -3, 1, 2, 3, 4, 100)

  # With na.rm = TRUE, should be identical to calculation without NA
  expect_equal(get_iqr_tails(x_na, na.rm = TRUE), get_iqr_tails(x))

  # With na.rm = FALSE (default), quantile() returns NA, which propagates
  expect_equal(get_iqr_tails(x_na, na.rm = FALSE), c(NA_real_, NA_real_))
})

test_that("handles vectors with all identical values", {
  x <- rep(7, 10)
  # Q1=7, Q3=7, IQR=0. Fences are both 7.
  # min(x[x >= 7]) is 7. max(x[x <= 7]) is 7.
  expect_equal(get_iqr_tails(x), c(7, 7))
})

test_that("handles edge case with negative k leading to empty subsets", {
  # This test confirms the current behavior (returning Inf/-Inf).
  # A more robust function might return c(NA, NA) instead.
  x <- 1:10
  # With k=-3, fences are inverted and outside the data range.
  # min(numeric(0)) -> Inf; max(numeric(0)) -> -Inf
  suppressWarnings({
    # Suppress "no non-missing arguments to min/max" warnings
    expect_equal(get_iqr_tails(x, k = -3), c(Inf, -Inf))
  })
})


test_that("calculates tails correctly for a standard case", {
  x <- c(-100, -1, -2, -3, 1, 2, 3, 4, 100)

  # Manual calculation:
  # sorted x = -100, -3, -2, -1, 1, 2, 3, 4, 100
  # Q1 = -2, Q3 = 3, IQR = 5
  # lower_fence = -2 - 1.5 * 5 = -9.5
  # upper_fence = 3 + 1.5 * 5 = 10.5
  # lo = min(x[x >= -9.5]) = -3
  # up = max(x[x <= 10.5]) = 4
  expect_equal(get_iqr_tails(x, k = 1.5), c(-3, 4))
})

test_that("calculates tails correctly with a custom k", {
  x <- c(1, 2, 3, 4, 5, 100)

  # With default k=1.5:
  # Q1=2, Q3=4, IQR=2. Fences are -1 and 7. Tails are 1 and 5.
  expect_equal(get_iqr_tails(x), c(1, 5))

  # With a large k that includes the outlier:
  # k=50 -> Fences are -98 and 104. Tails are 1 and 100.
  expect_equal(get_iqr_tails(x, k = 50), c(1, 100))
})


# ---------------------------------------------------------------------------
# Tests for find_closest
# ---------------------------------------------------------------------------

test_that("throws an error for an invalid method", {
  expect_error(
    find_closest(5, 1:10, method = "wrong"),
    "should be one of"
  )
})

test_that("handles empty available_numbers vector", {
  # which.min(numeric(0)) returns integer(0), subsetting returns numeric(0)
  expect_equal(length(find_closest(5, numeric(0), "absolute")), 0)
  # min(numeric(0)) returns Inf
  suppressWarnings({
    expect_equal(find_closest(5, numeric(0), "lower"), Inf)
    expect_equal(find_closest(5, numeric(0), "higher"), -Inf)
  })
})


test_that("absolute method finds the correct closest number", {
  nums <- c(1, 2, 8, 10)
  expect_equal(find_closest(5, nums, "absolute"), 2)
  expect_equal(find_closest(8.9, nums, "absolute"), 8)
  expect_equal(find_closest(9.1, nums, "absolute"), 10)
})

test_that("absolute method breaks ties by picking the first occurrence", {
  nums <- c(4, 6, 10)
  # The distance to 4 and 6 is identical (1). `which.min` picks the first one.
  expect_equal(find_closest(5, nums, "absolute"), 4)

  nums_rev <- c(6, 4, 10)
  # With reversed order, it should now pick 6.
  expect_equal(find_closest(5, nums_rev, "absolute"), 6)
})


test_that("lower method finds the largest number less than or equal to x", {
  nums <- c(1, 5, 10, 15)
  expect_equal(find_closest(11, nums, "lower"), 10)
  # Should include the number itself if it's a perfect match
  expect_equal(find_closest(10, nums, "lower"), 10)
})

test_that("lower method handles out-of-bounds case correctly", {
  nums <- c(10, 20, 30)
  # If x is smaller than all available numbers, it should return the min of available numbers.
  expect_equal(find_closest(5, nums, "lower"), 10)
})


test_that("higher method finds the smallest number greater than or equal to x", {
  nums <- c(1, 5, 10, 15)
  expect_equal(find_closest(6, nums, "higher"), 10)
  # Should include the number itself if it's a perfect match
  expect_equal(find_closest(5, nums, "higher"), 5)
})

test_that("higher method handles out-of-bounds case correctly", {
  nums <- c(10, 20, 30)
  # If x is larger than all available numbers, it should return the max of available numbers.
  expect_equal(find_closest(35, nums, "higher"), 30)
})
