test_that("some_na works", {
  expect_true(some_na(c(NA, 1, 2)))
  expect_false(some_na(c(NA, NA, NA)))
  expect_false(some_na(c(1, 2, 3)))
  expect_false(some_na(NA))
  expect_false(some_na(1))
  expect_false(some_na(NULL))
})

test_that("safe_min works", {
  expect_equal(safe_min(c(4, 3, 2, 1)), 1)
  expect_equal(safe_min(c(NA, NA, 3, 2, 1)), NA_real_)
  expect_equal(safe_min(c(NA, NA, 3, 2, 1), na.rm = TRUE), 1)
  expect_equal(safe_min(c(NaN, NaN, 3, 2, 1)), NA_real_)
  expect_equal(safe_min(c(NaN, NaN, 3, 2, 1), na.rm = TRUE), 1)
  expect_equal(safe_min(c(NA, NaN, 1, 2, 3), na.rm = FALSE), NA_real_)
  expect_equal(safe_min(c(NA, NaN, 1, 2, 3), na.rm = TRUE), 1)
  expect_equal(safe_min(c(NA, NaN, NA)), NA_real_)
  expect_equal(safe_min(c(NaN, NaN, NaN)), NA_real_)
})

test_that("safe_max works", {
  expect_equal(safe_max(c(4, 3, 2, 1)), 4)
  expect_equal(safe_max(c(NA, NA, 3, 2, 1)), NA_real_)
  expect_equal(safe_max(c(NA, NA, 3, 2, 1), na.rm = TRUE), 3)
  expect_equal(safe_max(c(NaN, NaN, 3, 2, 1)), NA_real_)
  expect_equal(safe_max(c(NaN, NaN, 3, 2, 1), na.rm = TRUE), 3)
  expect_equal(safe_max(c(NA, NaN, 1, 2, 3), na.rm = FALSE), NA_real_)
  expect_equal(safe_max(c(NA, NaN, 1, 2, 3), na.rm = TRUE), 3)
  expect_equal(safe_max(c(NA, NaN, NA)), NA_real_)
  expect_equal(safe_max(c(NaN, NaN, NaN)), NA_real_)
})


test_that("check_groupwise_identical_ids works", {
  df_identical <- tibble::tibble(
    group = c("A", "A", "A", "B", "B"),
    id = c(1, 1, 1, 2, 2))
  expect_true(check_groupwise_identical_ids(df_identical, group_col = group, id_col = id))

  df_non_identical <- tibble::tibble(
    group = c("A", "A", "A", "B", "B"),
    id = c(1, 2, 1, 2, 3))
  expect_false(check_groupwise_identical_ids(df_non_identical, group_col = group, id_col = id))

  df_missing <- tibble::tibble(
    group = c("A", "A", "B", "B"),
    id = c(1, NA, 2, 3))
  expect_false(check_groupwise_identical_ids(df_missing, group_col = group, id_col = id))

  df_single <- tibble::tibble(group = "A", id = 1)
  expect_true(check_groupwise_identical_ids(df_single, group_col = group, id_col = id))

  df_empty <- tibble::tibble(
    group = character(0),
    id = integer(0))
  expect_error(check_groupwise_identical_ids(df_empty, group_col = group, id_col = id),
               "data has no rows")
})

test_that("comp_val works", {
  tbl <- tibble::tibble(
    value1 = c(1, 2, NA, 4),
    value2 = c(5, NA, 7, 8)
  )

  expect_equal(comp_val(tbl, val = "non_existing_column", threshold = 5, operator = ">"), NA)

  tbl_with_na <- tibble::tibble(value1 = c(NA, NA, NA, NA))
  expect_equal(comp_val(tbl_with_na, val = "value1", threshold = NA, operator = ">"), c(NA, NA, NA, NA))

  expect_equal(comp_val(tbl, val = "value1", threshold = 3, operator = ">"), c(FALSE, FALSE, NA, TRUE))
  expect_equal(comp_val(tbl, val = "value1", threshold = 3, operator = "<"), c(TRUE, TRUE, NA, FALSE))
  expect_equal(comp_val(tbl, val = "value1", threshold = 2, operator = "=="), c(FALSE, TRUE, NA, FALSE))
  expect_equal(comp_val(tbl, val = "value2", threshold = 8, operator = "=="), c(FALSE, NA, FALSE, TRUE))

  expect_equal(comp_val(tbl, val = "noexist_col", threshold = 10, operator = ">"), NA)

  df_empty <- tibble::tibble(a = character(0),b = integer(0))
  expect_error(comp_val(df_empty, val = "value1", threshold = 3, operator = ">"),
               "tbl has no rows")
})


test_that("comp_lgl_vec works as it should", {
  expect_equal(comp_lgl_vec(list(c(TRUE, TRUE, TRUE), c(FALSE, TRUE, TRUE)),
                           .operator = "AND"), c(FALSE, TRUE, TRUE))

  expect_equal(comp_lgl_vec(list(c(TRUE, TRUE, TRUE), c(FALSE, TRUE, TRUE)),
                            .operator = "OR"), c(TRUE, TRUE, TRUE))

  expect_equal(comp_lgl_vec(list(c(TRUE, TRUE, TRUE), c(FALSE, TRUE, TRUE)),
                            .operator = "XOR"), c(TRUE, FALSE, FALSE))

  expect_equal(comp_lgl_vec(list(c(NA, NA, NA), c(NA, NA, NA)),
                            .operator = "AND"), c(NA, NA, NA))

  expect_null(comp_lgl_vec(list(c(TRUE, FALSE, TRUE), c(TRUE, TRUE, FALSE)), .operator = "XAND"))
  expect_null(comp_lgl_vec(list(), .operator = "AND"))

})

test_that("has_any_name works in assertr::verify as it should", {
  dt <- tibble(
    col_a = c(1, 2, 3, 4, 5),
    col_b = c(1, 2, 3, 4, 5),
    col_c = c(1, 2, 3, 4, 5)
  )
  expect_equal(dim (dt |> assertr::verify(has_any_name("col_a"), obligatory=TRUE, description = "")), c(5, 3))
  expect_equal(dim (dt |> assertr::verify(has_any_name("col_a", "col_b"), obligatory=TRUE, description = "")), c(5, 3))
  res <- dt |> assertr::verify(has_any_name("col_noexist"), obligatory=TRUE, description = "", error_fun = assertr::error_df_return)
  expect_equal(dim(res), c(1, 6)) # means it is an rrror deta frame
  res <- dt |> assertr::verify(has_any_name("col_a", "col_noexist"), obligatory=TRUE, description = "", error_fun = assertr::error_df_return)
  expect_equal(dim(res), c(5, 3))
})

test_that("add_missing_column works", {
  # Create a sample data frame without the target column
  dt <- tibble(A = 1:5, B = 6:10)

  result <- add_missing_column(dt, "c", 99, make_lowercase = FALSE)
  expect_equal(result$c, rep(99, 5))

  result <- add_missing_column(dt, "A", 99, make_lowercase = TRUE)
  expect_equal(result$a, 1:5)
  result <- add_missing_column(dt, "A", 99, make_lowercase = FALSE)
  expect_equal(result$A, 1:5)
})


test_that("get_conc_unit works as expected", {
  expect_equal(get_conc_unit("ul"), "\U003BCmol/L")
  expect_equal(get_conc_unit("mL"), "pmol/mL")
  expect_equal(get_conc_unit("L"), "pmol/L")
  expect_equal(get_conc_unit(c("ul", "ml")), "pmol/sample amount unit (multiple units)")
  expect_equal(get_conc_unit("mg"), "pmol/mg")
  expect_equal(get_conc_unit("Ul"), "\U003BCmol/L")
})
