# library(testthat)
# library(tibble)
# library(cli)
# library(pillar)

# Tests for `as_assertr_tibble`
test_that("as_assertr_tibble correctly converts data.frame to assertr_tibble", {
  df <- data.frame(a = 1:3, b = c("x", "y", "z")) |> as_tibble()
  assertr_tbl <- as_assertr_tibble(df)

  expect_s3_class(assertr_tbl, "assertr_tibble")
  expect_s3_class(assertr_tbl, "tbl_df")
  expect_s3_class(assertr_tbl, "data.frame")
})

test_that("as_assertr_tibble throws error if input is not a data.frame", {
  expect_error(as_assertr_tibble(123), "x must be a data.frame")
  expect_error(as_assertr_tibble("not a data.frame"), "x must be a data.frame")
})

# Tests for `tbl_sum.assertr_tibble`
test_that("tbl_sum.assertr_tibble returns the correct metadata summary", {
  df <- data.frame(a = 1:3, b = c("x", "y", "z"))  |> as_tibble()
  assertr_tbl <- as_assertr_tibble(df)

  expect_output(print(tbl_sum(assertr_tbl)), "-------------------------------------------------------------------------------------------")
})

test_that("tbl_format_header.assertr_tibble formats header correctly", {
  df <- data.frame(a = 1:3, b = c("x", "y", "z"))
  assertr_tbl <- as_assertr_tibble(df)

  setup <- list(tbl_sum = c(Metadata = "Errors and Warnings"))
  header <- tbl_format_header(assertr_tbl, setup)

  # Match the actual output, accounting for no spaces and ANSI strings
  expect_match(
    as.character(header),
    "Metadata-Errors and Warnings",
    fixed = TRUE
  )
})


test_that("tbl_format_footer.assertr_tibble formats footer correctly", {
  df <- data.frame(a = 1:3, b = c("x", "y", "z"))
  assertr_tbl <- as_assertr_tibble(df)

  setup <- list(rows_total = 3)
  footer <- tbl_format_footer(assertr_tbl, setup)

  # Matching only the content while ignoring the formatting (escape sequences)
  expect_output(
    print(footer),
    "E = Error, W = Warning, W\\* = Supressed Warning, N = Note"
  )
})



test_that("ctl_new_pillar.assertr_tibble creates a pillar correctly and prints it", {
  # Create a sample tibble with the assertr_tibble class
  df <- data.frame(ColumnA = 1:3, `Column B` = c("x", "y", "z"))
  assertr_tbl <- as_assertr_tibble(df)



  # Capture the output of the pillar
  captured_output <- capture_output(print(assertr_tbl))

  # Check if the title is present
  expect_true(any(grepl("Column.B", captured_output)))

  # Check if the data is included in the output
  expect_true(any(grepl("W", captured_output)))
  expect_true(any(grepl("Supressed Warning", captured_output)))
  expect_true(any(grepl("E = Error, W = Warning", captured_output)))
})


