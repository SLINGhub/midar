# library(testthat)
# library(dplyr)

mexp_orig <- lipidomics_dataset
mexp <- normalize_by_istd(mexp_orig)
mexp <- quantify_by_istd(mexp)
mexp_proc <- calc_qc_metrics(mexp,  use_batch_medians = FALSE)


test_that("get_response_curve_stats works", {
  res <- get_response_curve_stats(mexp)
  expect_s3_class(res, "tbl_df")
  expect_equal(dim(res), c(29,7))

  res <- get_response_curve_stats(mexp,
                                  limit_to_rqc = TRUE)
  expect_s3_class(res, "tbl_df")
  expect_equal(dim(res), c(29,7))

  res <- get_response_curve_stats(mexp,
                                  with_staturation_stats = TRUE)
  expect_equal(dim(res), c(29,17))

})

test_that("get_response_curve_stats missing data correctly handled", {
  mexp_temp <- mexp
  mexp_temp@annot_responsecurves <- mexp_temp@annot_responsecurves[0,]

  expect_error(get_response_curve_stats(mexp_temp),
               "No response curve metadata found")

  mexp_temp <- mexp
  mexp_temp@annot_responsecurves$analysis_id[1] <- "unknown1"
  mexp_temp@annot_responsecurves$analysis_id[3] <- "unknown2"

  expect_error(get_response_curve_stats(mexp_temp),
               "One or more analysis IDs")

  mexp_temp <- mexp
  mexp_temp@dataset <-  mexp_temp@dataset |>
    mutate(qc_type = if_else(qc_type == "RQC", "SPL", qc_type))

    expect_error(get_response_curve_stats(mexp_temp,
                                          limit_to_rqc = TRUE),
               "No analyses/samples of QC type \\`RQC")


    # check error handler of lm
    mexp_temp <- mexp
    mexp_temp@annot_responsecurves$analyzed_amount[1:6] <- Inf

    res <- get_response_curve_stats(mexp_temp,
                                          limit_to_rqc = TRUE)
    expect_true(all(is.na(res$r2_rqc_A)))
    expect_false(any(is.na(res$r2_rqc_B)))
    expect_true(all(is.na(res$slopenorm_rqc_A)))
    expect_false(any(is.na(res$slopenorm_rqc_B)))
    expect_true(all(is.na(res$y0norm_rqc_A)))
    expect_false(any(is.na(res$y0norm_rqc_B)))
})
