
#' Percent coefficient of variation (%CV)
#'
#' This function computes the percent coefficient of variation  of the values in x.
#' If na.rm is TRUE then missing values are removed before computation proceeds.
#'
#' @param x a numeric vector with untransformed data
#' @param na.rm logical, if TRUE then NA values are stripped from x before
#' computation takes place
#' @param use_robust_cv logical, if TRUE then the robust coefficient of variation
#'
#' @return a numeric value. If x contains a zero or is not numeric,
#' NA_real_ is returned

#' @export
#'
#' @examples
#' cv(c(5, 6, 3, 4, 5, NA), na.rm = TRUE)

cv <- function(x, na.rm = FALSE, use_robust_cv = FALSE){
  if(use_robust_cv){
    mad(x, na.rm = na.rm, constant = 1)/median(x, na.rm = na.rm) * 100
  } else {
    sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm) * 100
  }

}

#' Percent coefficient of variation (%CV) based on log-transformation
#'
#' Computes the percent coefficient of variation (CV) based on the log-transformation of
#' the input data. This can be more robust in certain cases, especially when the
#' data is not normally distributed.
#'
#' @param x a numeric vector with untransformed data
#' @param na.rm logical, if TRUE then NA values are stripped from x before
#' computation takes place
#'
#' @return a numeric value. If x contains a zero, then NaN is returned. If x is
#' not numeric, NA_real_ is returned
#' @references
#' Canchola et al. (2017) Correct use of percent coefficient of variation (% CV) formula for log-transformed data. MOJ Proteomics Bioinform. 2017;6(4):316‒317.
#' [DOI: 10.15406/mojpb.2017.06.00200](https://doi.org/10.15406/mojpb.2017.06.00200)
#' @export
#'
#' @examples
#' cv_log(c(5, 6, 3, 4, 5, NA), na.rm = TRUE)
#'
cv_log <- function(x, na.rm = FALSE){
  # sqrt(exp(1)^(sd(log(x, ...))^2)-1) *100
  sqrt(10^(log(10) * stats::sd(log(x, 10), na.rm)^2) - 1) * 100
}


#' Get MAD-based tails
#'
#' Computes the lower and upper boundaries based on Median Absolute Deviation
#' (MAD) from a numeric vector. Returns c(NA_real_, NA_real_) if the input vector
#' has less than 2 elements.
#'
#' @param x a numeric vector
#' @param k multiplier for the MAD
#' @param na.rm if TRUE then NA values are stripped from x before computation takes place.
#'
#' @return a numeric vector of length 2 with the lower and upper boundaries.
#'
#' @export
#'
#' @examples
#' x <- c(-100,1, 2, 3, 4, 100)
#' k <- 1.5
#' get_mad_tails(x, k)

# https://stackoverflow.com/questions/9843660/marking-the-very-end-of-the-two-whiskers-in-each-boxplot-in-ggplot2-in-r-statist
get_mad_tails <- function(x, k, na.rm = FALSE) {
  if (length(x) < 2) return(c(NA_real_, NA_real_))

  med <- median(x, na.rm = na.rm)
  MAD <- mad(x, na.rm = na.rm)
  upper <- med + k * MAD
  lower <- med - k * MAD

  ## Trim upper and lower
  up <- max(x[x < upper], na.rm = na.rm)
  lo <- min(x[x > lower], na.rm = na.rm)
  return(c(lo, up))
}


# # Flag outliers, based on Tukey’s IQR fences
# flag_outlier_iqr <- function(data, include_calibdata, limit_iqr = 1.5) {
#   data <- data |>
#     group_by(.data$ceramideName, .data$SampleType) |>
#     mutate(
#       IQR_sp = IQR(.data$C_SinglePoint_mean, na.rm = TRUE),
#       Q1_sp = quantile(.data$C_SinglePoint_mean, 0.25, na.rm = TRUE),
#       Q3_sp = quantile(.data$C_SinglePoint_mean, 0.75, na.rm = TRUE),
#       Outlier_sp = !dplyr::between(.data$C_SinglePoint_mean, (.data$Q1_sp - limit_iqr * .data$IQR_sp), (.data$Q3_sp + limit_iqr * .data$IQR_sp)),
#     ) |>
#     ungroup()
#   data
# }


