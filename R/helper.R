
# check if one or more but not all values are NA
some_na <- function(x){
  !all(is.na(x)) & any(is.na(x))
}

min_val <- function(x, na.rm = FALSE){
  if (all(is.na(x) | is.nan(x))) NA_real_ else min(x, na.rm = na.rm)
}

max_val <- function(x, na.rm = FALSE) {
  if (all(is.na(x) | is.nan(x))) NA_real_ else max(x, na.rm = na.rm)
}


# Custom assertr function to test if at least one of provided columns exists
has_any_name = function(...){
  check_this <- list(...)
  parent <- parent.frame()
  given_names <- rlang::env_names(parent$.top_env)
  given_names <- given_names[given_names != ".data"]
  any(check_this %in% given_names)
}

add_missing_column <- function(data, col_name, init_value, make_lowercase, all_na_replace = FALSE) {
  if (!tolower(col_name) %in% tolower(names(data))) {
    data |> tibble::add_column({{ col_name }} := init_value)
  } else {
    if (make_lowercase) data <- data |> dplyr::rename_with(tolower, dplyr::matches(c(col_name), ignore.case = TRUE))
    if (all_na_replace & all(is.na(data |> pull({{col_name}}))))  data <- data |> mutate({{col_name}} := init_value)
    data
  }
}


# https://stackoverflow.com/questions/9843660/marking-the-very-end-of-the-two-whiskers-in-each-boxplot-in-ggplot2-in-r-statist
get_tails <- function(x) {
  q1 <- quantile(x, na.rm = TRUE)[2]
  q3 <- quantile(x, na.rm = TRUE)[4]
  iqr <- q3 - q1
  upper <- q3 + 1.5 * iqr
  lower <- q1 - 1.5 * iqr
  if (length(x) == 1) {
    return(x)
  } # will deal with abnormal marks at the periphery of the plot if there is one value only
  ## Trim upper and lower
  up <- max(x[x < upper])
  lo <- min(x[x > lower])
  return(c(lo, up))
}


# Flag outliers, based on Tukey’s IQR fences
flag_outlier_iqr <- function(data, include_calibdata, limit_iqr = 1.5) {
  data <- data |>
    group_by(.data$ceramideName, .data$SampleType) |>
    mutate(
      IQR_sp = IQR(.data$C_SinglePoint_mean, na.rm = TRUE),
      Q1_sp = quantile(.data$C_SinglePoint_mean, 0.25, na.rm = TRUE),
      Q3_sp = quantile(.data$C_SinglePoint_mean, 0.75, na.rm = TRUE),
      Outlier_sp = !dplyr::between(.data$C_SinglePoint_mean, (.data$Q1_sp - limit_iqr * .data$IQR_sp), (.data$Q3_sp + limit_iqr * .data$IQR_sp)),
    ) |>
    ungroup()
  data
}


# https://dewey.dunnington.ca/post/2018/modifying-facet-scales-in-ggplot2/

FacetEqualWrap <- ggplot2::ggproto(
  "FacetEqualWrap", FacetWrap,
  train_scales = function(self, x_scales, y_scales, layout, data, params) {
    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      cli::cli_abort("X and Y scales required for facet_equal_wrap")
    }

    # regular training of scales
    ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)

    # switched training of scales (x and y and y on x)
    for (layer_data in data) {
      match_id <- match(layer_data$PANEL, layout$PANEL)

      x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
      y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))

      SCALE_X <- layout$SCALE_X[match_id]
      ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)

      SCALE_Y <- layout$SCALE_Y[match_id]
      ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
    }
  }
)
#' @importFrom ggplot2 ggproto
facet_wrap_equal <- function(...) {
  # take advantage of the sanitizing that happens in facet_wrap
  facet_super <- facet_wrap(...)

  ggplot2::ggproto(NULL, FacetEqualWrap,
    shrink = facet_super$shrink,
    params = facet_super$params
  )
}


cv <- function(x, ...){
  sd(x, ...)/mean(x, ...)
}
