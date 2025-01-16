
# check if one or more but not all values are NA
some_na <- function(x){
  !all(is.na(x)) & any(is.na(x))
}

safe_min <- function(x, na.rm = FALSE){
  if (all(is.na(x) | is.nan(x))) NA_real_ else min(x, na.rm = na.rm)
}

safe_max <- function(x, na.rm = FALSE) {
  if (all(is.na(x) | is.nan(x))) NA_real_ else max(x, na.rm = na.rm)
}

# checks if all values of a specific column are identical in a group
# return NA if data frame has no data
check_groupwise_identical_ids <- function(data, group_col, id_col) {
  if(nrow(data) == 0) stop("data has no rows")
  data |>
    summarise(all_identical = dplyr::n_distinct({{id_col}}) == 1, .by = {{group_col}}) |>
    pull(.data$all_identical) |> all()
}

# Used for qc filtering ####
# Function to used to compare qc values with criteria and deal with NA


#The comp_val function compares a column in a data frame to a threshold using a
#specified operator. It handles NA values by returning NA when both are NA, FALSE
#when the column is NA and the threshold is numeric, and applies the operator
#(e.g., >, <, ==) when both are numeric. If the column does not exist, it returns NA.

# Behaviour:
# value is NA , threshold NA -> NA
# value is num , threshold NA -> NA
# value is NA , threshold is Num -> FALSE
# value is num , threshold is num -> TRUE/FAL

# TODO: Add to function description,
# TODO: make this function public for user to build own?

comp_val <- function(tbl, val, threshold, operator) {
  if(nrow(tbl) == 0) stop("tbl has no rows")
  if (!val %in% names(tbl)) return(NA)
  v_val <- tbl[[val]]
  result <- get(operator)(v_val, threshold)
}


# performs element-wise logical operations (AND or OR) across multiple
# logical vectors in a list. It returns a vector of the results,
# where each element is the result of applying the specified operation to the
# corresponding elements of the input vectors.
# If any element is NA, the result for that position will also be NA.
# If list is empty NULL is returned
comp_lgl_vec <- function(lgl_list, .operator){
  if (.operator == "AND") {
    return(Reduce("&", lgl_list))
  } else if (.operator == "OR") {
    return(Reduce("|", lgl_list))
  } else if (.operator == "XOR") {
    return(Reduce(function(x, y) xor(x, y), lgl_list))
  } else {
    # Return NULL for unsupported operators
    return(NULL)
  }
}

# Custom assertr function to test if at least one of provided columns exists
has_any_name = function(...){
  check_this <- list(...)
  parent <- parent.frame()
  given_names <- rlang::env_names(parent$.top_env)
  given_names <- given_names[given_names != ".data"]
  any(check_this %in% given_names)
}


# Add a new column to a data frame if the specified column does not exist.
# If the column already exists, it can rename the column to
# lowercase (if `make_lowercase = TRUE`) and replace all `NA` values with a
# specified initial value (`init_value`) `all_na_replace = TRUE`

add_missing_column <- function(data, col_name, init_value, make_lowercase, all_na_replace = FALSE) {
  if (!tolower(col_name) %in% tolower(names(data))) {
    data |> tibble::add_column({{ col_name }} := init_value)
  } else {
    if (make_lowercase) data <- data |> dplyr::rename_with(tolower, dplyr::matches(col_name, ignore.case = TRUE))
    if (all_na_replace && all(is.na(data[[col_name]]))) data <- data |> mutate({{col_name}} := init_value)
    data
  }
}

#' Get Concentration Unit Based on Sample Amount Unit
#' internall analyte amount is pmol, thus when sample amount unit is uL,
#' pmol/uL equal to umol/L
#'
#' @param sample_amount_unit string with sample amount unit
#' @return string with feature_conc unit
#' @noRd

get_conc_unit <- function(sample_amount_unit) {
  units <- tolower(unique(sample_amount_unit))

  if (length(units) > 1) {
    conc_unit <- "pmol/sample amount unit (multiple units)"
  } else if (units == "ul" | units == "\U003BCl") {
    conc_unit <- "\U003BCmol/L"
  } else {
    conc_unit <- glue::glue("pmol/{sample_amount_unit}")
  }
  conc_unit
}


#
# # https://dewey.dunnington.ca/post/2018/modifying-facet-scales-in-ggplot2/
#
# FacetEqualWrap <- ggplot2::ggproto(
#   "FacetEqualWrap", FacetWrap,
#   train_scales = function(self, x_scales, y_scales, layout, data, params) {
#     # doesn't make sense if there is not an x *and* y scale
#     if (is.null(x_scales) || is.null(x_scales)) {
#       cli::cli_abort("X and Y scales required for facet_equal_wrap")
#     }
#
#     # regular training of scales
#     ggproto_parent(FacetWrap, self)$train_scales(x_scales, y_scales, layout, data, params)
#
#     # switched training of scales (x and y and y on x)
#     for (layer_data in data) {
#       match_id <- match(layer_data$PANEL, layout$PANEL)
#
#       x_vars <- intersect(x_scales[[1]]$aesthetics, names(layer_data))
#       y_vars <- intersect(y_scales[[1]]$aesthetics, names(layer_data))
#
#       SCALE_X <- layout$SCALE_X[match_id]
#       ggplot2:::scale_apply(layer_data, y_vars, "train", SCALE_X, x_scales)
#
#       SCALE_Y <- layout$SCALE_Y[match_id]
#       ggplot2:::scale_apply(layer_data, x_vars, "train", SCALE_Y, y_scales)
#     }
#   }
# )
#
# facet_wrap_equal <- function(...) {
#   # take advantage of the sanitizing that happens in facet_wrap
#   facet_super <- facet_wrap(...)
#
#   ggplot2::ggproto(NULL, FacetEqualWrap,
#     shrink = facet_super$shrink,
#     params = facet_super$params
#   )
# }

