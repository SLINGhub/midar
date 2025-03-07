
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


#The compare_values function compares a column in a data frame to a threshold using a
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

compare_values <- function(tbl, val, threshold, operator) {
  if(nrow(tbl) == 0) stop("tbl has no rows")
  if (!val %in% names(tbl)  && !is.na(threshold)) {
    var_name <- deparse(substitute(threshold))
    error_message <- case_when(
      str_detect(var_name, "conc") ~ "Cannot filter by `{var_name}` because concentration data is unavailable. Please quantify the data first using `quantify_by_*()` functions.",
      str_detect(var_name, "normint") ~ "Cannot filter by `{var_name}` because normalized data is unavailable. Please normalize the data first using `normalize_by_*()` functions.",
      str_detect(var_name, "response") ~ "Cannot filter by `{var_name}` because response curve data is unavailable. Please verify the corresponding data and metadata.",
      TRUE ~ "QC parameter is not available. Please verify the argument `val`."
    )
    cli_abort(col_red(error_message))
   }
  if (all(is.na(tbl[[val]])) && !is.na(threshold)) {
    var_name <- deparse(substitute(threshold))
    cli_abort("Underlying QC parameter to filter for {var_name} is not available. Please verify data if selected QC type is present and contains results.")
  }
  v_val <- tbl[[val]]
  if(is.null(v_val) && is.na(threshold)) return(rep(NA, nrow(tbl)))
  result <- get(operator)(v_val, threshold)
}


# performs element-wise logical operations (AND or OR) across multiple
# logical vectors in a list. It returns a vector of the results,
# where each element is the result of applying the specified operation to the
# corresponding elements of the input vectors.
# If any element is NA, the result for that position will also be NA.
# If list is empty NULL is returned
comp_lgl_vec <- function(lgl_list, .operator) {

  # Convert the list of vectors into a matrix
  matrix_data <- do.call(cbind, lgl_list)

  # Define the operation function based on the operator
  op_func <- switch(.operator,
                    "AND" = function(x) if (all(is.na(x))) NA else all(x, na.rm = TRUE),
                    "OR"  = function(x) if (all(is.na(x))) NA else any(x, na.rm = TRUE),
                    "XOR" = function(x) if (all(is.na(x))) NA else Reduce(xor, x[!is.na(x)]),
                    stop("Unsupported operator"))

  # Apply the operation function across the columns
  result <- apply(matrix_data, 1, op_func)

  return(result)
}



# comp_lgl_vec <- function(lgl_list, .operator){
#   browser()
#   if (.operator == "AND") {
#     return(Reduce("&", lgl_list))
#   } else if (.operator == "OR") {
#     return(Reduce("|", lgl_list))
#   } else if (.operator == "XOR") {
#     return(Reduce(function(x, y) xor(x, y), lgl_list))
#   } else {
#     # Return NULL for unsupported operators
#     return(NULL)
#   }
# }

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
#' @param analyte_amount_unit string with base analyre amount unit
#' @return string with feature_conc unit
#' @noRd

get_conc_unit <- function(sample_amount_unit, analyte_amount_unit) {
  units <- tolower(unique(sample_amount_unit))
  analyte_units <- tolower(unique(analyte_amount_unit))

  if (length(units) > 1) {
    conc_unit <- glue::glue("{analyte_amount_unit}/sample amount unit (multiple units)")
  } else if (analyte_amount_unit == "pmol"  && (units == "ul" | units == "\U003BCl")) {
    conc_unit <- "\U003BCmol/L"
  } else if (!str_detect(analyte_units, "\\/") && !str_detect(analyte_units, "\\-1")) {
    conc_unit <- glue::glue("{analyte_amount_unit}/{sample_amount_unit}")
  } else {
    conc_unit <- analyte_amount_unit
  }
  unique(conc_unit)
}

#' Reorder Data Frame based on a chain of linked values in two columns.
#'
#' This function orders rows of a data frame based on chained relationships defined by two columns.
#' It can also handle fully disconnected rows (i.e., rows where both `From` and `To` values
#' are not present in other rows). The behavior for disconnected rows is controlled via the
#' `disconnected_action` parameter.
#'
#' @param df A data frame containing the chain relationships.
#' @param from_col A string specifying the column name representing the starting point of the chain.
#' @param to_col A string specifying the column name representing the endpoint of the chain.
#' @param include_chain_id A logical indicating whether to include a `chain_id` column in the output.
#' @param disconnected_action A string indicating how to handle fully disconnected rows. Options are:
#'   \describe{
#'     \item{"exclude"}{Exclude disconnected rows from the output.}
#'     \item{"keep"}{Keep disconnected rows in the result.}
#'   }
#'
#' @return A data frame containing ordered chains with a `chain_id` column to distinguish between different chains.
#'   If disconnected rows are included, they will have their own `chain_id`.
#'
#' @examples
#' df_unordered <- data.frame(
#'   From = c("INSPECT", "VERIFY", "START", "NULL", "NEW", "CREATE", "MID", "DIFFERENT", "OUTLIER"),
#'   To = c("VERIFY", "PUBLISH", "MID", "NEW", "CREATE", "INSPECT", "END", "NOTSAME", "INSIDER"),
#' stringsAsFactors = FALSE
#' )
#'
#' # Order keeping disconnected rows
#' order_chained_columns_tbl(df_unordered, "From", "To", FALSE, "keep")
#'
#' # Ordr excluding disconnected rows
#' order_chained_columns_tbl(df_unordered, "From", "To", FALSE, "exclude")
#'
#'
#' @export
order_chained_columns_tbl <- function(df, from_col, to_col, include_chain_id, disconnected_action = "keep") {
  # Match the argument for disconnected_action
  disconnected_action <- match.arg(disconnected_action, c("keep", "exclude"))


  if(nrow(df) == 0) stop("Data frame has no rows")
  if(!all(c(from_col, to_col) %in% colnames(df))) stop("One or more columns are not present in the data frame.")

  # Step 1: Identify connected nodes (rows that are involved in a chain)
  df_initial <- df
  all_from <- unique(df[[from_col]])
  all_to <- unique(df[[to_col]])

  # Find rows where From or To is not in any other row's From or To
  unconnected_rows <- df[!(df[[from_col]] %in% all_to) & !(df[[to_col]] %in% all_from), ]

  # Step 2: Exclude unconnected rows if disconnected_action is "exclude"
  if (disconnected_action == "exclude") {
    df <- df[!(df[[from_col]] %in% unconnected_rows[[from_col]] | df[[to_col]] %in% unconnected_rows[[to_col]]), ]
  }

  # Step 3: Identify connected rows (after filtering disconnected if needed)
  all_from <- unique(df[[from_col]])
  all_to <- unique(df[[to_col]])

  # Step 4: Filter out disconnected rows based on the current df
  connected_df <- df[df[[from_col]] %in% all_from | df[[to_col]] %in% all_to, ]

  # Step 5: Build the chains from the connected rows only
  chain_map <- setNames(connected_df[[to_col]], connected_df[[from_col]])

  # Initialize ordered chains and visited set
  ordered_chains <- list()
  visited <- character(0)

  # Start the chain from any From node that is not a To node
  starts <- setdiff(connected_df[[from_col]], connected_df[[to_col]])

  if(length(starts) == 0) stop("Circular dependency detected. Please verify the input data.")
  # Traverse each chain from the starting node
  for (start in starts) {
    chain <- c(start)
    current <- start

    # Follow the chain from From -> To
    while (!is.null(chain_map[current]) && !is.na(chain_map[current])) {
      next_value <- chain_map[current]

      # Check for circular dependencies
      if (next_value %in% visited) {
        stop("Circular dependency detected in the chain.")
      }

      chain <- c(chain, next_value)
      visited <- union(visited, next_value)  # Mark the next node as visited
      current <- next_value
    }

    # Add the chain to the list
    ordered_chains[[length(ordered_chains) + 1]] <- chain
  }

  # Step 6: Convert chains into a data frame
  connected_chains_df <- do.call(rbind, lapply(seq_along(ordered_chains), function(i) {
    chain <- ordered_chains[[i]]
    data.frame(
      chain_id = i,
      From = chain[-length(chain)],  # All except last
      To = chain[-1],  # All except first
      stringsAsFactors = FALSE
    )
  }))

  names(connected_chains_df) <- c("chain_id", from_col, to_col)

  if(nrow(connected_chains_df) < nrow(df)) stop("Circular dependency detected. Please verify the input data.")

  # Cleanup and readd columns that were not part of the chain
  rownames(connected_chains_df) <- NULL


  cols_to_add <- setdiff(names(df_initial), c("chain_id", from_col, to_col))
  connected_chains_df$order <- seq_len(nrow(connected_chains_df))
  merged_df <- merge(connected_chains_df, df_initial[, c(from_col, to_col, cols_to_add)], by = c(from_col, to_col), all.x = TRUE)
  merged_df <- merged_df[order(merged_df$order), ]
  merged_df$order <- NULL

  # Remove chain_id if not specified
  if(!include_chain_id) merged_df$chain_id <- NULL

  # Return final data frame
  merged_df
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

