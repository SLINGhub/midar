add_missing_column <- function(data, col_name, init_value, make_caps) {
  if (!tolower(col_name) %in% tolower(names(data)))
    data |> tibble::add_column({{col_name}}:= init_value)
  else
  {
    if(make_caps) data <- data |>  dplyr::rename_with(toupper, dplyr::matches(c(col_name), ignore.case = TRUE))
    data
  }
}
