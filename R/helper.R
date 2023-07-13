add_missing_column <- function(data, col_name, init_value, make_caps) {
  if (!tolower(col_name) %in% tolower(names(data)))
    data |> tibble::add_column({{col_name}}:= init_value)
  else
  {
    if(make_caps) data <- data |>  dplyr::rename_with(toupper, dplyr::matches(c(col_name), ignore.case = TRUE))
    data
  }
}



# https://dewey.dunnington.ca/post/2018/modifying-facet-scales-in-ggplot2/

FacetEqualWrap <- ggplot2::ggproto(
  "FacetEqualWrap", FacetWrap,

  train_scales = function(self, x_scales, y_scales, layout, data, params) {

    # doesn't make sense if there is not an x *and* y scale
    if (is.null(x_scales) || is.null(x_scales)) {
      stop("X and Y scales required for facet_equal_wrap")
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
