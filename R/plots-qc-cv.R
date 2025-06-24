#' Compare %CV values before and after normalization
#'
#' This function compares the coefficient of variation (CV) of QC or study
#' samples before and after normalization. It preselects the relevant QC metrics
#' based on the chosen arguments and visualizes the comparison through a scatter
#' plot. The plot can be faceted by `feature_class`.
#'
#' @param data A `MidarExperiment` object
#' @param before_norm_var A string specifying the variable from the QC metrics
#'   table to be used for the x-axis (before normalization).
#' @param after_norm_var A string specifying the variable from the QC metrics
#'   table to be used for the y-axis (after normalization).
#' @param qc_types A character vector specifying the QC types to plot. It
#' must contain at least one element. The default is `NA`, which means any
#' of the non-blank QC types ("SPL", "TQC", "BQC", "HQC", "MQC", "LQC",
#' "NIST", "LTR") will be plotted if present in the dataset.
#' @param facet_by_class If `TRUE`, facets the plot by `feature_class`, as defined
#'   in the feature metadata.
#' @param filter_data Whether to use all data (default) or only
#'   QC-filtered data (filtered via [filter_features_qc()]).
#' @param include_qualifier Whether to include qualifier features
#'   (default is `TRUE`).
#' @param cv_threshold_value Numerical threshold value to be shown as dashed
#'   lines in the plot (default is `25`).
#' @param xlim Numeric vector of length 2 for x-axis limits. Use `NA` for
#'   auto-scaling (default is `c(0, NA)`).
#' @param ylim Numeric vector of length 2 for y-axis limits. Use `NA` for
#'   auto-scaling (default is `c(0, NA)`).
#' @param cols_page Number of facet columns per page, representing
#'   different feature classes (default is `5`). Only used if
#'   `facet_by_class = TRUE`.
#' @param point_size Size of points in millimeters (default is `1`).
#' @param point_alpha Transparency of points (default is `0.5`).
#' @param font_base_size Base font size in points (default is `8`).
#'
#' @return A `ggplot2` object representing the scatter plot comparing CV values
#'   before and after normalization.
#'
#' @details
#' The function preselects the corresponding variables from the QC metrics and uses
#' [plot_qcmetrics_comparison()] to visualize the results.
#'
#' - The data must be normalized before using [normalize_by_istd()] followed by
#' calculation of the QC metrics table via [calc_qc_metrics()] or
#' [filter_features_qc()], see examples below.
#'
#' - When `facet_by_class = TRUE`, then the `feature_class` must be defined in the
#'   metadata or retrieved via specific functions, e.g., [parse_lipid_feature_names()].
#'
#' @seealso
#' [plot_qcmetrics_comparison()], [calc_qc_metrics()], [filter_features_qc()],
#'
#' @examples
#' # Example usage:
#' mexp <- lipidomics_dataset
#' mexp <- normalize_by_istd(mexp)
#' mexp <- calc_qc_metrics(mexp)
#' plot_normalization_qc(
#'   data = mexp,
#'   before_norm_var = "intensity",
#'   after_norm_var = "norm_intensity",
#'   qc_type = "SPL",
#'   filter_data = FALSE,
#'   facet_by_class = TRUE,
#'   cv_threshold_value = 25
#' )
#'
#' @export
plot_normalization_qc <- function(data = NULL,
                                  before_norm_var = c("intensity", "norm_intensity", "conc_raw"),
                                  after_norm_var = c("norm_intensity", "conc", "conc_raw"),
                                  qc_types,
                                  facet_by_class = FALSE,
                                  filter_data = FALSE,
                                  include_qualifier = FALSE,
                                  cv_threshold_value = 25,
                                  xlim = c(0, NA),
                                  ylim = c(0, NA),
                                  cols_page = 5,
                                  point_size = 1,
                                  point_alpha = 0.5,
                                  font_base_size = 8) {

  check_data(data)

  # Match and clean variable names
  before_norm_var <- str_remove(before_norm_var, "feature_")
  rlang::arg_match(before_norm_var, c("intensity", "norm_intensity", "conc_raw"))

  after_norm_var <- str_remove(after_norm_var, "feature_")
  rlang::arg_match(after_norm_var, c("norm_intensity", "conc", "conc_raw"))

  # Match qc_type
  #rlang::arg_match(qc_type, c("SPL", "BQC", "TQC", "NIST", "LTR", "PQC", "SST", "RQC"))


  if(all(is.na(qc_types))){
    qc_types <- intersect(data$dataset$qc_type, c("SPL", "TQC", "BQC", "TQC", "HQC", "MQC", "LQC", "NIST", "LTR", "SBLK", "PBLK", "IBLK", "QC"))
  }

  if (!all(qc_types %in% unique(data$dataset$qc_type))) {
    cli::cli_abort(col_red("One or more specified `qc_types` are not present in the dataset. Please verify data or analysis metadata."))
  }

  start_regex <- paste(c(before_norm_var, after_norm_var), collapse = "|")
  middle_string <- "_cv_"
  end_regex <- paste(qc_types, collapse = "|")


  col_pattern <- paste0("^(", start_regex, ")", middle_string, "(", end_regex, ")$")



  # Generate variable names for CV metrics
  x_variable <- stringr::str_c(before_norm_var, "_cv")
  y_variable <- stringr::str_c(after_norm_var, "_cv")

 if (x_variable == y_variable) {
  cli::cli_abort(col_red("`before_norm_var` and `after_norm_var` cannot be the same."))
 }

  # Check if QC metrics table is available
  if (nrow(data@metrics_qc) == 0) {
    cli::cli_abort(col_red("No QC metrics available yet. Please run `calc_qc_metrics()`, or apply QC filters first."))
  }

  if(!data@is_istd_normalized) {
    cli::cli_abort(col_red("Data has not yet been normalized. Please run `normalize_by_istd()` first."))
  } else {
    if(!data@is_quantitated & after_norm_var == "conc") {
      cli::cli_abort(col_red("Data has not yet been quantitated. Please run one of the`quantitate_...()` functions first."))
    }
  }

   # Call plot_qcmetrics_comparison to generate the plot
  plot_qcmetrics_comparison(
    data = data,
    col_pattern = col_pattern,
    filter_data = filter_data,
    facet_by_class = facet_by_class,
    x_variable = x_variable,
    y_variable = y_variable,
    threshold_value = cv_threshold_value,
    equality_line = TRUE,
    include_qualifier = include_qualifier,
    xlim = xlim,
    ylim = ylim,
    cols_page = cols_page,
    point_size = point_size,
    point_alpha = point_alpha,
    font_base_size = font_base_size
  )
}


#' Comparison of two feature QC metrics variables
#'
#' This function generates scatter plots comparing two QC metrics variables
#' across feature classes. A list of available QC metrics is available from the
#' [calc_qc_metrics()] documentation.
#'
#' @param data A `MidarExperiment` object containing pre-calculated QC metrics.
#' @param x_variable The name of the QC metric variable to be plotted on the
#'   x-axis.
#' @param y_variable The name of the QC metric variable to be plotted on the
#'   y-axis.
#' @param col_pattern A string pattern to match the columns in the QC metrics
#' @param facet_by_class If `TRUE`, facets the plot by `feature_class`, as defined
#'   in the feature metadata.
#' @param filter_data Logical; whether to use all data (default) or only
#'   QC-filtered data (filtered via [filter_features_qc()]).
#' @param include_qualifier Logical; whether to include qualifier features
#'   (default is `TRUE`).
#' @param equality_line Logical; whether to show a line indicating
#'   identical values in both compared variables (default is `FALSE`).
#' @param threshold_value Numeric; threshold value to be shown as dashed lines
#'   from both axes on the plot (default is `NA`).
#' @param xlim Numeric vector of length 2 for x-axis limits. Use `NA` for
#'   auto-scaling (default is `c(0, NA)`).
#' @param ylim Numeric vector of length 2 for y-axis limits. Use `NA` for
#'   auto-scaling (default is `c(0, NA)`).
#' @param cols_page Integer; number of facet columns per page (default is `5`).
#' @param point_size Numeric; size of points in millimeters (default is `1`).
#' @param point_alpha Numeric; transparency of points (default is `0.5`).
#' @param font_base_size Numeric; base font size in points (default is `8`).
#'
#' @return A `ggplot2` object representing the scatter plot.
#'
#' @details
#'
#' - `x_variable` and `y_variable` must be available in the QC metrics table.
#' Please refer to the help page of [calc_qc_metrics()] for more information on
#' the available QC metric variables.
#'
#' - When `facet_by_class = TRUE`, then the `feature_class` must be defined in the
#'   metadata or retrieved via specific functions, e.g., [parse_lipid_feature_names()].
#'
#' @seealso
#' [calc_qc_metrics()], [filter_features_qc()], [plot_normalization_qc()], [normalize_by_istd()]
#'
#' @export
plot_qcmetrics_comparison <- function(data = NULL,
                                   x_variable,
                                   y_variable,
                                   col_pattern,
                                   facet_by_class = FALSE,
                                   filter_data = FALSE,
                                   include_qualifier = TRUE,
                                   equality_line = FALSE,
                                   threshold_value = NA,
                                   xlim = c(0, NA),
                                   ylim = c(0, NA),
                                   cols_page = 5,
                                   point_size = 1,
                                   point_alpha = 0.5,
                                   font_base_size = 8) {

  # Check if data is provided and valid
  check_data(data)

  # Check if QC metrics table is available
  if (nrow(data@metrics_qc) == 0) {
    cli::cli_abort(col_red("No QC metrics available yet. Please run
                           `calc_qc_metrics()`, or apply QC filters first."))
  }

  # if(!all(c(x_variable, y_variable) %in% colnames(data@metrics_qc))) {
  #   cli::cli_abort(col_red("One or both of the specified variables (`x_variable`
  #                          and `y_variable`) are not present in the QC metrics
  #                          table. Please verify the help page"))
  # }


  d_qc <- data@metrics_qc

  # Include or exclude qualifier features
  if (!include_qualifier) {
    d_qc <- d_qc |> filter(.data$is_quantifier, .data$valid_feature)
  }

  # Filter data based on valid features
  d_qc <- d_qc|>
    filter(.data$valid_feature) |>
    select("feature_id", "feature_class", dplyr::matches(col_pattern)) |>
    tidyr::pivot_longer(
      cols = dplyr::matches(col_pattern),
      names_to = c(".value", "qc_type"),
      names_pattern = "^(.*)_(.*)$"
    ) |>
    mutate(across(c(!!x_variable, !!y_variable), ~ ifelse(.x == 0, NA, .x))) |>
    drop_na()

  # Apply additional QC filtering if requested
  if (filter_data) {
    if (data@is_filtered) {
      d_qc <- d_qc |> filter(.data$valid_feature, .data$all_filter_pass)
    } else {
      cli_abort(cli::col_red("Data has not yet been QC-filtered. Apply filter, or set `filter_data = FALSE`."))
    }
  }



  # Check if feature_class exists and is valid
  if (facet_by_class && (!"feature_class" %in% names(d_qc) || all(is.na(d_qc$feature_class)))) {
    cli::cli_abort(col_red("`feature_class` to be defined in the metadata. Please ammend metadata or set `facet_by_class = FALSE`."))
  }



  # Begin ggplot object
  g <- ggplot(data = d_qc, aes(x = !!rlang::sym(x_variable),
                               y = !!rlang::sym(y_variable),
                               color = .data$qc_type,
                               fill = .data$qc_type,
                               shape = .data$qc_type))
  #g <- ggplot(data = d_qc, aes(x = x_var, y = y_var, color = .data$qc_type, fill = .data$qc_type))
  # Apply faceting if requested
  if (facet_by_class) {
    g <- g + facet_wrap(vars(.data$feature_class), scales = "free", ncol  = cols_page)
  }

  # Plot points
  g <- g + geom_point(size = point_size, alpha = point_alpha, stroke = 0.3, na.rm = TRUE)

  # Plot equality line if specified
  if (equality_line) {
    # Create a new column to hold the maximum x and y values per feature_class
    # used for hidden points that ensure both axis have the same scale
    d_qc <- d_qc |>
      group_by(.data$feature_class) |>
      mutate(xy_max = max(!!rlang::sym(x_variable), !!rlang::sym(y_variable), na.rm = TRUE)) |>
      ungroup()
    g <- g + geom_abline(intercept = 0, slope = 1, linewidth = 0.3, color = "orange")
  }

  # Add threshold lines if specified
  if (!is.na(threshold_value)) {
    g <- g +
      geom_hline(yintercept = threshold_value, linetype = "dashed", color = "darkgreen") +
      geom_vline(xintercept = threshold_value, linetype = "dashed", color = "darkgreen")
  }

  # Set axis limits
  g <- g +
    scale_x_continuous(limits = xlim, expand = expansion(mult = c(0, 0.2))) +
    scale_y_continuous(limits = ylim, expand = expansion(mult = c(0, 0.2))) +
    ggplot2::scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
    ggplot2::scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
    ggplot2::scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)

  # Customize plot theme
  g <- g + theme_bw(base_size = font_base_size) +
    theme(
      plot.title = element_text(size = font_base_size, face = "bold"),
      strip.text = ggplot2::element_text(size = font_base_size, face = "bold"),
      axis.text = element_text(size = font_base_size),
      axis.title = element_text(size = font_base_size, face = "bold"),
      panel.grid = element_line(linewidth = 0.001),
      strip.background = ggplot2::element_rect(linewidth = 0.0001, fill = "#00283d"),
      strip.text.x = ggplot2::element_text(color = "white"),
      #strip.switch.pad.wrap = ggplot2::unit(1, "mm"),
      panel.border = element_rect(linewidth = 0.5, color = "grey40"),
      legend.position = "right"
    )

  return(g)
}

