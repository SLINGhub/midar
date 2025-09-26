#' Compare Feature Variability Before and After Normalization
#'
#' @description
#' Evaluates the effectiveness of normalization by comparing feature variability
#' (measured as %CV) in QC and/or study samples before and after normalization.
#' The comparison is visualized through one of three plot types:
#' * Scatter plot: CV values before vs after normalization
#' * Difference plot: (CV after - CV before) vs mean CV
#' * Ratio plot: log2 of (CV after / CV before) vs mean CV
#'
#' Features can be grouped and visualized by their fature class using facets.
#'
#' The resulting visualization helps assess whether normalization improved measurement
#' precision across different features and sample/QC types.
#'
#' @param data A `MidarExperiment` object
#' @param before_norm_var A string specifying the variable from the QC metrics
#'   table to be used for the x-axis (before normalization).
#' @param after_norm_var A string specifying the variable from the QC metrics
#'   table to be used for the y-axis (after normalization).
#' @param plot_type A character string specifying the type of plot to generate.
#' Must be one of "scatter", "diff", or "ratio". Selecting "scatter" plots the before and after normalization CV values
#' as a scatter plot, "diff" plots the difference between the two CV values against the average CV, and "ratio" plots the log2 ratio of the two CV values against the average CV.
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
#' @param x_lim Numeric vector of length 2 for x-axis limits. Use `NA` for
#'   auto-scaling (default is `c(0, NA)`).
#' @param y_lim Numeric vector of length 2 for y-axis limits. Use `NA` for
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
#'   plot_type = "scatter",
#'   qc_type = "SPL",
#'   filter_data = FALSE,
#'   facet_by_class = TRUE,
#'   cv_threshold_value = 25
#' )
#'
#' @export
plot_normalization_qc <- function(
  data = NULL,
  before_norm_var = c("intensity", "norm_intensity", "conc_raw"),
  after_norm_var = c("norm_intensity", "conc", "conc_raw"),
  plot_type = c("scatter", "diff", "ratio"),
  qc_types,
  facet_by_class = FALSE,
  filter_data = FALSE,
  include_qualifier = FALSE,
  cv_threshold_value = 25,
  x_lim = c(0, NA),
  y_lim = c(0, NA),
  cols_page = 5,
  point_size = 1,
  point_alpha = 0.5,
  font_base_size = 8
) {
  check_data(data)

  # Match and clean variable names
  before_norm_var <- str_remove(before_norm_var, "feature_")
  rlang::arg_match(
    before_norm_var,
    c("intensity", "norm_intensity", "conc_raw")
  )

  after_norm_var <- str_remove(after_norm_var, "feature_")
  rlang::arg_match(after_norm_var, c("norm_intensity", "conc", "conc_raw"))

  rlang::arg_match(plot_type, c("scatter", "diff", "ratio"))

  # Match qc_type
  #rlang::arg_match(qc_type, c("SPL", "BQC", "TQC", "NIST", "LTR", "PQC", "SST", "RQC"))

  if (all(is.na(qc_types))) {
    qc_types <- intersect(
      data$dataset$qc_type,
      c(
        "SPL",
        "TQC",
        "BQC",
        "TQC",
        "HQC",
        "MQC",
        "LQC",
        "NIST",
        "LTR",
        "SBLK",
        "PBLK",
        "IBLK",
        "QC"
      )
    )
  }

  if (!all(qc_types %in% unique(data$dataset$qc_type))) {
    cli::cli_abort(col_red(
      "One or more specified `qc_types` are not present in the dataset. Please verify data or analysis metadata."
    ))
  }
  start_regex <- paste(c(before_norm_var, after_norm_var), collapse = "|")
  middle_string <- "_cv_"
  end_regex <- paste(qc_types, collapse = "|")

  #col_pattern <- paste0("^(", start_regex, ")", middle_string, "(", end_regex, ")$")

  # Generate variable names for CV metrics
  x_variable <- stringr::str_c(before_norm_var, "_cv")
  y_variable <- stringr::str_c(after_norm_var, "_cv")

  if (x_variable == y_variable) {
    cli::cli_abort(col_red(
      "`before_norm_var` and `after_norm_var` cannot be the same."
    ))
  }

  # Check if QC metrics table is available
  if (nrow(data@metrics_qc) == 0) {
    cli::cli_abort(col_red(
      "No QC metrics available yet. Please run `calc_qc_metrics()`, or apply QC filters first."
    ))
  }

  if (!data@is_istd_normalized) {
    cli::cli_abort(col_red(
      "Data has not yet been normalized. Please run `normalize_by_istd()` first."
    ))
  } else {
    if (!data@is_quantitated & after_norm_var == "conc") {
      cli::cli_abort(col_red(
        "Data has not yet been quantitated. Please run one of the`quantitate_...()` functions first."
      ))
    }
  }

  # Call plot_qcmetrics_comparison to generate the plot
  plot_qcmetrics_comparison(
    data = data,
    plot_type = plot_type,
    facet_type = "wrap",
    filter_data = filter_data,
    facet_by_class = facet_by_class,
    x_variable = x_variable,
    y_variable = y_variable,
    qc_types = qc_types,
    threshold_value = cv_threshold_value,
    equality_line = TRUE,
    include_qualifier = include_qualifier,
    x_lim = x_lim,
    y_lim = y_lim,
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
#' The comparison is visualized through one of three plot types:
#' * Scatter plot: Values of `y_variable` vs `x_variable`
#' * Difference plot: (`y_variable` - `x_variable``) vs mean of both values
#' * Ratio plot: log2(`y_variable` / `x_variable`) vs mean of both values
#'
#'
#' @param data A `MidarExperiment` object containing pre-calculated QC metrics.
#' @param plot_type A character string specifying the type of plot to generate.
#' Must be one of "scatter", "diff", or "ratio". Selecting "scatter" plots the "y_variable" against the "y_variable" values
#' as a scatter plot, "diff" plots the difference between the two values against the average value, and "ratio" plots the log2 ratio of the two values against the average value.
#' @param facet_type A character string specifying the type of faceting to apply. "grid" creates a grid of facets with shared axes, while "wrap" creates a wrapped layout of facets, each with its own axes. Default is "wrap".
#' @param x_variable The name of the QC metric variable to be plotted on the
#'   x-axis.
#' @param y_variable The name of the QC metric variable to be plotted on the
#'   y-axis.
#' @param qc_types A character vector specifying the QC types to plot.
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
#' @param log_scale Logical, whether to use a log10 scale the axes.
#' @param x_lim Numeric vector of length 2 for x-axis limits. Use `NA` for
#'   auto-scaling (default is `c(0, NA)`).
#' @param y_lim Numeric vector of length 2 for y-axis limits. Use `NA` for
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
plot_qcmetrics_comparison <- function(
  data = NULL,
  plot_type = c("scatter", "diff", "ratio"),
  facet_type = c("grid", "wrap"),
  x_variable,
  y_variable,
  qc_types = NA,
  facet_by_class = FALSE,
  filter_data = FALSE,
  include_qualifier = TRUE,
  equality_line = FALSE,
  threshold_value = NA,
  log_scale = FALSE,
  x_lim = c(0, NA),
  y_lim = c(0, NA),
  cols_page = 5,
  point_size = 1,
  point_alpha = 0.5,
  font_base_size = 8
) {
  # Check if data is provided and valid
  check_data(data)

  rlang::arg_match(plot_type, c("scatter", "diff", "ratio"))
  rlang::arg_match(facet_type, c("grid", "wrap"))

  # Check if QC metrics table is available
  if (nrow(data@metrics_qc) == 0) {
    cli::cli_abort(col_red(
      "No QC metrics available yet. Please run
                           `calc_qc_metrics()`, or apply QC filters first."
    ))
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
  col_pattern <- paste0(x_variable, "|", y_variable)
  # Filter data based on valid features

  d_qc <- d_qc |>
    filter(.data$valid_feature) |>
    select("feature_id", "feature_class", dplyr::matches(col_pattern))

  var_match <- str_remove(x_variable, "_(tqc|bqc|spl)$") ==
    str_remove(y_variable, "_(tqc|bqc|spl)$")

  if (!var_match) {
    d_qc <- d_qc |>
      tidyr::pivot_longer(
        cols = dplyr::matches(col_pattern),
        names_to = c(".value", "qc_type"),
        names_pattern = "^(.*)_(.*)$"
      )
  } else {
    d_qc$qc_type <- "spl"
  }

  d_qc <- d_qc |>
    mutate(across(c(!!x_variable, !!y_variable), ~ ifelse(.x == 0, NA, .x))) |>
    drop_na()

  if (all(!is.na(qc_types))) {
    d_qc <- d_qc |>
      filter(.data$qc_type %in% tolower(qc_types))
  }

  if (plot_type == "diff") {
    d_qc <- d_qc |>
      mutate(
        y_values = (!!rlang::sym(y_variable) - !!rlang::sym(x_variable)),
        x_values = (!!rlang::sym(x_variable) + !!rlang::sym(y_variable)) / 2
      )
  } else if (plot_type == "ratio") {
    d_qc <- d_qc |>
      mutate(
        y_values = log2(!!rlang::sym(y_variable) / !!rlang::sym(x_variable)),
        x_values = (!!rlang::sym(x_variable) + !!rlang::sym(y_variable)) / 2
      )
  }

  # Apply additional QC filtering if requested
  if (filter_data) {
    if (data@is_filtered) {
      d_qc <- d_qc |> filter(.data$all_filter_pass)
    } else {
      cli_abort(cli::col_red(
        "Data has not yet been QC-filtered. Apply filter, or set `filter_data = FALSE`."
      ))
    }
  }

  x_label <- case_when(
    plot_type == "scatter" ~ paste("QC metric:", x_variable),
    plot_type == "diff" ~ paste("Mean of", x_variable, "and", y_variable),
    plot_type == "ratio" ~ paste("Mean of", x_variable, "and", y_variable),
  )

  y_label <- case_when(
    plot_type == "scatter" ~ paste("QC metric:", y_variable),
    plot_type == "diff" ~ paste(y_variable, "-", x_variable),
    plot_type == "ratio" ~ paste("log2(", y_variable, "/", x_variable, ")"),
  )

  # Check if feature_class exists and is valid
  if (
    facet_by_class &&
      (!"feature_class" %in% names(d_qc) || all(is.na(d_qc$feature_class)))
  ) {
    cli::cli_abort(col_red(
      "`feature_class` to be defined in the metadata. Please ammend metadata or set `facet_by_class = FALSE`."
    ))
  }

  d_qc$qc_type <- toupper(d_qc$qc_type)

  if (plot_type == "scatter") {
    g <- ggplot(
      data = d_qc,
      aes(
        x = !!rlang::sym(x_variable),
        y = !!rlang::sym(y_variable),
        color = .data$qc_type,
        fill = .data$qc_type,
        shape = .data$qc_type
      )
    )
  } else {
    g <- ggplot(
      data = d_qc,
      aes(
        x = .data$x_values,
        y = .data$y_values,
        color = .data$qc_type,
        fill = .data$qc_type,
        shape = .data$qc_type
      )
    ) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "red")
  }

  # Apply faceting if requested
  if (facet_by_class) {
    if (facet_type == "wrap") {
      if (is.na(y_lim[1])) {
        scalemode <- "free_y"
      } else {
        scalemode <- "fixed"
      }
      g <- g +
        facet_wrap(
          vars(.data$feature_class),
          scales = scalemode,
          ncol = cols_page
        )
    } else if (facet_type == "grid") {
      g <- g +
        facet_wrap(vars(.data$feature_class), scales = "free", ncol = cols_page)
    } else {
      cli::cli_abort(col_red("Invalid `facet_type`. Use 'stacked' or 'grid'."))
    }
  }

  # Plot points
  g <- g +
    geom_point(
      size = point_size,
      alpha = point_alpha,
      stroke = 0.3,
      na.rm = TRUE
    )

  # Plot equality line if specified
  if (equality_line && plot_type == "scatter") {
    # Create a new column to hold the maximum x and y values per feature_class
    # used for hidden points that ensure both axis have the same scale
    d_qc <- d_qc |>
      group_by(.data$feature_class) |>
      mutate(
        xy_max = max(
          !!rlang::sym(x_variable),
          !!rlang::sym(y_variable),
          na.rm = TRUE
        )
      ) |>
      ungroup()
    g <- g +
      geom_abline(intercept = 0, slope = 1, linewidth = 0.3, color = "orange")
  }
  # Add threshold lines if specified
  if (!is.na(threshold_value)) {
    g <- g +
      geom_vline(
        xintercept = threshold_value,
        linetype = "dashed",
        color = "darkgreen"
      )
    if(plot_type == "scatter"){
        g <- g +
        geom_hline(
          yintercept = threshold_value,
          linetype = "dashed",
          color = "darkgreen"
        )
    }

  }


  # Set axis limits
  g <- g +

    ggplot2::scale_color_manual(
      values = (pkg.env$qc_type_annotation$qc_type_col),
      drop = TRUE
    ) +
    ggplot2::scale_fill_manual(
      values = (pkg.env$qc_type_annotation$qc_type_fillcol),
      drop = TRUE
    ) +
    ggplot2::scale_shape_manual(
      values = (pkg.env$qc_type_annotation$qc_type_shape),
      drop = TRUE
    )

  if (log_scale) {
    g <- g +
      ggplot2::scale_x_log10(
        #labels = scientific_format_end,
        #limits = c(0.1, NA),
        name = x_label,
        limits = x_lim,
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::scale_y_log10(
        #labels = scientific_format_end,
        limits = y_lim,
        name = y_label,
        expand = ggplot2::expansion(mult = c(0, 0.02))
      )
  } else {
    g <- g +
      ggplot2::scale_x_continuous(
        #labels = scientific_format_end,
        name = x_label,
        limits = x_lim,
        expand = ggplot2::expansion(mult = c(0, 0.02))
      ) +
      ggplot2::scale_y_continuous(
        #labels = scientific_format_end,
        limits = y_lim,
        name = y_label,
        expand = ggplot2::expansion(mult = c(0, 0.02))
      )
  }

  # Customize plot theme
  g <- g +
    theme_bw(base_size = font_base_size) +
    theme(
      plot.title = element_text(size = font_base_size, face = "bold"),
      strip.text = ggplot2::element_text(size = font_base_size, face = "bold"),
      axis.text = element_text(size = font_base_size),
      axis.title = element_text(size = font_base_size, face = "bold"),
      panel.grid = element_line(linewidth = 0.001),
      strip.background = ggplot2::element_rect(
        linewidth = 0.0001,
        fill = "#00283d"
      ),
      strip.text.x = ggplot2::element_text(color = "white"),
      #strip.switch.pad.wrap = ggplot2::unit(1, "mm"),
      panel.border = element_rect(linewidth = 0.5, color = "grey40"),
      legend.position = "right"
    )

  return(g)
}
