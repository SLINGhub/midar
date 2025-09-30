#' Plot Abundance Profile
#'
#' Creates a profile plot showing the abundance distribution of features across different classes.
#'
#' @param data A MidarExperiment object.
#' @param variable A character string indicating the variable to plot. For `use_qc_metrics = FALSE`,
#'   this must be a base name like "area" or "conc". For `use_qc_metrics = TRUE`, this is the
#'   base name of a metric in the `metrics_qc` table (e.g., "rt" for "rt_mean_SPL").
#' @param qc_types A character vector specifying the QC types to be averaged and plotted.
#'   If `use_qc_metrics` is `TRUE`, this must be a single character string (e.g., "SPL").
#' @param log_scale A logical value indicating whether to use a log10 scale for the x-axis.
#' @param use_qc_metrics A logical value. If `FALSE` (default), data is summarized on the fly from the main dataset.
#'   If `TRUE`, pre-calculated summary data is used from the `metrics_qc` table, which is much faster.
#'   When `TRUE`, `qc_types` must specify only one QC type.
#' @param show_sum A logical value indicating whether to plot a summary point (diamond shape)
#'   representing the mean of the summed abundances for each class. Defaults to `TRUE` for
#'   abundance-related variables and `FALSE` for others (e.g., "rt", "fwhm").
#' @param feature_map A named character vector specifying feature classes, their order, and color. The order of the vector determines the y-axis order.
#'   Alternatively, can be a single string: `lipidomics` to use a default theme comprising mapping for lipid classes frequently measured in lipidomics methods. When using
#'  the default theme, set `drop_empty_classes = TRUE` to avoid showing classes not present in the data.
#'   If `NA` (default), the default theme is used.
#' @param drop_empty_classes A logical value. If `TRUE` (default), feature classes from
#'   `feature_map` that are not found in the data are removed from the plot. If `FALSE`,
#'   they are kept as empty rows.
#' @param density_strip A logical value. If `TRUE`, a density strip showing the distribution of all
#'   features is added to the top of the plot.
#' @param exclude_classes A character vector of feature classes to be excluded from the plot.
#' @param filter_data A logical value indicating whether to use all data (`FALSE`, default)
#'   or only QC-filtered data (`TRUE`, via [filter_features_qc()]). This is ignored if `use_qc_metrics` is `TRUE`.
#' @param include_qualifier A logical value indicating whether to include qualifier features. Default is `FALSE`.
#' @param include_istd A logical value indicating whether to include internal standard (ISTD) features. Default is `FALSE`.
#' @param include_feature_filter A character or regex pattern to filter features by `feature_id`.
#'   If a vector is supplied, features matching any pattern will be included.
#' @param exclude_feature_filter A character or regex pattern to exclude features by `feature_id`.
#'   If a vector is supplied, features matching any pattern will be excluded.
#' @param analysis_range A numeric vector of length 2 specifying the inclusive range of `analysis_order`. If `NA` (default), all analyses are included.
#' @param scale_factor A numeric value to scale the selected variable. Default is 1.
#' @param x_lim A numeric vector of length 2 specifying the x-axis limits. If `log_scale` is `TRUE`,
#'   these are interpreted as exponents (e.g., `c(2, 6)` for `10^2` to `10^6`). If `NA`, limits are
#'   determined automatically from the data.
#' @param x_label A character string for the x-axis label. If `NA`, a default label is generated.
#' @param y_axis_position A character string specifying the y-axis position ("left" or "right"). Default is "right".
#' @param segment_width A numeric value for the linewidth of the feature segments. Default is 0.25.
#' @param font_base_size A numeric value for the base font size. Default is 8.
#' @param grid_major_color Color for major grid lines. Default is "grey50".
#' @param grid_major_linewidth Linewidth for major grid lines. Default is 0.1.
#' @param grid_minor_color Color for minor grid lines. Default is "grey50".
#' @param grid_minor_linewidth Linewidth for minor grid lines. Default is 0.1.
#'
#' @returns A `ggplot` object representing the abundance profile plot.
#'
#' @export
plot_abundanceprofile <- function(
  data,
  variable,
  qc_types,
  log_scale,
  use_qc_metrics = FALSE,
  feature_map = NA,
  show_sum = NA,
  exclude_classes = NA,
  filter_data = FALSE,
  include_qualifier = FALSE,
  include_istd = FALSE,
  include_feature_filter = NA,
  exclude_feature_filter = NA,
  analysis_range = NA,
  drop_empty_classes = TRUE,
  density_strip = FALSE,
  scale_factor = 1,
  x_lim = NA,
  x_label = NA,
  y_axis_position = "right",
  segment_width = 0.25,
  font_base_size = 8,
  grid_major_color = "grey50",
  grid_major_linewidth = 0.1,
  grid_minor_color = "grey50",
  grid_minor_linewidth = 0.1
) {
  # --- 1. Argument Validation and Setup ---
  variable_clean <- stringr::str_remove(variable, "feature_")
  rlang::arg_match(y_axis_position, c("left", "right"))

  ## FIX: The arg_match for `variable_clean` has been moved into the `else` block
  ## of the `use_qc_metrics` check to allow for flexible column names from the metrics table.

  if (is.na(show_sum)) {
    show_sum <- !(stringr::str_detect(
      variable_clean,
      "rt|fwhm|width|symmetry"
    ) ||
      (use_qc_metrics &&
        stringr::str_detect(
          variable_clean,
          "cv|rt|dratio|prop|signalblank|response"
        )))
  }

  ## Expanded logic to resolve feature_map from shortcuts, NA, or a custom vector.
  feature_map_resolved <- NULL
  if (
    is.character(feature_map) &&
      length(feature_map) == 1 &&
      is.null(names(feature_map))
  ) {
    # Case 1: A shortcut string like "lipidomics"
    if (feature_map == "lipidomics") {
      feature_map_resolved <- pkg.env$lipid_class_annotations$lipid_class_map
    } else if (feature_map == "metabolomics") {
      feature_map_resolved <- pkg.env$lipid_class_annotations$lipid_class_map
    }
  } else if (length(feature_map) == 1 && is.na(feature_map)) {
    # Case 2: Default NA, infer from data object
    feature_map_resolved <- pkg.env$lipid_class_annotations$lipid_class_map
  } else if (is.character(feature_map) && !is.null(names(feature_map))) {
    # Case 3: A user-provided named vector
    feature_map_resolved <- feature_map
  }

  # Final validation of the resolved map
  if (is.null(feature_map_resolved)) {
    stop(
      "Invalid `feature_map` argument. Must be NA, 'lipidomics', 'metabolomics', or a named character vector of colors."
    )
  }

  # --- 2. Data Source Selection ---
  if (use_qc_metrics) {
    # --- 2a. Use Pre-summarized QC Metrics ---
    if (length(qc_types) != 1) {
      warning("When `use_qc_metrics` is TRUE, `qc_types` will be ignored.")
    }

    if (!variable %in% names(data@metrics_qc)) {
      stop(paste(
        "Could not find a summary column in `metrics_qc` matching:",
        variable
      ))
    }

    d_features <- data@metrics_qc
    if (filter_data) {
      d_features <- d_features |> dplyr::filter(.data$all_filter_pass)
    }
    d_features <- d_features |>
      dplyr::select(
        "feature_id",
        "feature_class",
        "is_istd",
        "is_quantifier",
        abundance_mean = !!rlang::sym(variable)
      ) |>
      dplyr::mutate(abundance_mean = .data$abundance_mean * scale_factor)

    if (!include_istd) {
      d_features <- d_features |> dplyr::filter(!.data$is_istd)
    }
    if (!include_qualifier) {
      d_features <- d_features |> dplyr::filter(.data$is_quantifier)
    }
  } else {
    # --- 2b. Summarize Data from Raw Dataset (Original Logic) ---
    ## FIX: Validation is now performed here, only for this code path.
    valid_vars <- c(
      "area",
      "height",
      "intensity",
      "norm_intensity",
      "intensity_raw",
      "intensity_before",
      "norm_intensity_raw",
      "norm_intensity_before",
      "response",
      "conc",
      "conc_raw",
      "conc_before",
      "rt",
      "fwhm",
      "width",
      "symmetry"
    )
    rlang::arg_match(variable_clean, valid_vars)

    variable_with_prefix <- stringr::str_c("feature_", variable_clean)
    check_var_in_dataset(data@dataset, variable_with_prefix)

    d_filt <- get_dataset_subset(
      data,
      filter_data = filter_data,
      qc_types = qc_types,
      include_qualifier = include_qualifier,
      include_istd = include_istd,
      include_feature_filter = include_feature_filter,
      exclude_feature_filter = exclude_feature_filter
    ) |>
      dplyr::mutate(
        !!rlang::sym(variable_with_prefix) := !!rlang::sym(
          variable_with_prefix
        ) *
          scale_factor
      )
    if (!all(is.na(analysis_range))) {
      d_filt <- d_filt |>
        dplyr::filter(
          .data$analysis_order >= analysis_range[1] &
            .data$analysis_order <= analysis_range[2]
        )
    }

    d_features <- d_filt |>
      dplyr::group_by(.data$feature_id, .data$feature_class) |>
      dplyr::summarise(
        abundance_mean = mean(!!rlang::sym(variable_with_prefix), na.rm = TRUE),
        .groups = "drop"
      )
  }

  # --- 3. Common Data Preparation & Y-Axis Setup ---
  if (!is.na(exclude_classes)) {
    d_features <- d_features |>
      dplyr::filter(!(.data$feature_class %in% exclude_classes))
  }

  if (drop_empty_classes) {
    present_classes <- unique(d_features$feature_class)
    feature_map_plot <- feature_map_resolved[
      names(feature_map_resolved) %in% present_classes
    ]
  } else {
    feature_map_plot <- feature_map_resolved
  }
  feature_order_plot <- names(feature_map_plot)

  d_features$feature_class <- factor(
    d_features$feature_class,
    levels = rev(feature_order_plot)
  )
  d_features <- d_features |>
    mutate(
      abundance_mean = ifelse(
        is.infinite(.data$abundance_mean),
        NA,
        .data$abundance_mean
      )
    ) |>
    tidyr::drop_na("feature_class", "abundance_mean")

  d_summary <- d_features |>
    dplyr::group_by(.data$feature_class) |>
    dplyr::summarise(
      abundance_min = min(.data$abundance_mean, na.rm = TRUE),
      abundance_max = max(.data$abundance_mean, na.rm = TRUE),
      abundance_sum_mean = sum(.data$abundance_mean, na.rm = TRUE),
      .groups = "drop"
    )

  if (log_scale) {
    d_summary <- d_summary |>
      dplyr::mutate(
        abundance_min_padded = .data$abundance_min * 0.8,
        abundance_max_padded = .data$abundance_max * 1.2
      )
  } else {
    d_summary <- d_summary |>
      dplyr::mutate(
        abundance_min_padded = .data$abundance_min * 0.98,
        abundance_max_padded = .data$abundance_max * 1.02
      )
  }

  # --- 4. Plot Aesthetics and Limits ---
  if (is.na(x_label)) {
    if (use_qc_metrics) {
      if (stringr::str_detect(variable_clean, "cv")) {
        x_label <- "Coefficient of Variation (%)"
      } else if (stringr::str_detect(variable_clean, "rt")) {
        x_label <- "Retention Time"
      } else if (stringr::str_detect(variable_clean, "dratio")) {
        x_label <- "D-ratio"
      } else if (stringr::str_detect(variable_clean, "prop")) {
        x_label <- "Missingness (Proportion of Samples)"
      } else if (stringr::str_detect(variable_clean, "sb")) {
        x_label <- "Signal/Blank Ratio"
      } else if (stringr::str_detect(variable_clean, "r2")) {
        x_label <- "R2 of Response Curve"
      } else if (stringr::str_detect(variable_clean, "conc")) {
        if (
          data@is_quantitated &&
            data@status_processing == "Calibration-quantitated data"
        ) {
          conc_unit_origin <- unique(
            data@annot_qcconcentrations$concentration_unit
          )
        } else {
          conc_unit_origin <- "pmol"
        }
        d_analyses <- data@annot_analyses |>
          dplyr::filter(.data$qc_type %in% qc_types)
        unit <- get_conc_unit(d_analyses$sample_amount_unit, conc_unit_origin)
        x_label <- paste0("Concentration (", unit, ")")
      } else {
        x_label <- variable_clean
      }
    } else {
      if (
        data@is_quantitated &&
          data@status_processing == "Calibration-quantitated data"
      ) {
        conc_unit_origin <- unique(
          data@annot_qcconcentrations$concentration_unit
        )
      } else {
        conc_unit_origin <- "pmol"
      }
      d_analyses <- data@annot_analyses |>
        dplyr::filter(.data$qc_type %in% qc_types)
      if (stringr::str_detect(variable_clean, "conc")) {
        unit <- get_conc_unit(d_analyses$sample_amount_unit, conc_unit_origin)
        x_label <- paste0("Concentration (", unit, ")")
      } else {
        x_label <- stringr::str_replace_all(variable_clean, "_", " ") |>
          stringr::str_to_title()
      }
    }
  }

  auto_calculate_limits <- length(x_lim) == 1 && is.na(x_lim)
  if (auto_calculate_limits) {
    min_vals <- c(
      d_summary$abundance_min_padded,
      if (show_sum) d_summary$abundance_sum_mean else NA_real_
    )
    max_vals <- c(
      d_summary$abundance_max_padded,
      if (show_sum) d_summary$abundance_sum_mean else NA_real_
    )

    if (log_scale) {
      min_val <- min(min_vals[min_vals > 0], na.rm = TRUE)
      max_val <- max(max_vals, na.rm = TRUE)
      x_lim <- c(floor(log10(min_val)), ceiling(log10(max_val)))
    } else {
      x_lim <- c(min(min_vals, na.rm = TRUE), max(max_vals, na.rm = TRUE))
    }
  }

  if (log_scale) {
    breaks <- 10^(x_lim[1]:x_lim[2])
    minor_breaks <- rep(1:9, times = diff(x_lim) + 1) *
      10^rep(x_lim[1]:x_lim[2], each = 9)
    plot_limits <- c(10^x_lim[1], 10^x_lim[2])
  } else {
    plot_limits <- x_lim
  }

  # Proactively check for data outside user-specified limits and issue a clear warning.
  if (!auto_calculate_limits) {
    removed_features <- sum(
      d_features$abundance_mean < plot_limits[1] |
        d_features$abundance_mean > plot_limits[2],
      na.rm = TRUE
    )
    removed_sums <- 0
    if (show_sum) {
      removed_sums <- sum(
        d_summary$abundance_sum_mean < plot_limits[1] |
          d_summary$abundance_sum_mean > plot_limits[2],
        na.rm = TRUE
      )
    }

    if (removed_features > 0 || removed_sums > 0) {
      # Dynamically build the message parts
      msg_parts <- c()
      if (removed_features > 0) {
        msg_parts <- c(
          msg_parts,
          cli::pluralize("{removed_features} feature{?s}")
        )
      }
      if (removed_sums > 0) {
        msg_parts <- c(
          msg_parts,
          cli::pluralize("{removed_sums} class sum{?s}")
        )
      }

      # Combine parts into a single string and issue the warning
      msg_body <- paste(msg_parts, collapse = " and ")
      cli::cli_alert_warning(
        col_yellow(
          "Some data points fall outside the `x_lim` range and were removed: {msg_body}."
        ),
        .envir = environment()
      )
    }

    # Pre-filter the data to prevent ggplot2 warnings.
    d_features <- d_features |>
      dplyr::filter(
        .data$abundance_mean >= plot_limits[1] &
          .data$abundance_mean <= plot_limits[2]
      )

    if (show_sum) {
      d_summary <- d_summary |>
        dplyr::filter(
          .data$abundance_sum_mean >= plot_limits[1] &
            .data$abundance_sum_mean <= plot_limits[2]
        )
    }
  }

  # --- 5. Main Plot Construction ---
  plt <- ggplot2::ggplot(
    d_features,
    ggplot2::aes(x = .data$abundance_mean, y = as.numeric(.data$feature_class))
  ) +
    ggplot2::geom_rect(
      data = d_summary,
      ggplot2::aes(
        xmin = .data$abundance_min_padded,
        xmax = .data$abundance_max_padded,
        ymin = as.numeric(.data$feature_class) - 0.4,
        ymax = as.numeric(.data$feature_class) + 0.4,
        fill = .data$feature_class
      ),
      inherit.aes = FALSE,
      alpha = 0.5
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        xend = .data$abundance_mean,
        y = as.numeric(.data$feature_class) - 0.3,
        yend = as.numeric(.data$feature_class) + 0.3
      ),
      color = "black",
      linewidth = segment_width
    )

  if (show_sum) {
    plt <- plt +
      ggplot2::geom_point(
        data = d_summary,
        ggplot2::aes(
          x = .data$abundance_sum_mean,
          y = as.numeric(.data$feature_class),
          color = .data$feature_class
        ),
        size = 1.3,
        shape = 23, # Diamond
        stroke = 0.8
      )
  }

  plt <- plt +
    ggplot2::scale_y_continuous(
      position = y_axis_position,
      limits = c(0.5, length(levels(d_features$feature_class)) + 0.5),
      breaks = seq_along(levels(d_features$feature_class)),
      labels = levels(d_features$feature_class),
      expand = c(0.015, 0.01)
    ) +
    ggplot2::labs(x = x_label, y = NULL) +
    ggplot2::scale_fill_manual(
      values = feature_map_plot,
      drop = FALSE,
      guide = "none"
    ) +
    ggplot2::scale_color_manual(
      values = feature_map_plot,
      drop = FALSE,
      guide = "none"
    )

  if (log_scale) {
    plt <- plt +
      ggplot2::scale_x_log10(
        limits = plot_limits,
        breaks = breaks,
        minor_breaks = minor_breaks,
        labels = scales::trans_format("log10", scales::math_format()),
        expand = c(0.05, 0.002)
      ) +
      ggplot2::annotation_logticks(
        base = 10,
        sides = "b",
        linewidth = 0.3,
        colour = "grey80",
        long = ggplot2::unit(1, "mm"),
        mid = ggplot2::unit(0.5, "mm"),
        short = ggplot2::unit(0, "mm")
      )
  } else {
    plt <- plt +
      ggplot2::scale_x_continuous(
        limits = plot_limits,
        breaks = scales::pretty_breaks(n = 10),
        expand = c(0.05, 0.002)
      )
  }

  plt <- plt +
    ggplot2::theme_bw(base_size = font_base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.x = ggplot2::element_line(
        color = grid_minor_color,
        linetype = "dotted",
        linewidth = grid_minor_linewidth
      ),
      panel.grid.major.x = ggplot2::element_line(
        color = grid_major_color,
        linetype = "solid",
        linewidth = grid_major_linewidth
      ),
      axis.text.x = ggplot2::element_text(size = font_base_size - 1),
      axis.title.x = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = font_base_size - 1),
      axis.ticks = ggplot2::element_line(colour = "grey40", linewidth = 0.5),
      axis.ticks.length = ggplot2::unit(3, "pt"),
      panel.border = ggplot2::element_rect(
        color = "grey40",
        fill = NA,
        linewidth = 0.5
      ),
      plot.margin = ggplot2::unit(c(0, 0.0, 0.0, 0.0), "cm")
    )

  # --- 6. Optional Density Strip ---
  if (density_strip) {
    check_installed("patchwork")

    x_for_density <- if (log_scale) {
      log10(d_features$abundance_mean)
    } else {
      d_features$abundance_mean
    }

    x_range <- range(x_for_density)

    dens <- stats::density(x_for_density, n = 1000, na.rm = TRUE)
    dens_df <- data.frame(x = dens$x, density = dens$y) |>
      dplyr::filter(.data$x >= x_range[1] & .data$x <= x_range[2])
    dens_df$density <- dens_df$density / max(dens_df$density)
    tile_width <- diff(dens_df$x)[1]

    p_density <- ggplot2::ggplot(dens_df, ggplot2::aes(x = .data$x, y = 1)) +
      ggplot2::geom_tile(
        ggplot2::aes(fill = .data$density),
        height = 1,
        width = tile_width, na.rm = FALSE
      ) +
      ggplot2::scale_fill_gradientn(
        colors = grDevices::colorRampPalette(c(
          "white",
          "#edf3f4ff",
          "#b8e4f7",
          "#72b1e4"
        ))(200)
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0.5, 1.5),
        breaks = 1,
        labels = "All",
        position = y_axis_position,
        expand = c(0.002, 0.002)
      ) +
      geom_segment(
        data = d_features,
        aes(
          x = if (log_scale) {
            log10(.data$abundance_mean)
          } else {
            .data$abundance_mean
          },
          y = 1 - 0.3,
          xend = if (log_scale) {
            log10(.data$abundance_mean)
          } else {
            .data$abundance_mean
          },
          yend = 1 + 0.3,
        ),
        #fill = "black",
        color = "#226ca1",
        alpha = 0.90,
        linewidth = 0.25,
        show.legend = FALSE,
        inherit.aes = TRUE
      ) +
      ggplot2::theme_bw(base_size = font_base_size) +
      ggplot2::theme(
        legend.position = "none",
        panel.border = ggplot2::element_rect(
          color = "grey40",
          fill = NA,
          linewidth = 0.5
        ),
        axis.ticks.length = ggplot2::unit(3, "pt"),
        axis.text.y = ggplot2::element_text(size = font_base_size - 1),
        axis.ticks.y = ggplot2::element_line(
          colour = "grey40",
          linewidth = 0.5
        ),
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title = ggplot2::element_blank(),
        plot.margin = ggplot2::unit(c(0.0, 0.0, 0.0, 0.0), "cm"),
        panel.grid = ggplot2::element_blank()
      )

    if (log_scale) {
      p_density <- p_density +
        ggplot2::scale_x_continuous(
          limits = log10(plot_limits),
          expand = c(0.05, 0.002)
        )
    } else {
      p_density <- p_density +
        ggplot2::scale_x_continuous(
          limits = plot_limits,
          breaks = scales::pretty_breaks(n = 10),
          expand = c(0.05, 0.002)
        )
    }

    # Combine plots using patchwork
    plt <- p_density / plt + patchwork::plot_layout(heights = c(0.4, 9.6))
  }
  plt + theme(plot.margin = margin(0, 0, 0, 0))
}
