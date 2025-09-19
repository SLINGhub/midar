#' RunSequence Plot
#'
#' @description
#' The RunSequence plot provides an overview of the analysis design and
#' timelines, which can be useful for subsequent processing steps. The plot
#' illustrates the batch structure, the quality control (QC) samples included
#' with their respective positions, and additional information regarding the
#' date, duration, and run time of the analysis.
#'
#' Setting `show_timestamp = TRUE` allows you to check for any interruptions in
#' the analysis timeline.
#'
#' @param data MidarExperiment object
#' @param qc_types QC types to be plotted. Can be a vector of QC types or a
#' regular expression pattern. `NA` (default) displays all available QC/Sample
#' types.
#' @param show_batches Logical, whether to show batch separators in the plot.
#' @param show_timestamp Logical, whether to use the acquisition timestamp as
#' the x-axis instead of the run sequence number.
#' @param add_info_title Logical, whether to add a title with the experiment
#' title, analysis date, and analysis times.
#' @param single_row Logical, whether to show all QC types in a single row.
#' @param segment_linewidth Width of the segment lines, default is 0.5.
#' @param batch_zebra_stripe Logical, whether to show batches as shaded areas
#' instead of line separators.
#' @param batch_line_color Color of the batch separator lines.
#' @param batch_fill_color Color of the batch shaded areas.
#' @param base_font_size Numeric, base font size for the plot.
#' @return A ggplot object representing the run sequence plot.
#' @export
plot_runsequence <- function(
  data = NULL,
  qc_types = NA,
  show_batches = TRUE,
  show_timestamp = FALSE,
  add_info_title = TRUE,
  single_row = FALSE,
  segment_linewidth = 0.5,
  batch_zebra_stripe = FALSE,
  batch_line_color = "#b6f0c5",
  batch_fill_color = "grey90",
  base_font_size = 8
) {
  # Check if data is valid
  check_data(data)

  # Extract the required columns from the dataset
  d_filt <- data$dataset |>
    select(
      "analysis_order",
      "acquisition_time_stamp",
      "batch_id",
      "analysis_id",
      "qc_type"
    ) |>
    distinct()

  # Filter QC types if provided
  if (!all(is.na(qc_types)) && length(qc_types) > 0) {
    d_filt <- d_filt |>
      filter(
        if (is.vector(qc_types) && length(qc_types) > 1) {
          .data$qc_type %in% qc_types
        } else {
          str_detect(.data$qc_type, qc_types)
        }
      )
  }

  # Convert acquisition_time_stamp to POSIXct if using datetime
  if (show_timestamp) {
    d_filt$acquisition_time_stamp <- as.POSIXct(d_filt$acquisition_time_stamp)
  }

  # Convert qc_type to factor and create sample_category
  # d_filt$qc_type <- factor(d_filt$qc_type,
  #                          levels = pkg.env$qc_type_annotation$qc_type_levels) |>
  d_filt$qc_type <- factor(d_filt$qc_type) |>
    forcats::fct_drop()
  d_filt$sample_category <- as.character(d_filt$qc_type)

  # Count samples per analysis_order
  sample_counts <- d_filt |>
    group_by(.data$qc_type) |>
    summarise(sample_count = stringr::str_c("", dplyr::n()), .groups = 'drop')

  # Define QC colors
  qc_colors <- replace(pkg.env$qc_type_annotation$qc_type_col, "SPL", "grey35")

  # Initialize the ggplot object
  p <- ggplot(
    d_filt,
    aes(
      x = if (show_timestamp) {
        .data$acquisition_time_stamp
      } else {
        .data$analysis_order
      },
      y = rev(.data$qc_type),
      color = .data$sample_category
    )
  ) +
    labs(
      x = if (show_timestamp) "Acquisition Time" else "Analysis Order",
      y = "Sample Type"
    ) +
    theme_bw(base_size = base_font_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(
        colour = "grey80",
        linetype = "dotted",
        linewidth = 0.25
      ),
      panel.grid.minor.x = element_line(
        colour = "grey90",
        linetype = "dotted",
        linewidth = 0.25
      ),
      panel.border = element_rect(linewidth = 1),
      axis.title = element_text(face = "bold", size = base_font_size),
      axis.text.x = element_text(face = "plain", size = base_font_size),
      axis.text.y = element_text(face = "plain", size = base_font_size),
      axis.text.y.right = element_text(face = "plain", size = base_font_size),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_text(
        angle = 0,
        vjust = 0.5,
        size = base_font_size
      ),
      axis.ticks.x = element_line(
        colour = "grey80",
        linetype = "dotted",
        linewidth = 0.5
      ),
      plot.margin = unit(c(1, 1, 1, 1), "mm"),
      legend.position = if (single_row) "right" else "none" # Show legend if single_row
    )

  # Add batch shading if defined
  if (show_batches) {
    # Retrieve batch info
    d_batch_info <- data@annot_batches |>
      left_join(
        data$dataset |>
          select("analysis_order", "acquisition_time_stamp", "batch_id") |>
          distinct(),
        by = c("id_batch_start" = "analysis_order"),
        suffix = c("", "_start"),
        keep = FALSE
      ) |>
      left_join(
        data$dataset |>
          select("analysis_order", "acquisition_time_stamp", "batch_id") |>
          distinct(),
        by = c("id_batch_end" = "analysis_order"),
        suffix = c("", "_end"),
        keep = FALSE
      )

    if (batch_zebra_stripe) {
      d_batch_shading <- d_batch_info |>
        slice(-1) |>
        filter(.data$batch_no %% 2 != 1)
      p <- p +
        geom_rect(
          data = d_batch_shading,
          inherit.aes = FALSE,
          aes(
            xmin = if (show_timestamp) {
              .data$acquisition_time_stamp
            } else {
              .data$id_batch_start + 0.5
            },
            xmax = if (show_timestamp) {
              .data$acquisition_time_stamp_end
            } else {
              .data$id_batch_end - 0.5
            },
            ymin = -Inf,
            ymax = Inf
          ),
          fill = batch_fill_color,
          alpha = 1,
          color = NA
        )
    } else {
      p <- p +
        geom_vline(
          data = d_batch_info,
          aes(
            xintercept = if (show_timestamp) {
              .data$acquisition_time_stamp
            } else {
              (.data$id_batch_end + 0.5)
            }
          ),
          colour = batch_line_color,
          linewidth = segment_linewidth * 2
        )
    }
  }

  # Add segments for qc_type
  if (single_row) {
    p <- p +
      geom_segment(
        aes(
          x = if (show_timestamp) {
            .data$acquisition_time_stamp
          } else {
            .data$analysis_order
          },
          xend = if (show_timestamp) {
            .data$acquisition_time_stamp
          } else {
            .data$analysis_order
          },
          y = -1,
          yend = 1
        ),
        linewidth = segment_linewidth
      ) +
      scale_y_continuous(breaks = NULL) # Hide y-axis breaks
  } else {
    p <- p +
      geom_segment(
        aes(
          x = if (show_timestamp) {
            .data$acquisition_time_stamp
          } else {
            .data$analysis_order
          },
          xend = if (show_timestamp) {
            .data$acquisition_time_stamp
          } else {
            .data$analysis_order
          },
          y = as.integer(.data$qc_type) - 0.4,
          yend = as.integer(.data$qc_type) + 0.4
        ),
        linewidth = segment_linewidth
      )

    # Position sample counts outside the y-axis
    p <- p +
      scale_y_continuous(
        breaks = seq(1, nlevels(d_filt$qc_type), by = 1),
        labels = sample_counts$qc_type,
        expand = expansion(0.02, 0.02),
        sec.axis = ggplot2::sec_axis(
          ~.,
          name = "n",
          breaks = seq(1, nlevels(d_filt$qc_type), by = 1),
          labels = sample_counts$sample_count
        )
      )
  }

  # Format x-axis as date-time if using acquisition_time_stamp
  if (show_timestamp) {
    p <- p +
      ggplot2::scale_x_datetime(
        date_labels = "%Y-%m-%d",
        expand = expansion(0.02, 0.02),
        date_breaks = "day",
        date_minor_breaks = "hour"
      )
  } else {
    p <- p +
      scale_x_continuous(
        expand = expansion(0.02, 0.02),
        breaks = seq(
          0,
          max(d_filt$analysis_order),
          10^ceiling(log10(max(d_filt$analysis_order))) / 10
        )
      )
  }

  # Add additional information in the title
  if (add_info_title) {
    title_text <- if (data@title == "") "A" else glue::glue("{data@title} - A")

    p <- p +
      labs(
        title = glue::glue(
          "{title_text}nalysis time: {get_analysis_duration(data, estimate_sequence_end = TRUE) |> stringr::str_sub(end = -5)}  ({get_analyis_start(data) |> stringr::str_sub(end = -4)} - {get_analyis_end(data, estimate_sequence_end = TRUE) |> stringr::str_sub(end = -4)})
                                     Median run time: {get_runtime_median(data)@minute}: {get_runtime_median(data)@.Data} min - interruptions > 1 hour: {get_analysis_breaks(data, 60)}"
        )
      )
  }

  # Color mapping
  p <- p + scale_color_manual(values = qc_colors, name = "QC type")

  p
}


#' Relative Log Abundance (RLA) Plot
#'
#' @description
#' Relative log abundance (RLA) plots show standardized feature abundances across samples. Standardization is done by removing either the within-batch or across-batch median from each feature
#'
#' RLA plots are useful for visualizing technical effects that impact all features in a similar manner, such as batch effects due to changes in instrument response, pipetting errors, or sample spillage. Unlike plots of raw or normalized abundances, RLA plots are more robust to these types of effects.
#'
#' @param data MidarExperiment
#' @param rla_type_batch Character, must be either "within" or "across", defining whether to use within-batch or across-batch RLA
#' @param variable Variable to plot, must be one of "intensity", "norm_intensity", "conc", "area", "height", "fwhm", or one of
#' "intensity_raw", "intensity_before", "norm_intensity_raw", "norm_intensity_before", "conc_raw", "conc_before"
#' @param qc_types QC types to be plotted. Can be a vector of QC types or a regular expression pattern. `NA` (default) displays all available QC/Sample types.
#' @param plot_range Numeric vector of length 2, specifying the start and end indices of the analysis order to be plotted. `NA` plots all samples.
#' @param rla_limit_to_range Logical, whether to limit the RLA values to the specified `plot_range`. Default is `FALSE`, which means RLA values are calculated for all samples.
#' @param remove_gaps Logical, whether to remove gaps in the x-axis, occuring from QC types that were not selected. Default is `TRUE`.
#' is supplied, only features with exactly these names are excluded (applied individually as OR conditions).
#'
#' @param filter_data Logical, whether to use QC-filtered data based on criteria set via `filter_features_qc()`.

#' @param include_qualifier Logical, whether to include qualifier features. Default is `TRUE`.
#' @param include_istd Logical, whether to include internal standard (ISTD) features. Default is `TRUE`.
#' @param include_feature_filter A regex pattern or a vector of feature names used to filter features by `feature_id`.
#' If `NA` or an empty string (`""`) is provided, the filter is ignored. When a vector of length > 1 is supplied,
#' is supplied, only features with exactly these names are selected (applied individually as OR conditions).
#' @param exclude_feature_filter A regex pattern or a vector of feature names to exclude features by feature_id.
#' If `NA` or an empty string (`""`) is provided, the filter is ignored. When a vector of length > 1 is supplied.

#' @param show_timestamp Logical, whether to use the acquisition timestamp as
#' the x-axis instead of the run sequence number
#'
#'
#' @param outlier_detection Logical, whether to show outlier fences on the plot and return a table with detect outliers based on the method defined by `outlier_method`.
#' @param outlier_hide Logical, whether to exclude outlier values from the plot. Default is `FALSE`, which means outliers are shown.
#' @param outlier_method Character, method used for outlier detection. Default is "mad" (median absolute deviation).
#' Other possible values are "iqr", "sd", "z_normal", "z_robust", "quantile", and "fold". See get_outlier_bounds() for details.
#' @param outlier_qctypes Character vector, QC types to use for outlier detection. Default is `c("SPL", "TQC", "BQC")`.
#' @param outlier_k Numeric, multiplier for the outlier detection method. Default is `NULL`, which uses the default value for the selected method.
#' See get_outlier_bounds() for details. When using the "fold" method, either single numeric value or a vector with two values (lower and upper fences) can be supplied.
#'
#' @param min_feature_intensity Numeric, exclude features with overall median signal below this value
#' @param y_lim Numeric vector of length 2, specifying the lower and upper y-axis limits. Default is `NA`, which uses limits calculated based on `outlier_hide`.

#' @param show_batches Logical, whether to show batch separators in the plot
#' @param batch_zebra_stripe Logical, whether to show batches as shaded areas instead of line separators
#' @param batch_line_color Character, color of the batch separator lines
#' @param batch_fill_color Character, color of the batch shaded areas

#' @param x_gridlines Logical, whether to show major x-axis gridlines
#' @param linewidth Numeric, line width used for whiskers of the boxplot
#' @param base_font_size Numeric, base font size for the plot
#' @param relative_log_abundances Logical, whether to use relative log abundances (RLA) or just log-transformed values
#' @param show_plot Logical, whether to display the plot. Default is `TRUE`.
#' @return A list with the ggplot object representing the RLA plot and a table with detected outliers if `outlier_detection = TRUE`.
#' @references
#' De Livera et al. (2012) Normalizing and integrating metabolomics data. Analytical Chemistry 10768-10776
#' [DOI: 10.1021/ac302748b](https://doi.org/10.1021/ac302748b)
#' De Livera et al. (2015) Statistical Methods for Handling Unwanted Variation in Metabolomics Data. Analytical Chemistry 87(7):3606-3615
#' [DOI: 10.1021/ac502439y](https://doi.org/10.1021/ac502439y)
#' @export

# TODO: Add minor ticks to x-axis
plot_rla_boxplot <- function(
  data = NULL,
  rla_type_batch = c("within", "across"),
  variable = c(
    "intensity",
    "norm_intensity",
    "conc",
    "conc_raw",
    "area",
    "height",
    "fwhm",
    "width"
  ),
  qc_types = NA,
  plot_range = NA,
  rla_limit_to_range = FALSE,
  remove_gaps = TRUE,

  filter_data = FALSE,
  include_qualifier = TRUE,
  include_istd = TRUE,
  include_feature_filter = NA,
  exclude_feature_filter = NA,

  show_timestamp = FALSE,

  min_feature_intensity = 0,
  y_lim = NA,

  outlier_detection = TRUE,
  outlier_hide = FALSE,
  outlier_method = "mad",
  outlier_qctypes = c("SPL", "TQC", "BQC", "LTR", "NIST"),
  outlier_k = NULL,
  show_batches = TRUE,
  batch_zebra_stripe = FALSE,
  batch_line_color = "#b6f0c5",
  batch_fill_color = "grey93",

  x_gridlines = FALSE,
  linewidth = 0.2,
  base_font_size = 8,
  relative_log_abundances = TRUE,
  show_plot = TRUE
) {
  check_data(data)
  if (nrow(data@dataset) < 1) {
    cli::cli_abort("No data available. Please import data and metadata first.")
  }

  rlang::arg_match(rla_type_batch, c("within", "across"))

  # Check if selected variable is valid
  rlang::arg_match(
    variable,
    c(
      "intensity",
      "norm_intensity",
      "conc",
      "intensity_raw",
      "norm_intensity_raw",
      "conc_raw",
      "area",
      "height",
      "rt",
      "fwhm",
      "response",
      "width",
      "symmetry"
    )
  )
  variable <- str_remove(variable, "feature_")
  variable <- stringr::str_c("feature_", variable)

  # Check arguments are valid
  check_var_in_dataset(data@dataset, variable)
  if (
    str_detect(variable, "\\_raw") &&
      !any(data@var_drift_corrected) &&
      !any(data@var_batch_corrected)
  ) {
    cli::cli_abort(cli::col_red(
      "`{variable} is only available after drift or/and batch correction. Please run drift and/or batch corrections, or choose another variable."
    ))
  }
  variable_sym = rlang::sym(variable)

  # Subset dataset according to arguments
  d_filt <- get_dataset_subset(
    data,
    filter_data = filter_data,
    qc_types = qc_types,
    include_qualifier = include_qualifier,
    include_istd = include_istd,
    include_feature_filter = include_feature_filter,
    exclude_feature_filter = exclude_feature_filter
  )
  if (show_timestamp) {
    d_filt$acquisition_time_stamp <- as.POSIXct(d_filt$acquisition_time_stamp)
    x_axis_variable <- "acquisition_time_stamp"
  } else {
    x_axis_variable <- "analysis_order"
  }

  x_axis_variable_sym <- rlang::sym(x_axis_variable)

  if (!all(is.na(plot_range)) && rla_limit_to_range) {
    d_filt <- d_filt |>
      dplyr::filter(
        .data$analysis_order >= plot_range[1] &
          .data$analysis_order <= plot_range[2]
      ) |>
      droplevels()
  }

  d_filt <- d_filt |>
    mutate(value = ifelse(is.infinite(!!variable_sym), NA, !!variable_sym))

  d_filt <- d_filt |>
    dplyr::select(any_of(c(
      "analysis_id",
      "analysis_order",
      "acquisition_time_stamp",
      "qc_type",
      "batch_id",
      "feature_id",
      "feature_intensity",
      "feature_norm_intensity",
      "feature_conc"
    ))) |>
    group_by(.data$feature_id) |>
    filter(median(.data$feature_intensity) >= min_feature_intensity) |>
    droplevels() |>
    dplyr::arrange(.data$analysis_order)

  if (relative_log_abundances) {
    d_filt = d_filt |>
      mutate(val = log2(!!variable_sym))

    if (rla_type_batch == "within") {
      grp = c("feature_id", "batch_id")
    } else {
      grp = c("feature_id")
    }

    d_filt_medians <- d_filt |>
      group_by(across(all_of(grp))) |>
      summarise(val_median = median(.data$val, na.rm = TRUE)) |>
      ungroup()

    d_filt <- d_filt |>
      left_join(d_filt_medians, by = grp) |>
      mutate(val_res = .data$val - .data$val_median)
  } else {
    d_filt <- d_filt |> mutate(val_res = log2(!!variable_sym))
  }

  unique_orders <- sort(unique(d_filt$analysis_order))

  order_map <- tibble(
    analysis_order = unique_orders,
    analysis_order_index = seq_along(unique_orders)
  )

  d_filt <- d_filt |>
    left_join(order_map, by = "analysis_order")

  if (outlier_detection) {
    # Print outliers
    # outlier_qc_types <- rlang::arg_match(
    #   outlier_qctypes,
    #   pkg.env$qc_type_annotation$qc_type_levels
    # )

    d_sum <- d_filt |>
      filter(.data$qc_type %in% outlier_qctypes) |>
      group_by(.data$analysis_id, .data$analysis_order, .data$batch_id) |>
      summarise(
        val_res_median = median(.data$val_res, na.rm = TRUE),
        qc_type = dplyr::first(.data$qc_type),
        .groups = "drop"
      ) |>
      arrange(.data$analysis_order)

    outlier_bounds <- get_outlier_bounds(
      d_sum$val_res_median,
      method = outlier_method,
      k = outlier_k,
      na.rm = TRUE
    )
    if (outlier_method != "fold") {
      outlier_bounds <- get_outlier_bounds(
        d_sum$val_res_median,
        k = outlier_k,
        method = outlier_method,
        na.rm = TRUE
      )
    } else {
      if (length(outlier_k) == 1) {
        outlier_bounds <- c(-outlier_k, outlier_k)
      } else if (length(outlier_k) == 2) {
        outlier_bounds <- outlier_k
      } else {
        cli::cli_abort(
          "When using the 'fold' method, `outlier_k` must be a single numeric value or a vector of two values (lower and upper fences)."
        )
      }
    }

    d_outliers <- d_sum |>
      filter(
        .data$val_res_median < outlier_bounds[1] |
          .data$val_res_median > outlier_bounds[2]
      ) |>
      ungroup()

    if (nrow(d_outliers) > 0) {
      cli::cli_alert_info(
        "Found {nrow(d_outliers)} outliers in the {length(unique(d_filt$analysis_id))} shown analyses"
      )
    } else {
      cli::cli_alert_info("No outliers found.")
    }
  }

  # Get labels corresponding to the breaks. TODO: write it more elegant and clear
  if (x_axis_variable != "analysis_order") {
    labels <- unique(d_filt[[x_axis_variable]])[seq(
      1,
      length(unique(d_filt[[x_axis_variable]])),
      length.out = 10
    )]
    breaks <- d_filt |>
      filter(!!x_axis_variable_sym %in% labels) |>
      pull(.data$analysis_order) |>
      unique()
    p <- ggplot(
      d_filt,
      aes(
        x = .data$analysis_order,
        y = .data$val_res,
        group = .data$analysis_order
      )
    )
  } else {
    n_breaks <- 10
    breaks <- scales::breaks_pretty(n = n_breaks)(seq_along(unique_orders))
    breaks <- breaks[breaks > 0 & breaks <= length(unique_orders)]
    if (remove_gaps) {
      # If remove_gaps is TRUE, we use the analysis_order_index
      labels <- unique_orders[breaks] |> round(-2) |> format(big.mark = ",")
      p <- ggplot(
        d_filt,
        aes(
          x = .data$analysis_order_index,
          y = .data$val_res,
          group = .data$analysis_order_index
        )
      )
    } else {
      # If remove_gaps is FALSE, we use the analysis_order
      labels <- breaks |> round(-2) |> format(big.mark = ",")
      p <- ggplot(
        d_filt,
        aes(
          x = .data$analysis_order,
          y = .data$val_res,
          group = .data$analysis_order
        )
      )
    }
  }
  #TODO cleanup
  p <- p +
    scale_x_continuous(
      breaks = breaks,
      labels = if(show_timestamp) labels else breaks, #ToDo
      # limits = if (remove_gaps) {
      #   range(d_filt$analysis_order_index)
      # } else {
      #   range(d_filt$analysis_order)
      # },
      expand = c(0.02, 0.02)
    )

  if (show_batches) {
    d_batches <- data@annot_batches |>
      filter(
        .data$id_batch_start <= max(order_map$analysis_order) &
          .data$id_batch_start >= min(order_map$analysis_order)
      ) |>
      mutate(
        mapped_start = purrr::map_dbl(
          .data$id_batch_start,
          ~ find_closest(.x, order_map$analysis_order, method = "higher")
        ),
        mapped_end = purrr::map_dbl(
          .data$id_batch_end,
          ~ find_closest(.x, order_map$analysis_order, method = "lower")
        )
      ) |>
      left_join(order_map, by = c("mapped_start" = "analysis_order")) |>
      rename(id_batch_start_index = "analysis_order_index") |>
      left_join(order_map, by = c("mapped_end" = "analysis_order")) |>
      rename(id_batch_end_index = "analysis_order_index")

    if (!batch_zebra_stripe) {
      if (remove_gaps) {
        p <- p +
          geom_vline(
            data = d_batches |> slice(-1),
            aes(xintercept = .data$id_batch_start_index - 0.5),
            colour = batch_line_color,
            linetype = "solid",
            linewidth = 1
          )
      } else {
        p <- p +
          geom_vline(
            data = d_batches |> slice(-1),
            aes(xintercept = .data$id_batch_start - 0.5),
            colour = batch_line_color,
            linetype = "solid",
            linewidth = 1
          )
      }
    } else {
      d_batch_2nd <- d_batches |>
        slice(-1) |>
        filter(.data$batch_no %% 2 != 1)
      if (remove_gaps) {
        p <- p +
          geom_rect(
            data = d_batch_2nd,
            aes(
              xmin = .data$id_batch_start_index - 0.5,
              xmax = .data$id_batch_end_index + 0.5,
              ymin = -Inf,
              ymax = Inf
            ),
            inherit.aes = FALSE,
            fill = batch_fill_color,
            color = NA,
            alpha = 1,
            linetype = "solid",
            linewidth = 0.5,
            na.rm = TRUE
          )
      } else {
        p <- p +
          geom_rect(
            data = d_batch_2nd,
            aes(
              xmin = .data$id_batch_start - 0.5,
              xmax = .data$id_batch_end + 0.5,
              ymin = -Inf,
              ymax = Inf
            ),
            inherit.aes = FALSE,
            fill = batch_fill_color,
            color = NA,
            alpha = 1,
            linetype = "solid",
            linewidth = 0.5,
            na.rm = TRUE
          )
      }
    }
  }

  x_text_angle <- ifelse(x_axis_variable != "analysis_order", 90, 0)
  x_text_just <- ifelse(x_axis_variable != "analysis_order", 1, 0.5)

  p <- p +
    geom_boxplot(
      aes(fill = .data$qc_type, color = .data$qc_type),
      notch = FALSE,
      outlier.colour = NA,
      linewidth = linewidth,
      na.rm = TRUE
    ) +
    #stat_summary(mapping = aes(color = .data$qc_type), fun = median, geom = "line", color = "red", linetype = "solid") +
    scale_fill_manual(
      name = NULL,
      values = desaturate_colors(pkg.env$qc_type_annotation$qc_type_col, 0.4)
    ) +
    scale_color_manual(
      name = NULL,
      values = pkg.env$qc_type_annotation$qc_type_col
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(title = NULL, override.aes = list(size = 3))
    ) +
    theme_bw(base_size = base_font_size) +
    ylab(bquote(bold(log[2] ~ .(variable)))) +
    xlab("Analysis order") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(), #element_line(colour = "#bdbdbd", linetype = "dotted", size = .5),
      panel.grid.minor.x = element_blank(), #element_line(colour = "grey88", linetype = "dotted", size = .5),
      axis.text.y = element_text(size = base_font_size),
      axis.text.x = element_text(
        size = base_font_size,
        angle = x_text_angle,
        vjust = 0.5,
        hjust = x_text_just
      ),
      axis.title = element_text(size = base_font_size * 1, face = "bold"),
      panel.border = element_rect(linewidth = 1, color = "grey20"),
      legend.position = "inside",
      legend.direction = "horizontal", # vertical layout
      legend.text = element_text(size = base_font_size * 0.8), # text size
      legend.title = element_text(size = base_font_size * 0.8),
      legend.key.size = unit(base_font_size * 0.8, "pt"), # box size
      legend.position.inside = c(0.6, 0.1)
    )

  if (outlier_detection) {
    p <- p +
      geom_hline(
        data = tibble(y = outlier_bounds),
        mapping = aes(yintercept = .data$y),
        color = "red",
        linewidth = 0.5,
        alpha = 0.5
      )
  }

  if (x_gridlines) {
    p <- p +
      theme(
        panel.grid.major.x = element_line(
          colour = "#bdbdbd",
          linetype = "dotted",
          linewidth = .3
        )
      )
  } else {
    p <- p + theme(panel.grid.major.x = element_blank())
  }

  if (relative_log_abundances) {
    p <- p +
      geom_hline(
        yintercept = 0,
        colour = "#5fe3f5",
        linetype = "longdash",
        linewidth = 0.5
      ) +
      ylab(bquote(bold(
        log[2] ~ "( relative" ~ .(stringr::str_remove(variable, "feature\\_")) ~
          ")"
      )))
  }
  ylim = c(NA, NA)
  xlim = c(NA, NA)
  # Set y-axis limits
  if (!all(is.na(y_lim))) {
    ylim = y_lim
  } else if (outlier_hide) {
    tails <- get_outlier_bounds(
      d_filt$val_res,
      method = outlier_method,
      k = outlier_k,
      na.rm = TRUE
    )
    ylim = tails
  }

  if (!all(is.na(plot_range))) {
    if (remove_gaps) {
      xlim <- c(
        order_map[
          order_map$analysis_order ==
            find_closest(
              plot_range[1],
              order_map$analysis_order,
              method = "lower"
            ),
        ]$analysis_order_index,
        order_map[
          order_map$analysis_order ==
            find_closest(
              plot_range[2],
              order_map$analysis_order,
              method = "higher"
            ),
        ]$analysis_order_index
      )
    } else {
      xlim = plot_range
    }
  }

  p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  res <-   list(
    outliers = if (is.null(d_outliers)) NULL else d_outliers,
    plot = p
  )

  if (show_plot) {
    print(p)
    invisible(res)
  } else {
    res
  }


}
