#' Generate Correlation Matrix in Long Format
#'
#' @description
#' Creates a correlation matrix and transforms it to long format, filtering by correlation thresholds.
#'
#' @param tbl A data frame containing numeric columns for correlation analysis
#' @param cor_min_neg Numeric. Lower Pearson's correlation threshold
#' @param cor_min Numeric. Upper Pearson's correlation threshold
#'
#' @return A data frame in long format containing filtered correlations
#' @keywords internal
#'
#' @importFrom stats cor
get_feature_correlations <- function(tbl, cor_min_neg, cor_min) {

  mat <- tbl |>
    tibble::column_to_rownames("analysis_id") |>
    dplyr::select(-"qc_type") |>
    dplyr::select(where(is.numeric)) |>
    as.matrix()

  stats::cor(mat, method = "pearson") |>
    as.data.frame() |>
    tibble::rownames_to_column("var1") |>
    tidyr::pivot_longer(
      cols = -"var1",
      names_to = "var2",
      values_to = "value"
    ) |>
    dplyr::filter(.data$var1 < .data$var2) |>  # Keep only upper triangle
    dplyr::filter(.data$value <= cor_min_neg | .data$value >= cor_min)
}

#' Plot Highly Correlated Feature Pairs
#'
#' @description
#' Creates scatter plots for pairs of features that have correlations outside specified thresholds.
#' Each pair is displayed in a separate facet with its correlation coefficient.
#'
#' @param data A data frame containing numeric columns for correlation analysis
#'
#' @param variable A character string indicating the variable to use for PCA
#' analysis. Must be one of: "area", "height", "intensity", "norm_intensity", "response",
#' "conc", "conc_raw", "rt", "fwhm".
#'
#' @param qc_types A character vector specifying the QC types to plot. It
#' must contain at least one element. The default is `NA`, which means any
#' of the non-blank QC types ("SPL", "TQC", "BQC", "HQC", "MQC", "LQC",
#' "QC", "NIST", "LTR") will be plotted if present in the dataset.
#'
#' @param cor_min Numeric. Minimum correlation threshold. Only feature pairs with
#' positive correlations above this value will be shown. Set to Inf to exclude positive
#' corrections.
#'
#' @param cor_min_neg Numeric. Minimum nagative correlation threshold. Only feature pairs with
#' negative correlations above this value will be shown. Set to -Inf to exclude nagative
#' corrections.
#' @param log_scale A logical value indicating whether to use a log10 scale for
#' both axes. Default is `FALSE`.

#' @param filter_data A logical value indicating whether to use all data
#' (default) or only QC-filtered data (filtered via [filter_features_qc()]).
#' @param include_qualifier A logical value indicating whether to include
#' qualifier features. Default is `TRUE`.
#' @param include_istd A logical value indicating whether to include internal
#' standard (ISTD) features. Default is `TRUE`.
#' @param include_feature_filter A character or regex pattern used to filter
#' features by `feature_id`. If `NA` or an empty string (`""`) is provided,
#' the filter is ignored. When a vector of length > 1 is supplied, only
#' features with exactly these names are selected (applied individually as
#' OR conditions).
#' @param exclude_feature_filter A character or regex pattern used to exclude
#' features by `feature_id`. If `NA` or an empty string (`""`) is provided,
#' the filter is ignored. When a vector of length > 1 is supplied, only
#' features with exactly these names are excluded (applied individually as
#' OR conditions).
#' @param min_median_value Minimum median
#' feature value (as determined by the `variable`) across all samples from
#' selected QC types that must be met for a feature to be included in the
#' PCA analysis. `NA` (default) means no filtering will be applied. This
#' parameter provides an fast way to exclude noisy features from the
#' analysis. However, it is recommended to use `filter_data` with
#' [filter_features_qc()].
#'
#' @param point_size A numeric value indicating the size of points in
#' millimeters. Default is 1.
#' @param point_alpha A numeric value indicating the transparency of
#' points (0-1). Default is 0.8.
#' @param point_stroke A numeric value indicating the stroke width of the points. Default is 0.3.
#' @param line_size A numeric value indicating the size of the correlation line. Default is 0.5.
#' @param line_color A character string indicating the color of the correlation line. Default is orange.
#' @param line_alpha A numeric value indicating the transparency of the correlation line (0-1). Default is 0.5.
#'
#' @param font_base_size A numeric value indicating the base font size for
#' plot text elements. Default is 8.

#'
#' @return A ggplot object showing scatter plots of highly correlated feature pairs.
#' Returns NULL if no correlations meet the threshold criteria.
#'

plot_feature_correlations <- function(data,
                              variable,
                              qc_types = NA,
                              cor_min,
                              cor_min_neg = -0.99,
                              log_scale = FALSE,
                              filter_data = FALSE,
                              include_qualifier = FALSE,
                              include_istd = FALSE,
                              include_feature_filter = NA,
                              exclude_feature_filter = NA,
                              min_median_value = NA,
                              point_size = 1,
                              point_alpha = 0.8,
                              point_stroke = 0.3,
                              line_size = 0.5,
                              line_color = "orange",
                              line_alpha  = 0.5,
                              font_base_size = 8) {

  check_data(data)

  variable <- str_remove(variable, "feature_")
  rlang::arg_match(variable, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm", "width", "symmetry"))
  variable <- stringr::str_c("feature_", variable)
  check_var_in_dataset(data@dataset, variable)
  variable_sym = rlang::sym(variable)



  if(all(is.na(qc_types))){
    qc_types <- intersect(data$dataset$qc_type, c("SPL", "TQC", "BQC", "TQC", "HQC", "MQC", "LQC", "QC", "NIST", "LTR"))
  }

  # Subset dataset according to filter arguments
  # -------------------------------------
  d_filt <- get_dataset_subset(
    data,
    filter_data = filter_data,
    qc_types = qc_types,
    include_qualifier = include_qualifier,
    include_istd = include_istd,
    include_feature_filter = include_feature_filter,
    exclude_feature_filter = exclude_feature_filter
  )

  d_filt <- d_filt |>
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_id", {{ variable }})

  if(!is.na(min_median_value)){
    d_minsignal <- d_filt |>
      summarise(median_signal = median(!!variable_sym, na.rm = TRUE), .by = "feature_id") |>
      filter(.data$median_signal >= min_median_value)
    if(nrow(d_minsignal) == 0)
      cli_abort(col_red("No features passed the `min_median_value` filter. Please review the filter value, `variable` and data."))
    else if(nrow(d_minsignal) == 1)
      cli_abort(col_red("Only 1 feature passed the `min_median_value` filter. Please review the filter value, `variable`, and data."))

    d_filt <- d_filt |> semi_join(d_minsignal, by = "feature_id")
  }

  d_wide <- d_filt |>
    tidyr::pivot_wider(id_cols = c("analysis_id", "qc_type"), names_from = "feature_id", values_from = all_of(variable))


  if (cor_min_neg >= cor_min) {
    cli_abort(col_red("Lower correlation threshold must be less than upper threshold"))
  }

  # Generate correlation matrix
  cor_matrix <- get_feature_correlations(d_wide, cor_min_neg, cor_min)

  if (nrow(cor_matrix) == 0) {
    cli_alert_info("No correlations found exceeding thresholds (r < {.val {cor_min_neg}} or r > {.val {cor_min}})")
    return(NULL)
  }

  # Create plotting data
  plot_data <- purrr::map_df(1:nrow(cor_matrix), function(i) {
    var1 <- cor_matrix$var1[i]
    var2 <- cor_matrix$var2[i]
    cor_val <- round(cor_matrix$value[i], 3)

    data.frame(
      analysis_id = d_wide$analysis_id,
      qc_type = d_wide$qc_type,
      x = d_wide[[var1]],
      y = d_wide[[var2]],
      pair = sprintf("%s\n%s", var1, var2),
      r = sprintf("r = %.3f", cor_val),
      abs_cor = abs(cor_val),
      stringsAsFactors = FALSE
    )
  })  |>
    arrange(desc(.data$abs_cor)) |>
    dplyr::mutate(pair = factor(.data$pair, levels = unique(.data$pair)))



  plot_data$qc_type <- droplevels(factor(plot_data$qc_type, levels = c("SPL", "UBLK", "SBLK", "TQC", "BQC", "RQC", "LTR", "NIST", "PBLK")))
  plot_data <- plot_data |>
    dplyr::arrange(.data$qc_type)

  # Create plot
  p <- plot_data |>
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(size = point_size,
                        aes(color = .data$qc_type, shape = .data$qc_type, fill = .data$qc_type),
                        alpha = point_alpha,
                        stroke = point_stroke) +
    ggplot2::scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
    ggplot2::scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
    ggplot2::scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)
    if (log_scale) {
      p <- p +
        ggplot2::scale_x_log10(labels = scientific_format_end,
                               expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::scale_y_log10(labels = scientific_format_end,
                               expand = ggplot2::expansion(mult = c(0, 0.05)))
    } else {
      p <- p +
        ggplot2::scale_x_continuous(labels = scientific_format_end,
                                    limits = function(x) c(0, max(x)),
                                    expand = ggplot2::expansion(mult = c(0, 0.05))) +
        ggplot2::scale_y_continuous(labels = scientific_format_end,
                                    limits = function(x) c(0, max(x)),
                                    expand = ggplot2::expansion(mult = c(0, 0.05)))
    }

  p <- p +
    suppressWarnings({ggplot2::geom_smooth(method = "lm", formula = y ~ x, linewidth = line_size, color = line_color, alpha = line_alpha, se = FALSE, , na.rm = TRUE)}) +
    ggplot2::geom_text(
      data = plot_data |>
        dplyr::group_by(.data$pair, .data$r) |>
        dplyr::slice(1),
      ggplot2::aes(label = .data$r),
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      size = 2.5
    ) +
    ggplot2::facet_wrap(~pair, scales = "free") +
    theme_bw(base_size = font_base_size) +
    theme(
      plot.title = element_text(size = font_base_size, face = "bold"),
      strip.text = ggplot2::element_text(size = font_base_size, face = "bold"),
      axis.text = element_text(size = font_base_size),
      axis.title = element_text(size = font_base_size, face = "bold"),
      panel.grid = element_line(linewidth = 0.001),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2),
      strip.background = ggplot2::element_rect(linewidth = 0.0001, fill = "#00283d"),
      strip.text.x = ggplot2::element_text(color = "white"),
      #strip.switch.pad.wrap = ggplot2::unit(1, "mm"),
      panel.border = element_rect(linewidth = 0.5, color = "grey40"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      x = "Feature 1",
      y = "Feature 2"
    )

  p
}
