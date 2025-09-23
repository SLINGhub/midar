#' Plot standardized feature intensities grouped by QC type
#'
#' This function creates a grouped beeswarm plot of standardized feature intensities,
#' where the y-axis represents intensity standardized such that the mean across all
#' features is 100%. Points are grouped by `qc_type` and spread using quasirandom jitter.
#'
#' @param data A MidarExperiment object
#'
#' @param variable A character string indicating the variable to use for PCA
#' analysis. Must be one of: "area", "height", "intensity", "norm_intensity", "response",
#' "conc", "conc_raw", "rt", "fwhm".
#' @param qc_types A character vector specifying the QC types to plot. It
#' must contain at least one element. The default is `NA`, which means any
#' of the non-blank QC types ("SPL", "TQC", "BQC", "HQC", "MQC", "LQC",
#' "NIST", "LTR") will be plotted if present in the dataset.
#' @param batchwise_normalization A logical value indicating whether to normalize the signals by batch instead of globally.
#' @param include_qualifier A logical value indicating whether to include
#' qualifier features. Default is `TRUE`.
#' @param only_istd A logical value indicating whether to include features used as internal
#' standards (ISTD).  Default is `TRUE`. Set to `FALSE` in combination with feature_filter parameters to show other features.
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
#' @param y_lim A numeric vector of length 2 specifying the y-axis limits.
#' @param point_size A numeric value indicating the size of points in
#' millimeters. Default is 2.
#' @param dodge_width Numeric. Width used to dodge overlapping points by `qc_type`. Default is `0.6`.
#' @param point_alpha A numeric value indicating the transparency of
#' points (0-1). Default is 0.5.
#' @param font_base_size A numeric value indicating the base font size for
#' plot text elements. Default is 8.
#' @param point_alpha Numeric. Transparency of the plotted points. Default is `0.7`.
#' @param box_linewidth Numeric. Width of the boxplot lines. Default is `0.5`.
#' @param box_alpha Numeric. Transparency of the boxplot. Default is `0.3`.
#' @param angle_x Numeric. Angle of the x-axis text labels. Default is `45`.
#'
#' @return A `ggplot` object showing the grouped standardized beeswarm plot.
#' @export
#'

plot_qc_matrixeffects <- function(data,
                                  variable = "intensity",
                                  qc_types = c("SPL", "TQC", "PBLK", "BQC"),
                                  batchwise_normalization = TRUE,
                                  include_qualifier = FALSE,
                                  only_istd = TRUE,
                                  include_feature_filter = NA,
                                  exclude_feature_filter = NA,
                                  min_median_value = NA,
                                  y_lim = c(-NA, NA),
                                  point_size = 0.5,
                                  dodge_width = 0.6,
                                  point_alpha = 0.3,
                                  box_alpha = 0.3,
                                  box_linewidth = 0.5,
                                  font_base_size  = 8,
                                  angle_x = 45) {


  variable <- str_remove(variable, "feature_")
  rlang::arg_match(variable, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm", "width", "symmetry"))
  variable <- stringr::str_c("feature_", variable)
  check_var_in_dataset(data@dataset, variable)
  variable_sym = rlang::sym(variable)

  if(all(is.na(qc_types))){
    qc_types <- intersect(data$dataset$qc_type, c("SPL", "TQC", "BQC", "HQC", "MQC", "LQC", "QC", "NIST", "LTR", "PBLK", "SBLK"))
  }

  d_filt <- get_dataset_subset(
    data,
    filter_data = FALSE,
    qc_types = qc_types,
    include_qualifier = include_qualifier,
    include_istd = TRUE,
    include_feature_filter = include_feature_filter,
    exclude_feature_filter = exclude_feature_filter
  )


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

  df <- d_filt |>
    dplyr::select("feature_id", "qc_type", "batch_id", "is_istd", variable) 
  if(only_istd) df <- df |> filter(.data$is_istd)

  df$qc_type <- factor(df$qc_type, levels = c("PBLK", "TQC", "BQC", "LQC","MQC","HQC", "SPL", "NIST", "LTR"))


  grp <- if(batchwise_normalization) c("feature_id", "batch_id") else "feature_id"
  df_std <- df |>
    group_by(across(all_of(grp))) |>
    dplyr::mutate(scaled_intensity = !!variable_sym / mean(!!variable_sym, na.rm = TRUE) * 100) |> 
    drop_na(.data$scaled_intensity)

  ggplot2::ggplot(df_std, ggplot2::aes(
    x = .data$feature_id,
    y = .data$scaled_intensity,
    color = .data$qc_type,
    fill = .data$qc_type,
  )) +
    ggbeeswarm::geom_quasirandom(
      dodge.width = dodge_width,
      alpha = point_alpha,
      size = point_size,
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = .data$qc_type),
      position = ggplot2::position_dodge(width = dodge_width),
      alpha = box_alpha,
      width = 0.6,
      outlier.shape = NA,
      linewidth = box_linewidth,
    ) +
    ggplot2::geom_hline(yintercept = 100, linewidth = 0.5, color = "grey80", linetype = "dashed") +
    # ggplot2::labs(
    #   x = NULL,
    #   y = ""
    # ) +
    ggplot2::scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
    ggplot2::scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
    ggplot2::coord_cartesian(ylim = y_lim, expand = FALSE) +
    ggplot2::theme_bw(base_size = font_base_size) +
    ylab("Standardized Intensity (% of uncorrected)") +
    xlab("Internal Standard") +
    theme(
      axis.text.x = ggplot2::element_text(angle = angle_x, hjust = 1),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      #panel.grid.major.x = element_blank(), #element_line(colour = "#bdbdbd", linetype = "dotted", size = .5),
      panel.grid.minor.x = element_blank(), #element_line(colour = "grey88", linetype = "dotted", size = .5),
      axis.text.y = element_text(size = font_base_size),
      #axis.text.x = element_text(size = font_base_size, angle = x_text_angle, vjust = 0.5, hjust = x_text_just),
      axis.title = element_text(size = font_base_size * 1, face = "plain"),
      panel.border = element_rect(linewidth = 0.7, color = "grey40"),
      legend.position = "inside",
      legend.direction = "horizontal",         # vertical layout
      legend.text = element_text(size = font_base_size*0.8), # text size
      legend.title = element_text(size = font_base_size*0.8),
      legend.key.size = unit(font_base_size *0.8, "pt"),   # box size
      legend.position.inside = c(0.6, 0.1)
    )

}
