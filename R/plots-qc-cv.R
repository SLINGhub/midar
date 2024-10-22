#' Comparison of CV values before and after normalization
#' @param data MidarExperiment object
#' @param filter_data Use all (default) or qc-filtered data
#' @param qc_type QC type to use for plot.
#' @param var_before Variable of QC metrics table used for x axis
#' @param var_after Variable of QC metrics table used for y axis
#' @param cv_treshhold Threshold values shown in plot as dashed lines
#' @param only_quantifier Show only quantifier
#' @param xlim Lower and upper limits of the x axis as vector
#' @param ylim Lower and upper limits of the y axis as vector
#' @param ncol Number of facet columns
#' @param point_size size of points. Default is `2`
#' @param point_alpha transparency of points
#' @param scale_factor Scale factor for text labels
#' @param font_base_size Base font size for plot text elements
#' @return ggplot2 object
#' @export


qc_plot_normalization_cv <- function(data,
                           filter_data,
                           qc_type,
                           var_before = c("intensity"),
                           var_after = c("norm_intensity"),
                           cv_treshhold = 25,
                           only_quantifier = TRUE,
                           xlim = c(0, NA),
                           ylim = c(0, NA),
                           ncol = 5,
                           point_size = 1,
                           point_alpha = 0.5,
                           scale_factor = 1,
                           font_base_size = 8) {

  var_before <- str_remove(var_before, "feature_")
  rlang::arg_match(var_before, c("intensity", "norm_intensity", "conc"))

  var_after <- str_remove(var_after, "feature_")
  rlang::arg_match(var_after, c("intensity", "norm_intensity", "conc"))

  rlang::arg_match(qc_type,  c("SPL", "BQC", "TQC", "NIST", "LTR"))

  x_variable <- stringr::str_c(var_before, "_cv_", qc_type)
  y_variable <- stringr::str_c(var_after, "_cv_", qc_type)

  qc_plot_x_vs_y(data = data,
                             filter_data = filter_data,
                             x_variable = x_variable,
                             y_variable = y_variable,
                             cv_treshhold = cv_treshhold,
                             only_quantifier = only_quantifier,
                             xlim = xlim,
                             ylim = ylim,
                             ncol = ncol,
                             point_size = point_size,
                             point_alpha = point_alpha,
                             scale_factor = scale_factor,
                             font_base_size = font_base_size)

}

#' Contrast two variables from QC metrics table for all features per feature class
#' @param data MidarExperiment object
#' @param filter_data Use all (default) or qc-filtered data
#' @param x_variable Variable of QC metrics table used for x axis
#' @param y_variable Variable of QC metrics table used for y axis
#' @param cv_treshhold Treshhold values shown in plot as dashed lines
#' @param only_quantifier Show only quantifier
#' @param xlim Lower and upper limits of the x axis as vector
#' @param ylim Lower and upper limits of the y axis as vector
#' @param ncol Number of facet columns
#' @param point_size size of points. Default is `2`
#' @param point_alpha transparency of points
#' @param scale_factor Scale factor for text labels
#' @param font_base_size Base font size for plot text elements
#' @return ggplot2 object
#' @export


qc_plot_x_vs_y <- function(data,
                        filter_data,
                        x_variable,
                        y_variable,
                        cv_treshhold = 25,
                        only_quantifier = TRUE,
                        xlim = c(0, NA),
                        ylim = c(0, NA),
                        ncol = 5,
                        point_size = 1,
                        point_alpha = 0.5,
                        scale_factor = 1,
                        font_base_size = 8) {


  # Prepare data
  if(nrow(data@metrics_qc) == 0) {
    cli::cli_alert("No QC metrics available et. Please run `qc_calc_metrics()`, or set qc filters first.")
    return(NULL)
  }

  d_qc <- data@metrics_qc |> filter(.data$valid_feature)
  if (filter_data){
    if(data@is_filtered)
      d_qc <- d_qc |> filter(.data$valid_feature, .data$qc_pass)
    else {
      cli_abort(cli::col_red("Data has not yet been qc-filtered. Apply filter or set`use_filter_data = FALSE`."))
    }
  } else {
    d_qc <- d_qc |> filter(.data$valid_feature)
  }


  if (only_quantifier) d_qc <- d_qc |> filter(.data$is_quantifier, .data$valid_feature)

  x_sym <- rlang::sym(x_variable)
  y_sym <- rlang::sym(y_variable)

  if (!c("feature_class") %in% names(d_qc) || all(is.na(d_qc$feature_class)))
    cli::cli_abort("This function currently only works with data where `feature_class` has been defined.
                   Please define in metadata, or try `lipidomics_get_lipid_class_names()` to get the `feature_class`.")

  # get max value for the pair for each lipid class (so that 45deg line will be shown in the square plots)
  d_qc <- d_qc |>
    group_by(.data$feature_class) |>
    mutate(xy_max = max(!!sym(x_variable), !!sym(y_variable), na.rm = TRUE)) |>
    ungroup()


  # Use geom_polygon for coloring the two areas
  g <- ggplot(data = d_qc, aes(x = !!x_sym, y = !!y_sym)) +
    facet_wrap(vars(.data$feature_class), scales = "free", ncol = ncol) +
    geom_point(fill = "blue",
               size = point_size,
               shape = 21,
               alpha = 0.7,
               stroke = scale_factor/2,
               na.rm = TRUE) +
    geom_point(aes(x = .data$xy_max, y = .data$xy_max),
               fill = "darkblue",
               shape = ".",
               size = 1)

  g <- g +
    scale_x_continuous(limits = xlim, expand = expansion(mult = c(0, 0.2))) +
    scale_y_continuous(limits = ylim, expand = expansion(mult = c(0, 0.2))) +
    aes(ymin = 0, xmin = 0) +
    geom_abline(intercept = 0, slope = 1, linewidth = scale_factor * 0.3, color = "#ed5578") +
    # scale_colour_brewer(palette = "Set1") +
    #scale_shape_manual(na.value = NA, values = c(15, 0, 4), drop = FALSE, name = "Transition") +
    theme_bw(base_size = font_base_size * scale_factor) +
    theme(
      plot.title = element_text(size = 8 * scale_factor, face = "bold"),
      strip.text = element_text(size = 8 * scale_factor, margin = margin(1, 1, 0, 0), lineheight = 0),
      strip.background = element_rect(size = 0.0001),
      axis.text = element_text(size = font_base_size * scale_factor * 0.8, face = NULL),
      axis.title = element_text(size = font_base_size * scale_factor, face = "bold"),
      panel.grid = element_line(size = 0.001),
      strip.switch.pad.wrap = unit(1, "mm"),
      legend.position = "right"
    ) +
    geom_hline(yintercept = cv_treshhold, linetype = "dashed", color = "darkgreen") +
    geom_vline(xintercept = cv_treshhold, linetype = "dashed", color = "darkgreen")

  # if (with_histogram) {
  #   g <- g + theme(legend.position = c(0.9, 0.3))
  #   g <- ggExtra::ggMarginal(g, type = "histogram", margins = "y")
  # }
  return(g)
}
