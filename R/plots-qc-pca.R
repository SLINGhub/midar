#' PCA Plot for Quality Control
#'
#' @description Generates a Principal Component Analysis (PCA) plot for
#' visualizing samples including quality control (QC) samples. This function
#' provides options for filtering data, applying transformations, and customizing
#' visual elements to enhance the visualization.
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
#' @param ellipse_variable String specifying which sample variable to show
#' as ellipses. Must be one of: "none", "qc_type", "batch_id".
#' "none" omits ellipses.
#' @param ellipse_levels A character vector specifying the levels of
#' `ellipse_variable` to display as ellipses.
#' @param pca_dims A numeric vector of length 2 indicating the PCA dimensions
#' to plot. Default is c(1, 2).
#' @param log_transform A logical value indicating whether to log-transform
#' the data before the PCA. Default is `TRUE`.
#'
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
#' @param show_labels A logical value indicating whether to show analysis_id
#' labels for points outside k * MAD of the selected PCA dimensions. Default
#' is `TRUE`.
#' @param labels_threshold_mad A numeric value determining the threshold
#' for showing labels based on the median absolute deviation (MAD). Default
#' is 3. Set to `NULL` to suppress labels.
#' @param shared_labeltext_hide A character string representing text shared
#' across labels to be hidden (case-sensitive). If this results in
#' non-unique analysis_id's, an error will be raised.
#' @param label_font_size Number indicating the font size for labels in 'mm'.
#' Note the unit is different from font_base_size that is in 'pt'.
#'
#' @param point_size A numeric value indicating the size of points in
#' millimeters. Default is 2.
#' @param point_alpha A numeric value indicating the transparency of
#' points (0-1). Default is 0.5.
#' @param font_base_size A numeric value indicating the base font size for
#' plot text elements. Default is 8.
#' @param ellipse_confidence_level A numeric value indicating the confidence level
#' for the ellipses. Default is 0.95.
#' @param ellipse_linewidth A numeric value indicating the line width of the
#' ellipses. Default is 1.
#' @param ellipse_fill A logical value indicating whether to fill the ellipses.
#' @param ellipse_fillcolor A vector specifying the fill colors for ellipse corresponding
#'   to different `ellipse_variable` levels. This can be either an unnamed vector or a named
#'   vector, with names corresponding to leves in `ellipse_variable`. Unused fill colors will be ignored.
#'   Default is `NA` which corresponds to the default fill colors in case of
#'   `ellipse_variable = qc_type`, and to automatically generated colors otherwise.
#' @param ellipse_alpha A numeric value indicating the transparency of the
#' ellipse fill (0-1). Default is 0.3.

#'
#' @return A ggplot object with the PCA plot
#'
#' @export
plot_pca <- function(data = NULL,
                    variable,
                    qc_types = NA,
                    ellipse_variable = "qc_type",
                    ellipse_levels = NA,
                    pca_dims = c(1,2),

                    log_transform = TRUE,

                    filter_data = FALSE,
                    include_qualifier = FALSE,
                    include_istd = FALSE,
                    include_feature_filter = NA,
                    exclude_feature_filter = NA,
                    min_median_value = NA,

                    show_labels = TRUE,
                    labels_threshold_mad = 3,
                    shared_labeltext_hide = NA,
                    label_font_size = 3,

                    point_size = 2,
                    point_alpha = 0.8,
                    font_base_size = 8,

                    ellipse_confidence_level = 0.95,
                    ellipse_linewidth = 1,
                    ellipse_fill = TRUE,
                    ellipse_fillcolor = NA,
                    ellipse_alpha = 0.1) {


  # Check and define arguments

  check_data(data)
  variable <- str_remove(variable, "feature_")
  rlang::arg_match(variable, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt"))
  variable <- stringr::str_c("feature_", variable)
  check_var_in_dataset(data@dataset, variable)
  variable_sym = rlang::sym(variable)

  # Match qc_type TODO
  #rlang::arg_match(qc_types, c("SPL", "BQC", "TQC", "NIST", "LTR", "PQC", "SST", "RQC"))

  rlang::arg_match(ellipse_variable, c("none", "qc_type","batch_id"))
  ellipse_variable_sym = rlang::sym(ellipse_variable)

  PCx <- rlang::sym(paste0(".fittedPC", pca_dims[1]))
  PCy <- rlang::sym(paste0(".fittedPC", pca_dims[2]))


  if(show_labels) check_installed("ggrepel")
  if(ellipse_variable != "none") check_installed("ggnewscale")

  if(all(is.na(qc_types))){
    qc_types <- intersect(data$dataset$qc_type, c("SPL", "TQC", "BQC", "TQC", "HQC", "MQC", "LQC", "NIST", "LTR"))
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
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_id", values_from = all_of(variable))


  # if(!all(d_filt |> pull(analysis_id) == d_metadata |> pull(AnalyticalID))) cli::cli_abort("Data and Metadata missmatch")

  # ToDo: warning when rows/cols with NA
  d_clean <- d_wide |>
    filter(if_any(dplyr::where(is.numeric), ~ !is.na(.))) |>
    dplyr::select(where(~ !any(is.na(.) | is.nan(.) | is.infinite(.) | . <= 0)))

  n_removed <- ncol(d_wide) - ncol(d_clean)
  if(n_removed > 0)
    cli_alert_warning(col_yellow("{n_removed} features contained missing or non-numeric values and were exluded."))

  d_metadata <- d_filt |>
    dplyr::select("analysis_id", "qc_type", "batch_id") |>
    dplyr::distinct() |>
    dplyr::right_join(d_clean |> dplyr::select("analysis_id") |> distinct(), by = c("analysis_id"))

  m_raw <- d_clean |>
    tibble::column_to_rownames("analysis_id") |>
    as.matrix()



  if (log_transform) m_raw <- log2(m_raw)
  # get pca result with annotation
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)
  pca_annot <- pca_res |>
    broom::augment(d_metadata)

  pca_annot$qc_type <- droplevels(factor(pca_annot$qc_type, levels = c("SPL", "UBLK", "SBLK", "TQC", "BQC", "RQC", "LTR", "NIST", "PBLK")))
  pca_annot <- pca_annot |>
    dplyr::arrange(.data$qc_type)

  pca_contrib <- pca_res |> broom::tidy(matrix = "eigenvalues")



  if (!is.null(labels_threshold_mad) & !is.na(labels_threshold_mad)) {
   d_outlier <- pca_annot |> filter(abs(!!PCx) > (median(!!PCx) + labels_threshold_mad * mad(!!PCx)) |
                                      abs(!!PCy) > (median(!!PCy) + labels_threshold_mad * mad(!!PCy)))
   }

  pca_annot <- pca_annot |>
    select("analysis_id":stringr::str_c(".fittedPC", max(pca_dims))) |>
    left_join(d_outlier |> select("analysis_id"), by = "analysis_id", keep = TRUE, suffix = c("", "_outlier"))

  if (!is.na(shared_labeltext_hide)) {
    pca_annot <- pca_annot |> mutate(analysis_id_outlier = str_remove(.data$analysis_id_outlier, shared_labeltext_hide))
    if(any(duplicated(pca_annot$analysis_id_outlier, incomparables=NA)))
      cli_abort(cli::col_red("`shared_labeltext_hide` setting causes duplicate labels. Please adjust or set to `NA` to show full labels."))
  }

  #Check if ellipse_fillcolor is provided or is NA
  if(ellipse_variable != "none"){

    if(all(is.na(ellipse_levels))){
      pca_annot_ellipses <- pca_annot
    } else {
      if(length(setdiff(ellipse_levels, unique(pca_annot[[ellipse_variable]])))>0){
        cli_abort(col_red("One or more levels in `ellipse_levels` are not present in `{ellipse_variable}`. Please verify the levels and `ellipse_variable`."))
      }
      pca_annot_ellipses <- pca_annot |> filter(.data[[ellipse_variable]] %in% ellipse_levels)
    }

    if (is.null(ellipse_fillcolor) || length(ellipse_fillcolor) == 0 || all(is.na(ellipse_fillcolor))) {
      # If no ellipse_fillcolor is provided (NA or NULL), generate a discrete color scale
        n_col <- if(ellipse_variable == "qc_type") length(unique(d_filt$qc_type)) else length(unique(d_filt$batch_id))
        if (ellipse_variable == "qc_type") {
          ellipse_fillcolor <- pkg.env$qc_type_annotation$qc_type_col
          # Make the color for SPL lighter for fill
          if (!is.null(ellipse_fillcolor[["SPL"]]) && ellipse_fillcolor[["SPL"]] == "#8e9b9e") {
            ellipse_fillcolor[["SPL"]] <- "#cad6d9"
          }
        } else if (ellipse_variable == "batch_id") {
          ellipse_fillcolor <- scales::hue_pal()(n_col)
        } else {
          ellipse_fillcolor <-ellipse_fillcolor
        }
    } else {
      # If ellipse_fillcolor is provided, check if it has enough colors
      num_levels <- length(unique(pca_annot_ellipses[[ellipse_variable]]))
      if (length(ellipse_fillcolor) < num_levels) {
        cli::cli_abort(
          cli::col_red("Insufficient colors in `ellipse_fillcolor`. Provide at least {num_levels} unique colors for the number of {ellipse_variable}"))
      }
    }
}


  p <- ggplot(
    data = pca_annot,
    mapping = aes(
      x = !!sym(paste0(".fittedPC", pca_dims[1])),
      y = !!sym(paste0(".fittedPC", pca_dims[2])))
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.5, color = "grey80", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.5, color = "grey80", linetype = "dashed")

  if(show_labels){
        p <- p +
          ggrepel::geom_text_repel(aes(label = .data$analysis_id_outlier), size = label_font_size, na.rm = TRUE, seed = 1237)


  }

  if(ellipse_variable != "none"){
    p <- p +
      suppressWarnings(
        stat_ellipse(
          data = pca_annot_ellipses,
          aes(color = !!sym(ellipse_variable),
              fill = if (ellipse_fill) !!sym(ellipse_variable) else NA),
          geom = "polygon",
          level = ellipse_confidence_level,
          alpha = ellipse_alpha,
          linewidth = ellipse_linewidth,
          na.rm = TRUE)
        ) +
      scale_fill_manual(
        name = "Batch Id",
        values = ellipse_fillcolor,
        drop = TRUE,
        na.value = "transparent") +
      scale_color_manual(
        name = "Batch Id",
        values = ellipse_fillcolor,
        drop = TRUE) +
      ggplot2::guides(
        fill = if (ellipse_fill) {
          ggplot2::guide_legend(
            title = ellipse_variable,
            override.aes = list(size = 1, alpha = ellipse_alpha)
          )
        } else {
          "none"
        },
        color = if (ellipse_fill) {
          ggplot2::guide_legend(
            title = ellipse_variable
            #override.aes = list(size = 1, alpha = ellipse_alpha)
          )
        } else {
          ggplot2::guide_legend(
            title = ellipse_variable
            #override.aes = list(size = 1, alpha = 0.0)
          )
        }
      )




  }

  p <- p +
    ggnewscale::new_scale_fill() +
    ggnewscale::new_scale_color() +
    ggplot2::geom_point(size = point_size,
                        aes(color = .data$qc_type, shape = .data$qc_type, fill = .data$qc_type),
                        alpha = point_alpha) +
    ggplot2::scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
    ggplot2::scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
    ggplot2::scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)

  p <- p +
    ggplot2::theme_bw(base_size = font_base_size) +
    ggplot2::xlab(glue::glue("PC{pca_dims[1]} ({round(pca_contrib[[pca_dims[1],'percent']]*100,1)}%)")) +
    ggplot2::ylab(glue::glue("PC{pca_dims[2]} ({round(pca_contrib[[pca_dims[2],'percent']]*100,1)}%)")) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(linewidth = 1, color = "grey40"),
      axis.text.x = ggplot2::element_text(size = font_base_size, face = NULL),
      axis.text.y = ggplot2::element_text(size = font_base_size, face = NULL),
      axis.title.x = ggplot2::element_text(size = font_base_size * 1.2, face = NULL),
      axis.title.y = ggplot2::element_text(size = font_base_size * 1.2, face = NULL),
      aspect.ratio = 1
    )

  txt <- if(log_transform) "log2-transformed" else "raw"
  cli::cli_alert_info(col_green("The PCA was calculated based on `{variable}` values of {length(unique(d_filt$feature_id))} features."))

  p
}

#' Plot PCA loadings
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
#' @param pca_dims A numeric vector indicating for which PCA dimensions
#' to the loadings should be shown. Default is c(1, 2, 3, 4).
#' @param log_transform A logical value indicating whether to log-transform
#' the data before the PCA. Default is `TRUE`.
#' 
#' @param top_n Number of top features with highest absolute loading that will be to shown for each PC dimension. Default is 30.
#' @param vertical_bars Show vertical bars instead of horizontal bars in the plot. Default is `FALSE`.
#' @param abs_loading Show absolute loading values instead of signed loadings. Default is `TRUE`.
#'
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
#' @param font_base_size A numeric value indicating the base font size for
#' plot text elements. Default is 7.
#'
#' @returns ggplot object with PCA loadings plot
#'
#' @export
plot_pca_loading <- function(data = NULL,
                            variable,
                            qc_types = NA,
                            pca_dims = c(1, 2, 3, 4),
                            log_transform = TRUE,
                            top_n = 30,
                            vertical_bars = FALSE,
                            abs_loading = TRUE,
                            filter_data = FALSE,
                            include_qualifier = FALSE,
                            include_istd = FALSE,
                            include_feature_filter = NA,
                            exclude_feature_filter = NA,
                            min_median_value = NA,
                            font_base_size = 7) {
 # Check and define arguments

  check_data(data)
  variable <- str_remove(variable, "feature_")
  rlang::arg_match(variable, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt"))
  variable <- stringr::str_c("feature_", variable)
  check_var_in_dataset(data@dataset, variable)
  variable_sym = rlang::sym(variable)


  if(all(is.na(qc_types))){
    qc_types <- intersect(data$dataset$qc_type, c("SPL", "TQC", "BQC", "TQC", "HQC", "MQC", "LQC", "NIST", "LTR"))
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
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_id", {{ variable }}) |>
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_id", values_from = {{ variable }})

  d_filt <- d_filt |>
    tibble::column_to_rownames("analysis_id") |>
    dplyr::select(where(~ !any(is.na(.))))

  m_raw <- d_filt |>
    filter(if_any(dplyr::where(is.numeric), ~ !is.na(.))) |>
    dplyr::select(where(~ !any(is.na(.) | is.nan(.) | is.infinite(.) | . <= 0))) |>
    as.matrix()

  if (log_transform) m_raw <- log2(m_raw)
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)

  d_loading <- pca_res |>
    broom::tidy(matrix = "rotation") |>
    tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") |>
    dplyr::rename(feature_name = .data$column)

  d_loadings_selected <- d_loading |>
    tidyr::pivot_longer(cols = -.data$feature_name, names_to = "PC", values_to = "Value") |>
    dplyr::mutate(PC = as.numeric(stringr::str_remove(.data$PC, "PC"))) |>
    filter(.data$PC %in% pca_dims)


  d_loadings_selected <- d_loadings_selected |>
    dplyr::rowwise() |>
    dplyr::mutate(
      direction = if_else(.data$Value < 0, "neg", "pos"),
      Value = dplyr::if_else(abs_loading, abs(.data$Value), .data$Value),
      abs_value = abs(.data$Value)
    ) |>
    group_by(.data$PC)

  if (vertical_bars) {
    d_loadings_selected <- d_loadings_selected |> dplyr::arrange(.data$abs_value)
  } else {
    d_loadings_selected <- d_loadings_selected |> dplyr::arrange(desc(.data$abs_value))
  }

  d_loadings_selected <- d_loadings_selected |>
    dplyr::slice_max(order_by = .data$abs_value, n = top_n) |>
    ungroup() |>
    tidyr::unite("Feature", .data$feature_name, .data$PC, remove = FALSE) |>
    mutate(
      PC = as.factor(.data$PC),
      Feature = forcats::fct_reorder(.data$Feature, .data$abs_value)
    )


  p <- ggplot(d_loadings_selected, ggplot2::aes(x = .data$Feature, y = .data$Value, color = .data$direction, fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(ggplot2::vars(.data$PC), scales = "free", ncol = ifelse(vertical_bars, 1, length(pca_dims))) +
    ggplot2::scale_x_discrete(labels = d_loadings_selected$feature_name, breaks = d_loadings_selected$Feature) +
    ggplot2::scale_color_manual(values = c("neg" = "lightblue", "pos" = "#FF8C00")) +
    ggplot2::scale_fill_manual(values = c("neg" = "lightblue", "pos" = "#FF8C00")) +
    ggplot2::labs(x = "Feature", y = "Loading") +
    ggplot2::theme_bw(base_size = 8)

  if (!vertical_bars) {
    p <- p +
      ggplot2::coord_flip()
  } else {
    p <- p +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1)) +
      ggplot2::scale_x_discrete(limits = rev)
  }

  p +
    theme_bw(base_size = font_base_size) +
      theme(
        plot.title = element_text(size = font_base_size, face = "bold"),
        strip.text = ggplot2::element_text(size = font_base_size, face = "bold"),
        axis.text.x = element_text(size = font_base_size * 0.6),
        axis.text.y = element_text(size = font_base_size * 0.6),
        axis.title = element_text(size = font_base_size, face = "bold"),
        panel.grid = element_line(linewidth = 0.001),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(linewidth = 0.2),
        strip.background = ggplot2::element_rect(linewidth = 0.0001, fill = "#00283d"),
        strip.text.x = ggplot2::element_text(color = "white"),
        #strip.switch.pad.wrap = ggplot2::unit(1, "mm"),
        panel.border = element_rect(linewidth = 0.5, color = "grey40"),
        legend.position = "right"
      )

}


# plot_pca_pairs <- function(data = NULL, variable, dim_range = c(1, 8), use_filter_data, qc_types = c("BQC", "TQC", "NIST", "LTR", "SPL"), log_transform = TRUE, grouping = "qc_type", exclude_istds = TRUE, sliding = FALSE, ncol = 3,
#                            point_size = 0.5, fill_alpha = 0.1, legend_pos = "right") {
#   check_data(data)
#   if(use_filter_data)
#     if(data@is_filtered)
#       d_wide <- data@dataset_filtered
#   else
#     cli_abort(cli::col_red("Data has not yet been qc-filtered. Use `use_filter_data = FALSE` to use unfiltered data."))
#   else
#     d_wide <- data@dataset
#
#   # TODO: (IS as criteria for ISTD.. dangerous...
#   if (exclude_istds) d_wide <- d_wide |> filter(!.data$is_istd, !str_detect(.data$feature_id, "\\(IS\\)| ISTD")) # !stringr::str_detect(.data$feature_id, "\\(IS")
#
#   d_wide <- d_wide |>
#     filter(.data$qc_type %in% qc_types, .data$is_quantifier) |>
#     dplyr::select("analysis_id", "qc_type", "batch_id", "feature_id", {{ variable }})
#
#   d_filt <- d_wide |>
#     tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_id", values_from = {{ variable }})
#
#
#   # if(!all(d_filt |> pull(analysis_id) == d_metadata |> pull(AnalyticalID))) cli::cli_abort("Data and Metadata missmatch")
#
#   # ToDo: warning when rows/cols with NA
#   d_clean <- d_filt |>
#     filter(if_any(dplyr::where(is.numeric), ~ !is.na(.))) |>
#     dplyr::select(where(~ !any(is.na(.) | is.nan(.) | is.infinite(.) | . <= 0)))
#
#
#   d_metadata <- d_wide |>
#     dplyr::select("analysis_id", "qc_type", "batch_id") |>
#     dplyr::distinct() |>
#     dplyr::right_join(d_clean |> dplyr::select("analysis_id") |> distinct(), by = c("analysis_id"))
#
#
#   m_raw <- d_clean |>
#     tibble::column_to_rownames("analysis_id") |>
#     as.matrix()
#
#
#
#   if (log_transform) m_raw <- log2(m_raw)
#
#   # get pca result with annotation
#   pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)
#   pca_annot <- pca_res |> broom::augment(d_metadata)
#   pca_contrib <- pca_res |> broom::tidy(matrix = "eigenvalues")
#
#   dim_range <- seq(dim_range[1]:dim_range[2])
#   if (sliding) {
#     dim_range <- dim_range[-length(dim_range)]
#   } else {
#     dim_range <- dim_range[seq(1, length(dim_range), 2)]
#   }
#
#   plot_list <- list()
#
#   j <- 1
#
#   for (i in dim_range) {
#     p <- ggplot(data = pca_annot, aes(
#       x = !!sym(paste0(".fittedPC", i)),
#       y = !!sym(paste0(".fittedPC", i + 1)),
#       color = grouping,
#       shape = grouping,
#       fill = grouping
#     )) +
#       geom_hline(yintercept = 0, size = 0.2, color = "grey80") +
#       geom_vline(xintercept = 0, size = 0.2, color = "grey80") +
#       ggplot2::stat_ellipse(geom = "polygon", level = 0.95, alpha = fill_alpha, size = 0.2) +
#       geom_point(size = point_size)
#
#     p <- p +
#       scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
#       scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
#       scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)
#
#
#     p <- p +
#       ggplot2::theme_light(base_size = 6) +
#       ggplot2::xlab(glue::glue("PC{i} ({round(pca_contrib[[i,'percent']]*100,1)}%)")) +
#       ggplot2::ylab(glue::glue("PC{i+1} ({round(pca_contrib[[i+1,'percent']]*100,1)}%)")) +
#       ggplot2::theme(
#         panel.grid = ggplot2::element_line(size = 0.2, color = "grey95"),
#         panel.border = ggplot2::element_rect(size = .5, color = "grey70"),
#         aspect.ratio = 1
#       )
#
#     plot_list[[j]] <- p + ggplot2::theme(legend.position = "none")
#     j <- j + 1
#   }
#   # lg <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = legend_pos,
#   #                                                  legend.margin = ggplot2::margin(c(0,0,0,0))))
#   # pl1 <- cowplot::plot_grid( plotlist = plot_list, ncol = ncol)
#   # print(cowplot::plot_grid( pl1,lg, ncol = 1, rel_heights = c(1,0.2)))
# }
#
# plot_pca_loading_coord <- function(data = NULL, variable, log_transform, dim_x, dim_y, top_n, text_size = 1, fill_alpha = 0.1) {
#   check_data(data)
#   PCx <- rlang::sym(paste0("PC", dim_x))
#   PCy <- rlang::sym(paste0("PC", dim_y))
#
#   d_wide <- data@dataset |>
#     filter(.data$qc_type %in% c("BQC", "TQC", "NIST", "LTR", "SPL"), !stringr::str_detect(.data$feature_id, "\\(IS")) |>
#     dplyr::select("analysis_id", "qc_type", "batch_id", "feature_id", {{ variable }})
#
#   d_filt <- d_wide |>
#     tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_id", values_from = {{ variable }})
#
#   m_raw <- d_filt |>
#     dplyr::select(where(~ !any(is.na(.)))) |>
#     tibble::column_to_rownames("analysis_id") |>
#     as.matrix()
#
#   if (log_transform) m_raw <- log2(m_raw)
#   pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)
#
#   d_loading <- pca_res |>
#     broom::tidy(matrix = "rotation") |>
#     tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value")
#
#   d_top_loadings <- d_loading |>
#     dplyr::select("column", !!PCx, !!PCy) |>
#     mutate(vl = (!!PCx)^2 + (!!PCy)^2) |>
#     dplyr::arrange(dplyr::desc(.data$vl)) |>
#     head(top_n)
#   x_max <- max(abs(d_top_loadings[, glue::glue("PC{dim_x}")]))
#   y_max <- max(abs(d_top_loadings[, glue::glue("PC{dim_y}")]))
#
#   # define arrow style for plotting
#   arrow_style <- ggplot2::arrow(
#     angle = 20, ends = "first", type = "closed", length = grid::unit(6, "pt")
#   )
#
#   p <- ggplot(d_top_loadings, aes(
#     x = !!sym(glue("PC{dim_x}")),
#     y = !!sym(glue("PC{dim_y}"))
#   )) +
#     ggplot2::geom_segment(xend = 0, yend = 0, arrow = arrow_style, color = "grey60", size = 0.3) +
#     ggplot2::geom_text(
#       aes(label = .data$column),
#       size = text_size,
#       hjust = 0, nudge_x = -0.01,
#       color = "#904C2F"
#     ) +
#     ggplot2::scale_x_continuous(limits = c(-x_max * 1.2, x_max * 1.2), expand = ggplot2::expansion(mult = .2)) +
#     ggplot2::scale_y_continuous(limits = c(-y_max * 1.2, y_max * 1.2), expand = ggplot2::expansion(mult = .2)) +
#     # coord_fixed(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) + # fix aspect ratio to 1:1
#     ggplot2::theme_light(base_size = 7) +
#     ggplot2::theme(aspect.ratio = 1)
#   p
# }
#
