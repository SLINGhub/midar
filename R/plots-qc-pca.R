#' PCA plot for QC
#' @param data MidarExperiment object
#' @param variable Which variable to use for the PCA. Must be any of "area", "height", "intensity", "response", "conc", "conc_raw", "rt", "fwhm".
#' @param use_filtered_data Use all (default) or qc-filtered data
#' @param pca_dim PCA dimensions to plot as a vector of length 2. Default is `c(1,2)`
#' @param qc_types qc types to plot. Default is `c("SPL", "BQC", "TQC", "NIST", "LTR")`
#' @param label_k_mad Show analysis_id label for points outside k * mad of any of two defined PCA dimensions. Default is `3`. Set to `NULL` to supress labels.
#' @param log_transform log-transform data before calculating PCA. Default is `TRUE`
#' @param remove_istds Exlude ISTD features. Default is `TRUE`
#' @param min_median_signal Minimum median signal across all samples from all selected qc types. `NA` (default) will not filter any features based on signal.
#' @param point_size size of points. Default is `2`
#' @param point_alpha transparency of points
#' @param ellipse_alpha transparency of ellipse fill
#' @param font_base_size Base font size for plot text elements
#' @param hide_label_text Hide shared text (case-senstive) in labels. If it results in non-unique analysis_id's then an error will be raised.
#'
#' @return ggplot2 object
#' @export
plot_pca_qc <- function(data, variable, use_filtered_data, pca_dim = c(1,2), qc_types = c("SPL", "BQC", "TQC", "NIST", "LTR"),
                        label_k_mad = 3, log_transform = TRUE, remove_istds = TRUE, min_median_signal = NA, point_size = 2, point_alpha = 0.7,
                        ellipse_alpha = 0.8, font_base_size = 8, hide_label_text = NA) {

  variable <- str_remove(variable, "feature_")
  rlang::arg_match(variable, c("area", "height", "intensity", "response", "conc", "conc_raw", "rt", "fwhm"))

  variable <- stringr::str_c("feature_", variable)

  variable_sym = rlang::sym(variable)

  if(use_filtered_data)
    if(data@is_filtered)
      d_wide <- data@dataset_filtered
    else
      cli_abort(cli::col_red("Data has not yet been qc-filtered. Use `use_filtered_data = FALSE` to use unfiltered data."))
  else
    d_wide <- data@dataset


  PCx <- rlang::sym(paste0(".fittedPC", pca_dim[1]))
  PCy <- rlang::sym(paste0(".fittedPC", pca_dim[2]))

  # TODO: (IS as criteria for ISTD.. dangerous...
  if (remove_istds) d_wide <- d_wide |> filter(!.data$is_istd, !str_detect(.data$feature_id, "\\(IS\\)| ISTD")) # !stringr::str_detect(.data$feature_id, "\\(IS")

  d_wide <- d_wide |>
    filter(.data$qc_type %in% qc_types, .data$is_quantifier) |>
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_id", {{ variable }})

  if(!is.na(min_median_signal)){
    d_minsignal <- d_wide |>
      summarise(median_signal = median(!!variable_sym, na.rm = TRUE), .by = "feature_id") |>
      filter(.data$median_signal >= min_median_signal)
    d_wide <- d_wide |> semi_join(d_minsignal, by = "feature_id")
  }

  d_filt <- d_wide |>
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_id", values_from = variable )


  # if(!all(d_filt |> pull(analysis_id) == d_metadata |> pull(AnalyticalID))) cli::cli_abort("Data and Metadata missmatch")

  # ToDo: warning when rows/cols with NA
  d_clean <- d_filt |>
    filter(if_any(dplyr::where(is.numeric), ~ !is.na(.))) |>
    dplyr::select(where(~ !any(is.na(.) | is.nan(.) | is.infinite(.) | . <= 0)))


  d_metadata <- d_wide |>
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



  if (!is.null(label_k_mad) & !is.na(label_k_mad)) {
   d_outlier <- pca_annot |> filter(abs(!!PCx) > (median(!!PCx) + label_k_mad * mad(!!PCx)) |
                                      abs(!!PCy) > (median(!!PCy) + label_k_mad * mad(!!PCy)))
   }


  pca_annot <- pca_annot |>
    select("analysis_id":stringr::str_c(".fittedPC", max(pca_dim))) |>
    left_join(d_outlier |> select("analysis_id"), by = "analysis_id", keep = TRUE, suffix = c("", "_outlier"))

  if (!is.na(hide_label_text)) {
    pca_annot <- pca_annot |> mutate(analysis_id_outlier = str_remove(.data$analysis_id_outlier, hide_label_text))
    if(any(duplicated(pca_annot$analysis_id_outlier, incomparables=NA)))
      cli_abort(cli::col_red("Hiding text defined with `hide_label_text` would result in non-unique labels. Please modify or set to `NA` to show labels."))
  }

  p <- ggplot(data = pca_annot, aes_string(paste0(".fittedPC", pca_dim[1]),
    paste0(".fittedPC", pca_dim[2]),
    color = "qc_type",
    fill = "qc_type",
    shape = "qc_type",
    group = "qc_type"
  )) +
    ggplot2::geom_hline(yintercept = 0, size = 0.5, color = "grey80", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, size = 0.5, color = "grey80", linetype = "dashed") +
    suppressWarnings(ggplot2::stat_ellipse(data = pca_annot |> filter(.data$qc_type %in% c("BQC", "TQC", "SPL")), geom = "polygon", level = 0.95, alpha = ellipse_alpha, size = 0.3, na.rm = TRUE)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    #ggplot2::geom_text(hjust=0, vjust=0)
    ggrepel::geom_text_repel(aes(label = .data$analysis_id_outlier), size = 3, na.rm = TRUE)
  p <- p +
    ggplot2::scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
    ggplot2::scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
    ggplot2::scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)

  p <- p +
    ggplot2::theme_bw(base_size = font_base_size) +
    ggplot2::xlab(glue::glue("PC{pca_dim[1]} ({round(pca_contrib[[pca_dim[1],'percent']]*100,1)}%)")) +
    ggplot2::ylab(glue::glue("PC{pca_dim[2]} ({round(pca_contrib[[pca_dim[2],'percent']]*100,1)}%)")) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(size = 1, color = "grey40"),
      axis.text.x = ggplot2::element_text(size = font_base_size, face = NULL),
      axis.text.y = ggplot2::element_text(size = font_base_size, face = NULL),
      axis.title.x = ggplot2::element_text(size = font_base_size * 1.2, face = NULL),
      axis.title.y = ggplot2::element_text(size = font_base_size * 1.2, face = NULL),
      aspect.ratio = 1
    )

  p
}

plot_pca_pairs <- function(data, variable, dim_range = c(1, 8), use_filtered_data, qc_types = c("BQC", "TQC", "NIST", "LTR", "SPL"), log_transform = TRUE, grouping = "qc_type", remove_istds = TRUE, sliding = FALSE, ncol = 3,
                           point_size = 0.5, fill_alpha = 0.1, legend_pos = "right") {
  if(use_filtered_data)
    if(data@is_filtered)
      d_wide <- data@dataset_filtered
  else
    cli_abort(cli::col_red("Data has not yet been qc-filtered. Use `use_filtered_data = FALSE` to use unfiltered data."))
  else
    d_wide <- data@dataset

  # TODO: (IS as criteria for ISTD.. dangerous...
  if (remove_istds) d_wide <- d_wide |> filter(!.data$is_istd, !str_detect(.data$feature_id, "\\(IS\\)| ISTD")) # !stringr::str_detect(.data$feature_id, "\\(IS")

  d_wide <- d_wide |>
    filter(.data$qc_type %in% qc_types, .data$is_quantifier) |>
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_id", {{ variable }})

  d_filt <- d_wide |>
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_id", values_from = {{ variable }})


  # if(!all(d_filt |> pull(analysis_id) == d_metadata |> pull(AnalyticalID))) cli::cli_abort("Data and Metadata missmatch")

  # ToDo: warning when rows/cols with NA
  d_clean <- d_filt |>
    filter(if_any(dplyr::where(is.numeric), ~ !is.na(.))) |>
    dplyr::select(where(~ !any(is.na(.) | is.nan(.) | is.infinite(.) | . <= 0)))


  d_metadata <- d_wide |>
    dplyr::select("analysis_id", "qc_type", "batch_id") |>
    dplyr::distinct() |>
    dplyr::right_join(d_clean |> dplyr::select("analysis_id") |> distinct(), by = c("analysis_id"))


  m_raw <- d_clean |>
    tibble::column_to_rownames("analysis_id") |>
    as.matrix()



  if (log_transform) m_raw <- log2(m_raw)

  # get pca result with annotation
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)
  pca_annot <- pca_res |> broom::augment(d_metadata)
  pca_contrib <- pca_res |> broom::tidy(matrix = "eigenvalues")

  dim_range <- seq(dim_range[1]:dim_range[2])
  if (sliding) {
    dim_range <- dim_range[-length(dim_range)]
  } else {
    dim_range <- dim_range[seq(1, length(dim_range), 2)]
  }

browser()
  plot_list <- list()

  j <- 1

  for (i in dim_range) {
    p <- ggplot(data = pca_annot, aes_string(paste0(".fittedPC", i),
      paste0(".fittedPC", i + 1),
      color = grouping,
      shape = grouping,
      fill = grouping
    )) +
      geom_hline(yintercept = 0, size = 0.2, color = "grey80") +
      geom_vline(xintercept = 0, size = 0.2, color = "grey80") +
      ggplot2::stat_ellipse(geom = "polygon", level = 0.95, alpha = fill_alpha, size = 0.2) +
      geom_point(size = point_size)

    p <- p +
      scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
      scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
      scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)


    p <- p +
      ggplot2::theme_light(base_size = 6) +
      ggplot2::xlab(glue::glue("PC{i} ({round(pca_contrib[[i,'percent']]*100,1)}%)")) +
      ggplot2::ylab(glue::glue("PC{i+1} ({round(pca_contrib[[i+1,'percent']]*100,1)}%)")) +
      ggplot2::theme(
        panel.grid = ggplot2::element_line(size = 0.2, color = "grey95"),
        panel.border = ggplot2::element_rect(size = .5, color = "grey70"),
        aspect.ratio = 1
      )

    plot_list[[j]] <- p + ggplot2::theme(legend.position = "none")
    j <- j + 1
  }
  # lg <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = legend_pos,
  #                                                  legend.margin = ggplot2::margin(c(0,0,0,0))))
  # pl1 <- cowplot::plot_grid( plotlist = plot_list, ncol = ncol)
  # print(cowplot::plot_grid( pl1,lg, ncol = 1, rel_heights = c(1,0.2)))
}

plot_pca_loading_coord <- function(data, variable, log_transform, dim_x, dim_y, top_n, text_size = 1, fill_alpha = 0.1) {
  PCx <- rlang::sym(paste0("PC", dim_x))
  PCy <- rlang::sym(paste0("PC", dim_y))

  d_wide <- data@dataset |>
    filter(.data$qc_type %in% c("BQC", "TQC", "NIST", "LTR", "SPL"), !stringr::str_detect(.data$feature_id, "\\(IS")) |>
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_id", {{ variable }})

  d_filt <- d_wide |>
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_id", values_from = {{ variable }})

  m_raw <- d_filt |>
    dplyr::select(where(~ !any(is.na(.)))) |>
    tibble::column_to_rownames("analysis_id") |>
    as.matrix()

  if (log_transform) m_raw <- log2(m_raw)
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)

  d_loading <- pca_res |>
    broom::tidy(matrix = "rotation") |>
    tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value")

  d_top_loadings <- d_loading |>
    dplyr::select("column", !!PCx, !!PCy) |>
    mutate(vl = (!!PCx)^2 + (!!PCy)^2) |>
    dplyr::arrange(dplyr::desc(.data$vl)) |>
    head(top_n)
  x_max <- max(abs(d_top_loadings[, glue::glue("PC{dim_x}")]))
  y_max <- max(abs(d_top_loadings[, glue::glue("PC{dim_y}")]))

  # define arrow style for plotting
  arrow_style <- ggplot2::arrow(
    angle = 20, ends = "first", type = "closed", length = grid::unit(6, "pt")
  )

  p <- ggplot(d_top_loadings, aes_string(glue::glue("PC{dim_x}"), glue::glue("PC{dim_y}"))) +
    ggplot2::geom_segment(xend = 0, yend = 0, arrow = arrow_style, color = "grey60", size = 0.3) +
    ggplot2::geom_text(
      aes(label = .data$column),
      size = text_size,
      hjust = 0, nudge_x = -0.01,
      color = "#904C2F"
    ) +
    ggplot2::scale_x_continuous(limits = c(-x_max * 1.2, x_max * 1.2), expand = ggplot2::expansion(mult = .2)) +
    ggplot2::scale_y_continuous(limits = c(-y_max * 1.2, y_max * 1.2), expand = ggplot2::expansion(mult = .2)) +
    # coord_fixed(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) + # fix aspect ratio to 1:1
    ggplot2::theme_light(base_size = 7) +
    ggplot2::theme(aspect.ratio = 1)
  p
}

plot_pca_loading <- function(data, variable, log_transform, pc_dimensions, top_n, remove_istds, vertical_bars = FALSE, scale_pos_neg = FALSE, point_size = 2, fill_alpha = 0.1) {
  d_wide <- data@dataset_filtered |> filter(.data$qc_type %in% c("BQC", "TQC", "NIST", "LTR", "SPL"))

  if (remove_istds) d_wide <- d_wide |> filter(!.data$is_istd) # !stringr::str_detect(.data$feature_id, "\\(IS")

  d_filt <- d_wide |>
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_id", {{ variable }}) |>
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_id", values_from = {{ variable }})

  m_raw <- d_filt |>
    tibble::column_to_rownames("analysis_id") |>
    dplyr::select(where(~ !any(is.na(.)))) |>
    as.matrix()

  if (log_transform) m_raw <- log2(m_raw)
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)

  d_loading <- pca_res |>
    broom::tidy(matrix = "rotation") |>
    tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") |>
    dplyr::rename(feauture_name = .data$column)

  d_loadings_selected <- d_loading |>
    tidyr::pivot_longer(cols = -.data$feauture_name, names_to = "PC", values_to = "Value") |>
    dplyr::mutate(PC = as.numeric(stringr::str_remove(.data$PC, "PC"))) |>
    filter(.data$PC %in% pc_dimensions)


  d_loadings_selected <- d_loadings_selected |>
    dplyr::rowwise() |>
    dplyr::mutate(
      direction = if_else(.data$Value < 0, "neg", "pos"),
      Value = dplyr::if_else(!scale_pos_neg, abs(.data$Value), .data$Value),
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
    tidyr::unite("Feature", .data$feauture_name, .data$PC, remove = FALSE) |>
    mutate(
      PC = as.factor(.data$PC),
      Feature = forcats::fct_reorder(.data$Feature, .data$abs_value)
    )


  p <- ggplot(d_loadings_selected, ggplot2::aes(x = .data$Feature, y = .data$Value, color = .data$direction, fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(ggplot2::vars(.data$PC), scales = "free", ncol = ifelse(vertical_bars, 1, length(pc_dimensions))) +
    ggplot2::scale_x_discrete(labels = d_loadings_selected$feauture_name, breaks = d_loadings_selected$Feature) +
    ggplot2::scale_color_manual(values = c("neg" = "#75CEFF", "pos" = "#FFA166")) +
    ggplot2::scale_fill_manual(values = c("neg" = "#75CEFF", "pos" = "#FFA166")) +
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

  p
}
