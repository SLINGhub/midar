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
#' @param sort_by_corr A logical value indicating whether to sort the features in
#' the plot by correlation or alphabetically by feature ID. Default is `TRUE`.
#' @param rows_page Number of rows of plots per page.
#' @param cols_page Number of columns of plots per page.
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
#' @param output_pdf If `TRUE`, saves the generated plots as a PDF
#'   file. When `FALSE`, plots are directly plotted.
#' @param path The file path for saving the PDF. Must be defined if
#'   `output_pdf` is `TRUE`.
#' @param return_plots Logical. If `TRUE`, returns the plots as a list of
#'   `ggplot2` objects.
#' @param specific_page An integer specifying a specific page to plot. If
#'   `NA` (default), all pages are plotted.
#' @param page_orientation Orientation of the PDF paper: `"LANDSCAPE"` or
#'   `"PORTRAIT"`.
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
#' @param show_progress Logical. If `TRUE`, displays a progress bar during
#'   plot creation.
#'
#' @return A ggplot object showing scatter plots of highly correlated feature pairs.
#' Returns NULL if no correlations meet the threshold criteria.
#' @export

plot_feature_correlations <- function(data,
                              variable,
                              qc_types = NA,
                              cor_min,
                              cor_min_neg = -0.99,
                              log_scale = FALSE,
                              sort_by_corr = TRUE,
                              rows_page = 4,
                              cols_page = 5,
                              filter_data = FALSE,
                              include_qualifier = FALSE,
                              include_istd = FALSE,
                              include_feature_filter = NA,
                              exclude_feature_filter = NA,
                              min_median_value = NA,
                              output_pdf = FALSE,
                              path = NA,
                              specific_page = NA,
                              page_orientation = "LANDSCAPE",
                              return_plots = FALSE,
                              point_size = 1,
                              point_alpha = 0.8,
                              point_stroke = 0.3,
                              line_size = 0.5,
                              line_color = "orange",
                              line_alpha  = 0.5,
                              font_base_size = 8,
                              show_progress = TRUE) {

  check_data(data)

  variable <- str_remove(variable, "feature_")
  rlang::arg_match(variable, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm", "width", "symmetry"))
  variable <- stringr::str_c("feature_", variable)
  check_var_in_dataset(data@dataset, variable)
  variable_sym = rlang::sym(variable)



  if(all(is.na(qc_types))){
    qc_types <- intersect(data$dataset$qc_type, c("SPL", "TQC", "BQC", "HQC", "MQC", "LQC", "QC", "NIST", "LTR"))
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
  d_plot <- purrr::map_df(1:nrow(cor_matrix), function(i) {
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
  })



  d_plot$qc_type <- droplevels(factor(d_plot$qc_type, levels = c("SPL", "UBLK", "SBLK", "TQC", "BQC", "RQC", "LTR", "NIST", "PBLK")))
  d_plot <- d_plot |>
    dplyr::arrange(.data$qc_type)

  # Prepare PDF output
  if (output_pdf && !is.na(path)) { # nocov start
    path <- ifelse(stringr::str_detect(path, ".pdf"), path,
                   paste0(path, ".pdf"))
    if (page_orientation == "LANDSCAPE") {
      pdf(file = path, onefile = TRUE, paper = "A4r", useDingbats = FALSE,
          width = 28 / 2.54, height = 20 / 2.54)
    } else {
      pdf(file = path, onefile = TRUE, paper = "A4", useDingbats = FALSE,
          height = 28 / 2.54, width = 20 / 2.54)
    }
  } # nocov end

  # Determine the range of pages to generate
  if (!is.na(specific_page)) {
    total_pages <- ceiling(n_distinct(d_plot$pair) /
                             (cols_page * rows_page))
    if (specific_page > total_pages) {
      cli::cli_abort(col_red(
        "Selected page exceeds the total number of pages. Please select a page number between {.strong 1} and {.strong {total_pages}}."
      ))
    }
    page_range <- specific_page
  } else {
    page_range <- 1:ceiling(n_distinct(d_plot$pair) /
                              (cols_page * rows_page))
  }

  # Action text for progress output
  action_text <- if (output_pdf) "Saving plots to pdf" else "Generating plots"
  page_suffix <- if (max(page_range) > 1) {
    glue::glue("{max(page_range)} pages")
  } else {
    glue::glue("{max(page_range)} page")
  }
  progress_suffix <- if (show_progress) ":" else "..."

  # Output formatted message
  message(glue::glue("{action_text} ({page_suffix}){progress_suffix}"))

  # Initialize progress bar if requested
  if (show_progress) pb <- txtProgressBar(min = 0, max = max(page_range),
                                          width = 30, style = 3)

  p_list <- list()  # List to store plots for each page
  for (i in page_range) {
    p <- plot_feature_correlations_page(
      d_plot = d_plot,
      output_pdf = output_pdf,
      path = path,
      rows_page = rows_page,
      cols_page = cols_page,
      specific_page = i,
      sort_by_corr = sort_by_corr,
      log_scale = log_scale,
      point_size = point_size,
      point_alpha = point_alpha,
      point_stroke = point_stroke,
      line_size = line_size,
      line_color = line_color,
      line_alpha  = line_alpha,
      font_base_size = font_base_size
    )
    if(!return_plots) plot(p)
    dev.flush()  # Flush the plot
    flush.console()  # Ensure plot is rendered
    if (show_progress) setTxtProgressBar(pb, i)  # Update progress bar
    p_list[[i]] <- p
  }

  if (output_pdf) dev.off()  # Close PDF device
  message(" - done!")  # Completion message
  if (show_progress) close(pb)  # Close progress bar if open

  # Return plot list or invisible
  if (return_plots) {
    return(p_list[page_range])
  } else {
    invisible()
  }
}


plot_feature_correlations_page <- function(d_plot, ...){

  args <- base::list(...)
  # Subset dataset for current page
  n_samples <- length(unique(d_plot$analysis_id))
  row_start <- n_samples * args$cols_page * args$rows_page * (args$specific_page - 1) + 1
  row_end <- n_samples * args$cols_page * args$rows_page * args$specific_page

  if (args$sort_by_corr) {
    d_plot <- d_plot |>
      arrange(desc(.data$abs_cor),.data$analysis_id) |>
      dplyr::mutate(pair = factor(.data$pair, levels = unique(.data$pair)))
  } else {
    d_plot <- d_plot |>
      arrange(.data$pair, .data$analysis_id) |>
      dplyr::mutate(pair = factor(.data$pair, levels = unique(.data$pair)))
  }

  d_plot <- d_plot |>
    slice(row_start:row_end)

  # Create plot
  p <- d_plot |>
    ggplot2::ggplot(ggplot2::aes(x = .data$x, y = .data$y)) +
    ggplot2::geom_point(size = args$point_size,
                        aes(color = .data$qc_type, shape = .data$qc_type, fill = .data$qc_type),
                        alpha = args$point_alpha,
                        stroke = args$point_stroke) +
    ggplot2::scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
    ggplot2::scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
    ggplot2::scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)
  if (args$log_scale) {
    p <- p +
      # ggplot2::scale_x_log10(
      #                        expand = ggplot2::expansion(mult = c(0, 0.05))) +
      # ggplot2::scale_y_log10(
      #                        expand = ggplot2::expansion(mult = c(0, 0.05)))
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
    suppressWarnings({ggplot2::geom_smooth(method = "lm", formula = y ~ x, linewidth = args$line_size, color = args$line_color, alpha = args$line_alpha, se = FALSE,  na.rm = TRUE)}) +
    ggplot2::geom_text(
      data = d_plot |>
        dplyr::group_by(.data$pair, .data$r) |>
        dplyr::slice(1),
      ggplot2::aes(label = .data$r),
      x = -Inf, y = Inf,
      hjust = -0.1, vjust = 1.5,
      size = args$font_base_size / 3.5,
    ) +
    ggh4x::facet_wrap2(~pair, scales = "free", ncol = args$cols_page, nrow = args$rows_page, trim_blank = FALSE) +
    theme_bw(base_size = args$font_base_size) +
    theme(
      plot.title = element_text(size = args$font_base_size, face = "bold"),
      strip.text = ggplot2::element_text(size = args$font_base_size, face = "bold"),
      axis.text = element_text(size = args$font_base_size),
      axis.title = element_text(size = args$font_base_size, face = "bold"),
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
