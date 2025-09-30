#' RunScatter Plot
#'
#' @description
#' The `runscatter` function visualizes raw or processed feature signals across
#' different sample/QC types along the analysis sequence. It helps identify
#' trends, detect outliers, and assess analytical performance. Available
#' feature variables, such as retention time (RT) and full width at
#' half maximum (FWHM), can be plotted against analysis order or timestamps.
#'
#' By default, all QC types present in the dataset will be plotted. QC types that
#' predefined colors or shapes are assigned black shapes.
#' User-defined QC types that have no predefined colors or shapes in midar.
#' will be assigned black shapes.
#
#' have no predefined color and shape, will be assigned shapes in black. To show specific
#' QC types use the `qc_types` argument.
#'
#' To plot the feature values before the last applied drift/batch correction,
#' add `*_before` to the variable name, e.g., `intensity_before` or `conc_before`.
#' To plot the uncorrected feature values (before any drift/batch correction),
#' add `*_raw` to the variable name, e.g., `intensity_raw` or `conc_raw`. To show
#' corresponding fit curves, set `show_trend = TRUE`.
#'
#' The function also supports visualizing analysis batches, reference lines
#' (mean \eqn{\pm} SD), and trends. It offers customization options to display batch
#' separators, apply outlier capping, show smoothed trend curves, add reference
#' lines, and incorporate other features. Outlier capping is particularly useful
#' to focus on QC or study sample trends that might otherwise be obscured by
#' extreme values or high variability.
#'
#' The `runscatter` function serves as a central QC tool in the workflow,
#' providing critical insights into data quality.

#' @details
#' - The outlier capping feature (`cap_outliers`) allows you to cap upper outliers
#' based on median absolute deviation (MAD) fences of SPL and QC samples, or to
#' remove the top n points. This can help to focus on the trends of interest when
#' there are outlier or a high variability in the data, e.g. in the study samples.
#'
#' - When using log-scale (`log_scale = TRUE`), zero or negative values will
#' replaced with the minimum positive value divided by 5 to avoid log 0 errors
#'
#' - Reference lines/ranges corresponding to mean \eqn{\pm} k x SD can be shown across
#' or within batches as lines or shaded stripes.
#'
#' - Trend curves can be displayed before or after drift/batch correction. In
#' either case, a drift and/or batch correction must be applied to the data
#' to enable plotting of trend curves. To show trend curves used for the last drift or batch correction,
#' add "_before" to the variable name, e.g. `conc_before` or `intensity_before`
#' and set `show_trend = TRUE`.
#'
#'
#' @param data A `MidarExperiment` object containing the dataset and metadata.
#' @param variable The variable to plot on the y-axis, one of 'intensity',
#' 'norm_intensity', 'conc', 'conc', 'rt', 'fwhm', 'area', 'height', response'.
#' Add `_before` after the variable name to plot the feature values before the last
#' applied drift/batch correction, (e.g., `conc_before`). Add `_raw` after the
#' variable name to plot the raw uncorrected feature values (e.g., `conc_raw`).
#'
#' @param filter_data Logical, whether to use QC-filtered data based on criteria set via `filter_features_qc()`.
#' @param qc_types QC types to be plotted. Can be a vector of QC types or a regular expression pattern. `NA` (default) displays all available QC/Sample types.
#' @param include_qualifier Logical, whether to include qualifier features. Default is `TRUE`.
#' @param include_istd Logical, whether to include internal standard (ISTD) features. Default is `TRUE`.
#' @param include_feature_filter A regex pattern or a vector of feature names used to filter features by `feature_id`.
#' If `NA` or an empty string (`""`) is provided, the filter is ignored. When a vector of length > 1 is supplied,
#' is supplied, only features with exactly these names are selected (applied individually as OR conditions).
#' @param exclude_feature_filter A regex pattern or a vector of feature names to exclude features by feature_id.
#' If `NA` or an empty string (`""`) is provided, the filter is ignored. When a vector of length > 1 is supplied,
#' is supplied, only features with exactly these names are excluded (applied individually as OR conditions).
#' @param plot_range Numeric vector of length 2, specifying the start and end indices of the analysis order to be plotted. `NA` plots all samples.
#'
#' @param output_pdf Logical, whether to save the plot as a PDF file.
#' @param path File name for the PDF output.
#' @param multithreading Logical, whether to use parallel processing to speed up plot generation.
#' @param return_plots Logical, whether to return the list of ggplot objects.
#'
#' @param show_batches Logical, whether to show batch separators in the plot.
#' @param batch_zebra_stripe Logical, whether to display batches with alternating shaded and non-shaded areas.
#' @param batch_line_color Color of the batch separator lines.
#' @param batch_fill_color Color for the shaded areas representing batches.
#'
#' @param cap_outliers Logical, whether to cap upper outliers based on MAD fences of SPL and QC samples.
#' @param cap_sample_k_mad Numeric, k * MAD (median absolute deviation) for outlier capping of SPL samples.
#' @param cap_qc_k_mad Numeric, k * MAD (median absolute deviation) for outlier capping of QC samples.
#' @param cap_top_n_outliers Numeric, cap the top n outliers regardless of MAD fences. `NA` or `0` ignores this filter.
#'
#' @param show_reference_lines Whether to display reference lines (mean \eqn{\pm} n x SD).
#' @param reference_k_sd Multiplier for standard deviations to define SD reference lines.
#' @param ref_qc_types QC type for which the reference lines are calculated.
#' @param reference_batchwise Whether to calculate reference lines per batch.
#' @param reference_sd_shade `TRUE` plots a colored band indicating the \eqn{\pm} n x SD
#' reference range. `FALSE` (default) shows reference lines instead.
#' @param reference_line_color Color of the reference lines.
#' @param reference_fill_color Fill color of the batch-wise reference ranges.
#' If `NA` (default), the color assigned to the qc_type is used.
#' @param reference_linewidth Width of the reference lines.
#'
#' @param show_trend If `TRUE` trend curves before or after drift/batch correction are shown.
#' @param trend_color Color of the trend curve.
#'
#' @param point_size Size of the data points. Default is `1.5`.
#' @param point_transparency Alpha transparency of the data points.
#' @param point_border_width Width of the data point borders.
#' @param y_label_text Override the default y-axis label text.
#' @param pages_per_core Number of pages to be plotted by core when multithreading is enabled. Default is `6`. Changing this number may improve performance.
#' @param show_gridlines Whether to show major x and y gridlines.
#' @param y_lim Numeric vector of length 2, specifying the lower and upper y-axis limits. Default is `c(0,NA)``, which sets the lower limit to 0 and the upper limit automatically.
#' @param log_scale Logical, whether to use a log10 scale for the y-axis.
#'
#' @param rows_page Number of rows per page.
#' @param cols_page Number of columns per page.
#' @param specific_page Show/save a specific page number only. `NA` plots/saves all pages.
#' @param page_orientation Page orientation, "LANDSCAPE" or "PORTRAIT".
#' @param base_font_size Base font size for the plot.
#' @param use_dingbats Logical, whether to use Dingbats font in the PDF output for improved plotting speed. Default is `TRUE`. Set to `FALSE` if your PDF viewer does not show points correctly.
#' @param show_progress Logical, whether to show a progress bar. Default is `TRUE`.
#'
#' @return A list of ggplot2 plots, or `NULL` if `return
#' @export

plot_runscatter <- function(
  data = NULL,
  variable = c(
    "intensity",
    "norm_intensity",
    "conc",
    "rt",
    "area",
    "height",
    "fwhm",
    "width",
    "symmetry",
    "intensity_raw",
    "intensity_before",
    "norm_intensity_raw",
    "norm_intensity_before",
    "conc_raw",
    "conc_before"
  ),
  # Data and filtering arguments
  filter_data = FALSE,
  qc_types = NA,
  include_qualifier = TRUE,
  include_istd = TRUE,
  include_feature_filter = NA,
  exclude_feature_filter = NA,

  plot_range = NA,
  # Output settings
  output_pdf = FALSE,
  path = NA,
  multithreading = FALSE,
  return_plots = FALSE,

  # Display of batches
  show_batches = TRUE,
  batch_zebra_stripe = FALSE,
  batch_line_color = "#cdf7d9",
  batch_fill_color = "grey93",

  # Outlier capping
  cap_outliers = FALSE,
  cap_sample_k_mad = 4,
  cap_qc_k_mad = 4,
  cap_top_n_outliers = NA,

  # Control reference lines (mean and SD lines)
  show_reference_lines = FALSE,
  ref_qc_types = NA,
  reference_k_sd = 2,
  reference_batchwise = FALSE,
  reference_line_color = "#04bf9a",
  reference_sd_shade = FALSE,
  reference_fill_color = NA,
  reference_linewidth = 0.75,

  # Smoothed trend curve for selected QC sample type
  show_trend = FALSE,
  trend_color = "#22e06b",

  # Plot customization
  y_lim = c(0, NA),
  log_scale = FALSE,
  show_gridlines = FALSE,
  point_size = 1.5,
  point_transparency = 1,
  point_border_width = NA,
  base_font_size = 10,

  # Layout settings
  rows_page = 3,
  cols_page = 3,
  specific_page = NA,
  page_orientation = "LANDSCAPE",

  # Others
  y_label_text = NA,
  pages_per_core = 1,

  use_dingbats = TRUE,

  # Progress bar settings
  show_progress = TRUE
) {
  # Check the validity of input data
  check_data(data)

  if (nrow(data@dataset) < 1) {
    cli::cli_abort(col_red("No data available. Please import data and metadata first."))
  }

  # Handle saving output to PDF
  if (output_pdf && (!is.character(path) || path == "")) {
    cli::cli_abort(
      "No valid path defined. When using `output_pdf = TRUE`, please set the output path
                   for the PDF via `path = `."
    )
  }

  if (multithreading) {
    res <- check_installed(
      "mirai",
      reason = "to use multithreading for plot generation."
    )
    if (rlang::is_installed("mirai") && mirai::daemons()$connections == 0) {
      cli::cli_alert_warning(cli::col_yellow(
        "To use multithreading for plot generation, please set `mirai::daemon(number_of_cores)` where `number_of_cores` is the number of CPU cores to use."
      ))
    }
  }

  # Match the selected variable with predefined options
  variable <- str_remove(variable, "feature_")
  rlang::arg_match(
    variable,
    c(
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
  )
  variable <- stringr::str_c("feature_", variable)
  variable_sym = rlang::sym(variable)

  # Check arguments are valid
  check_var_in_dataset(data@dataset, variable)


  if (!all(is.numeric(y_lim) || is.na(y_lim))) {
    cli::cli_abort(cli::col_red("`y_lim` must have numeric values or `NA`s."))
  }

  if (show_reference_lines && is.na(ref_qc_types)) {
    cli::cli_abort(
      "Please define a QC to show reference lines, via the `ref_qc_types` argument or set `show_reference_lines = FALSE`."
    )
  }
  if (
    show_reference_lines &&
      (!any(ref_qc_types %in% unique(data@dataset$qc_type)))
  ) {
    cli::cli_abort("Selected `ref_qc_types` not present in the dataset.")
  }

  if (str_detect(variable, "_before|_raw")) {
    if (!any(data@var_drift_corrected) && !any(data@var_batch_corrected)) {
      cli::cli_abort(col_red(
        "Variables `_before` and `_raw` after only available after drift/batch corrections. Please set chose an other variable, or first apply drift/batch correction."
      ))
    }
  }

  if (show_trend | str_detect(variable, "_before|_raw")) {
    if (!any(data@var_drift_corrected) && !any(data@var_batch_corrected)) {
      cli::cli_abort(col_red(
        "Drift or batch correction is currently required to show trend lines. Please set `show_trend = FALSE`, or apply drift/batch correction first."
      ))
    }
  }

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

  # Cleanup the data

  # Replace infinite values and store the variable values
  d_filt <- d_filt |>
    mutate(value = ifelse(is.infinite(!!variable_sym), NA, !!variable_sym))

  # Set the y-axis label text
  y_label <- dplyr::if_else(
    cap_outliers,
    paste0(
      ifelse(
        is.na(y_label_text),
        stringr::str_remove(variable, "feature\\_"),
        y_label_text
      ),
      " (capped by MAD outlier filter) "
    ),
    stringr::str_remove(variable, "feature\\_")
  )

  # Reorder QC types and assign values
  d_filt$qc_type <- d_filt$qc_type |>
    factor() |>
    forcats::fct_expand(pkg.env$qc_type_annotation$qc_type_levels) |>
    forcats::fct_relevel(pkg.env$qc_type_annotation$qc_type_levels)

  d_filt <- d_filt |>
    dplyr::mutate(value = !!variable_sym)

  # Cap outliers if the option is selected
  if (cap_outliers) {
    outlier_offset_ratio <- 1.03

    if (
      is.na(cap_top_n_outliers) &&
        is.na(cap_sample_k_mad) &&
        is.na(cap_qc_k_mad)
    ) {
      cli::cli_abort(col_red(
        "One or more of `cap_sample_k_mad`, `cap_qc_k_mad`, and  `cap_top_n_outliers` must be a positive number when `cap_outlier = TRUE`, o."
      ))
    }

    if (is.na(cap_sample_k_mad)) {
      cap_sample_k_mad <- Inf
    }
    if (is.na(cap_qc_k_mad)) {
      cap_qc_k_mad <- Inf
    }

    if (!is.na(cap_top_n_outliers) && cap_top_n_outliers > 0) {
      d_filt <- d_filt |>
        dplyr::group_by(.data$feature_id) |>
        dplyr::arrange(desc(.data$value)) |>
        mutate(
          value = ifelse(
            dplyr::row_number() <= cap_top_n_outliers,
            .data$value[row_number() > cap_top_n_outliers] *
              outlier_offset_ratio,
            .data$value
          )
        ) |>
        dplyr::arrange(.data$feature_id, .data$analysis_order) |>
        dplyr::ungroup()
    }

    thresholds <- d_filt |>
      group_by(.data$feature_id) |>
      summarise(
        max_spl = median(.data$value[.data$qc_type == "SPL"], na.rm = TRUE) +
          cap_sample_k_mad *
            mad(.data$value[.data$qc_type == "SPL"], na.rm = TRUE),
        max_tqc = median(.data$value[.data$qc_type == "TQC"], na.rm = TRUE) +
          cap_qc_k_mad * mad(.data$value[.data$qc_type == "TQC"], na.rm = TRUE),
        max_bqc = median(.data$value[.data$qc_type == "BQC"], na.rm = TRUE) +
          cap_qc_k_mad * mad(.data$value[.data$qc_type == "BQC"], na.rm = TRUE),
        .groups = "drop"
      ) |>
      mutate(
        value_max = pmax(
          .data$max_spl,
          .data$max_tqc,
          .data$max_bqc,
          na.rm = TRUE
        )
      )

    d_filt <- d_filt |>
      left_join(thresholds, by = "feature_id") |>
      mutate(
        value_max = ifelse(
          is.infinite(.data$value_max),
          .data$value,
          .data$value_max
        ),
        value_mod = if_else(
          .data$value > .data$value_max,
          .data$value_max * outlier_offset_ratio,
          .data$value
        )
      )
  } else {
    d_filt <- d_filt |>
      dplyr::mutate(value_mod = .data$value)
  }
  if (log_scale) {
    # check if value_mod contains any negative or zero
    if (any(d_filt$value_mod <= 0)) {
      cli::cli_alert_warning(cli::col_yellow(
        "Zero or negative values were replaced with the minimum positive value divided by 5 to avoid log(0) errors."
      ))

      d_filt <- d_filt |>
        dplyr::mutate(
          value_mod = if_else(
            .data$value_mod <= 0,
            min(.data$value_mod[.data$value_mod > 0]) / 5,
            .data$value_mod
          )
        )
    }
  }

  # Determine page range for the plots
  if (!is.numeric(specific_page)) {
    page_range <- 1:ceiling(
      dplyr::n_distinct(d_filt$feature_id) /
        (cols_page * rows_page)
    )
  } else {
    page_range <- specific_page
  }

  arglist <- list(
    multithreading = multithreading,
    return_plots = return_plots,
    y_var = variable,
    d_batches = data@annot_batches,
    cols_page = cols_page,
    rows_page = rows_page,
    show_trend = show_trend,
    output_pdf = output_pdf,
    page_orientation = page_orientation,
    point_size = point_size,
    cap_outliers = cap_outliers,
    point_transparency = point_transparency,
    show_batches = show_batches,
    batch_zebra_stripe = batch_zebra_stripe,
    batch_line_color = batch_line_color,
    batch_fill_color = batch_fill_color,
    y_label = y_label,
    base_font_size = base_font_size,
    point_border_width = point_border_width,
    show_grid = show_gridlines,
    y_min = y_lim[1],
    y_max = y_lim[2],
    log_scale = log_scale,
    plot_range = plot_range,
    show_reference_lines = show_reference_lines,
    ref_qc_types = ref_qc_types,
    reference_k_sd = reference_k_sd,
    reference_batchwise = reference_batchwise,
    reference_line_color = reference_line_color,
    reference_sd_shade = reference_sd_shade,
    reference_fill_color = reference_fill_color,
    reference_linewidth = reference_linewidth,
    trend_color = trend_color,
    show_progress = show_progress,
    use_dingbats = use_dingbats
  )
  # subset the dataset with only the rows used for plotting the facets of the selected page
  n_samples <- length(unique(d_filt$analysis_id))

  # Arrange the data first
  d_arranged <- d_filt |>
    arrange(.data$feature_id, .data$analysis_order)

  # Calculate page size
  n_samples <- length(unique(d_arranged$analysis_id))
  page_size <- n_samples * cols_page * rows_page

  if (multithreading) {
    page_group_size <- page_size * pages_per_core
  } else {
    page_group_size <- ceiling(nrow(d_arranged))
    pages_per_core = ceiling(nrow(d_arranged) / page_size)
  }

  # Add page_id
  d_with_page <- d_arranged |>
    mutate(
      page_id = ceiling(row_number() / page_size),
      page_group = ceiling(.data$page_id / pages_per_core)
    )

  # Split into list of page groups for parallel processing if multithreading is enabled otherwise include all pages in one group
  page_group_list <- split(d_with_page, d_with_page$page_group)

  if (output_pdf) {
    path <- ifelse(
      stringr::str_detect(path, ".pdf"),
      path,
      paste0(path, ".pdf")
    )
    if (multithreading) {
      tmp_dir <- fs::dir_create(fs::path_temp("mrmhub_plotpages"))
      page_group_files <- paste0(
        tmp_dir,
        "runscatter_page_",
        seq_along(page_group_list),
        ".pdf"
      )
    } else {
      page_group_files <- path
    }
    action_text <- "Saving plots to pdf"
  } else {
    tmp_dir <- NULL
    page_group_files <- seq_along(page_group_list)
    action_text <- "Generating plots"
  }
  message(
    cli::col_green(glue::glue(
      "{action_text} ({max(page_range)} {ifelse(max(page_range) > 1, 'pages', 'page')}){ifelse(show_progress, '...', '')}"
    )),
    appendLF = FALSE
  )
  # if (show_progress) {
  #   pb <- txtProgressBar(min = 0, max = max(page_range), width = 30, style = 3)
  # }
  # flush.console()

  if (multithreading) {
    p_list <- purrr::map2(
      page_group_list,
      page_group_files,
      purrr::in_parallel(
        .f = function(pg_group, file) {
          args <- c(list(d_subset = pg_group, file = file), arglist)
          do.call(runscatter_plot_pages, args)
        },
        runscatter_plot_pages = runscatter_plot_pages,
        arglist = arglist
      ),
      .progress = show_progress
    )
  } else {
    p_list <- purrr::map2(
      page_group_list,
      page_group_files,
      .f = function(pg_group, file) {
        args <- c(list(d_subset = pg_group, file = file), arglist)
        do.call(runscatter_plot_pages, args)
      },
      .progress = show_progress
    )
  }

  p_list <- purrr::flatten(p_list)
  flush.console()


  if(output_pdf) message(cli::col_green("  done!"))

  flush.console()
  if (multithreading && output_pdf) {
    check_installed("qpdf")
    # Combine individual page PDFs into a single PDF
    qpdf::pdf_combine(page_group_files, output = path)
    fs::dir_delete(tmp_dir)
  }

  # if (show_progress) {
  #   close(pb)
  # }

  # flush.console()

  if (return_plots) {
    return(p_list[page_range])
  } else {
    invisible()
  }
}


runscatter_plot_pages <- function(
  d_subset,
  file,
  multithreading,
  return_plots,
  y_var,
  d_batches,
  cols_page,
  rows_page,
  show_trend,
  output_pdf,
  page_orientation,
  point_size,
  cap_outliers,
  point_transparency,
  show_batches,
  batch_zebra_stripe,
  batch_line_color,
  batch_fill_color,
  y_label,
  base_font_size,
  point_border_width,
  show_grid,
  y_min,
  y_max,
  log_scale,
  plot_range,
  show_reference_lines,
  ref_qc_types,
  reference_k_sd,
  reference_batchwise,
  reference_line_color,
  reference_sd_shade,
  reference_fill_color,
  reference_linewidth,
  trend_color,
  show_progress,
  use_dingbats
) {
  runscatter_one_page <- function(d_subset) {
    # For debugging
    # p <- ggplot(data = data.frame(speed = 1:4, dist = cumsum(runif(4, 0, 22))), aes(x = speed, y = dist)) + geom_point()
    # plot(p)
    # return(p)

    point_size <- ifelse(is.na(point_size), 2, point_size)

    if (is.na(point_border_width)) {
      point_border_width <- dplyr::if_else(output_pdf, .1, .4)
    }

    #d_subset$qc_type <- forcats::fct_relevel(d_subset$qc_type, pkg.env$qc_type_annotation$qc_type_levels)

    # Reorder QC types and assign values

    d_subset$qc_type <- d_subset$qc_type |>
      factor() |>
      forcats::fct_expand(pkg.env$qc_type_annotation$qc_type_levels) |>
      forcats::fct_relevel(rev(pkg.env$qc_type_annotation$qc_type_levels)) |>
      forcats::fct_rev()

    d_subset <- d_subset |> arrange(desc(.data$qc_type))

    defined_qctypes <- levels(d_subset$qc_type)

    # Handle QC types that are not defined in midar

    # Assign "grey50" to any undefined qc types
    qc_types_color <- setNames(
      rep("grey0", length(defined_qctypes)),
      defined_qctypes
    )
    qc_types_color[names(
      pkg.env$qc_type_annotation$qc_type_col
    )] <- pkg.env$qc_type_annotation$qc_type_col # Override known groups with defined colors
    qc_types_fill <- setNames(
      rep("grey40", length(defined_qctypes)),
      defined_qctypes
    )
    qc_types_fill[names(
      pkg.env$qc_type_annotation$qc_type_fillcol
    )] <- pkg.env$qc_type_annotation$qc_type_fillcol # Override known groups with defined colors

    undefined_qctypes <- setdiff(
      defined_qctypes,
      names(pkg.env$qc_type_annotation$qc_type_shape)
    )

    if (length(undefined_qctypes) > 0) {
      cli_alert_warning(col_yellow(
        "The QC types '{glue::glue_collapse(undefined_qctypes, sep = ', ')}' not predefined in Midar and will be displayed in black with auto-assigned shapes."
      ))
    }

    extra_shapes <- c(8, 3, 4, 9, 21, 22, 23, 24, 25, 7, 10, 12, 13, 14, 11) # Extra distinct shapes
    new_shapes <- setNames(
      extra_shapes[seq_along(undefined_qctypes)],
      undefined_qctypes
    )

    qc_types_shape <- c(pkg.env$qc_type_annotation$qc_type_shape, new_shapes)

    d_subset <- d_subset |>
      dplyr::arrange(rev(.data$qc_type))

    # https://stackoverflow.com/questions/46327431/facet-wrap-add-geom-hline
    if (nrow(d_subset) > 0) {
      dMax <- d_subset |>
        dplyr::group_by(.data$feature_id) |>
        dplyr::summarise(
          y_max = if (!all(is.na(.data$value_mod))) {
            max(.data$value_mod, na.rm = TRUE) * 1
          } else {
            NA_real_
          },
          y_min = if (!all(is.na(.data$value_mod))) {
            min(.data$value_mod, na.rm = TRUE) * 1
          } else {
            NA_real_
          }
        )
    }
    d_batches <- d_batches
    d_batch_data <- d_batches |>
      dplyr::slice(rep(1:dplyr::n(), each = nrow(dMax)))
    d_batch_data$feature_id <- rep(dMax$feature_id, times = nrow(d_batches))
    d_batch_data <- d_batch_data |> dplyr::left_join(dMax, by = c("feature_id"))

    p <- ggplot2::ggplot(d_subset, aes(x = !!sym("analysis_order")))

    # browser()
    if (show_batches) {
      if (!batch_zebra_stripe) {
        d_batches_temp <- d_batch_data |> filter(.data$id_batch_start != 1)
        p <- p +
          ggplot2::geom_vline(
            data = d_batches_temp,
            ggplot2::aes(xintercept = .data$id_batch_start - 0.5),
            colour = batch_line_color,
            linetype = "solid",
            linewidth = .5,
            na.rm = TRUE
          )
      } else {
        d_batches_temp <- d_batch_data |>
          dplyr::filter(.data$batch_no %% 2 != 1)
        p <- p +
          ggplot2::geom_rect(
            data = d_batches_temp,
            ggplot2::aes(
              xmin = .data$id_batch_start - 0.5,
              xmax = .data$id_batch_end + 0.5,
              ymin = .data$y_min,
              ymax = .data$y_max
            ),
            inherit.aes = FALSE,
            fill = batch_fill_color,
            color = NA,
            alpha = 1,
            linetype = "solid",
            linewidth = 0.3,
            na.rm = TRUE
          )
      }
    }

    # Plot reference range as shaded background
    if (show_reference_lines) {
      if (reference_batchwise) {
        grp = c("feature_id", "batch_id")
      } else {
        grp = c("feature_id")
      }
      d_subset_stats <- d_subset |>
        left_join(d_batches, by = c("batch_id")) |>
        filter(.data$qc_type == ref_qc_types) |>
        group_by(across(all_of(grp))) |>
        # TODO: could be cleaned up
        summarise(
          mean = mean(.data$value_mod, na.rm = TRUE),
          sd = if (!is.na(reference_k_sd)) {
            reference_k_sd * sd(.data$value_mod, na.rm = TRUE)
          } else {
            0
          },
          y_min = .data$mean - .data$sd,
          y_max = .data$mean + .data$sd,
          y_min_cap = if_else(.data$y_min < 0, 0, .data$y_min),
          y_max_cap = if_else(
            .data$y_max > safe_max(.data$value_mod, na.rm = TRUE),
            safe_max(.data$value_mod, na.rm = TRUE),
            .data$y_max
          ),
          batch_start = min(.data$id_batch_start),
          batch_end = max(.data$id_batch_end),
          batch_id = min(.data$batch_id),
          .groups = 'drop'
        )

      if (reference_sd_shade) {
        if (is.na(reference_fill_color)) {
          reference_fill_color <- pkg.env$qc_type_annotation$qc_type_col[
            ref_qc_types
          ]
        }
        p <- p +
          ggplot2::geom_rect(
            data = d_subset_stats,
            inherit.aes = FALSE,
            aes(
              xmin = .data$batch_start,
              xmax = .data$batch_end,
              ymin = .data$y_min_cap,
              ymax = .data$y_max_cap,
              group = .data$batch_id
            ),
            fill = reference_fill_color,
            linewidth = reference_linewidth,
            alpha = .15,
            na.rm = TRUE
          )
      }
    }

    if (cap_outliers) {
      p <- p +
        ggplot2::geom_hline(
          data = dMax,
          ggplot2::aes(yintercept = .data$y_max),
          color = "#ffefbf",
          linewidth = 3,
          alpha = 1,
          na.rm = TRUE
        )
    }

    p <- p +
      ggplot2::geom_point(
        aes(
          x = !!sym("analysis_order"),
          y = !!sym("value_mod"),
          color = .data$qc_type,
          fill = .data$qc_type,
          shape = .data$qc_type,
          group = .data$batch_id
        ),
        size = point_size,
        alpha = point_transparency,
        stroke = point_border_width,
        na.rm = TRUE
      )

    if (show_trend) {
      #browser()
      y_var_trend <- if_else(
        str_detect(y_var, "\\_before|\\_raw"),
        paste0(y_var, "_fit"),
        paste0(y_var, "_fit_after")
      )
      p <- p +
        ggplot2::geom_line(
          aes(
            x = !!sym("analysis_order"),
            y = !!sym(y_var_trend),
            group = .data$batch_id
          ),
          color = trend_color,
          linewidth = 1,
          na.rm = TRUE
        )
    }

    # Plot reference lines
    if (show_reference_lines) {
      if (reference_batchwise) {
        p <- p +
          ggplot2::geom_segment(
            data = d_subset_stats,
            inherit.aes = FALSE,
            aes(
              x = .data$batch_start,
              xend = .data$batch_end,
              y = .data$mean,
              yend = .data$mean,
              group = .data$batch_id
            ),
            color = reference_line_color,
            linewidth = reference_linewidth,
            linetype = "solid",
            alpha = 1
          )
      } else {
        p <- p +
          ggplot2::geom_hline(
            data = d_subset_stats,
            aes(yintercept = .data$mean),
            color = reference_line_color,
            linewidth = reference_linewidth,
            alpha = 1,
            linetype = "longdash",
            na.rm = TRUE
          )
      }

      if (!is.na(reference_k_sd) && !reference_sd_shade) {
        if (reference_batchwise) {
          p <- p +
            ggplot2::geom_segment(
              data = d_subset_stats,
              inherit.aes = FALSE,
              aes(
                x = .data$batch_start,
                xend = .data$batch_end,
                y = .data$y_min_cap,
                yend = .data$y_min_cap,
                group = .data$batch_id
              ),
              color = reference_line_color,
              linewidth = reference_linewidth,
              linetype = "dashed",
              alpha = 1
            ) +
            ggplot2::geom_segment(
              data = d_subset_stats,
              inherit.aes = FALSE,
              aes(
                x = .data$batch_start,
                xend = .data$batch_end,
                y = .data$y_max_cap,
                yend = .data$y_max_cap,
                group = .data$batch_id
              ),
              color = reference_line_color,
              linewidth = reference_linewidth,
              linetype = "dashed",
              alpha = 1
            )
        } else {
          p <- p +
            ggplot2::geom_hline(
              data = d_subset_stats,
              aes(yintercept = .data$y_min_cap),
              color = reference_line_color,
              linewidth = reference_linewidth,
              alpha = 1,
              linetype = "dashed",
              na.rm = TRUE
            ) +
            ggplot2::geom_hline(
              data = d_subset_stats,
              aes(yintercept = .data$y_max_cap),
              color = reference_line_color,
              linewidth = reference_linewidth,
              alpha = 1,
              linetype = "dashed",
              na.rm = TRUE
            )
        }
      }
    }

    p <- p +
      ggh4x::facet_wrap2(
        ggplot2::vars(.data$feature_id),
        scales = "free_y",
        ncol = cols_page,
        nrow = rows_page,
        trim_blank = FALSE
      ) +
      ggplot2::scale_color_manual(
        name = NULL,
        values = qc_types_color,
        drop = TRUE,
        na.value = "yellow"
      ) +
      ggplot2::scale_fill_manual(
        name = NULL,
        values = qc_types_fill,
        drop = TRUE,
        na.value = "yellow"
      ) +
      ggplot2::scale_shape_manual(
        name = NULL,
        values = qc_types_shape,
        drop = TRUE,
        na.value = 4
      )

    p <- p +
      # aes(ymin=0) +
      ggplot2::xlab("Analysis order") +
      ggplot2::ylab(label = y_label)

    if (log_scale) {
      p <- p + scale_y_log10(expand = ggplot2::expansion(mult = c(0.02, 0.03)))
    } else {
      p <- p +
        scale_y_continuous(
          limits = c(y_min, y_max),
          expand = ggplot2::expansion(mult = c(0.02, 0.03))
        ) +
        ggplot2::expand_limits(y = 0)
    }

    p <- p +
      ggplot2::theme_bw(base_size = base_font_size) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          size = base_font_size * 1,
          face = "bold"
        ),
        strip.text = ggplot2::element_text(
          size = base_font_size * 1,
          face = "bold"
        ),
        strip.background = ggplot2::element_rect(
          linewidth = 0.0001,
          fill = "#00283d"
        ),
        strip.text.x = ggplot2::element_text(color = "white"),
        axis.text.x = ggplot2::element_text(
          size = base_font_size * 0.8,
          ,
          face = NULL
        ),
        axis.text.y = ggplot2::element_text(
          size = base_font_size * 0.8,
          face = NULL
        ),
        axis.title = ggplot2::element_text(
          size = base_font_size,
          face = NULL
        ),
        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        strip.switch.pad.wrap = ggplot2::unit(-1, "mm"),
        panel.border = element_rect(linewidth = 0.5, color = "grey40")
      )

    if (show_grid) {
      p <- p +
        ggplot2::theme(
          panel.grid.major = ggplot2::element_line(
            linewidth = 0.3,
            colour = "grey88",
            linetype = "dashed"
          )
        )
    } else {
      p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
    }

    if (!all(is.na(plot_range))) {
      p <- p + ggplot2::coord_cartesian(xlim = plot_range)
    }

    plot(p)
    if (return_plots) return(p)
  }

  # Split into list of page groups for parallel processing if multithreading is enabled otherwise include all pages in one group
  #tick <- Sys.time()
  page_list <- split(d_subset, d_subset$page_id)
  #tock <- Sys.time()

  ##rint("Start PDF FILE")

  if (output_pdf) {
    pdf(
      file = file,
      onefile = !multithreading,
      #paper = "A4r",
      useDingbats = use_dingbats,
      useKerning = TRUE,
      width = ifelse(page_orientation == "LANDSCAPE", 28 / 2.54, 20 / 2.54),
      height = ifelse(page_orientation == "LANDSCAPE", 20 / 2.54, 28 / 2.54)
    )
  }

  # # add a ggplot test plot with penguins
  #p_list <- ggplot(data = data.frame(speed = 1:4, dist = cumsum(runif(4, 0, 22))), aes(x = speed, y = dist)) + geom_point()
  #plot(p_list)
  # return(p)

  #p_list <- runscatter_one_page(d_subset = d_subset)

  p_list <- purrr::map(
    page_list,
    #function(pg) {
    # Combine page-specific arguments with arglist
    #args <- c(list())
    ~ runscatter_one_page(d_subset = .x),
    .progress = show_progress
  )

  if (output_pdf) {
    dev.off()
  }

  return(p_list)
}
