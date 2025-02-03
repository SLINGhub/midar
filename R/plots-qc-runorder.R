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
#' @param show_batches Logical, whether to show batch separators in the plot
#' @param show_timestamp Logical, whether to use the acquisition timestamp as
#' the x-axis instead of the run sequence number
#' @param add_info_title Logical, whether to add a title with the experiment
#' title, analysis date, and analysis times
#' @param single_row Logical, whether to show all QC types in a single row
#' @param batch_zebra_stripe Logical, whether to show batches as shaded areas
#' instead of line separators
#' @param batch_line_color Color of the batch separator lines
#' @param batch_fill_color Color of the batch shaded areas
#' @param base_font_size Numeric, base font size for the plot
#' @return A ggplot object representing the run sequence plot
#' @export
plot_runsequence <- function(data = NULL,
                             qc_types = NA,
                             show_batches = TRUE,
                             show_timestamp = FALSE,
                             add_info_title = TRUE,
                             single_row = FALSE,
                             segment_width = 0.5,
                             batch_zebra_stripe = FALSE,
                             batch_line_color = "#b6f0c5",
                             batch_fill_color = "grey90",
                             base_font_size = 8) {

  # Check if data is valid
  check_data(data)

  # Extract the required columns from the dataset
  d_filt <- data$dataset |>
    select("run_seq_num", "acquisition_time_stamp", "batch_id", "analysis_id",
           "qc_type") |>
    distinct()

  # Filter QC types if provided
  if (!all(is.na(qc_types)) && length(qc_types) > 0) {
    d_filt <- d_filt |>
      filter(if (is.vector(qc_types) && length(qc_types) > 1) {
        .data$qc_type %in% qc_types
      } else {
        str_detect(.data$qc_type, qc_types)
      })
  }

  # Convert acquisition_time_stamp to POSIXct if using datetime
  if (show_timestamp) {
    d_filt$acquisition_time_stamp <- as.POSIXct(d_filt$acquisition_time_stamp)
  }

  # Convert qc_type to factor and create sample_category
  d_filt$qc_type <- factor(d_filt$qc_type,
                           levels = pkg.env$qc_type_annotation$qc_type_levels) |>
    forcats::fct_drop()
  d_filt$sample_category <- as.character(d_filt$qc_type)

  # Count samples per run_seq_num
  sample_counts <- d_filt |>
    group_by(.data$qc_type) |>
    summarise(sample_count = stringr::str_c("", dplyr::n()), .groups = 'drop')

  # Define QC colors
  qc_colors <- replace(pkg.env$qc_type_annotation$qc_type_col, "SPL", "grey35")

  # Initialize the ggplot object
  p <- ggplot(d_filt, aes(x = if (show_timestamp) .data$acquisition_time_stamp
                          else .data$run_seq_num,
                          y = rev(.data$qc_type),
                          color = .data$sample_category)) +
    labs(x = if (show_timestamp) "Acquisition Time" else "Analysis Order",
         y = "Sample Type") +
    theme_bw(base_size = base_font_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey80", linetype = "dotted",
                                        linewidth = 0.25),
      panel.grid.minor.x = element_line(colour = "grey90", linetype = "dotted",
                                        linewidth = 0.25),
      panel.border = element_rect(linewidth = 1),
      axis.title = element_text(face = "bold", size = base_font_size *0.8),
      axis.text.x = element_text(face = "plain", size = base_font_size),
      axis.text.y = element_text(face = "bold", size = base_font_size),
      axis.text.y.right = element_text(face = "plain", size = base_font_size),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_text(angle = 0, vjust = 0.5, size = base_font_size ),
      axis.ticks.x = element_line(colour = "grey80", linetype = "dotted", linewidth = 0.5),
      plot.margin = unit(c(1, 5, 1, 1), "lines"),
      legend.position = if (single_row) "right" else "none"  # Show legend if single_row
    )

  # Add batch shading if defined
  if (show_batches) {
    # Retrieve batch info
    d_batch_info <- data@annot_batches |>
      left_join(data$dataset |>
                  select("run_seq_num", "acquisition_time_stamp", "batch_id") |>
                  distinct(), by = c("id_batch_start" = "run_seq_num"),
                suffix = c("", "_start"), keep = FALSE) |>
      left_join(data$dataset |>
                  select("run_seq_num", "acquisition_time_stamp", "batch_id") |>
                  distinct(), by = c("id_batch_end" = "run_seq_num"),
                suffix = c("", "_end"), keep = FALSE)

    if (batch_zebra_stripe) {
      d_batch_shading <- d_batch_info |>
        slice(-1) |>
        filter(.data$batch_no %% 2 != 1)
      p <- p + geom_rect(data = d_batch_shading, inherit.aes = FALSE,
                         aes(xmin = if (show_timestamp) .data$acquisition_time_stamp
                             else .data$id_batch_start + 0.5,
                             xmax = if (show_timestamp) .data$acquisition_time_stamp_end
                             else .data$id_batch_end - 0.5,
                             ymin = -Inf, ymax = Inf),
                         fill = batch_fill_color, alpha = 1, color = NA)
    } else {
      p <- p + geom_vline(data = d_batch_info,
                          aes(xintercept = if (show_timestamp) .data$acquisition_time_stamp
                              else (.data$id_batch_end + 0.5)),
                          colour = batch_line_color, linewidth = segment_width * 2)
    }
  }

  # Add segments for qc_type
  if (single_row) {
    p <- p + geom_segment(aes(
      x = if (show_timestamp) .data$acquisition_time_stamp else .data$run_seq_num,
      xend = if (show_timestamp) .data$acquisition_time_stamp else .data$run_seq_num,
      y = -1, yend = 1
    ), linewidth = segment_width) +
      scale_y_continuous(breaks = NULL)  # Hide y-axis breaks
  } else {
    p <- p + geom_segment(aes(
      x = if (show_timestamp) .data$acquisition_time_stamp else .data$run_seq_num,
      xend = if (show_timestamp) .data$acquisition_time_stamp else .data$run_seq_num,
      y = as.integer(.data$qc_type) - 0.4,
      yend = as.integer(.data$qc_type) + 0.4
    ), linewidth = segment_width)

    # Position sample counts outside the y-axis
    p <- p +
      scale_y_continuous(breaks = seq(1, nlevels(d_filt$qc_type), by = 1),
                         labels = sample_counts$qc_type,
                         expand = expansion(0.02, 0.02),
                         sec.axis = ggplot2::sec_axis(~ ., name = "n",
                                                      breaks = seq(1, nlevels(d_filt$qc_type), by = 1),
                                                      labels = sample_counts$sample_count))
  }

  # Format x-axis as date-time if using acquisition_time_stamp
  if (show_timestamp) {
    p <- p + ggplot2::scale_x_datetime(date_labels = "%Y-%m-%d",
                                       expand = expansion(0.02, 0.02),
                                       date_breaks = "day", date_minor_breaks = "hour")
  } else {
    p <- p + scale_x_continuous(expand = expansion(0.02, 0.02),
                                breaks = seq(0, max(d_filt$run_seq_num),
                                             10^ceiling(log10(max(d_filt$run_seq_num))) / 10))
  }

  # Add additional information in the title
  if (add_info_title) {

    title_text <- if(data@title == "") "A" else glue::glue("{data@title} — A" )

    p <- p + labs(title = glue::glue("{title_text}nalysis time: {get_analysis_duration(data, estimate_sequence_end = TRUE) |> stringr::str_sub(end = -5)}  ({get_analyis_start(data) |> stringr::str_sub(end = -4)} - {get_analyis_end(data, estimate_sequence_end = TRUE) |> stringr::str_sub(end = -4)}) — median run time: {get_runtime_median(data)@minute}: {get_runtime_median(data)@.Data} min — interruptions > 1 hour: {get_analysis_breaks(data, 60)}"))
  }

  # Color mapping
  p <- p + scale_color_manual(values = qc_colors, name = "QC type")

  p
}


#' RunScatter Plot
#'
#' @description
#' The `runscatter` function visualizes raw or processed feature signals across
#' different sample/QC types along the analysis sequence. It helps identify
#' trends, detect outliers, and assess analytical performance. Available
#' feature variables, such as retention time (RT) and full width at
#' half maximum (FWHM), can be plotted against analysis order or timestamps.
#' \n
#' The function also supports visualizing analysis batches, reference lines
#' (mean ± SD), and trends. It offers customization options to display batch
#' separators, apply outlier capping, show smoothed trend curves, add reference
#' lines, and incorporate other features. Outlier capping is particularly useful
#' to focus on QC or study sample trends that might otherwise be obscured by
#' extreme values or high variability.
#' \n
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
#' - Reference lines/ranges corresponding to mean ± k x SD can be shown across
#' or within batches as lines or shaded stripes.
#'
#' - Trend curves can be displayed before or after drift/batch correction. In
#' either case, a drift and/or batch correction must be applied to the data
#' to enable plotting of trend curves.
#'
#'
#' @param data A `MidarExperiment` object containing the dataset and metadata.
#' @param variable The variable to plot on the y-axis, one of 'area', 'height', 'intensity', 'norm_intensity', 'intensity_raw', 'norm_intensity_raw', 'response', 'conc', 'conc_raw', 'rt', 'fwhm.'
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
#' @param analysis_order_range Numeric vector of length 2, specifying the start and end indices of the analysis order to be plotted. `NA` plots all samples.
#'
#' @param output_pdf Logical, whether to save the plot as a PDF file.
#' @param path File name for the PDF output.
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
#' @param show_reference_lines Whether to display reference lines (mean ± n x SD).
#' @param reference_k_sd Multiplier for standard deviations to define SD reference lines.
#' @param reference_qc_type QC type for which the reference lines are calculated.
#' @param reference_batchwise Whether to calculate reference lines per batch.
#' @param reference_sd_shade `TRUE` plots a colored band indicating the ± n x SD
#' reference range. `FALSE` (default) shows reference lines instead.
#' @param reference_line_color Color of the reference lines.
#' @param reference_fill_color Fill color of the batch-wise reference ranges.
#' If `NA` (default), the color assigned to the qc_type is used.
#' @param reference_linewidth Width of the reference lines.
#'
#' @param show_trend If `TRUE` trend curves before or after drift/batch correction are shown.
#' @param trend_color Color of the trend curve.
#'
#' @param point_size Size of the data points.
#' @param point_transparency Alpha transparency of the data points.
#' @param point_border_width Width of the data point borders.
#' @param y_label_text Override the default y-axis label text.
#' @param show_gridlines Whether to show major x and y gridlines.
#' @param log_scale Logical, whether to use a log10 scale for the y-axis.
#'
#' @param rows_page Number of rows per page.
#' @param cols_page Number of columns per page.
#' @param specific_page Show/save a specific page number only. `NA` plots/saves all pages.
#' @param page_orientation Page orientation, "LANDSCAPE" or "PORTRAIT".
#' @param base_font_size Base font size for the plot.

#' @param show_progress Logical, whether to show a progress bar.
#'
#' @return A list of ggplot2 plots, or `NULL` if `return

plot_runscatter <- function(data = NULL,
                            variable = c(
                              "intensity", "norm_intensity", "conc",
                              "conc_raw", "rt", "area", "height", "fwhm"
                            ),
                            # Data and filtering arguments
                            filter_data = FALSE,
                            qc_types = NA,
                            include_qualifier = TRUE,
                            include_istd = TRUE,
                            include_feature_filter = NA,
                            exclude_feature_filter = NA,

                            analysis_order_range = NA,

                            # Output settings
                            output_pdf = FALSE,
                            path = NA,
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
                            reference_qc_type = NA,
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
                            log_scale = FALSE,
                            show_gridlines = FALSE,
                            point_size = 1.5,
                            point_transparency = 1,
                            point_border_width = 1,
                            base_font_size = 11,

                            # Layout settings
                            rows_page = 3,
                            cols_page = 3,
                            specific_page = NA,
                            page_orientation = "LANDSCAPE",

                            # Others
                            y_label_text = NA,

                            # Progress bar settings
                            show_progress = FALSE) {
  # Check the validity of input data
  check_data(data)
  if (nrow(data@dataset) < 1)
    cli::cli_abort("No data available. Please import data and metadata first.")

  # Handle saving output to PDF
  if (output_pdf && (!is.character(path) || path == "")){
    cli::cli_abort("No valid path defined. When using `output_pdf = TRUE`, please set the output path
                   for the PDF via `path = `.")
  }

  # Match the selected variable with predefined options
  variable <- str_remove(variable, "feature_")
  rlang::arg_match(variable, c("area", "height", "intensity", "norm_intensity",
                               "intensity_raw", "norm_intensity_raw", "response",
                               "conc", "conc_raw", "rt", "fwhm"))
  variable <- stringr::str_c("feature_", variable)
  variable_sym = rlang::sym(variable)

  # Check arguments are valid
  check_var_in_dataset(data@dataset, variable)
  if (str_detect(variable, "\\_raw") && !any(data@var_drift_corrected) && !any(data@var_batch_corrected)){
      cli::cli_abort(cli::col_red("`{variable} is only available after drift or/and batch correction. Please run drift and/or batch corrections, or choose another variable."))
  }



  if(show_reference_lines && is.na(reference_qc_type)) {
    cli::cli_abort("Please define a QC to show reference lines, via the `reference_qc_type` argument or set `show_reference_lines = FALSE`.")
  }
  if(show_reference_lines && (!any(reference_qc_type %in% unique(data@dataset$qc_type)))) {
    cli::cli_abort("Selected `reference_qc_type` not present in the dataset.")
  }

  if (show_trend) {
    if (!any(data@var_drift_corrected) && !any(data@var_batch_corrected)) {
      cli::cli_abort(col_red("Drift or batch correction is currently required to show trend lines. Please apply drift/batch correction first."))
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
  y_label <- dplyr::if_else(cap_outliers,
                            paste0(ifelse(is.na(y_label_text),
                                          stringr::str_remove(variable, "feature\\_"),
                                          y_label_text),
                                   " (capped by MAD outlier filter) "),
                            stringr::str_remove(variable, "feature\\_"))


  # Reorder QC types and assign values
  d_filt$qc_type <- factor(as.character(d_filt$qc_type), pkg.env$qc_type_annotation$qc_type_levels)
  d_filt <- d_filt |>
    dplyr::mutate(value = !!variable_sym)

  # Cap outliers if the option is selected
  if (cap_outliers) {

    outlier_offset_ratio <- 1.03


    if(is.na(cap_top_n_outliers) && is.na(cap_sample_k_mad) && is.na(cap_qc_k_mad))
      cli::cli_abort(col_red("One or more of `cap_sample_k_mad`, `cap_qc_k_mad`, and  `cap_top_n_outliers` must be a positive number when `cap_outlier = TRUE`, o."))

    if (is.na(cap_sample_k_mad)) cap_sample_k_mad <- Inf
    if (is.na(cap_qc_k_mad)) cap_qc_k_mad <- Inf

    if (!is.na(cap_top_n_outliers) & cap_top_n_outliers > 0) {


      d_filt <- d_filt |>
        dplyr::group_by(.data$feature_id) |>
        dplyr::arrange(desc(.data$value)) |>
        mutate(value = ifelse(dplyr::row_number() <= cap_top_n_outliers,
                              .data$value[row_number() > cap_top_n_outliers] * outlier_offset_ratio, .data$value)) |>
        dplyr::arrange(.data$feature_id, .data$run_seq_num) |>
        dplyr::ungroup()
    }

    d_filt <- d_filt |>
      dplyr::group_by(.data$feature_id) |>
      dplyr::mutate(
        value_max_spl = median(.data$value[.data$qc_type == "SPL"], na.rm = TRUE) +
          cap_sample_k_mad * mad(.data$value[.data$qc_type == "SPL"],
                                 na.rm = TRUE),
        value_max_tqc = median(.data$value[.data$qc_type == "TQC"], na.rm = TRUE) +
          cap_qc_k_mad * mad(.data$value[.data$qc_type == "TQC"],
                             na.rm = TRUE),
        value_max_bqc = median(.data$value[.data$qc_type == "BQC"], na.rm = TRUE) +
          cap_qc_k_mad * mad(.data$value[.data$qc_type == "BQC"],
                             na.rm = TRUE),
        value_max = pmax(.data$value_max_spl, .data$value_max_tqc, .data$value_max_bqc,
                         na.rm = TRUE),
        value_max = ifelse(is.infinite(.data$value_max), .data$value, .data$value_max),
        value_mod = dplyr::if_else(
          .data$value > .data$value_max,
          if (!all(is.na(.data$value))) {
            max(.data$value[.data$value <= .data$value_max], na.rm = TRUE) * outlier_offset_ratio
          } else {
            NA_real_
          },
          .data$value
        )
      ) |>
      dplyr::ungroup()
  } else {
    d_filt <- d_filt |>
      dplyr::mutate(value_mod = .data$value)
  }

  if (log_scale){
    # check if value_mod contains any negative or zero
    if(any(d_filt$value_mod <= 0)){
      cli::cli_alert_warning(cli::col_yellow("Zero or negative values were replaced with the minimum positive value divided by 5 to avoid log(0) errors."))

      d_filt <- d_filt |>
        dplyr::mutate(value_mod = if_else(.data$value_mod <= 0, min(.data$value_mod[.data$value_mod > 0])/5, .data$value_mod))
    }
  }

  if (output_pdf) { # nocov start
    path <- ifelse(stringr::str_detect(path, ".pdf"), path, paste0(path, ".pdf"))
    if (page_orientation == "LANDSCAPE") {
      pdf(file = path, onefile = T, paper = "A4r", useDingbats = FALSE,
          width = 28 / 2.54, height = 20 / 2.54)
    } else {
      pdf(file = path, onefile = T, paper = "A4", useDingbats = FALSE,
          height = 28 / 2.54, width = 20 / 2.54)
    }
  } # nocov end

  # Determine page range for the plots
  if (!is.numeric(specific_page)) {
    page_range <- 1:ceiling(dplyr::n_distinct(d_filt$feature_id) /
                              (cols_page * rows_page))
  } else {
    page_range <- specific_page
  }

  # Prepare and render the plots for each page
  if(output_pdf)
    action_text = "Saving plots to pdf"
  else
    action_text = "Generating plots"

  message(cli::col_green(glue::glue("{action_text} ({max(page_range)} {ifelse(max(page_range) > 1, 'pages', 'page')}){ifelse(show_progress, ':', '...')}")))
  if(show_progress) pb <- txtProgressBar(min = 0, max = max(page_range),
                                         width = 50, style = 3)

  p_list <- list()
  for (i in page_range) {
    p <- runscatter_one_page(
      d_filt = d_filt, data = data, y_var = variable, d_batches = data@annot_batches,
      cols_page = cols_page, rows_page = rows_page, show_trend = show_trend,
      fit_qc_type = fit_qc_type, output_pdf = output_pdf, page_no = i,
      point_size = point_size, cap_outliers = cap_outliers, point_transparency =
        point_transparency,
      show_batches = show_batches, batch_zebra_stripe = batch_zebra_stripe,
      batch_line_color = batch_line_color, batch_fill_color =
        batch_fill_color, y_label = y_label, base_font_size = base_font_size,
      point_border_width = point_border_width, show_grid = show_gridlines,
      log_scale = log_scale, analysis_order_range = analysis_order_range,
      show_reference_lines = show_reference_lines, reference_qc_type = reference_qc_type, reference_k_sd = reference_k_sd,
      reference_batchwise = reference_batchwise, reference_line_color = reference_line_color, reference_sd_shade = reference_sd_shade, reference_fill_color = reference_fill_color,
      reference_linewidth = reference_linewidth, trend_color = trend_color
    )

    plot(p)
    dev.flush()
    flush.console()
    if(show_progress) setTxtProgressBar(pb, i)
    p_list[[i]] <- p
  }

  if (output_pdf) dev.off() # nocov
  message(cli::col_green(" - done!"))
  if(show_progress) close(pb)

  flush.console()

  if (return_plots) {
    return(p_list[page_range])
  }
}



runscatter_one_page <- function(d_filt, data, y_var, d_batches, cols_page, rows_page, page_no,
                                show_trend, fit_qc_type, cap_outliers,
                                show_batches, batch_zebra_stripe, batch_line_color, batch_fill_color,
                                output_pdf, point_transparency, point_size = point_size, y_label, base_font_size, point_border_width,
                                show_grid, log_scale, analysis_order_range, show_reference_lines,reference_qc_type, reference_k_sd, reference_batchwise, reference_line_color, reference_sd_shade, reference_fill_color,
                                reference_linewidth, trend_color) {
  point_size <- ifelse(missing(point_size), 2, point_size)
  point_border_width <- dplyr::if_else(output_pdf, .3, .5)


  # subset the dataset with only the rows used for plotting the facets of the selected page
  n_cmpd <- length(unique(d_filt$analysis_id))
  row_start <- n_cmpd * cols_page * rows_page * (page_no - 1) + 1
  row_end <- n_cmpd * cols_page * rows_page * page_no


  d_subset <- d_filt |>
    dplyr::arrange(.data$feature_id, .data$run_seq_num) |>
    dplyr::slice(row_start:row_end)

  d_subset$qc_type <- forcats::fct_relevel(d_subset$qc_type, pkg.env$qc_type_annotation$qc_type_levels)
  d_subset <- d_subset |>
    dplyr::arrange(rev(.data$qc_type))

  # https://stackoverflow.com/questions/46327431/facet-wrap-add-geom-hline
  if (nrow(d_subset) > 0) {
    dMax <- d_subset |>
      dplyr::group_by(.data$feature_id) |>
      dplyr::summarise(
        y_max =
          if (!all(is.na(.data$value_mod))) {
            max(.data$value_mod, na.rm = TRUE) * 1.05
          } else {
            NA_real_
          },
        y_min =
          if (!all(is.na(.data$value_mod))) {
            min(.data$value_mod, na.rm = TRUE) * 0.95
          } else {
            NA_real_
          }
      )
  }


  d_batch_data <- d_batches |> dplyr::slice(rep(1:dplyr::n(), each = nrow(dMax)))
  d_batch_data$feature_id <- rep(dMax$feature_id, times = nrow(d_batches))
  d_batch_data <- d_batch_data |> dplyr::left_join(dMax, by = c("feature_id"))

  p <- ggplot2::ggplot(d_subset, aes(x = !!sym("run_seq_num")))

  # browser()
  if (show_batches) {
    if (!batch_zebra_stripe) {
      d_batches_temp <- d_batch_data |> filter(.data$id_batch_start != 1)
      p <- p + ggplot2::geom_vline(data = d_batches_temp, ggplot2::aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, linetype = "solid", linewidth = .5, na.rm = TRUE)
    } else {
      d_batches_temp <- d_batch_data |> dplyr::filter(.data$batch_no %% 2 != 1)
      p <- p + ggplot2::geom_rect(
        data = d_batches_temp, ggplot2::aes(xmin = .data$id_batch_start - 0.5, xmax = .data$id_batch_end + 0.5, ymin = .data$y_min, ymax = .data$y_max),
        inherit.aes = FALSE, fill = batch_fill_color, color = NA, alpha = 1, linetype = "solid", linewidth = 0.3, na.rm = TRUE
      )
    }
  }

  # Plot reference range as shaded background
  if(show_reference_lines) {
    if(reference_batchwise) grp = c("feature_id", "batch_id") else grp = c("feature_id")
    d_subset_stats <- d_subset |>
      left_join(d_batches, by = c("batch_id")) |>
      filter(.data$qc_type == reference_qc_type) |>
      group_by(across(all_of(grp))) |>
      # TODO: could be cleaned up
      summarise(mean = mean(.data$value_mod, na.rm = TRUE),
                sd = if(!is.na(reference_k_sd)) reference_k_sd * sd(.data$value_mod, na.rm = TRUE) else 0,
                y_min = .data$mean - .data$sd,
                y_max = .data$mean + .data$sd,
                y_min_cap = if(y_min < 0) 0 else y_min,
                y_max_cap = if(y_max > max(.data$value_mod, na.rm = TRUE)) max(.data$value_mod, na.rm = TRUE) else .data$y_max,
                batch_start = min(.data$id_batch_start),
                batch_end = max(.data$id_batch_end),
                batch_id = min(.data$batch_id),
                .groups = 'drop')


    if(reference_sd_shade){
      if(is.na(reference_fill_color)) reference_fill_color <- pkg.env$qc_type_annotation$qc_type_col[reference_qc_type]
      p <- p +
        ggplot2::geom_rect(data = d_subset_stats, inherit.aes = FALSE, aes(xmin = .data$batch_start, xmax = .data$batch_end, ymin = .data$y_min_cap , ymax = .data$y_max_cap, group = .data$batch_id), fill = reference_fill_color, linewidth = reference_linewidth, alpha = .15)
    }
  }

  if (cap_outliers) {
    p <- p + ggplot2::geom_hline(data = dMax, ggplot2::aes(yintercept = .data$y_max), color = "#ffefbf", linewidth = 3, alpha = 1)
  }

  p <- p +
    ggplot2::geom_point(
      aes(
        x = !!sym("run_seq_num"),
        y = !!sym("value_mod"),
        color = qc_type,
        fill = qc_type,
        shape = qc_type,
        group = batch_id
      ),
      size = point_size,
      alpha = point_transparency,
      stroke = point_border_width,
      na.rm = TRUE)

  if (show_trend) {
    #browser()
    y_var_trend <- if_else(str_detect(y_var, "\\_raw"), paste0(y_var, "_fit"), paste0(y_var, "_fit_after"))
    p <- p +
      ggplot2::geom_line(aes(
        x = !!sym("run_seq_num"),
        y = !!sym(y_var_trend),
        group = batch_id
      ),
      color = trend_color,
      linewidth = 1,
      na.rm = TRUE)
  }

  # Plot reference lines
  if(show_reference_lines) {
    if(reference_batchwise){
      p <- p +
        ggplot2::geom_segment(data = d_subset_stats, inherit.aes = FALSE, aes(x = .data$batch_start, xend = .data$batch_end, y = .data$mean, yend = .data$mean, group = .data$batch_id), color = reference_line_color, linewidth = reference_linewidth, linetype = "solid", alpha = 1)
      } else{
      p <- p +
        ggplot2::geom_hline(data = d_subset_stats, aes(yintercept = .data$mean), color = reference_line_color, linewidth = reference_linewidth, alpha = 1, linetype = "longdash")
      }
    if(!is.na(reference_k_sd) && !reference_sd_shade) {
      if(reference_batchwise){
        p <- p +
          ggplot2::geom_segment(data = d_subset_stats, inherit.aes = FALSE, aes(x = .data$batch_start, xend = .data$batch_end, y = .data$y_min_cap, yend = .data$y_min_cap, group = .data$batch_id), color = reference_line_color, linewidth = reference_linewidth, linetype = "dashed", alpha = 1) +
          ggplot2::geom_segment(data = d_subset_stats, inherit.aes = FALSE, aes(x = .data$batch_start, xend = .data$batch_end, y = .data$y_max_cap, yend = .data$y_max_cap, group = .data$batch_id), color = reference_line_color, linewidth = reference_linewidth, linetype = "dashed", alpha = 1)
       } else {
        p <- p +
          ggplot2::geom_hline(data = d_subset_stats, aes(yintercept = .data$y_min_cap), color = reference_line_color, linewidth = reference_linewidth, alpha = 1, linetype = "dashed") +
          ggplot2::geom_hline(data = d_subset_stats, aes(yintercept = .data$y_max_cap), color = reference_line_color, linewidth = reference_linewidth, alpha = 1, linetype = "dashed")
       }
    }
  }

  p <- p +
    ggh4x::facet_wrap2(ggplot2::vars(.data$feature_id), scales = "free_y", ncol = cols_page, nrow = rows_page, trim_blank = FALSE) +
    ggplot2::scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
    ggplot2::scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
    ggplot2::scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)



  p <- p +
    # aes(ymin=0) +
    ggplot2::xlab("Analysis order") +
    ggplot2::ylab(label = y_label)

  if (log_scale) {
    p <- p + scale_y_log10(expand = ggplot2::expansion(mult = c(0.02, 0.03)))

  } else{
    p <- p + scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.02, 0.03))) +
      ggplot2::expand_limits(y = 0)

  }

  p <- p +
    ggplot2::theme_bw(base_size = base_font_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = base_font_size * 1, face = "bold"),
      strip.text = ggplot2::element_text(size = base_font_size * 1, face = "bold"),
      strip.background = ggplot2::element_rect(linewidth = 0.0001, fill = "#00283d"),
      strip.text.x = ggplot2::element_text(color = "white"),
      axis.text.x = ggplot2::element_text(size = base_font_size * 0.8, , face = NULL),
      axis.text.y = ggplot2::element_text(size = base_font_size * 0.8 , face = NULL),
      axis.title = ggplot2::element_text(size = base_font_size, face = NULL),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.switch.pad.wrap = ggplot2::unit(-1, "mm"),
      panel.border = element_rect(linewidth = 0.5, color = "grey40")
    )

  if (show_grid) {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_line(linewidth = 0.3, colour = "grey88", linetype = "dashed"))
  } else {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  }




  if(!all(is.na(analysis_order_range))) {
    p <- p + ggplot2::coord_cartesian(xlim = analysis_order_range)
  }

  return(p)
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
#' @param variable Variable to plot, must be one of "intensity", "norm_intensity", "conc", "conc_raw", "area", "height", "fwhm".

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
#' @param analysis_order_range Numeric vector of length 2, specifying the start and end indices of the analysis order to be plotted. `NA` plots all samples.
#' @param show_timestamp Logical, whether to use the acquisition timestamp as
#' the x-axis instead of the run sequence number
#'
#' @param min_feature_intensity Numeric, exclude features with overall median signal below this value
#' @param y_lim Numeric vector of length 2, specifying the lower and upper y-axis limits. Default is `NA`, which uses limits calculated based on `ignore_outliers`.
#' @param ignore_outliers Logical, whether to exclude outlier values based on 4x MAD (median absolute deviation) fences

#' @param show_batches Logical, whether to show batch separators in the plot
#' @param batch_zebra_stripe Logical, whether to show batches as shaded areas instead of line separators
#' @param batch_line_color Character, color of the batch separator lines
#' @param batch_fill_color Character, color of the batch shaded areas

#' @param x_gridlines Logical, whether to show major x-axis gridlines
#' @param linewidth Numeric, line width used for whiskers of the boxplot
#' @param base_font_size Numeric, base font size for the plot
#' @param relative_log_abundances Logical, whether to use relative log abundances (RLA) or just log-transformed values
#' @return A ggplot object representing the RLA plot
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
                                variable = c("intensity", "norm_intensity", "conc", "conc_raw", "area", "height", "fwhm"),
                                filter_data = FALSE,
                                qc_types = NA,

                                include_qualifier = TRUE,
                                include_istd = TRUE,
                                include_feature_filter = NA,
                                exclude_feature_filter = NA,

                                analysis_order_range = NA,
                                show_timestamp = FALSE,

                                min_feature_intensity = 0,
                                y_lim = NA,
                                ignore_outliers = FALSE,

                                show_batches = TRUE,
                                batch_zebra_stripe = FALSE,
                                batch_line_color = "#b6f0c5",
                                batch_fill_color = "grey93",

                                x_gridlines = FALSE,
                                linewidth = 0.2,
                                base_font_size = 8,
                                relative_log_abundances = TRUE) {
  check_data(data)
  if (nrow(data@dataset) < 1) cli::cli_abort("No data available. Please import data and metadata first.")

  rlang::arg_match(rla_type_batch, c("within", "across"))


  # Check if selected variable is valid
  rlang::arg_match(variable, c("area", "height", "intensity", "response", "conc", "conc_raw", "rt", "fwhm"))
  variable <- str_remove(variable, "feature_")
  variable <- stringr::str_c("feature_", variable)

  # Check arguments are valid
  check_var_in_dataset(data@dataset, variable)
  if (str_detect(variable, "\\_raw") && !any(data@var_drift_corrected) && !any(data@var_batch_corrected)){
    cli::cli_abort(cli::col_red("`{variable} is only available after drift or/and batch correction. Please run drift and/or batch corrections, or choose another variable."))
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
  if(show_timestamp) {
    d_filt$acquisition_time_stamp <- as.POSIXct(d_filt$acquisition_time_stamp)
    x_axis_variable <- "acquisition_time_stamp"
  } else {
    x_axis_variable <- "run_seq_num"
  }

  x_axis_variable_sym <- rlang::sym(x_axis_variable)


  if (!all(is.na(analysis_order_range))) {
    d_filt <- d_filt |> dplyr::filter(.data$run_seq_num >= analysis_order_range[1] & .data$run_seq_num <= analysis_order_range[2]) |> droplevels()
  }



  d_filt <- d_filt |>
    mutate(value = ifelse(is.infinite(!!variable_sym), NA, !!variable_sym))

  d_filt <- d_filt |>
    dplyr::select(any_of(c("analysis_id", "run_seq_num", "acquisition_time_stamp", "qc_type", "batch_id", "feature_id", "feature_intensity", "feature_norm_intensity", "feature_conc"))) |>
    group_by(.data$feature_id) |>
    filter(median(.data$feature_intensity) >= min_feature_intensity) |>
    droplevels() |>
    dplyr::arrange(.data$run_seq_num)


  if (relative_log_abundances) {

    d_filt = d_filt |>
      mutate(val = log2(!!variable_sym))

    if(rla_type_batch == "within") grp = c("feature_id","batch_id") else grp = c("feature_id")

    d_filt_medians <- d_filt |>
      group_by(across(all_of(grp))) |>
      summarise(val_median = median(.data$val, na.rm = TRUE)) |> ungroup()

    d_filt <- d_filt |>
      left_join(d_filt_medians, by = grp) |>
      mutate(val_res = .data$val - .data$val_median)

  } else {
    d_filt <- d_filt |> mutate(val_res = log2(!!variable_sym))
  }


  # Get labels corresponding to the breaks. TODO: write it more elegant and clear
  if(x_axis_variable != "run_seq_num"){
    labels <- unique(d_filt[[x_axis_variable]])[seq(1, length(unique(d_filt[[x_axis_variable]])), length.out = 20)]
    breaks <- d_filt |> filter(!!x_axis_variable_sym %in% labels) |> pull(.data$run_seq_num) |> unique()
  } else {
    breaks <- scales::breaks_pretty(n = 20)(range(d_filt$run_seq_num))
    labels = breaks
  }
  p <- ggplot(d_filt, aes(x = .data$run_seq_num, y = .data$val_res, group = .data$run_seq_num))

  p <- p + scale_x_continuous(breaks = breaks, labels = labels)

  if (show_batches) {
    if (!batch_zebra_stripe) {
      p <- p + geom_vline(data = data@annot_batches |> slice(-1), aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, linetype = "solid", linewidth = 1)
    } else {

      d_batch_2nd <- data@annot_batches |>
        slice(-1) |>
        filter(.data$batch_no %% 2 != 1)
      p <- p + geom_rect(
        data = d_batch_2nd, aes(xmin = .data$id_batch_start - 0.5, xmax = .data$id_batch_end + 0.5, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = batch_fill_color, color = NA, alpha = 1, linetype = "solid", linewidth = 0.5, na.rm = TRUE
      )
    }
  }


  x_text_angle <- ifelse(x_axis_variable != "run_seq_num", 90, 0)
  x_text_just <- ifelse(x_axis_variable != "run_seq_num", 1, 0.5)

  p <- p +
    geom_boxplot(aes(fill = .data$qc_type, color = .data$qc_type), notch = FALSE, outlier.colour = NA, linewidth = linewidth, na.rm = TRUE) +
    scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_col) +
    scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col) +
    theme_bw(base_size = base_font_size) +
    ylab(bquote(bold(log[2] ~ .(variable)))) +
    xlab("Analysis order") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(), #element_line(colour = "#bdbdbd", linetype = "dotted", size = .5),
      panel.grid.minor.x = element_blank(), #element_line(colour = "grey88", linetype = "dotted", size = .5),
      axis.text.y = element_text(size = base_font_size),
      axis.text.x = element_text(size = base_font_size, angle = x_text_angle, vjust = 0.5, hjust = x_text_just),
      axis.title = element_text(size = base_font_size * 1, face = "bold"),
      panel.border = element_rect(linewidth = 1, color = "grey20")
    )

  if(x_gridlines)
    p <- p + theme(panel.grid.major.x = element_line(colour = "#bdbdbd", linetype = "dotted", linewidth = .3))
  else
    p <- p + theme(panel.grid.major.x = element_blank())

  if (relative_log_abundances) {
    p <- p + geom_hline(yintercept = 0, colour = "#5fe3f5", linetype = "longdash", linewidth = 0.5) +
      ylab(bquote(bold(log[2] ~ "( relative" ~ .(stringr::str_remove(variable, "feature\\_")) ~ ")")))
  }
  ylim = c(NA, NA)
  xlim = c(NA, NA)
  # Set y-axis limits
  if (!all(is.na(y_lim))) {
    ylim = y_lim
  } else if(ignore_outliers) {
    tails <- get_mad_tails(d_filt$val_res, k = 4, na.rm = TRUE)
    ylim = tails
  }

  if(!all(is.na(analysis_order_range))) {
    xlim = analysis_order_range
  }

  p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  return(p)
}
