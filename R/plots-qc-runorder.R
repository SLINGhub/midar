#' Run Sequence Plot to Visualize Analysis Design and Timelines
#'
#' @description
#' The run sequence plot provides an overview of the analysis design and timelines, which can be useful for subsequent processing steps. The plot illustrates the batch structure, the quality control (QC) samples included with their respective positions, and additional information regarding the date, duration, and run time of the analysis.
#'
#' Setting `show_timestamp = TRUE` allows you to check for any interruptions in the analysis timeline.
#'
#' @param data MidarExperiment object
#' @param qc_types QC types to be shown. Can be a vector of QC types or a regular expression pattern. `NA` (default) displays all available QC/Sample types.
#' @param show_batches Logical, whether to show batch separators in the plot
#' @param show_timestamp Logical, whether to use the acquisition timestamp as the x-axis instead of the run sequence number
#' @param add_info_title Logical, whether to add a title with the experiment title, analysis date, and analysis times
#' @param single_row Logical, whether to show all QC types in a single row
#' @param batches_as_shades Logical, whether to show batches as shaded areas instead of line separators
#' @param batch_line_color Color of the batch separator lines
#' @param batch_shading_color Color of the batch shaded areas
#' @param scale_factor Numeric, overall plot scale factor
# #' @param show_excluded Logical, whether to show excluded analyses in a separate row
#' @param segment_width Numeric, line width of the segments
#' @param base_font_size Numeric, base font size for the plot
#' @return A ggplot object representing the run sequence plot
#' @export
qc_plot_runsequence <- function(data = NULL,
                                qc_types = NA,
                                show_batches = TRUE,
                                show_timestamp = FALSE,
                                add_info_title = TRUE,
                                single_row = FALSE,
                                segment_width = 0.5,
                                batches_as_shades = FALSE,
                                batch_line_color = "#b6f0c5",
                                batch_shading_color = "grey85",
                                scale_factor = 1,
                                #show_excluded = TRUE,
                                base_font_size = 8) {


  check_data(data)

  # Extract columns needed for plotting
  d_filt <- data$dataset |>
    select("run_seq_num", "acquisition_time_stamp", "batch_id", "analysis_id", "qc_type") |>
    distinct()

  # Filter qc_types if provided
  if (all(!is.na(qc_types)) & all(qc_types != "")) {
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
  d_filt$qc_type <- factor(d_filt$qc_type, levels = pkg.env$qc_type_annotation$qc_type_levels) |> forcats::fct_drop()
  d_filt$sample_category <- as.character(d_filt$qc_type)

  # Retrieve batch info
  d_batch_info <- data@annot_batches |>
    left_join(data$dataset |>
                select("run_seq_num", "acquisition_time_stamp", "batch_id") |>
                distinct(), by = c("id_batch_start" = "run_seq_num"), suffix = c("", "_start"), keep = FALSE) |>
    left_join(data$dataset |>
                select("run_seq_num", "acquisition_time_stamp", "batch_id") |>
                distinct(), by = c("id_batch_end" = "run_seq_num"), suffix = c("", "_end"), keep = FALSE)

  # Count samples per run_seq_num
  sample_counts <- d_filt |>
    group_by(.data$qc_type) |>
    summarise(sample_count = stringr::str_c("", dplyr::n()), .groups = 'drop')


  qc_colors <- replace(pkg.env$qc_type_annotation$qc_type_col, "SPL", "grey35")

  # Initialize ggplot
  p <- ggplot(d_filt, aes(x = if (show_timestamp) .data$acquisition_time_stamp else .data$run_seq_num,
                          y = rev(.data$qc_type),
                          color = .data$sample_category)) +
    labs(x = if (show_timestamp) "Acquisition Time" else "Analysis Order", y = "Sample Type") +
    theme_bw(base_size = base_font_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey80", linetype = "dotted", size = 0.25),
      panel.grid.minor.x = element_line(colour = "grey90", linetype = "dotted", size = 0.25),
      panel.border = element_rect(size = 1),
      axis.title = element_text(face = "bold", size = base_font_size - 1),
      axis.text.x = element_text(face = "plain", size = base_font_size),
      axis.text.y = element_text(face = "bold", size = base_font_size),
      axis.text.y.right = element_text(face = "plain", size = base_font_size),
      axis.ticks.y.right = element_blank(),
      axis.title.y.right = element_text(angle = 0, vjust = 0.5, size = 6 * scale_factor),
      axis.ticks.x = element_line(colour = "grey80", linetype = "dotted", size = 0.25),
      plot.margin = unit(c(1, 5, 1, 1), "lines"),
      legend.position = if (single_row) "right" else "none"  # Show legend if single_row
    )

  # Add batch shading if defined
  if (show_batches) {
    if (batches_as_shades) {
      d_batch_shading <- d_batch_info |>
        slice(-1) |>
        filter(.data$batch_no %% 2 != 1)
      p <- p + geom_rect(data = d_batch_shading,inherit.aes = FALSE,
                         aes(xmin = if(show_timestamp) .data$acquisition_time_stamp else .data$id_batch_start + 0.5,
                             xmax = if(show_timestamp) .data$acquisition_time_stamp_end else .data$id_batch_end - 0.5,
                             ymin = -Inf, ymax = Inf),
                         fill = batch_shading_color, alpha = 1, color = NA)
    } else {

      p <- p + geom_vline(data = d_batch_info |>  slice(-1),
                          aes(xintercept = if(show_timestamp) .data$acquisition_time_stamp else (.data$id_batch_end - .data$id_batch_start)/2),
                          colour = batch_line_color, size = 0.5)
    }
  }

  # Add segments for qc_type
  if (single_row) {
    p <- p + geom_segment(aes(
      x = if (show_timestamp) .data$acquisition_time_stamp else .data$run_seq_num,
      xend = if (show_timestamp) .data$acquisition_time_stamp else .data$run_seq_num,
      y = -1, yend = 1
    ), size = segment_width) +
      scale_y_continuous(breaks = NULL)  # Hide y-axis breaks
  } else {
    p <- p + geom_segment(aes(
      x = if (show_timestamp) .data$acquisition_time_stamp else .data$run_seq_num,
      xend = if (show_timestamp) .data$acquisition_time_stamp else .data$run_seq_num,
      y = as.integer(.data$qc_type) - 0.4,
      yend = as.integer(.data$qc_type) + 0.4
    ), size = segment_width)

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
    p <- p + ggplot2::scale_x_datetime(date_labels = "%Y-%m-%d", expand = expansion(0.02, 0.02), date_breaks = "day", date_minor_breaks = "hour")
  } else {
    p <- p + scale_x_continuous(expand = expansion(0.02, 0.02), breaks = seq(0, max(d_filt$run_seq_num), 10^ceiling(log10(max(d_filt$run_seq_num))) / 10))
  }

 if(add_info_title) {
    p <- p + labs(title = glue::glue("{data@title} \u2014 analysis time: {get_analysis_duration(data)|> stringr::str_sub(end = -5)}  ({get_analyis_start(data)|> stringr::str_sub(end = -4)} - {get_analyis_end(data)|> stringr::str_sub(end = -4)}) \u2014 median run time: {get_run_time(data)@minute}:{get_run_time(data)@.Data} min \u2014 interruptions > 1 hour: {get_analysis_interruptions(data, 3600)}"))
  }

  # Color mapping
  p <- p + scale_color_manual(values = qc_colors, name = "QC type")

  p
}

#' Create a Run Scatter Plot to Visualize Analyte Distributions
#'
#' @description
#' The run scatter plot allows you to visualize the sequence of measurements and distribution of QC (sample) types across the analysis run sequence.
#' This can be useful for identifying trends, outliers, or other patterns in the data.
#'
#' The function provides a range of options for customizing the plot, such as:
#' - Filtering the data (e.g., by QC type, feature ID)
#' - Applying outlier capping
#' - Displaying batch information and control limits
#' - Showing trend lines for drift correction
#' - Controlling the layout and appearance of the plot
#'
#' @param data MidarExperiment object
#' @param variable Variable to plot, one of 'area', 'height', 'intensity', 'norm_intensity', 'intensity_raw', 'norm_intensity_raw', 'response', 'conc', 'conc_raw', 'rt', 'fwhm.'
#' @param filter_data Logical, whether to use QC-filtered data based on criteria set via `qc_apply_feature_filter()`
#' @param qc_types QC types to be shown. Can be a vector of QC types or a regular expression pattern. `NA` (default) displays all available QC/Sample types.
#' @param include_qualifier Logical, whether to include qualifier features
#' @param include_feature_filter Regex pattern to select features by feature_id. `NA` or `""` ignores the filter.
#' @param exclude_feature_filter Regex pattern to exclude features by feature_id. `NA` or `""` ignores the filter.
#' @param plot_range_indices Numeric vector of length 2, specifying the start and end indices of the sequence to be plotted. `NA` plots all samples.
#' @param output_pdf Logical, whether to save the plot as a PDF file
#' @param path File name for the PDF output
#' @param log_scale Logical, whether to use a log10 scale for the y-axis
#' @param cap_outliers Logical, whether to cap upper outliers based on MAD fences of SPL and QC samples
#' @param cap_spl_k_mad Numeric, k * MAD (median absolute deviation) for outlier capping of SPL samples
#' @param cap_qc_k_mad Numeric, k * MAD (median absolute deviation) for outlier capping of QC samples
#' @param cap_top_n_values Numeric, cap the top n values regardless of MAD fences. `NA` or `0` ignores this filter.
#' @param show_batches Logical, whether to show batch separators in the plot
#' @param show_control_limits Logical, whether to show the mean line for the specified QC types. `NA` ignores the mean line.
#' @param set_limits_n_sd Numeric, show +/- n x SD lines for the QC types defined in `show_control_limits`. `NA` ignores the SD lines.
#' @param limits_batchwise Logical, whether to calculate the mean and SD lines per batch
#' @param limits_linecolor Color of the mean and SD lines
#' @param limits_linewidth Width of the mean and SD lines
#' @param batches_as_shades Logical, whether to show batches as shaded areas instead of line separators
#' @param batch_line_color Color of the batch separator lines
#' @param batch_shading_color Color of the batch shaded areas
#' @param show_trend Logical, whether to show a drift correction trend line
#' @param trend_linecolor Color of the trend line
#' @param fit_qc_type QC type used for smoothing the trend line
#' @param cols_page Number of columns per page
#' @param rows_page Number of rows per page
#' @param base_font_size Base font size for the plot
#' @param annot_scale Scale factor for text elements
#' @param page_orientation Page orientation, "LANDSCAPE" or "PORTRAIT"
#' @param point_transparency Alpha transparency of the data points
#' @param point_size Size of the data points
#' @param point_border_width Width of the data point borders
#' @param specific_page Integer, show/save a specific page number. `NA` plots/saves all pages.
#' @param y_label_text Override the default y-axis label text
#' @param return_plots Logical, whether to return the list of ggplot objects
#' @param show_gridlines Logical, whether to show major x and y gridlines
#' @param show_progress Logical, whether to show a progress bar
#' @return A list of ggplot2 plots, or `NULL` if `return_plots = FALSE`
#' @export

qc_plot_runscatter <- function(data = NULL,
                            variable = c("intensity", "norm_intensity", "conc", "conc_raw", "area", "height", "fwhm"),
                            filter_data = FALSE,
                            qc_types = NA,
                            include_qualifier = TRUE,
                            include_feature_filter = NA,
                            exclude_feature_filter = NA,
                            plot_range_indices = NA,
                            output_pdf = FALSE,
                            path = "",
                            log_scale = FALSE,
                            cap_outliers = FALSE,
                            cap_qc_k_mad = 4,
                            cap_spl_k_mad = 4,
                            cap_top_n_values = NA,
                            show_batches = TRUE,
                            show_control_limits = NA,
                            set_limits_n_sd = NA,
                            limits_batchwise = FALSE,
                            limits_linecolor = "#38dff5",
                            limits_linewidth = 0.75,
                            batches_as_shades = FALSE,
                            batch_line_color = "#b6f0c5",
                            batch_shading_color = "grey93",
                            cols_page = 3,
                            rows_page = 3,
                            base_font_size = 11,
                            annot_scale = 1.0,
                            show_trend = FALSE,
                            trend_linecolor = "#22e06b",
                            fit_qc_type = "BQC",
                            #show_trend_samples = FALSE,
                            #trend_samples_function = "loess",
                            #trend_samples_col = "darkred",
                            #show_before_after_smooth = FALSE,
                            #plot_other_qc = FALSE,
                            page_orientation = "LANDSCAPE",
                            point_transparency = 1,
                            point_size = 2,
                            point_border_width = 1,
                            specific_page = NA,
                            y_label_text = NA,
                            return_plots = FALSE,
                            show_gridlines = FALSE,
                            show_progress = FALSE) {

  check_data(data)
  if (nrow(data@dataset) < 1) cli::cli_abort("No data available. Please import data and metadata first.")


  # Check if selected variable is valid
  variable <- str_remove(variable, "feature_")
  rlang::arg_match(variable, c("area", "height", "intensity", "norm_intensity", "intensity_raw", "norm_intensity_raw", "response", "conc", "conc_raw", "rt", "fwhm"))
  variable <- stringr::str_c("feature_", variable)
  check_var_in_dataset(data@dataset, variable)
  variable_sym = rlang::sym(variable)

  # if (show_trend) {
  #   if(!data@var_drift_corrected[[str_remove(variable, "_raw")]])
  #     cli::cli_abort("This option is only available for drift corrected datasets. Please apply `correct_drift_...()` first.")
  # }

  # Filter data if filter_data is TRUE
  if (filter_data) {
    d_filt <- data@dataset_filtered |> dplyr::ungroup()
    if (!data@is_filtered) cli::cli_abort("Data has not been qc filtered, or has changed. Please run `qc_apply_feature_filter` first.")
  } else {
    d_filt <- data@dataset |> dplyr::ungroup()
  }

  if(!include_qualifier){
    d_filt <- d_filt |> filter(!.data$is_qualifier)
  }

  if (!all(is.na(plot_range_indices))) {
    d_filt <- d_filt |> dplyr::filter(.data$run_seq_num >= plot_range_indices[1] & .data$run_seq_num <= plot_range_indices[2]) |> droplevels()
  }


  if (!variable %in% names(data@dataset)){
    if(str_detect(variable, "_raw"))
      cli::cli_abort(cli::col_red("`{str_remove(variable, '_raw')}` values have not yet been drift and/or batch corrected. Please apply drift and/or batch corrections,  or choose another variable."))
    else
      cli::cli_abort(cli::col_red("`{variable}` values are not available. Please process accordingly or choose another variable."))
  }

  d_filt <- d_filt |>
    mutate(value = ifelse(is.infinite(!!variable_sym), NA, !!variable_sym))

  # Define y axis label based on if cap_outlier was selected
  y_label <- dplyr::if_else(cap_outliers,
                            paste0(ifelse(is.na(y_label_text), stringr::str_remove(variable, "feature\\_"), y_label_text), " (capped by MAD outlier filter) "),
                            stringr::str_remove(variable, "feature\\_"))

  # Subset data based on incl and excl argument values ----

  if (all(!is.na(include_feature_filter)) & all(include_feature_filter != "")) {
    if (length(include_feature_filter) == 1) {
      d_filt <- d_filt |> dplyr::filter(stringr::str_detect(.data$feature_id, include_feature_filter))
    } else {
      d_filt <- d_filt |> dplyr::filter(.data$feature_id %in% include_feature_filter)
    }
  }

  if (all(!is.na(exclude_feature_filter)) & all(exclude_feature_filter != "")) {
    if (length(exclude_feature_filter) == 1) {
      d_filt <- d_filt |> dplyr::filter(!stringr::str_detect(.data$feature_id, exclude_feature_filter))
    } else {
      d_filt <- d_filt |> dplyr::filter(!.data$feature_id %in% exclude_feature_filter)
    }
  }

  # Subset data based on qc_types argument ----
  if (all(!is.na(qc_types)) & all(qc_types != "")) {
    if (length(qc_types) == 1) {
      d_filt <- d_filt |> dplyr::filter(stringr::str_detect(.data$qc_type, qc_types))
    } else {
      d_filt <- d_filt |> dplyr::filter(.data$qc_type %in% qc_types)
    }
  }

  # Error in case the set filters excluded all features

  if (nrow(d_filt) < 1) cli::cli_abort("None of the feature id's meet the filter criteria. Please check `include_feature_filter` and `exclude_feature_filter` parameters.")

  # Re-order qc_type levels to define plot layers, i.e.  QCs are plotted over StudySamples
  d_filt$qc_type <- factor(as.character(d_filt$qc_type), pkg.env$qc_type_annotation$qc_type_levels)

  # Cap upper outliers  ----

  d_filt <- d_filt |>
    dplyr::mutate(value = !!variable_sym)

  if (cap_outliers) {

    if (!is.na(cap_top_n_values) & cap_top_n_values > 0) {
      d_filt <- d_filt |>
        dplyr::group_by(.data$feature_id) |>
        dplyr::arrange(.data$value) |>
        mutate(
          value = ifelse(dplyr::row_number() < cap_top_n_values, .data$value[cap_top_n_values], .data$value)
        ) |>
        dplyr::arrange(.data$feature_id, .data$run_seq_num) |>
        dplyr::ungroup()
    }

    d_filt <- d_filt |>
      dplyr::group_by(.data$feature_id) |>
      dplyr::mutate(
        value_max_spl = median(.data$value[.data$qc_type=="SPL"], na.rm=T) + cap_spl_k_mad * mad(.data$value[.data$qc_type=="SPL"], na.rm=T),
        value_max_tqc = median(.data$value[.data$qc_type=="TQC"], na.rm=T) + cap_qc_k_mad * mad(.data$value[.data$qc_type=="TQC"], na.rm=T),
        value_max_bqc = median(.data$value[.data$qc_type=="BQC"], na.rm=T) + cap_qc_k_mad * mad(.data$value[.data$qc_type=="BQC"], na.rm=T),
        value_max = pmax(.data$value_max_spl, .data$value_max_tqc, .data$value_max_bqc, na.rm = TRUE),
        value_max = ifelse(is.infinite(.data$value_max), NA, .data$value_max),
        value_mod = dplyr::if_else(
         .data$value > .data$value_max,
          if (!all(is.na(.data$value))) {
            max(.data$value[.data$value <= .data$value_max], na.rm = TRUE)
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


  if (output_pdf & (is.na(path) | path == "")) cli::cli_abort("Save to PDF selected, but no valid path defined. Please set path via `path`.")

  if (output_pdf & !is.na(path)) {
    path <- ifelse(stringr::str_detect(path, ".pdf"), path, paste0(path, ".pdf"))
    if (page_orientation == "LANDSCAPE") {
      pdf(file = path, onefile = T, paper = "A4r", useDingbats = FALSE, width = 28 / 2.54, height = 20 / 2.54)
    } else {
      pdf(file = path, onefile = T, paper = "A4", useDingbats = FALSE, height = 28 / 2.54, width = 20 / 2.54)
    }
  }

  if (is.na(specific_page)) {
    page_range <- 1:ceiling(dplyr::n_distinct(d_filt$feature_id) / (cols_page * rows_page))
  } else {
    page_range <- specific_page
  }

  if(output_pdf) action_text = "Saving plots to pdf" else action_text = "Generating plots"
  #if(show_progress) cli::cli_progress_bar(action_text, total = max(page_range))

  cat(cli::col_green(glue::glue("{action_text} ({max(page_range)} {ifelse(max(page_range) > 1, 'pages', 'page')}){ifelse(show_progress, ':', '...')}")))
  if(show_progress) pb <- txtProgressBar( min = 0, max = max(page_range), width = 50, style = 3)

   p_list <- list()  # p_list <- vector("list", length(page_range))
  for (i in page_range) {

    p <- runscatter_one_page(
      d_filt = d_filt, data = data, y_var = variable, d_batches = data@annot_batches, cols_page = cols_page, rows_page = rows_page, show_trend = show_trend,
      fit_qc_type = fit_qc_type, output_pdf = output_pdf, page_no = i,
      point_size = point_size, cap_outliers = cap_outliers, point_transparency = point_transparency, annot_scale = annot_scale,
      show_batches = show_batches, batches_as_shades = batches_as_shades, batch_line_color = batch_line_color,
      batch_shading_color = batch_shading_color, y_label = y_label, base_font_size = base_font_size, point_border_width = point_border_width, show_grid = show_gridlines,
      log_scale = log_scale, plot_range_indices = plot_range_indices, show_control_limits = show_control_limits, set_limits_n_sd = set_limits_n_sd, limits_batchwise = limits_batchwise, limits_linecolor = limits_linecolor,
      limits_linewidth = limits_linewidth, trend_linecolor = trend_linecolor
    )

    plot(p)
    dev.flush()
    flush.console()
    if(show_progress) setTxtProgressBar(pb, i)
    p_list[[i]] <- p

  }
  if (output_pdf) dev.off()
  cat(cli::col_green(" - done!"))
  if(show_progress) close(pb)

  flush.console()
  # on.exit(
  #   dev.off())

  if (return_plots) {
    return(p_list)
  }
}



runscatter_one_page <- function(d_filt, data, y_var, d_batches, cols_page, rows_page, page_no,
                                show_trend, fit_qc_type, cap_outliers,
                                show_batches, batches_as_shades, batch_line_color, batch_shading_color,
                                output_pdf, annot_scale, point_transparency, point_size = point_size, y_label, base_font_size, point_border_width,
                                show_grid, log_scale, plot_range_indices, show_control_limits, set_limits_n_sd, limits_batchwise, limits_linecolor,
                                limits_linewidth, trend_linecolor) {
  point_size <- ifelse(missing(point_size), 2, point_size)
  point_border_width <- dplyr::if_else(output_pdf, .3, .2 * (1 + annot_scale / 5))


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
            max(.data$value_mod, na.rm = TRUE) * 1.0
          } else {
            NA_real_
          }
      )
  }


  d_batch_data <- d_batches |> dplyr::slice(rep(1:dplyr::n(), each = nrow(dMax)))
  d_batch_data$feature_id <- rep(dMax$feature_id, times = nrow(d_batches))
  d_batch_data <- d_batch_data |> dplyr::left_join(dMax, by = c("feature_id"))

  p <- ggplot2::ggplot(d_subset, ggplot2::aes_string(x = "run_seq_num"))

  # browser()
  if (show_batches) {
    if (!batches_as_shades) {
      d_batches_temp <- d_batch_data |> filter(.data$id_batch_start != 1)
      p <- p + ggplot2::geom_vline(data = d_batches_temp, ggplot2::aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, linetype = "solid", size = .5, na.rm = TRUE)
    } else {
      d_batches_temp <- d_batch_data |> dplyr::filter(.data$batch_no %% 2 != 1)
      p <- p + ggplot2::geom_rect(
        data = d_batches_temp, ggplot2::aes(xmin = .data$id_batch_start - 0.5, xmax = .data$id_batch_end + 0.5, ymin = 0, ymax = .data$y_max),
        inherit.aes = FALSE, fill = batch_shading_color, color = NA, alpha = 1, linetype = "solid", size = 0.3, na.rm = TRUE
      )
    }
  }

  if (cap_outliers) {
    p <- p + ggplot2::geom_hline(data = dMax, ggplot2::aes(yintercept = .data$y_max), color = "#ffefbf", size = 3, alpha = 1)
  }

  p <- p +
    ggplot2::geom_point(aes_string(x = "run_seq_num", y = "value_mod", color = "qc_type", fill = "qc_type", shape = "qc_type", group = "batch_id"), size = point_size, alpha = point_transparency, stroke = point_border_width, na.rm = TRUE)


  if (show_trend) {
    #browser()
    y_var_trend <- if_else(str_detect(y_var, "\\_raw"), paste0(y_var, "_fit"), paste0(y_var, "_fit_after"))
    p <- p +
      ggplot2::geom_line(aes_string(x = "run_seq_num", y = y_var_trend, group = "batch_id"), color = trend_linecolor, size = 1, na.rm = TRUE)
      # pkg.env$qc_type_annotation$qc_type_fillcol[fit_qc_type]
  }

  p <- p +
    ggh4x::facet_wrap2(ggplot2::vars(.data$feature_id), scales = "free_y", ncol = cols_page, nrow = rows_page, trim_blank = FALSE) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
    ggplot2::scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
    ggplot2::scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)



  # if (show_trend) {
  #   if (show_before_after_smooth) {
  #     p <- p +
  #       ggplot2::geom_smooth(
  #         data = filter(d_subset, .data$qc_type == fit_qc_type), ggplot2::aes_string(x = "run_seq_num", y = "y_fit", group = "batch_id"), se = TRUE, na.rm = TRUE,
  #         colour = pkg.env$qc_type_annotation$qc_type_fillcol[fit_qc_type], fill = pkg.env$qc_type_annotation$qc_type_fillcol[fit_qc_type],
  #         method = "loess", alpha = 0.35, size = 0.8 # TODO before used MASS::rlm
  #       )
  #
  #     if (show_trend_samples) {
  #       p <- p + ggplot2::geom_smooth(
  #         data = filter(d_subset, .data$qc_type == "SPL"), ggplot2::aes_string(x = "run_seq_num", y = "value", group = "batch_id"), colour = trend_samples_col, fill = trend_samples_col, linetype = "dashed",
  #         method = trend_samples_function, se = FALSE, alpha = 0.2, size = .8, na.rm = TRUE
  #       )
  #     }
  #
  #     if (plot_other_qc) {
  #       other_qc <- dplyr::if_else(fit_qc_type == "BQC", "TQC", "BQC")
  #       other_qc_col <- pkg.env$qc_type_annotation$qc_type_fillcol[other_qc]
  #       p <- p +
  #         ggplot2::geom_smooth(
  #           data = dplyr::filter(d_subset, .data$qc_type == other_qc), ggplot2::aes_string(x = "run_seq_num", y = "value", group = "batch_id"), colour = other_qc_col, fill = other_qc_col,
  #           method = trend_samples_function, se = TRUE, alpha = 0.2, size = .4, na.rm = FALSE
  #         )
  #     }
  #   } else {
  #     p <- p +
  #       geom_smooth(
  #         data = dplyr::filter(d_subset, .data$qc_type == "SPL"), ggplot2::aes_string(x = "run_seq_num", y = "y_fit", group = "batch_id"), colour = trend_samples_col, fill = trend_samples_col,
  #         method = trend_samples_function, se = TRUE, alpha = 0.2, size = .4, na.rm = FALSE
  #       )
  #
  #     if (plot_other_qc) {
  #       other_qc <- dplyr::if_else(.data$fit_qc_type == "BQC", "TQC", "BQC")
  #       other_qc_col <- pkg.env$qc_type_annotation$qc_type_fillcol[other_qc]
  #       p <- p +
  #         ggplot2::geom_smooth(
  #           data = dplyr::filter(d_subset, .data$qc_type == other_qc), ggplot2::aes_string(x = "run_seq_num", y = "value", group = "batch_id"), colour = other_qc_col, fill = other_qc_col,
  #           method = trend_samples_function, se = TRUE, alpha = 0.2, size = .4, na.rm = FALSE
  #         )
  #     }
  #   }
  # }

  if(!all(is.na(show_control_limits))) {
    if(limits_batchwise) grp = c("feature_id", "batch_id") else grp = c("feature_id")
    d_subset_stats <- d_subset |>
      left_join(d_batches, by = c("batch_id")) |>
      filter(.data$qc_type %in% show_control_limits) |>
       group_by(dplyr::pick(grp)) |>
       summarise(mean = mean(.data$value_mod, na.rm = TRUE),
                 sd = set_limits_n_sd * sd(.data$value_mod, na.rm = TRUE),
                 batch_start = min(.data$id_batch_start),
                 batch_end = max(.data$id_batch_end),
                 .groups = 'drop')
    if(limits_batchwise){
      p <- p +
        ggplot2::geom_segment(data = d_subset_stats, inherit.aes = FALSE, aes(x = .data$batch_start, xend = .data$batch_end, y = .data$mean, yend = .data$mean, group = .data$batch_id), color = limits_linecolor, size = limits_linewidth, linetype = "solid", alpha = 1) +
        ggplot2::geom_rect(data = d_subset_stats, inherit.aes = FALSE, aes(xmin = .data$batch_start, xmax = .data$batch_end, ymin = .data$mean - .data$sd , ymax = .data$mean + .data$sd, group = .data$batch_id), fill = limits_linecolor, size = limits_linewidth, alpha = .25)
        #ggplot2::geom_segment(data = d_subset_stats, inherit.aes = FALSE, aes(x = .data$batch_start, xend = .data$batch_end, y = .data$mean - .data$sd , yend = .data$mean - .data$sd, group = .data$batch_id), color = limits_linecolor, linetype = "solid", size = limits_linewidth, alpha = 1)
    } else{
      p <- p +
        ggplot2::geom_hline(data = d_subset_stats, aes(yintercept = .data$mean), color = limits_linecolor, size = limits_linewidth, alpha = 1, linetype = "longdash") +
        ggplot2::geom_hline(data = d_subset_stats, aes(yintercept = .data$mean + .data$sd ), color = limits_linecolor, size = limits_linewidth, alpha = 1, linetype = "dashed") +
        ggplot2::geom_hline(data = d_subset_stats, aes(yintercept = .data$mean - .data$sd), color = limits_linecolor, size = limits_linewidth, alpha = 1, linetype = "dashed")
    }

  }

  p <- p +
    # aes(ymin=0) +
    ggplot2::xlab("Analysis order") +
    ggplot2::ylab(label = y_label) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.02, 0.03))) +
    # expand_limits(y = 0) +
    ggplot2::theme_bw(base_size = base_font_size) +

    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 1, face = "bold"),
      strip.text = ggplot2::element_text(size = 9 * annot_scale, face = "bold"),
      strip.background = ggplot2::element_rect(size = 0.0001, fill = "#00283d"),
      strip.text.x = ggplot2::element_text(color = "white"),
      axis.text.x = ggplot2::element_text(size = 7 * annot_scale, face = NULL),
      axis.text.y = ggplot2::element_text(size = 7 * annot_scale, face = NULL),
      axis.title = ggplot2::element_text(size = 8 * annot_scale, face = NULL),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      strip.switch.pad.wrap = ggplot2::unit(-1, "mm"),
      panel.border = element_rect(linewidth = 0.5, color = "grey40")
    )

  if (show_grid) {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_line(size = 0.3, colour = "grey88", linetype = "dashed"))
  } else {
    p <- p + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  }


  if (log_scale) {
    p <- p + scale_y_log10()
  }

  if(!all(is.na(plot_range_indices))) {
    p <- p + ggplot2::coord_cartesian(xlim = plot_range_indices)
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
#' @param filter_data Logical, whether to use QC-filtered data
#' @param qc_types QC type(s) to plot. Can be a vector of QC types or `NA`/`NULL` to plot all available QC types.
#' @param x_axis_variable Character, variable to use for the x-axis. Must be one of: "run_seq_num", "run_no", "analysis_id", or "timestamp"
#' @param include_feature_filter Regex pattern to select specific features
#' @param exclude_feature_filter Regex pattern to exclude specific features
#' @param plot_range_indices Numeric vector of length 2, specifying the start and end indices of the sequence to be plotted. `NA` plots all samples.
#' @param min_feature_intensity Numeric, exclude features with overall median signal below this value
#' @param y_lim Numeric vector of length 2, specifying the lower and upper y-axis limits. Default is `NA`, which uses limits calculated based on `ignore_outliers`.
#' @param ignore_outliers Logical, whether to exclude outlier values based on 4x MAD (median absolute deviation) fences
#' @param show_batches Logical, whether to show batch separators in the plot
#' @param batches_as_shades Logical, whether to show batches as shaded areas instead of line separators
#' @param batch_line_color Character, color of the batch separator lines
#' @param batch_shading_color Character, color of the batch shaded areas
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
qc_plot_rla_boxplot <- function(
                                data = NULL,
                                rla_type_batch = c("within", "across"),
                                variable = c("intensity", "norm_intensity", "conc", "conc_raw", "area", "height", "fwhm"),
                                filter_data = FALSE,
                                qc_types = NA,
                                x_axis_variable = c("run_seq_num"),
                                include_feature_filter = "",
                                exclude_feature_filter = "",
                                plot_range_indices = NA,
                                min_feature_intensity = 0,
                                y_lim = NA,
                                ignore_outliers = FALSE,
                                show_batches = TRUE,
                                batches_as_shades = FALSE,
                                batch_line_color = "#b6f0c5",
                                batch_shading_color = "grey93",
                                x_gridlines = FALSE,
                                linewidth = 0.2,
                                base_font_size = 8,
                                relative_log_abundances = TRUE) {
  check_data(data)
  if (nrow(data@dataset) < 1) cli::cli_abort("No data available. Please import data and metadata first.")

  # Check if selected variable is valid
  rlang::arg_match(variable, c("area", "height", "intensity", "response", "conc", "conc_raw", "rt", "fwhm"))
  variable <- str_remove(variable, "feature_")
  variable <- stringr::str_c("feature_", variable)
  check_var_in_dataset(data@dataset, variable)
  variable_sym = rlang::sym(variable)
  rlang::arg_match(x_axis_variable, c("run_seq_num", "run_no", "analysis_id", "timestamp"))
  if(x_axis_variable == "run_no") x_axis_variable <- "run_seq_num"
  if(x_axis_variable == "timestamp") x_axis_variable <- "acquisition_time_stamp"
  x_axis_variable_sym <- rlang::sym(x_axis_variable)
  rlang::arg_match(rla_type_batch, c("within", "across"))

  # Filter data if filter_data is TRUE
  if (filter_data) {
    d_filt <- data@dataset_filtered |> dplyr::ungroup()
    if (!data@is_filtered) cli::cli_abort(cli::col_red("Data has not been qc filtered, or has changed. Please run `qc_apply_feature_filter` first."))
  } else {
    d_filt <- data@dataset |> dplyr::ungroup()
  }

  # Subset data based on qc_types argument ----
  if (all(!is.na(qc_types)) & all(qc_types != "")) {
    if (length(qc_types) == 1) {
      d_filt <- d_filt |> dplyr::filter(stringr::str_detect(.data$qc_type, qc_types)) |> droplevels()
    } else {
      d_filt <- d_filt |> dplyr::filter(.data$qc_type %in% qc_types) |> droplevels()
    }
  }

  if (!all(is.na(plot_range_indices))) {
    d_filt <- d_filt |> dplyr::filter(.data$run_seq_num >= plot_range_indices[1] & .data$run_seq_num <= plot_range_indices[2]) |> droplevels()
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
      group_by(dplyr::pick(grp)) |>
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
    if (!batches_as_shades) {
      p <- p + geom_vline(data = data@annot_batches |> slice(-1), aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, linetype = "solid", size = 1)
    } else {
      d_batch_2nd <- data@annot_batches |>
        slice(-1) |>
        filter(.data$batch_no %% 2 != 1)
      p <- p + geom_rect(
        data = d_batch_2nd, aes(xmin = .data$id_batch_start - 0.5, xmax = .data$id_batch_end + 0.5, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = batch_shading_color, color = NA, alpha = 1, linetype = "solid", size = 0.5, na.rm = TRUE
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
    p <- p + theme(panel.grid.major.x = element_line(colour = "#bdbdbd", linetype = "dotted", size = .3))
  else
    p <- p + theme(panel.grid.major.x = element_blank())

  if (relative_log_abundances) {
    p <- p + geom_hline(yintercept = 0, colour = "#5fe3f5", linetype = "longdash", size = 0.5) +
      ylab(bquote(bold(log[2] ~ "( relative" ~ .(stringr::str_remove(variable, "feature\\_")) ~ ")")))
  }
  ylim = c(NA, NA)
  xlim = c(NA, NA)
  # Set y-axis limits
  if (!all(is.na(y_lim))) {
    ylim = y_lim
  } else if(ignore_outliers) {
    tails <- get_mad_tails(d_filt$val_res, k = 4)
    ylim = tails
  }

  if(!all(is.na(plot_range_indices))) {
    xlim = plot_range_indices
  }

  p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)

  return(p)
}
