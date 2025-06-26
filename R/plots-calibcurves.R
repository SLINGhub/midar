#' Plot Calibration Curves
#'
#' This function plots calibration curves of each feature where defined
#' and displays QC samples with defined concentrations within the plot.
#' Users can select a regression model (`linear` or `quadratic`) and apply
#' weighting (`none`, `"1/x"`, or `"1/x^2"`), either through function arguments
#' or feature metadata.
#'
#' Features for plotting can be filtered using QC filters defined via
#' [filter_features_qc()] or through `include_feature_filter` and
#' `exclude_feature_filter` arguments. The resulting plots offer extensive
#' customization options, including point size, line width, point color, point
#' fill, point shape, line color, ribbon fill, and font base size.
#'
#' Plots will be divided into multiple pages if the number of features exceeds
#' the product of `rows_page` and `cols_page` settings. The function supports
#' both direct plotting within R and saving plots as PDF files. Additionally,
#' plots can be returned as a list of ggplot2 objects for further manipulation
#' or integration into other analyses.
#'
#' @param data A `MidarExperiment` object containing the dataset.
#' @param variable Variable to plot on the y-axis, usually intensity. Default
#'   is `"intensity"`.
#' @param qc_types A character vector specifying the QC types to plot. It must contain
#' at least `CAL`, which represents calibration curve samples. Other QC types will be
#' plotted as points when they have assigned concentrations (see QC-concentration metadata).
#' These QC types need to be present in the data and defined in the analysis metadata.
#' The default is `NA`, which means any of the QC types "CAL", "HQC", "MQC", "LQC",
#' "EQA", "QC", will be plotted if present and have assigned concentrations.
#' @param fit_overwrite If `TRUE`,
#'   the function will use the provided `fit_model` and `fit_weighting` values
#'   for all analytes and ignore any fit method and weighting settings defined in
#'   the metadata .
#' @param fit_model A character string specifying the default regression fit
#'   method to use for the calibration curve. Must be one of `"linear"` or
#'   `"quadratic"`. This method will be applied if no specific fit method is
#'   defined for a feature in the metadata, or
#'   when `fit_overwrite = TRUE`.
#' @param fit_weighting A character string specifying the default weighting
#'   method for the regression points in the calibration curve. Must be one of
#'   `"none"`, `"1/x"`, or `"1/x^2"`. This method will be applied if no
#'   specific weighting method is defined for a feature in the metadata, or
#'   when `fit_overwrite = TRUE`.
#' @param show_confidence_interval Logical, if `TRUE`, displays the confidence interval as ribbon.
#' Default is `NA`, in which case confidence intervals are plotted in a linear
#' scale and ommitted in log-log scale.
#' @param log_axes Logical. Determines whether the x and y axes are displayed in a logarithmic scale (log-log scale).
#'   Set to `TRUE` to enable logarithmic scaling; otherwise, set to `FALSE` for a linear scale.
#'   Note: If `TRUE`, any regression curves or standard error regions with negative
#'   values will be omitted from display.
#' @param filter_data Logical, if `TRUE`, uses QC filtered data; otherwise uses
#'   raw data. Default is `FALSE`.
#' @param include_qualifier Logical, whether to include qualifier features. Default is `TRUE`.
#' @param include_istd Logical, whether to include internal standard (ISTD) features. Default is `TRUE`.
#' @param include_feature_filter Regex pattern to include features. If omitted,
#'   considers all features.
#' @param exclude_feature_filter Regex pattern to exclude features. If omitted,
#'   excludes none.
#' @param output_pdf Logical, if `TRUE`, saves plots as a PDF file. Default is
#'   `FALSE`.
#' @param path File path for saving the PDF. Default is an empty string.
#' @param return_plots Logical, if `TRUE`, returns plots as a list of ggplot2
#'   objects. Default is `FALSE`.
#' @param point_size Size of points in the plot. Default is 1.5.
#' @param line_width Width of regression lines. Default is 0.7.
#' @param point_color A vector specifying the colors for points corresponding
#'   to different QC types. This can be either an unnamed vector or a named
#'   vector, with names corresponding to QC types. Unused colors will be ignored.
#'   Default is `NA` which corresponds to the default colors for QC types defined in the package.
#' @param point_fill A vector specifying the fill colors for points corresponding
#'   to different QC types. This can be either an unnamed vector or a named
#'   vector, with names corresponding to QC types. Unused fill colors will be ignored.
#'   Default is `NA` which corresponds to the default fill colors for QC types defined in the package.
#' @param point_shape A vector specifying the shapes for points corresponding
#'   to different QC types. This can be either an unnamed vector or a named
#'   vector, with names corresponding to QC types. Unused shapes will be ignored.
#'   Default is `NA` which corresponds to the default shapes for QC types defined in the package.

#'
#' @param line_color Color of the regression line. Default is `"#4575b4"`.
#' @param ribbon_fill Color for the confidence interval ribbon. Default is
#'   `"#91bfdb40"`.
#' @param font_base_size Base font size for text in plots. Default is 7.
#' @param rows_page Number of plot rows. Default is 4.
#' @param cols_page Number of plot columns. Default is 5.
#' @param specific_page Show/save a specific page number only. `NA` plots/saves all pages.
#' @param page_orientation Orientation of PDF, either `"LANDSCAPE"` or
#'   `"PORTRAIT"`. Default is `"LANDSCAPE
#' @param show_progress Logical. If `TRUE`, displays a progress bar during
#'   plot creation.
#'
#' @export
plot_calibrationcurves <- function(data = NULL,
                                   variable = "norm_intensity",
                                   qc_types = NA,
                                   fit_overwrite,
                                   fit_model = c("linear", "quadratic"),
                                   fit_weighting = c(NA, "none", "1/x", "1/x^2"),
                                   show_confidence_interval = NA,
                                   log_axes = FALSE,
                                   filter_data = FALSE,
                                   include_qualifier = TRUE,
                                   include_istd = FALSE,
                                   include_feature_filter = NA,
                                   exclude_feature_filter = NA,
                                   output_pdf = FALSE,
                                   path = NA,
                                   return_plots = FALSE,
                                   point_size = 1.5,
                                   line_width = 0.7,
                                   point_color = NA,
                                   point_fill = NA,
                                   point_shape = NA,
                                   line_color = "#4575b4",
                                   ribbon_fill = "#91bfdb40",
                                   font_base_size = 7,
                                   rows_page = 4,
                                   cols_page = 5,
                                   specific_page = NA,
                                   page_orientation = "LANDSCAPE",
                                   # Progress bar settings
                                   show_progress = TRUE) {
  check_data(data)

  variable_strip <- str_remove(variable, "feature_")
  rlang::arg_match(
    variable_strip,
    c(
      "area",
      "height",
      "intensity",
      "norm_intensity",
      "response",
      "conc",
      "conc_raw",
      "rt",
      "fwhm"
    )
  )
  variable <- stringr::str_c("feature_", variable_strip)
  variable_sym <- rlang::sym(variable)
  rlang::arg_match(fit_model, c(NA, "linear", "quadratic"))
  rlang::arg_match(fit_weighting, c(NA, "none", "1/x", "1/x^2"))
  # rlang::arg_match(fit_overwrite, c("TRUE", "FALSE"))

  plot_var <- rlang::sym(variable)

  if (is.na(show_confidence_interval)) {
    show_confidence_interval <- !log_axes
  }


  if (all(is.na(point_color))) {
    point_color <- pkg.env$qc_type_annotation$qc_type_col
  }
  if (all(is.na(point_fill))) {
    point_fill <- pkg.env$qc_type_annotation$qc_type_fillcol
  }
  if (all(is.na(point_shape))) {
    point_shape <- pkg.env$qc_type_annotation$qc_type_shape
  }


  check_var_in_dataset(data@dataset, variable)

  # Ensure path is defined when saving PDF
  if (output_pdf && (is.na(path) || path == "")) {
    cli::cli_abort(
      "The argument {.strong `path`} must be defined when {.strong output_pdf} is {.strong TRUE}."
    )
  }
  if (nrow(data@dataset) < 1) {
    cli::cli_abort(
      col_red(
        "No data available in the dataset. Please import the necessary data and metadata before proceeding."
      )
    )
  }

  if (nrow(data@annot_qcconcentrations) < 1) {
    cli::cli_abort(
      col_red(
        "No QC-concentration metadata is available. Please import the corresponding metadata to proceed."
      )
    )
  }
  if (!"CAL" %in% unique(data@dataset$qc_type)) {
    cli::cli_abort(
      col_red(
        "No QC type {.strong 'CAL'} defined in the data. Please assign `CAL` as `qc_type` to corresponding calibration analyses/samples in the analysis metadata."
      )
    )
  }

  # Subset dataset according to filter arguments

  d_filt <- get_dataset_subset(
    data,
    filter_data = filter_data,
    qc_types = qc_types,
    include_qualifier = include_qualifier,
    include_istd = include_istd,
    include_feature_filter = include_feature_filter,
    exclude_feature_filter = exclude_feature_filter
  )

  d_calib <- d_filt |>
    dplyr::select(any_of(
      c(
        "analysis_id",
        "sample_id",
        "qc_type",
        "feature_id",
        "analyte_id",
        variable
      )
    )) |>
    # REMOVE filter(str_detect(.data$qc_type, "CAL|[MLH]QC|^QC|EQA")) |>
    dplyr::right_join(
      data@annot_qcconcentrations,
      by = c("sample_id" = "sample_id", "analyte_id" = "analyte_id")
    ) |> drop_na("concentration")

  if (all(is.na(qc_types))) {
    qc_types <- unique(d_calib$qc_type)
  } else {
    if (length(setdiff(qc_types, unique(d_calib$qc_type))) > 0) {
      cli::cli_abort(cli::col_red(
        paste(
          "One or more selected `qc_types` have no defined analyte concentrations. Please verify the feature and QC-concentration metadata, or select other `qc_types`."
        )
      ))
    }
  }

  # If color_curves is provided, check if it has enough colors
  num_levels <- length(unique(d_calib$qc_type))
  if (length(point_color) < num_levels) {
    cli::cli_abort(cli::col_red(
      paste(
        "Insufficient colors in `point_colors`. Provide at least",
        num_levels,
        "unique colors for the number of selected `qc_types`"
      )
    ))
  }
  if (length(point_fill) < num_levels) {
    cli::cli_abort(cli::col_red(
      paste(
        "Insufficient fill colors in `point_fill`. Provide at least",
        num_levels,
        "unique colors for the number of selected `qc_types`"
      )
    ))
  }
  if (length(point_shape) < num_levels) {
    cli::cli_abort(cli::col_red(
      paste(
        "Insufficient shape codes in `point_shape`. Provide at least",
        num_levels,
        "unique shape codes for the number of selected `qc_types`"
      )
    ))
  }
  qc_types <- unique(d_calib$qc_type)
  if (!is.null(names(point_color)) &&
    length(setdiff(qc_types, names(point_color))) > 0) {
    cli::cli_abort(cli::col_red(
      paste(
        "The names in `point_color` must match the `qc_types` in the dataset. Please verify, or provide only colors."
      )
    ))
  }
  if (!is.null(names(point_fill)) &&
    length(setdiff(qc_types, names(point_fill))) > 0) {
    cli::cli_abort(cli::col_red(
      paste(
        "The names in `point_fill` must match the `qc_types` in the dataset. Please verify, or provide only fill colors."
      )
    ))
  }
  if (!is.null(names(point_shape)) &&
    length(setdiff(qc_types, names(point_shape))) > 0) {
    cli::cli_abort(cli::col_red(
      paste(
        "The names in `point_shape` must match the `qc_types` in the dataset. Please verify, or provide only shapes codes."
      )
    ))
  }

  d_calib$curve_id <- 1


  d_calib <- d_calib |>
    left_join(
      data@annot_features |> select("feature_id", "curve_fit_model", "curve_fit_weighting"),
      by = "feature_id"
    )


  d_calib <- d_calib |> mutate(
    curve_fit_model = ifelse(
      is.na(.data$curve_fit_model) |
        fit_overwrite,
      fit_model,
      .data$curve_fit_model
    ),
    fit_weighting = ifelse(
      is.na(.data$curve_fit_weighting) |
        fit_overwrite,
      fit_weighting,
      .data$curve_fit_weighting
    )
  )



  d_calib <- d_calib |> mutate(qc_type_cat = if_else(str_detect(.data$qc_type, "CAL"), "CAL", "QC"))



  # Verify if unit is the same for all data points/curves
  x_axis_unit <- unique(d_calib$concentration_unit)
  if (length(x_axis_unit) > 1) {
    cli::cli_abort(
      cli::col_red(
        "The `concentration_unit` (x axis) must be identical for selected curves and data points. Please change selection of curves or update calibration curve metadata."
      )
    )
  }


  # add calibration metrics to data
  data <- calc_calibration_results(
    data = data,
    variable = variable,
    include_qualifier = include_qualifier,
    fit_overwrite = fit_overwrite,
    fit_model = fit_model,
    fit_weighting = fit_weighting,
    include_fit_object = TRUE
  )



  count_regfailed <- sum(data@metrics_calibration$reg_failed_cal_1)
  if (count_regfailed > 0) {
    cli::cli_alert_info(
      col_yellow(
        "Regression failed for {count_regfailed} features, no curves shown for these."
      )
    )
  }

  get_predictions <- function(stats, log_scale) {
    if (!stats$reg_failed_cal_1) {
      fit <- stats$fit_cal_1[[1]]
      conc_orig <- model.matrix(fit)[, 2]
      if (log_scale) {
        concs <- 10^(seq(log10(min(conc_orig)), log10(max(conc_orig)), length.out = 100))
      } else {
        concs <- seq(min(conc_orig), max(conc_orig), length.out = 100)
      }
      predictions <- predict(fit,
        newdata = data.frame(concentration = concs),
        interval = "confidence"
      )
      prediction_data <- tibble(
        feature_id = stats$feature_id,
        curve_id = "1",
        concentration = concs,
        y_pred = predictions[, "fit"],
        lwr = predictions[, "lwr"],
        upr = predictions[, "upr"]
      )
    } else {
      prediction_data <- tibble(
        feature_id = stats$feature_id,
        curve_id = "1",
        concentration = NA_real_,
        y_pred = NA_real_,
        lwr = NA_real_,
        upr = NA_real_
      )
    }
    prediction_data
  }

  d_calib_stats <- data@metrics_calibration
  d_calib_stats_grp <- d_calib_stats |>
    dplyr::group_split(.data$feature_id) # TODO .data$curve_id

  d_pred <- map(d_calib_stats_grp, function(x) {
    get_predictions(x, log_axes)
  }) |> bind_rows()



  # Prepare PDF output
  if (output_pdf && !is.na(path)) {
    # nocov start
    path <- ifelse(stringr::str_detect(path, ".pdf"),
      path,
      paste0(path, ".pdf")
    )
    if (page_orientation == "LANDSCAPE") {
      pdf(
        file = path,
        onefile = TRUE,
        paper = "A4r",
        useDingbats = FALSE,
        width = 28 / 2.54,
        height = 20 / 2.54
      )
    } else {
      pdf(
        file = path,
        onefile = TRUE,
        paper = "A4",
        useDingbats = FALSE,
        height = 28 / 2.54,
        width = 20 / 2.54
      )
    }
  } # nocov end


  # Determine the range of pages to generate
  if (!rlang::is_na(specific_page)) {
    total_pages <- ceiling(n_distinct(d_calib$feature_id) /
      (cols_page * rows_page))
    if (specific_page > total_pages) {
      cli::cli_abort(
        col_red(
          "Selected page exceeds the total number of pages. Please select a page number between {.strong 1} and {.strong {total_pages}}."
        )
      )
    }
    page_range <- specific_page
  } else {
    page_range <- 1:ceiling(n_distinct(d_calib$feature_id) /
      (cols_page * rows_page))
  }

  a <- !all(is.na(d_pred$concentration))
  log_flag <- log_axes && a

  #TODO: COVR still does not detect this part, even confirm being tested
  if (log_flag) {
    txt <- "" # nocov start
    txt2 <- ""
    if (show_confidence_interval) {
      txt2 <- ifelse(
        show_confidence_interval,
        "Consider set `show_confidence_interval = FALSE`",
        ""
      )
      if (any(d_pred$y_pred <= 0, na.rm = T)) {
        txt <- "regression curve and confidence intervals"
      } else if (any(d_pred$lwr <= 0, na.rm = T)) {
        txt <- "regression confidence intervals"
      }
    } else {
      if (any(d_pred$y_pred <= 0, na.rm = T)) {
        txt <- "regression curve"
      }
    }

    if (txt != "") {
      cli::cli_alert_info(
        cli::col_yellow(
          "Regions of the {txt} are partially <= 0 and will be omitted. {txt2}."
        )
      )
    }

    d_pred <- d_pred |>
      group_by(.data$feature_id) |>
      mutate(
        y_pred = if_else(.data$y_pred < 0, NA_real_, .data$y_pred),
        lwr = if_else(.data$lwr < 0, NA_real_, .data$lwr)
      ) |>
      ungroup() # nocov end
  }

  # Action text for progress output
  action_text <- if (output_pdf) {
    "Saving plots to pdf"
  } else {
    "Generating plots"
  }
  page_suffix <- if (max(page_range) > 1) {
    glue::glue("{max(page_range)} pages")
  } else {
    glue::glue("{max(page_range)} page")
  }
  progress_suffix <- if (show_progress) {
    ":"
  } else {
    "..."
  }

  # Output formatted message
  message(glue::glue("{action_text} ({page_suffix}){progress_suffix}"))

  # Initialize progress bar if requested
  if (show_progress) {
    pb <- txtProgressBar(
      min = 0,
      max = max(page_range),
      width = 30,
      style = 3
    )
  }

  p_list <- list() # p_list <- vector("list", length(page_range))

  for (i in page_range) {
    p <- plot_calibcurves_page(
      d_pred = d_pred,
      d_calib = d_calib,
      d_calib_stats = d_calib_stats,
      output_pdf = output_pdf,
      response_variable = variable,
      include_qualifier = include_qualifier,
      path = path,
      rows_page = rows_page,
      cols_page = cols_page,
      specific_page = i,
      point_size = point_size,
      line_width = line_width,
      point_color = point_color,
      point_fill = point_fill,
      point_shape = point_shape,
      line_color = line_color,
      ribbon_fill = ribbon_fill,
      font_base_size = font_base_size,
      x_axis_title = x_axis_unit,
      fit_model = fit_model,
      fit_weighting = fit_weighting,
      log_axes = log_axes,
      show_confidence_interval = show_confidence_interval,
      fit_overwrite = fit_overwrite
    )
    plot(p)
    dev.flush()
    flush.console()
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
    p_list[[i]] <- p
  }

  if (output_pdf) {
    dev.off()
  } # Close PDF device
  message(" - done!") # Completion message
  if (show_progress) {
    close(pb)
  } # Close progress bar if open

  # Return plot list or invisible
  if (return_plots) {
    return(p_list[page_range])
  } else {
    invisible()
  }
}






# Define function to plot 1 page
plot_calibcurves_page <- function(d_pred,
                                  d_calib,
                                  d_calib_stats,
                                  output_pdf,
                                  response_variable,
                                  include_qualifier,
                                  path,
                                  rows_page,
                                  cols_page,
                                  specific_page,
                                  point_size,
                                  line_width,
                                  point_color,
                                  point_fill,
                                  point_shape,
                                  line_color,
                                  ribbon_fill,
                                  font_base_size,
                                  x_axis_title,
                                  fit_model,
                                  fit_weighting,
                                  log_axes,
                                  show_confidence_interval,
                                  fit_overwrite) {
  plot_var <- rlang::sym(response_variable)
  d_calib$curve_id <- as.character(d_calib$curve_id)

  # subset the dataset with only the rows used for plotting the facets of the selected page
  n_samples <- length(unique(d_calib$analysis_id))
  row_start <- n_samples * cols_page * rows_page * (specific_page - 1) + 1
  row_end <- n_samples * cols_page * rows_page * specific_page

  dat_subset <- d_calib |>
    dplyr::arrange(.data$feature_id, .data$curve_id) |>
    dplyr::slice(row_start:row_end) |>
    mutate(!!plot_var := if_else(is.nan(!!plot_var), NA_real_, !!plot_var)) |>
    drop_na((!!plot_var))

  d_pred_filt <- d_pred |>
    dplyr::semi_join(dat_subset, by = c("feature_id")) |>
    dplyr::arrange(.data$feature_id, .data$curve_id) |>
    group_by(.data$feature_id) %>%
    filter(!(all(is.na(lwr)) | all(is.na(upr)))) |>
    ungroup()

  d_pred_sum <- d_pred_filt |>
    dplyr::group_by(.data$feature_id) |>
    dplyr::summarise(
      x_min = safe_min(.data$concentration, na.rm = TRUE),
      y_max = if (show_confidence_interval) {
        safe_max(.data$upr, na.rm = TRUE)
      } else {
        safe_max(.data$y_pred, na.rm = TRUE)
      },
      y_min = safe_min(.data$y_pred, na.rm = TRUE)
    )

  d_calib_stats <- d_calib_stats |>
    dplyr::semi_join(dat_subset, by = c("feature_id")) |>
    dplyr::inner_join(d_pred_sum, by = c("feature_id"))

  p <- ggplot(data = dat_subset, aes(x = .data$concentration, y = !!plot_var))

  if (show_confidence_interval &&
    nrow(d_pred_filt |> filter(!is.na(.data$concentration))) > 0) {
    # Confidence interval
    p <- p +
      ggplot2::geom_ribbon(
        data = d_pred_filt,
        aes(
          x = .data$concentration,
          ymin = .data$y_pred,
          ymax = .data$upr
        ),
        fill = ribbon_fill,
        inherit.aes = FALSE,
        na.rm = TRUE
      ) +
      ggplot2::geom_ribbon(
        data = d_pred_filt,
        aes(
          x = .data$concentration,
          ymin = .data$lwr,
          ymax = .data$y_pred
        ),
        fill = ribbon_fill,
        inherit.aes = FALSE,
        na.rm = TRUE
      )
  }
  p <- p +
    ggplot2::geom_line(
      data = d_pred_filt,
      aes(x = .data$concentration, y = .data$y_pred),
      color = line_color,
      linewidth = line_width,
      inherit.aes = FALSE,
      na.rm = TRUE
    ) +


    # color = ifelse(after_stat(r.squared) < 0.80, "red", "darkgreen")), size = 1.4) +
    scale_color_manual(values = point_color) +
    scale_fill_manual(values = point_fill) +
    scale_shape_manual(values = point_shape) +
    # scale_y_continuous(limits = c(0, NA)) +
    # scale_x_continuous(limits = c(0, NA), breaks = scales::breaks_extended(6)) +
    ggh4x::facet_wrap2(
      vars(.data$feature_id),
      scales = "free",
      nrow = rows_page,
      ncol = cols_page,
      trim_blank = FALSE
    ) +
    geom_point(
      aes(
        fill = .data$qc_type,
        color = .data$qc_type,
        shape = .data$qc_type
      ),
      size = point_size,
      na.rm = TRUE
    ) +

    labs(
      x = x_axis_title,
      y = stringr::str_remove(response_variable, "feature\\_")
    ) +
    theme_light(base_size = font_base_size) +
    theme(
      strip.text = element_text(size = font_base_size, face = "bold"),
      strip.background = element_rect(linewidth = 0.0001, fill = "#496875"),
      panel.grid.major = element_line(
        color = "grey70",
        linewidth = 0.2,
        linetype = "dotted"
      ),
      # Light and dotted major gridlines
      panel.grid.minor = element_line(
        color = "grey90",
        linewidth = 0.1,
        linetype = "dotted"
      ) # Lighter minor gridlines
    )



  if (log_axes) {
    p <- p +
      scale_x_continuous(trans = "log10") +
      scale_y_continuous(trans = "log10")
  }

  if (!log_axes) {
    p <- p + ggplot2::coord_cartesian(xlim = c(0, NA), ylim = c(0, NA))
  }

  if (nrow(d_pred |> filter(!is.na(.data$concentration))) > 0) {
    p <- p +
      ggplot2::geom_text(
        data = d_calib_stats,
        aes(
          x = if (!log_axes) 0 else .data$x_min,
          y = Inf,
          label = glue::glue(
            "{stringr::str_to_title(fit_model)}\nR\u00B2 = {round(.data$r2_cal_1, 4)}"
          )
        ),
        inherit.aes = FALSE,
        hjust = if (!log_axes) -0.1 else -0.2,
        vjust = 1.5,
        #vjust = 0,
        #hjust = 0,
        size = 2,
        color = "grey36",
        fontface = "italic",
        parse = FALSE
      )
  }
  p
}
