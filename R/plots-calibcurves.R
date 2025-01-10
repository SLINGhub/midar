' Plot Calibration Curves
#'
#' Plots calibration curves for each measured feature.
#'
#' @param data A `MidarExperiment` object .
#' @param variable The variable name to plot on the y-axis, usually a measure of intensity. Defaults to `"intensity"`
#' @param overwrite_metadata A logical value (`TRUE` or `FALSE`). If `TRUE`, the function will ignore any fit method and weighting settings defined in the metadata and use the provided `fit_method` and `fit_weighting` values for all analytes.
#' @param fit_method A character string specifying the default regression fit method to use for the calibration curve. Must be one of `"linear"` or `"quadratic"`. This method will be applied if no specific fit method is defined for a feature in the metadata.
#' @param fit_weighting A character string specifying the default weighting method for the regression points in the calibration curve. Must be one of `"none"`, `"1/x"`, or `"1/x^2"`. This method will be applied if no specific weighting method is defined for a feature in the metadata.
#' @param log_axes Logical. If `TRUE`, the x and y axes are displayed in logarithmic scale. Defaults to `FALSE`.
#' @param filter_data Logical. Indicates whether to use quality control (QC) filtered data (`TRUE`) or the raw data (`FALSE`). Defaults to `FALSE`.
#' @param output_pdf Logical. If `TRUE`, saves the generated plots as a PDF file. Defaults to `FALSE`.
#' @param include_feature_filter A regex pattern to filter and include features that match the criteria. If omitted, all features are considered.
#' @param exclude_feature_filter A regex pattern to filter and exclude features that match the criteria. If omitted, no features are excluded.
#' @param max_regression_value The maximum x value (analyzed amount) for which the regression line is fitted. Defaults to `NA`, considering all data points.
#' @param path The file path for saving the PDF. Defaults to an empty string, meaning no file is saved.
#' @param rows_page The number of rows of plots per page for the PDF output. Defaults to 4.
#' @param cols_page The number of columns of plots per page for the PDF output. Defaults to 5.
#' @param specific_page An integer specifying a particular page to plot. If `NA`, all pages are plotted.
#' @param point_size Numeric value specifying the size of points in the plot. Defaults to 1.5.
#' @param line_width Numeric value specifying the width of regression lines. Defaults to 0.7.
#' @param scaling_factor A universal scaling factor for fonts, symbols, and lines. Defaults to 1.
#' @param page_orientation The orientation of the PDF paper, either `"LANDSCAPE"` or `"PORTRAIT"`. Defaults to `"LANDSCAPE"`.
#' @param font_base_size The base font size for text in the plots. Defaults to 7.
#' @param return_plots Logical. If `TRUE`, returns the plots as a list of `ggplot2` objects. Each item represents a page of plots. Defaults to `FALSE`.
#' @param show_progress Logical. Show a progress bar during plot creation if set to `TRUE`. Defaults to `TRUE`.
#' @return If `return_plots` is `TRUE`, a list of `ggplot2` objects is returned. Otherwise, the function may return a plot output and/or save a PDF, depending on the `output_pdf` parameter.

#' @export
plot_calibrationcurves <- function(data = NULL,
                                   variable = "intensity",
                                   overwrite_metadata = FALSE,
                                   fit_method = c("linear", "quadratic"),
                                   fit_weighting = c(NA, "none", "1/x", "1/x^2"),
                                   log_axes = FALSE,
                                   filter_data = FALSE,
                                   output_pdf = FALSE,
                                   include_feature_filter = "",
                                   exclude_feature_filter = "",
                                   max_regression_value = NA,
                                   path = "",
                                   rows_page = 4,
                                   cols_page = 5,
                                   specific_page = NA,
                                   point_size = 1.5,
                                   line_width = 0.7,
                                   scaling_factor = 1,
                                   page_orientation = "LANDSCAPE",
                                   font_base_size = 7,
                                   show_progress = TRUE,
                                   return_plots = FALSE) {
  check_data(data)
  variable_strip <- str_remove(variable, "feature_")
  rlang::arg_match(variable_strip, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm"))
  variable <- stringr::str_c("feature_", variable_strip)
  variable_sym = rlang::sym(variable)
  rlang::arg_match(fit_method, c(NA, "linear", "quadratic"))
  rlang::arg_match(fit_weighting, c(NA, "none", "1/x", "1/x^2"))
  #rlang::arg_match(overwrite_metadata, c("TRUE", "FALSE"))
  check_var_in_dataset(data@dataset, variable)



  if (output_pdf & path == "") cli::cli_abort("Please define parameter `path`")
  if (nrow(data@dataset) < 1) cli::cli_abort("No data available. Please import data and metadata first.")

  # Filter data if filter_data is TRUE
  if (filter_data) {
    dat_filt <- data@dataset_filtered |> dplyr::ungroup()
    if (!data@is_filtered) cli::cli_abort("Data has not been qc filtered, or has changed. Please run `filter_features_qc` first.")
  } else {
    dat_filt <- data@dataset |> dplyr::ungroup()
  }



  dat_filt <- dat_filt |> dplyr::semi_join(data@annot_qcconcentrations, by = c("sample_id", "feature_id"))

  if (all(!is.na(include_feature_filter)) & all(include_feature_filter != "")) {
    if (length(include_feature_filter) == 1) {
      dat_filt <- dat_filt |> dplyr::filter(stringr::str_detect(.data$feature_id, include_feature_filter))
    } else {
      dat_filt <- dat_filt |> dplyr::filter(.data$feature_id %in% include_feature_filter)
    }
  }

  if (all(!is.na(exclude_feature_filter)) & all(exclude_feature_filter != "")) {
    if (length(exclude_feature_filter) == 1) {
      dat_filt <- dat_filt |> dplyr::filter(!stringr::str_detect(.data$feature_id, exclude_feature_filter))
    } else {
      dat_filt <- dat_filt |> dplyr::filter(!.data$feature_id %in% exclude_feature_filter)
    }
  }


  d_calib <- dat_filt |>
    dplyr::select(tidyselect::any_of(
      c("analysis_id", "sample_id", "qc_type", "feature_id", variable)
    )) |>
    filter(str_detect(.data$qc_type, "CAL|[MLH]QC|^QC")) |>
    dplyr::right_join(data@annot_qcconcentrations, by = c("sample_id" = "sample_id", "feature_id" = "feature_id"))

  d_calib$curve_id = 1


  d_calib <-  d_calib |>
    left_join(data@annot_features |> select("feature_id", "curve_fit_method", "curve_fit_weighting"), by = "feature_id")


  d_calib <- d_calib |> mutate(curve_fit_method = ifelse(is.na(.data$curve_fit_method) | overwrite_metadata, fit_method, .data$curve_fit_method),
                               fit_weighting = ifelse(is.na(.data$curve_fit_weighting) | overwrite_metadata, fit_weighting, .data$curve_fit_weighting))




  # Verify if unit is the same for all data points/curves
  x_axis_unit <- unique(d_calib$concentration_unit)
  if (length(x_axis_unit) > 1)
    cli::cli_abort(cli::col_red("The `concentration_unit` (x axis) must be identical for selected curves and data points. Please change selection of curves or update calibration curve metadata."))


  max_regression_value <-
    ifelse(
      is.na(max_regression_value),
      max(d_calib$concentration, na.rm = TRUE),
      max_regression_value
    )

  if (output_pdf & !is.na(path)) {
    path <- ifelse(stringr::str_detect(path, ".pdf"), path, paste0(path, ".pdf"))
    if (page_orientation == "LANDSCAPE") {
      pdf(file = path, onefile = T, paper = "A4r", useDingbats = FALSE, width = 28 / 2.54, height = 20 / 2.54)
    } else {
      pdf(file = path, onefile = T, paper = "A4", useDingbats = FALSE, height = 28 / 2.54, width = 20 / 2.54)
    }
  }

  if (is.na(specific_page)) {
    page_range <- 1:ceiling(dplyr::n_distinct(d_calib$feature_id) / (cols_page * rows_page))
  } else {
    page_range <- specific_page
  }

  if(output_pdf) action_text = "Saving plots to pdf" else action_text = "Generating plots"
  cat(cli::col_green(glue::glue("{action_text} ({max(page_range)} {ifelse(max(page_range) > 1, 'pages', 'page')}){ifelse(show_progress, ':', '...')}")))
  if(show_progress) pb <- txtProgressBar( min = 0, max = max(page_range), width = 50, style = 3)

  p_list <- list()  # p_list <- vector("list", length(page_range))

  for (i in page_range) {


    p <-  plot_calibcurves_page(
      data = data,
      dataset = d_calib,
      output_pdf = output_pdf,
      response_variable = variable,
      max_regression_value = max_regression_value,
      path = path,
      rows_page = rows_page,
      cols_page = cols_page,
      specific_page = i,
      point_size = point_size,
      line_width = line_width,
      scaling_factor = scaling_factor,
      font_base_size = font_base_size,
      x_axis_title = x_axis_unit,
      fit_method = fit_method,
      fit_weighting = fit_weighting,
      log_axes = log_axes,
      overwrite_metadata = overwrite_metadata
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






# Define function to plot 1 page
plot_calibcurves_page <- function(data,
                                  dataset,
                                     output_pdf,
                                     response_variable,
                                     max_regression_value,
                                     path,
                                     rows_page,
                                     cols_page,
                                     specific_page,
                                     point_size,
                                     line_width,
                                     scaling_factor,
                                     font_base_size,
                                     x_axis_title,
                                  fit_method,
                                  fit_weighting,
                                  log_axes,
                                  overwrite_metadata
                                  ) {
  plot_var <- rlang::sym(response_variable)
  dataset$curve_id <- as.character(dataset$curve_id)

  # subset the dataset with only the rows used for plotting the facets of the selected page
  n_cmpd <- length(unique(dataset$analysis_id))
  row_start <- n_cmpd * cols_page * rows_page * (specific_page - 1) + 1
  row_end <- n_cmpd * cols_page * rows_page * specific_page

  dat_subset <- dataset |>
    dplyr::arrange(.data$feature_id, .data$curve_id) |>
    dplyr::slice(row_start:row_end)

  dat_subset <- dat_subset |> mutate(qc_type_cat = if_else(str_detect(.data$qc_type, "CAL"), "CAL", "QC"))

  data <- calc_calibration_results(data = data,
                                   overwrite_metadata = overwrite_metadata,
                                   fit_method = fit_method,
                                   fit_weighting = fit_weighting,
                                   include_fit_object = TRUE)

  d_calib_stats <- data@metrics_calibration



 get_predictions <- function(stats){

    fit <- stats$fit_cal_1[[1]]
    conc_orig <- model.matrix(fit)[,2]
    concs <- seq(min(conc_orig), max(conc_orig), length.out = 100)
    predictions <- predict(fit, newdata = data.frame(concentration = concs), interval = "confidence")
    prediction_data <- tibble(
      feature_id = stats$feature_id,
      curve_id = "1",
      concentration = concs,
      y_pred = predictions[, "fit"],
      lwr = predictions[, "lwr"],
      upr = predictions[, "upr"]
    )
    prediction_data
 }


 d_calib_stats_grp <- d_calib_stats |>
   dplyr::semi_join(dat_subset , by = c("feature_id")) |>
   dplyr::group_split(.data$feature_id) #TODO .data$curve_id

 d_pred <- map(d_calib_stats_grp, function(x) get_predictions(x)) |> bind_rows()

 d_pred <- d_pred |>
   dplyr::arrange(.data$feature_id, .data$curve_id)


    p <- ggplot(
        data = dat_subset,
        aes(
          x = .data$concentration,
          y = !!plot_var,
          color = .data$qc_type,
          shape = .data$qc_type
        )
      ) +
      ggplot2::geom_line(data = d_pred,
                aes(x = .data$concentration,
                    y = .data$y_pred,
                    shape = "CAL"),
                color = "#4575b4",
                fill = "#4575b4",
                linewidth = 1, inherit.aes = TRUE) +
      # # Confidence interval
      ggplot2::geom_ribbon(data = d_pred,
                  aes(x = .data$concentration, ymin = .data$lwr, ymax = .data$upr, shape = "CAL"),
                  fill = "#91bfdb", alpha = 0.25, inherit.aes = FALSE) +

        #color = ifelse(after_stat(r.squared) < 0.80, "red", "darkgreen")), size = 1.4) +
        scale_color_manual(values = c("CAL" ="#254f6e", "HQC" = "#f27507", "LQC" = "#f27507")) +
        scale_fill_manual(values = c("CAL" ="white", "HQC" = "#f5c969", "LQC" = "#f5c969")) +
        scale_shape_manual(values = c("CAL" = 21, "HQC" = 24, "LQC" = 25)) +
        scale_y_continuous(limits = c(0, NA)) +
        scale_x_continuous(limits = c(0, NA), breaks = scales::breaks_extended(6)) +
        ggh4x::facet_wrap2(
          vars(.data$feature_id),
          scales = "free",
          nrow = rows_page,
          ncol = cols_page,
          trim_blank = FALSE
        ) +
        geom_point(aes(fill = .data$qc_type), size = point_size) +

        labs( x= x_axis_title, y = stringr::str_remove(response_variable, "feature\\_")) +
        theme_light(base_size = font_base_size) +
        theme(
          strip.text = element_text(size = font_base_size * scaling_factor, face = "bold"),
          strip.background = element_rect(size = 0.0001, fill = "#496875"),
          panel.grid.major = element_line(color = "grey70", size = 0.2, linetype = "dotted"),  # Light and dotted major gridlines
          panel.grid.minor = element_line(color = "grey90", size = 0.1, linetype = "dotted")   # Lighter minor gridlines

           )


    if(log_axes){
      p <- p +
        scale_x_continuous(trans = "log10") +
        scale_y_continuous(trans = "log10")
    }

    p <- p +
      ggplot2::geom_text(data = d_calib_stats,
                aes(x = 0, y = Inf, label = paste0(stringr::str_to_title(.data$fit_model), "\nR^2 = ", round(.data$r2_cal_1, 4))),
                inherit.aes = FALSE,
                hjust = 0, vjust = 1.3, size = 2, color = "grey36", fontface = "italic", parse = TRUE)


  p
}


