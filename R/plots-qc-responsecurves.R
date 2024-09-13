# Define function to plot 1 page
plot_responsecurves_page <- function(dataset,
                                     output_pdf,
                                     response_variable,
                                     regr_max_percent,
                                     path,
                                     rows_page,
                                     columns_page,
                                     point_size,
                                     line_width,
                                     text_scale_factor,
                                     base_size) {
  plot_var <- rlang::sym(response_variable)

  dataset$rqc_series_id <- as.character(dataset$rqc_series_id)

  ggplot2::ggplot(
    data = dataset,
    ggplot2::aes(
      x = .data$relative_sample_amount,
      y = !!plot_var,
      color = .data$rqc_series_id
    )
  ) +
    ggpmisc::stat_poly_line(
      data = subset(dataset, dataset$relative_sample_amount <= (regr_max_percent / 100)),
      ggplot2::aes(
        x = .data$relative_sample_amount,
        y = !!plot_var,
        color = .data$rqc_series_id
      ),
      se = FALSE, na.rm = TRUE, size = line_width, inherit.aes = FALSE
    ) +
    ggpmisc::stat_poly_eq(
      ggplot2::aes(group = .data$rqc_series_id, label = ggplot2::after_stat(.data$rr.label)),
      size = 2 * text_scale_factor, rr.digits = 3, vstep = .1
    ) +
    # color = ifelse(after_stat(r.squared) < 0.80, "red", "darkgreen")), size = 1.4) +
    ggplot2::scale_color_manual(values = c("#4575b4", "#91bfdb", "#fc8d59", "#d73027")) +
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::scale_x_continuous(
      limits = c(0, NA),
      breaks = c(0, 0.5, 1, 1.5, 2, 4),
      labels = scales::percent_format(accuracy = NULL)
    ) +
    ggh4x::facet_wrap2(
      ggplot2::vars(.data$feature_id),
      scales = "free",
      nrow = rows_page,
      ncol = columns_page,
      trim_blank = FALSE
    ) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::xlab("Sample Amount (Relative to BQC/TQC)") +
    ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(
      strip.text = ggplot2::element_text(size = 9 * text_scale_factor, face = "bold"),
      strip.background = ggplot2::element_rect(size = 0.0001, fill = "#8C8C8C")
    )
}


#' Response curves plot
#' @param data MidarExperiment object
#' @param use_filt_data Use QC-filtered data
#' @param output_pdf Save as PDF
#' @param response_variable Variable to plot
#' @param feature_incl_filt Filter features names matching the criteria (regex). When empty, `NA` or `NULL` all available features are included.
#' @param feature_excl_filt Exclude features names matching the criteria (regex).  When empty, `NA` or `NULL` all available features are included.
#' @param regr_max_percent Max relative sample amount to use in regressionb
#' @param path file name of pdf file
#' @param rows_page rows per page
#' @param columns_page columns per page
#' @param point_size point size
#' @param line_width regression line width
#' @param text_scale_factor text scale factor
#' @param return_plot_list return plot as list
#' @param base_size base font size
#' @return A list of ggplot2 objects
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @export
plot_responsecurves <- function(data,
                                use_filt_data,
                                output_pdf,
                                response_variable = "feature_intensity",
                                feature_incl_filt = "",
                                feature_excl_filt = "",
                                regr_max_percent = NA,
                                path = "",
                                rows_page = 4,
                                columns_page = 5,
                                point_size = 2,
                                line_width = 1,
                                text_scale_factor = 1,
                                return_plot_list = FALSE, base_size = 7) {
  if (output_pdf & path == "") cli::cli_abort("Please define parameter `path`")

  rows_page <- rows_page
  columns_page <- columns_page


  if (nrow(data@dataset) < 1) cli::cli_abort("No data available. Please import data and metadata first.")

  if (use_filt_data) {
    dat_filt <- data@dataset_filtered %>% dplyr::ungroup()
    if (nrow(dat_filt) < 1) cli::cli_abort("Data has not been qc filtered. Please apply `apply_qc_filter` first.")
  } else {
    dat_filt <- data@dataset %>% dplyr::ungroup()
  }

  dat_filt <- dat_filt |> dplyr::semi_join(data@annot_responsecurves, by = c("analysis_id"))

  if (all(!is.na(feature_incl_filt)) & all(feature_incl_filt != "")) {
    if (length(feature_incl_filt) == 1) {
      dat_filt <- dat_filt |> dplyr::filter(stringr::str_detect(.data$feature_id, feature_incl_filt))
    } else {
      dat_filt <- dat_filt |> dplyr::filter(.data$feature_id %in% feature_incl_filt)
    }
  }

  if (all(!is.na(feature_excl_filt)) & all(feature_excl_filt != "")) {
    if (length(feature_excl_filt) == 1) {
      dat_filt <- dat_filt |> dplyr::filter(!stringr::str_detect(.data$feature_id, feature_excl_filt))
    } else {
      dat_filt <- dat_filt |> dplyr::filter(!.data$feature_id %in% feature_excl_filt)
    }
  }

  d_rqc <- dat_filt |>
    dplyr::select(tidyselect::any_of(
      c("analysis_id", "feature_id", "feature_intensity", "feature_norm_intensity")
    )) |>
    dplyr::right_join(data@annot_responsecurves, by = c("analysis_id" = "analysis_id"))


  regr_max_percent <-
    ifelse(
      is.na(regr_max_percent),
      max(d_rqc$relative_sample_amount * 100),
      regr_max_percent
    )
  # browser()
  d_rqc_grp <- d_rqc %>%
    dplyr::left_join(tibble::tibble(feature_id = unique(d_rqc$feature_id)) |>
      mutate(grp = ceiling(
        row_number() / (rows_page * columns_page)
      )), by = "feature_id") %>%
    dplyr::group_by(.data$grp) %>%
    tidyr::nest() %>%
    dplyr::mutate(plt = purrr::map(
      data,
      function(x) {
        plot_responsecurves_page(
          dataset = x,
          output_pdf = output_pdf,
          response_variable = response_variable,
          regr_max_percent = regr_max_percent,
          path = path,
          rows_page = rows_page,
          columns_page = columns_page,
          point_size = point_size,
          line_width = line_width,
          text_scale_factor = text_scale_factor,
          base_size = base_size
        )
      }
    ))

  # Print pages
  # browser()
  if (!output_pdf) {
    if (!return_plot_list) {
    print(d_rqc_grp$plt)
    }
  } else {
    pdf(
      file = path,
      onefile = TRUE,
      paper = "A4r",
      width = 11
    )
    print(d_rqc_grp$plt)
    dev.off()
  }
  if (return_plot_list) {
    d_rqc_grp
  }
}
