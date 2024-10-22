# Define function to plot 1 page
qc_plot_responsecurves_page <- function(dataset,
                                     save_pdf,
                                     response_variable,
                                     regr_max_percent,
                                     path,
                                     rows_page,
                                     columns_page,
                                     page_no,
                                     point_size,
                                     line_width,
                                     text_scale_factor,
                                     base_size) {
  plot_var <- rlang::sym(response_variable)
  dataset$rqc_series_id <- as.character(dataset$rqc_series_id)

  # subset the dataset with only the rows used for plotting the facets of the selected page
  n_cmpd <- length(unique(dataset$analysis_id))
  row_start <- n_cmpd * columns_page * rows_page * (page_no - 1) + 1
  row_end <- n_cmpd * columns_page * rows_page * page_no


  dat_subset <- dataset |>
    dplyr::arrange(.data$feature_id, .data$rqc_series_id) |>
    dplyr::slice(row_start:row_end)

  p <- ggplot(
      data = dat_subset,
      aes(
        x = .data$relative_sample_amount,
        y = !!plot_var,
        color = .data$rqc_series_id
      )
    ) +
      ggpmisc::stat_poly_line(
        data = subset(dat_subset, dat_subset$relative_sample_amount <= (regr_max_percent / 100)),
        aes(
          x = .data$relative_sample_amount,
          y = !!plot_var,
          color = .data$rqc_series_id
        ),
        se = FALSE, na.rm = TRUE, size = line_width, inherit.aes = FALSE
      ) +
      ggpmisc::stat_poly_eq(
        aes(group = .data$rqc_series_id, label = ggplot2::after_stat(.data$rr.label)),
        size = 2 * text_scale_factor, rr.digits = 3, vstep = .1
      ) +
      # color = ifelse(after_stat(r.squared) < 0.80, "red", "darkgreen")), size = 1.4) +
      scale_color_manual(values = c("#4575b4", "#91bfdb", "#fc8d59", "#d73027")) +
      scale_y_continuous(limits = c(0, NA)) +
      scale_x_continuous(
        limits = c(0, NA),
        breaks = c(0, 0.5, 1, 1.5, 2, 4),
        labels = scales::percent_format(accuracy = NULL)
      ) +
      ggh4x::facet_wrap2(
        vars(.data$feature_id),
        scales = "free",
        nrow = rows_page,
        ncol = columns_page,
        trim_blank = FALSE
      ) +
      geom_point(size = point_size) +
      xlab("Sample Amount (Relative to BQC/TQC)") +
      theme_light(base_size = base_size) +
      theme(
        strip.text = element_text(size = 9 * text_scale_factor, face = "bold"),
        strip.background = element_rect(size = 0.0001, fill = "#8C8C8C")
      )
  p
}


#' Response curves plot
#' @param data MidarExperiment object
#' @param filter_data Use QC-filtered data
#' @param save_pdf Save as PDF
#' @param response_variable Variable to plot
#' @param feature_incl_filt Filter features names matching the criteria (regex). When empty, `NA` or `NULL` all available features are included.
#' @param feature_excl_filt Exclude features names matching the criteria (regex).  When empty, `NA` or `NULL` all available features are included.
#' @param regr_max_percent Max relative sample amount to use in regressionb
#' @param path file name of pdf file
#' @param rows_page rows per page
#' @param columns_page columns per page
#' @param page_no Specific page to plot. Default `NA`, meaning all pages are plotted
#' @param point_size point size
#' @param line_width regression line width
#' @param text_scale_factor text scale factor
#' @param paper_orientation Landscape/Portrait
#' @param return_plot_list return plot as list
#' @param base_size base font size
#' @param show_progress show progress bar
#' @return A list of ggplot2 objects
#' @export
qc_plot_responsecurves <- function(data,
                                response_variable = "feature_intensity",
                                filter_data,
                                feature_incl_filt = "",
                                feature_excl_filt = "",
                                save_pdf = FALSE,
                                path = "",
                                regr_max_percent = NA,
                                rows_page = 4,
                                columns_page = 5,
                                page_no = NA,
                                point_size = 2,
                                line_width = 1,
                                text_scale_factor = 1,
                                paper_orientation = "LANDSCAPE",
                                base_size = 7,
                                show_progress = TRUE,
                                return_plot_list = FALSE) {
  if (save_pdf & path == "") cli::cli_abort("Please define parameter `path`")



  if (nrow(data@dataset) < 1) cli::cli_abort("No data available. Please import data and metadata first.")

  if (filter_data) {
    dat_filt <- data@dataset_filtered |> dplyr::ungroup()
    if (nrow(dat_filt) < 1) cli::cli_abort("Data has not been qc filtered. Please apply `qc_set_feature_filters` first.")
  } else {
    dat_filt <- data@dataset |> dplyr::ungroup()
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

  if (save_pdf & !is.na(path)) {
    path <- ifelse(stringr::str_detect(path, ".pdf"), path, paste0(path, ".pdf"))
    if (paper_orientation == "LANDSCAPE") {
      pdf(file = path, onefile = T, paper = "A4r", useDingbats = FALSE, width = 28 / 2.54, height = 20 / 2.54)
    } else {
      pdf(file = path, onefile = T, paper = "A4", useDingbats = FALSE, height = 28 / 2.54, width = 20 / 2.54)
    }
  }

  if (is.na(page_no)) {
    page_range <- 1:ceiling(dplyr::n_distinct(d_rqc$feature_id) / (columns_page * rows_page))
  } else {
    page_range <- page_no
  }

  if(save_pdf) action_text = "Saving plots to pdf" else action_text = "Generating plots"
  message(cli::col_green(glue::glue("{action_text} ({max(page_range)} {ifelse(max(page_range) > 1, 'pages', 'page')}){ifelse(show_progress, ':', '...')}")))
  if(show_progress) pb <- txtProgressBar( min = 1, max = max(page_range), width = 50, style = 3)

  p_list <- list()  # p_list <- vector("list", length(page_range))
  for (i in page_range) {


    p <-  qc_plot_responsecurves_page(
        dataset = d_rqc,
        save_pdf = save_pdf,
        response_variable = response_variable,
        regr_max_percent = regr_max_percent,
        path = path,
        rows_page = rows_page,
        columns_page = columns_page,
        page_no = i,
        point_size = point_size,
        line_width = line_width,
        text_scale_factor = text_scale_factor,
        base_size = base_size
    )
    plot(p)
    dev.flush()
    flush.console()
    if(show_progress) setTxtProgressBar(pb, i)
    p_list[[i]] <- p
  }
  if (save_pdf) dev.off()
  cat(cli::col_green(" - done!"))
  if(show_progress) close(pb)

  flush.console()
  # on.exit(
  #   dev.off())

  if (return_plot_list) {
    return(p_list)
  }
}


