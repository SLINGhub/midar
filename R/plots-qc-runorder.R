#' RunSequence plot
#'
#' @param data MidarExperiment object
#' @param qc_type_subet Select QC types to be show, NA will select all available QC/Sample types
#' @param show_qc_dataset_only Show only available QC types
#' @param show_batches Show batches
#' @param batches_as_shades Show batches as shades
#' @param batch_line_color batch separator color
#' @param batch_shading_color batch shade color
#' @param scale_factor Overall plot scale factor
#' @param show_outlier Should samples defined as outlier as separate row
#' @param segment_width Linewidth of the segments
#' @param base_size base font size
#' @param factor_a First factor referring to a field in the study sample metadata used to show randomization
#' @param factor_b First factor referring to a field in the study sample metadata used to show randomization
#' @return ggplot object
#' @export

plot_runsequence <- function(data,
                             qc_type_subet,
                             show_qc_dataset_only = FALSE,
                             show_batches = TRUE,
                             segment_width = 0.5,
                             batches_as_shades = FALSE,
                             batch_line_color = "darkred",
                             batch_shading_color = "grey60",
                             scale_factor = 1,
                             show_outlier = TRUE,
                             base_size = 12,
                             factor_a = NA,
                             factor_b = NA) {
  d_temp <- data$dataset %>%
    dplyr::select(.data$run_id, .data$batch_id, .data$analysis_id, .data$qc_type, .data$sample_id) %>%
    distinct()

  d_temp$qc_type <- factor(d_temp$qc_type, c("EQC", "SST", "MBLK", "SBLK", "UBLK", "PBLK", "RQC", "LTR", "NIST", "TQC", "BQC", "SPL"))

  d_temp$sample_category <- as.character(d_temp$qc_type)



  # if (factor_a != "" & factor_b !="" & nrow(data@annot_studysamples)> 0){
  #   fac_a <- sym(factor_a)
  #   fac_b <- sym(factor_b)
  #
  #   data@annot_studysamples <-  data@annot_studysamples %>%
  #     mutate(factor_a = !!fac_a,
  #            factor_b = !!fac_b)
  #
  #   d_temp <- bind_rows(d_temp, d_temp %>% filter(sample_category == "SPL") %>%
  #                         mutate(qc_type = ifelse(qc_type == "SPL", "SPL_B", qc_type)))
  #
  #
  #   d_temp <- d_temp %>% left_join(data$annot_studysamples %>% dplyr::select(sample_id, factor_a, factor_b), by=(c("sample_id"="sample_id")))
  #   d_temp <- d_temp %>% mutate(sample_category = ifelse(qc_type=="SPL", paste0("SPL_", factor_a), sample_category),
  #                               sample_category = ifelse(qc_type=="SPL_B", paste0("SPL_", factor_b), sample_category))
  #   d_temp <- d_temp %>% mutate(sample_category = fct_rev(as_factor(sample_category))) %>% droplevels()
  # }

  # d_temp <- d_temp %>% mutate(sample_category = forcats::fct_rev(forcats::as_factor(sample_category)))
  d_temp <- d_temp %>%
    filter(str_detect(as.character(.data$sample_category), qc_type_subet)) %>%
    {
      if (show_qc_dataset_only | qc_type_subet != "") droplevels(.) else .
    }

  d_batch_info <- data@annot_batches

  # round to next 10 and divide by number of breaks will be showb
  scale_dataset_size <- 10^ceiling(log10(max(d_temp$run_id))) / 10
  p <- ggplot(d_temp, aes(x = .data$run_id, y = rev(.data$qc_type), color = .data$sample_category))
  if (show_batches) {
    if (batches_as_shades) {
      d_batch_2nd <- d_batch_info %>%
        dplyr::slice(-1) %>%
        filter(.data$batch_no %% 2 != 1)
      p <- p + geom_rect(
        data = d_batch_2nd, aes(xmin = .data$id_batch_start - 0.5, xmax = .data$id_batch_end + 0.5, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = batch_shading_color, alpha = 0.1, color = NA, linetype = "solid", size = 0.5
      )
    }
  }
  p <- p + geom_segment(
    aes(
      x = .data$run_id,
      xend = .data$run_id,
      y = as.integer(.data$qc_type) - 0.4,
      yend = as.integer(.data$qc_type) + 0.4
    ),
    size = segment_width
  ) +
    labs(
      x = "Analysis order",
      y = "Sample Type"
    ) +
    scale_x_continuous(breaks = seq(0, max(d_temp$run_id), scale_dataset_size)) +
    scale_y_discrete(
      limits = c(levels(d_temp$qc_type)),
      expand = expansion(0.05, 0.05)
    ) +
    # scale_color_manual(values =  data$sample_category) +
    theme_light(base_size = base_size) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey80", linetype = "solid", size = .5),
      panel.grid.minor.x = element_line(colour = "grey90", linetype = "dotted", size = 0.25),
      panel.border = element_rect(size = 2),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text.x = element_text(face = "plain", size = base_size),
      axis.text.y = element_text(face = "bold", size = base_size),
      legend.position = "none"
    )


  if (show_batches) {
    if (!batches_as_shades) {
      p <- p + geom_vline(data = d_batch_info %>% dplyr::slice(-1), aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, size = 0.5)
    }
  }
  # if (factor_a != "" & factor_b !=""){

  #  p <- p + scale_color_brewer(palette="Dark2")
  # } else {
  p <- p + scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col)
  # }
  return(p)
}




#' RunScatter plot
#'
#' @param data MidarExperiment object
#' @param use_filt_data Use QC-filtered data
#' @param plot_var Variable to plot
#' @param qc_types QC type to plot. When qc_types us NA or NULL, all available QC types are plotted.
#' @param feature_incl_filt Filter features names matching the criteria (regex). When empty, `NA` or `NULL` all available features are included.
#' @param feature_excl_filt Exclude features names matching the criteria (regex).  When empty, `NA` or `NULL` all available features are included.
#' @param log_scale Use log10 scale for the y axis
#' @param cap_outliers Cap outliers to a specific range in the plot
#' @param cap_spl_iqr_factor Multiplicator for the upper Tukey's IQR fence use for outlier capping of samples
#' @param cap_qc_iqr_factor Multiplicator for the upper Tukey's IQR fence use for outlier capping of QCs
#' @param cap_top_n_values Cap top n values
#' @param qc_type_fit QC TYPE used for loess fit
#' @param show_driftcorrection Show drift correction
#' @param show_trend_samples Fit trend line for study samples, if show_driftcorrection is TRUE
#' @param trend_samples_fun Function used for drift correction. Default 'loess'
#' @param trend_samples_col Color of drift line
#' @param after_correction Show before/after correction
#' @param plot_other_qc Plot all QCS
#' @param show_batches Show batches
#' @param batches_as_shades Show batches as shades
#' @param batch_line_color batch separator color
#' @param batch_shading_color batch shade color
#' @param save_pdf save as PDF
#' @param path file name of PDF
#' @param cols_page columns per page
#' @param rows_page rows per page
#' @param annot_scale scale factor of text elements
#' @param paper_orientation Landscape/Portrait
#' @param point_transparency Alpha of points
#' @param point_size point size
#' @param page_no Show page number
#' @param y_label_text Overwrite y label with this text
#' @param silent Verbose or silent
#' @param point_stroke_width point stroke width
#' @param base_size base font size of plot
#' @param return_plot_list return list with plots
#' @param show_gridlines show x and y major gridlines
#' @return A list of ggplot2 plots or NULL
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom MASS rlm
#' @importFrom scales percent_format
#' @importFrom ggpmisc stat_poly_line stat_poly_eq
#' @importFrom grDevices dev.off pdf
#' @importFrom stats prcomp
#' @importFrom utils head
#' @export

plot_runscatter <- function(data,
                            plot_var = c("feature_intensity", "feature_norm_intensity", "feature_conc"),
                            use_filt_data = FALSE,
                            qc_types = NA,
                            feature_incl_filt = "",
                            feature_excl_filt = "",
                            log_scale = FALSE,
                            cap_outliers = FALSE,
                            cap_qc_iqr_factor = 3,
                            cap_spl_iqr_factor = 3,
                            cap_top_n_values = NA,
                            show_driftcorrection = FALSE,
                            qc_type_fit = "BQC",
                            show_trend_samples = FALSE,
                            trend_samples_fun = "loess",
                            trend_samples_col = "",
                            after_correction = FALSE,
                            plot_other_qc = TRUE,
                            show_batches = TRUE,
                            batches_as_shades = FALSE,
                            batch_line_color = "#9dbecf",
                            batch_shading_color = "grey90",
                            save_pdf = FALSE,
                            path = "",
                            cols_page = 3,
                            rows_page = 3,
                            annot_scale = 1.0,
                            paper_orientation = "LANDSCAPE",
                            point_transparency = 1,
                            point_size = 2,
                            point_stroke_width = 1,
                            page_no = NA,
                            y_label_text = NA,
                            silent = TRUE,
                            return_plot_list = FALSE,
                            base_size = 12,
                            show_gridlines = FALSE) {
  if (nrow(data@dataset) < 1) cli::cli_abort("No data available. Please import data and metadata first.")

  if (use_filt_data) {
    dat_filt <- data@dataset_filtered %>% dplyr::ungroup()
    if (nrow(dat_filt) < 1) cli::cli_abort("Data has not been qc filtered. Please apply `apply_qc_filter` first.")
  } else {
    dat_filt <- data@dataset %>% dplyr::ungroup()
  }

  plot_var <- rlang::arg_match(plot_var)
  plot_var_s <- rlang::sym(plot_var)
  y_label <- dplyr::if_else(cap_outliers, paste0(ifelse(is.na(y_label_text), plot_var, y_label_text), " (capped min(", cap_spl_iqr_factor, "x IQR+Q3[SPL]) ,", cap_qc_iqr_factor, "x IQR+Q3[QC]"), stringr::str_remove(plot_var, "feature\\_"))




  # Subset data ----

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

  if (all(!is.na(qc_types)) & all(qc_types != "")) {
    if (length(qc_types) == 1) {
      dat_filt <- dat_filt |> dplyr::filter(stringr::str_detect(.data$qc_type, qc_types))
    } else {
      dat_filt <- dat_filt |> dplyr::filter(.data$qc_type %in% qc_types)
    }
  }

  if (nrow(dat_filt) < 1) cli::cli_abort("None of the feature names meet the filter criteria. Please check feature_filter_include and feature_filter_exclude parameters.")

  # Re-order qc_type levels to define plot layers, i.e.  QCs are plotted over StudySamples
  dat_filt$qc_type <- factor(as.character(dat_filt$qc_type), pkg.env$qc_type_annotation$qc_type_levels)

  # Cap upper outliers  ----

  dat_filt <- dat_filt %>% dplyr::mutate(value = !!plot_var_s)


  dat_filt <- dat_filt %>%
    dplyr::group_by(.data$feature_id) %>%
    dplyr::mutate(
      # value_max_spl = mean(.data$value[.data$qc_type=="SPL"], na.rm=T) + cap_SPL_SD * sd(.data$value[.data$qc_type=="SPL"], na.rm=T),
      # value_max_qc = mean(.data$value[.data$qc_type==qc_type_fit], na.rm=T) + cap_QC_SD * sd(.data$value[.data$qc_type==qc_type_fit]), na.rm=T,
      value_max_spl = quantile(.data$value[.data$qc_type == "SPL"], 0.75, names = FALSE, na.rm = TRUE) + cap_spl_iqr_factor * IQR(.data$value[.data$qc_type == "SPL"], na.rm = TRUE),
      value_max_qc = quantile(.data$value[.data$qc_type == qc_type_fit], 0.75, names = FALSE, na.rm = TRUE) + cap_qc_iqr_factor * IQR(.data$value[.data$qc_type == qc_type_fit], na.rm = TRUE),
      value_max = pmax(.data$value_max_spl, .data$value_max_qc, na.rm = TRUE),
      value_mod = dplyr::if_else(cap_outliers & .data$value > .data$value_max,
        if (!all(is.na(.data$value))) {
          max(.data$value[.data$value < .data$value_max], na.rm = TRUE)
        } else {
          NA_real_
        },
        .data$value
      )
    ) %>%
    dplyr::ungroup()

  dat_filt <- dat_filt %>%
    dplyr::group_by(.data$feature_id) %>%
    dplyr::arrange(.data$value) %>%
    mutate(
      value = ifelse(dplyr::row_number() < cap_top_n_values, .data$value[cap_top_n_values], .data$value)
    ) |>
    dplyr::arrange(.data$feature_id, .data$run_id) %>%
    dplyr::ungroup()


  if (save_pdf & !is.na(path)) {
    path <- ifelse(stringr::str_detect(path, ".pdf"), path, paste0(path, ".pdf"))
    if (paper_orientation == "LANDSCAPE") {
      pdf(file = path, onefile = T, paper = "A4r", useDingbats = FALSE, width = 28 / 2.54, height = 20 / 2.54)
    } else {
      pdf(file = path, onefile = T, paper = "A4", useDingbats = FALSE, height = 28 / 2.54, width = 20 / 2.54)
    }
  }

  if (is.na(page_no)) {
    page_range <- 1:ceiling(dplyr::n_distinct(dat_filt$feature_id) / (cols_page * rows_page))
  } else {
    page_range <- page_no
  }

  # TODO if(!silent) print(paste0("Plotting ", max(page_range), " pages..."))

  # p_list <- vector("list", length(page_range))
  p_list <- list()
  for (i in page_range) {
    if (!silent) print(paste0("page ", i))
    p <- runscatter_one_page(
      dat_filt = dat_filt, data = data, d_batches = data@annot_batches, cols_page = cols_page, rows_page = rows_page, show_driftcorrection = show_driftcorrection,
      show_trend_samples, trend_samples_fun, trend_samples_col, after_correction = after_correction, qc_type_fit = qc_type_fit, save_pdf = save_pdf, page_no = i,
      point_size = point_size, cap_outliers = cap_outliers, point_transparency = point_transparency, annot_scale = annot_scale,
      show_batches = show_batches, batches_as_shades = batches_as_shades, batch_line_color = batch_line_color, plot_other_qc,
      batch_shading_color = batch_shading_color, y_label = y_label, base_size = base_size, point_stroke_width = point_stroke_width, show_grid = show_gridlines, log_scale = log_scale
    )
    plot(p)
    p_list[[i]] <- p
  }
  on.exit(if (save_pdf) {
    dev.off()
  })
  if (return_plot_list) {
    return(p_list)
  }
}
#' @importFrom ggplot2 Stat
runscatter_one_page <- function(dat_filt, data, d_batches, cols_page, rows_page, page_no,
                                show_driftcorrection, after_correction = FALSE, qc_type_fit, cap_outliers,
                                show_batches, batches_as_shades, batch_line_color, batch_shading_color, show_trend_samples, trend_samples_fun, trend_samples_col, plot_other_qc,
                                save_pdf, annot_scale, point_transparency, point_size = point_size, y_label, base_size, point_stroke_width, show_grid, log_scale) {
  point_size <- ifelse(missing(point_size), 2, point_size)
  point_stroke_width <- dplyr::if_else(save_pdf, .3, .2 * (1 + annot_scale / 5))


  # subset the dataset with only the rows used for plotting the facets of the selected page
  n_cmpd <- length(unique(dat_filt$analysis_id))
  row_start <- n_cmpd * cols_page * rows_page * (page_no - 1) + 1
  row_end <- n_cmpd * cols_page * rows_page * page_no


  dat_subset <- dat_filt %>%
    dplyr::arrange(.data$feature_id, .data$run_id) %>%
    dplyr::slice(row_start:row_end)

  dat_subset$qc_type <- forcats::fct_relevel(dat_subset$qc_type, c("SPL", "UBLK", "SBLK", "TQC", "BQC", "RQC", "LTR", "NIST", "PBLK"))
  dat_subset <- dat_subset %>%
    dplyr::arrange(.data$qc_type)




  # https://stackoverflow.com/questions/46327431/facet-wrap-add-geom-hline
  if (nrow(dat_subset) > 0) {
    dMax <- dat_subset %>%
      dplyr::group_by(.data$feature_id) %>%
      dplyr::summarise(
        y_max =
          if (!all(is.na(.data$value_mod))) {
            max(.data$value_mod, na.rm = TRUE) * 1.0
          } else {
            NA_real_
          }
      )
  }


  d_batch_data <- d_batches %>% dplyr::slice(rep(1:dplyr::n(), each = nrow(dMax)))
  d_batch_data$feature_id <- rep(dMax$feature_id, times = nrow(d_batches))
  d_batch_data <- d_batch_data %>% dplyr::left_join(dMax, by = c("feature_id"))
  p <- ggplot2::ggplot(dat_subset, ggplot2::aes_string(x = "run_id"))

  # browser()
  if (show_batches) {
    if (!batches_as_shades) {
      d_batches_temp <- d_batch_data |> filter(.data$id_batch_start != 1)
      p <- p + ggplot2::geom_vline(data = d_batches_temp, ggplot2::aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, linetype = "solid", size = .5, na.rm = TRUE)
    } else {
      d_batches_temp <- d_batch_data %>% dplyr::filter(.data$batch_no %% 2 != 1)
      p <- p + ggplot2::geom_rect(
        data = d_batches_temp, ggplot2::aes(xmin = .data$id_batch_start - 0.5, xmax = .data$id_batch_end + 0.5, ymin = 0, ymax = .data$y_max, label = .data$batch_id),
        inherit.aes = FALSE, fill = batch_shading_color, color = NA, alpha = 0.5, linetype = "solid", size = 0.3, na.rm = TRUE
      )
    }
  }

  if (cap_outliers) {
    p <- p + ggplot2::geom_hline(data = dMax, ggplot2::aes(yintercept = .data$y_max), color = "#fa9b9b", size = 3, alpha = .4)
  }

  p <- p +
    ggplot2::geom_point(aes_string(x = "run_id", y = "value_mod", color = "qc_type", fill = "qc_type", shape = "qc_type", group = "batch_id"), size = point_size, alpha = point_transparency, stroke = point_stroke_width, na.rm = TRUE)



  if (after_correction & show_driftcorrection) {
    p <- p +
      ggplot2::geom_line(aes_string(x = "run_id", y = "CURVE_Y_PREDICTED", group = "batch_id"), color = pkg.env$qc_type_annotation$qc_type_fillcol[qc_type_fit], size = .5, na.rm = TRUE)
  }

  p <- p +
    ggh4x::facet_wrap2(ggplot2::vars(.data$feature_id), scales = "free_y", ncol = cols_page, nrow = rows_page, trim_blank = FALSE) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col, drop = TRUE) +
    ggplot2::scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_fillcol, drop = TRUE) +
    ggplot2::scale_shape_manual(values = pkg.env$qc_type_annotation$qc_type_shape, drop = TRUE)



  if (show_driftcorrection) {
    if (after_correction) {
      p <- p +
        ggplot2::geom_smooth(
          data = filter(dat_subset, .data$qc_type == qc_type_fit), ggplot2::aes_string(x = "run_id", y = "value", group = "batch_id"), se = TRUE, na.rm = TRUE,
          colour = pkg.env$qc_type_annotation$qc_type_fillcol[qc_type_fit], fill = pkg.env$qc_type_annotation$qc_type_fillcol[qc_type_fit],
          method = MASS::rlm, alpha = 0.35, size = 0.8
        )

      if (show_trend_samples) {
        p <- p + ggplot2::geom_smooth(
          data = filter(dat_subset, .data$qc_type == "SPL"), ggplot2::aes_string(x = "run_id", y = "value", group = "batch_id"), colour = trend_samples_col, fill = trend_samples_col, linetype = "dashed",
          method = trend_samples_fun, se = FALSE, alpha = 0.2, size = .8, na.rm = TRUE
        )
      }

      if (plot_other_qc) {
        other_qc <- dplyr::if_else(qc_type_fit == "BQC", "TQC", "BQC")
        other_qc_col <- pkg.env$qc_type_annotation$qc_type_fillcol[other_qc]
        p <- p +
          ggplot2::geom_smooth(
            data = dplyr::filter(dat_subset, .data$qc_type == other_qc), ggplot2::aes_string(x = "run_id", y = "value", group = "batch_id"), colour = other_qc_col, fill = other_qc_col,
            method = trend_samples_fun, se = TRUE, alpha = 0.2, size = .4, na.rm = FALSE
          )
      }
    } else {
      p <- p +
        geom_smooth(
          data = dplyr::filter(dat_subset, .data$qc_type == "SPL"), ggplot2::aes_string(x = "run_id", y = "value", group = "batch_id"), colour = trend_samples_col, fill = trend_samples_col,
          method = trend_samples_fun, se = TRUE, alpha = 0.2, size = .4, na.rm = FALSE
        )

      if (plot_other_qc) {
        other_qc <- dplyr::if_else(.data$qc_type_fit == "BQC", "TQC", "BQC")
        other_qc_col <- pkg.env$qc_type_annotation$qc_type_fillcol[other_qc]
        p <- p +
          ggplot2::geom_smooth(
            data = dplyr::filter(dat_subset, .data$qc_type == other_qc), ggplot2::aes_string(x = "run_id", y = "value", group = "batch_id"), colour = other_qc_col, fill = other_qc_col,
            method = trend_samples_fun, se = TRUE, alpha = 0.2, size = .4, na.rm = FALSE
          )
      }
    }
  }

  p <- p +
    # aes(ymin=0) +
    ggplot2::xlab("Analysis order") +
    ggplot2::ylab(label = y_label) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.02, 0.03))) +
    # expand_limits(y = 0) +
    ggplot2::theme_bw(base_size = base_size) +

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
  }

  if (log_scale) {
    p <- p + scale_y_log10()
  }


  return(p)
}




#' Sample sequence Boxplots
#'
#' @param data MidarExperiment object
#' @param relative_log_abundances Plot as RLA (Relative Log Abundnance Plot)
#' @param plot_var Variable to plot
#' @param use_qc_filtered_data Use QC-filtered data
#' @param min_feature_intensity Exclude features with median signal below this value
#' @param qc_types QC types to be plotted
#' @param ignore_outliers Exclude outlier values based on 3x IQR fences
#' @param feature_filter Filter text to select specific features (regex string)
#' @param show_batches Show batches
#' @param batches_as_shades Show batches as shades
#' @param batch_line_color batch separator color
#' @param batch_shading_color batch shade color
#' @param base_size base font size for plots (default is 8)
#' @return ggplot object
#'
#' @export
#'

plot_runboxplots <- function(data,
                             relative_log_abundances,
                             plot_var,
                             use_qc_filtered_data,
                             min_feature_intensity,
                             qc_types = c("BQC|TQC|SPL|NIST|LTR"),
                             ignore_outliers = TRUE,
                             feature_filter = "",
                             show_batches = TRUE,
                             batches_as_shades = FALSE,
                             batch_line_color = "red",
                             batch_shading_color = "grey70",
                             base_size = 8) {
  plot_var_sym <- sym(plot_var)





  if (!use_qc_filtered_data) {
    d_temp <- data@dataset
  } else {
    d_temp <- data@dataset_filtered
  }

  d_temp <- d_temp %>%
    dplyr::select(.data$analysis_id, .data$run_id, .data$qc_type, .data$batch_id, .data$feature_id, .data$feature_intensity, .data$feature_norm_intensity, .data$feature_conc) %>%
    filter(.data$feature_intensity > min_feature_intensity) %>%
    filter(str_detect(.data$qc_type, qc_types)) %>%
    droplevels()

  if (relative_log_abundances) {
    d_temp <- d_temp %>%
      group_by(.data$feature_id) %>%
      mutate(val = !!plot_var_sym) %>%
      mutate(val = .data$val / mean(.data$val[.data$qc_type == "BQC" | .data$qc_type == "TQC" | .data$qc_type == "SPL"], na.rm = TRUE))
  } else {
    d_temp <- d_temp %>% mutate(val = !!plot_var_sym)
  }

  breaks <- data$dataset %>%
    dplyr::select(.data$run_id) %>%
    distinct() %>%
    mutate(ticks_to_plot = .data$run_id %% 10 == 0) %>%
    pull(.data$run_id)


  # d_temp$run_id <- as_factor(d_temp$run_id)

  p <- ggplot(d_temp, aes(x = .data$run_id, y = log2(.data$val), group = .data$run_id))

  if (show_batches) {
    if (!batches_as_shades) {
      p <- p + geom_vline(data = data@annot_batches %>% slice(-1), aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, linetype = "solid", size = 1)
    } else {
      d_batch_2nd <- data$batch_info %>%
        slice(-1) %>%
        filter(.data$batch_id %% 2 != 1)
      p <- p + geom_rect(
        data = d_batch_2nd, aes(xmin = .data$id_batch_start - 0.5, xmax = .data$id_batch_end + 0.5, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE, fill = batch_shading_color, color = NA, alpha = 0.1, linetype = "solid", size = 0.5, na.rm = TRUE
      )
    }
  }

  scale_dataset_size <- 2^ceiling(log2(max(data$dataset$run_id))) / 100

  # geom_point(size=3, color = "#0053a8",alpha=0.6) +
  # stat_summary(aes(x=lipidClass, y=BQC_normIntensity_CV),fun.data="plot.median", geom="errorbar", colour="#fc0000", width=0.8, size=2, inherit.aes=FALSE,na.rm = TRUE) +
  p <- p +
    geom_boxplot(aes(fill = .data$qc_type, color = .data$qc_type), notch = FALSE, outlier.colour = NA, linewidth = 0.2, na.rm = TRUE) +
    # scale_colour_gradient(low = "white", high = "#004489") +
    scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_col) +
    scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col) +
    # scale_y_continuous(limits = c(0, 1e6)) +
    # scale_x_continuous(breaks = seq(0, max(data$dataset$run_id)+1, scale_dataset_size)) +
    #  scale_x_discrete(breaks = breaks) +
    theme_bw(base_size = base_size) +
    ylab(bquote(bold(log[2] ~ .(plot_var)))) +
    xlab("Analysis order") +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_line(colour = "#bdbdbd", linetype = "dotted", size = .5),
      panel.grid.minor.x = element_line(colour = "#bdbdbd", linetype = "dotted", size = .25),
      axis.text.y = element_text(size = base_size),
      axis.text.x = element_text(size = base_size),
      axis.title = element_text(size = base_size * 1, face = "bold"),
      panel.border = element_rect(linewidth = 1, color = "grey20")
    )

  if (ignore_outliers) {
    tails <- get_tails(log2(d_temp$val))
    p <- p + scale_y_continuous(limits = c(tails[1] * 2, tails[2] * 2))
  }

  if (relative_log_abundances) {
    p <- p + geom_hline(yintercept = 0, colour = "#666666", linetype = "dashed", size = 0.8) +
      ylab(bquote(bold("Rel. " ~ log[2] ~ .(plot_var))))
  }


  return(p)
}
