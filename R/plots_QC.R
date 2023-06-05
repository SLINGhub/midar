#' RunScatter plot
#'
#' @param data MidarExperiment object
#' @param y_var Variable to plot
#' @param transition_filter Filter features containing
#' @param filter_exclude Exclude or include transition_filter
#' @param cap_values Cap y axis to ignore outliers
#' @param cap_SPL_SD Minimum s.d. of samples
#' @param cap_QC_SD Minimum s.d. of QCs
#' @param cap_top_n Cap top n values
#' @param QC_TYPE_fit QC TYPE used for loess fit
#' @param show_driftcorrection Show drift correction
#' @param trend_samples_fun Function used for drift correction. Default 'loess'
#' @param trend_samples_col Color of drift line
#' @param after_correction Show before/after correction
#' @param plot_other_qc Plot all QCS
#' @param show_batches Show batches
#' @param batches_as_shades Show batches as shades
#' @param batch_line_color batch separator color
#' @param batch_shading_color batch shade color
#' @param outputPDF save as PDF
#' @param filename file name of PDF
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

plot_runscatter <- function(data, y_var, transition_filter, filter_exclude = FALSE,
                            cap_values, cap_SPL_SD = 4, cap_QC_SD = 4, cap_top_n = 10, QC_TYPE_fit = "TQC",
                            show_driftcorrection = FALSE, trend_samples_fun = "loess", trend_samples_col ="" , after_correction = FALSE,  plot_other_qc = TRUE,
                            show_batches = FALSE, batches_as_shades = FALSE, batch_line_color = "red1", batch_shading_color = "grey85",
                            outputPDF, filename = "", cols_page, rows_page, annot_scale = 1, paper_orientation = "LANDSCAPE" ,
                            point_transparency=1, point_size=2, point_stroke_width = .8, page_no = NA, y_label_text=NA, silent = FALSE, return_plot_list = TRUE, base_size = 7) {

  y_var_s <- rlang::sym(y_var)
  y_label <- dplyr::if_else(cap_values, paste0(ifelse(is.na(y_label_text), y_var, y_label_text), " (capped at min(", cap_SPL_SD, "x SD[SPL]) ,", cap_QC_SD, "x SD[QC]"), y_var)
  # Re-order QC_TYPE levels to define plot layers, e.g. that QCs are plotted over StudySamples
  data@dataset$QC_TYPE <- factor(as.character(data@dataset$QC_TYPE), pkg.env$qc_type_annotation$qc_type_levels)

  #  filter data
  dat_filt <- data@dataset %>% dplyr::ungroup() %>%
    dplyr::arrange(.data$FEATURE_NAME, .data$RUN_ID) %>%
    dplyr::filter(stringr::str_detect(.data$FEATURE_NAME, paste0("^$|", transition_filter), negate = filter_exclude))

  # cap upper range of dataset to avoid skewness
  dat_filt <- dat_filt %>%
    dplyr::mutate(value =  !!y_var_s)

  dat_filt <- dat_filt %>%
    dplyr::group_by(.data$FEATURE_NAME) %>%
    dplyr::mutate(
      value_max_spl = mean(.data$value[.data$QC_TYPE=="SPL"], na.rm=T) + cap_SPL_SD * sd(.data$value[.data$QC_TYPE=="SPL"], na.rm=T),
      value_max_qc = mean(.data$value[.data$QC_TYPE==QC_TYPE_fit], na.rm=T) + cap_QC_SD * sd(.data$value[.data$QC_TYPE==QC_TYPE_fit]), na.rm=T,
      value_max = max(.data$value_max_spl, .data$value_max_qc, na.rm=T),
      value_mod = dplyr::if_else(.data$value > .data$value_max & cap_values, .data$value_max, .data$value)
    ) %>%
    dplyr::ungroup()

  dat_filt <- dat_filt %>%
    dplyr::group_by(.data$FEATURE_NAME) %>%
    dplyr::arrange(.data$value) %>%
    mutate(
      value = ifelse(dplyr::row_number() < cap_top_n, .data$value[cap_top_n], .data$value)
    ) |>
    dplyr::arrange(.data$FEATURE_NAME, .data$RUN_ID) %>%
    dplyr::ungroup()


  if (outputPDF & !is.na(filename)){
    filename = ifelse(stringr::str_detect(filename, ".pdf"), filename, paste0(filename, ".pdf"))
    if(paper_orientation == "LANDSCAPE")
      pdf(file=filename , onefile=T, paper="A4r", useDingbats=FALSE, width=28/2.54, height=20/2.54)
    else
      pdf(file=filename , onefile=T, paper="A4", useDingbats=FALSE, height=28/2.54, width=20/2.54)
  }

  if(is.na(page_no))
    page_range <- 1:ceiling(dplyr::n_distinct(dat_filt$FEATURE_NAME)/(cols_page * rows_page))
  else
    page_range <- page_no

  if(!silent) print(paste0("Plotting ", max(page_range), " pages..."))

  #p_list <- vector("list", length(page_range))
  p_list <- list()
  for (i in page_range){
    if(!silent) print(paste0("page ", i))
    p <- runscatter_one_page(dat_filt = dat_filt, data= data, d_batches = data@annot_batch_info, cols_page = cols_page, rows_page = rows_page, show_driftcorrection = show_driftcorrection,
                             trend_samples_fun, trend_samples_col, after_correction = after_correction, QC_TYPE_fit = QC_TYPE_fit, outputPDF = outputPDF, page_no = i,
                             point_size = point_size, cap_values = cap_values, point_transparency = point_transparency, annot_scale = annot_scale,
                             show_batches = show_batches, batches_as_shades = batches_as_shades, batch_line_color = batch_line_color, plot_other_qc,
                             batch_shading_color = batch_shading_color, y_label=y_label, base_size=base_size, point_stroke_width=point_stroke_width)
    plot(p)
    p_list[[i]] <- p
  }
  on.exit(if(outputPDF) {dev.off()})
  if(return_plot_list) return(p_list)
}
#' @importFrom ggplot2
runscatter_one_page <- function(dat_filt, data, d_batches, cols_page, rows_page, page_no,
                                show_driftcorrection, after_correction = FALSE, QC_TYPE_fit,cap_values,
                                show_batches, batches_as_shades, batch_line_color, batch_shading_color, trend_samples_fun, trend_samples_col, plot_other_qc,
                                outputPDF, annot_scale, point_transparency, point_size=2, y_label, base_size, point_stroke_width){

  point_size = ifelse(missing(point_size), 2, point_size)
  point_stroke_width <- dplyr::if_else(outputPDF, .3, .2 * (1 + annot_scale/5))


  # subset the dataset with only the rows used for plotting the facets of the selected page
  n_cmpd <- length(unique(dat_filt$ANALYSIS_ID))
  row_start <- n_cmpd * cols_page * rows_page * (page_no - 1) + 1
  row_end <- n_cmpd * cols_page * rows_page * page_no



  dat_subset <- dat_filt %>%
    dplyr::arrange(.data$FEATURE_NAME, .data$RUN_ID) %>%
    dplyr::slice(row_start:row_end)

  dat_subset$QC_TYPE <- forcats::fct_relevel(dat_subset$QC_TYPE, c("SPL", "UBLK", "SBLK", "TQC", "BQC", "RQC", "LTR", "NIST", "PBLK"))
  dat_subset <- dat_subset %>%
    dplyr::arrange(.data$QC_TYPE)

  # https://stackoverflow.com/questions/46327431/facet-wrap-add-geom-hline
  dMax <- dat_subset %>%
    dplyr::group_by(.data$FEATURE_NAME) %>%
    dplyr::summarise(y_max = max(.data$value_mod, na.rm = TRUE)*1.0)

  d_batch_data <- d_batches %>% dplyr::slice(rep(1:dplyr::n(), each = nrow(dMax)))
  d_batch_data$FEATURE_NAME <- rep(dMax$FEATURE_NAME, times = nrow(d_batches))
  d_batch_data <- d_batch_data %>% dplyr::left_join(dMax, by=c("FEATURE_NAME"))


  p <- ggplot2::ggplot(dat_subset, ggplot2::aes_string(x="RUN_ID", label = "ANALYSIS_ID"))

  if (show_batches) {
    if (!batches_as_shades) {
      p <- p + ggplot2::geom_vline(data = d_batch_data %>% dplyr::slice(-1), ggplot2::aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, linetype = "solid", size = .5)
    }
    else {
      d_batches_temp <- d_batch_data %>% dplyr::slice(-1) %>% dplyr::filter(.data$BATCH_NO %% 2 != 1)
      p <- p + ggplot2::geom_rect(data = d_batches_temp, ggplot2::aes(xmin = .data$id_batch_start - 0.5 , xmax = .data$id_batch_end + 0.5, ymin = 0, ymax = .data$y_max, label = .data$BATCH_ID),
                         inherit.aes = FALSE, fill = batch_shading_color, color = NA, alpha = 0.5, linetype = "solid", size = 0.3)
    }
  }

  p <- p +
    ggplot2::geom_point(aes_string(x = "RUN_ID", y= "value_mod", color="QC_TYPE", fill="QC_TYPE", shape="QC_TYPE", group="BATCH_ID"), size=point_size, alpha=point_transparency, stroke = point_stroke_width)



  if(after_correction & show_driftcorrection){
    p <- p +
      ggplot2::geom_line(aes_string(x = "RUN_ID", y= "value_fitted", group = "BATCH_ID"), color = pkg.env$qc_type_annotation$qc_type_fillcol[QC_TYPE_fit], size = .5)
  }

  p <- p +
    ggh4x::facet_wrap2(ggplot2::vars(.data$FEATURE_NAME), scales = "free_y", ncol = cols_page, nrow = rows_page,trim_blank = FALSE) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_color_manual(values=pkg.env$qc_type_annotation$qc_type_col, drop=TRUE) +
    ggplot2::scale_fill_manual(values=pkg.env$qc_type_annotation$qc_type_fillcol, drop=TRUE)+
    ggplot2::scale_shape_manual(values=pkg.env$qc_type_annotation$qc_type_shape, drop=TRUE)
  if(show_driftcorrection){
    if(after_correction) {
      p <- p +
        ggplot2::geom_smooth(data = filter(dat_subset, .data$QC_TYPE == QC_TYPE_fit), ggplot2::aes_string(x = "RUN_ID", y= "value", group = "BATCH_ID"), se=TRUE,
                    colour=pkg.env$qc_type_annotation$qc_type_fillcol[QC_TYPE_fit], fill = pkg.env$qc_type_annotation$qc_type_fillcol[QC_TYPE_fit],
                    method = MASS::rlm, alpha = 0.2, size=0.4) +
        ggplot2::geom_smooth(data = filter(dat_subset, .data$QC_TYPE == "SPL"), ggplot2::aes_string(x = "RUN_ID", y= "value", group = "BATCH_ID"), colour=trend_samples_col, fill = trend_samples_col,
                    method = trend_samples_fun, se=TRUE, alpha = 0.2, size=.4, na.rm = FALSE)

      if(plot_other_qc){
        other_qc <- dplyr::if_else(QC_TYPE_fit == "BQC", "TQC", "BQC")
        other_qc_col <- pkg.env$qc_type_annotation$qc_type_fillcol[other_qc]
        p <- p +
          ggplot2::geom_smooth(data = dplyr::filter(dat_subset, .data$QC_TYPE == other_qc), ggplot2::aes_string(x = "RUN_ID", y= "value", group = "BATCH_ID"), colour=other_qc_col, fill = other_qc_col,
                      method = trend_samples_fun, se=TRUE, alpha = 0.2, size=.4, na.rm = FALSE)
      }
    }
    else {

      p <- p +
        geom_smooth(data = dplyr::filter(dat_subset, .data$QC_TYPE == "SPL"), ggplot2::aes_string(x = "RUN_ID", y= "value", group = "BATCH_ID"), colour=trend_samples_col, fill = trend_samples_col,
                    method = trend_samples_fun, se=TRUE, alpha = 0.2, size=.4, na.rm = FALSE)

      if(plot_other_qc){
        other_qc <- dplyr::if_else(.data$QC_TYPE_fit == "BQC", "TQC", "BQC")
        other_qc_col <- pkg.env$qc_type_annotation$qc_type_fillcol[other_qc]
        p <- p +
          ggplot2::geom_smooth(data = dplyr::filter(dat_subset, .data$QC_TYPE == other_qc), ggplot2::aes_string(x = "RUN_ID", y= "value", group = "BATCH_ID"), colour=other_qc_col, fill = other_qc_col,
                      method = trend_samples_fun, se=TRUE, alpha = 0.2, size=.4, na.rm = FALSE)
      }
    }
  }
  if(cap_values)
    p <- p + ggplot2::geom_hline(data = dMax, ggplot2::aes(yintercept = .data$y_max), color = "#C7C5BF", size = 3, alpha = .3)

  p <- p  +
    #aes(ymin=0) +
    ggplot2::xlab("Injection number") +
    ggplot2::ylab(label = y_label) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.02,0.03))) +
    #expand_limits(y = 0) +
    ggplot2::theme_light(base_size = base_size) +

    ggplot2::theme(plot.title = ggplot2::element_text(size=1, face="bold"),
          strip.text = ggplot2::element_text(size=10*annot_scale, face="bold"),
          strip.background = ggplot2::element_rect(size=0.0001,fill="#8C8C8C"),
          axis.text.x = ggplot2::element_text( size=9*annot_scale, face=NULL),
          axis.text.y = ggplot2::element_text( size=7*annot_scale, face=NULL),
          panel.grid =  ggplot2::element_line(size=0.00001,colour = "#DEDEDE",linetype = "dotted"),
          strip.switch.pad.wrap = ggplot2::unit(0,"mm"))


  return(p)

}



# Define function to plot 1 page
plot_responsecurves_page <- function(dataset,
                                     include_features_containing,
                                     exclude_features_containing,
                                     output_PDF,
                                     response_variable,
                                     regr_max_percent,
                                     pdf_file_name,
                                     rows_page,
                                     columns_page,
                                     point_size,
                                     line_width,
                                     text_scale_factor,
                                     base_size){

  y_var <- rlang::sym(response_variable)
  ggplot2::ggplot(data = dataset,
                  ggplot2::aes(x = .data$RELATIVE_SAMPLE_AMOUNT ,
             y = !!y_var,
             color = .data$RQC_SERIES_ID)) +
    ggpmisc::stat_poly_line(data = subset(dataset, dataset$RELATIVE_SAMPLE_AMOUNT<= (regr_max_percent/100)),
                            ggplot2::aes(x = .data$RELATIVE_SAMPLE_AMOUNT ,
                       y  = !!y_var,
                       color = .data$RQC_SERIES_ID),
                   se = FALSE, na.rm = TRUE, size = line_width, inherit.aes = FALSE) +
    # ggpmisc::stat_poly_eq(
    #   ggplot2::aes(group = .data$RQC_SERIES_ID, label = ggplot2::after_stat(.data$rr.label)),
    #   size = 2* text_scale_factor, rr.digits = 3, vstep = .1) +
    #color = ifelse(after_stat(r.squared) < 0.80, "red", "darkgreen")), size = 1.4) +
    ggplot2::scale_color_manual(values = c("#4575b4", "#91bfdb","#fc8d59",  "#d73027"))+
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::scale_x_continuous(limits = c(0, NA),
                       breaks = c(0,0.5,1,1.5,2,4),
                       labels = scales::percent_format(accuracy = NULL))+
    ggh4x::facet_wrap2(
      ggplot2::vars(.data$FEATURE_NAME),
      scales = "free",
      nrow = rows_page,
      ncol = columns_page,
      trim_blank = FALSE) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::xlab("Sample Amount (Relative to BQC/TQC)")+
    ggplot2::theme_light(base_size = base_size) +
    ggplot2::theme(strip.text = ggplot2::element_text(size=9*text_scale_factor, face="bold"),
                    strip.background = ggplot2::element_rect(size=0.0001,fill="#8C8C8C"))
}


#' Response curves plot
#' @param data MidarExperiment object
#' @param output_PDF Save as PDF
#' @param response_variable Variable to plot
#' @param include_features_containing Filter features containing
#' @param exclude_features_containing Exclude features containing
#' @param regr_max_percent Max relative sample amount to use in regressionb
#' @param pdf_file_name file name of pdf file
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
                                output_PDF,
                                response_variable = "Intensity",
                                include_features_containing = "",
                                exclude_features_containing = "",
                                regr_max_percent = NA,
                                pdf_file_name = "",
                                rows_page = 4,
                                columns_page = 5,
                                point_size = 2,
                                line_width = 1,
                                text_scale_factor = 1,
                                return_plot_list = FALSE, base_size = 7) {

  if (output_PDF & pdf_file_name == "") stop("Please set 'pdf_file_name'")

  rows_page = rows_page
  columns_page = columns_page

  #
  d_rqc <- data@dataset |>
    dplyr::select(tidyselect::any_of(
      c("ANALYSIS_ID", "FEATURE_NAME", "Intensity", "normIntensity")
    )) |>
    dplyr::filter(stringr::str_detect(.data$FEATURE_NAME, paste0("^$|", include_features_containing))) |>
    dplyr::filter(!stringr::str_detect(.data$FEATURE_NAME, paste0("^$|", exclude_features_containing, negate = TRUE))) |>
    dplyr::right_join(data@annot_responsecurves, by = c("ANALYSIS_ID" = "ANALYSIS_ID"))


  regr_max_percent <-
    ifelse(
      is.na(regr_max_percent),
      max(d_rqc$RELATIVE_SAMPLE_AMOUNT * 100),
      regr_max_percent
    )
  #browser()
  d_rqc_grp <- d_rqc %>%
    dplyr::left_join(tibble::tibble(FEATURE_NAME = unique(d_rqc$FEATURE_NAME)) |>
                       mutate(grp = ceiling(
                         row_number() / (rows_page * columns_page)
                       )), by = "FEATURE_NAME") %>%
    dplyr::group_by(.data$grp) %>%
    tidyr::nest() %>%
    dplyr::mutate(plt = purrr::map(
      data,
      function(x)
        plot_responsecurves_page(
          dataset = x,
          include_features_containing = include_features_containing,
          exclude_features_containing = exclude_features_containing,
          output_PDF = output_PDF,
          response_variable = response_variable,
          regr_max_percent = regr_max_percent,
          pdf_file_name = pdf_file_name,
          rows_page = rows_page,
          columns_page = columns_page,
          point_size = point_size,
          line_width = line_width,
          text_scale_factor = text_scale_factor,
          base_size = base_size
        )
    ))

  # Print pages
  #browser()
  if (!output_PDF) {
    if (!return_plot_list)
      d_rqc_grp$plt
  } else{
    pdf(
      file = pdf_file_name,
      onefile = TRUE,
      paper = "A4r",
      width = 11
    )
    print(d_rqc_grp$plt)
    dev.off()
  }
  if (return_plot_list)
    d_rqc_grp
}

#' PCA plot for QC
#' @param data MidarExperiment object
#' @param variable which variable to use for plot
#' @param log_transform log transform data for plot
#' @param dim_x PCA dimension on x axis
#' @param dim_y PCA dimension on y axis
#' @param grouping field used for ellipses
#' @param point_size size of points
#' @param fill_alpha transparency of points
#'
#' @return ggplot2 object
#' @export
plot_pca_qc <- function(data, variable, log_transform, dim_x, dim_y, grouping, point_size = 2, fill_alpha = 0.1) {

  d_wide = data@dataset_QC_filtered  %>%
    filter(.data$QC_TYPE %in% c("BQC", "TQC", "NIST", "LTR", "SPL"), !stringr::str_detect(.data$FEATURE_NAME, "\\(IS"), .data$isQUANTIFIER ) %>%
    dplyr::select("ANALYSIS_ID", "QC_TYPE", "BATCH_ID", "FEATURE_NAME", {{variable}})

  d_filt <- d_wide %>%
    tidyr::pivot_wider(id_cols = "ANALYSIS_ID", names_from = "FEATURE_NAME", values_from = {{variable}})


  d_metadata <- d_wide  %>% dplyr::select("ANALYSIS_ID", "QC_TYPE", "BATCH_ID") |> dplyr::distinct()
  #if(!all(d_filt |> pull(ANALYSIS_ID) == d_metadata |> pull(AnalyticalID))) stop("Data and Metadata missmatch")

  m_raw <- d_filt  |>
    dplyr::select(where(~!any(is.na(.)))) |>
    tibble::column_to_rownames("ANALYSIS_ID") |>
    as.matrix()

  if(log_transform) m_raw <- log2(m_raw)


  # get pca result with annotation
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)
  pca_annot <- pca_res |> broom::augment(d_metadata)
  pca_contrib <- pca_res |> broom::tidy(matrix = "eigenvalues")

  p <- ggplot(data = pca_annot, aes_string(paste0(".fittedPC", dim_x),
                                           paste0(".fittedPC", dim_y),
                                           color = "QC_TYPE",
                                           fill = "QC_TYPE",
                                           shape = "QC_TYPE",
                                           label = "ANALYSIS_ID"
  )) +
    ggplot2::geom_hline(yintercept = 0, size = 0.4, color = "grey80") +
    ggplot2::geom_vline(xintercept = 0, size = 0.4, color = "grey80") +
    ggplot2::stat_ellipse(geom = "polygon", level = 0.95,alpha = fill_alpha, size = 0.3) +
    ggplot2::geom_point(size = point_size)

    p <- p +
      ggplot2::scale_color_manual(values=pkg.env$qc_type_annotation$qc_type_col, drop=TRUE) +
      ggplot2::scale_fill_manual(values=pkg.env$qc_type_annotation$qc_type_fillcol, drop=TRUE)+
      ggplot2::scale_shape_manual(values=pkg.env$qc_type_annotation$qc_type_shape, drop=TRUE)

  p <- p +
    ggplot2::theme_light(base_size = 8) +
    ggplot2::xlab(glue::glue("PC{dim_x} ({round(pca_contrib[[dim_x,'percent']]*100,1)}%)"))+
    ggplot2::ylab(glue::glue("PC{dim_y} ({round(pca_contrib[[dim_y,'percent']]*100,1)}%)"))+
    ggplot2::theme(
      panel.grid = ggplot2::element_line(size = 0.3, color = "grey95"),
      panel.border = ggplot2::element_rect(size = 1, color = "grey70"),
      aspect.ratio=1)

  p
}

plot_pca_pairs <- function(data, variable, dim_range = c(1,8), log_transform = TRUE, grouping = "QC_TYPE", sliding = FALSE, ncol = 3,
                           point_size = 0.5, fill_alpha = 0.1, legend_pos = "right"){


  d_wide = data@dataset  %>% filter(.data$QC_TYPE %in% c("BQC", "TQC", "NIST", "LTR", "SPL"), !stringr::str_detect(.data$FEATURE_NAME, "\\(IS") ) %>%
    dplyr::select("ANALYSIS_ID", "QC_TYPE", "BATCH_ID", "FEATURE_NAME", {{variable}})

  d_filt <- d_wide %>%
    tidyr::pivot_wider(id_cols = "ANALYSIS_ID", names_from = "FEATURE_NAME", values_from = {{variable}})


  d_metadata <- d_wide  %>% dplyr::select("ANALYSIS_ID", "QC_TYPE", "BATCH_ID") |> dplyr::distinct()
  #if(!all(d_filt |> pull(ANALYSIS_ID) == d_metadata |> pull(AnalyticalID))) stop("Data and Metadata missmatch")

  m_raw <- d_filt  |>
    dplyr::select(where(~!any(is.na(.)))) |>
    tibble::column_to_rownames("ANALYSIS_ID") |>
    as.matrix()

  if(log_transform) m_raw <- log2(m_raw)

  # get pca result with annotation
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)
  pca_annot <- pca_res |> broom::augment(d_metadata)
  pca_contrib <- pca_res |> broom::tidy(matrix = "eigenvalues")

  dim_range <- seq(dim_range[1]:dim_range[2])
  if(sliding) dim_range <- dim_range[-length(dim_range)] else
    dim_range <- dim_range[seq(1, length(dim_range), 2)]


  plot_list <- list()

  j = 1

  for(i in dim_range) {

    p <- ggplot(data = pca_annot, aes_string(paste0(".fittedPC", i),
                                             paste0(".fittedPC", i+1),
                                             color = grouping,
                                             shape = grouping,
                                             fill = grouping )) +
      geom_hline(yintercept = 0, size = 0.2, color = "grey80") +
      geom_vline(xintercept = 0, size = 0.2, color = "grey80") +
      ggplot2::stat_ellipse(geom = "polygon", level = 0.95,alpha = fill_alpha, size = 0.2) +
      geom_point(size = point_size)

    p <- p +
      scale_color_manual(values=pkg.env$qc_type_annotation$qc_type_col, drop=TRUE) +
      scale_fill_manual(values=pkg.env$qc_type_annotation$qc_type_fillcol, drop=TRUE)+
      scale_shape_manual(values=pkg.env$qc_type_annotation$qc_type_shape, drop=TRUE)


    p <- p +
      ggplot2::theme_light(base_size = 6) +
      ggplot2::xlab(glue::glue("PC{i} ({round(pca_contrib[[i,'percent']]*100,1)}%)"))+
      ggplot2::ylab(glue::glue("PC{i+1} ({round(pca_contrib[[i+1,'percent']]*100,1)}%)"))+
      ggplot2::theme(
        panel.grid = ggplot2::element_line(size = 0.2, color = "grey95"),
        panel.border = ggplot2::element_rect(size = .5, color = "grey70"),
        aspect.ratio=1)

    plot_list[[j]] <- p + ggplot2::theme(legend.position = "none")
    j = j + 1
  }
  lg <- cowplot::get_legend(plot_list[[1]] + theme(legend.position = legend_pos,
                                                   legend.margin = ggplot2::margin(c(0,0,0,0))))
  pl1 <- cowplot::plot_grid( plotlist = plot_list, ncol = ncol)
  print(cowplot::plot_grid( pl1,lg, ncol = 1, rel_heights = c(1,0.2)))
}

plot_pca_loading_coord <- function(data, variable, log_transform, dim_x, dim_y, top_n, text_size = 1, fill_alpha = 0.1){

  PCx = rlang::sym(paste0("PC",dim_x))
  PCy = rlang::sym(paste0("PC",dim_y))

  d_wide = data@dataset  %>% filter(.data$QC_TYPE %in% c("BQC", "TQC", "NIST", "LTR", "SPL"), !stringr::str_detect(.data$FEATURE_NAME, "\\(IS") ) %>%
    dplyr::select("ANALYSIS_ID", "QC_TYPE", "BATCH_ID", "FEATURE_NAME", {{variable}})

  d_filt <- d_wide %>%
    tidyr::pivot_wider(id_cols = "ANALYSIS_ID", names_from = "FEATURE_NAME", values_from = {{variable}})

  m_raw <- d_filt  |>
    dplyr::select(where(~!any(is.na(.)))) |>
    tibble::column_to_rownames("ANALYSIS_ID") |>
    as.matrix()

  if(log_transform) m_raw <- log2(m_raw)
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)

  d_loading <- pca_res %>%
    broom::tidy(matrix = "rotation") %>%
    tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value")

  d_top_loadings <- d_loading |> dplyr::select("column", !!PCx, !!PCy) |> mutate(vl = (!!PCx)^2 + (!!PCy)^2) |> dplyr::arrange(dplyr::desc(.data$vl)) |> head(top_n)
  x_max <- max(abs(d_top_loadings[, glue::glue("PC{dim_x}")]))
  y_max <- max(abs(d_top_loadings[, glue::glue("PC{dim_y}")]))

  # define arrow style for plotting
  arrow_style <- ggplot2::arrow(
    angle = 20, ends = "first", type = "closed", length = grid::unit(6, "pt")
  )

  p <-  ggplot(d_top_loadings, aes_string(glue::glue("PC{dim_x}"), glue::glue("PC{dim_y}"))) +
    ggplot2::geom_segment(xend = 0, yend = 0, arrow = arrow_style, color = "grey60", size = 0.3) +
    ggplot2::geom_text(
      aes(label = .data$column),
      size = text_size,
      hjust = 0, nudge_x = -0.01,
      color = "#904C2F"
    ) +
    ggplot2::scale_x_continuous(limits = c(-x_max*1.2, x_max*1.2), expand = ggplot2::expansion(mult = .2))+
    ggplot2::scale_y_continuous(limits = c(-y_max*1.2, y_max*1.2), expand = ggplot2::expansion(mult = .2))+
    #coord_fixed(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) + # fix aspect ratio to 1:1
    ggplot2::theme_light(base_size = 7) +
    ggplot2::theme(aspect.ratio=1)
  p
}

plot_pca_loading <- function(data, variable, log_transform, pc_dimensions, top_n, scale_pos_neg = FALSE, point_size = 2, fill_alpha = 0.1){

  d_wide = data@dataset  %>% filter(.data$QC_TYPE %in% c("BQC", "TQC", "NIST", "LTR", "SPL"), !stringr::str_detect(.data$FEATURE_NAME, "\\(IS") ) %>%
    dplyr::select("ANALYSIS_ID", "QC_TYPE", "BATCH_ID", "FEATURE_NAME", {{variable}})

  d_filt <- d_wide %>%
    tidyr::pivot_wider(id_cols = "ANALYSIS_ID", names_from = "FEATURE_NAME", values_from = {{variable}})

  m_raw <- d_filt  |>
    tibble::column_to_rownames("ANALYSIS_ID") |>
    dplyr::select(where(~!any(is.na(.)))) |>
    as.matrix()

  if(log_transform) m_raw <- log2(m_raw)
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)

  d_loading <- pca_res %>%
    broom::tidy(matrix = "rotation") %>%
    tidyr::pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") |>
    dplyr::rename(Compound = .data$column)

  d_loadings_selected <- d_loading |>
    tidyr::pivot_longer(cols = -.data$Compound,  names_to = "PC", values_to = "Value") |>
    dplyr::mutate(PC = as.numeric(stringr::str_remove(.data$PC, "PC"))) |>
    filter(.data$PC %in% pc_dimensions)


  d_loadings_selected <- d_loadings_selected |>
    dplyr::rowwise() |>
    dplyr::mutate(direction = if_else(.data$Value < 0, "neg", "pos"),
           Value = dplyr::if_else(!scale_pos_neg, abs(.data$Value), .data$Value),
           abs_value = abs(.data$Value)) |>
    group_by(.data$PC) |>
    dplyr::arrange(.data$abs_value) |>
    dplyr::slice_max(order_by = .data$abs_value, n = .data$top_n) |>
    ungroup() |>
    tidyr::unite("E", .data$Compound, .data$PC, remove = FALSE) |>
    mutate(PC = as.factor(.data$PC),
           E = forcats::fct_reorder(.data$E, .data$abs_value))


  p <- ggplot(d_loadings_selected, ggplot2::aes(x = .data$E, y = .data$Value, color = .data$direction, fill = .data$direction)) +
    ggplot2::geom_col()+
    ggplot2::facet_wrap(ggplot2::vars(.data$PC), nrow=1,scales = "free") +
    ggplot2::scale_x_discrete(labels=d_loadings_selected$Compound, breaks=d_loadings_selected$E) +
    ggplot2::scale_color_manual(values = c("neg" = "blue", "pos" = "red")) +
    ggplot2::scale_fill_manual(values = c("neg" = "blue", "pos" = "red")) +
    ggplot2::coord_flip() +
    ggplot2::theme_light(base_size = 8)

  p

}

