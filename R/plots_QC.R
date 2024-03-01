#' RunScatter plot
#'
#' @param data MidarExperiment object
#' @param plot_variable Variable to plot
#' @param feature_filter Filter features containing
#' @param filter_exclude Exclude or include feature_filter
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

plot_runscatter <- function(data, plot_variable = c("feature_intensity", "feature_norm_intensity", "feature_conc"), feature_filter = "", filter_exclude = FALSE,
                            cap_outliers = FALSE, cap_qc_iqr_factor = 1.5, cap_spl_iqr_factor = 1.5, cap_top_n_values = NA, qc_type_fit = "BQC",
                            show_driftcorrection = FALSE, show_trend_samples, trend_samples_fun = "loess", trend_samples_col ="" , after_correction = FALSE,  plot_other_qc = TRUE,
                            show_batches = FALSE, batches_as_shades = TRUE, batch_line_color = "#9dbecf", batch_shading_color = "grey90",
                            outputPDF = FALSE, filename = "", cols_page = 4, rows_page = 3, annot_scale = 1, paper_orientation = "LANDSCAPE" ,
                            point_transparency=1, point_size=2, point_stroke_width = .8, page_no = NA, y_label_text=NA, silent = FALSE, return_plot_list = FALSE, base_size = 7, show_gridlines = FALSE) {

  plot_variable <- rlang::arg_match(plot_variable)
  plot_variable_s <- rlang::sym(plot_variable)
  y_label <- dplyr::if_else(cap_outliers, paste0(ifelse(is.na(y_label_text), plot_variable, y_label_text), " (capped at min(", cap_spl_iqr_factor, "x IQR+Q3[SPL]) ,", cap_qc_iqr_factor, "x IQR+Q3[QC]"), plot_variable)

  # Re-order qc_type levels to define plot layers, e.g. that QCs are plotted over StudySamples
  data@dataset$qc_type <- factor(as.character(data@dataset$qc_type), pkg.env$qc_type_annotation$qc_type_levels)

  #  filter data
  dat_filt <- data@dataset %>% dplyr::ungroup() %>%
    dplyr::arrange(.data$feature_name, .data$run_id) %>%
    dplyr::filter(stringr::str_detect(.data$feature_name, paste0("^$|", feature_filter), negate = filter_exclude))

  # Cap upper outliers to reduce skewness
  dat_filt <- dat_filt %>%
    dplyr::mutate(value =  !!plot_variable_s)


  dat_filt <- dat_filt %>%
    dplyr::group_by(.data$feature_name) %>%
    dplyr::mutate(
      #value_max_spl = mean(.data$value[.data$qc_type=="SPL"], na.rm=T) + cap_SPL_SD * sd(.data$value[.data$qc_type=="SPL"], na.rm=T),
      #value_max_qc = mean(.data$value[.data$qc_type==qc_type_fit], na.rm=T) + cap_QC_SD * sd(.data$value[.data$qc_type==qc_type_fit]), na.rm=T,
      value_max_spl = quantile(.data$value[.data$qc_type == "SPL"], 0.75, na.rm = TRUE) + cap_spl_iqr_factor * IQR(.data$value[.data$qc_type=="SPL"], na.rm = TRUE),
      value_max_qc = quantile(.data$value[.data$qc_type == qc_type_fit], 0.75, na.rm = TRUE) + cap_qc_iqr_factor * IQR(.data$value[.data$qc_type == qc_type_fit], na.rm = TRUE),

      value_max = max(.data$value_max_spl, .data$value_max_qc, na.rm = TRUE),
      value_mod = dplyr::if_else(cap_outliers & .data$value > .data$value_max, suppressWarnings(max(.data$value[.data$value < .data$value_max], na.rm = TRUE)), .data$value)
    ) %>%
    dplyr::ungroup()

  dat_filt <- dat_filt %>%
    dplyr::group_by(.data$feature_name) %>%
    dplyr::arrange(.data$value) %>%
    mutate(
      value = ifelse(dplyr::row_number() < cap_top_n_values, .data$value[cap_top_n_values], .data$value)
    ) |>
    dplyr::arrange(.data$feature_name, .data$run_id) %>%
    dplyr::ungroup()


  if (outputPDF & !is.na(filename)){
    filename = ifelse(stringr::str_detect(filename, ".pdf"), filename, paste0(filename, ".pdf"))
    if(paper_orientation == "LANDSCAPE")
      pdf(file=filename , onefile=T, paper="A4r", useDingbats=FALSE, width=28/2.54, height=20/2.54)
    else
      pdf(file=filename , onefile=T, paper="A4", useDingbats=FALSE, height=28/2.54, width=20/2.54)
  }

  if(is.na(page_no))
    page_range <- 1:ceiling(dplyr::n_distinct(dat_filt$feature_name)/(cols_page * rows_page))
  else
    page_range <- page_no

  if(!silent) print(paste0("Plotting ", max(page_range), " pages..."))

  #p_list <- vector("list", length(page_range))
  p_list <- list()
  for (i in page_range){
    if(!silent) print(paste0("page ", i))
    p <- runscatter_one_page(dat_filt = dat_filt, data= data, d_batches = data@annot_batch_info, cols_page = cols_page, rows_page = rows_page, show_driftcorrection = show_driftcorrection,
                             show_trend_samples, trend_samples_fun, trend_samples_col, after_correction = after_correction, qc_type_fit = qc_type_fit, outputPDF = outputPDF, page_no = i,
                             point_size = point_size, cap_outliers = cap_outliers, point_transparency = point_transparency, annot_scale = annot_scale,
                             show_batches = show_batches, batches_as_shades = batches_as_shades, batch_line_color = batch_line_color, plot_other_qc,
                             batch_shading_color = batch_shading_color, y_label=y_label, base_size=base_size, point_stroke_width=point_stroke_width, show_grid = show_gridlines)
    plot(p)
    p_list[[i]] <- p
  }
  on.exit(if(outputPDF) {dev.off()})
  if(return_plot_list) return(p_list)
}
#' @importFrom ggplot2 Stat
runscatter_one_page <- function(dat_filt, data, d_batches, cols_page, rows_page, page_no,
                                show_driftcorrection, after_correction = FALSE, qc_type_fit,cap_outliers,
                                show_batches, batches_as_shades, batch_line_color, batch_shading_color, show_trend_samples, trend_samples_fun, trend_samples_col, plot_other_qc,
                                outputPDF, annot_scale, point_transparency, point_size=2, y_label, base_size, point_stroke_width, show_grid){

  point_size = ifelse(missing(point_size), 2, point_size)
  point_stroke_width <- dplyr::if_else(outputPDF, .3, .2 * (1 + annot_scale/5))


  # subset the dataset with only the rows used for plotting the facets of the selected page
  n_cmpd <- length(unique(dat_filt$analysis_id))
  row_start <- n_cmpd * cols_page * rows_page * (page_no - 1) + 1
  row_end <- n_cmpd * cols_page * rows_page * page_no



  dat_subset <- dat_filt %>%
    dplyr::arrange(.data$feature_name, .data$run_id) %>%
    dplyr::slice(row_start:row_end)

  dat_subset$qc_type <- forcats::fct_relevel(dat_subset$qc_type, c("SPL", "UBLK", "SBLK", "TQC", "BQC", "RQC", "LTR", "NIST", "PBLK"))
  dat_subset <- dat_subset %>%
    dplyr::arrange(.data$qc_type)




  # https://stackoverflow.com/questions/46327431/facet-wrap-add-geom-hline
  dMax <- dat_subset %>%
    dplyr::group_by(.data$feature_name) %>%
    dplyr::summarise(y_max = max(.data$value_mod, na.rm = TRUE)*1.0)


  d_batch_data <- d_batches %>% dplyr::slice(rep(1:dplyr::n(), each = nrow(dMax)))
  d_batch_data$feature_name <- rep(dMax$feature_name, times = nrow(d_batches))
  d_batch_data <- d_batch_data %>% dplyr::left_join(dMax, by=c("feature_name"))
  p <- ggplot2::ggplot(dat_subset, ggplot2::aes_string(x="run_id", label = "analysis_id"))

  #browser()
  if (show_batches) {

    if (!batches_as_shades) {
      d_batches_temp <- d_batch_data  |> filter(.data$id_batch_start != 1)
      p <- p + ggplot2::geom_vline(data = d_batches_temp, ggplot2::aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, linetype = "solid", size = .5)
    }
    else {
      d_batches_temp <- d_batch_data  %>% dplyr::filter(.data$batch_no %% 2 != 1)
      p <- p + ggplot2::geom_rect(data = d_batches_temp, ggplot2::aes(xmin = .data$id_batch_start - 0.5 , xmax = .data$id_batch_end + 0.5, ymin = 0, ymax = .data$y_max, label = .data$batch_id),
                         inherit.aes = FALSE, fill = batch_shading_color, color = NA, alpha = 0.5, linetype = "solid", size = 0.3)
    }
  }

  if(cap_outliers)
    p <- p + ggplot2::geom_hline(data = dMax, ggplot2::aes(yintercept = .data$y_max), color = "#fa9b9b", size = 3, alpha = .4)

  p <- p +
    ggplot2::geom_point(aes_string(x = "run_id", y= "value_mod", color="qc_type", fill="qc_type", shape="qc_type", group="batch_id"), size=point_size, alpha=point_transparency, stroke = point_stroke_width)



  if(after_correction & show_driftcorrection){
    p <- p +
      ggplot2::geom_line(aes_string(x = "run_id", y= "CURVE_Y_PREDICTED", group = "batch_id"), color = pkg.env$qc_type_annotation$qc_type_fillcol[qc_type_fit], size = .5)
  }

  p <- p +
    ggh4x::facet_wrap2(ggplot2::vars(.data$feature_name), scales = "free_y", ncol = cols_page, nrow = rows_page,trim_blank = FALSE) +
    ggplot2::expand_limits(y = 0) +
    ggplot2::scale_color_manual(values=pkg.env$qc_type_annotation$qc_type_col, drop=TRUE) +
    ggplot2::scale_fill_manual(values=pkg.env$qc_type_annotation$qc_type_fillcol, drop=TRUE)+
    ggplot2::scale_shape_manual(values=pkg.env$qc_type_annotation$qc_type_shape, drop=TRUE)

  if(show_driftcorrection){
    if(after_correction) {
      p <- p +
        ggplot2::geom_smooth(data = filter(dat_subset, .data$qc_type == qc_type_fit), ggplot2::aes_string(x = "run_id", y= "value", group = "batch_id"), se=TRUE,
                    colour=pkg.env$qc_type_annotation$qc_type_fillcol[qc_type_fit], fill = pkg.env$qc_type_annotation$qc_type_fillcol[qc_type_fit],
                    method = MASS::rlm, alpha = 0.35, size=0.8)

      if(show_trend_samples){
        p <- p + ggplot2::geom_smooth(data = filter(dat_subset, .data$qc_type == "SPL"), ggplot2::aes_string(x = "run_id", y= "value", group = "batch_id"), colour=trend_samples_col, fill = trend_samples_col, linetype = "dashed",
                    method = trend_samples_fun, se=FALSE, alpha = 0.2, size=.8, na.rm = TRUE)
      }

      if(plot_other_qc){
        other_qc <- dplyr::if_else(qc_type_fit == "BQC", "TQC", "BQC")
        other_qc_col <- pkg.env$qc_type_annotation$qc_type_fillcol[other_qc]
        p <- p +
          ggplot2::geom_smooth(data = dplyr::filter(dat_subset, .data$qc_type == other_qc), ggplot2::aes_string(x = "run_id", y= "value", group = "batch_id"), colour=other_qc_col, fill = other_qc_col,
                      method = trend_samples_fun, se=TRUE, alpha = 0.2, size=.4, na.rm = FALSE)
      }
    }
    else {

      p <- p +
        geom_smooth(data = dplyr::filter(dat_subset, .data$qc_type == "SPL"), ggplot2::aes_string(x = "run_id", y= "value", group = "batch_id"), colour=trend_samples_col, fill = trend_samples_col,
                    method = trend_samples_fun, se=TRUE, alpha = 0.2, size=.4, na.rm = FALSE)

      if(plot_other_qc){
        other_qc <- dplyr::if_else(.data$qc_type_fit == "BQC", "TQC", "BQC")
        other_qc_col <- pkg.env$qc_type_annotation$qc_type_fillcol[other_qc]
        p <- p +
          ggplot2::geom_smooth(data = dplyr::filter(dat_subset, .data$qc_type == other_qc), ggplot2::aes_string(x = "run_id", y= "value", group = "batch_id"), colour=other_qc_col, fill = other_qc_col,
                      method = trend_samples_fun, se=TRUE, alpha = 0.2, size=.4, na.rm = FALSE)
      }
    }
  }


  p <- p  +
    #aes(ymin=0) +
    ggplot2::xlab("Injection number") +
    ggplot2::ylab(label = y_label) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.02,0.03))) +
    #expand_limits(y = 0) +
    ggplot2::theme_bw(base_size = base_size) +

    ggplot2::theme(plot.title = ggplot2::element_text(size=1, face="bold"),
          strip.text = ggplot2::element_text(size=10*annot_scale, face="bold"),
            strip.background = ggplot2::element_rect(size=0.0001,fill="#00283d"),
          strip.text.x = ggplot2::element_text(color = "white"),
          axis.text.x = ggplot2::element_text( size=9*annot_scale, face=NULL),
          axis.text.y = ggplot2::element_text( size=7*annot_scale, face=NULL),
          panel.grid.major =  ggplot2::element_blank(),
          panel.grid.minor =  ggplot2::element_blank(),
          strip.switch.pad.wrap = ggplot2::unit(0,"mm"))

  if(show_grid)
    p <- p  + ggplot2::theme(panel.grid.major =  ggplot2::element_line(size=0.3,colour = "grey88",linetype = "dashed"))



  return(p)

}



# Define function to plot 1 page
plot_responsecurves_page <- function(dataset,
                                     include_features_containing,
                                     exclude_features_containing,
                                     output_PDF,
                                     response_variable,
                                     regr_max_percent,
                                     pdf_filename,
                                     rows_page,
                                     columns_page,
                                     point_size,
                                     line_width,
                                     text_scale_factor,
                                     base_size){

  plot_variable <- rlang::sym(response_variable)
  ggplot2::ggplot(data = dataset,
                  ggplot2::aes(x = .data$relative_sample_amount ,
             y = !!plot_variable,
             color = .data$rqc_series_id)) +
    ggpmisc::stat_poly_line(data = subset(dataset, dataset$relative_sample_amount<= (regr_max_percent/100)),
                            ggplot2::aes(x = .data$relative_sample_amount ,
                       y  = !!plot_variable,
                       color = .data$rqc_series_id),
                   se = FALSE, na.rm = TRUE, size = line_width, inherit.aes = FALSE ) +
    ggpmisc::stat_poly_eq(
      ggplot2::aes(group = .data$rqc_series_id, label = ggplot2::after_stat(.data$rr.label)),
      size = 2* text_scale_factor, rr.digits = 3, vstep = .1) +
    #color = ifelse(after_stat(r.squared) < 0.80, "red", "darkgreen")), size = 1.4) +
    ggplot2::scale_color_manual(values = c("#4575b4", "#91bfdb","#fc8d59",  "#d73027"))+
    ggplot2::scale_y_continuous(limits = c(0, NA)) +
    ggplot2::scale_x_continuous(limits = c(0, NA),
                       breaks = c(0,0.5,1,1.5,2,4),
                       labels = scales::percent_format(accuracy = NULL))+
    ggh4x::facet_wrap2(
      ggplot2::vars(.data$feature_name),
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
#' @param pdf_filename file name of pdf file
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
                                response_variable = "feature_intensity",
                                include_features_containing = "",
                                exclude_features_containing = "",
                                regr_max_percent = NA,
                                pdf_filename = "",
                                rows_page = 4,
                                columns_page = 5,
                                point_size = 2,
                                line_width = 1,
                                text_scale_factor = 1,
                                return_plot_list = FALSE, base_size = 7) {

  if (output_PDF & pdf_filename == "") stop("Please set 'pdf_filename'")

  rows_page = rows_page
  columns_page = columns_page

  #
  d_rqc <- data@dataset |>
    dplyr::select(tidyselect::any_of(
      c("analysis_id", "feature_name", "feature_intensity", "feature_norm_intensity")
    )) |>
    dplyr::filter(stringr::str_detect(.data$feature_name, paste0("^$|", include_features_containing))) |>
    dplyr::filter(!stringr::str_detect(.data$feature_name, paste0("^$|", exclude_features_containing, negate = TRUE))) |>
    dplyr::right_join(data@annot_responsecurves, by = c("analysis_id" = "analysis_id"))


  regr_max_percent <-
    ifelse(
      is.na(regr_max_percent),
      max(d_rqc$relative_sample_amount * 100),
      regr_max_percent
    )
  #browser()
  d_rqc_grp <- d_rqc %>%
    dplyr::left_join(tibble::tibble(feature_name = unique(d_rqc$feature_name)) |>
                       mutate(grp = ceiling(
                         row_number() / (rows_page * columns_page)
                       )), by = "feature_name") %>%
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
          pdf_filename = pdf_filename,
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
      file = pdf_filename,
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


#' Plot comparison of analytical CV before and after ISTD-normalization
#' @param data MidarExperiment object
#' @param x Variable of QC metrics table used for x axis
#' @param y Variable of QC metrics table used for y axis
#' @param only_quantifier Show only quantifier
#' @param xlim Lower and upper limits of the x axis as vector
#' @param ylim Lower and upper limits of the y axis as vector
#' @param ncol Numer of columns with facets
#' @param scale_factor Scale factor for text labels
#' @param point_size Point size
#' @param with_histogram Show side histograms indicating distribution
#' @return ggplot2 object
#' @export


plot_x_vs_y <- function(data, x, y, only_quantifier = TRUE, xlim=c(0,NA), ylim=c(0,NA), ncol = 5, scale_factor = 1, point_size = 3, with_histogram=FALSE){


  d_QC <- data@metrics_qc |> filter(.data$valid_integration)

  if(only_quantifier) d_QC <- d_QC |> filter(.data$is_quantifier, .data$valid_integration)

  x_sym <- rlang::ensym(x)
  y_sym <- rlang::ensym(y)

  if(!c("LipidClass") %in% names(d_QC)) stop("This function currently only works with lipidomics data. Please add lipid names/class with the function `add_lipid_class_transition` before calling this function.")
  # get max value for the pair for each lipid class (so that 45deg line will be shown in the square plots)
  d_QC <- d_QC %>% group_by(.data$LipidClass) %>% mutate(xy_max = max(!!sym(x),!!sym(y),na.rm = TRUE)) %>% ungroup()

  point_size = 0.75 * point_size
  point_linewidth =1 * scale_factor/2

  # Use geom_polygon for coloring the two areas
  g <- ggplot(data=d_QC, aes(x=!!x_sym, y=!!y_sym)) +
    facet_wrap(vars(.data$LipidClass), scales = "free",ncol = ncol) +
    geom_point(aes(colour=as.factor(.data$TransitionGroup),fill=as.factor(.data$TransitionGroup)), size = point_size, alpha=0.7, stroke = point_linewidth) +
    geom_point(aes(x=.data$xy_max,y=.data$xy_max), fill="#5cf442", shape = ".") +
    labs(fill = "Transition Group", colour = "Transition Group")

  g <- g +
    scale_x_continuous(limits = xlim) +
    scale_y_continuous(limits = ylim) +
    aes(ymin=0,xmin=0) +
    geom_abline(intercept = 0, slope = 1) +
    #scale_colour_brewer(palette = "Set1") +
    scale_shape_manual(na.value=NA, values = c(15,0,4), drop=FALSE, name="Transition") +
    theme(plot.title = element_text(size=01, face="bold"),
          strip.text = element_text(size=11*scale_factor, margin = margin(1,1,0,0), lineheight=0),
          strip.background = element_rect(size=0.0001),
          axis.text = element_text(size=09*scale_factor, face=NULL),
          axis.title = element_text(size=12*scale_factor, face="bold"),
          panel.grid =  element_line(size=0.001),
          strip.switch.pad.wrap = unit(1,"mm"),
          legend.position="right") +
    geom_hline(yintercept=20, linetype="dashed", color = "#ed5578") +
    geom_vline(xintercept=20, linetype="dashed", color = "#ed5578")

  if(with_histogram){
    g <- g + theme(legend.position= c(0.9, 0.3))
    g <- ggExtra::ggMarginal(g, type = "histogram",margins = "y")
  }
  return(g)
}




#' PCA plot for QC
#' @param data MidarExperiment object
#' @param variable which variable to use for plot
#' @param log_transform log transform data for plot
#' @param dim_x PCA dimension on x axis
#' @param dim_y PCA dimension on y axis
#' @param point_size size of points
#' @param point_alpha transparency of points
#' @param ellipse_alpha transparency of ellipse fill
#' @param font_base_size Base font size for plot text elements
#' @param remove_istds exclude internal standards
#'
#' @return ggplot2 object
#' @export
plot_pca_qc <- function(data, variable, dim_x, dim_y, log_transform, remove_istds, point_size = 2, point_alpha = 0.5, ellipse_alpha = 0.8, font_base_size = 8) {


  d_wide <- data@dataset_QC_filtered

  #TODO: (IS as criteria for ISTD.. dangerous...
  if(remove_istds)  d_wide <- d_wide |> filter(!is_istd) # !stringr::str_detect(.data$feature_name, "\\(IS")

  d_wide <- d_wide |>  filter(.data$qc_type %in% c("BQC", "TQC", "NIST", "LTR", "SPL"), .data$is_quantifier ) |>
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_name", {{variable}})

  d_filt <- d_wide %>%
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_name", values_from = {{variable}})


  #if(!all(d_filt |> pull(analysis_id) == d_metadata |> pull(AnalyticalID))) stop("Data and Metadata missmatch")

  #ToDo: warning when rows/cols with NA
  d_clean <- d_filt  |>
    filter(if_any(dplyr::where(is.numeric), ~ !is.na(.))) |>
    dplyr::select(where(~!any(is.na(.) | is.nan(.) | is.infinite(.) | . <= 0)))


  d_metadata <- d_wide  %>%
    dplyr::select("analysis_id", "qc_type", "batch_id") |>
    dplyr::distinct() |>
    dplyr::right_join(d_clean |> dplyr::select("analysis_id") |> distinct(), by = c("analysis_id"))


  m_raw <- d_clean |>
    tibble::column_to_rownames("analysis_id") |>
    as.matrix()



  if(log_transform) m_raw <- log2(m_raw)
  # get pca result with annotation
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)
  pca_annot <- pca_res |>
    broom::augment(d_metadata)

  pca_annot$qc_type <- forcats::fct_relevel(pca_annot$qc_type, c("SPL", "UBLK", "SBLK", "TQC", "BQC", "RQC", "LTR", "NIST", "PBLK"))
  pca_annot <- pca_annot %>%
    dplyr::arrange(.data$qc_type)

  pca_contrib <- pca_res |> broom::tidy(matrix = "eigenvalues")

  p <- ggplot(data = pca_annot, aes_string(paste0(".fittedPC", dim_x),
                                           paste0(".fittedPC", dim_y),
                                           color = "qc_type",
                                           fill = "qc_type",
                                           shape = "qc_type",
                                           label = "analysis_id"
  )) +
    ggplot2::geom_hline(yintercept = 0, size = 0.5, color = "grey80", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, size = 0.5, color = "grey80", linetype = "dashed") +
    ggplot2::stat_ellipse(geom = "polygon", level = 0.95,alpha = ellipse_alpha, size = 0.3) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha)

    p <- p +
      ggplot2::scale_color_manual(values=pkg.env$qc_type_annotation$qc_type_col, drop=TRUE) +
      ggplot2::scale_fill_manual(values=pkg.env$qc_type_annotation$qc_type_fillcol, drop=TRUE)+
      ggplot2::scale_shape_manual(values=pkg.env$qc_type_annotation$qc_type_shape, drop=TRUE)

  p <- p +
    ggplot2::theme_bw(base_size = font_base_size) +
    ggplot2::xlab(glue::glue("PC{dim_x} ({round(pca_contrib[[dim_x,'percent']]*100,1)}%)"))+
    ggplot2::ylab(glue::glue("PC{dim_y} ({round(pca_contrib[[dim_y,'percent']]*100,1)}%)"))+
    ggplot2::theme(
      panel.grid.major =  ggplot2::element_blank(),
      panel.grid.minor =  ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(size = 1, color = "grey40"),
      axis.text.x = ggplot2::element_text( size=font_base_size, face=NULL),
      axis.text.y = ggplot2::element_text( size=font_base_size, face=NULL),
      axis.title.x = ggplot2::element_text( size=font_base_size*1.2, face=NULL),
      axis.title.y = ggplot2::element_text( size=font_base_size*1.2, face=NULL),
      aspect.ratio=1)

  p
}

plot_pca_pairs <- function(data, variable, dim_range = c(1,8), log_transform = TRUE, grouping = "qc_type", sliding = FALSE, ncol = 3,
                           point_size = 0.5, fill_alpha = 0.1, legend_pos = "right"){


  d_wide = data@dataset  %>% filter(.data$qc_type %in% c("BQC", "TQC", "NIST", "LTR", "SPL"), !stringr::str_detect(.data$feature_name, "\\(IS") ) %>%
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_name", {{variable}})

  d_filt <- d_wide %>%
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_name", values_from = {{variable}})


  d_metadata <- d_wide  %>% dplyr::select("analysis_id", "qc_type", "batch_id") |> dplyr::distinct()
  #if(!all(d_filt |> pull(analysis_id) == d_metadata |> pull(AnalyticalID))) stop("Data and Metadata missmatch")

  m_raw <- d_filt  |>
    dplyr::select(where(~!any(is.na(.)))) |>
    tibble::column_to_rownames("analysis_id") |>
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

  d_wide = data@dataset  %>% filter(.data$qc_type %in% c("BQC", "TQC", "NIST", "LTR", "SPL"), !stringr::str_detect(.data$feature_name, "\\(IS") ) %>%
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_name", {{variable}})

  d_filt <- d_wide %>%
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_name", values_from = {{variable}})

  m_raw <- d_filt  |>
    dplyr::select(where(~!any(is.na(.)))) |>
    tibble::column_to_rownames("analysis_id") |>
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

plot_pca_loading <- function(data, variable, log_transform, pc_dimensions, top_n, remove_istds, vertical_bars = FALSE, scale_pos_neg = FALSE, point_size = 2, fill_alpha = 0.1){

  d_wide = data@dataset_QC_filtered  %>% filter(.data$qc_type %in% c("BQC", "TQC", "NIST", "LTR", "SPL"))

  if(remove_istds)  d_wide <- d_wide |>  filter(!is_istd) # !stringr::str_detect(.data$feature_name, "\\(IS")

  d_filt <- d_wide |>
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_name", {{variable}}) |>
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_name", values_from = {{variable}})

  m_raw <- d_filt  |>
    tibble::column_to_rownames("analysis_id") |>
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
    group_by(.data$PC)

  if(vertical_bars)
    d_loadings_selected <- d_loadings_selected |> dplyr::arrange(.data$abs_value)
  else
    d_loadings_selected <- d_loadings_selected |> dplyr::arrange(desc(.data$abs_value))

  d_loadings_selected <- d_loadings_selected  |>
    dplyr::slice_max(order_by = .data$abs_value, n = top_n) |>
    ungroup() |>
    tidyr::unite("Feature", .data$Compound, .data$PC, remove = FALSE) |>
    mutate(PC = as.factor(.data$PC),
           Feature = forcats::fct_reorder(.data$Feature, .data$abs_value))


  p <- ggplot(d_loadings_selected, ggplot2::aes(x = .data$Feature, y = .data$Value, color = .data$direction, fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(ggplot2::vars(.data$PC), scales = "free", ncol = ifelse(vertical_bars, 1, length(pc_dimensions))) +
    ggplot2::scale_x_discrete(labels=d_loadings_selected$Compound, breaks=d_loadings_selected$Feature) +
    ggplot2::scale_color_manual(values = c("neg" = "#75CEFF", "pos" = "#FFA166")) +
    ggplot2::scale_fill_manual(values = c("neg" = "#75CEFF", "pos" = "#FFA166")) +
    ggplot2::labs(x = "Feature", y = "Loading") +
    ggplot2::theme_bw(base_size = 8)

  if (!vertical_bars) {
    p <- p +
      ggplot2::coord_flip()
  } else
  {
    p <- p +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust=0.5, hjust=1)) +
      scale_x_discrete(limits = rev)
  }

  p

}

