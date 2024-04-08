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
                              factor_a=NA,
                              factor_b=NA) {

  d_temp <- data$dataset %>% dplyr::select(run_id, batch_id, analysis_id, qc_type, sample_id) %>% distinct()

  d_temp$qc_type <- factor(d_temp$qc_type, c("EQC", "SST", "MBLK", "SBLK", "UBLK", "PBLK", "RQC",   "LTR", "NIST", "TQC" ,"BQC", "SPL"))

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

  #d_temp <- d_temp %>% mutate(sample_category = forcats::fct_rev(forcats::as_factor(sample_category)))
  d_temp <- d_temp %>% filter(str_detect(as.character(sample_category), qc_type_subet)) %>%
    {if (show_qc_dataset_only|qc_type_subet != "") droplevels(.) else .}

  d_batch_info <- data@annot_batches

  # round to next 10 and divide by number of breaks will be showb
  scale_dataset_size = 10^ceiling(log10(max(d_temp$run_id)))/10
  p <- ggplot(d_temp, aes(x = run_id, y = rev(qc_type), color = sample_category))
  if(show_batches){
    if (batches_as_shades){
      d_batch_2nd <- d_batch_info %>% dplyr::slice(-1) %>% filter(batch_no %% 2 != 1)
      p <- p + geom_rect(data = d_batch_2nd, aes(xmin = id_batch_start-0.5 , xmax = id_batch_end+0.5, ymin = -Inf, ymax = Inf),
                         inherit.aes=FALSE, fill=batch_shading_color, alpha = 0.1, color= NA,linetype="solid", size=0.5)
    }
  }
  p <- p + geom_segment(
    aes(
      x = run_id,
      xend = run_id,
      y = as.integer(qc_type) - 0.4,
      yend = as.integer(qc_type) + 0.4
    ),
    size = segment_width) +
    labs(x="Analysis order",
         y = "Sample Type") +
    scale_x_continuous(breaks = seq(0, max(d_temp$run_id), scale_dataset_size)) +
    scale_y_discrete(
      limits = c(levels(d_temp$qc_type)),
      expand = expansion(0.05,0.05)
    ) +
    #scale_color_manual(values =  data$sample_category) +
    theme_light(base_size = base_size) +
    theme( panel.grid.major.y = element_blank(),
           panel.grid.major.x =  element_line(colour="grey80", linetype="solid", size=.5),
           panel.grid.minor.x =  element_line(colour="grey90", linetype="dotted", size=0.25),
           panel.border = element_rect(size = 2),
           axis.title = element_text(face = "bold",size = base_size),
           axis.text.x = element_text(face = "plain",size = base_size),
           axis.text.y = element_text(face = "bold",size = base_size),
           legend.position = "none")


  if(show_batches){
    if (!batches_as_shades){
      p <- p + geom_vline(data = d_batch_info %>% dplyr::slice(-1), aes(xintercept=id_batch_start-0.5), colour=batch_line_color, size=0.5)
    }
  }
 # if (factor_a != "" & factor_b !=""){

  #  p <- p + scale_color_brewer(palette="Dark2")
  #} else {
    p <- p + scale_color_manual(values =  pkg.env$qc_type_annotation$qc_type_col)
  #}
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
                            trend_samples_col ="" ,
                            after_correction = FALSE,
                            plot_other_qc = TRUE,
                            show_batches = TRUE,
                            batches_as_shades = FALSE,
                            batch_line_color = "#9dbecf",
                            batch_shading_color = "grey90",
                            outputPDF = FALSE,
                            filename = "",
                            cols_page = 3,
                            rows_page = 3,
                            annot_scale = 1.0,
                            paper_orientation = "LANDSCAPE" ,
                            point_transparency= 1,
                            point_size=2,
                            point_stroke_width = 1,
                            page_no = NA,
                            y_label_text=NA,
                            silent = TRUE,
                            return_plot_list = FALSE,
                            base_size = 12,
                            show_gridlines = FALSE) {

  if(nrow(data@dataset) < 1) stop("No data available. Please import data and metadata first.")

  if (use_filt_data){
    dat_filt <- data@dataset_filtered %>% dplyr::ungroup()
    if(nrow(dat_filt) < 1) stop("Data has not been qc filtered. Please apply `apply_qc_filter` first.")
  } else {
    dat_filt <- data@dataset %>% dplyr::ungroup()
  }

  plot_var <- rlang::arg_match(plot_var)
  plot_var_s <- rlang::sym(plot_var)
  y_label <- dplyr::if_else(cap_outliers, paste0(ifelse(is.na(y_label_text), plot_var, y_label_text), " (capped min(", cap_spl_iqr_factor, "x IQR+Q3[SPL]) ,", cap_qc_iqr_factor, "x IQR+Q3[QC]"), stringr::str_remove(plot_var, "feature\\_"))




  # Subset data ----

  if(all(!is.na(feature_incl_filt)) & all(feature_incl_filt != "")){
    if(length(feature_incl_filt) == 1)
      dat_filt <- dat_filt |> dplyr::filter(stringr::str_detect(.data$feature_name, feature_incl_filt))
    else
      dat_filt <- dat_filt |> dplyr::filter(.data$feature_name %in% feature_incl_filt)
  }

  if(all(!is.na(feature_excl_filt)) & all(feature_excl_filt != "")){
    if(length(feature_excl_filt) == 1)
      dat_filt <- dat_filt |> dplyr::filter(!stringr::str_detect(.data$feature_name, feature_excl_filt))
    else
      dat_filt <- dat_filt |> dplyr::filter(!.data$feature_name %in% feature_excl_filt)
  }

  if(all(!is.na(qc_types)) & all(qc_types != "")){
    if(length(qc_types) == 1)
      dat_filt <- dat_filt |> dplyr::filter(stringr::str_detect(.data$qc_type, qc_types))
    else
      dat_filt <- dat_filt |> dplyr::filter(.data$qc_type %in% qc_types)
  }

  if(nrow(dat_filt) < 1) stop("None of the feature names meet the filter criteria. Please check feature_filter_include and feature_filter_exclude parameters.")

  # Re-order qc_type levels to define plot layers, i.e.  QCs are plotted over StudySamples
  dat_filt$qc_type <- factor(as.character(dat_filt$qc_type), pkg.env$qc_type_annotation$qc_type_levels)

  # Cap upper outliers  ----

  dat_filt <- dat_filt %>% dplyr::mutate(value =  !!plot_var_s)


  dat_filt <- dat_filt %>%
    dplyr::group_by(.data$feature_name) %>%
    dplyr::mutate(
      #value_max_spl = mean(.data$value[.data$qc_type=="SPL"], na.rm=T) + cap_SPL_SD * sd(.data$value[.data$qc_type=="SPL"], na.rm=T),
      #value_max_qc = mean(.data$value[.data$qc_type==qc_type_fit], na.rm=T) + cap_QC_SD * sd(.data$value[.data$qc_type==qc_type_fit]), na.rm=T,
      value_max_spl = quantile(.data$value[.data$qc_type == "SPL"], 0.75, names = FALSE, na.rm = TRUE) + cap_spl_iqr_factor * IQR(.data$value[.data$qc_type=="SPL"], na.rm = TRUE),
      value_max_qc = quantile(.data$value[.data$qc_type == qc_type_fit], 0.75, names = FALSE, na.rm = TRUE) + cap_qc_iqr_factor * IQR(.data$value[.data$qc_type == qc_type_fit], na.rm = TRUE),
      value_max =  pmax(.data$value_max_spl, .data$value_max_qc, na.rm = TRUE),
      value_mod = dplyr::if_else(cap_outliers & .data$value > .data$value_max,
                                 if(!all(is.na(.data$value)))
                                   max(.data$value[.data$value < .data$value_max], na.rm = TRUE)
                                 else
                                   NA_real_,
                                 .data$value)) %>%
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

  #TODO if(!silent) print(paste0("Plotting ", max(page_range), " pages..."))

  #p_list <- vector("list", length(page_range))
  p_list <- list()
  for (i in page_range){
    if(!silent) print(paste0("page ", i))
    p <- runscatter_one_page(dat_filt = dat_filt, data= data, d_batches = data@annot_batches, cols_page = cols_page, rows_page = rows_page, show_driftcorrection = show_driftcorrection,
                             show_trend_samples, trend_samples_fun, trend_samples_col, after_correction = after_correction, qc_type_fit = qc_type_fit, outputPDF = outputPDF, page_no = i,
                             point_size = point_size, cap_outliers = cap_outliers, point_transparency = point_transparency, annot_scale = annot_scale,
                             show_batches = show_batches, batches_as_shades = batches_as_shades, batch_line_color = batch_line_color, plot_other_qc,
                             batch_shading_color = batch_shading_color, y_label=y_label, base_size=base_size, point_stroke_width=point_stroke_width, show_grid = show_gridlines, log_scale = log_scale)
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
                                outputPDF, annot_scale, point_transparency, point_size=point_size, y_label, base_size, point_stroke_width, show_grid, log_scale){

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
  if (nrow(dat_subset) > 0){
    dMax <- dat_subset %>%
      dplyr::group_by(.data$feature_name) %>%
      dplyr::summarise(y_max =
                         if(!all(is.na(.data$value_mod)))
                           max(.data$value_mod, na.rm = TRUE)*1.0
                       else
                         NA_real_)
  }


  d_batch_data <- d_batches %>% dplyr::slice(rep(1:dplyr::n(), each = nrow(dMax)))
  d_batch_data$feature_name <- rep(dMax$feature_name, times = nrow(d_batches))
  d_batch_data <- d_batch_data %>% dplyr::left_join(dMax, by=c("feature_name"))
  p <- ggplot2::ggplot(dat_subset, ggplot2::aes_string(x="run_id"))

  #browser()
  if (show_batches) {

    if (!batches_as_shades) {
      d_batches_temp <- d_batch_data  |> filter(.data$id_batch_start != 1)
      p <- p + ggplot2::geom_vline(data = d_batches_temp, ggplot2::aes(xintercept = .data$id_batch_start - 0.5), colour = batch_line_color, linetype = "solid", size = .5,  na.rm = TRUE)
    }
    else {
      d_batches_temp <- d_batch_data  %>% dplyr::filter(.data$batch_no %% 2 != 1)
      p <- p + ggplot2::geom_rect(data = d_batches_temp, ggplot2::aes(xmin = .data$id_batch_start - 0.5 , xmax = .data$id_batch_end + 0.5, ymin = 0, ymax = .data$y_max, label = .data$batch_id),
                         inherit.aes = FALSE, fill = batch_shading_color, color = NA, alpha = 0.5, linetype = "solid", size = 0.3,  na.rm = TRUE)
    }
  }

  if(cap_outliers)
    p <- p + ggplot2::geom_hline(data = dMax, ggplot2::aes(yintercept = .data$y_max), color = "#fa9b9b", size = 3, alpha = .4)

  p <- p +
    ggplot2::geom_point(aes_string(x = "run_id", y= "value_mod", color="qc_type", fill="qc_type", shape="qc_type", group="batch_id"), size=point_size, alpha=point_transparency, stroke = point_stroke_width, na.rm = TRUE)



  if(after_correction & show_driftcorrection){
    p <- p +
      ggplot2::geom_line(aes_string(x = "run_id", y= "CURVE_Y_PREDICTED", group = "batch_id"), color = pkg.env$qc_type_annotation$qc_type_fillcol[qc_type_fit], size = .5,  na.rm = TRUE)
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
        ggplot2::geom_smooth(data = filter(dat_subset, .data$qc_type == qc_type_fit), ggplot2::aes_string(x = "run_id", y= "value", group = "batch_id"), se=TRUE,na.rm = TRUE,
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
    ggplot2::xlab("Analysis order") +
    ggplot2::ylab(label = y_label) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0.02,0.03))) +
    #expand_limits(y = 0) +
    ggplot2::theme_bw(base_size = base_size) +

    ggplot2::theme(plot.title = ggplot2::element_text(size=1, face="bold"),
          strip.text = ggplot2::element_text(size=9*annot_scale, face="bold"),
          strip.background = ggplot2::element_rect(size=0.0001,fill="#00283d"),
          strip.text.x = ggplot2::element_text(color = "white"),
          axis.text.x = ggplot2::element_text( size=7*annot_scale, face=NULL),
          axis.text.y = ggplot2::element_text( size=7*annot_scale, face=NULL),
          axis.title = ggplot2::element_text( size=8*annot_scale, face=NULL),
          panel.grid.major =  ggplot2::element_blank(),
          panel.grid.minor =  ggplot2::element_blank(),
          strip.switch.pad.wrap = ggplot2::unit(-1,"mm"),
          panel.border = element_rect(linewidth = 0.5, color = "grey40"))

  if(show_grid)
    p <- p  + ggplot2::theme(panel.grid.major =  ggplot2::element_line(size=0.3,colour = "grey88",linetype = "dashed"))

  if(log_scale)
      p <- p + scale_y_log10()


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

  plot_var_sym = sym(plot_var)

  #https://stackoverflow.com/questions/9843660/marking-the-very-end-of-the-two-whiskers-in-each-boxplot-in-ggplot2-in-r-statist
  get_tails = function(x) {
    q1 = quantile(x,na.rm = TRUE)[2]
    q3 = quantile(x,na.rm = TRUE)[4]
    iqr = q3 -q1
    upper = q3+1.5*iqr
    lower = q1-1.5*iqr
    if(length(x) == 1){return(x)} # will deal with abnormal marks at the periphery of the plot if there is one value only
    ##Trim upper and lower
    up = max(x[x < upper])
    lo = min(x[x > lower])
    return(c(lo, up))
  }



  if(!use_qc_filtered_data)
    d_temp <- data@dataset
  else
    d_temp <- data@dataset_filtered

  d_temp <- d_temp %>%
    dplyr::select(analysis_id, run_id, qc_type, batch_id, feature_name, feature_intensity, feature_norm_intensity, feature_conc) %>%
    filter(feature_intensity > min_feature_intensity) %>%
    filter(str_detect(qc_type, qc_types)) %>%
    droplevels()

  if(relative_log_abundances){
    d_temp <- d_temp %>%
      group_by(feature_name) %>%
      mutate(val = !!plot_var_sym) %>%
      mutate(val = val/ mean(val[qc_type == "BQC"|qc_type == "TQC"|qc_type == "SPL"],na.rm = TRUE))
  } else
  {
    d_temp <- d_temp %>% mutate(val = !!plot_var_sym)
  }

  breaks <- data$dataset %>%
    dplyr::select(run_id) %>% distinct() %>%
    mutate(ticks_to_plot = run_id %% 10 == 0) %>%
    pull(run_id)


  #d_temp$run_id <- as_factor(d_temp$run_id)

  p <- ggplot(d_temp, aes(x=run_id, y=log2(val), group=run_id))

  if(show_batches){
    if (!batches_as_shades){
      p <- p + geom_vline(data = data@annot_batches %>% slice(-1), aes(xintercept=id_batch_start-0.5), colour=batch_line_color, linetype="solid", size=1)
    }
    else {
      d_batch_2nd <- data$batch_info %>% slice(-1) %>% filter(batch_id %% 2 != 1)
      p <- p + geom_rect(data = d_batch_2nd, aes(xmin = id_batch_start-0.5 , xmax = id_batch_end+0.5, ymin = -Inf, ymax = Inf),
                         inherit.aes=FALSE, fill=batch_shading_color, color= NA,alpha= 0.1, linetype="solid", size=0.5, na.rm = TRUE)
    }
  }

  scale_dataset_size = 2^ceiling(log2(max(data$dataset$run_id)))/100

  #geom_point(size=3, color = "#0053a8",alpha=0.6) +
  #stat_summary(aes(x=lipidClass, y=BQC_normIntensity_CV),fun.data="plot.median", geom="errorbar", colour="#fc0000", width=0.8, size=2, inherit.aes=FALSE,na.rm = TRUE) +
  p <- p +
    geom_boxplot(aes(fill = qc_type, color = qc_type), notch=FALSE, outlier.colour = NA, linewidth = 0.2, na.rm = TRUE) +
    #scale_colour_gradient(low = "white", high = "#004489") +
    scale_fill_manual(values = pkg.env$qc_type_annotation$qc_type_col) +
    scale_color_manual(values = pkg.env$qc_type_annotation$qc_type_col) +
    #scale_y_continuous(limits = c(0, 1e6)) +
    #scale_x_continuous(breaks = seq(0, max(data$dataset$run_id)+1, scale_dataset_size)) +
    #  scale_x_discrete(breaks = breaks) +
    theme_bw(base_size = base_size) +
    ylab(bquote(bold(log[2] ~ .(plot_var)))) +
    xlab("Analysis order") +
    theme( panel.grid.major.y = element_blank(),
           panel.grid.minor.y = element_blank(),
           panel.grid.major.x =  element_line(colour="#bdbdbd", linetype="dotted", size=.5),
           panel.grid.minor.x =  element_line(colour="#bdbdbd", linetype="dotted", size=.25),
           axis.text.y = element_text(size = base_size),
           axis.text.x = element_text(size = base_size),
           axis.title = element_text(size = base_size*1, face = "bold"),
           panel.border = element_rect(linewidth =1, color = "grey20"))

  if(ignore_outliers){
    tails <- get_tails(log2(d_temp$val))
    p <- p + scale_y_continuous(limits = c(tails[1]*2,tails[2]*2))
  }

  if(relative_log_abundances){
    p <- p + geom_hline(yintercept = 0, colour="#666666", linetype="dashed", size=0.8) +
      ylab(bquote(bold( 'Rel. ' ~ log[2] ~ .(plot_var))))
  }


  return(p)
}


# Define function to plot 1 page
plot_responsecurves_page <- function(dataset,
                                     output_pdf,
                                     response_variable,
                                     regr_max_percent,
                                     pdf_filename,
                                     rows_page,
                                     columns_page,
                                     point_size,
                                     line_width,
                                     text_scale_factor,
                                     base_size){

  plot_var <- rlang::sym(response_variable)
  ggplot2::ggplot(data = dataset,
                  ggplot2::aes(x = .data$relative_sample_amount ,
             y = !!plot_var,
             color = .data$rqc_series_id)) +
    ggpmisc::stat_poly_line(data = subset(dataset, dataset$relative_sample_amount<= (regr_max_percent/100)),
                            ggplot2::aes(x = .data$relative_sample_amount ,
                       y  = !!plot_var,
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
#' @param use_filt_data Use QC-filtered data
#' @param output_pdf Save as PDF
#' @param response_variable Variable to plot
#' @param feature_incl_filt Filter features names matching the criteria (regex). When empty, `NA` or `NULL` all available features are included.
#' @param feature_excl_filt Exclude features names matching the criteria (regex).  When empty, `NA` or `NULL` all available features are included.
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
                                use_filt_data,
                                output_pdf,
                                response_variable = "feature_intensity",
                                feature_incl_filt = "",
                                feature_excl_filt = "",
                                regr_max_percent = NA,
                                pdf_filename = "",
                                rows_page = 4,
                                columns_page = 5,
                                point_size = 2,
                                line_width = 1,
                                text_scale_factor = 1,
                                return_plot_list = FALSE, base_size = 7) {

  if (output_pdf & pdf_filename == "") stop("Please define parameter `pdf_filename`")

  rows_page = rows_page
  columns_page = columns_page


  if(nrow(data@dataset) < 1) stop("No data available. Please import data and metadata first.")

  if (use_filt_data){
    dat_filt <- data@dataset_filtered %>% dplyr::ungroup()
    if(nrow(dat_filt) < 1) stop("Data has not been qc filtered. Please apply `apply_qc_filter` first.")
  } else {
    dat_filt <- data@dataset %>% dplyr::ungroup()
  }



  if(all(!is.na(feature_incl_filt)) & all(feature_incl_filt != "")){
    if(length(feature_incl_filt) == 1)
      dat_filt <- dat_filt |> dplyr::filter(stringr::str_detect(.data$feature_name, feature_incl_filt))
    else
      dat_filt <- dat_filt |> dplyr::filter(.data$feature_name %in% feature_incl_filt)
  }

  if(all(!is.na(feature_excl_filt)) & all(feature_excl_filt != "")){
    if(length(feature_excl_filt) == 1)
      dat_filt <- dat_filt |> dplyr::filter(!stringr::str_detect(.data$feature_name, feature_excl_filt))
    else
      dat_filt <- dat_filt |> dplyr::filter(!.data$feature_name %in% feature_excl_filt)
  }

  d_rqc <- dat_filt |>
    filter(.data$qc_type == "RQC") |>
    dplyr::select(tidyselect::any_of(
      c("analysis_id", "feature_name", "feature_intensity", "feature_norm_intensity")
    )) |>
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
          output_pdf = output_pdf,
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
  if (!output_pdf) {
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

  if(!c("LipidClass") %in% names(d_QC)) stop("This function currently only works with lipidomics data. Please add lipid names/class with the function `lipidomics_get_lipid_class_names` before calling this function.")
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
plot_pca_qc <- function(data, variable, dim_x, dim_y, log_transform, remove_istds, point_size = 5, point_alpha = 0.8, ellipse_alpha = 0.8, font_base_size = 12) {


  d_wide <- data@dataset_filtered #|> filter(qc_types %in% c("BQC", "TQC", "SPL", "NIST"))

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

  pca_annot$qc_type <-  droplevels(factor(pca_annot$qc_type, levels = c("SPL", "UBLK", "SBLK", "TQC", "BQC", "RQC", "LTR", "NIST", "PBLK")))
  pca_annot <- pca_annot %>%
    dplyr::arrange(.data$qc_type)

  pca_contrib <- pca_res |> broom::tidy(matrix = "eigenvalues")

  p <- ggplot(data = pca_annot, aes_string(paste0(".fittedPC", dim_x),
                                           paste0(".fittedPC", dim_y),
                                           color = "qc_type",
                                           fill = "qc_type",
                                           shape = "qc_type",
                                           group = "qc_type"
  )) +
    ggplot2::geom_hline(yintercept = 0, size = 0.5, color = "grey80", linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, size = 0.5, color = "grey80", linetype = "dashed") +
    suppressWarnings(ggplot2::stat_ellipse(data = pca_annot |> filter(qc_type %in% c("BQC", "TQC", "SPL")),  geom = "polygon", level = 0.95,alpha = ellipse_alpha, size = 0.3, na.rm = TRUE)) +
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

  d_wide = data@dataset_filtered  %>% filter(.data$qc_type %in% c("BQC", "TQC", "NIST", "LTR", "SPL"))

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
    dplyr::rename(feauture_name = .data$column)

  d_loadings_selected <- d_loading |>
    tidyr::pivot_longer(cols = -.data$feauture_name,  names_to = "PC", values_to = "Value") |>
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
    tidyr::unite("Feature", .data$feauture_name, .data$PC, remove = FALSE) |>
    mutate(PC = as.factor(.data$PC),
           Feature = forcats::fct_reorder(.data$Feature, .data$abs_value))


  p <- ggplot(d_loadings_selected, ggplot2::aes(x = .data$Feature, y = .data$Value, color = .data$direction, fill = .data$direction)) +
    ggplot2::geom_col() +
    ggplot2::facet_wrap(ggplot2::vars(.data$PC), scales = "free", ncol = ifelse(vertical_bars, 1, length(pc_dimensions))) +
    ggplot2::scale_x_discrete(labels=d_loadings_selected$feauture_name, breaks=d_loadings_selected$Feature) +
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


#' Plot results of QC filter per analyte class
#' @param data MidarExperiment object
#' @param user_defined_keeper Include user-defined feature inclusion list, even they did not pass QC filtering
#' @param base_size font size of plots
#' @return ggplot2 object
#' @export


plot_qc_summary_classes <- function(data, user_defined_keeper = FALSE, base_size=8) {

  if(user_defined_keeper) stop("user_defined_keeper = TRUE not yet supported")

  d_qc <- data@metrics_qc |>
    filter(valid_integration, !is_istd) |>
    mutate(feature_class = tidyr::replace_na(feature_class, "Undefined"))


  #TODO: cleanup feature/lipidclasses
  #if(!all(is.na(d_qc$feature_class)) & any(is.na(d_qc$lipid_class))) d_qc$feature_class <- d_qc$lipid_class

  if(all(is.na(d_qc$feature_class))) stop("This plot required feature_class to be defined . Please define feature classes in the metadata or retrieve via corresponding {midar} functions.")


  d_qc$feature_class <- forcats::fct(d_qc$feature_class)

  d_QC_sum <- d_qc  %>%
    group_by(feature_class) %>%
    summarise(
      has_only_na = sum(!pass_no_na),
      below_lod = sum(pass_no_na & !pass_lod),
      below_sb = sum(pass_lod & pass_no_na & !pass_sb),
      above_cva = sum(pass_lod & pass_no_na & pass_sb & !pass_cva),
      bad_linearity = sum(pass_lod & pass_no_na & pass_sb & pass_cva & !pass_linearity),
      above_dratio = sum(pass_lod & pass_no_na & pass_sb & pass_cva & pass_linearity & !pass_dratio),
      qc_pass = sum(qc_pass)
    ) %>%
    tidyr::pivot_longer(-feature_class, names_to = "qc_criteria", values_to  = "count_pass") %>%
    ungroup() %>%
    mutate(qc_criteria = factor(qc_criteria, c("below_lod", "has_only_na", "below_sb", "above_cva", "above_dratio",  "bad_linearity", "qc_pass")))



  ggplot(d_QC_sum, aes(forcats::fct_rev(feature_class), count_pass)) +
    geom_bar(aes(fill = qc_criteria), stat="identity", na.rm = TRUE) +
    scale_fill_manual(values=  c(below_lod = "#c7c7c7", has_only_na = "#5e555e", below_sb = "#8f8f8f", above_cva = "#870331", bad_linearity = "#009ec9", above_dratio = "#ffabab", qc_pass = "#02bd62")) +
    #facet_wrap(~Tissue) +
    #guides(fill = guide_legend(override.aes = list(size = 6))) +
    coord_flip() +
    labs(y = "Number of features", x = "Feature class") +
    scale_x_discrete(limits = rev(levels(d_QC_sum$feature_class))) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = c(0.8,0.8),
      panel.grid.major.x = element_line(color = "grey80", linewidth = .1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = .1),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = base_size-1 ), # Legend text size
      legend.title = element_text(size = base_size ), # Legend title size
      legend.key.size = unit(2, "mm")) # Legend key size
}

#' Plot summary of QC filtering
#' @param data MidarExperiment object
#' @param user_defined_keeper Include user-defined feature inclusion list, even they did not pass QC filtering
#' @param base_size font size of plots
#' @return ggplot2 object
#' @export


plot_qc_summary_venn <- function(data, user_defined_keeper, base_size= 12) {

  if(user_defined_keeper) stop("user_defined_keeper = TRUE not yet supported")

  d_qc <- data@metrics_qc |>
    filter(valid_integration, !is_istd) |>
    mutate(feature_class = tidyr::replace_na(feature_class, "Undefined"))


  d_QC_sum_total <- d_qc  %>%
    ungroup() |>
    summarise(
      has_only_na = sum(!pass_no_na),
      below_lod = sum(pass_no_na & !pass_lod),
      below_sb = sum(pass_lod & pass_no_na & !pass_sb),
      above_cva = sum(pass_lod & pass_no_na & pass_sb & !pass_cva),
      bad_linearity = sum(pass_lod & pass_no_na & pass_sb & pass_cva & !pass_linearity),
      above_dratio = sum(pass_lod & pass_no_na & pass_sb & pass_cva & pass_linearity & !pass_dratio),
      qc_pass = sum(qc_pass)
    ) %>%
    tidyr::pivot_longer(names_to = "qc_criteria", values_to  = "count_pass", cols = everything()) %>%
    ungroup() %>%
    mutate(qc_criteria = factor(qc_criteria, c("below_lod", "has_only_na", "below_sb", "above_cva", "above_dratio", "bad_linearity", "qc_pass")))



  # d_QC_sum_total <- d_QC_sum %>%
  #   group_by(qc_criteria) %>%
  #   summarise(Count = sum(Count,na.rm = TRUE)) %>%
  #   mutate(totalCount = sum(Count,na.rm = TRUE),
  #          cumCount = cumsum(Count),
  #          centres = totalCount - (cumCount - Count / 2)) %>% ungroup()

  p_bar <- ggplot(d_QC_sum_total, aes(x=reorder(qc_criteria, count_pass), y=count_pass, fill=qc_criteria))+
    geom_bar(width = 1, stat = "identity") +
    coord_flip() +
    scale_fill_manual(values=  c(below_lod = "#c7c7c7", has_only_na = "#5e555e", below_sb = "#8f8f8f", above_cva = "#870331", bad_linearity = "#009ec9", above_dratio = "#ffabab", qc_pass = "#02bd62")) +
    #geom_text(aes(label = Count), size=4 ) +
    geom_text(
      aes(
        label= count_pass,
        hjust=ifelse(count_pass < max(count_pass) / 1.5, -2, 2) # <- Here lies the magic
      )) +
    labs(x="", "Number of analytes")
    #facet_wrap(~Tissue) +
    theme_bw(base_size = base_size) +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          legend.text = element_text(size = base_size - 2), # Legend text size
          legend.title = element_text(size = base_size - 4), # Legend title size
          legend.key.size = unit(3, "mm")) # Legend key size

  # prevent creating log file

    d_qc_venn <- d_qc

    cva_failed = d_qc_venn$feature_name[!d_qc_venn$pass_cva & d_qc_venn$pass_no_na]
    lod_sb_failed = d_qc_venn$feature_name[(!d_qc_venn$pass_lod | !d_qc_venn$pass_sb) & d_qc_venn$pass_no_na]
    #sb_failed = d_qc_venn$feature_name[!d_qc_venn$pass_sb  & d_qc_venn$pass_no_na]
    lin_failed = d_qc_venn$feature_name[!d_qc_venn$pass_linearity & d_qc_venn$pass_no_na]

    lod_sb_label <- "below S/B or LoD" #paste0('LoD < ', parameter_processing$)
    #sb_label <- "below S/B or LOD" #paste0('S/B < ', parameter_processing$)
    cva_label <- "above CVa" #paste0('CV > ', percent(MAX_CV_NORM/100))
    lin_label <- "below R^2" #paste0('RQC r^2 < ', MIN_LINEARITY_RSQUARE, ' OR rel y0 > ', REL_Y_INTERSECT)o

    x2 = list(lod_sb_failed, cva_failed, lin_failed)
    names(x2) <- c(lod_sb_label, cva_label, lin_label)

    p_venn <- ggvenn::ggvenn(x2, c(lod_sb_label ,cva_label, lin_label),
                     show_percentage = FALSE,
                     fill_color = c("#c7c7c7", "#ff8080", "#009ec9"),
                     fill_alpha = 0.5,
                     stroke_size = 0.0,
                     text_size = 5,
                     set_name_size = 6) # + ggplot2::coord_cartesian(clip="off")



  #plt <- arrangeGrob(p_bar, gTree(NULL),gTree(children=p_venn), ncol=3, widths=c(1.6, 0, 1.3), padding = 0)
  plt <- p_bar +  p_venn + patchwork::plot_layout(ncol = 2, widths = c(1,1)) + patchwork::plot_annotation(tag_levels = c("A", "B"))

  return(plt)
}


