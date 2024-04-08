

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


