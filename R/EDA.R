
#' Get FC and p values from a t test contrast
#'
#' @param d_FC  data.frame with FEATURE_NAME, p and log2FC
#' @param p_adjust calculate and use FDR
#' @param sig_FC_min |FC| treshold
#' @param sig_p_value_min P value treshold
#' @param symmetric_x make x axis symetric around 0
#' @param x_min x min
#' @param x_max x max
#' @param point_size point size
#' @param point_transparency point alpha
#' @param scale_factor scale factor for fonts
#' @param scale_factor_species_label scale factor point labels
#' @param squared make plot a square
#' @param axistick axis ticks
#' @return tibble with p values, fdr and log2FC
#' @export
volcano_plot <- function(d_FC, p_adjust, sig_FC_min, sig_p_value_min, symmetric_x, x_min, x_max, point_size, point_transparency, scale_factor, scale_factor_species_label, squared = TRUE,  axistick =1){

  if (p_adjust){
    d_FC$p_value_mod <- p.adjust(d_FC$p_value, method = "BH")
    y_lab = " FDR-adj. "}
  else {
    d_FC$p_value_mod <- d_FC$p_value
    y_lab = "Unadj. "}

  d_FC$Significant <- ifelse((d_FC$p_value_mod < sig_p_value_min) & (abs(d_FC$log2FC) > log2(sig_FC_min)), "sign", "ns")

  FC_max_positive <- abs(max(d_FC$log2FC))
  FC_max_negative <- abs(min(d_FC$log2FC))
  FC_abs_max <- max(FC_max_negative, FC_max_positive)

  FC_seq <- seq(from = -(floor(FC_abs_max)), to=floor(FC_abs_max), by=1)

  y_max <- -log10(min(d_FC$p_value_mod))
  y_break_max = ceiling(y_max*1.2)
  y_breaks <- seq(0, y_break_max)
  y_labels <- 10^(-seq(0, y_break_max))


  p <- ggplot(d_FC, aes(x = .data$log2FC, y = -log10(.data$p_value_mod))) +
    #ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
    #geom_hline(yintercept = -log10(0.001), colour="#990000", linetype="dashed", size=0.1) +
    #geom_hline(yintercept = -log10(0.01), colour="#990000", linetype="dashed", size=0.1) +
    geom_hline(yintercept = -log10(0.05), colour="#d99898", linetype="dashed", size=0.3) +
    geom_vline(xintercept = -log2(sig_FC_min), colour="#d99898", linetype="dashed", size=0.3) +
    geom_vline(xintercept = log2(sig_FC_min), colour="#d99898", linetype="dashed", size=0.3) +
    geom_vline(xintercept = 0, colour="#8fb0c4", size=.3) +
    geom_point(size = .data$point_size, alpha = .data$point_transparency, na.rm = T, aes(colour = .data$Significant)) +
    scale_color_manual(values = c("grey", "red")) +
    theme_bw(base_size = 14) +
    labs(
      y = bquote(.(y_lab) ~ italic('P') ~ ' value (log scale)'),
      x = bquote(log[2] ~ 'FC')
    ) +
    scale_y_continuous(
      limits = c(0, y_max *1.2),
      breaks = y_breaks,
      labels = y_labels) +
   #scale_y_continuous(
  #    breaks = c()
    #   breaks = c()
    #   #trans  = scales::compose_trans("log10", "reverse"),
    #   #labels = scales::label_log(base = 10),
    #   #breaks  = log_breaks_125
    # ) +
    #ggplot2::annotation_logticks(sides = 'tb', outside = TRUE, size = .2) +
    #scale_y_continuous(trans = "log1p", limits=c(0,5), breaks = c(1, 2, 3,4,5),  expand = c(0,0)) +
    ggrepel::geom_text_repel(data = d_FC %>% filter(.data$Significant =="sign"),
                    aes(label = .data$FEATURE_NAME),
                    size = 3 * scale_factor_species_label,
                    box.padding = unit(0.15, "lines"),
                    point.padding = unit(0.1, "lines")) +
    theme(axis.text = element_text(size=8 * scale_factor),
          axis.title = element_text(size=9 * scale_factor,face = "bold"),
          plot.title = element_text(size=12, face = "bold")) +
    theme(legend.position="none")

  if(symmetric_x) {
    FC_seq <- seq(from = -(ceiling(FC_abs_max)), to=ceiling(FC_abs_max), by=axistick)
    #print(FC_seq)
    p <- p + scale_x_continuous( limits=c(-(FC_abs_max*1.1), FC_abs_max*1.1), breaks = FC_seq,  expand = c(0.02,0.02))
  } else
  {
    min_x <- ifelse(!is.na(x_min), x_min,-(FC_max_negative))
    max_x <- ifelse(!is.na(x_max), x_max,(FC_max_positive))
    FC_seq <- seq(from = min_x, to=max_x, by=axistick)
    #print(FC_seq)
    p <- p + scale_x_continuous(limits=c(min_x*1.1, max_x*1.1), breaks = FC_seq,  expand = c(0.02,0.02), oob = scales::oob_squish)
  }

  p <- p + theme(
    panel.grid.major =  element_line(size=0.3),
    panel.grid.minor =  element_line(size=0.15))

  if (squared) p <- p + theme(aspect.ratio=1)

  p
}

get_CompoundName <- function(transition_name){
  compound_name <- str_trim(str_replace(transition_name, "\\[.*?\\]",""))
  compound_name <- str_trim(str_replace(compound_name, "  "," "))
  compound_name <- str_trim(str_replace(compound_name, "\\+",""))
  compound_name <- str_trim(str_replace(compound_name, "\\[.*?\\]",""))
  compound_name <- str_trim(str_replace(compound_name, "  "," "))
  return(compound_name)
}

add_lipid_transition_classnames <- function(datLong){
  #cat("Retrieving lipid class/transition names...", fill = FALSE)
  datLong_temp <- datLong %>%  mutate(lipidClassBase = (str_trim(str_extract(.data$Compound, "[A-z0-9]+[[:blank:]]*"))))
  datLong_temp <- datLong_temp %>%  mutate(lipidClass = (str_trim(str_extract(.data$Compound, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m])"))))

  # add a "-", except for between name and sphingoid base
  datLong_temp <- datLong_temp %>% mutate(lipidClass = str_replace(.data$lipidClass, "([[:blank:]]+)([^d|t|m]{1})", "-\\2"))
  datLong_temp <- datLong_temp %>% mutate(lipidClassSL = str_extract(.data$Compound, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m][0-9\\:]{4})"))
  # Add  transition name  (defined by flanking "[]")
  datLong_temp <- datLong_temp %>%  mutate(transitionName = str_extract(.data$Compound, "(?<=\\[).*?(?=\\])"))
  datLong_temp <- datLong_temp %>%  mutate(transitionName = ifelse(is.na(.data$transitionName),"",.data$transitionName))

  # Add  transition class (adds measured product class in square brackets after the compound class name). In case of NL indicated with [-...] (e.g.  TG 48:2 [-16:1]), [NL] is indicated (TG [NL]). Maybe useful for normalization
  datLong_temp <- datLong_temp %>%  mutate(transitionClass =  ifelse(grepl("\\[\\-|\\[NL", .data$Compound),"M>M-NL",str_replace(paste("", str_trim(str_extract(.data$Compound, "\\[[a-zA-Z0-9\\-\\: ]*\\]")))," NA", .data$transitionName)))

  datLong_temp <- datLong_temp %>%  mutate(CompoundName = get_CompoundName(.data$Compound))

  return(datLong_temp %>% ungroup())
}




library(ComplexHeatmap)

plot_heatmap <- function(data, d_metadata, annot_color, log_transform, split_variable = NULL, clust_rows = TRUE, clust_columns = FALSE, row_names_size = 8, col_names_size = 8) {

  annot_col = ComplexHeatmap::columnAnnotation(
    df = d_metadata |> dplyr::select(!.data$ANALYSIS_ID) |> as.data.frame(),
    col = annot_color
  )

 d_wide <- data %>%
    pivot_wider(names_from = "FEATURE_NAME", values_from = "Concentration")

  d_filt = left_join(d_metadata[,"ANALYSIS_ID"], d_wide )

  m_raw <- d_filt %>%
    #dplyr::select(where(~!any(is.na(.)))) |>
    column_to_rownames("ANALYSIS_ID") |>
    as.matrix()

  if(log_transform) m_raw <- log2(m_raw)
  m_scaled <- scale(m_raw, center = TRUE, scale = TRUE)


  # cluster samples and compoundss
  # d1 <- dist(t(m), method = "manhattan", diag = TRUE, upper = TRUE)
  # d2 <- dist(m,method = "maximum", diag = TRUE, upper = FALSE)
  # c1 <- hclust(d1, method = "ward.D2", members = NULL)
  # c2 <- hclust(d2, method = "ward.D2", members = NULL)

  z_score_seq <-
    seq(ceiling(min(t(m_scaled),na.rm = T)), ceiling(max(t(m_scaled),na.rm = T)), by = 1)

  lgd = Legend(title = "Row z score", at = z_score_seq)



  lipidCat_color <- c(
    "AcylC" = "#fff263",
    "CAR" = "#fff263",
    "CE" = "#87044e",
    "Cer1P d" = "#bdb782",
    "Cer d" = "#fa9c0f",
    "Cer m" = "#fa9c0f",
    "COH" = "#cf98d6",
    "DG" = "#fab6b6",
    "GM3 d" = "#bae012",
    "Hex1Cer d" = "#53e012",
    "Hex2Cer d" = "#3a9410",
    "Hex3Cer d" = "#3a9410",
    "Hex4Cer d" = "#3a9410",
    "C1P" = "#291b47",
    "LPC" = "#06c9b9",
    "LPC-O" = "#d4eb8f",
    "LPC-P" = "#d4eb8f",
    "LPE" = "#017a70",
    "LPE-P" = "#d4eb8f",
    "PC" = "#0493cc",
    "PC-O" = "yellow",
    "PC-P" = "yellow",
    "PA" = "#291b47",
    "PE" = "#015170" ,
    "PE-O" = "yellow" ,
    "PE-P" = "yellow" ,
    "PG" = "#626682",
    "PI" = "#1d00ff",
    "LPI" = "#5ebdf7",
    "PS" = "#bfc7ff",
    "SM" = "#a41db3",
    "Sph d" = "#946599",
    "TG" = "#99022a",
    "TG-O" = "#f03768"
  )

  df_compounds <- data.frame(Compound = colnames(m_scaled))  %>% add_lipid_transition_classnames(.data$.)

  ha2 = rowAnnotation(
    df = data.frame(lipidClass = df_compounds$lipidClass),
    col = list(lipidClass = lipidCat_color),
    gap = unit(10, "mm"),
    width = unit(10, "mm"),show_annotation_name = FALSE
  )

  Heatmap(
    t(m_scaled),
    row_dend_side = "left",
    cluster_columns = clust_columns,
    cluster_rows = clust_rows,
    #column_split = d_metadata[[split_variable]],
    cluster_column_slices = FALSE,
    column_gap = unit(3, "mm"),
    # reverse order (rev) makes the heatmap look as created with heatmap.2
    color_space = "RGB",
    top_annotation = annot_col,
    left_annotation = ha2,
    row_dend_reorder = TRUE,
    column_dend_reorder = TRUE,
    #column_dend_height = unit(30, "mm"),
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = row_names_size),
    column_names_gp =  gpar(fontsize = col_names_size),
    heatmap_legend_param = list(title = "Row z score", color_bar = "discrete", at = z_score_seq)
  )
}


plot_pca_sling2 <- function(data, d_metadata, annot_color = NULL,log_transform, dim_x, dim_y, grouping, point_size = 2, fill_alpha = 0.1, ellipse = TRUE, mark_ellipse = FALSE, show_labels = FALSE, label_size = 6, max_label_overlaps = Inf) {

  d_wide <- data %>%
    pivot_wider(names_from = "FEATURE_NAME", values_from = "Concentration")

  d_filt = left_join(d_metadata[,"ANALYSIS_ID"], d_wide )


  if(!all(d_filt |> pull(.data$ANALYSIS_ID) == d_metadata |> pull(.data$ANALYSIS_ID))) stop("Data and Metadata missmatch")

  m_raw <- d_filt  |>
    dplyr::select(where(~!any(is.na(.)))) |>
    column_to_rownames("ANALYSIS_ID") |>
    as.matrix()

  if(log_transform) m_raw <- log2(m_raw)


  # get pca result with annotation
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)
  pca_annot <- pca_res |> broom::augment(d_metadata)
  pca_contrib <- pca_res |> broom::tidy(matrix = "eigenvalues")

  p <- ggplot(data = pca_annot, aes_string(paste0(".fittedPC", dim_x),
                                           paste0(".fittedPC", dim_y),
                                           color = grouping,

                                           fill = grouping, label = "ANALYSIS_ID"
  )) +
    geom_hline(yintercept = 0, size = 0.4, color = "grey80") +
    geom_vline(xintercept = 0, size = 0.4, color = "grey80")

    if (ellipse){
      if (!mark_ellipse){
        p <- p +stat_ellipse(geom = "polygon", level = 0.95,alpha = fill_alpha, size = 0.3)
      } else {
        p <- p + ggforce::geom_mark_ellipse(aes(fill = !!rlang::ensym(grouping), color = !!rlang::ensym(grouping),  label = "" ), alpha = fill_alpha, label.colour = c("white", "white"), con.size = 0, linewidth  = 0.01)
      }
    }
  p <- p +  geom_point(size = point_size, shape = 1, stroke = 1.5)

  if(!is.null(annot_color)){
    p <- p +
      scale_color_manual(values = annot_color) +
      scale_fill_manual(values = annot_color)
  }

  if(show_labels){
    p <- p +
      ggrepel::geom_text_repel(
        size = label_size,
        max.overlaps = max_label_overlaps,
        min.segment.length = unit(2, "mm"), seed = 42, box.padding = 0.5, segment.linewidth = .2,segment.size = .2,segment.alpha = .5


      )
  }

  p <- p +
    theme_light(base_size = 8) +
    xlab(glue::glue("PC{dim_x} ({round(pca_contrib[[dim_x,'percent']]*100,1)}%)"))+
    ylab(glue::glue("PC{dim_y} ({round(pca_contrib[[dim_y,'percent']]*100,1)}%)"))+
    theme(
      panel.grid = element_line(size = 0.3, color = "grey95"),
      panel.border = element_rect(size = 1, color = "grey70"),
      aspect.ratio=1)

  p
}

# ToDo: check error behaviour (returns P = 1 currently)
mod_t_test <- function(...) {
  obj<-try(t.test(...), silent=TRUE)
  if (is(obj, "try-error")) return(t.test(x = c(1,2,3), y= c(1,2,3))) else return(obj)
}

get_pval <- function(df, test_col, contrasts, paired) {
  df$`_group` <- df |> pull({{test_col}})
  map_dfr(.x = contrasts,
          .f = ~ broom::tidy(
            mod_t_test(formula = Concentration ~ factor(`_group`), paired = paired, data = subset(df, `_group` %in% .x))) |>
            mutate(id = paste0(.x, collapse = "-"),
                   grp1 = .x[1],
                   grp2 = .x[2],
                   y_max_1 = mean(df$Concentration[df$`_group` == .x[1]])+sd(df$Concentration[df$`_group` == .x[1]]),
                   y_max_2 = mean(df$Concentration[df$`_group` == .x[2]])+sd(df$Concentration[df$`_group` == .x[2]]),
                   y_max = max(y_max_1, y_max_2)))
}

plot_dotboxplus_onepage <- function(data, outer_inner, GroupOuter=NULL, GroupInner, n_row, n_col, scale_colors, contrasts, scale_text,scale_signf,  paired, point_size){

  if(outer_inner) {
    group_levels <- levels(data |> pull({{GroupInner}}))
    n_groups <- length(group_levels)
    print(group_levels)
    d_stat <- data |>
      arrange({{GroupInner}}) |>
      group_by({{GroupOuter}}, .data$FEATURE_NAME) |>
      nest() |>
      mutate(res = map(.x = data, .f = \(x) res = get_pval(x, {{GroupInner}}, contrasts = contrasts, paired=paired))) |>
      unnest(-data) |>
      mutate(
        p_val_sym = case_when(
          p.value > 0.05 ~ "", p.value > 0.01 ~ "*",
          p.value > 0.001 ~ "**", TRUE ~ "***"))|>
      mutate(y_max = max(.data$y_max)) |>
      group_by(.data$FEATURE_NAME) |>
      mutate(y_max_all = max(.data$y_max)) |>
      group_by({{GroupOuter}}) |>
      mutate(group_outer_id = cur_group_id()) |>
      group_by(.data$FEATURE_NAME, {{GroupOuter}}) |>
      mutate(
        x_min = (.data$group_outer_id - 0.4) + 0.8/n_groups/2 + 0.8*1/(n_groups) * (match(.data$grp1, group_levels)-1),
        x_max = (.data$group_outer_id - 0.4) + 0.8/n_groups/2 + 0.8*1/(n_groups) * (match(.data$grp2, group_levels)-1),
        #x_min = (group_outer_id - 0.8/n_groups) + 0.8/1 * (match(grp1, group_levels)-1),
        #x_max = (group_outer_id - 0.8/n_groups) + 0.8/1 * (match(grp2, group_levels)-1),
        x_min_U = match(.data$grp1, group_levels)-1,
        x_max_U = match(.data$grp2, group_levels)-1,
        y_max = .data$y_max_all + (.data$y_max_all/15 * row_number() )
      ) |>
      drop_na(.data$y_max)

    d_stat$group <- 1:nrow(d_stat)
    pos <- position_jitter(width = 0.3, seed = 2)

    plt <- ggplot(data, aes(x = {{GroupOuter}}, y = .data$Concentration, group = .data$SUBJECT_ID)) +
      ggtitle(label = paste(unique(data$lipidClass), collapse = " | " ))

    plt <-  plt +
      geom_rect(xmin=2.5, xmax=6.7, ymin = -Inf, ymax=Inf, fill = "grey94") +
      geom_boxplot(width=0.7, aes( color = {{GroupInner}},fill = {{GroupInner}}),
                   alpha = 0.1, lwd=0.03, outlier.shape = NA,
                   position = position_dodge(width=0.8)) +
      geom_point(aes( color = {{GroupInner}},fill = {{GroupInner}}),size = point_size,position = pos, width = .4, shape = 1, stroke = 0.7) +
      facet_wrap2(vars(.data$FEATURE_NAME), scales = "free", nrow = n_row, ncol = n_col, drop = FALSE, trim_blank = FALSE) +
      scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, .15)))+
      scale_color_manual(values = scale_colors) +
      scale_fill_manual(values = scale_colors) +
      ggsignif::geom_signif(
        mapping = aes(y_position = .data$y_max,
                      xmin = .data$x_min,
                      xmax = .data$x_max,
                      annotations = .data$p_val_sym,
                      group=.data$group),
        data = d_stat,
        manual = TRUE,
        inherit.aes = FALSE,
        margin_top = 0.05,
        step_increase = 0.1,
        size = 0.2 * scale_signf,
        vjust = 0.5,
        tip_length = 0.005,
        textsize = 5,
        color = "grey60") +
      theme_bw(base_size = 10 * scale_text) +
      theme(
        axis.text.x = element_text( size=10 , angle = 45,vjust = 1, hjust = 1),
        axis.text.y = element_text( size=7),
        panel.grid = element_blank(),
        legend.position="right",
        plot.title = element_text(size=10, face = "bold")
      )
  } else {

    #pos <- ggbeeswarm::position_beeswarm(dodge.width = 0.08, method = "center")

    sig_annot = function(x){
      if(x < 0.001) {"***"}
      else if(x < 0.01){"**"}
      else if(x < 0.05){"*"}
      else{NA}}
    pos <- position_jitter(width = 0.3, seed = 2)
    plt <- ggplot(data, aes(x = {{GroupInner}}, y = .data$Concentration)) +
      ggtitle(label = paste(unique(data$lipidClass), collapse = " | " ))

    if (paired) plt <-  plt + geom_line(aes(group = .data$SUBJECT_ID), linewidth = 0.15, color = "grey50", alpha = 0.3)


    plt <-  plt +
      #geom_rect(xmin=2.5, xmax=6.7, ymin = -Inf, ymax=Inf, fill = "grey94") +
      geom_boxplot(width=0.7, aes( color = {{GroupInner}},fill = {{GroupInner}}),
                   alpha = 0.2, lwd=0.1, outlier.shape = NA) +
      ggbeeswarm::geom_quasirandom(aes( color = {{GroupInner}},fill = {{GroupInner}}), size = point_size,na.rm = TRUE, shape = 1, stroke = 0.4,width = 0.2) +
      facet_wrap2(vars(.data$FEATURE_NAME), scales = "free", nrow = n_row, ncol = n_col, drop = FALSE, trim_blank = FALSE) +
      scale_y_continuous(limits = c(0,NA), expand = expansion(mult = c(0, .25)))+
      scale_color_manual(values = scale_colors) +
      scale_fill_manual(values = scale_colors) +
      ggsignif::geom_signif(
        #mapping = aes(y_position = y_max,
        #              xmin = x_min,
        #              xmax = x_max,
        #              annotations = p_val_sym,
        #              group=group),
        #data = d_stat,
        #manual = TRUE,
        #inherit.aes = FALSE,
        comparisons = contrasts,
        test = "t.test",
        test.args = list(paired = paired, alternative = "two.sided", var.equal = FALSE),
        margin_top = 0.05,
        step_increase = 0.2,
        color = "darkblue",
        size = 0.2 *scale_signf,
        vjust = 0.3,
        tip_length = 0.005,
        textsize = 5,
        map_signif_level = sig_annot,
        color = "grey60") +
      theme_bw(base_size = 10* scale_text) +
      theme(
        axis.text.x = element_text( size=10 *  scale_text, angle = 45,vjust = 1, hjust = 1),
        axis.text.y = element_text( size= 7* scale_text),
        panel.grid = element_blank(),
        strip.text = element_text(size = 10* scale_text),
        legend.position="right",
        plot.title = element_text(size=12 * scale_text, face = "bold")
      )
  }
  plt
}


plot_dotboxplus <- function(data, d_metadata, inner_group, outer_group, contrasts, paired, inner_group_colors = NULL,
                            rows_page,columns_page, class_per_page = TRUE, save_pdf = FALSE, filename = NULL, scale_text =1, scale_signf = 1, point_size = 1) {

  d_long_full <- data |>
    right_join(d_metadata |> rename(ANALYSIS_ID = .data$SAMPLE_ID))

  d_lipidnames <- d_long_full |>
    dplyr::select(Compound = .data$FEATURE_NAME) |> distinct() %>%
    add_lipid_transition_classnames(.data$.) |>
    dplyr::select(FEATURE_NAME = "Compound", "lipidClass", "lipidClassSL")

  df <- d_long_full |>
    left_join(d_lipidnames)


  page_list <- df |>
    group_by(.data$lipidClass) |>
    mutate(lipidclass_id = cur_group_id()) |>
    group_by(.data$lipidClass, .data$ANALYSIS_ID) |>
    mutate(lipid_class_page_no =  max(ceiling(which(unique(.data$FEATURE_NAME)==.data$FEATURE_NAME)/ (rows_page * columns_page)))) |>
    group_by(.data$FEATURE_NAME) |>
    mutate(page_no =  ceiling(cur_group_id()/ (rows_page * columns_page)))

  if (class_per_page) {
    plt_list <- page_list |>
      group_by(.data$lipidclass_id, .data$lipid_class_page_no) %>%
      nest() %>%
      mutate(plt = map(data, ~ plot_dotboxplus_onepage(.,outer_inner = FALSE, GroupOuter = {{outer_group}},GroupInner = {{inner_group}},
                                                       n_row = rows_page, n_col = columns_page, scale_colors = inner_group_colors, contrasts = contrasts, paired=paired, scale_text, scale_signf, point_size=point_size)))
  } else {
    plt_list <- page_list |>
      group_by(.data$page_no) %>%
      nest() %>%
      mutate(plt = map(data, ~ plot_dotboxplus_onepage(.,outer_inner = FALSE, GroupOuter = {{outer_group}},GroupInner = {{inner_group}},
                                                       n_row = rows_page, n_col = columns_page, scale_colors = inner_group_colors, contrasts = contrasts, paired=paired, scale_text, scale_signf, point_size=point_size)))
  }


  if(save_pdf){
    pdf(file = filename, onefile = TRUE, paper = "A4r", width = 10, height = 8)
    print(plt_list$plt)
    dev.off()
  } else {
    plt_list$plt
  }
  #plt_list
}




