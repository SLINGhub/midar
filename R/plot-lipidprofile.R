
plot_abundanceprofile <- function(data,
                                  variable, 
                                  qc_types,
                                  exclude_classes = NA,
                                  filter_data = FALSE,
                                  include_qualifier = FALSE,
                                  include_istd = FALSE,
                                  include_feature_filter = NA,
                                  exclude_feature_filter = NA,
                                  lipiomics = NA,
                                  feature_cat_color = NA,
                                  scalef,
                                  xlim = NA,
                                  unit,
                                  legend_pos = "right",
                                  d_QC = NULL) {
  

  # Load required packages  # Match the selected variable with predefined options
  variable <- str_remove(variable, "feature_")
  rlang::arg_match(variable, c("area", "height", "intensity", "norm_intensity",
                               "intensity_raw", "intensity_before", "norm_intensity_raw", "norm_intensity_before", "response",
                               "conc", "conc_raw", "conc_before", "rt", "fwhm", "width", "symmetry"))
  variable <- stringr::str_c("feature_", variable)
  variable_sym = rlang::sym(variable)     
  check_var_in_dataset(data@dataset, variable)  
  


  d_filt <- get_dataset_subset(
    data,
    filter_data = filter_data,
    qc_types = qc_types,
    include_qualifier = include_qualifier,
    include_istd = include_istd,
    include_feature_filter = include_feature_filter,
    exclude_feature_filter = exclude_feature_filter
  )


  d_filt <-d_filt |>  
    select("feature_id", "analysis_id", "qc_type", "feature_class", !!variable_sym) |> 
    distinct()
    
  if(data@analysis_type == "lipidomics"){
    d_filt <-d_filt |> 
      mutate(feature_class = factor(.data$feature_class, levels = rev(pkg.env$lipid_class_themes$lipid_class_order))) |> 
      arrange(.data$feature_class)
  }


  feature_classes <- d_filt$feature_class |> unique()
  feature_classes_indices <- setNames(feature_classes, (1:length(feature_classes)))  

  d_filt$feature_class_index <- d_filt$feature_class
  levels(dat$feature_class_index) <- as.list(feature_classes_indices)
    
  d_filt <- d_filt |>
    mutate(!!variable_sym := !!variable_sym * scalef)


  if(!is.na(exclude_classes)){
    d_filt <- d_filt |>
      filter(
        !(.data$class_index %in% exclude_classes
        )
      )
  }

  d_filt <- d_filt |> 
    group_by(.data$feature_id, .data$feature_class_index, .data$feature_class) |> 
    summarise(abundance_mean = mean(!!variable_sym, na.rm = TRUE))

  d_filt_minmax <- d_filt |>
    group_by(.data$class_index, .data$feature_class) |>
    summarise(
      abundance_mim = min(.data$abundance_mean, na.rm = TRUE),
      abundance_max = max(.data$abundance_mean, na.rm = TRUE)
    )

  d_sum <- d_filt |>
    group_by(.data$analysis_id, .data$feature_class_index, .data$feature_class) |>
    summarise(
      abundance_sum = sum(!!variable_sym, na.rm = T)
    ) |>
    group_by(.data$class_name) |>
    summarise(
      abundance_sum_mean = mean(.data$abundance_sum, na.rm = TRUE),
      CV = sd(.data$abundance_sum, na.rm = TRUE) / .data$abundance_sum_mean * 100
    ) |>
    left_join(d_filt_minmax)


  breaks <- 10^(xlim[1]:xlim[2])
  minor_breaks <- rep(1:10, 1) * (10^rep(xlim[1]:xlim[2], each = 9))

  #  print(1:length(levels(d_temp$class_index)))

  d_blank <- tibble(x = 1:65, y = 0)

  browser()
  plt <- ggplot(
    d_filt,
    aes(
      x = .data$abundance_mean,
      y = as.numeric(as.character(.data$feature_class_index))
    )
  ) +
    ggplot2::scale_y_continuous(
      position = legend_pos,
      limits = c(0.0, length(levels(d_class$feature_class_index)) + 1),
      breaks = 1:length(levels(d_class$feature_class_index)),
      labels = feature_classes,
      expand = c(0.002, 0.002)
    ) +
    xlab(unit) +
    ylab("") +
    ggplot2::scale_x_log10(
      limits = c(10^xlim[1], 10^xlim[2]),
      breaks = breaks,
      minor_breaks = minor_breaks,
      labels = scales::trans_format("log10", function(x) scales::math_format()(log10(x)))
    ) +
    ggplot2::annotation_logticks(
      base = 10,
      sides = "b",
      size = 0.3,
      long = unit(1, "mm"),
      short = unit(0, "mm"),
      mid = unit(.5, "mm"),
      colour = "grey80"
    ) +
    geom_rect(
      data = d_class,
      aes(
        xmin = .data$abundance_min * 0.8,
        xmax = .data$abundance_max * 1.2,
        ymin = as.numeric(as.character(.data$feature_class_index)) - 0.4,
        ymax = as.numeric(as.character(.data$feature_class_index)) + 0.4,
        group = .data$feature_class,
        fill = .data$feature_class
      ),
      inherit.aes = FALSE,
      show.legend = FALSE,
      alpha = 0.5
    ) +

    # does not work (yet) with non-linear scales
    #geom_rect_pattern(data = d_temp_class,
    #                  aes(xmin = LOD*0.8 , xmax = LOD, ymin =  as.numeric(lipidClass)-0.4, ymax = as.numeric(lipidClass)+0.4, group = lipidClass,  pattern_fill2 = feature_cat),
    #                  inherit.aes = FALSE, show.legend = F, pattern = "gradient",  pattern_orientation = "horizontal", pattern_density = .3, colour = NA) +
    geom_segment(
      aes(
        x = .data$conc_mean,
        y = as.numeric(as.character(.data$feature_class_index)) - 0.3,
        xend = .data$conc_mean,
        yend = as.numeric(as.character(.data$feature_class_index)) + 0.3
      ),
      color = "black",
      size = 0.25,
      show.legend = T,
      inherit.aes = TRUE
    ) +

    geom_point(
      data = d_class,
      aes(
        x = .data$conc_sum_mean,
        y = as.numeric(as.character(.data$class_index)),
        color = .data$feature_class
      ),
      size = 1.3,
      shape = 23,
      stroke = 0.8,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = feature_cat_color, drop = FALSE) +
    scale_color_manual(values = feature_cat_color, drop = FALSE) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor.x = element_line(
        color = "grey68",
        linetype = "dotted",
        size = .01
      ),
      panel.grid.major.x = element_line(
        color = "grey70",
        linetype = "solid",
        size = .03
      ),
      axis.text.y = element_text(size = 7),
      axis.text.x = element_text(size = 7),
      axis.title.x = element_text(size = 8, face = "bold"),
      axis.title.y = element_text(size = 8, face = "bold"),
      axis.ticks = element_line(colour = "grey40", size = .5),
      legend.position = "none",
      panel.border = element_rect(color = "grey40", size = 0.5),
      plot.margin = unit(c(0.75, 0.10, 0.0, 0.1), "cm"),
      plot.title = element_text(size = 8, hjust = 0.5)
    ) +
    ggplot2::guides(fill = FALSE)

  plt
  }