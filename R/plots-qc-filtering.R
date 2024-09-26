#' Plot results of QC filter per analyte class
#' @param data MidarExperiment object
#' @param use_batches How batches should be used, either across batches (ignoring batches), plot batch individually, or summarize batches (median)
#' @param user_defined_keeper Include user-defined feature inclusion list, even they did not pass QC filtering
#' @param base_size font size of plots
#' @return ggplot2 object
#' @export

# TODO: handling of features with (many) missing values, in SPL, in QC
plot_qc_summary_classes <- function(data, use_batches = c("across", "individual", "summarise"), user_defined_keeper = FALSE, base_size = 8) {
  if (user_defined_keeper) cli::cli_abort("user_defined_keeper = TRUE not yet supported")

  rlang::arg_match(use_batches)

  if(!"pass_missingval" %in% names(data@metrics_qc))
    cli_abort(col_red("QC filter has not yet been applied. Please use `apply_qc_filter()` to filter the data."))


  if(use_batches != "summarise") stop("Currently only `summarise` supported for parameter `batches`")

  d_qc <- data@metrics_qc |>
    filter(.data$valid_feature, !.data$is_istd) |>
    mutate(feature_class = tidyr::replace_na(.data$feature_class, "Undefined"))


  # TODO: cleanup feature/lipidclasses
  # if(!all(is.na(d_qc$feature_class)) & any(is.na(d_qc$lipid_class))) d_qc$feature_class <- d_qc$lipid_class

  if (all(is.na(d_qc$feature_class))) cli::cli_abort("This plot requires `feature_class` to be defined. Please define classes in the feature metadata or retrieve via corresponding {midar} functions.")

  d_qc$feature_class <- forcats::fct(d_qc$feature_class)

# Count how many features failed qc criteria, excluding features that failed before tested criteria (lower hiarchy)
# TODO: can surely be better implemented
  d_qc_sum <- d_qc |> ungroup() |>
    group_by(.data$feature_class) |>
    summarise(
      has_only_na = sum(.data$na_in_all_spl, na.rm = TRUE),
      exceed_missingness = sum((!replace_na(.data$na_in_all_spl, TRUE) & !.data$pass_missingval) | all(is.na(.data$pass_missingval)), na.rm = TRUE),
      below_lod = sum((!.data$na_in_all_spl & replace_na(.data$pass_missingval, TRUE)) & !replace_na(.data$pass_lod, TRUE) | all(is.na(.data$pass_lod)), na.rm = TRUE),
      below_sb = sum((!.data$na_in_all_spl & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & !.data$pass_sb) | all(is.na(.data$pass_sb)), na.rm = TRUE),
      above_cva = sum((!.data$na_in_all_spl & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & replace_na(.data$pass_sb, TRUE) & !.data$pass_cva) | all(is.na(.data$pass_cva)), na.rm = TRUE),
      bad_linearity = sum((!.data$na_in_all_spl & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & replace_na(.data$pass_sb, TRUE) & replace_na(.data$pass_cva, TRUE) & !.data$pass_linearity) | all(is.na(.data$pass_linearity)), na.rm = TRUE),
      above_dratio = sum((!.data$na_in_all_spl & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & replace_na(.data$pass_sb, TRUE) & replace_na(.data$pass_cva, TRUE) & replace_na(.data$pass_linearity, TRUE) & !.data$pass_dratio) | all(is.na(.data$pass_dratio)), na.rm = TRUE),
      qc_pass = sum(.data$qc_pass, na.rm = TRUE)
    ) |>
    tidyr::pivot_longer(-.data$feature_class, names_to = "qc_criteria", values_to = "count_pass") |>
    ungroup()

  qc_colors <- c(qc_pass = "#4bcc7f", above_dratio = "#ffcf57", above_cva = "#bf0426", bad_linearity = "#1fc1ed", below_sb = "#cccaca", below_lod = "#919191", exceed_missingness = "yellow", has_only_na = "#111111")
  d_qc_sum$qc_criteria <- forcats::fct_relevel(d_qc_sum$qc_criteria, rev(names(qc_colors)))


  # Remove levels/qc criteria for which was not filtered for
  if(all(is.na(d_qc$na_in_all_spl))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "has_only_na")
  if(all(is.na(d_qc$pass_missingval))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "exceed_missingness")
  if(all(is.na(d_qc$pass_lod))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "below_lod")
  if(all(is.na(d_qc$pass_sb))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "below_sb")
  if(all(is.na(d_qc$pass_cva))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "above_cva")
  if(all(is.na(d_qc$pass_linearity))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "bad_linearity")
  if(all(is.na(d_qc$pass_dratio))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "above_dratio")

  d_qc_sum <- d_qc_sum |> drop_na(.data$qc_criteria)

  ggplot(d_qc_sum, aes(forcats::fct_rev(.data$feature_class), .data$count_pass)) +
    ggplot2::geom_bar(aes(fill = .data$qc_criteria), stat = "identity", na.rm = TRUE) +
    scale_fill_manual(values = qc_colors, na.value = "purple", drop = FALSE) +
    # facet_wrap(~Tissue) +
    # guides(fill = guide_legend(override.aes = list(size = 6))) +
    ggplot2::coord_flip() +
    labs(y = "Number of features", x = "Feature class") +
    #facet_wrap(~batch_id, scales = "free") +
    ggplot2::scale_x_discrete(limits = rev(levels(d_qc_sum$feature_class))) +
    theme_bw(base_size = base_size) +
    theme(
      legend.position = c(0.8, 0.8),
      panel.grid.major.x = element_line(color = "grey80", linewidth = .1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = .1),
      panel.grid.minor = element_blank(),
      legend.text = element_text(size = base_size - 1), # Legend text size
      legend.title = element_text(size = base_size), # Legend title size
      legend.key.size = unit(2, "mm")
    ) # Legend key size
}

#' Plot summary of QC filtering
#' @param data MidarExperiment object
#' @param user_defined_keeper Include user-defined feature inclusion list, even they did not pass QC filtering
#' @param base_size font size of plots
#' @return ggplot2 object
#' @export


plot_qc_summary_venn <- function(data, user_defined_keeper, base_size = 12) {
  if (user_defined_keeper) cli::cli_abort("user_defined_keeper = TRUE not yet supported")

  d_qc <- data@metrics_qc |>
    filter(.data$valid_feature, !.data$is_istd) |>
    mutate(feature_class = tidyr::replace_na(.data$feature_class, "Undefined"))


  d_qc_sum_total <- d_qc |>
    ungroup() |>
    summarise(
      has_only_na = sum(!.data$pass_no_na),
      below_lod = sum(.data$pass_no_na & !.data$pass_lod),
      below_sb = sum(.data$pass_lod & .data$pass_no_na & !.data$pass_sb),
      above_cva = sum(.data$pass_lod & .data$pass_no_na & .data$pass_sb & !.data$pass_cva),
      bad_linearity = sum(.data$pass_lod & .data$pass_no_na & .data$pass_sb & .data$pass_cva & !.data$pass_linearity),
      above_dratio = sum(.data$pass_lod & .data$pass_no_na & .data$pass_sb & .data$pass_cva & .data$pass_linearity & !.data$pass_dratio),
      qc_pass = sum(.data$qc_pass)
    ) |>
    tidyr::pivot_longer(names_to = "qc_criteria", values_to = "count_pass", cols = everything()) |>
    ungroup() |>
    mutate(qc_criteria = factor(.data$qc_criteria, c("below_lod", "has_only_na", "below_sb", "above_cva", "above_dratio", "bad_linearity", "qc_pass")))



  # d_qc_sum_total <- d_qc_sum |>
  #   group_by(qc_criteria) |>
  #   summarise(Count = sum(Count,na.rm = TRUE)) |>
  #   mutate(totalCount = sum(Count,na.rm = TRUE),
  #          cumCount = cumsum(Count),
  #          centres = totalCount - (cumCount - Count / 2)) |> ungroup()

  p_bar <- ggplot(d_qc_sum_total, aes(x = reorder(.data$qc_criteria, .data$count_pass), y = .data$count_pass, fill = .data$qc_criteria)) +
    geom_bar(width = 1, stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c(below_lod = "#c7c7c7", has_only_na = "#5e555e", below_sb = "#8f8f8f", above_cva = "#870331", bad_linearity = "#009ec9", above_dratio = "#ffabab", qc_pass = "#02bd62")) +
    # geom_text(aes(label = Count), size=4 ) +
    geom_text(
      aes(
        label = .data$count_pass,
        hjust = ifelse(.data$count_pass < max(.data$count_pass) / 1.5, -2, 2) # <- Here lies the magic
      )
    ) +
    labs(x = "", "Number of analytes")
  # facet_wrap(~Tissue) +
  theme_bw(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      legend.text = element_text(size = base_size - 2), # Legend text size
      legend.title = element_text(size = base_size - 4), # Legend title size
      legend.key.size = unit(3, "mm")
    ) # Legend key size

  # prevent creating log file

  d_qc_venn <- d_qc

  cva_failed <- d_qc_venn$feature_id[!d_qc_venn$pass_cva & d_qc_venn$pass_no_na]
  lod_sb_failed <- d_qc_venn$feature_id[(!d_qc_venn$pass_lod | !d_qc_venn$pass_sb) & d_qc_venn$pass_no_na]
  # sb_failed = d_qc_venn$feature_id[!d_qc_venn$pass_sb  & d_qc_venn$pass_no_na]
  lin_failed <- d_qc_venn$feature_id[!d_qc_venn$pass_linearity & d_qc_venn$pass_no_na]

  lod_sb_label <- "below S/B or LoD" # paste0('LoD < ', parameter_processing$)
  # sb_label <- "below S/B or LOD" #paste0('S/B < ', parameter_processing$)
  cva_label <- "above CVa" # paste0('CV > ', percent(MAX_CV_NORM/100))
  lin_label <- "below R^2" # paste0('RQC r^2 < ', MIN_LINEARITY_RSQUARE, ' OR rel y0 > ', REL_Y_INTERSECT)o

  x2 <- list(lod_sb_failed, cva_failed, lin_failed)
  names(x2) <- c(lod_sb_label, cva_label, lin_label)

  p_venn <- ggvenn::ggvenn(x2, c(lod_sb_label, cva_label, lin_label),
    show_percentage = FALSE,
    fill_color = c("#c7c7c7", "#ff8080", "#009ec9"),
    fill_alpha = 0.5,
    stroke_size = 0.0,
    text_size = 5,
    set_name_size = 6
  ) # + ggplot2::coord_cartesian(clip="off")



  # plt <- arrangeGrob(p_bar, gTree(NULL),gTree(children=p_venn), ncol=3, widths=c(1.6, 0, 1.3), padding = 0)
  plt <- p_bar + p_venn + patchwork::plot_layout(ncol = 2, widths = c(1, 1)) + patchwork::plot_annotation(tag_levels = c("A", "B"))

  return(plt)
}
