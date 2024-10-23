#' Plot summary of feature QC filtering per class
#' @param data MidarExperiment object
#' @param include_qualifier Include qualifier features. Is not used when `filter_data = TRUE` was applied.
#' @param include_istd Include internal standard features. Default is `TRUE`. Is not used when `filter_data = TRUE` was applied.
# #' @param use_batches How batches should be used, either across batches (ignoring batches), plot batch individually, or summarize batches (median)
#' @param user_defined_keeper Include user-defined feature inclusion list, even they did not pass QC filtering
#' @param font_base_size font size of plots
#' @return ggplot2 object
#' @export

# TODO: handling of features with (many) missing values, in SPL, in QC
qc_plot_summary_classes <- function(data, include_qualifier = FALSE, include_istd = FALSE, user_defined_keeper = FALSE, font_base_size = 8) {
  if (user_defined_keeper) cli::cli_abort("user_defined_keeper = TRUE not yet supported")

  #rlang::arg_match(use_batches, c("across", "individual", "summarise"))
  #if(use_batches != "summarise") stop("Currently only `summarise` supported for parameter `batches`")

  if(!data@is_filtered)
    cli_abort(col_red("QC filter has not yet been applied, or data has changed. Please run `qc_apply_feature_filter()` first."))

  d_qc <- data@metrics_qc |>
    filter(.data$valid_feature, .data$in_data) |>
    mutate(feature_class = tidyr::replace_na(.data$feature_class, "Undefined"))

  if(!include_qualifier){
    d_qc <- d_qc |> filter(.data$is_quantifier)
  }

  if(!include_istd){
    d_qc <- d_qc |> filter(!.data$is_istd)
  }


  # TODO: cleanup feature/lipidclasses
  # if(!all(is.na(d_qc$feature_class)) & any(is.na(d_qc$lipid_class))) d_qc$feature_class <- d_qc$lipid_class

  if (all(is.na(d_qc$feature_class))) cli::cli_abort("This plot requires `feature_class` to be defined. Please define classes in the feature metadata or retrieve via corresponding {midar} functions.")

  d_qc$feature_class <- forcats::fct(d_qc$feature_class)

# Count how many features failed qc criteria, excluding features that failed before tested criteria (lower hiarchy)
# TODO: can surely be better implemented
  d_qc_sum <- d_qc |> ungroup() |>
    group_by(.data$feature_class) |>
    summarise(
      has_only_na = sum(.data$na_in_all, na.rm = TRUE),
      exceed_missingness = sum((!replace_na(.data$na_in_all, TRUE) & !replace_na(.data$pass_missingval, TRUE)), na.rm = TRUE),
      below_lod = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE)) & !replace_na(.data$pass_lod, TRUE), na.rm = TRUE),
      below_sb = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE)) & !replace_na(.data$pass_sb, TRUE), na.rm = TRUE),
      above_cva = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & replace_na(.data$pass_sb, TRUE)) & !replace_na(.data$pass_cva, TRUE), na.rm = TRUE),
      bad_linearity = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & replace_na(.data$pass_sb, TRUE) & replace_na(.data$pass_cva, TRUE)) & !replace_na(.data$pass_linearity, TRUE), na.rm = TRUE),
      above_dratio = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & replace_na(.data$pass_sb, TRUE) & replace_na(.data$pass_cva, TRUE) & replace_na(.data$pass_linearity, TRUE)) & !replace_na(.data$pass_dratio, TRUE), na.rm = TRUE),
      qc_pass = sum(.data$qc_pass, na.rm = TRUE)
    ) |>
    tidyr::pivot_longer(-.data$feature_class, names_to = "qc_criteria", values_to = "count_pass") |>
    ungroup() |>
    group_by(.data$feature_class) |>
    mutate(percent_pass = .data$count_pass / sum(.data$count_pass, na.rm = TRUE) * 100)

  qc_colors <- c(qc_pass = "#02bf83", above_dratio = "#b5a2f5", bad_linearity = "#abdeed", above_cva = "#F44336", below_sb = "#d9d5b6", below_lod = "#ada3a3", exceed_missingness = "yellow", has_only_na = "#111111")
  d_qc_sum$qc_criteria <- forcats::fct_relevel(d_qc_sum$qc_criteria, rev(names(qc_colors)))


  # Remove levels/qc criteria for which was not filtered for
  if(all(is.na(d_qc$na_in_all))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "has_only_na")
  if(all(is.na(d_qc$pass_missingval))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "exceed_missingness")
  if(all(is.na(d_qc$pass_lod))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "below_lod")
  if(all(is.na(d_qc$pass_sb))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "below_sb")
  if(all(is.na(d_qc$pass_cva))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "above_cva")
  if(all(is.na(d_qc$pass_linearity))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "bad_linearity")
  if(all(is.na(d_qc$pass_dratio))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "above_dratio")

  d_qc_sum <- d_qc_sum |> drop_na(.data$qc_criteria)
  ggplot(d_qc_sum, aes(x = as.numeric(rev(.data$feature_class)), y = .data$count_pass)) +
    ggplot2::geom_bar(aes(fill = .data$qc_criteria), stat = "identity", na.rm = TRUE) +
    scale_fill_manual(values = qc_colors, na.value = "purple", drop = FALSE) +
    # facet_wrap(~Tissue) +
    # guides(fill = guide_legend(override.aes = list(size = 6))) +
    ggplot2::coord_flip() +
    labs(y = "Number of features", x = "Feature class") +
    #facet_wrap(~batch_id, scales = "free") +
    ggplot2::scale_y_continuous(expand = expansion(0.02, 0.03)) +
    ggplot2::scale_x_continuous(breaks = seq(1, nlevels(d_qc_sum$feature_class), by = 1),
                                labels = rev(levels(d_qc_sum$feature_class)),
                                expand = expansion(0.02, 0.02),
                                sec.axis = ggplot2::sec_axis(~ ., name = "Features passed QC",
                                                  breaks = seq(1, nlevels(d_qc_sum$feature_class), by = 1),
                                                  labels = rev(d_qc_sum |> filter(.data$qc_criteria == "qc_pass") |> mutate(txt = (paste0(.data$count_pass, " (", round(.data$percent_pass,0), "%)"))) |> pull(.data$txt)))) +
    theme_bw(base_size = font_base_size) +
    theme(
      legend.position = c(0.8, 0.8),
      panel.grid.major.x = element_line(color = "grey80", linewidth = .1),
      panel.grid.major.y = element_line(color = "grey90", linewidth = .1),
      panel.grid.minor = element_blank(),
      axis.title.y.right = element_text(angle = 90, vjust = 0.5, hjust = 0.5),
      legend.text = element_text(size = font_base_size - 1), # Legend text size
      legend.title = element_text(size = font_base_size), # Legend title size
      legend.key.size = unit(2, "mm")
    ) # Legend key size
}

#' Plot summary of feature QC filtering
#' @param data MidarExperiment object
# #' @param use_batches How batches should be used, either across batches (ignoring batches), plot batch individually, or summarize batches (median)
#' @param include_qualifier Include qualifier features. Is not used when `filter_data = TRUE` was applied.
#' @param include_istd Include internal standard features. Default is `TRUE`. Is not used when `filter_data = TRUE` was applied.
#' @param with_venn_diag Include Venn diagram of features failing S/B, CVa, and linearity
#' @param user_defined_keeper Include user-defined feature inclusion list, even they did not pass QC filtering
#' @param font_base_size font size of plots
#' @return ggplot2 object
#' @export



qc_plot_summary <- function(data, include_qualifier = FALSE, include_istd = FALSE, with_venn_diag = TRUE, user_defined_keeper = FALSE, font_base_size = 8) {

  if (user_defined_keeper) cli::cli_abort("user_defined_keeper = TRUE not yet supported")
  #rlang::arg_match(use_batches, c("across", "individual", "summarise"))
  #if(use_batches != "summarise") stop("Currently only `summarise` supported for parameter `batches`")

  if(!data@is_filtered)
    cli_abort(col_red("QC filter has not yet been applied, or data has changed. Please run `qc_apply_feature_filter()` first."))

  d_qc <- data@metrics_qc |>
    filter(.data$valid_feature, .data$in_data) |>
    mutate(feature_class = tidyr::replace_na(.data$feature_class, "Undefined"))


  if(!include_qualifier){
    d_qc <- d_qc |> filter(.data$is_quantifier)
  }

  if(!include_istd){
    d_qc <- d_qc |> filter(!.data$is_istd)
  }

  d_qc_sum <- d_qc |>
    ungroup() |>
    summarise(
      has_only_na = sum(.data$na_in_all, na.rm = TRUE),
      exceed_missingness = sum((!replace_na(.data$na_in_all, TRUE) & !replace_na(.data$pass_missingval, TRUE)), na.rm = TRUE),
      below_lod = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE)) & !replace_na(.data$pass_lod, TRUE), na.rm = TRUE),
      below_sb = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE)) & !replace_na(.data$pass_sb, TRUE), na.rm = TRUE),
      above_cva = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & replace_na(.data$pass_sb, TRUE)) & !replace_na(.data$pass_cva, TRUE), na.rm = TRUE),
      bad_linearity = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & replace_na(.data$pass_sb, TRUE) & replace_na(.data$pass_cva, TRUE)) & !replace_na(.data$pass_linearity, TRUE), na.rm = TRUE),
      above_dratio = sum((!replace_na(.data$na_in_all, TRUE) & replace_na(.data$pass_missingval, TRUE) & replace_na(.data$pass_lod, TRUE) & replace_na(.data$pass_sb, TRUE) & replace_na(.data$pass_cva, TRUE) & replace_na(.data$pass_linearity, TRUE)) & !replace_na(.data$pass_dratio, TRUE), na.rm = TRUE),
      qc_pass = sum(.data$qc_pass, na.rm = TRUE)
    ) |>
    tidyr::pivot_longer(names_to = "qc_criteria", values_to = "count_pass", cols = everything()) |>
    ungroup() |>
    mutate(percent_pass = .data$count_pass / sum(.data$count_pass, na.rm = TRUE) * 100) |>
    ungroup() |>
    mutate(qc_criteria = factor(.data$qc_criteria, c("exceed_missingness", "below_lod", "has_only_na", "below_sb", "above_cva", "above_dratio", "bad_linearity", "qc_pass")))

  qc_colors <- c(qc_pass = "#02bf83", above_dratio = "#b5a2f5", bad_linearity = "#abdeed", above_cva = "#F44336", below_sb = "#d9d5b6", below_lod = "#ada3a3", exceed_missingness = "yellow", has_only_na = "#111111")
  d_qc_sum$qc_criteria <- forcats::fct_relevel(d_qc_sum$qc_criteria, rev(names(qc_colors)))


  # Remove levels/qc criteria for which was not filtered for
  if(all(is.na(d_qc$na_in_all)) | d_qc_sum$count_pass[d_qc_sum$qc_criteria == "has_only_na"] == 0) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "has_only_na")
  if(all(is.na(d_qc$pass_missingval))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "exceed_missingness")
  if(all(is.na(d_qc$pass_lod))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "below_lod")
  if(all(is.na(d_qc$pass_sb))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "below_sb")
  if(all(is.na(d_qc$pass_cva))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "above_cva")
  if(all(is.na(d_qc$pass_linearity))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "bad_linearity")
  if(all(is.na(d_qc$pass_dratio))) d_qc_sum$qc_criteria <- forcats::fct_recode(d_qc_sum$qc_criteria, NULL = "above_dratio")

  p_bar <- ggplot(d_qc_sum |> drop_na(.data$qc_criteria), aes(x = .data$qc_criteria, y = .data$count_pass, fill = .data$qc_criteria)) +
    geom_bar(width = 1, stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = qc_colors) +
    ggplot2::scale_y_continuous(expand = expansion(mult = c(0.02, 0.1))) +
    ggplot2::scale_x_discrete(expand = expansion(0.12, 0.12)) +
    # geom_text(aes(label = Count), size=4 ) +
    geom_text(
      aes(
        label = .data$count_pass,
        hjust = ifelse(.data$count_pass < max(.data$count_pass) / 1.5, -2, 2), # <- Here lies the magic
      ),
      size = 2
    ) +
    labs(x = "", y= "Number of analytes") +
  # facet_wrap(~Tissue) +
  theme_bw(base_size = font_base_size) +
    theme(
      legend.position = "none",
      #axis.text.x = element_blank(),
      panel.grid.major.y = element_blank(), #element_line(color = "grey80", linewidth = .1),
      panel.grid.major.x = element_line(color = "grey80", linewidth = .2, linetype = "dotted"),
      panel.grid.minor = element_blank()
    ) # Legend key size

  # prevent creating log file
  if(with_venn_diag) {
    d_qc_venn <- d_qc

    sb_failed <- d_qc_venn$feature_id[!replace_na(d_qc_venn$na_in_all, TRUE) & replace_na(d_qc_venn$pass_lod, TRUE) & !replace_na(d_qc_venn$pass_sb, TRUE)]
    cva_failed <- d_qc_venn$feature_id[!replace_na(d_qc_venn$na_in_all, TRUE) & replace_na(d_qc_venn$pass_lod, TRUE) & !replace_na(d_qc_venn$pass_cva, TRUE)]
    lin_failed <- d_qc_venn$feature_id[!replace_na(d_qc_venn$na_in_all, TRUE) & replace_na(d_qc_venn$pass_lod, TRUE) & !replace_na(d_qc_venn$pass_linearity, TRUE)]

    sb_label <- "below S/B"
    cva_label <- "above CV(A)" # paste0('CV > ', percent(MAX_CV_NORM/100))
    lin_label <- "bad linearity" # paste0('RQC r^2 < ', MIN_LINEARITY_RSQUARE, ' OR rel y0 > ', REL_Y_INTERSECT)o

    x2 <- list(sb_failed, cva_failed, lin_failed)
    names(x2) <- c(sb_label, cva_label, lin_label)

    p_venn <- ggvenn::ggvenn(x2, c(sb_label, cva_label, lin_label),
      show_percentage = FALSE,
      fill_color = c("#d9d5b6", "#F44336", "#abdeed"),
      fill_alpha = 0.5,
      stroke_size = 0.0,
      text_size = 2,
      set_name_size = 2
    ) # + ggplot2::coord_cartesian(clip="off")


    plt <- patchwork::wrap_plots(p_bar, p_venn) +
      patchwork::plot_layout(ncol = 2, widths = c(1.3, 1)) +
      patchwork::plot_annotation(tag_levels = c("A", "B")) &
      theme(plot.tag = element_text(face = 'bold', size = font_base_size, color = 'black', hjust = 0, vjust = 0, margin = margin(10, 20, 10, 10)))
  } else {
    plt <- p_bar
  }

  return(plt)
}
