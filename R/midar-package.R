#' @keywords internal
"_PACKAGE"

## usethis namespace: start
# #' @importFrom future plan multisession sequential
#' @import rlang
#' @importFrom assertr chain_start chain_end verify assert has_all_names in_set is_uniq not_na
#' @importFrom cli cli_abort cli_alert_success cli_alert_warning cli_alert_info cli_alert_danger cli_alert col_green col_red col_yellow
#' @importFrom dplyr n select filter group_by desc ungroup mutate summarise slice if_else rowwise arrange left_join right_join inner_join full_join anti_join semi_join join_by row_number cur_group_id case_when distinct rename relocate pull across all_of any_of if_any ends_with if_all case_match case_when bind_rows group_split pick
#' @importFrom fs is_dir path_tidy file_exists dir_ls
#' @importFrom ggplot2 ggplot aes Stat unit geom_point geom_line geom_abline margin coord_flip geom_text geom_bar geom_segment scale_y_discrete scale_x_discrete scale_y_log10 xlab ylab vars theme facet_wrap position_jitter scale_color_manual labs scale_fill_manual theme_light geom_vline geom_rect geom_vline aes_string scale_fill_manual scale_shape_manual expand_limits geom_smooth geom_hline scale_y_continuous element_blank element_text theme_bw stat_ellipse element_rect element_line expansion scale_x_continuous geom_boxplot ggtitle position_dodge
#' @importFrom glue glue
#' @importFrom grDevices pdf dev.off dev.flush
#' @importFrom methods is validObject
#' @importFrom pillar ctl_new_pillar new_pillar tbl_format_footer tbl_format_header tbl_sum
#' @importFrom purrr map_dfr map
#' @importFrom rlang .data ensym sym := arg_match is_na are_na caller_env
#' @importFrom stats as.formula lm median na.exclude sd quantile na.omit setNames p.adjust t.test IQR dnorm mad reorder prcomp model.matrix predict
#' @importFrom stringr str_remove str_replace str_replace_all str_detect str_trim str_extract str_squish fixed
#' @importFrom tibble column_to_rownames as_tibble tibble
#' @importFrom tidyr unite drop_na pivot_wider nest unnest replace_na
#' @importFrom tidyselect vars_select_helpers everything starts_with
#' @importFrom utils tail txtProgressBar setTxtProgressBar head flush.console
## usethis namespace: end
NULL


.datatable.aware <- TRUE
