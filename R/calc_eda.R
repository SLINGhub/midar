#' Get FC and p values from a t test contrast
#'
#' @param data data.frame or tibble
#' @param feature feature
#' @param grouping grouping
#' @param group_case case group
#' @param group_ref reference group
#' @param paired paired data
#' @param min_groupsize min group size for t test
#' @param log_transform log transform before t test
#'
#' @return tibble with p values, fdr and log2FC
#' @noRd
get_summary_stats <- function(data, feature, grouping, group_case, group_ref, paired, min_groupsize = 3, log_transform = FALSE) {
  if (min_groupsize <= 2) cli::cli_abort("group size must be >= 2 for the statistical tests")

  feature <- rlang::ensym(feature)
  grouping <- rlang::ensym(grouping)
  f <- as.formula(paste0("feature_conc ~ ", grouping))

  d_filt <- data |>
    filter(!!grouping == group_case | !!grouping == group_ref) |>
    mutate(!!grouping := forcats::fct_relevel(!!grouping, c(group_case, group_ref)))

  d_incompletegroups <- d_filt |>
    group_by(!!feature) |>
    summarise(
      grp_size_fail =
        sum(!is.na(.data$feature_conc[!!grouping == group_case])) < min_groupsize |
        sum(!is.na(.data$feature_conc[!!grouping == group_ref])) < min_groupsize
    ) |>
    filter(.data$grp_size_fail)

  incomplete <- paste(d_incompletegroups |> pull(!!feature), collapse = ", ")

  if (nrow(d_incompletegroups) > 0) {
    warning(glue::glue("Following features were ignored as they have group size(s) of less than {min_groupsize}:\n {incomplete}"))
  }

  d_stats <- d_filt |>
    rowwise() |>
    mutate(feature_conc = if_else(log_transform, log2(.data$feature_conc), .data$feature_conc)) |>
    ungroup() |>
    anti_join(d_incompletegroups, by = join_by(!!feature)) |>
    group_by(!!feature) |>
    # drop_na(feature_conc) |>
    nest() |>
    mutate(res = purrr::map(data, \(x) broom::tidy(
      t.test(f, paired = paired, var.equal = FALSE, alternative = "two.sided", data = x)
    ))) |>
    unnest(.data$res) |>
    ungroup() |>
    mutate(FDR = p.adjust(.data$p.value, method = "BH")) |>
    select(!!feature, p_value = .data$p.value, .data$FDR)

  if (!paired) {
    d_fc <- d_filt |>
      group_by(!!feature) |>
      summarise(log2FC = log2(mean(.data$feature_conc[!!grouping == group_case]) / mean(.data$feature_conc[!!grouping == group_ref])))
  } else {
    d_fc <- d_filt |>
      group_by(!!feature) |>
      summarise(log2FC = mean(log2(.data$feature_conc[!!grouping == group_case] / .data$feature_conc[!!grouping == group_ref])))
  }

  d_stats |> inner_join(d_fc, by = "feature_id")
}
