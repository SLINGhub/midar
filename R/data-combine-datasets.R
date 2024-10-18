#' Combines a list of MidarExperiments into one
#'
#' @param ... MidarExperiment objects
#' @param ordered_by_runsequence Boolean if list of provided MidarExperiment objects is in the run order
#' @export
#'
#' @importFrom glue glue
#' @importFrom openxlsx write.xlsx
#' @importFrom lubridate now
#' @importFrom tibble tribble
#' @importFrom utils packageVersion
#'
#'


combine_experiments <- function(..., ordered_by_runsequence) {
  exp_list <- list(...)
  if (is.null(attr(exp_list, which = "class")[[1]])) exp_list <- exp_list[[1]]

  mexp <- MidarExperiment()
  mexp@dataset_orig <- purrr::map_dfr(.x = exp_list, .f = \(x) x@dataset_orig) |> dplyr::distinct()
  mexp@dataset <- purrr::map_dfr(.x = exp_list, .f = \(x) x@dataset) |> dplyr::distinct()
  mexp@annot_analyses <- purrr::map_dfr(.x = exp_list, .f = \(x) x@annot_analyses) |>
    dplyr::distinct() |>
    mutate(analysis_seq_num = dplyr::row_number())
  mexp@annot_istd <- purrr::map_dfr(.x = exp_list, .f = \(x) x@annot_istd) |> dplyr::distinct()
  mexp@annot_features <- purrr::map_dfr(.x = exp_list, .f = \(x) x@annot_features) |> dplyr::distinct()
  # ToDo: Combine batch and curve id to give unique curve id over the combined experiment
  mexp@annot_responsecurves <- purrr::map_dfr(.x = exp_list, .f = \(x) x@annot_responsecurves) |> dplyr::distinct()

  mexp@dataset <- mexp@dataset |>
    dplyr::rename(batch_analysis_seq_num = .data$analysis_seq_num) |>
    dplyr::group_by(.data$feature_id) |>
    dplyr::mutate(analysis_seq_num = dplyr::row_number(), .before = .data$batch_analysis_seq_num) |>
    dplyr::ungroup()

  mexp@annot_batches <- mexp@annot_analyses |>
    dplyr::group_by(.data$batch_id) |>
    dplyr::summarise(
      batch_id = .data$batch_id[1],
      batch_no = .data$batch_no[1],
      id_batch_start = dplyr::first(.data$analysis_seq_num),
      id_batch_end = dplyr::last(.data$analysis_seq_num)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$id_batch_start) |>
    dplyr::mutate(batch_no = dplyr::row_number()) |>
    dplyr::bind_rows(pkg.env$table_templates$annot_batch_info_template)
  mexp
}
