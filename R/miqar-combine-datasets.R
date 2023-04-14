
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


combine_experiments <- function(..., ordered_by_runsequence){

  exp_list <- list(...)

  mexp <- MidarExperiment()


  mexp@dataset <- purrr::map_dfr(.x = exp_list,  .f = \(x) x@dataset)  |> dplyr::distinct()
  mexp@annot_analyses <- purrr::map_dfr(.x = exp_list,  .f = \(x) x@annot_analyses) |> mutate(RUN_ID_ANNOT = dplyr::row_number())
  mexp@annot_istd <- purrr::map_dfr(.x = exp_list,  .f = \(x) x@annot_istd)  |> dplyr::distinct()
  mexp@annot_features <- purrr::map_dfr(.x = exp_list,  .f = \(x) x@annot_features) |> dplyr::distinct()

  mexp@dataset <- mexp@dataset %>%
    dplyr::rename(BATCH_RUN_ID = .data$RUN_ID) %>%
    dplyr::group_by(.data$FEATURE_NAME) %>%
    dplyr::mutate(RUN_ID = dplyr::row_number(), .before = .data$BATCH_RUN_ID) %>%
    dplyr::ungroup()

  mexp@annot_batch_info <- mexp@annot_analyses %>%
    dplyr::group_by(.data$BATCH_ID) %>%
    dplyr::mutate(BATCH_NO = dplyr::cur_group_id()) %>%
    dplyr::summarise(
      BATCH_ID = .data$BATCH_ID[1],
      BATCH_NO = .data$BATCH_NO[1],
      id_batch_start = dplyr::first(.data$RUN_ID_ANNOT),
      id_batch_end = dplyr::last(.data$RUN_ID_ANNOT)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$id_batch_start) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_batch_info_template)

  mexp

}









