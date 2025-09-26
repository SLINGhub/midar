#' Sum up feature intensities per analyte
#'
#' @description
#' This function sums up feature intensities per analyte_id.
#'
#' THis is is useful when you have multiple features (e.g. adducts, isotopes, in-source fragments) or isomers that you
#' want to combine into a single analyte intensity value, such as LPC sn1 and sn2 species.
#'
#' NOTE: This is still an experimental function! It will overwrite the feature_id in the dataset and analysis metadata
#' of featured that share same analyte_id. Currently the original feature_id is not backed up anywhere. Use with caution
#' and check results carefully!
#'
#'
#' @param data MidarExperiment object
#' @param qualifier_action Character. How to handle qualifier features. To sum them up separately select "separate",
#' @param feature_classes Select feature classes to sum up. Default (NA) is to sum up all feature classes according to their analyte_id.
#' to include them in the sum if quantifier select "include", to not sum them up select "exclude".
#'
#' @return MidarExperiment object
#' @export
#'
# The main function structure remains the same
data_sum_features <- function(
  data,
  qualifier_action = "include",
  feature_classes = NA
) {
  qualifier_action <- match.arg(
    qualifier_action,
    c("separate", "include", "exclude")
  )
  ds <- data@dataset

  ds_na_analytes <- ds |> dplyr::filter(is.na(.data$analyte_id))
  ds_to_process <- ds |> dplyr::filter(!is.na(.data$analyte_id))

  aggregate_quant <- function(df) {
    if (nrow(df) == 0) {
      return(df)
    }

    feature_cols <- setdiff(
      names(df)[sapply(df, is.numeric) & grepl("^feature_", names(df))],
      c("feature_id", "feature_class")
    )
    metadata_cols <- c("feature_class", "feature_label")
    other_cols <- setdiff(names(df), c(feature_cols, "feature_id"))

    agg_numeric <- df |>
      dplyr::group_by(.data$analysis_id, .data$analyte_id) |>
      dplyr::summarise(
        dplyr::across(
          any_of(c("feature_intensity", "feature_height", "feature_fwhm")),
          ~ sum(.x, na.rm = TRUE)
        ),
        dplyr::across(any_of(c("feature_rt")), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )

    agg_meta <- df |>
      dplyr::select(
        "analysis_id",
        "analyte_id",
        all_of(metadata_cols),
        all_of(setdiff(other_cols, metadata_cols))
      ) |>
      dplyr::distinct(.data$analysis_id, .data$analyte_id, .keep_all = TRUE)

    agg <- dplyr::left_join(
      agg_numeric,
      agg_meta,
      by = c("analysis_id", "analyte_id")
    )
    agg$feature_id <- agg$analyte_id
    agg
  }

  # --- MODIFICATION: Add a 'suffix' argument to the helper function ---
  process_and_aggregate_subset <- function(df, suffix = "") {
    if (nrow(df) == 0) {
      return(df)
    }

    df_tagged <- df |>
      dplyr::group_by(.data$analysis_id, .data$analyte_id) |>
      dplyr::mutate(.group_size = n()) |>
      dplyr::ungroup()

    rows_to_aggregate <- df_tagged |> dplyr::filter(.data$.group_size > 1)
    rows_to_keep <- df_tagged |>
      dplyr::filter(.data$.group_size <= 1) |>
      dplyr::select(-".group_size")

    aggregated_data <- aggregate_quant(rows_to_aggregate)

    # --- MODIFICATION: Apply the suffix if provided ---
    # This targets ONLY the newly aggregated rows.
    if (nchar(suffix) > 0 && nrow(aggregated_data) > 0) {
      aggregated_data <- aggregated_data |>
        dplyr::mutate(feature_id = paste0(.data$feature_id, suffix))
    }

    dplyr::bind_rows(aggregated_data, rows_to_keep)
  }

  if (qualifier_action == "include") {
    ds_res <- process_and_aggregate_subset(ds_to_process) |>
      dplyr::mutate(is_quantifier = TRUE)
  } else if (qualifier_action == "exclude") {
    quant_rows <- dplyr::filter(ds_to_process, .data$is_quantifier)
    ds_res <- process_and_aggregate_subset(quant_rows) |>
      dplyr::mutate(is_quantifier = TRUE)
  } else {
    # separate
    quant_rows <- dplyr::filter(ds_to_process, .data$is_quantifier)
    qual_rows <- dplyr::filter(ds_to_process, !.data$is_quantifier)

    processed_quants <- process_and_aggregate_subset(quant_rows)

    # --- MODIFICATION: Pass the suffix when processing qualifiers ---
    processed_quals <- process_and_aggregate_subset(qual_rows, suffix = "_qual")

    ds_res <- dplyr::bind_rows(processed_quants, processed_quals)
  }

  ds_res <- dplyr::bind_rows(ds_res, ds_na_analytes)

  ds_res <- ds_res |>
    dplyr::select(any_of(names(data@dataset)))

  data@dataset <- ds_res

  annot <- data@annot_features |>
    mutate(
      is_duplicate = duplicated(.data$analyte_id) |
        duplicated(.data$analyte_id, fromLast = TRUE)
    ) |>
    mutate(
      feature_id = if_else(
        !is.na(.data$analyte_id) & .data$analyte_id != "" & .data$is_duplicate,
        .data$analyte_id,
        .data$feature_id
      )
    ) |>
    distinct()

  data@annot_features <- annot
  data
}
