#' Sum up feature intensities per analyte
#'
#' @description
#' Retrieves analysis IDs of data outliers based  the principal components PCA
#' with SD or MAD fences
#'
#' @param data MidarExperiment object
#' @param qualifier_action Character. How to handle qualifier features. To sum them up separately select "separate", 
#' @param feature_classes Select feature classes to sum up. Default (NA) is to sum up all feature classes according to their analyte_id.
#' to include them in the sum if quantifier select "include", to not sum them up select "exclude".
#' 
#' @return MidarExperiment object
#' @export
#' 

data_sum_features <- function(
  data,
  qualifier_action = "separate",
  feature_classes = NA
) {
  qualifier_action <- match.arg(qualifier_action, c("separate", "include", "exclude"))
  ds <- data@dataset

  # Numeric feature columns to sum
  feature_cols <- setdiff(
    names(ds)[sapply(ds, is.numeric) & grepl("^feature_", names(ds))],
    c("feature_id", "feature_class")
  )

  # Metadata columns to keep
  metadata_cols <- c("feature_class", "feature_label")

  # Other columns (non-numeric, not metadata)
  other_cols <- setdiff(names(ds), c(feature_cols, "feature_id"))

  aggregate_quant <- function(df) {
    # Sum numeric feature columns per group
    agg_numeric <- df |>
      dplyr::group_by(.data$analysis_id, .data$analyte_id) |>
      dplyr::summarise(
        dplyr::across(any_of(c("feature_intensity", "feature_height", "feature_fwhm")), ~ sum(.x, na.rm = TRUE)),
        dplyr::across(any_of(c("feature_rt")), ~ mean(.x, na.rm = TRUE)),
        .groups = "drop"
      )

    # Get metadata columns from first row per group
    agg_meta <- df |>
      dplyr::select("analysis_id", "analyte_id", all_of(metadata_cols), all_of(setdiff(other_cols, metadata_cols))) |>
      dplyr::distinct(.data$analysis_id, .data$analyte_id, .keep_all = TRUE)

    # Merge numeric sums and metadata
    agg <- dplyr::left_join(agg_numeric, agg_meta, by = c("analysis_id", "analyte_id"))

    # Replace feature_id with analysis_id in aggregated rows
    agg$feature_id <- agg$analyte_id
    agg
  }

  if (qualifier_action == "include") {
    ds_res <- aggregate_quant(ds) |>
      dplyr::mutate(is_quantifier = TRUE)

  } else if (qualifier_action == "exclude") {
    ds_res <- aggregate_quant(dplyr::filter(ds, .data$is_quantifier)) |>
      dplyr::mutate(is_quantifier = TRUE)

  } else { # separate
    quant <- aggregate_quant(dplyr::filter(ds, .data$is_quantifier)) |>
      dplyr::mutate(is_quantifier = TRUE)
    qual <- dplyr::filter(ds, !.data$is_quantifier)
    ds_res <- dplyr::bind_rows(quant, qual)
  }

  ds_res <- ds_res |> 
  select(any_of(names(data@dataset)))
  
  data@dataset <- ds_res

  data@annot_features <- data@annot_features |> 
    mutate(feature_id = .data$analyte_id) |>
    distinct()

  data
}


# data_sum_features <- function(
#   data,
#   feature_classes = NA,
#   qualifier_action = c("separate", "include", "exclude"),
#   include_feature_filter = NA,
#   exclude_feature_filter = NA
# ) {
#   qualifier_action <- match.arg(qualifier_action)
#   ds <- data@dataset

#   if (qualifier_action == "include") {
#     # combine quant + qual
#     ds_res <- ds |>
#       dplyr::group_by(analysis_id, analyte_id) |>
#       dplyr::summarise(
#         feature_intensity = sum(.data$feature_intensity, na.rm = TRUE),
#         .groups = "drop"
#       ) |>
#       dplyr::mutate(is_quantifier  = TRUE)

#   } else if (qualifier_action == "exclude") {
#     # only quantifiers
#     ds_res <- ds |>
#       dplyr::filter(is_quantifier) |>
#       dplyr::group_by(analysis_id, analyte_id) |>
#       dplyr::summarise(
#         feature_intensity = sum(.data$feature_intensity, na.rm = TRUE),
#         .groups = "drop"
#       ) |>
#       dplyr::mutate(is_quantifier = TRUE)

#   } else { # "separate"
#     quant <- ds |>
#       dplyr::filter(is_quantifier) |>
#       dplyr::group_by(analysis_id, analyte_id) |>
#       dplyr::summarise(
#         dplyr::across(any_of(c("feature_intensity", "feature_height", "feature_fwhm")), ~ sum(.x, na.rm = TRUE)),
#         dplyr::across(any_of(c("feature_class", "qc_type")), ~ first(.x)),
#         .groups = "drop"
#       ) |>
      
#       dplyr::mutate(is_quantifier = TRUE)

#     qual <- ds |> dplyr::filter(!is_quantifier)

#     common_cols <- intersect(names(quant), names(qual))
#     ds_res <- dplyr::bind_rows(
#       dplyr::select(quant, dplyr::all_of(common_cols)),
#       dplyr::select(qual, dplyr::all_of(common_cols))
#     )
#   }

#   #data@dataset <- ds_res
#   ds_res
# }
