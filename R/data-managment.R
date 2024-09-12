
check_rawdata_present <- function(data){nrow(data@dataset_orig) > 0}
check_dataset_present <- function(data){nrow(data@dataset) > 0}

#' @title Set the analysis order
#' @description
#' Sets the analysis order (sequence), based on either (i) analysis timestamp, if available, (ii) the order in which analysis appeared in the imported raw data file, or (iii) the order in which analyses were defined in the Analysis metadata
#' @param data MidarExperiment object
#' @param analysis_sequence Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
#' @return MidarExperiment object
#' @examples
#' file_path <- system.file("extdata", "sPerfect_MRMkit.tsv", package = "midar")
#' mexp <- MidarExperiment()
#' mexp <- rawdata_import_mrmkit(data = mexp, path = file_path, use_metadata = TRUE)
#' mexp <- set_analysis_order(data, "timestamp")

#' @export

#c("timestamp", "resultfile", "metadata")
set_analysis_order <- function(data, analysis_sequence = "default"){

  analysis_sequence <- rlang::arg_match(arg = analysis_sequence, c("timestamp", "resultfile", "metadata", "default"), multiple = FALSE)

  d_temp <- data@dataset_orig |>
    select("analysis_id", any_of("acquisition_time_stamp")) |>
    distinct()

  if (analysis_sequence == "default"){
    if ("acquisition_time_stamp" %in% names(d_temp))
        analysis_sequence <- "timestamp"
    else {
      analysis_sequence <- "resultfile"
      cli::cli_alert_info(cli::col_grey(glue::glue("Analysis order was based on sequence of analysis results, as no timestamps were found. Use `set_analysis_order` to define alternative analysis orders.")))
    }
  }

  if (analysis_sequence %in% c("timestamp")) {
    if("acquisition_time_stamp" %in% names(d_temp)){
      d_temp <- d_temp |>
        arrange(.data$acquisition_time_stamp) |>
        mutate(run_id = row_number()) |>
        select(-"acquisition_time_stamp")


      data@annot_analyses <- data@annot_analyses |>
        left_join(d_temp,  by = "analysis_id")
    } else {
      cli::cli_abort(call. = FALSE, "No acquisition timestamp field present in analysis results, please set parameter `analysis_sequence` to `resultfile` or `metadata`.")
    }
  } else if (analysis_sequence == "resultfile") {
    d_temp <- d_temp |> mutate(run_id = row_number())
    data@annot_analyses <- data@annot_analyses |>
      left_join(d_temp,  by = "analysis_id")
  } else if (analysis_sequence == "metadata") {
      data@annot_analyses <- data@annot_analyses |>  mutate(run_id = row_number())
  }
  data
}


# # ##### TODO TODO =====================
# d_dataset <- data@dataset_orig %>%
#   dplyr::inner_join(
#     data@annot_analyses %>%
#       dplyr::select("run_id", "analysis_id", "qc_type", "specimen", "sample_id", "replicate_no", "valid_analysis", "batch_id"),
#     by = c("analysis_id")) %>%
#   dplyr::inner_join(
#     metadata$annot_features %>%
#       filter(.data$valid_feature) |>
#       dplyr::select(dplyr::any_of(c("feature_id", "feature_id", "feature_class", "norm_istd_feature_id", "quant_istd_feature_id", "is_istd", "feature_id", "is_quantifier", "valid_feature", "table_templates", "interference_feature_id", "interference_proportion"))),
#     by = c("feature_id"), keep = FALSE
#   )



# Link DATA with METADATA and create DATASET table. =================
## - Only valid analyses and features will be added
## - Only key information will be added
link_data_metadata <- function(data, minimal_info = TRUE){
  data@dataset <- data@dataset_orig |>
    select(
      "analysis_id",
      "acquisition_time_stamp",
      "feature_id",
      starts_with("method_"),
      starts_with("feature_")
    ) |>
    inner_join(data@annot_analyses, by = "analysis_id") |>
    inner_join(data@annot_features, by = "feature_id") |>
    filter(valid_analysis, valid_feature) |>
    select(
      "run_id",
      "analysis_id",
      "acquisition_time_stamp",
      "qc_type",
      "batch_id",
      "feature_id",
      "feature_class",
      "is_istd",
      starts_with("method_"),
      starts_with("feature_")
    )
  # NOTE: To adjust
  #data@dataset <-data@dataset |> select(all_of("feature_intensity" ))

  if(minimal_info)
    data@dataset <- data@dataset |> select(-starts_with("method_"), -starts_with("feature_int"))



  data@dataset <- data@dataset |>
    dplyr::bind_rows(pkg.env$table_templates$dataset_template)

  # NOTE: To adjust
    # mutate(
    #   corrected_interference = FALSE,
    #   outlier_technical = FALSE
    # )

  # Arrange run_id and then by feature_id, as they appear in the metadata
  # TODO:  check if as intended
  data@dataset <- data@dataset |>
    dplyr::arrange(match(.data$feature_id, data@annot_features$feature_id)) |>
    arrange(.data$run_id)

  check_integrity(data, excl_unannotated_analyses = excl_unannotated_analyses)

  data
}





