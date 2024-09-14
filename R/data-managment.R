
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
      "is_quantifier",
      starts_with("method_"),
      starts_with("feature_")
    )



  if(minimal_info)
    data@dataset <- data@dataset |> select(-starts_with("method_"), -starts_with("feature_int_"))

  data@dataset <- data@dataset |>
    dplyr::bind_rows(pkg.env$table_templates$dataset_template)

  data@dataset <- data@dataset |> mutate(feature_intensity = !!(sym(data@feature_intensity_var)))


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

  # TODOTODO: DECIDE WHEN WHERE TO RUN THIS
  # check_integrity(data, excl_unannotated_analyses = excl_unannotated_analyses)

  data
}

#' @title Set default variable to be used as feature raw signal value
#' @description
#' Sets the raw signal variable used for calculations starting from raw signal
#' values (i.e., normalization)
#' @param data MidarExperiment object
#' @param analysis_sequence Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
#' @param analysis_sequence Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
#' @return MidarExperiment object
#' @export

set_intensity_var <- function(data, variable_name, auto_select = FALSE, ...){
  if (auto_select) {
    var_list <- unlist(rlang::list2(...), use.names = FALSE)
    idx = match(var_list,names(data@dataset_orig))[1]
    if (!is.na(idx)) {
      data@feature_intensity_var = var_list[1]
      cli_alert_info(text = cli::col_grey("{.var {var_list[1]}} selected as default raw signal intensity. Use {.fn set_intensity_var} to modify."))
      variable_name <- var_list[1]
      } else {
      cli_alert_warning(text = cli::col_yellow("No typical raw signal intensity variable found in the data. Use {.fn set_intensity_var to set it.}}"))
      return(data)
      }
  } else {
    #TODO: Double check behavior if there a feature_intensity in the raw data file
    if (! variable_name %in% names(data@dataset_orig))
      cli_abort(c("x" = "{.var variable_name} is not present in the raw data."))

    if (! variable_name %in% c("feature_area", "feature_height", "feature_response", "feature_intensity"))
      cli_alert_warning(text = "{.var {variable_name}} is not a typically used raw signal name (i.e., area, response, intensity, height).")
  }
  data@feature_intensity_var <- variable_name

  if (check_dataset_present(data)) {
    calc_cols <- c("featue_norm_intensity", "feature_conc", "feature_amount", "feature_raw_conc")
    if (any(calc_cols %in% names(data@dataset))){
      data@dataset <- data@dataset |> select(-any_of(calc_cols))
      cli_alert_warning(text = "New feature variable set as default raw signal, please re-process data")
    }
    data <- link_data_metadata(data)
    }
  data
}


