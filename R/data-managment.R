

# TODO: export functions
check_data_present <- function(data){nrow(data@dataset_orig) > 0}
check_dataset_present <- function(data){nrow(data@dataset) > 0}

get_analysis_count <- function(data, qc_types = NULL) {
  if (is.null(qc_types))
    return(data@dataset|> select("analysis_id") |> distinct() |> nrow())
  else
    return(data@dataset|> filter(.data$qc_type %in% qc_types) |>
             select("analysis_id") |> distinct() |> nrow())
}

get_feature_count <- function(data, isistd = NULL, isquantifier = NULL) {
  d <- data@dataset
  if (!is.null(isistd)) d <- d |> filter(.data$is_istd == isistd)
  if (!is.null(isquantifier)) d <- d |> filter(.data$is_quantifier == isquantifier)

  d |> select("feature_id") |> distinct() |> nrow()
}


get_analyis_start <- function(data){
  if (check_data_present(data))
    return(min(data@dataset$acquisition_time_stamp))
  else
    return(NA)
}

get_analyis_end <- function(data){
  if (check_data_present(data))
    return(max(data@dataset$acquisition_time_stamp))
  else
    return(NA)
}

get_run_time <- function(data){
  if (check_data_present(data))
    median(diff(unique(data@dataset$acquisition_time_stamp), units = "secs")) |>  lubridate::seconds_to_period()
  else
    return(NA)
}

get_analysis_duration <- function(data){
  if (check_data_present(data))
    difftime(max(unique(data@dataset$acquisition_time_stamp)), min(unique(data@dataset$acquisition_time_stamp)), units = "secs") |> lubridate::seconds_to_period()
  else
    return(NA)
}


get_analysis_interruptions <- function(data, break_mins){
  if (check_data_present(data))
    sum(diff(unique(data@dataset$acquisition_time_stamp), units = "secs") > break_mins)
  else
    return(NA)
}

check_var_in_dataset <- function(table, variable) {
  if(variable == "feature_conc" & !"feature_conc" %in% names(table)) cli_abort("Concentrations not available, please process data or choose another variable.", show = "none", parent = NULL, call= NULL)
  if(variable == "feature_area" & !"feature_area" %in% names(table)) cli_abort("Area is not available, please choose another variable.", show = "none", parent = NULL, call= NULL)
  if(variable == "response" & !"response" %in% names(table)) cli_abort("Response is not available, please choose another variable.", show = "none", parent = NULL, call= NULL)
  if(variable == "feature_norm_intensity" & !"feature_norm_intensity" %in% names(table)) cli_abort("Normalized intensities not available, please process data, or choose another variable.", show = "none", parent = NULL, call= NULL)
  if(variable == "height" & !"height" %in% names(table)) cli_abort("Heights not available, please choose another variable.", show = "none", parent = NULL, call= NULL)
  if(variable == "conc_raw" & !"conc_raw" %in% names(table)) cli_abort("Concentrations not available, please process data, or choose another variable.", show = "none", parent = NULL, call= NULL)
}

#' @title Get the start and end analysis numbers of specified batches
#' @description
#' Sets the analysis order (sequence) based on either (i) analysis timestamp, if available, (ii) the order in which analysis appeared in the imported raw data file, or (iii) the order in which analyses were defined in the Analysis metadata.
#' @param data MidarExperiment object
#' @param batch_ids A vector with one or two elements: the first and/or last batch ID. If NULL or invalid, the function will abort.
#' @return A vector with two elements: the lower and upper analysis number for the specified batch(es).
#' @export
get_batch_boundaries <- function(data = NULL, batch_ids = NULL) {
  check_data(data)
  if (nrow(data@annot_batches) == 0) {
    cli::cli_abort("No batches defined in the dataset.")
  }

  if (is.null(batch_ids) || length(batch_ids) == 0) {
    cli::cli_abort("No batch IDs provided.")
  }
  # Ensure batch_ids has at least one value
  if (length(batch_ids) == 1) {
    first_batch_id <- last_batch_id <- batch_ids[1]
  } else if (length(batch_ids) == 2) {
    first_batch_id <- batch_ids[1]
    last_batch_id <- batch_ids[2]
  } else {
    cli::cli_abort("Please provide a vector with one or two batch IDs.")
  }

  # Check if the specified batch IDs exist in the dataset
  existing_batches <- data@annot_batches$batch_id
  if (!all(c(first_batch_id, last_batch_id) %in% existing_batches)) {
    cli::cli_abort("One or more of the specified batch IDs do not exist in the dataset.")
  }

  # Filter the dataset for the given range of batch IDs
  d <- data@annot_batches |>
    filter(.data$batch_id >= first_batch_id & .data$batch_id <= last_batch_id)

  if (nrow(d) == 0) {
    cli::cli_abort("No batches found for the provided range of batch IDs.")
  }

  # Extract lower and upper analysis numbers
  lower_bound <- min(d$id_batch_start)
  upper_bound <- max(d$id_batch_end)

  return(c(lower_bound, upper_bound))
}


#' @title Set the analysis order
#' @description
#' Sets the analysis order (sequence), based on either (i) analysis timestamp, if available, (ii) the order in which analysis appeared in the imported raw data file, or (iii) the order in which analyses were defined in the Analysis metadata
#' @param data MidarExperiment object
#' @param analysis_sequence Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
#' @return MidarExperiment object
#' @examples
#' file_path <- system.file("extdata", "sPerfect_MRMkit.tsv", package = "midar")
#' mexp <- MidarExperiment()
#' mexp <- data_import_mrmkit(mexp, path = file_path, use_metadata = TRUE)
#' mexp <- set_analysis_order(mexp, "timestamp")

#' @export

#c("timestamp", "resultfile", "metadata")
set_analysis_order <- function(data = NULL, analysis_sequence = "default"){
  check_data(data)

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
        mutate(run_seq_num = row_number()) |>
        select(-"acquisition_time_stamp")

      data@annot_analyses <- data@annot_analyses |>
        inner_join(d_temp,  by = "analysis_id")
    } else {
      cli::cli_abort(call. = FALSE, "No acquisition timestamp field present in analysis results, please set parameter `analysis_sequence` to `resultfile` or `metadata`.")
    }
  } else if (analysis_sequence == "resultfile") {
    d_temp <- d_temp |> mutate(run_seq_num = row_number())
    data@annot_analyses <- data@annot_analyses |>
      inner_join(d_temp,  by = "analysis_id")
  } else if (analysis_sequence == "metadata") {
      data@annot_analyses <- data@annot_analyses |>  mutate(run_seq_num = row_number())
  }
  data
}


# # ##### TODO TODO =====================
# d_dataset <- data@dataset_orig |>
#   dplyr::inner_join(
#     data@annot_analyses |>
#       dplyr::select("run_seq_num", "analysis_id", "qc_type", "specimen", "sample_id", "replicate_no", "valid_analysis", "batch_id"),
#     by = c("analysis_id")) |>
#   dplyr::inner_join(
#     metadata$annot_features |>
#       filter(.data$valid_feature) |>
#       dplyr::select(dplyr::any_of(c("feature_id", "feature_id", "feature_class", "istd_feature_id", "quant_istd_feature_id", "is_istd", "feature_id", "is_quantifier", "valid_feature", "table_templates", "interference_feature_id", "interference_proportion"))),
#     by = c("feature_id"), keep = FALSE
#   )



# Link DATA with METADATA and create DATASET table. =================
## - Only valid analyses and features will be added
## - Only key information will be added
link_data_metadata <- function(data = NULL, minimal_info = TRUE){
  check_data(data)
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
    filter(.data$valid_analysis, .data$valid_feature) |>
    select(
      "run_seq_num",
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

  data@is_isotope_corr <- FALSE
  #data@is_istd_normalized <- FALSE
  #data@is_quantitated <- FALSE
  data@is_drift_corrected <- FALSE
  data@is_batch_corrected <- FALSE
  data@is_filtered <- FALSE


  data@status_processing <- glue::glue("Annotated raw {toupper(str_remove(data@feature_intensity_var, 'feature_'))} values")

  data@metrics_qc <- data@metrics_qc[FALSE,]

  # Arrange run_seq_num and then by feature_id, as they appear in the metadata
  # TODO:  check if as intended
  data@dataset <- data@dataset |>
    dplyr::arrange(match(.data$feature_id, data@annot_features$feature_id)) |>
    arrange(.data$run_seq_num)

  # TODOTODO: DECIDE WHEN WHERE TO RUN THIS
  # check_integrity(data, excl_unannotated_analyses = excl_unannotated_analyses)

  data
}

#' @title Set default variable to be used as feature raw signal value
#' @description
#' Sets the raw signal variable used for calculations starting from raw signal
#' values (i.e., normalization)
#' @param data MidarExperiment object
#' @param variable_name Feature variable to be used as default feature intensity for downstream processing.
#' @param auto_select If `TRUE` then the first available of these will be used as default: "feature_intensity", "feature_response", "feature_area", "feature_height".
#' @param warnings Suppress warnings
#' @param ... Feature variables to best search for one-by-one when `auto-detect = TRUE`
#' @return MidarExperiment object
#' @export

data_set_intensity_var <- function(data = NULL, variable_name, auto_select = FALSE, warnings = TRUE, ...){
  check_data(data)

  if (auto_select) {
    var_list <- unlist(rlang::list2(...), use.names = FALSE)
    id_all <- match(var_list,names(data@dataset_orig))
    idx <- which(!is.na(id_all))[1]
    if (!is.na(idx)) {
      data@feature_intensity_var = var_list[idx]
      cli_alert_info(text = cli::col_grey("{.var {var_list[idx]}} selected as default feature intensity. Modify with {.fn data_set_intensity_var}."))
      variable_name <- var_list[idx]
      } else {
      cli_alert_warning(text = cli::col_yellow("No typical feature intensity variable found in the data. Use {.fn data_set_intensity_var} to set it.}}."))
      return(data)
      }
  } else {
    #TODO: Double check behavior if there a feature_intensity in the raw data file
    if (! variable_name %in% names(data@dataset_orig))
      cli_abort(c("x" = "{.var variable_name} is not present in the raw data."))

    if (! variable_name %in% c("feature_intensity", "feature_response", "feature_area", "feature_height"))
      if(warnings) cli_alert_warning(text = "{.var {variable_name}} is not a typically used raw signal name (i.e., area, response, intensity, height).")
  }
  data@feature_intensity_var <- variable_name

  if (check_dataset_present(data)) {
    calc_cols <- c("featue_norm_intensity", "feature_conc", "feature_amount", "feature_raw_conc")
    if (any(calc_cols %in% names(data@dataset))){
      data@dataset <- data@dataset |> select(-any_of(calc_cols))
      cli_alert_info(cli::col_yellow("New feature intensity variable defined, please reprocess data."))
    } else
    {
      cli_alert_success(cli::col_green("`{variable_name}` was set as default feature intensity variable for downstream processing."))
    }
    data <- link_data_metadata(data)
  }

  if(variable_name == "feature_conc") data@is_quantitated <- TRUE
  data
}


#'  @title Exclude analyses from the dataset
#' @param data MidarExperiment object
#' @param analyses_exlude Vector of analysis IDs (case-sensitive) to exclude from the dataset.
#' @param overwrite If `TRUE` then existing valid_analysis flags will be overwritten, otherwise appended
#' @return `MidarExperiment` object
#' @export

data_exclude_analyses <- function(data = NULL, analyses_exlude, overwrite ){
  check_data(data)
  if (all(is.na(analyses_exlude)) | length(analyses_exlude) == 0) {
    if(!overwrite){
      cli_abort(cli::col_red("No `analysis id`s provided. To (re)include all analyses, use `analysis_ids_exlude = NA` and `overwrite = TRUE`."))
    } else{
      cli::cli_alert_info(cli::col_green("Exclusions were removed and all analyses are now included for subsequent steps. Please reprocess data."))
      data@analyses_excluded <- FALSE
      data@annot_analyses <- data@annot_analyses |> mutate(valid_analysis = TRUE)
      data <- link_data_metadata(data)
      return(data)
      }
  }
  if (any(!c(analyses_exlude) %in% data@annot_analyses$analysis_id) > 0) {
    cli_abort(cli::col_red("One or more provided `analysis id`s to exclude are not present. Please check the analysis metadata."))
  }
  if(!overwrite){
    data@annot_analyses <- data@annot_analyses |>
      mutate(valid_analysis = !(.data$analysis_id %in% analyses_exlude) & .data$valid_analysis)
    cli_alert_info(cli::col_green("A total of {data@annot_analyses |> filter(!.data$valid_analysis) |> nrow()} analyses were now excluded for downstream processing. Please reprocess data."))
    data@analyses_excluded <- TRUE
  } else {
    data@annot_analyses <- data@annot_analyses |>
      mutate(valid_analysis = !(.data$analysis_id %in% analyses_exlude))
    cli_alert_info(cli::col_green("{data@annot_analyses |> filter(!.data$valid_analysis) |> nrow()} analyses were excluded for downstream processing. Please reprocess data."))
    data@analyses_excluded <- TRUE
    }

  data <- link_data_metadata(data)

  data
}



#' @title Get the annotated or the originally imported analytical data
#' @param data MidarExperiment object
#' @param original Boolean indicating whether to return the original imported data (`TRUE`) or the annotated data (`FALSE`)
#' @param overwrite If `TRUE` then existing valid_feature flags will be overwritten, otherwise appended
#' @return A tibble with the analytical data in the long format
#' @export

data_get_analyticaldata <- function(data = NULL, original = FALSE){
  check_data(data)
  if (original) {
    return(data@dataset_orig)
  } else {
    return(data@dataset)
  }
}

#' @title Exclude features from the dataset
#' @param data MidarExperiment object
#' @param features_exlude Vector of feature IDs (case-sensitive) to exclude from the dataset
#' @param overwrite If `TRUE` then existing valid_feature flags will be overwritten, otherwise appended
#' @return `MidarExperiment` object
#' @export

data_exclude_features <- function(data = NULL, features_exlude, overwrite ){
  check_data(data)

  if (all(is.na(features_exlude)) | length(features_exlude) == 0) {
    if(!overwrite){
      cli_abort(cli::col_red("No `feature id`s provided. To include all analyses, use `feature_ids_exlude = NA` and `overwrite = TRUE`."))
    } else{
      cli::cli_alert_info(cli::col_green("All exlusions were removed, i.e. all analyses are included. Please reprocess data."))
      return(data)
    }
  }
  if (any(!c(features_exlude) %in% data@annot_features$feature_id) > 0) {
    cli_abort(cli::col_red("One or more provided `feature id`s to exclude are not present. Please check the feature metadata."))
  }
  if(!overwrite){
    data@annot_features <- data@annot_features |>
      mutate(valid_feature = !(.data$feature_id %in% features_exlude) & .data$valid_feature)
    cli_alert_info(cli::col_green("A total of {data@annot_features |> filter(!.data$valid_feature) |> nrow()} features were now excluded for downstream processing. Please reprocess data."))
  }
  else {
    data@annot_features <- data@annot_features |>
      mutate(valid_feature = !(.data$feature_id %in% features_exlude))
    cli_alert_info(cli::col_green("{data@annot_features |> filter(!.data$valid_feature) |> nrow()} analyses were excluded for downstream processing. Please reprocess data."))
  }

  data <- link_data_metadata(data)

  data
}



