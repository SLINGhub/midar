
# TODO: This below has to be changed, mapping of metadata to data should be a distinct function
#' @title Import metadata from the MSOrganizer template (.XLM) and associates it with the analysis data
#' @description Requires version 1.9.1 of the template
#' @param data MidarExperiment object
#' @param path file name and path
#' @param analysis_sequence Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
#' @param excl_unannotated_analyses Exclude analyses (samples) that have no matching metadata
#' @return MidarExperiment object
#' @export
#'

import_metadata_msorganizer <- function(data, path, analysis_sequence = "default", excl_unannotated_analyses = FALSE) {

  analysis_sequence <- rlang::arg_match(arg = analysis_sequence, c("timestamp", "resultfile", "metadata", "default"), multiple = FALSE)

  d_annot <- read_msorganizer_xlm(path)

  #browser()
  data@annot_analyses <- data@dataset_orig %>%
    dplyr::select("raw_data_filename", dplyr::any_of("acquisition_time_stamp") ) %>%
    dplyr::distinct() %>%
    dplyr::right_join(d_annot$annot_analyses, by = c("raw_data_filename" = "raw_data_filename"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_analyses_template)

  if("acquisition_time_stamp" %in% names(data@annot_analyses) & (analysis_sequence %in% c("timestamp", "default"))){
        data@annot_analyses <- data@annot_analyses |> arrange(.data$acquisition_time_stamp)
  } else {
    if(analysis_sequence == "resultfile")
      data@annot_analyses <- data@annot_analyses |> dplyr::arrange(match(.data$analysis_id, d_annot$dataset_orig$analysis_id |> unique()))
    else if(analysis_sequence == "metadata")
      data@annot_analyses <- data@annot_analyses |> dplyr::arrange(match(.data$analysis_id, d_annot$annot_analyses$analysis_id))
    else if(analysis_sequence == "timestamp")
      stop(call. = FALSE, "No acquisition timestamp field present in analysis results, please set parameter `analysis_sequence` to `resultfile` or `metadata` to define analysis order.")
    else
      writeLines(crayon::yellow(glue::glue("\u2713 No acquisition timestamps present in results, order was therefore based on analysis results sequence. Set parameter `analysis_sequence` to `metadata` to use this sequence as analysis order.")))

    }

  data@annot_analyses <- data@annot_analyses |>
    mutate(run_id = row_number(), .before = 1)

  data@annot_batches <- data@annot_analyses %>%
    dplyr::group_by(.data$batch_no) %>%
    dplyr::summarise(
      batch_id = .data$batch_id[1],
      id_batch_start = dplyr::first(.data$run_id),
      id_batch_end = dplyr::last(.data$run_id)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$id_batch_start)%>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_batch_info_template)

  data@annot_features <- data@dataset_orig %>%
    dplyr::select("feature_name") %>%
    dplyr::distinct() %>%
    dplyr::right_join(d_annot$annot_features, by = c("feature_name"="feature_name"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_features_template)

  data@annot_istd <- data@annot_features %>%
    dplyr::select("quant_istd_feature_name") %>%
    dplyr::distinct() %>%
    dplyr::left_join(d_annot$annot_istd, by = c("quant_istd_feature_name"="quant_istd_feature_name"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_istd_template)


  data@annot_responsecurves <- data@annot_analyses %>%
    dplyr::filter(.data$qc_type == "RQC") %>%
    dplyr::select("analysis_id", "raw_data_filename") %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::right_join(d_annot$annot_responsecurves, by = c("raw_data_filename"="raw_data_filename"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_responsecurves_template)


  data@dataset_orig <- data@dataset_orig %>%
    dplyr::bind_rows(pkg.env$dataset_templates$dataset_orig_template)


  d_dataset <- data@dataset_orig  %>%
    dplyr::inner_join(data@annot_analyses  %>%
                        dplyr::select("run_id", "analysis_id", "raw_data_filename", "qc_type", "specimen" ,"sample_id", "replicate_no", "valid_analysis", "batch_id"), by = c("raw_data_filename")) %>%
    dplyr::inner_join(d_annot$annot_features %>%
                        filter(.data$valid_integration) |>
                        dplyr::select(dplyr::any_of(c("feature_name", "feature_name", "feature_class", "norm_istd_feature_name", "quant_istd_feature_name", "is_istd", "feature_name", "is_quantifier", "valid_integration", "feature_response_factor", "interference_feature_name", "interference_proportion"))),
                      by = c("feature_name"), keep = FALSE)

  data@dataset <-
    dplyr::bind_rows(pkg.env$dataset_templates$dataset_template, d_dataset) |>
    mutate(corrected_interference = FALSE,
           outlier_technical = FALSE) |>
    dplyr::arrange(match(.data$feature_name, d_annot$annot_features$feature_name))

  data@dataset <- data@dataset |> arrange(.data$run_id)



  #stopifnot(methods::validObject(data, excl_nonannotated_analyses))
  check_integrity(data, excl_unannotated_analyses = excl_unannotated_analyses)
  data@status_processing <- "Annotated Raw Data"


  writeLines(crayon::green(glue::glue("\u2713 Metadata successfully associated with {length(data@dataset$analysis_id %>% unique())} samples and {length(data@dataset$feature_name %>% unique())} features.")))
  data
}



#' @title Reads and parses metadata provided by the MSOrganizer EXCEL  template.
#' @description Requires version 1.9.1 of the template
#'
#' @param path File path of the MSOrganizer EXCEL template (*.xlm)
#' @param trim_ws Trim all white spaces and double spaces
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom tidyselect vars_select_helpers
#' @importFrom dplyr select mutate filter group_by row_number
#' @importFrom stringr regex
#' @return A list of tibbles with different metadata
read_msorganizer_xlm <- function(path, trim_ws = TRUE){
  d_annot <- list()

  # ANALYSIS/SAMPLE annotation
  # ToDo: Make note if feature names are not original

  d_temp_analyses <- readxl::read_excel(path, sheet = "Analyses (Samples)", trim_ws = TRUE)
  names(d_temp_analyses) <- tolower(names(d_temp_analyses))

  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "analysis_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "valid_analysis", init_value = TRUE, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "replicate_no", init_value = 1L, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "specimen", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "panel_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "sample_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)

  # TODO: If analysis_id is defined, then it will overwrite the raw_data_filename defined in the raw data files
  # TODO: if user-defined analysis_id names (=remapping if IDs) are provided, then it should be reported somewhere,  possible source of user-error!

  d_annot$annot_analyses <- d_temp_analyses |>
    dplyr::mutate(
      valid_analysis = TRUE,
      batch_id = as.character(.data$batch_id),
      analysis_no = dplyr::row_number()) |>
    dplyr::select(
      "analysis_no",
      "analysis_id",
      "raw_data_filename" ,
      qc_type = "sample_type",
      "sample_amount",
      "sample_amount_unit",
      istd_volume ="istd_mixture_volume_[ul]",
      "batch_id" ,
      "replicate_no" ,
      "valid_analysis",
      "specimen",
      "sample_id",
      "remarks"
    ) |>
    dplyr::mutate(batch_no = dplyr::cur_group_id(), .by = c("batch_id")) |>
    dplyr::mutate(
      raw_data_filename = stringr::str_squish(as.character(.data$raw_data_filename)),
      raw_data_filename = stringr::str_remove(.data$raw_data_filename, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)) ,
      analysis_id = stringr::str_remove(.data$analysis_id, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)) ,
      analysis_id = stringr::str_squish(as.character(.data$analysis_id)),
      analysis_id = if_else(is.na(.data$analysis_id), .data$raw_data_filename, .data$analysis_id),
      specimen = stringr::str_squish(as.character(.data$specimen)),
      valid_analysis = as.logical(.data$valid_analysis),
      qc_type = if_else(.data$qc_type == "Sample" | is.na(.data$qc_type), "SPL", .data$qc_type))|>
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))


  # FEATURE annotation

  # ToDo: Make note if feature names are not original
  d_temp_features <- readxl::read_excel(path, sheet = "Features (Analytes)", trim_ws = TRUE)
  names(d_temp_features) <- tolower(names(d_temp_features))

  d_temp_features <- d_temp_features |> add_missing_column(col_name = "feature_class", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "quantifier", init_value = TRUE, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "valid_integration", init_value = TRUE, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "response_factor", init_value = 1, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "new_feature_name", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "interference_feature_name", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "interference_proportion", init_value = NA_real_, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)



  # NOTE: If feature_name is defined, then it will overwrite the feature name defined in the raw data files
  # Todo: if user-defined feature names are provided, then it should be reported somewhere,  possible source of user-error!

  d_annot$annot_features <- d_temp_features |>
    dplyr::mutate(
      new_feature_name = stringr::str_squish(.data$new_feature_name),
      feature_name = stringr::str_squish(.data$feature_name),
      feature_name = if_else(is.na(.data$new_feature_name), .data$feature_name, .data$new_feature_name),
      feature_class = stringr::str_squish(.data$feature_class),
      norm_istd_feature_name	= stringr::str_squish(.data$istd_feature_name),
      quant_istd_feature_name = stringr::str_squish(.data$istd_feature_name),
      is_istd = (.data$feature_name == .data$norm_istd_feature_name),
      is_quantifier = if_else(tolower(.data$quantifier) %in% c("yes","true"), TRUE, FALSE),
      is_quantifier = as.logical(.data$is_quantifier),
      interference_feature_name = stringr::str_squish(.data$interference_feature_name),
      remarks = NA_character_) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish)) %>%
    dplyr::select(
      "feature_name",
      "feature_class",
      "is_istd",
      "norm_istd_feature_name",
      "quant_istd_feature_name",
      feature_response_factor = "response_factor",
      "is_quantifier",
      "valid_integration",
      "interference_feature_name",
      "interference_proportion",
      "remarks")

  #ToDo: Merged cell in template
  annot_istd <- readxl::read_excel(path,
                           sheet = "Internal Standards",
                           trim_ws = TRUE, .name_repair = ~ ifelse(nzchar(.x), .x,LETTERS [seq_along(.x)]))

  names(annot_istd) <- tolower(names(annot_istd))
  names(annot_istd)[1] <- "istd_feature_name"

  d_annot$annot_istd <- annot_istd |>
    #dplyr::mutate(istd_compound_name = na_character_) |>
    dplyr::select(
      quant_istd_feature_name = "istd_feature_name",
      istd_conc_nmolar = "istd_conc_[nm]") %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  names(annot_istd) <- tolower(names(annot_istd))

  d_annot_responsecurves <- readxl::read_excel(path, sheet = "Response Curves")

  names(d_annot_responsecurves) <- tolower(names(d_annot_responsecurves))



  d_annot$annot_responsecurves <- d_annot_responsecurves |>
    dplyr::select(
      "raw_data_filename",
      rqc_series_id = "response_curve_name",
      "relative_sample_amount" = "relative_sample_amount_[%]",
      "injection_volumne" = "injection_volume_[ul]") |>
    dplyr::mutate(
      raw_data_filename = stringr::str_remove(.data$raw_data_filename, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)),
      raw_data_filename = stringr::str_squish(as.character(.data$raw_data_filename)),
      rqc_series_id = stringr::str_squish(as.character(.data$rqc_series_id)),
      relative_sample_amount = .data$relative_sample_amount/100) %>%
      dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  return(d_annot)
}
