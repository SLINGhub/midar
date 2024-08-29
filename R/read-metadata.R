
#' @title Retrieve available metadata from imported analysis data
#' @description Available information will depend on the format of imported raw analysis data. See [rawdata_import_agilent()] and [rawdata_import_mrmkit()]
#' @param data MidarExperiment object
#' @param analysis_sequence Must by any of: "resultfile" or "timestamp". Defines how the analysis order is according the sequence (default) in the loaded analysis data or the timestamp, if available, in the analysis data.
#' @return MidarExperiment object
#' @export
#'

metadata_from_data<- function(data, analysis_sequence = "resultfile") {
  analysis_sequence <- rlang::arg_match(arg = analysis_sequence, c("timestamp", "resultfile"))

  # get analysis metadata
   data@annot_analyses <- data@dataset_orig %>%
    dplyr::select("analysis_id", dplyr::any_of(c("acquisition_time_stamp", "qc_type", "batch_no"))) %>%
    dplyr::distinct() %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_analyses_template)

  if ("acquisition_time_stamp" %in% names(data@dataset_orig) & (analysis_sequence %in% c("timestamp", "default"))) {
    data@annot_analyses <- data@annot_analyses |> arrange(.data$acquisition_time_stamp)
  } else {
    if (analysis_sequence == "timestamp") {
      cli::cli_abort(call. = FALSE, "No acquisition timestamp field present in analysis results, please set parameter `analysis_sequence` to `resultfile`.")
    }
  }

  data@annot_analyses <- data@annot_analyses |>
    mutate(run_id = row_number(), .before = 1)

  # retrieve batches
  if ("batch_id" %in% names(data@dataset_orig)) data <- load_metadata_batches(data)



  # get feature metadata

  data@annot_features <- data@dataset_orig %>%
    dplyr::select("feature_id", dplyr::any_of(c("feature_class", "is_istd", "precursor_mz", "product_mz", "norm_istd_feature"))) %>%
    dplyr::distinct() %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_features_template)

  # TODO
  data@dataset <- data@dataset_orig


  data
}

# Retreive batch info from analysis metadata

load_metadata_batches <- function(data){
  data@annot_batches <- data@annot_analyses %>%
    dplyr::group_by(.data$batch_no) %>%
    dplyr::summarise(
      batch_id = .data$batch_id[1],
      id_batch_start = dplyr::first(.data$run_id),
      id_batch_end = dplyr::last(.data$run_id)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$id_batch_start) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_batch_info_template)
  data
}


# TODO: This below has to be changed, mapping of metadata to data should be a distinct function
#' @title Import metadata from the midar template (.XLM) and associates it with the analysis data
#' @description Requires version 1.9.1 of the template. Template is based on the MSorganizer template (see TODO)
#' @param data MidarExperiment object
#' @param path file name and path
#' @param analysis_sequence Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
#' @param excl_unannotated_analyses Exclude analyses (samples) that have no matching metadata
#' @return MidarExperiment object
#' @export
#'

metadata_import_midarxlm<- function(data, path, analysis_sequence = "default", excl_unannotated_analyses = FALSE) {

  analysis_sequence <- rlang::arg_match(arg = analysis_sequence, c("timestamp", "resultfile", "metadata", "default"), multiple = FALSE)
  d_annot <- read_metadata_midarxlm(path)

  add_metadata(data, d_annot, analysis_sequence, excl_unannotated_analyses)

}


assert_warnings_metadata <- function(list_of_errors, data=NULL, warn = TRUE, ...) {
  #browser()
  res <- as_tibble(do.call(rbind, list_of_errors)) |>
    select(message, description, num.violations) |>
    mutate(message = str_replace(message, "analysis_id", "'analysis_id'")) |>
    mutate(Field = str_extract(message, "(?<=\\').+?(?=\\')")) |>
    separate(col = description, into = c("Issue", "Table"), sep = ";", remove = TRUE) |>
    select(Table, Field, Issue, Count = num.violations)

  res$Count <- res$Count |> unlist()
  res
}


#' @title Add metadata an MidarExperiment object
#' @description Metadata provided as a list of tibbles will validates for consistency again loaded analysis data of the provided MidarExperiment object and then transfered.
#' @param data MidarExperiment object
#' @param path List of tibbles or data.frames containing analysis, feature, istd, response curve tables
#' @param analysis_sequence Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
#' @param excl_unannotated_analyses Exclude analyses (samples) that have no matching metadata
#' @return MidarExperiment object
#' @export
#'
 add_metadata <- function(data, d_annot, analysis_sequence, excl_unannotated_analyses, allow_na = FALSE) {
  analysis_sequence <- rlang::arg_match(arg = analysis_sequence, c("timestamp", "resultfile", "metadata", "default"), multiple = FALSE)


  t_defects <- tibble()
  t_warnings <- tibble()

  if (!is.null(d_annot$annot_analyses) && nrow(d_annot$annot_analyses) > 0){
    d_analyses <- d_annot$annot_analyses

    #browser()

    d_analyses <- d_analyses |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::verify(assertr::has_all_names("analysis_id"), obligatory=TRUE, description = "Column missing;Analyses/Samples") |>
      assertr::verify(assertr::has_all_names("qc_type"), obligatory=TRUE, description = "Column missing;Analyses/Samples") |>
      assertr::chain_end(error_fun = error_append)

    if(rlang::is_true(attributes(d_analyses)$assertr_data_defective))
      cli::cli_abort(message = c(
        "x" = "Analyses metadata must contain the columns `analysis_id` and `qc_type`",
        "Ensure both columns are defined in the Analyses metadata and re-import"))


    d_analyses <- d_analyses |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert({\(x) all(is_na(x))}, all_of(c("qc_type")), description = "Missing value(s);Analyses/Samples") |>
      assertr::verify(all(assertr::is_uniq(analysis_id)), description = "Analysis IDs duplicated;Analyses/Samples") |>
      assertr::verify((analysis_id %in% unique(mexp@dataset_orig$analysis_id)), description = "Analysis not present in analysis data;Analyses/Samples") |>
      assertr::verify((unique(mexp@dataset_orig$analysis_id) %in% analysis_id), description = "Unannotated analyses found analysis data;Analyses/Samples") |>
      assertr::assert(\(x){not_na(x)}, any_of(c("valid_analysis", "sample_amount", "sample_amount_unit", "istd_volume", "batch_id", "replicate_no")), description = "Missing value;Analyses/Samples") |>
      assertr::chain_end(error_fun = warning_append)

    if (!is.null(attr(d_analyses, "assertr_errors"))) {
      cli::cli_alert_danger(text = cli::col_yellow("Analysis metadata contains following issues"))
      print(assert_warnings_metadata(attr(d_analyses, "assertr_errors")))
      cli::cli_abort(message = cli::col_red(c("i" = "Please check the analysis metadata or set the parameter `ignore_warnings` to TRUE")))
    }

    data@annot_analyses <- data@dataset_orig %>%
      dplyr::select("analysis_id",  dplyr::any_of("acquisition_time_stamp")) %>%
      dplyr::distinct() %>%
      dplyr::right_join(d_analyses, by = c("analysis_id"), keep = FALSE) %>%
      dplyr::bind_rows(pkg.env$dataset_templates$annot_analyses_template)



    if ("acquisition_time_stamp" %in% names(data@annot_analyses) && (analysis_sequence %in% c("timestamp", "default"))) {
      data@annot_analyses <- data@annot_analyses |> arrange(.data$acquisition_time_stamp)
    } else {
      if (analysis_sequence == "resultfile") {
        data@annot_analyses <- data@annot_analyses |> dplyr::arrange(match(.data$analysis_id, d_annot$dataset_orig$analysis_id |> unique()))
      } else if (analysis_sequence == "metadata") {
        data@annot_analyses <- data@annot_analyses |> dplyr::arrange(match(.data$analysis_id, d_annot$annot_analyses$analysis_id))
      } else if (analysis_sequence == "timestamp") {
        cli::cli_abort(call. = FALSE, "No acquisition timestamp field present in analysis results, please set parameter `analysis_sequence` to `resultfile` or `metadata` to define analysis order.")
      } else {
        cli::cli_alert_warning(cli::col_yellow(glue::glue("No acquisition timestamps present in results, order was therefore based on analysis results sequence. Set parameter `analysis_sequence` to `metadata` to use this sequence as analysis order.")))
      }
    }
  }

  data@annot_analyses <- data@annot_analyses |>
    mutate(run_id = row_number(), .before = 1)

  data@annot_batches <- data@annot_analyses %>%
    dplyr::group_by(.data$batch_no) %>%
    dplyr::summarise(
      batch_id = .data$batch_id[1],
      id_batch_start = dplyr::first(.data$run_id),
      id_batch_end = dplyr::last(.data$run_id)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$id_batch_start) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_batch_info_template)

  data@annot_features <- data@dataset_orig %>%
    dplyr::select("feature_id") %>%
    dplyr::distinct() %>%
    dplyr::right_join(d_annot$annot_features, by = c("feature_id" = "feature_id"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_features_template)

  data@annot_istd <- data@annot_features %>%
    dplyr::select("quant_istd_feature_id") %>%
    dplyr::distinct() %>%
    dplyr::left_join(d_annot$annot_istd, by = c("quant_istd_feature_id" = "quant_istd_feature_id"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_istd_template)


  data@annot_responsecurves <- data@annot_analyses %>%
    dplyr::filter(.data$qc_type == "RQC") %>%
    dplyr::select("analysis_id") %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::right_join(d_annot$annot_responsecurves, by = c("analysis_id" = "analysis_id"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_responsecurves_template)


  data@dataset_orig <- data@dataset_orig %>%
    dplyr::bind_rows(pkg.env$dataset_templates$dataset_orig_template)

  #TODO
  d_dataset <- data@dataset_orig %>%
    dplyr::inner_join(
      data@annot_analyses %>%
        dplyr::select("run_id", "analysis_id", "qc_type", "specimen", "sample_id", "replicate_no", "valid_analysis", "batch_id"),
      by = c("analysis_id")) %>%
    dplyr::inner_join(
      d_annot$annot_features %>%
        filter(.data$valid_integration) |>
        dplyr::select(dplyr::any_of(c("feature_id", "feature_id", "feature_class", "norm_istd_feature_id", "quant_istd_feature_id", "is_istd", "feature_id", "is_quantifier", "valid_integration", "feature_response_factor", "interference_feature_id", "interference_proportion"))),
      by = c("feature_id"), keep = FALSE
    )

  data@dataset <-
    dplyr::bind_rows(pkg.env$dataset_templates$dataset_template, d_dataset) |>
    mutate(
      corrected_interference = FALSE,
      outlier_technical = FALSE
    ) |>
    dplyr::arrange(match(.data$feature_id, d_annot$annot_features$feature_id))

  data@dataset <- data@dataset |> arrange(.data$run_id)



  # stopifnot(methods::validObject(data, excl_nonannotated_analyses))
  check_integrity(data, excl_unannotated_analyses = excl_unannotated_analyses)
  data@status_processing <- "Annotated Raw Data"


  cli_alert_success(col_green(glue::glue("Metadata successfully associated with {length(data@dataset$analysis_id %>% unique())} samples and {length(data@dataset$feature_id %>% unique())} features.")))
  data
}



#' @title Reads and parses metadata provided by the MSOrganizer EXCEL  template.
#' @description Requires version 1.9.1 of the template
#' NOTES
#' - if no sample_type is defined then SPL will be assigned
#' - if valid_analysis is left blank for all analyses then samples then it will be replace by TRUE
#' - if valid_analysis is undefined for one or more, but all all samples, then an error will be returned
#' @param path File path of the MSOrganizer EXCEL template (*.xlm)
#' @param trim_ws Trim all white spaces and double spaces
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom tidyselect vars_select_helpers
#' @importFrom dplyr select mutate filter group_by row_number
#' @importFrom stringr regex
#' @return A list with tibbles containing different metadata
read_metadata_midarxlm <- function(path, trim_ws = TRUE) {

  d_annot <- list()

  # ANALYSIS/SAMPLE annotation --------
  d_temp_analyses<- readxl::read_excel(path, sheet = "Analyses (Samples)", trim_ws = TRUE)
  names(d_temp_analyses) <- tolower(names(d_temp_analyses))

  #d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "analysis_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "valid_analysis", init_value = TRUE, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "replicate_no", init_value = 1L, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "specimen", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "panel_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "sample_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)

  d_annot$annot_analyses <- d_temp_analyses |>
    select(
      "analysis_id" = "raw_data_filename",
      qc_type = "sample_type",
      "sample_amount",
      "sample_amount_unit",
      istd_volume = "istd_mixture_volume_[ul]",
      "batch_id",
      "replicate_no",
      "specimen",
      "sample_id",
      "valid_analysis",
      "remarks"
    ) |>
    mutate(analysis_no = dplyr::row_number(), .before = 1) |>
    mutate(
      batch_id = as.character(.data$batch_id),
      analysis_id = stringr::str_squish(as.character(.data$analysis_id)),
      analysis_id = stringr::str_remove(.data$analysis_id, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)),
      valid_analysis = as.logical(case_match(tolower(.data$valid_analysis),
                                  "yes" ~ TRUE,
                                  "no"~ FALSE,
                                  "true" ~ TRUE,
                                  "false" ~ FALSE,
                                  .default = NA)),
      qc_type = if_else(.data$qc_type == "Sample" | is.na(.data$qc_type), "SPL", .data$qc_type)) |>
    mutate(batch_no = dplyr::cur_group_id(), .by = c("batch_id"), .before = batch_id) |>
    mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish)) |>
    ungroup()

    # Handle the non-mandatory field Valid_Analysis
    if (all(is.na(d_annot$annot_analyses$valid_analysis)))
      d_annot$annot_analyses$valid_analysis <- TRUE
    else
      if (any(is.na(d_annot$annot_analyses$valid_analysis)))
        cli::cli_abort("`Valid_Analysis` is not defined for one or more analyses/samples. Please check sheet 'Analyses (Samples)'")


  # FEATURE annotation  -------------------------

  # ToDo: Make note if feature names are not original
  d_temp_features <- readxl::read_excel(path, sheet = "Features (Analytes)", trim_ws = TRUE)
  names(d_temp_features) <- tolower(names(d_temp_features))

  d_temp_features <- d_temp_features |> add_missing_column(col_name = "feature_class", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "quantifier", init_value = TRUE, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "valid_integration", init_value = TRUE, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "response_factor", init_value = 1, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "new_feature_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "interference_feature_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "interference_proportion", init_value = NA_real_, make_lowercase = FALSE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)



  # TODO: feauture_id and (new)_feature_name...find a clear way
  # feature_label... idea from

  d_annot$annot_features <- d_temp_features |>
    dplyr::mutate(
      feature_id = stringr::str_squish(.data$feature_name),
      feature_label = stringr::str_squish(.data$new_feature_name),
      feature_class = stringr::str_squish(.data$feature_class),
      norm_istd_feature_id = stringr::str_squish(.data$istd_feature_name),
      quant_istd_feature_id = stringr::str_squish(.data$istd_feature_name),
      is_istd = (.data$feature_id == .data$norm_istd_feature_id),
      is_quantifier = as.logical(case_match(tolower(.data$quantifier),
                                             "yes" ~ TRUE,
                                             "no"~ FALSE,
                                             "true" ~ TRUE,
                                             "false" ~ FALSE,
                                             .default = NA)),
      valid_integration = as.logical(case_match(tolower(.data$valid_integration),
                                            "yes" ~ TRUE,
                                            "no"~ FALSE,
                                            "true" ~ TRUE,
                                            "false" ~ FALSE,
                                            .default = NA)),
      interference_feature_id = stringr::str_squish(.data$interference_feature_id),
      remarks = NA_character_) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish)) |>
    dplyr::select(
      "feature_id",
      "feature_class",
      "is_istd",
      "norm_istd_feature_id",
      "quant_istd_feature_id",
      feature_response_factor = "response_factor",
      "is_quantifier",
      "valid_integration",
      "interference_feature_id",
      "interference_proportion",
      "remarks"
    )

  # Handle the non-mandatory field Valid_Analysis
  if (all(is.na(d_annot$annot_features$valid_integration)))
    d_annot$annot_features$valid_integration <- TRUE
  else
    if (any(is.na(d_annot$annot_features$valid_integration)))
      cli::cli_abort("`Valid_Integration` is not defined for one or more features/analytes. Please check sheet 'Features (Analytes)'")

  # Handle the non-mandatory field Valid_Analysis
  if (all(is.na(d_annot$annot_features$is_quantifier)))
    d_annot$annot_features$is_quantifier <- TRUE
  else
    if (any(is.na(d_annot$annot_features$is_quantifier)))
      cli::cli_abort("`Quantifier` is not defined for one or more features/analytes. Please check sheet 'Features (Analytes)'")



  # ISTD annotation -------------------------

  annot_istd <- readxl::read_excel(path,
    sheet = "Internal Standards",
    trim_ws = TRUE,
    .name_repair = ~ if_else(nzchar(.x), .x, LETTERS[seq_along(.x)])
  )

  names(annot_istd) <- tolower(names(annot_istd))
  names(annot_istd)[1] <- "istd_feature_id"

  d_annot$annot_istd <- annot_istd |>
    # dplyr::mutate(istd_compound_name = na_character_) |>
    dplyr::select(
      quant_istd_feature_id = "istd_feature_id",
      istd_conc_nmolar = "istd_conc_[nm]"
    ) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  names(annot_istd) <- tolower(names(annot_istd))

  # RESPONSE CURVE annotation -------------------------

  d_annot_responsecurves <- readxl::read_excel(path, sheet = "Response Curves")
  names(d_annot_responsecurves) <- tolower(names(d_annot_responsecurves))

  d_annot$annot_responsecurves <- d_annot_responsecurves |>
    dplyr::select(
      analysis_id = "raw_data_filename",
      rqc_series_id = "response_curve_name",
      "relative_sample_amount" = "relative_sample_amount_[%]",
      "injection_volume" = "injection_volume_[ul]"
    ) |>
    dplyr::mutate(
      analysis_id = stringr::str_remove(.data$analysis_id, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)),
      analysis_id = stringr::str_squish(as.character(.data$analysis_id)),
      rqc_series_id = stringr::str_squish(as.character(.data$rqc_series_id)),
      relative_sample_amount = .data$relative_sample_amount / 100
    ) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  return(d_annot)
}
