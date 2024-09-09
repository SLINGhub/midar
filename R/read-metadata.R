
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
  if ("batch_id" %in% names(data@dataset_orig)) data <- loametadata_batches(data)



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

loametadata_batches <- function(data){
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

metadata_import_midarxlm<- function(data, path, analysis_sequence = "default", ignore_warnings, excl_unannotated_analyses = FALSE) {

  analysis_sequence <- rlang::arg_match(arg = analysis_sequence, c("timestamp", "resultfile", "metadata", "default"), multiple = FALSE)
  metadata <- read_metadata_midarxlm(path)

  metadata  <- assert_metadata(data = data, metadata = metadata, excl_unannotated_analyses = excl_unannotated_analyses, ignore_warnings = ignore_warnings)
  data  <- include_metadata(data = data, metadata = metadata, analysis_sequence = analysis_sequence,excl_unannotated_analyses = excl_unannotated_analyses)
  data
}


get_assert_summary_table <- function(list_of_errors, data=NULL, warn = TRUE, ...) {
  if (is.null(list_of_errors)) return(NULL)
  res <- as_tibble(do.call(rbind, list_of_errors))
  res <- res |>
    select(message, description, num.violations) |>
    mutate(message = str_replace_all(message, fixed('"'), "'")) |>
    mutate(Field = paste0(str_extract(unlist(message), "(?<=\\').+?(?=\\')")) ) |>
    separate(col = description, into = c("Type", "Issue", "Table"), sep = ";", remove = TRUE) |>
    filter(Type != "DX") |>
    select(Type, Table, Issue, Column = Field, Count = num.violations) |> ungroup()

  res$Count <- res$Count |> unlist()
  res
}


alert_assertions <- function(data, data_label, assert_type = c("defect", "warning"), ignore_warnings) {

  t1 <- get_assert_summary_table(attr(data$annot_analyses, "assertr_errors"))
  t2 <- get_assert_summary_table(attr(data$annot_features, "assertr_errors"))
  t3 <- get_assert_summary_table(attr(data$annot_istd, "assertr_errors"))
  t4 <- get_assert_summary_table(attr(data$annot_responsecurves, "assertr_errors"))

  t_all <- bind_rows(t1, t2, t3, t4) |>
    dplyr::mutate(Table = stringr::str_remove(Table, "metadata")) |>
    arrange(Type, Table, Count)



  if(any(t_all$Type == "D")){
    cli::cli_alert_warning(text = cli::col_red(glue::glue("Metadata is invalid with following defects:")))
    print(as_assertr_tibble(t_all))
    cli::cli_abort(message = cli::col_yellow(" Please check corresponding metadata tables and try again."), trace = NULL, call = caller_env())
  } else if(any(t_all$Type == "E")){
    cli::cli_alert_warning(text = cli::col_red(glue::glue("Metadata has following errors{ifelse(any(t_all$Type == 'W'), ' and warnings', '')}:")))
    print(as_assertr_tibble(t_all))
    cli::cli_abort(message = cli::col_yellow("Please check corresponding metadata and try again."), trace = NULL, call = caller_env())
  } else if (any(t_all$Type == "W" | t_all$Type == "N")){
    cli::cli_alert_warning(text = cli::col_yellow(glue::glue("Metadata has following warnings and notifications:")))
    print(as_assertr_tibble(t_all))
    if (!ignore_warnings)
      cli::cli_abort(message = cli::col_red("Please check corresponding metadata. Use `ignore_warnings`= TRUE to ignore these warnings."), trace = NULL, call = caller_env())
    else
      cli::cli_alert_warning(text = cli::col_yellow("Ignoring metadata warnings (as 'ignore_warnings' was set to TRUE)"))
  }
}


#' @title Add metadata an MidarExperiment object
#' @description Metadata provided as a list of tibbles will validates for consistency again loaded analysis data of the provided MidarExperiment object and then transfered.
#' @param data MidarExperiment object
#' @param path List of tibbles or data.frames containing analysis, feature, istd, response curve tables
#' @param analysis_sequence Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
#' @param excl_unannotated_analyses Exclude analyses (samples) that have no matching metadata
#' @return metadata list
#' @export
#'
#'
# Verify/assert metadata consistency with analysis data
 assert_metadata <- function(data, metadata, excl_unannotated_analyses, ignore_warnings) {

   #Note: to check for multiple missing columns defects, first each column will be check
   # for presence and error is raised for reporting, then a defect is raised to
   # disable further check, but without reporting it (using flag "DX")

  # ANALYSES METADATA ====================

  if (!is.null(metadata$annot_analyses) && nrow(metadata$annot_analyses) > 0){

    ## Check for data defects ----
    metadata$annot_analyses <- metadata$annot_analyses |>
      assertr::verify(has_any_name("analysis_id"), obligatory=FALSE, description = "D;Column missing;Analyses", defect_fun = defect_append) |>
      assertr::verify(has_any_name("qc_type"), obligatory=FALSE, description = "D;Column missing;Analyses", defect_fun = defect_append) |>
      assertr::verify(has_any_name("analysis_id", "qc_type"), obligatory=TRUE, description = "DX;Column missing;Analyses", defect_fun = assertr::defect_append)

    # Check data integrity ----
    metadata$annot_analyses <- metadata$annot_analyses|>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert(\(x){not_na(x)}, any_of(c("analysis_id", "qc_type")), obligatory=FALSE, description = "E;Missing value(s);Analyses") |>
      assertr::verify(all(assertr::is_uniq(analysis_id)), obligatory=FALSE, description = "E;Duplicated analysis IDs;Analyses") |>
      assertr::assert({\(x) all(is_na(x))}, all_of(c("qc_type")), description = "W;Missing value(s);Analyses") |>
      assertr::verify((analysis_id %in% unique(data@dataset_orig$analysis_id)), description = "W;Analyses not in analysis data;Analyses") |>
      assertr::verify((unique(data@dataset_orig$analysis_id) %in% analysis_id), description = "W;Analyses without metadata;Analyses") |>
      assertr::assert(\(x){not_na(x)}, any_of(c("valid_analysis", "sample_amount", "sample_amount_unit", "istd_volume", "batch_id", "replicate_no")), description = "W;Missing value(s);Analyses") |>
      assertr::chain_end(error_fun = error_append)
  }

  # FEATURES METADATA ====================
  if (!is.null(metadata$annot_features) && nrow(metadata$annot_features) > 0){

    ## Check for data defects ----
    ## TODO: interference_feature_id
    metadata$annot_features <- metadata$annot_features |>
      assertr::verify(assertr::has_all_names("feature_id"), obligatory=TRUE, description = "D;Column missing;Features", defect_fun = defect_append)

    ## Check data integrity ------
    metadata$annot_features <- metadata$annot_features |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert(\(x){not_na(x)}, any_of(c("feature_id")), obligatory=FALSE, description = "E;Missing value(s);Features") |>
      assertr::verify(has_any_name("feature_class","norm_istd_feature_id","quant_istd_feature_id","feature_response_factor","is_quantifier","valid_integration","interference_feature_id"), obligatory=FALSE, description = "E;No metadata field(s) provided;Features") |>
      assertr::assert(assertr::is_uniq, feature_id, obligatory=FALSE, description = "E;IDs duplicated;Features") |>
      assertr::verify(unique(norm_istd_feature_id) %in% feature_id, description = "E;ISTD(s) not defined as feature;Features") |>
      #assertr::verify(unique(quant_istd_feature_id) %in% feature_id, description = "E;ISTD(s) not defined as feature;Features") |>
      assertr::assert(\(x) {any(metadata$annot_istd$quant_istd_feature_id %in% (x))},quant_istd_feature_id, obligatory=FALSE, description = "W;ISTD(s) not defined;ISTDs") |>
      #assertr::verify(unique(interference_feature_id) %in% feature_id, description = "Interfering feature(s) not defined under 'feature_id';Features") |>
      assertr::verify((feature_id %in% unique(data@dataset_orig$feature_id)), description = "W;Feature(s) not in analysis data;Features") |>
      assertr::verify((unique(data@dataset_orig$feature_id) %in% feature_id), description = "W;Feature(s) without metadata;Features") |>
      #assertr::assert(assertr::within_bounds(lower.bound = 0, upper.bound = Inf, include.lower = FALSE, include.upper = FALSE), any_of(c("feature_response_factor", "interference_proportion")), description = "W;Values 0 or negative;Features") |>
      assertr::chain_end(error_fun = error_append)
   }

  # ISTD METADATA ====================

  if (!is.null(metadata$annot_istd) && nrow(metadata$annot_istd) > 0){
    ## Check for data defects ----
    metadata$annot_istd  <- metadata$annot_istd |>
      assertr::verify(assertr::has_all_names("quant_istd_feature_id"), obligatory=FALSE, description = "D;Column missing;ISTDs", defect_fun = defect_append) |>
      assertr::verify(assertr::has_all_names("istd_conc_nmolar"), obligatory=FALSE, description = "D;Column missing;ISTDs", defect_fun = defect_append) |>
      assertr::verify(assertr::has_all_names("quant_istd_feature_id", "istd_conc_nmolar"), obligatory=TRUE, description = "DX;Column missing;ISTDs", defect_fun = defect_append)

    ## Check data integrity =====
    metadata$annot_istd  <- metadata$annot_istd  |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert(\(x){not_na(x)}, all_of(c("quant_istd_feature_id", "istd_conc_nmolar")), obligatory=TRUE, description = "E;Missing value(s);ISTDs") |>
      #assertr::assert(\(x) {unique(x) %in% metadata$annot_istd$quant_istd_feature_id},quant_istd_feature_id, obligatory=TRUE, description = "W;Internal standard(s) not defined;ISTDs") |>
      assertr::verify(all(assertr::is_uniq(quant_istd_feature_id)), obligatory=TRUE, description = "E;Internal standard(s) duplicated;ISTDs") |>
      assertr::assert(\(x) {x %in% metadata$annot_features$feature_id}, quant_istd_feature_id, description = "W;Internal standard(s) not used;ISTDs") |>
      assertr::assert(\(x){x > 0}, any_of(c("istd_conc_nmolar")), description = "W;Values 0 or negative;ISTDs") |>
      assertr::chain_end(error_fun = warning_append)
  }

  # RESPONSE CURVE METADATA ====================

  if (!is.null(metadata$annot_responsecurves) && nrow(metadata$annot_responsecurves) > 0){
    ## Check for data defects ----
    metadata$annot_responsecurves <- metadata$annot_responsecurves |>
      assertr::verify(assertr::has_all_names("analysis_id", "rqc_series_id", "relative_sample_amount"), obligatory=TRUE, description = "D;Column missing;Response Curves", defect_fun = defect_append)

    ## Check data integrity ----
    ### TODO: check for qc_type RQC
    metadata$annot_responsecurves <- metadata$annot_responsecurves |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert(\(x){not_na(x)}, any_of(c("analysis_id", "rqc_series_id")), obligatory=TRUE, description = "E;Missing value(s);Response Curves") |>
      assertr::verify(all(assertr::is_uniq(analysis_id)), obligatory=TRUE, description = "E;Duplicated analysis IDs;Response Curves") |>
      assertr::verify((analysis_id %in% unique(data@dataset_orig$analysis_id)), description = "E;Analysis not present in analysis data;Response Curves") |>
      #assertr::verify((data@annot_analyses |> filter(qc_type == "RQC") |> pull(analysis_id) %in% analysis_id), description = "W;Analyses of QC type 'RQC' not defined;Response Curves") |>
      assertr::verify((analysis_id %in% unique(data@annot_analyses$analysis_id)), description = "W;Analysis not present in analysis metadata;Response Curves") |>
      assertr::assert(\(x){not_na(x)}, any_of(c("relative_sample_amount", "injection_volume")), description = "W;Missing value(s);Response Curves") |>
      assertr::chain_end(error_fun = warning_append)
  }

  alert_assertions(metadata, "Response Curves",assert_type = "defect", ignore_warnings)

  metadata
}



 #' @title Add metadata an MidarExperiment object
 #' @description Metadata provided as a list of tibbles will validates for consistency again loaded analysis data of the provided MidarExperiment object and then transfered.
 #' @param data MidarExperiment object
 #' @param path List of tibbles or data.frames containing analysis, feature, istd, response curve tables
 #' @param analysis_sequence Must by any of: "timestamp", "resultfile" or "metadata". Defines how the analysis order is determined. Default is "timestamp", when not available the sequence in the analysis results are used.
 #' @param excl_unannotated_analyses Exclude analyses (samples) that have no matching metadata
 #' @return metadata list
 #' @export
 #'
 #'

# Add verified metadata to the MidarExperiment object
include_metadata <- function(data, metadata, analysis_sequence = c("resultfile", "timestamp", "metadata", "default"), excl_unannotated_analyses = FALSE) {

  # ANALYSES METADATA ====================



  if (!is.null(metadata$annot_analyses) && nrow(metadata$annot_analyses) > 0){
    data@annot_analyses <- metadata$annot_analyses |>
      dplyr::semi_join(data@dataset_orig, by = "analysis_id") |>
      dplyr::bind_rows(pkg.env$dataset_templates$annot_analyses_template)

    # Define the analysis order

    if ("acquisition_time_stamp" %in% names(data@dataset_orig) && (analysis_sequence %in% c("timestamp", "default"))) {
      data@dataset_orig <- data@dataset_orig |> arrange(.data$acquisition_time_stamp)
      data@annot_analyses <- data@annot_analyses |> dplyr::arrange(match(.data$analysis_id, metadata$dataset_orig$analysis_id |> unique()))
    } else {
      if (analysis_sequence == "resultfile") {
        data@annot_analyses <- data@annot_analyses |> dplyr::arrange(match(.data$analysis_id, metadata$dataset_orig$analysis_id |> unique()))
      } else if (analysis_sequence == "metadata") {
        data@annot_analyses <- data@annot_analyses |> dplyr::arrange(match(.data$analysis_id, metadata$annot_analyses$analysis_id))
      } else if (analysis_sequence == "timestamp") {
        cli::cli_abort(call. = FALSE, "No acquisition timestamp field present in analysis results, please set parameter `analysis_sequence` to `resultfile` or `metadata` to define analysis order.")
      } else {
        cli::cli_alert_warning(cli::col_yellow(glue::glue("No acquisition timestamps present in results, order therefore based on analysis results sequence. Set parameter `analysis_sequence` to `metadata` to use this sequence as analysis order.")))
      }
    }


    # Set order/run id based on previous sorting
    # TODO: part of metadata or analysis data?
    data@annot_analyses <- data@annot_analyses |>
      mutate(run_id = row_number(), .before = 1)

    #TODO
    d_dataset <- data@dataset_orig %>%
      dplyr::inner_join(
        data@annot_analyses %>%
          dplyr::select("run_id", "analysis_id", "qc_type", "specimen", "sample_id", "replicate_no", "valid_analysis", "batch_id"),
        by = c("analysis_id")) %>%
      dplyr::inner_join(
        metadata$annot_features %>%
          filter(.data$valid_integration) |>
          dplyr::select(dplyr::any_of(c("feature_id", "feature_id", "feature_class", "norm_istd_feature_id", "quant_istd_feature_id", "is_istd", "feature_id", "is_quantifier", "valid_integration", "feature_response_factor", "interference_feature_id", "interference_proportion"))),
        by = c("feature_id"), keep = FALSE
      )

    cli_alert_success(col_green(glue::glue("Analysis metadata successfully associated with {length(data@dataset$analysis_id %>% unique())} samples.")))
  }

  # FEATURE METADATA ====================
  if (!is.null(metadata$annot_features) && nrow(metadata$annot_features) > 0){
    data@annot_features <- metadata$annot_features |>
      dplyr::semi_join(data@dataset_orig, by = "feature_id") |>
      dplyr::bind_rows(pkg.env$dataset_templates$annot_features_template)

    data@dataset_orig <- data@dataset_orig %>%
      dplyr::bind_rows(pkg.env$dataset_templates$dataset_orig_template)

    cli_alert_success(col_green(glue::glue("Feature metadata successfully associated with {length(data@dataset$feature_id %>% unique())} features.")))
  }




  # ISTD METADATA ====================

  if (!is.null(metadata$annot_istd) && nrow(metadata$annot_istd) > 0){
    data@annot_istd <- metadata$annot_istd  |>
      dplyr::semi_join(data@annot_features, by = join_by(quant_istd_feature_id == feature_id)) |>
      dplyr::bind_rows(pkg.env$dataset_templates$annot_istd_template)

    cli_alert_success(col_green(glue::glue("Internal Standard metadata successfully associated with {length(data@annot_istd$quant_istd_feature_id %>% unique())} ISTDs")))
  }


  # RQC METADATA ====================
  if (!is.null(metadata$annot_responsecurves) && nrow(metadata$annot_responsecurves) > 0){
    data@annot_responsecurves <- metadata$annot_responsecurves |>
      dplyr::semi_join(data@dataset_orig, by = "analysis_id") |>
      dplyr::bind_rows(pkg.env$dataset_templates$annot_responsecurves_template)

    cli_alert_success(col_green(glue::glue("Response Curve metadata successfully associated with {length(data@annot_responsecurves$analysis_id %>% unique())} analyses")))
  }

  # BATCHES METADATA ------------

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

  data@dataset <-
    dplyr::bind_rows(pkg.env$dataset_templates$dataset_template, d_dataset) |>
    mutate(
      corrected_interference = FALSE,
      outlier_technical = FALSE
    ) |>
    dplyr::arrange(match(.data$feature_id, metadata$annot_features$feature_id))

  data@dataset <- data@dataset |> arrange(.data$run_id)



  # stopifnot(methods::validObject(data, excl_nonannotated_analyses))
  check_integrity(data, excl_unannotated_analyses = excl_unannotated_analyses)
  data@status_processing <- "Annotated Raw Data"

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

  metadata <- list()

  # ANALYSIS/SAMPLE annotation --------
  d_temp_analyses<- readxl::read_excel(path, sheet = "Analyses (Samples)", trim_ws = TRUE)
  names(d_temp_analyses) <- tolower(names(d_temp_analyses))

  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "valid_analysis", init_value = TRUE, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "replicate_no", init_value = 1L, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "specimen", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "panel_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "sample_id", init_value = NA_character_, make_lowercase = FALSE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)

  metadata$annot_analyses <- d_temp_analyses |>
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
    if (all(is.na(metadata$annot_analyses$valid_analysis)))
      metadata$annot_analyses$valid_analysis <- TRUE
    else
      if (any(is.na(metadata$annot_analyses$valid_analysis)))
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

  metadata$annot_features <- d_temp_features |>
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
  if (all(is.na(metadata$annot_features$valid_integration)))
    metadata$annot_features$valid_integration <- TRUE
  else
    if (any(is.na(metadata$annot_features$valid_integration)))
      cli::cli_abort("`Valid_Integration` is not defined for one or more features/analytes. Please check sheet 'Features (Analytes)'")

  # Handle the non-mandatory field Valid_Analysis
  if (all(is.na(metadata$annot_features$is_quantifier)))
    metadata$annot_features$is_quantifier <- TRUE
  else
    if (any(is.na(metadata$annot_features$is_quantifier)))
      cli::cli_abort("`Quantifier` is not defined for one or more features/analytes. Please check sheet 'Features (Analytes)'")



  # ISTD annotation -------------------------

  annot_istd <- readxl::read_excel(path,
    sheet = "Internal Standards",
    trim_ws = TRUE,
    .name_repair = ~ if_else(nzchar(.x), .x, LETTERS[seq_along(.x)])
  )

  names(annot_istd) <- tolower(names(annot_istd))
  names(annot_istd)[1] <- "istd_feature_id"

  metadata$annot_istd <- annot_istd |>
    # dplyr::mutate(istd_compound_name = na_character_) |>
    dplyr::select(
      quant_istd_feature_id = "istd_feature_id",
      istd_conc_nmolar = "istd_conc_[nm]"
    ) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  names(annot_istd) <- tolower(names(annot_istd))

  # RESPONSE CURVE annotation -------------------------

  metadata_responsecurves <- readxl::read_excel(path, sheet = "Response Curves")
  names(metadata_responsecurves) <- tolower(names(metadata_responsecurves))

  metadata$annot_responsecurves <- metadata_responsecurves |>
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

  return(metadata)
}
