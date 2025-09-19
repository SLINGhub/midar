#' @title Retrieve Metadata from Imported Analysis Data
#' @description Retrieves available metadata from the imported analysis data and
#'   associates it with the provided MidarExperiment object.
#' @param data A `MidarExperiment` object
#' @param qc_type_column_name Column name in the imported raw data representing
#'   the `qc_type`
#' @return An updated `MidarExperiment` object
#'

import_metadata_from_data<- function(data = NULL, qc_type_column_name = "qc_type") {
  check_data(data)

  # get analysis metadata
   annot_analyses <- data@dataset_orig |>
    dplyr::select("analysis_id", dplyr::any_of(c("analysis_order", qc_type_column_name, "batch_id"))) |>
    dplyr::distinct() |>
    dplyr::rename(any_of(c(qc_type = qc_type_column_name)))

   if(!"analysis_order" %in% names(annot_analyses))
     annot_analyses <- annot_analyses |>  dplyr::mutate(analysis_order = row_number(), before = 1)

   annot_analyses <- clean_analysis_metadata(annot_analyses)

  # get feature metadata

  annot_features <- data@dataset_orig |>
    dplyr::select("feature_id", dplyr::any_of(c("feature_class",  "feature_sum_formula", "molecular_weight", "precursor_mz", "product_mz", "istd_feature_id", "is_quantifier"))) |>
    dplyr::distinct()
  annot_features <- clean_feature_metadata(annot_features)

  # check and add metadata
  metadata <- assert_metadata(data = data, metadata = list(annot_analyses = annot_analyses, annot_features = annot_features), ignore_warnings = FALSE, excl_unmatched_analyses = FALSE)
  data <- add_metadata(data = data, metadata = metadata, excl_unmatched_analyses = FALSE)
  data <- link_data_metadata(data)

  data
}

# Retrieve batch info from analysis metadata

get_metadata_batches <- function(annot_analyses){

  annot_batches <- annot_analyses |>
    mutate(batch_no = dplyr::cur_group_id(), .by = c("batch_id"), .before = "batch_id") |>
    dplyr::group_by(.data$batch_no) |>
    dplyr::summarise(
      batch_id = .data$batch_id[1],
      id_batch_start = dplyr::first(.data$analysis_order),
      id_batch_end = dplyr::last(.data$analysis_order)
    ) |>
    dplyr::ungroup() |>
    dplyr::arrange(.data$id_batch_start)

  annot_batches <-  dplyr::bind_rows(pkg.env$table_templates$annot_batch_info_template,annot_batches)
  annot_batches
}


#' @title Import Metadata from a MIDAR Metadata Organizer file
#' @description Imports metadata from a 'MIDAR Metadata Organizer' file (.xlsm/.xlsx) file and associates it with analysis data.
#' @param data A `MidarExperiment` object
#' @param path File name and path of the 'MIDAR Metadata Organizer' file (.xlsm/.xlsx) file
#' @param excl_unmatched_analyses Exclude analyses (samples) that have no matching metadata
#' @param ignore_warnings Ignore warnings from data validation and proceed with importing metadata
#' @return An updated `MidarExperiment` object

#' @examples
#' mexp <- MidarExperiment()
#'
#' mexp <- import_data_mrmkit(
#'   data = mexp,
#'   path = system.file("extdata", "MRMkit_demo.tsv", package = "midar"),
#'   import_metadata = TRUE)
#'
#' mexp <- import_metadata_msorganiser(
#'  data = mexp,
#'  path = system.file("extdata", "Example_Metadata_1.xlsm", package = "midar"),
#'  excl_unmatched_analyses = FALSE,
#'  ignore_warnings = TRUE)
#'
#' print(mexp)
#'
#' @export
import_metadata_msorganiser <- function(data = NULL, path,  ignore_warnings = FALSE, excl_unmatched_analyses = FALSE) {
  check_data(data)

  metadata <- read_metadata_msorganiser(path)
  metadata  <- assert_metadata(data, metadata = metadata, ignore_warnings = ignore_warnings, excl_unmatched_analyses = excl_unmatched_analyses)
  data  <- add_metadata(data, metadata = metadata, excl_unmatched_analyses = excl_unmatched_analyses)
  data <- link_data_metadata(data)
  data
}

#' @title Import analysis metadata
#' @description Imports analysis metadata (annotation) from a preloaded data frame or tibble via the `data` argument, or from data from a file (CSV or Excel) via the `path` argument.
#' The analysis metadata must contain following columns: `analysis_id` and `qc_type`. Additional analysis metadata columns are described under details below.
#' @param data A `MidarExperiment` object
#' @param table A data frame or tibble with analysis (sample) metadata. If `path` is also provided, an error will be raised.
#' @param path A character string specifying the path to a CSV (.csv) or Excel (.xlsx) file. If `table` is also provided, an error will be raised.
#' @param sheet Defines the sheet name in case an Excel file is provided.
#' @param excl_unmatched_analyses Exclude analyses (samples) that have no matching metadata
#' @param ignore_warnings Ignore warnings from data validation and proceed with importing metadata
#' @return An updated `MidarExperiment` object
#' @examples
#' mexp <- MidarExperiment()

#' file_path = system.file("extdata", "MHQuant_demo.csv", package = "midar")
#' mexp <- import_data_masshunter(
#'   data = mexp,
#'   path = file_path,
#'   import_metadata = FALSE)
#'
#' meta_path = system.file("extdata", "MHQuant_demo_metadata_analyses.csv", package = "midar")
#'
#' mexp <- import_metadata_analyses(
#'   data = mexp,
#'   path = meta_path,
#'   excl_unmatched_analyses = TRUE)
#'
#' print(mexp)
#' @export
#'
import_metadata_analyses <- function(data = NULL, table = NULL, path = NULL,  sheet = NULL, ignore_warnings = FALSE, excl_unmatched_analyses = FALSE) {
  check_data(data)
  tbl_metadata <- get_metadata_table(table, path, sheet)
  tbl_metadata <- clean_analysis_metadata(tbl_metadata )
  metadata  <- assert_metadata(data, metadata = list(annot_analyses = tbl_metadata), ignore_warnings, excl_unmatched_analyses)
  data  <- add_metadata(data, metadata = metadata, excl_unmatched_analyses = excl_unmatched_analyses)
  data <- link_data_metadata(data)
  data
}

#' @title Import feature metadata
#' @description Imports analysis metadata (annotation) from a preloaded data frame or tibble via the `data` argument,  or from data from a file (CSV or Excel) via the `path` argument.
#' The analysis metadata must contain following columns: `analysis_id` and `qc_type`. Additional analysis metadata columns are described under details below.
#' @param data A `MidarExperiment` object
#' @param table A data frame or tibble with analysis (sample) metadata. If `path` is also provided, an error will be raised.
#' @param path A character string specifying the path to a CSV (.csv) or Excel (.xlsx) file. If `table` is also provided, an error will be raised.
#' @param sheet Defines the sheet name in case an Excel file is provided.
# #' @param excl_unmatched_analyses Exclude analyses (samples) that have no matching metadata
#' @param ignore_warnings Ignore warnings from data validation and proceed with importing metadata
#' @return An updated `MidarExperiment` object
#' @export
#'
import_metadata_features <- function(data = NULL, table = NULL, path = NULL,  sheet = NULL, ignore_warnings = FALSE) {
  check_data(data)
  tbl_metadata <- get_metadata_table(table, path, sheet)
  tbl_metadata <- clean_feature_metadata(tbl_metadata)
  metadata  <- assert_metadata(data, metadata = list(annot_features = tbl_metadata), ignore_warnings, excl_unmatched_analyses = FALSE)
  data  <- add_metadata(data, metadata = metadata, excl_unmatched_analyses = FALSE)
  data <- link_data_metadata(data)
  data
}


#' @title Import Internal Standards (ISTD) metadata
#' @description Imports ISTD metadata (annotation) from a preloaded data frame or tibble via the `data` argument, or from data from a file (CSV or Excel) via the `path` argument.
#' The analysis metadata must contain following columns: `istd_feature_id` and one of `istd_conc_nmolar` or `istd_conc_ngml`.
#' @param data A `MidarExperiment` object
#' @param table A data frame or tibble with analysis (sample) metadata. If `path` is also provided, an error will be raised.
#' @param path A character string specifying the path to a CSV (.csv) or Excel (.xlsx) file. If `table` is also provided, an error will be raised.
#' @param sheet Defines the sheet name in case an Excel file is provided.
# #' @param excl_unmatched_analyses Exclude analyses (samples) that have no matching metadata
#' @param ignore_warnings Ignore warnings from data validation and proceed with importing metadata
#' @return An updated `MidarExperiment` object
#' @export
#'
import_metadata_istds <- function(data = NULL, table = NULL, path = NULL,  sheet = NULL, ignore_warnings = FALSE) {
  check_data(data)
  tbl_metadata <- get_metadata_table(table, path, sheet)
  tbl_metadata <- clean_istd_metadata(tbl_metadata)
  metadata  <- assert_metadata(data, metadata = list(annot_istds = tbl_metadata), ignore_warnings, excl_unmatched_analyses = FALSE)
  data  <- add_metadata(data, metadata = metadata)
  data <- link_data_metadata(data)
  data
}


#' @title Import response curves metadata
#' @description Imports response curve metadata (annotation) from a preloaded data frame or tibble via the `data` argument, or from data from a file (CSV or Excel) via the `path` argument.
#' The analysis metadata must contain following columns: `analysis_id`, `curve_id`, `analyzed_amount` and `analyzed_amount_unit`.
#' @param data A `MidarExperiment` object
#' @param table A data frame or tibble with response curve metadata. If `path` is also provided, an error will be raised.
#' @param path A character string specifying the path to a CSV (.csv) or Excel (.xlsx) file. If `table` is also provided, an error will be raised.
#' @param sheet Defines the sheet name in case an Excel file is provided.
#' @param ignore_warnings Ignore warnings from data validation and proceed with importing metadata
#' @return An updated `MidarExperiment` object
#' @export
#'
import_metadata_responsecurves <- function(data = NULL, table = NULL, path = NULL,  sheet = NULL, ignore_warnings = FALSE) {
  check_data(data)
  tbl_metadata <- get_metadata_table(table, path, sheet)
  tbl_metadata <- clean_response_metadata(tbl_metadata )
  metadata  <- assert_metadata(data, metadata = list(annot_responsecurves = tbl_metadata), ignore_warnings, excl_unmatched_analyses = FALSE)
  data  <- add_metadata(data, metadata = metadata)
  data <- link_data_metadata(data)
  data
}

#' @title Import calibration curves metadata
#' @description Imports calibration curve metadata (annotation) from a preloaded data frame or tibble via the `data` argument, or from data from a file (CSV or Excel) via the `path` argument.
#' The analysis metadata must contain following columns: `analysis_id`, `curve_id`, `feature_id`, `concentration`, and `concentration_unit`.
#' @param data A `MidarExperiment` object
#' @param table A data frame or tibble with calibration curve metadata. If `path` is also provided, an error will be raised.
#' @param path A character string specifying the path to a CSV (.csv) or Excel (.xlsx) file. If `table` is also provided, an error will be raised.
#' @param sheet Defines the sheet name in case an Excel file is provided.
#' @param ignore_warnings Ignore warnings from data validation and proceed with importing metadata
#' @return An updated `MidarExperiment` object
#' @export
#'
import_metadata_qcconcentrations <- function(data = NULL, table = NULL, path = NULL,  sheet = NULL, ignore_warnings = FALSE) {
  check_data(data)
  tbl_metadata <- get_metadata_table(table, path, sheet)
  tbl_metadata <- clean_qcconc_metadata(tbl_metadata)
  metadata  <- assert_metadata(data, metadata = list(annot_qcconcentrations = tbl_metadata), ignore_warnings, excl_unmatched_analyses = FALSE)
  data  <- add_metadata(data, metadata = metadata)
  data <- link_data_metadata(data)
  data
}





get_assert_summary_table <- function(list_of_errors, data=NULL, warn = TRUE, ...) {

  if (is.null(list_of_errors)) return(NULL)
  res <- as_tibble(do.call(rbind, list_of_errors))

  res <- res |>
    rowwise() |>
    select("message", "description", "num.violations") |>
    mutate(message = str_replace_all(unlist(.data$message)[1], fixed('"'), "'")) |>
    mutate(Field = paste0(str_extract(unlist(.data$message), "(?<=\\').+?(?=\\')")) ) |>
    tidyr::separate(col = "description", into = c("Type", "Issue", "Table", "TargetField"), sep = ";", remove = TRUE) |>
    filter(.data$Type != "DX") |>
    mutate(Field = if_else(.data$Field == "NA", .data$TargetField, .data$Field)) |>
    #mutate(num.violations = if_else(.data$verb == "verify", "", .data$num.violations)) |>
    select("Type", "Table", Column = "Field", "Issue", Count = "num.violations") |> ungroup()

  res$Count <- res$Count |> unlist()
  res
}


print_assertion_summary <- function(data, metadata_new, data_label, assert_type = c("defect", "warning"), ignore_warnings, excl_unmatched_analyses) {
  t1 <- get_assert_summary_table(attr(data$annot_analyses, "assertr_errors"))
  t2 <- get_assert_summary_table(attr(data$annot_features, "assertr_errors"))
  t3 <- get_assert_summary_table(attr(data$annot_istds, "assertr_errors"))
  t4 <- get_assert_summary_table(attr(data$annot_responsecurves, "assertr_errors"))
  t5 <- get_assert_summary_table(attr(data$annot_qcconcentrations, "assertr_errors"))

  if(!is.null(t1)) t1 <- t1 |> mutate(ignore_warn_flag = attr(data$annot_analyses, "ignore_warnings"))
  if(!is.null(t2)) t2 <- t2 |> mutate(ignore_warn_flag = attr(data$annot_features, "ignore_warnings"))
  if(!is.null(t3)) t3 <- t3 |> mutate(ignore_warn_flag = attr(data$annot_istds, "ignore_warnings"))
  if(!is.null(t4)) t4 <- t4 |> mutate(ignore_warn_flag = attr(data$annot_responsecurves, "ignore_warnings"))
  if(!is.null(t5)) t5 <- t5 |> mutate(ignore_warn_flag = attr(data$annot_qcconcentrations, "ignore_warnings"))

  t_all <- bind_rows(t1, t2, t3, t4, t5)

  if(is.null(t_all) | nrow(t_all) == 0)
    return(NULL)
  else
    t_all <- t_all |> arrange("Type", "Table", "Count")

  # Ensure ignore_warnings only applies to the current metadata import
  if(!is.null(t1)) {
    t1 <- t1 |> mutate(ignore_warn_flag = attr(data$annot_analyses, "ignore_warnings"))
    if(attr(data$annot_analyses, "excl_unmatched_analyses") & "Analyses without metadata" %in% t_all$Issue)
      t_all$Type[t_all$Issue == "Analyses without metadata"] <- "N"

    if ("Analyses without metadata" %in% (t_all |> filter(.data$Type == "W") |> pull(.data$Issue))  & !ignore_warnings)
      if ("annot_analyses" %in% names(metadata_new))
        cli::cli_abort(message = cli::col_red("Not all analyses listed in metadata are present in the data. Please verify data or set `excl_unmatched_analyses = TRUE`"), trace = NULL, call = caller_env())
  }

  #else
  #cli::cli_alert_warning(text = cli::col_yellow("Ignoring metadata warnings (as 'ignore_warnings' was set to TRUE)"))
  t_all_print <- t_all |>
    mutate(Type = if_else(.data$Type == "W" & .data$ignore_warn_flag, "W*", .data$Type)) |>
    select(-"ignore_warn_flag")
  if(any(t_all$Type == "D")){
    cli::cli_alert_warning(text = cli::col_red(glue::glue("Metadata is invalid with following defects:")))
    print(as_assertr_tibble(t_all_print))
    cli::cli_abort(message = cli::col_yellow(" Please verify corresponding metadata tables and try again."), trace = NULL, call = caller_env())
  } else if(any(t_all$Type == "E")){
    cli::cli_alert_warning(text = cli::col_red(glue::glue("Metadata has following errors{ifelse(any(t_all$Type == 'W'), ' and warnings', '')}:")))
    print(as_assertr_tibble(t_all_print))
    cli::cli_abort(message = cli::col_yellow("Please verify corresponding metadata and try again."), trace = NULL, call = caller_env())
  } else if (any(t_all$Type == "W" | t_all$Type == "N")){
    cli::cli_alert_warning(text = cli::col_yellow(glue::glue("Metadata has following warnings and notifications:")))
    print(as_assertr_tibble(t_all_print))
  }

  if (!ignore_warnings & any(t_all |> filter(!.data$ignore_warn_flag) |>  pull(.data$Type) == "W"))
    cli::cli_abort(message = cli::col_red("Please verify warnings in corresponding metadata. Use `ignore_warnings`= TRUE to ignore warnings."), trace = NULL, call = caller_env())

}


#' @title Add metadata an MidarExperiment object
#' @description Metadata provided as a list of tibbles will validates for consistency again loaded analysis data of the provided MidarExperiment object and then transfered.
#' @param data MidarExperiment object
#' @param metadata List of tibbles or data.frames containing analysis, feature, istd, response curve tables
#' @param excl_unmatched_analyses Exclude analyses (samples) that have no matching metadata
#' @param ignore_warnings Ignore data validation warnings and proceed with adding metadata
#' @return metadata list
#'
#'
# Verify/assert metadata consistency with analysis data

#TODO: align with metadata assertions and when it is check_integrity called
 assert_metadata <- function(data = NULL, metadata, ignore_warnings, excl_unmatched_analyses ) {
   check_data(data)
   #TODO: NOTE to check for multiple missing columns defects, first each column will be check
   # for presence and error is raised for reporting, then a defect is raised to
   # disable further check, but without reporting it (using flag "DX")
   # Columns with all values =  NA will be ignored

  metadata_new <- metadata
   
# retrieve vom data if available and if not provided as metadata (for single annot adding)
 if(is.null(metadata$annot_analyses)) metadata$annot_analyses <- data@annot_analyses
 if(is.null(metadata$annot_features)) metadata$annot_features <- data@annot_features
 if(is.null(metadata$annot_istds)) metadata$annot_istds <- data@annot_istds
 if(is.null(metadata$annot_responsecurves)) metadata$annot_responsecurves <- data@annot_responsecurves
 if(is.null(metadata$annot_qcconcentrations)) metadata$annot_qcconcentrations <- data@annot_qcconcentrations

  # ANALYSES METADATA ====================
  if (!is.null(metadata$annot_analyses) && nrow(metadata$annot_analyses) > 0){

    ## Check for data defects ----
    #TODO remove: metadata$annot_analyses$qc_type[3] <- NA

    metadata$annot_analyses <- metadata$annot_analyses |>
      assertr::verify(has_any_name("analysis_id"), obligatory=FALSE, description = "D;Column missing;Analyses;analysis_id", defect_fun = assertr::defect_append) |>
      assertr::verify(has_any_name("qc_type"), obligatory=FALSE, description = "D;Column missing;Analyses;qc_type", defect_fun = assertr::defect_append) |>
      assertr::verify(has_any_name("analysis_id", "qc_type"), obligatory=TRUE, description = "DX;Column missing;Analyses;analysis_id|qc_type", defect_fun = assertr::defect_append)

    # Check data integrity ----
    metadata$annot_analyses <- metadata$annot_analyses |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert(assertr::is_uniq, "analysis_id", obligatory=FALSE, description = "E;IDs duplicated;Analyses;analysis_id") |>
      assertr::assert(\(x){not_na(x)}, any_of(c("analysis_id", "qc_type")), obligatory=FALSE, description = "E;Missing value(s);Analyses;analysis_id|qc_type") |>
      assertr::assert(\(x){not_na(x)}, where(\(x){!all(is.na(x))}) & dplyr::any_of(c("sample_id")), description = "N;Not defined for all analyses;Analyses;sample_id") |>
      assertr::assert(\(x){not_na(x)}, where(\(x){!all(is.na(x))}) & dplyr::any_of(c("sample_amount", "sample_amount_unit")), description = "W;Incomplete value(s);Analyses;sample_amount|sample_amount_unit") |>
      assertr::assert(\(x){not_na(x)}, where(\(x){!all(is.na(x))}) & dplyr::any_of(c("istd_volume")), description = "W;Not defined for all analyses;Analyses;istd_volume") |>
      assertr::assert(\(x){not_na(x)}, where(\(x){!all(is.na(x))}) & dplyr::any_of(c("batch_id")), description = "W;Not defined for all analyses;Analyses;batch_id") |>
      assertr::assert(\(x){not_na(x)}, where(\(x){!all(is.na(x))}) & dplyr::any_of(c("replicate_no")), description = "W;Not defined for all analyses;Analyses;replicate_no") |>
      assertr::assert(\(x){not_na(x)}, where(\(x){!all(is.na(x))}) & dplyr::any_of(c("valid_analysis")), description = "E;Not defined for all analyses;Analyses;valid_analysis")
    #assertr::verify(all(assertr::is_uniq(analysis_id)), obligatory=FALSE, description = "E;Duplicated analysis IDs;Analyses;analysis_id")

    if(!is.null(data)){
      metadata$annot_analyses <- metadata$annot_analyses |>
        assertr::verify((.data$analysis_id %in% unique(data@dataset_orig$analysis_id)), description = "W;Analyses not in analysis data;Analyses;analysis_id") |>
        assertr::verify((unique(data@dataset_orig$analysis_id) %in% .data$analysis_id), description = "W;Analyses without metadata;Analyses;analysis_id")
    }
    metadata$annot_analyses <- metadata$annot_analyses |> assertr::chain_end(error_fun = assertr::error_append)

    if("annot_analyses" %in% names(metadata_new)) {
      attr(metadata$annot_analyses, "excl_unmatched_analyses") <- excl_unmatched_analyses
      attr(metadata$annot_analyses, "ignore_warnings") <- ignore_warnings
    }
  }

  # FEATURES METADATA ====================
  if (!is.null(metadata$annot_features) && nrow(metadata$annot_features) > 0){

    ## Check for data defects ----
    ## TODO: interference_feature_id
    metadata$annot_features <- metadata$annot_features |>
      assertr::verify(assertr::has_all_names("feature_id"), obligatory=TRUE, description = "D;Column missing;Features", defect_fun = assertr::defect_append)

    ## Check data integrity ------
    #TODO: check if ISTD is also present in ISTD conc table. THis will need a seperate check to as the same ISTD occurs multiple times.

    metadata$annot_features <- metadata$annot_features |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert(\(x){not_na(x)}, any_of(c("feature_id")), obligatory=FALSE, description = "E;Missing value(s);Features") |>
      assertr::verify(has_any_name("feature_class","feature_sum", "istd_feature_id","quant_istd_feature_id","response_factor","is_quantifier","valid_feature","interference_feature_id"), obligatory=FALSE, description = "E;No metadata field(s) provided;Features; ") |>
      assertr::assert(assertr::is_uniq, "feature_id", obligatory=FALSE, description = "E;IDs duplicated;Features;feature_id") |>
      assertr::assert(assertr::in_set(unique(metadata$annot_features$feature_id)), "istd_feature_id", description = "E;ISTD(s) not defined as feature;Features;feature_id") |>
      assertr::assert(assertr::in_set(NA, "linear", "quadratic"), "curve_fit_model", description = "E; Must be NA, 'linear' or 'quadratic';Features;curve_fit_model") |>
      assertr::assert(assertr::in_set(NA, "1/x", "1/x^2"), "curve_fit_weighting", description = "E; Must be NA, '1/x' or '1/x^2';Features;curve_fit_weighting") |>
      assertr::assert(assertr::in_set(unique(metadata$annot_features$feature_id)), "interference_feature_id", description = "E;Interfering feature(s) not defined as feature;Features;interference_feature_id") |>
      assertr::verify(any(!xor(is.na(metadata$annot_features$interference_contribution), is.na(metadata$annot_features$interference_feature_id))), "interference_feature_id", obligatory=FALSE, description = "E;Incomplete interference info;Features;interference_contribution") |>
      assertr::assert(\(x){not_na(x)}, where(\(x){!all(is.na(x))}) & dplyr::any_of(c("analyte_id")), description = "N;Not defined for all features;Features;analyte_id")
      #assertr::assert(\(x){any(xor(is.na(x), is.na(metadata$annot_features$interference_feature_id)))}, "interference_feature_id", obligatory=FALSE, description = "E;Missing interference proportion(s);Features;interference_contribution") |>


    if(!is.null(data)){
      metadata$annot_features <- metadata$annot_features |>
        assertr::verify((.data$feature_id %in% unique(data@dataset_orig$feature_id)), description = "W;Feature(s) not in analysis data;Features;feature_id") |>
        assertr::verify((unique(data@dataset_orig$feature_id) %in% .data$feature_id), description = "W;Feature(s) without metadata;Features;feature_id")
    }
      #assertr::assert(assertr::within_bounds(lower.bound = 0, upper.bound = Inf, include.lower = FALSE, include.upper = FALSE), any_of(c("response_factor", "interference_contribution")), description = "W;Values 0 or negative;Features") |>
      metadata$annot_features <- metadata$annot_features |> assertr::chain_end(error_fun = assertr::error_append)

      if("annot_features" %in% names(metadata_new)) {
        attr(metadata$annot_features, "ignore_warnings") <- ignore_warnings
      }

   }

  # ISTD METADATA ====================

  if (!is.null(metadata$annot_istds) && nrow(metadata$annot_istds) > 0){
    ## Check for data defects ----
    metadata$annot_istds  <- metadata$annot_istds |>
      assertr::verify(assertr::has_all_names("quant_istd_feature_id"), obligatory=FALSE, description = "D;Column missing;ISTDs;quant_istd_feature_id", defect_fun = assertr::defect_append) |>
      # TODO ngml / assertr::verify(assertr::has_all_names("istd_conc_nmolar"), obligatory=FALSE, description = "D;Column missing;ISTDs;istd_conc_nmolar", defect_fun = assertr::defect_append) |>
      assertr::verify(assertr::has_all_names("quant_istd_feature_id", "istd_conc_nmolar"), obligatory=TRUE, description = "DX;Column missing;ISTDs;quant_istd_feature_id", defect_fun = assertr::defect_append)

    ## Check data integrity =====
    ### TODO: istd_conc_nmolar check, make optinal istd_conc_nmolar and istd_conc_ngml
    metadata$annot_istds  <- metadata$annot_istds  |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert(\(x){not_na(x)}, all_of(c("quant_istd_feature_id")), obligatory=TRUE, description = "E;Missing value(s);ISTDs; ") |>
      #assertr::assert(\(x) {unique(x) %in% metadata$annot_istds$quant_istd_feature_id},quant_istd_feature_id, obligatory=TRUE, description = "W;Internal standard(s) not defined;ISTDs;quant_istd_feature_id") |>
      assertr::verify(all(assertr::is_uniq("quant_istd_feature_id")), obligatory=TRUE, description = "E;Internal standard(s) duplicated;ISTDs;quant_istd_feature_id") |>
      assertr::assert(\(x) {x %in% metadata$annot_features$feature_id}, "quant_istd_feature_id", description = "W;Internal standard(s) not used;ISTDs;feature_id") |>
      #assertr::assert(\(x){x > 0}, any_of(c("istd_conc_nmolar")), description = "W;Values 0 or negative;ISTDs;istd_conc_nmolar") |>
      assertr::chain_end(error_fun = assertr::error_append)

    if("annot_istds" %in% names(metadata_new)) {
      attr(metadata$annot_istds, "ignore_warnings") <- ignore_warnings
    }

  }

  # RESPONSE CURVE METADATA ====================

  if (!is.null(metadata$annot_responsecurves) && nrow(metadata$annot_responsecurves) > 0){
    ## Check for data defects ----
    metadata$annot_responsecurves <- metadata$annot_responsecurves |>
      assertr::verify(assertr::has_all_names("analysis_id", "curve_id", "analyzed_amount"), obligatory=TRUE, description = "D;Column missing;Response Curves; ", defect_fun = assertr::defect_append)

    ## Check data integrity ----
    ### TODO: check for qc_type RQC
    metadata$annot_responsecurves <- metadata$annot_responsecurves |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert(\(x){not_na(x)}, any_of(c("analysis_id", "curve_id")), obligatory=TRUE, description = "E;Missing value(s);Response Curves;analysis_id|curve_id") |>
      assertr::verify(all(assertr::is_uniq(.data$analysis_id)), obligatory=TRUE, description = "E;Duplicated analysis IDs;Response Curves;analysis_id") |>
      assertr::verify(check_groupwise_identical_ids(metadata$annot_responsecurves , group_col = "curve_id", id_col = "analyzed_amount_unit"), obligatory=FALSE, description = "W;Units not identical in at least one group;Response Curves;analyzed_amount_unit") |>
      #assertr::verify((data@annot_analyses |> filter(qc_type == "RQC") |> pull(analysis_id) %in% analysis_id), description = "W;Analyses of QC type 'RQC' not defined;Response Curves;analysis_id") |>
      assertr::assert(\(x){not_na(x)}, any_of(c("analyzed_amount")), description = "W;Missing value(s);Response Curves;analyzed_amount")
    if(!is.null(data)){
      metadata$annot_responsecurves <- metadata$annot_responsecurves |>
        assertr::verify((.data$analysis_id %in% unique(data@dataset_orig$analysis_id)), description = "E;Analysis not present in analysis data;Response Curves;analysis_id")
    }
      metadata$annot_responsecurves <- metadata$annot_responsecurves |> assertr::chain_end(error_fun = assertr::error_append)

    if("annot_responsecurves" %in% names(metadata_new)) {
      attr(metadata$annot_responsecurves, "ignore_warnings") <- ignore_warnings
    }

  }


  # QC CONCENTRATION METADATA ====================

  if (!is.null(metadata$annot_qcconcentrations) && nrow(metadata$annot_qcconcentrations) > 0){
    ## Check for data defects ----
    metadata$annot_qcconcentrations <- metadata$annot_qcconcentrations |>
      assertr::verify(assertr::has_all_names("sample_id", "analyte_id", "concentration", "concentration_unit"), obligatory=TRUE, description = "D;Column missing;QC concentrations; ", defect_fun = assertr::defect_append)

    ## Check data integrity ----
    ### TODO: check for qc_type RQC
    metadata$annot_qcconcentrations <- metadata$annot_qcconcentrations |>
      assertr::chain_start(store_success = FALSE) |>
      assertr::assert(\(x){not_na(x)}, any_of(c("sample_id")), obligatory=TRUE, description = "E;Missing value(s);QC concentrations;sample_id") |>
      assertr::verify(check_groupwise_identical_ids(metadata$annot_qcconcentrations , group_col = "sample_id", id_col = .data$concentration_unit), obligatory=FALSE, description = "W;Units not identical in at least one group;QC concentrations;analyzed_amount_unit") |>
      assertr::assert(\(x){not_na(x)}, any_of(c("concentration")), description = "W;Missing value(s);QC concentrations;analyzed_amount") |>
      assertr::assert(\(x){not_na(x)}, where(\(x){!all(is.na(x))}) & dplyr::any_of(c("include_in_analysis")), description = "E;Not defined for all entries;QC concentrations;include_in_analysis")

    if(!is.null(data)){
      metadata$annot_qcconcentrations <- metadata$annot_qcconcentrations |>
        assertr::verify((.data$sample_id %in% unique(metadata$annot_analyses$sample_id)), description = "N;Samples not defined in analysis data;QC concentrations;sample_id") |>
        assertr::verify((.data$analyte_id %in% unique(metadata$annot_features$analyte_id)), description = "N;Analytes not defined in analysis data;QC concentrations;analyte_id")
    }
    metadata$annot_qcconcentrations <- metadata$annot_qcconcentrations |>
      assertr::verify(all(metadata$annot_features$analyte_id[!(metadata$annot_features$feature_id %in% metadata$annot_istds$quant_istd_feature_id)] %in% unique(.data$analyte_id)), description = "N;Non-ISTD analytes missing from QC concentrations;QC concentrations;analyte_id")
    metadata$annot_qcconcentrations <- metadata$annot_qcconcentrations |> assertr::chain_end(error_fun = assertr::error_append)

    if("annot_qcconcentrations" %in% names(metadata_new)) {
      attr(metadata$annot_qcconcentrations, "ignore_warnings") <- ignore_warnings
    }

  }


  # retrieve vom data if available and if not provide as metadata (for single annot adding)
  if(is.null(metadata_new$annot_analyses)) metadata$annot_analyses <- data@annot_analyses
  if(is.null(metadata_new$annot_features)) metadata$annot_features <- data@annot_features
  if(is.null(metadata_new$annot_istds)) metadata$annot_istds <- data@annot_istds
  if(is.null(metadata_new$annot_responsecurves)) metadata$annot_responsecurves <- data@annot_responsecurves
  if(is.null(metadata_new$annot_qcconcentrations)) metadata$annot_qcconcentrations <- data@annot_qcconcentrations

  print_assertion_summary(data = metadata, metadata_new = metadata_new, data_label = "Response Curves",assert_type = "defect", ignore_warnings = ignore_warnings, excl_unmatched_analyses = excl_unmatched_analyses)


  metadata
 }


 #' @title Add metadata an MidarExperiment object
 #' @description Metadata provided as a list of tibbles will validates for consistency again loaded analysis data of the provided MidarExperiment object and then transfered.
 #' @param data MidarExperiment object
 #' @param metadata List of tibbles or data.frames containing analysis, feature, istd, response curve tables
 #' @param excl_unmatched_analyses Exclude analyses (samples) that have no matching metadata
 #' @return metadata list
 #' @export


# Add verified metadata to the MidarExperiment object
add_metadata <- function(data = NULL, metadata, excl_unmatched_analyses = FALSE) {
  check_data(data)
  # ANALYSES METADATA ====================
  if (!is.null(metadata$annot_analyses) & nrow(metadata$annot_analyses) > 0){

    data@annot_analyses <- metadata$annot_analyses

    # Get info for cli output
    n_match <- intersect(data@dataset_orig$analysis_id, data@annot_analyses |> filter(.data$valid_analysis) |> pull(.data$analysis_id)) |> length()
    cli_alert_success(col_green(glue::glue("Analysis metadata associated with {n_match} analyses.")))

    n_invalid <- data@annot_analyses |> filter(!.data$valid_analysis) |> nrow()
    if (n_invalid > 0){
      cli_alert_info(col_yellow("{n_invalid} invalid analyses (as defined in the metadata) were excluded"))
    }

  }
  # FEATURE METADATA ====================
  if (!is.null(metadata$annot_features) & nrow(metadata$annot_features) > 0){
    # Determines if a feature is used as an ISTD
    all_istds <- unique(metadata$annot_features$istd_feature_id) |> na.omit()
    metadata$annot_features <- metadata$annot_features |>
      mutate(is_istd = .data$feature_id %in% all_istds, .after = "feature_class")

    data@annot_features <- metadata$annot_features

    # Get info for cli output
    n_match <- intersect(data@dataset_orig$feature_id, data@annot_features |> filter(.data$valid_feature) |> pull(.data$feature_id)) |> length()
    cli_alert_success(col_green(glue::glue("Feature metadata associated with {n_match} features.")))

    n_invalid <- data@annot_features |> filter(!.data$valid_feature) |> nrow()
    if (n_invalid > 0){
      cli_alert_info(col_yellow("{n_invalid} invalid features (as defined in the metadata) were excluded"))
    }
  }

  # lists of analyses (analyses_id) and features (feature_id) in the dataset that are annotated
  annot_analyses <- data@annot_analyses  |> semi_join(data@dataset_orig, by = "analysis_id") |> pull(.data$analysis_id) |> unique()
  annot_features <- data@annot_features  |> semi_join(data@dataset_orig, by = "feature_id") |> pull(.data$feature_id) |> unique()


  # ISTD METADATA ====================
  if (!is.null(metadata$annot_istds) & nrow(metadata$annot_istds) > 0){
    # Add template table structure
    data@annot_istds <- metadata$annot_istds

    # Get info for cli output
    n_match <- intersect(data@annot_istds$quant_istd_feature_id, data@annot_features$quant_istd_feature_id) |> length()
    cli_alert_success(col_green(glue::glue("Internal Standard metadata associated with {n_match} ISTDs.")))
  }



  # RQC METADATA ====================
  if (!is.null(metadata$annot_responsecurves) & nrow(metadata$annot_responsecurves) > 0){
    # Add template table structure
    data@annot_responsecurves <- metadata$annot_responsecurves

    # Get info for cli output
    n_match <- intersect(annot_analyses, data@annot_responsecurves$analysis_id) |> length()
    cli_alert_success(col_green(glue::glue("Response curve metadata associated with {n_match} annotated analyses.")))
  }

  # lists of samples (sample_id) and analytes (analyte_id) in the dataset that are annotated
  annot_samples <- data@annot_analyses  |> semi_join(data@dataset_orig, by = "analysis_id") |> pull(.data$sample_id) |> unique()
  annot_analytes <- data@annot_features  |> semi_join(data@dataset_orig, by = "feature_id") |> pull(.data$analyte_id) |> unique()


  # QC CONCENTRATION METADATA ====================
  if (!is.null(metadata$annot_qcconcentrations) & nrow(metadata$annot_qcconcentrations) > 0){
    # Add template table structure
    data@annot_qcconcentrations <- metadata$annot_qcconcentrations

    # Get info for cli output
    n_match_samples <- intersect(annot_samples, unique(data@annot_qcconcentrations$sample_id)) |> length()
    n_match_analytes <- intersect(annot_analytes, unique(data@annot_qcconcentrations$analyte_id)) |> length()
    cli_alert_success(col_green(glue::glue("QC concentration metadata associated with {n_match_samples} annotated samples and {n_match_analytes} annotated analytes")))
  }


  # ANALYSIS ORDER ---------
  # Select timestamp if available, otherwise use order in raw data
  data <- set_analysis_order_analysismetadata(data)

  # BATCH METADATA ------------

  data@annot_batches <- get_metadata_batches(data@annot_analyses )

  # FINALIZE =================
  data@status_processing <- "Raw and metadata imported and associated"
  data
}



# stopifnot(methods::validObject(data, excl_nonannotated_analyses))

#' @title Reads and parses metadata provided by the msorganiser EXCEL  template.
#' @description Requires version 1.9.1 of the template
#' NOTES
#' - if no sample_type is defined then SPL will be assigned
#' - if valid_analysis is left blank for all analyses then samples then it will be replace by TRUE
#' - if valid_analysis is undefined for one or more, but all all samples, then an error will be returned
#' @param path File path of the msorganiser EXCEL template (*.xlm)
#' @param trim_ws Trim all white spaces and double spaces
#' @return A list with tibbles containing different metadata
#' @noRd
read_metadata_msorganiser <- function(path, trim_ws = TRUE) {

  metadata <- list()

  # check if path exist
  if (!(fs::path_ext(path) %in% c("xlsm", "xlsx"))){
    cli::cli_abort(col_red("Invalid file type not. A MSOrganiser template file (*.xlsm or *.xlsx) must be provided."))
  } else {
    if (!fs::file_exists(path))
      cli::cli_abort(col_red("File not found. Please verify path."))
  }


  w_xlm <- openxlsx2::wb_load(path)

  if (!"About" %in% openxlsx2::wb_get_sheet_names(w_xlm))
    cli::cli_abort(col_red("This appears to be an invalid or unsupported MSOrganiser template file, without `About` sheet. Please verify format and version."))

  d_about <- openxlsx2::wb_to_df(w_xlm,
                           sheet = "About",
                           skip_empty_rows = FALSE,
                           skip_empty_cols = FALSE,
                           skip_hidden_rows = FALSE,
                           skip_hidden_cols = FALSE,
                           convert = TRUE,
                           col_names = FALSE)

  if(!(str_detect(d_about[[1,2]], "MSOrganiser") && str_detect(d_about[[3,2]], "Version")))
      cli::cli_abort(col_red("This appears to be an invalid or unsupported MSOrganiser template file. Please verify format and version."))

  version <- str_replace(d_about[[3,3]], "(\\.[^.]*)\\.", "\\1")
  version <- suppressWarnings(as.numeric(version))
  if(is.na(version))
    cli::cli_abort(col_red("Invalid version number found in the template. Please use an MSOrganiser template v0.2 or higher."))


  if(version < 0.2 || version >= 0.3)
    cli::cli_abort(col_red("Unsupported MSOrganiser template version. Please use an MSOrganiser template v0.2 or higher."))

  w_xlm <- openxlsx2::wb_load(path)

  # ANALYSIS/SAMPLE annotation --------
  d_analyses <- openxlsx2::wb_to_df(w_xlm,
                              sheet = "Analyses (Samples)",
                              skip_empty_rows = FALSE,
                              skip_empty_cols = FALSE,
                              skip_hidden_rows = FALSE,
                              skip_hidden_cols = FALSE,
                              na.strings = c("", "''", "NA", "N/A", "#N/A", "n/a", "NaN", "nan"),
                              convert = TRUE,
                              col_names = TRUE) |>
    mutate(across(where(is.character), str_trim)) |>
    as_tibble() |>
    filter(!if_all(everything(), is.na))

  names(d_analyses) <- tolower(names(d_analyses))

  d_analyses <- d_analyses |>
    rename(
      qc_type = "sample_type",
      istd_volume = "istd_mixture_volume_[ul]"
    )

  # Ensure compatibility with older template versions (< 1.9.2)

  if("raw_data_filename" %in% names(d_analyses))
    d_analyses <- d_analyses |> rename("analysis_id" = "raw_data_filename")

  metadata$annot_analyses <- clean_analysis_metadata(d_analyses)


  # FEATURE annotation  -------------------------

  # ToDo: Make note if feature names are not original
  d_features <- openxlsx2::wb_to_df(file = w_xlm,
                                         sheet = "Features (Analytes)",
                                         skip_empty_rows = FALSE,
                                         skip_empty_cols = FALSE,
                                         skip_hidden_rows = FALSE,
                                         skip_hidden_cols = FALSE,
                                         na.strings = c("", "''", "NA", "N/A", "#N/A", "n/a", "NaN", "nan"),
                                         convert = TRUE,
                                         col_names = TRUE) |>
    mutate(across(where(is.character), str_trim)) |>
    as_tibble() |>
    filter(!if_all(everything(), is.na))
  names(d_features) <- tolower(names(d_features))

  # Ensure compatibility with older template versions (< 1.9.2) ---
  if("feature_name" %in% names(d_features))
    d_features <- d_features |> rename("feature_id" = "feature_name")

  if("istd_feature_name" %in% names(d_features))
    d_features <- d_features |> rename("istd_feature_id" = "istd_feature_name")

  if("interference_feature_name" %in% names(d_features))
    d_features <- d_features |> rename("interference_feature_id" = "interference_feature_name")

  if("interference_proportion" %in% names(d_features))
    d_features <- d_features |> rename("interference_contribution" = "interference_proportion")

  if("new_feature_name" %in% names(d_features))
    d_features <- d_features |> rename("feature_label" = "new_feature_name")

  if("quantifier" %in% names(d_features))
    d_features <- d_features |> rename("is_quantifier" = "quantifier")

  if("valid_integration" %in% names(d_features))
    d_features <- d_features |> rename("valid_feature" = "valid_integration")

  metadata$annot_features <- clean_feature_metadata(d_features)


  # ISTD annotation -------------------------

  d_istds <- openxlsx2::wb_to_df(file = w_xlm,
                                         sheet = "Internal Standards",
                                         cols = c(1,2,4),
                                         skip_empty_rows = FALSE,
                                         skip_empty_cols = FALSE,
                                         skip_hidden_rows = FALSE,
                                         skip_hidden_cols = FALSE,
                                         na.strings = c("", "''", "NA", "N/A", "#N/A", "n/a", "NaN", "nan"),
                                         col_names = TRUE) |>
    mutate(across(where(is.character), str_trim)) |>
    as_tibble() |>
    filter(!if_all(everything(), is.na))

  names(d_istds) <- tolower(names(d_istds))
  names(d_istds)[1] <- "istd_feature_id"

  # Ensure compatibility with older template versions (< 1.9.2) ---

  d_istds <- d_istds |>
    dplyr::select(any_of(c(
      istd_feature_id = "istd_feature_id",
      istd_conc_nmolar = "istd_conc_[nm]",
      istd_conc_ngml = "istd_conc_[ng/ml]")
      )
    ) |>
    dplyr::mutate(dplyr::across(where(is.character), stringr::str_squish))

  metadata$annot_istds <- clean_istd_metadata(d_istds)

  # RESPONSE CURVE annotation -------------------------

    d_rqc <- openxlsx2::wb_to_df(file = w_xlm,
                                         sheet = "Response Curves",
                                         skip_empty_rows = FALSE,
                                         skip_empty_cols = FALSE,
                                         skip_hidden_rows = FALSE,
                                         skip_hidden_cols = FALSE,
                                         na.strings = c("", "''", "NA", "N/A", "#N/A", "n/a", "NaN", "nan"),
                                         col_names = TRUE) |>
    mutate(across(where(is.character), str_trim)) |>
    as_tibble() |>
    filter(!if_all(everything(), is.na))

  names(d_rqc) <- tolower(names(d_rqc))

  # Ensure compatibility with older template versions (< 1.9.2)
  if("raw_data_filename" %in% names(d_rqc))
    d_rqc <- d_rqc |> rename("analysis_id" = "raw_data_filename")
  if("response_curve_name" %in% names(d_rqc))
    d_rqc <- d_rqc |> rename("curve_id" = "response_curve_name")

  if("relative_sample_amount_[%]" %in% names(d_rqc)){
    d_rqc <- d_rqc |> rename("analyzed_amount" = "relative_sample_amount_[%]")
    d_rqc <- d_rqc |> mutate("analyzed_amount_unit" = "%")
  }

  metadata$annot_responsecurves <- clean_response_metadata(d_rqc)

  # QC CONCENTRATION annotation -------------------------

  if ("QC Concentrations" %in% openxlsx2::wb_get_sheet_names((w_xlm))) {
    d_cal <- openxlsx2::wb_to_df(file = w_xlm,
                                 sheet = "QC Concentrations",
                                 skip_empty_rows = FALSE,
                                 skip_empty_cols = FALSE,
                                 skip_hidden_rows = FALSE,
                                 skip_hidden_cols = FALSE,
                                 na.strings = c("", "''", "NA", "N/A", "#N/A", "n/a", "NaN", "nan"),
                                 col_names = TRUE)
      d_cal <- d_cal |>
      mutate(across(where(is.character), str_trim)) |>
      as_tibble() |>
      filter(!if_all(everything(), is.na))

    names(d_cal) <- tolower(names(d_cal))

    metadata$annot_qcconcentrations <- clean_qcconc_metadata(d_cal)
  } else {
    metadata$annot_qcconcentrations <- NULL
  }



  # FINALIZE METADATA-------------------------
  metadata
}


# Read CSV or XLSX sheet
get_metadata_table <- function(dataset = NULL, path = NULL, sheet = NULL) {

  # Check if both `metadata` and `path` are provided
  if (!is.null(dataset) && !is.null(path)) {
    cli::cli_abort("Both 'dataset' and 'path' cannot be specified at the same time. Please provide either a data frame object or a file path.")  }

  if (is.null(dataset) && is.null(path)) {
    cli::cli_abort("Please provide either an existing data frame/tibble or a path to data file (SC. ")
  }

  # If `metadata` is provided, use it directly
  if (!is.null(dataset)) {
    if (!is.data.frame(dataset)) {
      cli::cli_abort("Dataset must be a data frame or tibble. Please provide a corresponding variable or define the path to a data file.")
    }
    d <- dataset
  } else {

    # Import from data file
    if (!fs::is_dir(path)) {
      file_path <- fs::path_tidy(path)
    } else {
      cli::cli_abort("`path` must be a file, not directory")
    }

    if(length(file_path) > 1) cli::cli_abort("Only one file can be imported.")

    if (!all(fs::file_exists(file_path))) cli::cli_abort("Files do not exist. Please verify file path.")

    if (fs::path_ext(file_path) == "csv") {

      d  <-  readr::read_csv(file_path, col_names = TRUE, trim_ws = TRUE, locale = readr::locale(encoding = "UTF-8"), col_types = "c", progress = FALSE)

    } else if (fs::path_ext(file_path) == "xlsx"){
      w_xlm <- openxlsx2::wb_load(file_path)
      if (!sheet %in% openxlsx2::wb_get_sheet_names(w_xlm))
        cli::cli_abort("Sheet `{sheet}` not found in the given file.")

      d <- openxlsx2::wb_to_df(file = w_xlm,
                               sheet = sheet,
                               skip_empty_rows = FALSE,
                               skip_empty_cols = FALSE,
                               skip_hidden_rows = FALSE,
                               skip_hidden_cols = FALSE,
                               convert = TRUE,
                               col_names = TRUE) |>
        mutate(dplyr::across(where(is.character), str_trim)) |>
        as_tibble() |>
        filter(!if_all(everything(), is.na))

    } else {
      cli::cli_abort("Only files with extension `csv` and `xlsx` are supported.")
    }
  d
  }
}


clean_analysis_metadata <- function(d_analyses) {

  names(d_analyses) <- tolower(names(d_analyses))

  if(!all(c("analysis_id") %in% names(d_analyses))){
    cli::cli_abort(cli::col_red("Analysis (Sample) metadata must have the `analysis_id` columnn. Please verify the input data. "))
  }

  # Fill missing columns
  d_analyses <- d_analyses |> add_missing_column(col_name = "analysis_order", init_value = NA_real_, make_lowercase = FALSE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "qc_type", init_value = NA_character_, make_lowercase = FALSE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "batch_id", init_value = 1, make_lowercase = FALSE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "sample_amount", init_value = NA_real_, make_lowercase = FALSE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "sample_amount_unit", init_value = NA_character_, make_lowercase = FALSE, all_na_replace = TRUE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "istd_volume", init_value = NA_real_, make_lowercase = FALSE, all_na_replace = TRUE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "valid_analysis", init_value = TRUE, make_lowercase = FALSE, all_na_replace = TRUE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "replicate_no", init_value = 1L, make_lowercase = FALSE, all_na_replace = TRUE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "specimen", init_value = NA_character_, make_lowercase = FALSE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "panel_id", init_value = NA_character_, make_lowercase = FALSE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "sample_id", init_value = NA_character_, make_lowercase = FALSE)
  d_analyses <- d_analyses |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)

  d_analyses <- d_analyses |>
    dplyr::select(
      "analysis_order",
      "analysis_id",
      "qc_type",
      "sample_amount",
      "sample_amount_unit",
      "istd_volume",
      "batch_id",
      "replicate_no",
      "specimen",
      "sample_id",
      "valid_analysis",
      "remarks"
    ) |>
    mutate(
      batch_id = str_squish(as.character(.data$batch_id)),
      remarks = str_squish(as.character(.data$remarks)),
      analysis_id = stringr::str_squish(as.character(.data$analysis_id)),
      analysis_id = stringr::str_remove(.data$analysis_id, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)),
      sample_id = stringr::str_squish(as.character(.data$sample_id)),
      replicate_no = as.integer(stringr::str_squish(as.character(.data$replicate_no))),
      specimen = stringr::str_squish(as.character(.data$specimen)),
      valid_analysis = as.logical(case_match(tolower(.data$valid_analysis),
                                             "yes" ~ TRUE,
                                             "no"~ FALSE,
                                             "true" ~ TRUE,
                                             "false" ~ FALSE,
                                             .default = NA)),
      qc_type = if_else(.data$qc_type == "Sample" | is.na(.data$qc_type), "SPL", .data$qc_type),
      annot_order_num = dplyr::row_number()) |>
    mutate(across(where(is.character), str_trim)) |>
    ungroup()
  # Handle the non-mandatory field Valid_Analysis
  if (all(is.na(d_analyses$valid_analysis)))
    d_analyses$valid_analysis <- TRUE
  else
    if (any(is.na(d_analyses$valid_analysis)) & !all(is.na(d_analyses$valid_analysis))){
        cli::cli_abort("`valid_analysis` is inconsistently defined, i.e., not for one or more analyses. Please verify imported analysis metadata.")
  }
  # Add template table structure
  d_analyses <- dplyr::bind_rows(pkg.env$table_templates$annot_analyses_template, d_analyses)
  d_analyses
  }

clean_feature_metadata <- function(d_features) {

  names(d_features) <- tolower(names(d_features))

  if("quantifier" %in% names(d_features))
    d_features <- d_features |> rename("is_quantifier" = "quantifier")


  if(!all(c("feature_id") %in% names(d_features))){
    cli::cli_abort(cli::col_red("Feature metadata must have column: `feature_id`. Please verify the input data. "))
  }
  d_features <- d_features |> add_missing_column(col_name = "feature_class", init_value = NA_character_, make_lowercase = FALSE, all_na_replace = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "chem_formula", init_value = NA_character_, make_lowercase = FALSE, all_na_replace = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "molecular_weight", init_value = NA_real_, make_lowercase = FALSE, all_na_replace = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "is_quantifier", init_value = TRUE, make_lowercase = FALSE, all_na_replace = TRUE)
  d_features <- d_features |> add_missing_column(col_name = "valid_feature", init_value = TRUE, make_lowercase = FALSE, all_na_replace = TRUE)
  d_features <- d_features |> add_missing_column(col_name = "istd_feature_id", init_value = NA_character_, make_lowercase = FALSE, all_na_replace = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "quant_istd_feature_id", init_value = NA_character_, make_lowercase = FALSE, all_na_replace = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "response_factor", init_value = 1.0, make_lowercase = FALSE, all_na_replace = TRUE)
  d_features <- d_features |> add_missing_column(col_name = "feature_label", init_value = NA_character_, make_lowercase = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "analyte_id", init_value = NA_character_, make_lowercase = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "interference_feature_id", init_value = NA_character_, make_lowercase = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "interference_contribution", init_value = NA_real_, make_lowercase = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "curve_fit_model", init_value = NA_character_, make_lowercase = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "curve_fit_weighting", init_value = NA_character_, make_lowercase = FALSE)
  d_features <- d_features |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)

  d_features <- d_features |>
    dplyr::mutate(
      feature_id = stringr::str_squish(.data$feature_id),
      feature_label = stringr::str_squish(.data$feature_label),
      feature_class = stringr::str_squish(.data$feature_class),
      chem_formula = stringr::str_squish(.data$chem_formula),
      analyte_id = stringr::str_squish(.data$analyte_id),
      istd_feature_id = stringr::str_squish(.data$istd_feature_id),
      quant_istd_feature_id = stringr::str_squish(.data$istd_feature_id),
      is_quantifier = as.logical(case_match(tolower(.data$is_quantifier),
                                            "yes" ~ TRUE,
                                            "no"~ FALSE,
                                            "true" ~ TRUE,
                                            "false" ~ FALSE,
                                            .default = NA)),
      valid_feature = as.logical(case_match(tolower(.data$valid_feature),
                                            "yes" ~ TRUE,
                                            "no"~ FALSE,
                                            "true" ~ TRUE,
                                            "false" ~ FALSE,
                                            .default = NA)),
      interference_feature_id = stringr::str_squish(.data$interference_feature_id),
      curve_fit_model = stringr::str_squish(.data$curve_fit_model),
      curve_fit_weighting = stringr::str_squish(.data$curve_fit_weighting),
      remarks = as.character(.data$remarks)
    ) |>
    mutate(across(where(is.character), str_trim)) |>
    dplyr::select(
      "feature_id",
      "feature_class",
      "feature_label",
      "chem_formula",
      "molecular_weight",
      "analyte_id",
      "istd_feature_id",
      "quant_istd_feature_id",
      "response_factor",
      "is_quantifier",
      "valid_feature",
      "interference_feature_id",
      "interference_contribution",
      "curve_fit_model",
      "curve_fit_weighting",
      "remarks"
    )


  # Handle the non-mandatory field Valid_Analysis  TODO: still needed?
  if (all(is.na(d_features$valid_feature)))
    d_features$valid_feature <- TRUE
  else
    if (any(is.na(d_features$valid_feature)))
      cli::cli_abort("`valid_feature` is inconsistently defined, i.e., not for one or more features. Please verify imported feature metadata.")

  # Handle the non-mandatory field  Quantifier TODO: still needed?
  if (all(is.na(d_features$is_quantifier)))
    d_features$is_quantifier <- TRUE
  else
    if (any(is.na(d_features$is_quantifier)))
      cli::cli_abort("`is_quantifier` is inconsistently defined, i.e., not for one or more features. Please verify imported feature metadata.")
  # Add template table structure
  d_features <- dplyr::bind_rows(pkg.env$table_templates$annot_features_template, d_features)
  d_features
}

clean_istd_metadata <- function(d_istds) {

  if(!all(c("istd_feature_id") %in% names(d_istds)) | !any(c("istd_conc_nmolar", "istd_conc_ngml") %in% names(d_istds))){
    cli::cli_abort(cli::col_red("ISTD metadata must have following columns: `istd_feature_id`, and `istd_conc_nmolar` or `istd_conc_ngml`. Please verify the input data. "))
  }

  d_istds <- d_istds |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)
  d_istds <- d_istds |>
    mutate(across(where(is.character), str_trim)) |>
    dplyr::mutate(
      istd_feature_id = as.character(.data$istd_feature_id),,
      remarks = as.character(.data$remarks)
    ) |>
    mutate(across(starts_with("feature_conc_"), ~ as.numeric(.))) |>
    mutate(across(where(is.character), str_squish)) |>
    dplyr::select( any_of(
      c("quant_istd_feature_id" = "istd_feature_id",
      "istd_conc_nmolar",
      "istd_conc_ngml",
      "remarks"))
    )

  # Add template table structure
  d_istds <- dplyr::bind_rows(pkg.env$table_templates$annot_istds_template, d_istds)
  d_istds
}



clean_response_metadata <- function(d_rqc) {
  if(!all(c("analysis_id", "curve_id", "analyzed_amount", "analyzed_amount_unit") %in% names(d_rqc))){
    cli::cli_abort(cli::col_red("Response curves metadata must have following columns: `analysis_id`, `curve_id`, `analyzed_amount` and `analyzed_amount_unit`. Please verify the input data. "))
  }

  d_rqc <- d_rqc |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)

  d_rqc <- d_rqc |>
    dplyr::mutate(
      analysis_id = stringr::str_remove(.data$analysis_id, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)),
      analysis_id = stringr::str_squish(as.character(.data$analysis_id)),
      curve_id = stringr::str_squish(as.character(.data$curve_id)),
      analyzed_amount = as.numeric(stringr::str_squish(.data$analyzed_amount)),
      analyzed_amount_unit = stringr::str_squish(.data$analyzed_amount_unit)
    ) |>
    dplyr::select(
      "analysis_id",
      "curve_id",
      "analyzed_amount",
      "analyzed_amount_unit",
      "remarks"
    )
  d_rqc <- dplyr::bind_rows(pkg.env$table_templates$annot_responsecurves_template,
                            d_rqc)
  d_rqc
}

clean_qcconc_metadata <- function(d_cal) {
  if(!all(c("sample_id", "analyte_id", "concentration", "concentration_unit") %in% names(d_cal))){
    cli::cli_abort(cli::col_red("QC concentration metadata must have following columns: `sample_id`, `analyte_id`, `concentration` and `concentration_unit`. Please verify the input data. "))
  }

  d_cal <- d_cal |> add_missing_column(col_name = "remarks", init_value = NA_character_, make_lowercase = FALSE)
  d_cal <- d_cal |> add_missing_column(col_name = "include_in_analysis", init_value = TRUE, make_lowercase = FALSE, all_na_replace = TRUE)


  d_cal <- d_cal |>
    dplyr::mutate(
      sample_id = stringr::str_remove(.data$sample_id, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)),
      sample_id = stringr::str_squish(as.character(.data$sample_id)),
      analyte_id = stringr::str_squish(as.character(.data$analyte_id)),
      concentration_unit = stringr::str_squish(.data$concentration_unit),
      include_in_analysis = as.logical(case_match(tolower(.data$include_in_analysis),
                                             "yes" ~ TRUE,
                                             "no"~ FALSE,
                                             "true" ~ TRUE,
                                             "false" ~ FALSE,
                                             NA ~ TRUE,
                                             .default = NA)),
    ) |>
    dplyr::select(
      "sample_id",
      "analyte_id",
      "concentration",
      "concentration_unit",
      "include_in_analysis",
      "remarks"
    )

  if (all(is.na(d_cal$include_in_analysis)))
    d_cal$include_in_analysis <- TRUE
  else
    if (any(is.na(d_cal$include_in_analysis))){
      cli::cli_abort("Invalid value(s) detected in `include_in_analysis`. Please check the QC concentration metadata, allowed values are 'true', 'false', 'yes', 'no', or empty (case-insensitive).")
    }
  d_cal <- dplyr::bind_rows(pkg.env$table_templates$annot_qcconcentrations_template,
                            d_cal)
  d_cal
}

