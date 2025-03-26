#' Write Data Processing Report (EXCEL)
#'
#' Generates a data processing report from a `MidarExperiment` object and writes it to an Excel file.
#' The report includes information on the data processing steps, quality control metrics, feature concentrations, and metadata.
#' Following tables will be created as sheets in the EXCEL file:
#'
#' - Info: General  information including date, author, and MiDAR version, processing status and feature concentration unit.
#' - Feature_QC_metrics: Quality control metrics of all features.
#' - QCfilt_x_StudySamples: Feature (QC)-filtered data (variable defiend via `filtered_variable`) in study samples ('SPL'). Filter have to be set via [filter_features_qc()]. The _x_ corresponds to the `filtered_variable` argument.
#' - QCfilt_x_AllSamples: Feature (QC)-filtered data (variable defiend via `filtered_variable`) in all samples. Filter have to be set via [filter_features_qc()]. The _x_ corresponds to the `filtered_variable` argument.
#' - Conc_FullDataset: Final feature concentrations from the full, non-filtered dataset.
#' - Raw_Intensity_FullDataset: Raw feature intensities from the full, non-filtered dataset.
#' - Norm_Intensity_FullDataset: Normalized feature intensities from the full, non-filtered dataset.
#' - SampleMetadata:  Analysis metadata that was imported and used for processing steps
#' - FeatureMetadata: Feature metadata that was imported and used for processing steps
#' - InternalStandards: Internal standards metadata with concentrations
#' - BatchInfo: Information on batches and positions of first and last analysis/sample
#' in each batch
#'
#'
#' @param data A `MidarExperiment` object containing original and processed data and metadata.
#' @param path A character string specifying the file name and path for the Excel file.
#' If the path does not include an `.xlsx` extension, it is added automatically.
#' @param filtered_variable A character string specifying the variable name in the
#' filtered data to be exported. It must be one of "conc", "intensity", "norm_intensity",
#' "response", "area", "height", "conc_raw", "rt", or "fwhm". The defined variable
#' name will be included in the sheet name. Default is "conc".
#'
#' @details
#' #' If certain data sets are not available, the function includes empty tables for the corresponding dataset.
#'
#' Concentration corresponds to the final concentration values after applying isotope correction, and drift and batch correction, if applicable.
#' If any corrections, such as drift or batch correction, were applied to raw or normalized intensities, the exported values will reflect these corrections.
#'
#' @return The function does not return a value. It writes the report to the specified Excel file.
#'
#' @examples
#' \dontrun{
#' # Assuming `midarexp` is a MidarExperiment object and `output_path` is a valid path
#' save_report_xlsx(data = midarexp, path = "output_path/report.xlsx")
#' }
#'
#' @export

save_report_xlsx <- function(data = NULL, path, filtered_variable = "conc") {
  check_data(data)

  filtered_variable <- str_remove(filtered_variable, "feature_")
  filtered_variable_strip <- filtered_variable
  rlang::arg_match(filtered_variable, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm"))
  filtered_variable <- stringr::str_c("feature_", filtered_variable)
  if(data@is_filtered) check_var_in_dataset(data@dataset, filtered_variable) #TODO dataset_filt?

  if (!stringr::str_detect(path, ".xlsx")) path <- paste0(path, ".xlsx")

  if(nrow(data@dataset) > 0 ){
    d_intensity_wide <- data@dataset |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "PBLK", "SBLK", "UBLK", "MBLK")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_intensity"))) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_intensity")
  } else {
    d_intensity_wide <- tibble("No ISTD-normalized intensities available." = NA) |> tibble::add_row()
  }

  if(data@is_istd_normalized){
    d_norm_intensity_wide <- data@dataset |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "PBLK", "SBLK", "UBLK", "MBLK")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_norm_intensity"))) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_norm_intensity")
  } else {
    d_norm_intensity_wide <- tibble("No ISTD-normalized intensities available." = NA) |> tibble::add_row()
  }

  if(data@is_quantitated){
      d_conc_wide <- data@dataset |>
        dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) |>
        dplyr::filter(!str_detect(.data$feature_id, "\\(IS")) |>
        dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_conc"))) |>
        tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_conc")
  } else{
    d_conc_wide <- tibble("No concentration data available." = NA) |> tibble::add_row()
  }

  if (data@is_filtered) {
    d_conc_wide_QC_SPL <- data@dataset_filtered |>
      #dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "STD", "CTRL")) |>
      dplyr::filter(.data$qc_type %in% c("SPL")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "feature_id", filtered_variable))) |>
      dplyr::filter(!str_detect(.data$feature_id, "\\(IS")) |>
      #dplyr::filter(.data$is_quantifier) |>
      #dplyr::filter(!.data$is_istd) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = any_of(filtered_variable))


    d_conc_wide_QC_all <- data@dataset_filtered |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "STD", "CTRL")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "feature_id", "feature_norm_intensity"))) |>
      dplyr::filter(!str_detect(.data$feature_id, "\\(IS")) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_norm_intensity")

    d_conc_wide_QC_all <- d_conc_wide_QC_all |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "STD", "CTRL"))


  } else {
    filtered_variable_strip <- ""
    d_conc_wide_QC_SPL <- tibble("No qc-filtered data available." = NA) |> tibble::add_row()
    d_conc_wide_QC_all <- tibble("No qc-filtered data available." = NA) |> tibble::add_row()
  }

  if(data@is_quantitated && data@status_processing == "Calibration-quantitated data")
    conc_unit_origin <- unique(data@annot_qcconcentrations$concentration_unit)
  else
    conc_unit_origin <- "pmpol"

  d_info <- tibble::tribble(
    ~Info, ~Value,
    "Date Report", as.character(lubridate::now()),
    "Author", Sys.info()[["user"]],
    "MiDAR Version", as.character(utils::packageVersion("midar")[[1]]),
    "", "",
    "feature_conc Unit", get_conc_unit(data@annot_analyses$sample_amount_unit, conc_unit_origin)
  )


  if(nrow(data@metrics_qc) == 0)
    qc_metrics <- tibble(tibble("Feature qc metrics has not been calculated." = NA)) |> tibble::add_row()
  else
    qc_metrics <- data@metrics_qc

  if(nrow(data@metrics_calibration) == 0)
    metrics_calibration <- tibble(tibble("Calibration metrics has not been calculated." = NA)) |> tibble::add_row()
  else
    metrics_calibration <- data@metrics_calibration

  if(filtered_variable_strip == "norm_intensity") filtered_variable_strip <- "normInt"
  if(filtered_variable_strip != "") {
    name_filt<- paste0("_", paste0(toupper(substr(filtered_variable_strip, 1, 1)), substr(filtered_variable_strip, 2, nchar(filtered_variable_strip))) )
  } else {
    name_filt <- ""
  }
  name_filt_spl <- paste0("QCfilt",name_filt,"_StudySamples")
  name_filt_all <- paste0("QCfilt",name_filt,"_AllSamples")


  table_list <- list(
    "Info" = d_info,
    "Feature_QC_metrics" = qc_metrics,
    "Calibration_metrics" = metrics_calibration,
    name_filt_spl = d_conc_wide_QC_SPL,
    name_filt_all = d_conc_wide_QC_all,
    "Raw_Intensity_FullDataset" = d_intensity_wide,
    "Norm_Intensity_FullDataset" = d_norm_intensity_wide,
    "Conc_FullDataset" = d_conc_wide,
    "SampleMetadata" = if(nrow(data@annot_analyses) == 0) data@annot_analyses |> tibble::add_row() else data@annot_analyses,
    "FeatureMetadata" = if(nrow(data@annot_features) == 0) data@annot_features |> tibble::add_row() else data@annot_features,
    "InternalStandards" = if(nrow(data@annot_istds) == 0) data@annot_istds |> tibble::add_row() else data@annot_istds,
    "BatchInfo" = if(nrow(data@annot_batches) == 0) tibble("No batches defined" = NA) |> tibble::add_row() else data@annot_batches
  )

  names(table_list)[4:5] <- c(name_filt_spl, name_filt_all)

  message("\rSaving report to disk - please wait...")
  openxlsx2::write_xlsx(x = table_list,
                        file = path,
                        na.strings = "",
                        as_table = TRUE,
                        overwrite = TRUE,
                        col_names = TRUE,
                        grid_lines = FALSE,
                        col_widths = "auto",
                        first_col = c(FALSE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                        first_row = c(FALSE, TRUE, TRUE,TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                        with_filter = c(FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
                        tab_color = c("#d7fc5d", "#34fac5","#34fac5","#ff170f", "#9e0233", "#0A83ad", "#0313ad","#7113ad", "#c9c9c9", "#c9c9c9", "#c9c9c9", "#c9c9c9")
  )

  txtitle <- if(data@title != "") glue::glue(" of experiment '{data@title}' ") else " "
  cli_alert_success(col_green(glue::glue("\rThe data processing report{txtitle}has been saved to '{path}'.")))
}


#' Export Data to CSV file
#'
#' This function exports specific unprocessed or pr ocessed feature variable
#' (e.g. intensities or concentrations) from a `MidarExperiment` object to a CSV file.
#' Allows selection of features and optional QC filtering.
#' @param data MidarExperiment object
#' @param path File name with path of exported CSV file
#' @param variable Variable to be exported, must be present in the data and any of "area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm".
#' @param filter_data A logical value indicating whether to use all data
#' (default) or only QC-filtered data (filtered via [filter_features_qc()]).
#' @param qc_types QC types to be plotted. Can be a vector of QC types or a regular expression pattern. `NA` (default) displays all available QC/Sample types.
#' @param include_qualifier A logical value indicating whether to include
#' qualifier features. Default is `NA`, which will be automatically set to `FALSE`
#' if `variable` is `conc` or `conc_raw`, and `FALSE` otherwise.
#' @param include_istd A logical value indicating whether to include internal
#' standard (ISTD) features. Default is `NA`, which will be automatically set to `FALSE`
#' if `variable` is ''norm_intensity`, `conc` or `conc_raw`, and `TRUE` otherwise.
#' @param include_feature_filter A character or regex pattern used to filter
#' features by `feature_id`. If `NA` or an empty string (`""`) is provided,
#' the filter is ignored. When a vector of length > 1 is supplied, only
#' features with exactly these names are selected (applied individually as
#' OR conditions).
#' @param exclude_feature_filter A character or regex pattern used to exclude
#' features by `feature_id`. If `NA` or an empty string (`""`) is provided,
#' the filter is ignored. When a vector of length > 1 is supplied, only
#' features with exactly these names are excluded (applied individually as
#' OR conditions).
#' @param add_qctype Add the QC type as column
#' @export
save_dataset_csv <- function(data = NULL,
                             path,
                             variable,
                             filter_data,
                             qc_types = NA,
                             include_qualifier = NA,
                             include_istd = NA,
                             include_feature_filter = NA,
                             exclude_feature_filter = NA,
                             add_qctype = NA
                             ) {
  check_data(data)
  variable <- str_remove(variable, "feature_")
  variable_strip <- variable
  rlang::arg_match(variable, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm"))
  variable <- stringr::str_c("feature_", variable)
  check_var_in_dataset(data@dataset, variable)
  variable_sym = rlang::sym(variable)

  # Auto-choose some arg values if user does not define

  if(is.na(include_qualifier)) {
    if(variable %in% c("feature_conc", "feature_conc_raw")) include_qualifier <- FALSE else include_qualifier <- TRUE
  }

  if(is.na(include_istd)) {
    if(variable %in% c("feature_conc", "feature_conc_raw", "feature_norm_intensity")) include_istd <- FALSE else include_istd <- TRUE
  }

  if(is.na(add_qctype))
    add_qctype <- !(length(qc_types) == 1)

  if (!(variable %in% names(data@dataset))) cli::cli_abort("Variable '{variable}' has not yet been calculated. Please process data or choose other variable.")

  if(all(is.na(qc_types))){
    qc_types <- unique(data$dataset$qc_type)
  }


  # Subset dataset according to arguments
  d_filt <- get_dataset_subset(
    data,
    filter_data = filter_data,
    qc_types = qc_types,
    include_qualifier = include_qualifier,
    include_istd = include_istd,
    include_feature_filter = include_feature_filter,
    exclude_feature_filter = exclude_feature_filter
  )

  if(add_qctype)
    flds <- c("analysis_id", "qc_type", "feature_id")
  else
    flds <- c("analysis_id", "feature_id")

  ds <- d_filt |>
    dplyr::select(all_of(c(flds, variable))) |>
    tidyr::pivot_wider(names_from = "feature_id", values_from = !!variable_sym)

  readr::write_csv(ds, file = path, col_names = TRUE)
  if(variable_strip ==  "conc") variable_strip <- "concentration"
  cli_alert_success(col_green(glue::glue("{stringr::str_to_title(variable_strip)} values for {nrow(ds)} analyses and {length(unique(d_filt$feature_id))} features have been exported to '{path}'.")))
}

#' Save Feature QC Metrics to CSV
#'
#' This function exports the feature information and QC (Quality Control) metrics
#' from a MidarExperiment object to a CSV file.
#'
#' @param data A MidarExperiment object containing the QC metrics.
#' @param path A string specifying the file path where the CSV file will be saved.
#' @return A tibble with the QC metrics that have been exported.
#' @export
#'
save_feature_qc_metrics <- function(data = NULL, path) {
  check_data(data)

  # Verify that the QC metrics have been calculated
  if (nrow(data@metrics_qc) == 0) {
    cli::cli_abort(col_red("Feature QC metrics has not yet been calculated. Please run 'calc_qc_metrics()' first."))
  }

  # Write the QC metrics to a CSV file
  readr::write_csv(data@metrics_qc, file = path, col_names = TRUE)

  cli_alert_success(col_green("Feature QC metrics table was saved to '{path}'."))

  # Return the QC metrics invisibly as a side-effect
  invisible(data@metrics_qc)
}

#' Saves a Excel (xlsx) file with metadata templates
#'
#' This function saves a XLSX file with metadata template to the specified location.
#'
#' @param path File path where the XLSX file with templates will be saved.
#' If left empty (default), the file will be saved in the current working directory
#' under the file "metadata_template.xlsx"
#' @export
#'

save_metadata_templates <- function(path = "metadata_template.xlsx") {
  # Locate the template file inside the package
  template_path <- system.file("extdata", "midar_metadata_templates.xlsx", package = "midar")

  if (fs::file_exists(path)) {
    cli_abort(col_red("A file with this name already exists at the specified location. Please delete it or choose a different filename or location."))
  }

  if (template_path == "") {
    cli_abort(col_red("Template file not found in package data. Please re-install `midar`."))
  }

  # Copy the template to the desired location
  fs::file_copy(template_path, path, overwrite = TRUE)
  cli_alert_success(col_green("Metadata table templates were saved to '{path}'."))
}


#' Saves a MiDAR Metadata Organizer template
#'
#' This function saves a XLSX file with metadata template to the specified location.
#'
#' @param path File path where the MiDAR Metadata Organizer file will be saved.
#' If left empty (default), the file will be saved in the current working directory
#' under the file "metadata_msorganiser_template.xlsm"
#' @export
#'

save_metadata_msorganiser_template <- function(path = "metadata_msorganiser_template.xlsm") {
  # Locate the template file inside the package
  template_path <- system.file("extdata", "metadata_msorganiser_template.xlsm", package = "midar")

  if (fs::file_exists(path)) {
    cli_abort(col_red("A file with this name already exists at the specified location. Please delete it or choose a different filename or location."))
  }

  if (template_path == "") {
    cli_abort(col_red("Template file not found in package data. Please re-install `midar`."))
  }

  # Copy the template to the desired location
  fs::file_copy(template_path, path, overwrite = TRUE)
  cli_alert_success(col_green("A MiDAR Metadata Organizer template was saved to '{path}'."))
}
