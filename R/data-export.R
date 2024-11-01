#' Writes all a data processing report to an EXCEL file
#'
#' @param data MidarExperiment object
#' @param path File name and path of the Excel file
#' @export
#'
#' @importFrom lubridate now
#' @importFrom tibble tribble
#' @importFrom utils packageVersion


# TODO: filtering of names containing "(IS"
report_write_xlsx <- function(data = NULL, path) {
  check_data(data)

  if (!stringr::str_detect(path, ".xlsx")) path <- paste0(path, ".xlsx")

  if(nrow(data@dataset) > 0 ){
    d_intensity_wide <- data@dataset |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "PBLK", "SBLK", "UBLK", "MBLK")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_intensity"))) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_intensity")
  } else {
    d_intensity_wide <- tibble("No annotated raw data available." = NA)
  }

  if(data@is_istd_normalized){
    d_norm_intensity_wide <- data@dataset |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "PBLK", "SBLK", "UBLK", "MBLK")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_norm_intensity"))) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_norm_intensity")
  } else {
    d_norm_intensity_wide <- tibble("No ISTD-normalized intensities available." = NA)
  }

  if(data@is_quantitated){
      d_conc_wide <- data@dataset |>
        dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) |>
        dplyr::filter(!str_detect(.data$feature_id, "\\(IS")) |>
        dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_conc"))) |>
        tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_conc")
  } else{
    d_conc_wide <- tibble("No concentration data available." = NA)
  }

  if (data@is_filtered) {
    d_conc_wide_QC_SPL <- data@dataset_filtered |>
      #dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "STD", "CTRL")) |>
      dplyr::filter(.data$qc_type %in% c("SPL")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "is_quantifier", "is_istd", "acquisition_time_stamp", "feature_id", "feature_conc"))) |>
      dplyr::filter(!str_detect(.data$feature_id, "\\(IS")) |>
      dplyr::filter(.data$is_quantifier) |>
      dplyr::filter(!.data$is_istd) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_conc")


    d_conc_wide_QC_all <- data@dataset_filtered |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "STD", "CTRL")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_norm_intensity"))) |>
      dplyr::filter(!str_detect(.data$feature_id, "\\(IS")) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_norm_intensity")

    d_conc_wide_QC_all <- d_conc_wide_QC_all |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "STD", "CTRL"))


  } else {
    d_conc_wide_QC_SPL <- tibble("No concentration data available." = NA)
    d_conc_wide_QC_all <- tibble("No qc-filtered concentration data available." = NA)
  }

  d_info <- tibble::tribble(
    ~Info, ~Value,
    "Date Report", as.character(lubridate::now()),
    "Author", Sys.info()[["user"]],
    "MiDAR Version", as.character(packageVersion("midar")[[1]]),
    "", "",
    "feature_conc Unit", get_conc_unit(data@annot_analyses$sample_amount_unit)
  )


  if(nrow(data@metrics_qc) == 0)
    qc_metrics <- tibble(tibble("Feature qc metrics available." = NA))
  else
    qc_metrics <- data@metrics_qc

  table_list <- list(
    "Info" = d_info,
    "Feature_QC_metrics" = qc_metrics,
    "Conc_QCfilt_StudySamples" = d_conc_wide_QC_SPL,
    "Conc_QCfilt_AllSamples" = d_conc_wide_QC_all,
    "Conc_FullDataset" = d_conc_wide,
    "Raw_Intensity_FullDataset" = d_intensity_wide,
    "Norm_Intensity_FullDataset" = d_norm_intensity_wide,
    "SampleMetadata" = data@annot_analyses,
    "FeatureMetadata" = data@annot_features,
    "InternalStandards" = data@annot_istd,
    "BatchInfo" = data@annot_batches
  )

  cat("\rSaving report to disk - please wait...")
  openxlsx2::write_xlsx(x = table_list,
                        file = path,
                        na.strings = "",
                        as_table = TRUE,
                        overwrite = TRUE,
                        col_names = TRUE,
                        grid_lines = FALSE,
                        col_widths = "auto",
                        first_col = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                        first_row = c(FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
                        with_filter = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE),
                        tab_color = c("#d7fc5d", "#34fac5","#ff170f", "#9e0233", "#0A83ad", "#0313ad","#7113ad", "#c9c9c9", "#c9c9c9", "#c9c9c9", "#c9c9c9")
  )
  cli_alert_success(col_green(glue::glue("\rThe data processing report of analysis '{data@title}' has been saved.")))
}



#' Export any parameter to a wide-format table
#'
#' @param data MidarExperiment object
#' @param path File name with path of exported CSV fil
#' @param variable Variable to be exported
#' @param filter_data Use QC-filtered data, based on criteria set via `qc_apply_feature_filter()`. Overwrites `include_qualifier` and `include_istd`.
#' @param qc_types QC type to plot. When qc_types us NA or NULL, all available QC types are plotted.
#' @param include_qualifier Include qualifier features. Is not used when `filter_data = TRUE` was applied.
#' @param include_istd Include internal standard features. Default is `TRUE`. Is not used when `filter_data = TRUE` was applied.
#' @param add_qctype Add the QC type as column
#' @export
report_write_csv <- function(data = NULL,
                             path,
                             variable = c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm"),
                             filter_data,
                             qc_types = c("SPL", "BQC", "TQC", "NIST", "LTR", "PBLK", "SBLK", "UBLK", "MBLK"),
                             include_qualifier = NA,
                             include_istd = NA,
                             add_qctype = NA
                             ) {
  check_data(data)
  variable_strip <- str_remove(variable, "feature_")
  rlang::arg_match(variable_strip, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm"))
  variable <- stringr::str_c("feature_", variable_strip)
  check_var_in_dataset(data@dataset, variable)
  variable_sym = rlang::sym(variable)

  # Auto-choose some arg values if user does not define

  if(is.na(include_qualifier)) {
    if(variable %in% c("feature_conc", "feature_conc_raw")) include_qualifier <- FALSE else include_qualifier <- TRUE
  }

  if(is.na(include_istd)) {
    if(variable %in% c("feature_conc", "feature_conc_raw", "feature_norm_intensity")) include_istd <- FALSE else include_istd <- TRUE
  }

  if(length(qc_types) == 1) add_qctype <- FALSE else add_qctype <- TRUE

  if (!(variable %in% names(data@dataset))) cli::cli_abort("Variable '{variable}' has not yet been calculated. Please process data or choose other variable.")



  # Filter data if filter_data is TRUE
  if (filter_data) {
    dat_filt <- data@dataset_filtered |> dplyr::ungroup()
    if (!data@is_filtered) cli::cli_abort("Data has not been qc filtered. Please apply `qc_apply_feature_filter` first.")
  } else {
    dat_filt <- data@dataset |> dplyr::ungroup()
  }

  if(!include_qualifier){
    dat_filt <- dat_filt |> filter(.data$is_quantifier)
  }

  if(!include_istd){
    dat_filt <- dat_filt |> filter(!.data$is_istd)
  }


  # Subset data based on qc_types argument ----
  if (all(!is.na(qc_types)) & all(qc_types != "")) {
    if (length(qc_types) == 1) {
      dat_filt <- dat_filt |> dplyr::filter(stringr::str_detect(.data$qc_type, qc_types))
    } else {
      dat_filt <- dat_filt |> dplyr::filter(.data$qc_type %in% qc_types)
    }
  }

  if(add_qctype)
    flds <- c("analysis_id", "qc_type", "feature_id")
  else
    flds <- c("analysis_id", "feature_id")


  ds <- dat_filt |>
    dplyr::select(all_of(c(flds, variable))) |>
    tidyr::pivot_wider(names_from = .data$feature_id, values_from = !!variable_sym)

  readr::write_csv(ds, file = path, num_threads = 4, col_names = TRUE)
  if(variable_strip ==  "conc") variable_strip <- "concentration"
  cli_alert_success(col_green(glue::glue("{stringr::str_to_title(variable_strip)} values of {nrow(ds)} analyses and {length(unique(dat_filt$feature_id))} features have been exported.")))
}


#' Save the QC table to a CSV file
#'
#' @param data MidarExperiment object
#' @param path File name with path of exported CSV file
#' @return A tibble with the exported dataset
#' @export

report_write_qc_metrics <- function(data = NULL, path) {
  check_data(data)
  if (nrow(data@metrics_qc) == 0) cli::cli_abort("QC info has not yet been calculated. Please apply 'qc_calc_metrics' first.")

  readr::write_csv(data@metrics_qc, file = path, num_threads = 4, col_names = TRUE)
  invisible(data@metrics_qc)
}
