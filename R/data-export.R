#' Writes all a data processing report to an EXCEL file
#'
#' @param data MidarExperiment object
#' @param path File name and path of the Excel file
#' @export
#'
#' @importFrom glue glue
#' @importFrom openxlsx write.xlsx
#' @importFrom lubridate now
#' @importFrom tibble tribble
#' @importFrom utils packageVersion

#'


# TODO: filtering of names containing "(IS"
report_write_xlsx <- function(data, path) {
  if (!("feature_conc" %in% names(data@dataset))) cli::cli_abort("Variable '", "feature_conc", "' does not (yet) exist in dataset")
  if (!stringr::str_detect(path, ".xlsx")) path <- paste0(path, ".xlsx")

  d_intensity_wide <- data@dataset |>
    dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) |>
    dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_intensity"))) |>
    tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_intensity")

  d_conc_wide <- data@dataset |>
    dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) |>
    dplyr::filter(!str_detect(.data$feature_id, "\\(IS")) |>
    dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_conc"))) |>
    tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_conc")

  if ("feature_id" %in% names(data@dataset_filtered)) {
    d_conc_wide_QC <- data@dataset_filtered |>
      #dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "STD", "CTRL")) |>
      dplyr::filter(.data$qc_type %in% c("SPL")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "is_quantifier", "is_istd", "acquisition_time_stamp", "feature_id", "feature_conc"))) |>
      dplyr::filter(!str_detect(.data$feature_id, "\\(IS")) |>
      dplyr::filter(.data$is_quantifier) |>
      dplyr::filter(!.data$is_istd) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_conc")

    d_conc_wide_QC_SPL <- d_conc_wide_QC |>
      dplyr::filter(.data$qc_type == "SPL") |>
      dplyr::select(!"qc_type":"acquisition_time_stamp")
  } else {
    d_conc_wide_QC <- NULL
  }


  if ("feature_id" %in% names(data@dataset_filtered)) {
    d_conc_wide_QC_all <- data@dataset_filtered |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "STD", "CTRL")) |>
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", "feature_norm_intensity"))) |>
      dplyr::filter(!str_detect(.data$feature_id, "\\(IS")) |>
      tidyr::pivot_wider(names_from = "feature_id", values_from = "feature_norm_intensity")

    d_conc_wide_QC_all <- d_conc_wide_QC_all |>
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR", "STD", "CTRL"))
      #dplyr::select(!"qc_type":"acquisition_time_stamp")
  } else {
    d_conc_wide_QC_all <- NULL
  }

  d_info <- tibble::tribble(
    ~Info, ~Value,
    "Date Report", as.character(lubridate::now()),
    "Author", Sys.info()[["user"]],
    "MiDAR Version", as.character(packageVersion("midar")[[1]]),
    "", "",
    "feature_conc Unit", get_conc_unit(data@annot_analyses$sample_amount_unit)
  )


  table_list <- list(
    "Info" = d_info,
    "Feature_QC_metrics" = data@metrics_qc,
    "Conc_QCfilt_StudySamples" = d_conc_wide_QC,
    "Conc_QCfilt_AllSamples" = d_conc_wide_QC_all,
    "Conc_FullDataset" = d_conc_wide,
    "RawIntensiy_FullDataset" = d_intensity_wide,
    "SampleMetadata" = data@annot_analyses,
    "FeatureMetadata" = data@annot_features,
    "InternalStandards" = data@annot_istd,
    "BatchInfo" = data@annot_batches
  )

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
                        tab_color = c("#d7fc5d", "#34fac5", "#0383ad", "#0383ad", "#9e0233", "#ff170f", "#c9c9c9", "#c9c9c9", "#c9c9c9", "#c9c9c9")
  )
}



#' Export any parameter to a wide-format table
#'
#' @param data MidarExperiment object
#' @param variable Variable to be exported
#' @param path File name with path of exported CSV file
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider
#' @export
exportWideCSV <- function(data, variable, path) {
  var <- dplyr::sym(variable)

  if (!(variable %in% names(data@dataset))) cli::cli_abort("Variable '", variable, "' does not (yet) exist in dataset.")

  ds <- data@dataset |>
    dplyr::select("analysis_id", "qc_type", "acquisition_time_stamp", "feature_id", !!var) |>
    tidyr::pivot_wider(names_from = .data$feature_id, values_from = !!var)

  readr::write_csv(ds, file = path, num_threads = 4, col_names = TRUE)
  invisible(ds)
}


#' Save the QC table to a CSV file
#'
#' @param data MidarExperiment object
#' @param path File name with path of exported CSV file
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider
#' @export

saveQCinfo <- function(data, path) {
  if (nrow(data@metrics_qc) == 0) cli::cli_abort("QC info has not yet been calculated. Please apply 'qc_calc_metrics' first.")

  readr::write_csv(data@metrics_qc, file = path, num_threads = 4, col_names = TRUE)
  invisible(data@metrics_qc)
}
