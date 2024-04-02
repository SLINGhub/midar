#' Writes all a data processing report to an EXCEL file
#'
#' @param data MidarExperiment object
#' @param filename File name and path of the Excel file
#' @export
#'
#' @importFrom glue glue
#' @importFrom openxlsx write.xlsx
#' @importFrom lubridate now
#' @importFrom tibble tribble
#' @importFrom utils packageVersion

#'
writeReportXLS <- function(data, filename) {

  if (!("feature_conc" %in% names(data@dataset))) stop("Variable '", "feature_conc",  "' does not (yet) exist in dataset")
  if (!stringr::str_detect(filename, ".xlsx")) filename = paste0(filename, ".xlsx")

  d_intensity_wide <- data@dataset %>%
    dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) %>%
    dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_name",  "feature_intensity"))) %>%
    tidyr::pivot_wider(names_from = "feature_name", values_from = "feature_intensity")

  d_conc_wide <- data@dataset %>%
    dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) %>%
    dplyr::filter(!str_detect(.data$feature_name, "\\(IS")) %>%
    dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "acquisition_time_stamp", "feature_name",  "feature_conc"))) %>%
    tidyr::pivot_wider(names_from = "feature_name", values_from = "feature_conc")

  if("feature_name" %in% names(data@dataset_filtered)) {

    d_conc_wide_QC <- data@dataset_filtered %>%
      dplyr::filter(.data$qc_type %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) %>%
      dplyr::select(dplyr::any_of(c("analysis_id", "qc_type", "is_istd.x", "is_quantifier", "acquisition_time_stamp", "feature_name",  "feature_conc"))) %>%
      dplyr::filter(!str_detect(.data$feature_name, "\\(IS")) %>%
      tidyr::pivot_wider(names_from = "feature_name", values_from = "feature_conc")

    d_conc_wide_QC_SPL <- d_conc_wide_QC |> dplyr::filter(.data$qc_type == "SPL") |> dplyr::select(!"qc_type":"is_quantifier")

  } else {
    d_conc_wide_QC <- data@dataset_filtered
  }

  d_info <- tibble::tribble(
    ~Info, ~Value,
    "Date Report", lubridate::now(),
    "Author", Sys.info()[["user"]],
    "LIDAR Version", packageVersion("midar") ,
    "", "",
    "feature_conc Unit", get_conc_unit(data@annot_analyses$sample_amount_unit))


 table_list <- list(
            "Intensities_All" = d_intensity_wide,
            "Conc_All" = d_conc_wide,
            "Conc_QCfilt" = d_conc_wide_QC,
            "QC" = data@metrics_qc,
            "Info" = d_info,
            "SampleMetadata" = data@annot_analyses,
            "FeatureMetadata" = data@annot_features,
            "InternalStandards" = data@annot_istd,
            "BatchInfo" = data@annot_batches)

  openxlsx::write.xlsx(x = table_list, file = filename, overwrite = TRUE )
}
