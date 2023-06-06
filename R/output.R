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

  if (!("Concentration" %in% names(data@dataset))) stop("Variable '", "Concentration",  "' does not (yet) exist in dataset")
  if (!stringr::str_detect(filename, ".xlsx")) filename = paste0(filename, ".xlsx")

  d_intensity_wide <- data@dataset %>%
    dplyr::filter(.data$QC_TYPE %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) %>%
    dplyr::select(dplyr::any_of(c("ANALYSIS_ID", "QC_TYPE", "AcqTimeStamp", "FEATURE_NAME", "Intensity"))) %>%
    tidyr::pivot_wider(names_from = "FEATURE_NAME", values_from = "Intensity")

  d_conc_wide <- data@dataset %>%
    dplyr::filter(.data$QC_TYPE %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) %>%
    dplyr::filter(!str_detect(.data$FEATURE_NAME, "\\(IS")) %>%
    dplyr::select(dplyr::any_of(c("ANALYSIS_ID", "QC_TYPE", "AcqTimeStamp", "FEATURE_NAME", "Concentration"))) %>%
    tidyr::pivot_wider(names_from = "FEATURE_NAME", values_from = "Concentration")

  if("FEATURE_NAME" %in% names(data@dataset_QC_filtered)) {

    d_conc_wide_QC <- data@dataset_QC_filtered %>%
      dplyr::filter(.data$QC_TYPE %in% c("SPL", "TQC", "BQC", "NIST", "LTR")) %>%
      dplyr::select(dplyr::any_of(c("ANALYSIS_ID", "QC_TYPE", "isISTD.x", "isQUANTIFIER", "AcqTimeStamp", "FEATURE_NAME", "Concentration"))) %>%
      dplyr::filter(!str_detect(.data$FEATURE_NAME, "\\(IS")) %>%
      tidyr::pivot_wider(names_from = "FEATURE_NAME", values_from = "Concentration")

    d_conc_wide_QC_SPL <- d_conc_wide_QC |> dplyr::filter(.data$QC_TYPE == "SPL") |> dplyr::select(!"QC_TYPE":"isQUANTIFIER")

  } else {
    d_conc_wide_QC <- data@dataset_QC_filtered
  }

  d_info <- tibble::tribble(
    ~Info, ~Value,
    "Date Report", lubridate::now(),
    "Author", Sys.info()[["user"]],
    "LIDAR Version", packageVersion("midar") ,
    "", "",
    "Concentration Unit", get_conc_unit(data@annot_analyses$SAMPLE_AMOUNT_UNIT))


 table_list <- list(
            "Intensities_All" = d_intensity_wide,
            "Conc_All" = d_conc_wide,
            "Conc_QCfilt" = d_conc_wide_QC,
            "QC" = data@metrics_qc,
            "Info" = d_info,
            "SampleMetadata" = data@annot_analyses,
            "FeatureMetadata" = data@annot_features,
            "InternalStandards" = data@annot_istd,
            "BatchInfo" = data@annot_batch_info)

  openxlsx::write.xlsx(x = table_list, file = filename, overwrite = TRUE )
}
