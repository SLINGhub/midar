#' Imports metadata provided by an MSOrganizer EXCEL template
#'
#' @param filename File path of the MSOrganizer EXCEL template (*.xlm)
#' @param trim_ws Trim all white spaces and double spaces
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom tidyselect vars_select_helpers
#' @importFrom dplyr select mutate filter group_by row_number
#' @importFrom dplyr recode
#' @return A list of tibbles with different metadata
#' @export
import_msorganizer_xlm <- function(filename, trim_ws = TRUE){
  d_annot <- list()
  d_annot$annot_analyses <- readxl::read_excel(filename, sheet = "Sample_Annot") |>
    dplyr::mutate(
      VALID_ANALYSIS = TRUE,
      BATCH_ID = as.character(.data$BATCH_ID),
      RUN_ID_ANNOT = dplyr::row_number()) |>
    dplyr::select(
      "RUN_ID_ANNOT",
      ANALYSIS_ID = "Sample_Name",
      DATAFILE_NAME = "Sample_Name",
      QC_TYPE = "Sample_Type",
      SAMPLE_AMOUNT =	"Sample_Amount",
      SAMPLE_AMOUNT_UNIT = "Sample_Amount_Unit",
      ISTD_VOL ="ISTD_Mixture_Volume_[uL]",
      "BATCH_ID",
      "VALID_ANALYSIS"
    ) |>
    dplyr::group_by(.data$BATCH_ID) |>
    dplyr::mutate(BATCH_NO = dplyr::cur_group_id(),
                  ANALYSIS_ID = stringr::str_squish(as.character(.data$ANALYSIS_ID)),
                  ANALYSIS_ID = stringr::str_remove(.data$ANALYSIS_ID, "\\.mzML|\\.d"),
                  QC_TYPE = if_else(QC_TYPE == "Sample" | is.na(QC_TYPE), "SPL", QC_TYPE)
                  ) |>
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))


  d_annot$annot_features <- readxl::read_excel(filename, sheet = "Transition_Name_Annot", trim_ws = TRUE)|>
    dplyr::mutate(
      FEATURE_NAME = stringr::str_squish(.data$Transition_Name),
      NORM_ISTD_FEATURE_NAME	= stringr::str_squish(.data$Transition_Name_ISTD),
      QUANT_ISTD_FEATURE_NAME = stringr::str_squish(.data$Transition_Name_ISTD),
      isISTD = (.data$FEATURE_NAME == .data$NORM_ISTD_FEATURE_NAME),
      FEATURE_RESPONSE_FACTOR	= 1,
      isQUANTIFIER = dplyr::recode(tolower(.data$Quantifier), "yes" = TRUE, "no"=FALSE),
      isINTEGRATED = TRUE,
      REMARKS = NA_character_) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish)) %>%
    dplyr::select(dplyr::any_of(c("FEATURE_ID",
      "FEATURE_NAME",
      "isISTD",
      "NORM_ISTD_FEATURE_NAME",
      "QUANT_ISTD_FEATURE_NAME",
      "FEATURE_RESPONSE_FACTOR",
      "isQUANTIFIER",
      "isINTEGRATED",
      "REMARKS")))

  #ToDo: Merged cell in template
  annot_istd <- readxl::read_excel(filename,
                           sheet = "ISTD_Annot",
                           skip = 2,
                           trim_ws = TRUE,
                           .name_repair = ~ ifelse(nzchar(.x), .x, LETTERS[seq_along(.x)]))
  names(annot_istd)[1] <- "Transition_Name_ISTD"

  d_annot$annot_istd <- annot_istd |>
    dplyr::mutate(ISTD_COMPOUND_NAME = NA_character_) |>
    dplyr::select(
      QUANT_ISTD_FEATURE_NAME = "Transition_Name_ISTD",
      ISTD_CONC_nM = "ISTD_Conc_[nM]") %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  d_annot$annot_responsecurves <- readxl::read_excel(filename, sheet = "Dilution_Annot") |>
    dplyr::select(
      ANALYSIS_ID = "Sample_Name",
      RQC_SERIES_ID = "Dilution_Batch_Name",
      RELATIVE_SAMPLE_AMOUNT = "Relative_Sample_Amount_[%]",
      INJECTION_VOL = "Injection_Volume_[uL]") |>
    dplyr::mutate(
      ANALYSIS_ID = stringr::str_remove(.data$ANALYSIS_ID, "\\.mzML|\\.d"),
      ANALYSIS_ID = stringr::str_squish(as.character(.data$ANALYSIS_ID)),
      RQC_SERIES_ID = stringr::str_squish(as.character(.data$RQC_SERIES_ID)),
      RELATIVE_SAMPLE_AMOUNT = .data$RELATIVE_SAMPLE_AMOUNT/100) %>%
      dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  return(d_annot)
}
