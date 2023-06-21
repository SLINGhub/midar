#' Imports metadata provided by an MSOrganizer EXCEL template
#'
#' @param filename File path of the MSOrganizer EXCEL template (*.xlm)
#' @param trim_ws Trim all white spaces and double spaces
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom tidyselect vars_select_helpers
#' @importFrom dplyr select mutate filter group_by row_number
#' @importFrom stringr regex
#' @return A list of tibbles with different metadata
#' @export
import_msorganizer_xlm <- function(filename, trim_ws = TRUE){
  d_annot <- list()

  # ANALYSIS/SAMPLE annotation
  # -----------------
  # ToDo: Make note if feature names are not original

  d_temp_analyses <- readxl::read_excel(filename, sheet = "Sample_Annot")
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "ANALYSIS_ID", init_value = NA_character_, make_caps = TRUE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "VALID_ANALYSIS", init_value = TRUE, make_caps = TRUE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "REPLICATE", init_value = 1L, make_caps = TRUE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "SPECIMEN", init_value = NA_character_, make_caps = TRUE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "PANEL_ID", init_value = NA_character_, make_caps = TRUE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "SAMPLE_ID", init_value = NA_character_, make_caps = TRUE)
  d_temp_analyses <- d_temp_analyses |> add_missing_column(col_name = "REMARKS", init_value = NA_character_, make_caps = TRUE)

  # NOTE: If ANALYSIS_ID is defined, then it will overwrite the DATAFILE_NAME defined in the raw data files
  # Todo: if user-defined ANALYSIS_ID names (=remapping if IDs) are provided, then it should be reported somewhere,  possible source of user-error!

  d_annot$annot_analyses <- d_temp_analyses |>
    dplyr::mutate(
      VALID_ANALYSIS = TRUE,
      BATCH_ID = as.character(.data$BATCH_ID),
      RUN_ID_ANNOT = dplyr::row_number()) |>
    dplyr::select(
      "RUN_ID_ANNOT",
      "ANALYSIS_ID",
      DATAFILE_NAME = "Sample_Name",
      QC_TYPE = "Sample_Type",
      SAMPLE_AMOUNT =	"Sample_Amount",
      SAMPLE_AMOUNT_UNIT = "Sample_Amount_Unit",
      ISTD_VOL ="ISTD_Mixture_Volume_[uL]",
      "BATCH_ID",
      "VALID_ANALYSIS",
      "SPECIMEN",
      "SAMPLE_ID",
      "REMARKS"
    ) |>
    dplyr::mutate(BATCH_NO = dplyr::cur_group_id(), .by = c("BATCH_ID")) |>
    dplyr::mutate(
      DATAFILE_NAME = stringr::str_squish(as.character(.data$DATAFILE_NAME)),
      DATAFILE_NAME = stringr::str_remove(.data$DATAFILE_NAME, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)) ,
      ANALYSIS_ID = stringr::str_remove(.data$ANALYSIS_ID, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)) ,
      ANALYSIS_ID = stringr::str_squish(as.character(.data$ANALYSIS_ID)),
      ANALYSIS_ID = if_else(is.na(.data$ANALYSIS_ID), .data$DATAFILE_NAME, .data$ANALYSIS_ID),
      QC_TYPE = if_else(.data$QC_TYPE == "Sample" | is.na(.data$QC_TYPE), "SPL", .data$QC_TYPE))|>
    dplyr::ungroup() %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  # FEATURE annotation
  # -----------------
  # ToDo: Make note if feature names are not original
  d_temp_features <- readxl::read_excel(filename, sheet = "Transition_Name_Annot", trim_ws = TRUE)

  d_temp_features <- d_temp_features |> add_missing_column(col_name = "QUANTIFIER", init_value = TRUE, make_caps = TRUE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "VALID_INTEGRATION", init_value = TRUE, make_caps = TRUE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "RESPONSE_FACTOR", init_value = 1, make_caps = TRUE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "SOURCE_FEATURE_NAME", init_value = NA_character_, make_caps = TRUE)
  d_temp_features <- d_temp_features |> add_missing_column(col_name = "REMARKS", init_value = NA_character_, make_caps = TRUE)

  # NOTE: If FEATURE_NAME is defined, then it will overwrite the feature name defined in the raw data files
  # Todo: if user-defined feature names are provided, then it should be reported somewhere,  possible source of user-error!

  d_annot$annot_features <- d_temp_features |>
    dplyr::mutate(
      FEATURE_NAME = stringr::str_squish(.data$Transition_Name),
      SOURCE_FEATURE_NAME = stringr::str_squish(.data$SOURCE_FEATURE_NAME),
      SOURCE_FEATURE_NAME = if_else(is.na(.data$SOURCE_FEATURE_NAME), .data$FEATURE_NAME, .data$SOURCE_FEATURE_NAME),
      NORM_ISTD_FEATURE_NAME	= stringr::str_squish(.data$Transition_Name_ISTD),
      QUANT_ISTD_FEATURE_NAME = stringr::str_squish(.data$Transition_Name_ISTD),
      isISTD = (.data$FEATURE_NAME == .data$NORM_ISTD_FEATURE_NAME),
      QUANTIFIER = if_else(tolower(.data$QUANTIFIER) %in% c("yes","true"), TRUE, FALSE),
      REMARKS = NA_character_) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish)) %>%
    dplyr::select(dplyr::any_of(c(
      "SOURCE_FEATURE_NAME",
      "FEATURE_ID",
      "FEATURE_NAME",
      "isISTD",
      "NORM_ISTD_FEATURE_NAME",
      "QUANT_ISTD_FEATURE_NAME",
      FEATURE_RESPONSE_FACTOR = "RESPONSE_FACTOR",
      isQUANTIFIER = "QUANTIFIER",
      "VALID_INTEGRATION",
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
      DATAFILE_NAME = "Sample_Name",
      RQC_SERIES_ID = "Dilution_Batch_Name",
      RELATIVE_SAMPLE_AMOUNT = "Relative_Sample_Amount_[%]",
      INJECTION_VOL = "Injection_Volume_[uL]") |>
    dplyr::mutate(
      DATAFILE_NAME = stringr::str_remove(.data$DATAFILE_NAME, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE)),
      DATAFILE_NAME = stringr::str_squish(as.character(.data$DATAFILE_NAME)),
      RQC_SERIES_ID = stringr::str_squish(as.character(.data$RQC_SERIES_ID)),
      RELATIVE_SAMPLE_AMOUNT = .data$RELATIVE_SAMPLE_AMOUNT/100) %>%
      dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  return(d_annot)
}
