pkg.env <- new.env()

# Data structure templates
pkg.env$dataset_templates <- list(
  dataset_orig_template = tibble::tibble(
    "ANALYSIS_ID" = character(),
    "FEATURE_ID"= character(),
    "Intensity" = numeric(),
    "normIntensity" = numeric(),
    "Concentration" = numeric()
  ),

  annot_analyses_template = tibble::tibble(
    "ANALYSIS_ID" = character(),
    "DATAFILE_NAME"= character(),
    "SAMPLE_ID" = character(),
    "QC_TYPE" = factor(),
    "REPLICATE" = character(),
    "SAMPLE_NAME" = character(),
    "PANEL_ID" = character(),
    "BATCH_ID" = character(),
    "SPECIMEN" = character(),
    "SAMPLE_AMOUNT" = numeric(),
    "SAMPLE_AMOUNT_UNIT" = character(),
    "ISTD_VOL" = numeric(),
    "VALID_ANALYSIS" = logical(),
    "REMARKS" = character()
  ),
  annot_features_template = tibble::tibble(
    "FEATURE_ID" = character(),
    "FEATURE_NAME" = character(),
    "isISTD" = logical(),
    "NORM_ISTD_FEATURE_NAME" = character(),
    "QUANT_ISTD_FEATURE_NAME" = character(),
    "FEATURE_RESPONSE_FACTOR" = numeric(),
    "isQUANTIFIER" = logical(),
    "isINTEGRATED" = logical(),
    "REMARKS" = character()
  ),
  "annot_istd_template" = tibble::tibble(
    "ISTD_COMPOUND_NAME" = character(),
    "QUANT_ISTD_FEATURE_NAME" = character(),
    "ISTD_CONC_nM" = numeric()
  ),
  annot_responsecurves_template = tibble::tibble(
    "ANALYSIS_ID" = character(),
    "RQC_SERIES_ID" = character(),
    "RELATIVE_SAMPLE_AMOUNT" = numeric(),
    "INJECTION_VOL" = numeric()
  ),
  annot_batch_info_template = tibble::tibble(
    "BATCH_ID" = character(),
    "BATCH_NO" = numeric(),
    "id_batch_start" = numeric(),
    "id_batch_end" = numeric()
  )
)

pkg.env$qc_type_annotation <- list(
  qc_type_levels = c("SBLK", "TBLK", "UBLK", "PQC", "TQC", "BQC", "RQC", "EQC", "NIST",
                     "LTR", "PBLK", "SPL", "SST", "MBLK"),

  qc_type_col = c(
    "SBLK" = "#1854f9",
    "TBLK" = "#db0202",
    "UBLK" = "#de21de",
    "BQC" = "#db0202",
    "TQC" = "#1854f9",
    "PQC" = "#f99f18",
    "RQC" = "#96a4ff",
    "EQC" = "#513c3c",
    "NIST" = "#002e6b",
    "LTR" = "#880391",
    "PBLK" = "#08c105",
    "SPL" = "#aaaeaf",
    "SST" = "#bafc03",
    "MBLK" = "black"),

  qc_type_fillcol = c(
    "SBLK" = "#f891ff",
    "TBLK" = "#fffb03",
    "UBLK" = "#c1bd04",
    "BQC" = "#db0202",
    "TQC" = "#1854f9",
    "PQC" = "#f99f18",
    "RQC" = "#688ff9",
    "EQC" = "NA",
    "NIST" = "#cce2ff",
    "LTR" = "#880391",
    "PBLK" = "#08c105",
    "SPL" = "NA",
    "SST" = "#aaaeaf",
    "MBLK" = "black"),

  qc_type_shape = c(
    "SBLK" = 23,
    "TBLK" = 23,
    "UBLK" = 23,
    "BQC" = 16,
    "TQC" = 25,
    "PQC" = 25,
    "RQC" = 6,
    "EQC" = 24,
    "NIST" = 23,
    "LTR" = 23,
    "PBLK" = 23,
    "SPL" = 21,
    "SST" = 10,
    "MBLK" = 10)
)

# allow S4 to see the class tbl_df
setOldClass(c("tbl_df", "tbl", "data.frame"))

#' S4 Class Representing the MIDAR Dataset
#'
#' @description
#' Core class to store raw and processed data, annotations and metadata.
#'
#' @details
#' Provides also functions to validate data integrity, to perfom normalization and quantitation based on internal standards and to calculate the QC parameters of the experiment
#'
#' @docType class
#'
#' @slot dataset_orig Original imported analysis data. Required fields:
#' @slot dataset Processed analysis data. Required fields:
#' @slot dataset_QC_filtered Processed analysis data. Required fields:
#' @slot annot_analyses Annotation of analyses/runs
#' @slot annot_features Annotation of measured features.
#' @slot annot_istd Annotation of Internal Standard concentrations.
#' @slot annot_responsecurves Annotation of  Response curves (RQC). Required fields
#' @slot annot_studysamples Annotation of study samples. Required fields:
#' @slot annot_batch_info Annotation of batches. Required fields:
#' @slot d_QC QC information for each measured feature
#' @slot status_processing Status within the data processing workflow
#'
#' @export
#'
#' @importFrom utils tail
#' @importFrom tibble tibble

setClass("MidarExperiment",
         slots = c(
           dataset_orig = "tbl_df",
           dataset = "tbl_df",
           dataset_QC_filtered = "tbl_df",
           annot_analyses = "tbl_df",
           annot_features = "tbl_df",
           annot_istd = "tbl_df",
           annot_responsecurves= "tbl_df",
           annot_studysamples = "tbl_df",
           annot_batch_info = "tbl_df",
           d_QC = "tbl_df",
           status_processing = "character"
         ),
         prototype = list(
           dataset_orig = pkg.env$dataset_templates$annot_analyses_template,
           dataset = pkg.env$dataset_templates$annot_analyses_template,
           dataset_QC_filtered = pkg.env$dataset_templates$annot_analyses_template,
           annot_analyses = pkg.env$dataset_templates$annot_analyses_template,
           annot_features = pkg.env$dataset_templates$annot_features_template,
           annot_istd = pkg.env$dataset_templates$annot_istd_template,
           annot_responsecurves= pkg.env$dataset_templates$annot_responsecurves_template,
           annot_studysamples = tibble::tibble(),
           annot_batch_info = tibble::tibble(),
           d_QC = tibble::tibble(),
           status_processing = "No Data"
         )
)

#' Constructor for the MidarExperiment object.
#' @importFrom methods new

#' @return `MidarExperiment` object
#' @export
MidarExperiment <- function() {
  methods::new("MidarExperiment")
}


check_integrity <-  function(object) {
  #browser()
  if (nrow(object@dataset_orig) > 0 & nrow(object@annot_analyses) > 0) {
    d_xy <- length(setdiff(object@dataset_orig$ANALYSIS_ID %>% unique(), object@annot_analyses$ANALYSIS_ID))
    d_yx <- length(setdiff(object@annot_analyses$ANALYSIS_ID,object@dataset_orig$ANALYSIS_ID %>% unique()))
    if(d_xy > 0){
      print(glue::glue("{d_xy} of {object@dataset_orig$ANALYSIS_ID %>% unique() %>% length()} measurements have no matching sample metadata"))
      if (d_xy < 100) print(paste0(setdiff(object@dataset_orig$ANALYSIS_ID %>% unique(), object@annot_analyses$ANALYSIS_ID), collapse = ", "))
      else {
        print("too many to display")
        #rint(paste0(setdiff(object@dataset_orig$ANALYSIS_ID %>% unique(), object@annot_analyses$ANALYSIS_ID[1:100]), collapse = ", ")
        }
      FALSE
      } else if(d_yx > 0) {
      warning(glue::glue("{d_yx} of {object@annot_analyses$ANALYSIS_ID %>% length()} sample metadata are not found in the measurement data"))
      object@status_processing <- "Data and Sample Metadata loaded"
      TRUE
      } else {
      object@status_processing <- "Data and Sample Metadata loaded"
      TRUE
    }

  } else {

    TRUE
  }
}
#' @importFrom methods setValidity
#'
methods::setValidity("MidarExperiment", check_integrity)





#' loadMasshunterCSV
#'
#' @param data MidarExperiment object
#' @param filename file name of the MH CSV file
#'
#' @importFrom methods validObject
#'
#' @return MidarExperiment object
#' @export
#'
#'
loadMasshunterCSV <- function(data, filename) {
  data@dataset_orig <- read_MassHunterCSV(filename, silent = FALSE)
  data@dataset_orig <- data@dataset_orig %>% dplyr::rename(Intensity = "Area")
  stopifnot(methods::validObject(data))
  data
}

#' #' load_MRMkit_csv
#' #'
#' #' @param data MidarExperiment object
#' #' @param filename file name of the MH CSV file
#' #'
#' #' @importFrom methods validObject
#' #'
#' #' @return MidarExperiment object
#' #' @export
#' #'
#' #'
#' load_MRMkit_csv <- function(data, filename) {
#'   data@dataset_orig <- read_MRMkitCSV(filename, silent = FALSE)
#'   data@dataset <- data@dataset_orig  %>%
#'     dplyr::left_join(data@annot_analyses, by = c("ANALYSIS_ID"="ANALYSIS_ID"))
#'   stopifnot(methods::validObject(data))
#'   data
#' }


#' Import metadata from the MSOrganizer template (.XLM)
#' @param data MidarExperiment object
#' @param filename file name and path


#' @return MidarExperiment object
#' @export
#'

loadMSOrganizerXLM <- function(data, filename) {
  d_annot <- import_MSOrganizerXLM(filename)

  data@annot_analyses <- data@dataset_orig %>%
    dplyr::select("ANALYSIS_ID") %>%
    dplyr::distinct() %>%
    dplyr::right_join(d_annot$annot_analyses, by = c("ANALYSIS_ID" = "ANALYSIS_ID"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_analyses_template)

  data@annot_features <- data@dataset_orig %>%
    dplyr::select("FEATURE_NAME") %>%
    dplyr::distinct() %>%
    dplyr::right_join(d_annot$annot_features, by = c("FEATURE_NAME"="FEATURE_NAME"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_features_template)

  data@annot_istd <- data@annot_features %>%
    dplyr::select("QUANT_ISTD_FEATURE_NAME") %>%
    dplyr::distinct() %>%
    dplyr::left_join(d_annot$annot_istd, by = c("QUANT_ISTD_FEATURE_NAME"="QUANT_ISTD_FEATURE_NAME"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_istd_template)

  data@annot_responsecurves <- data@annot_analyses %>%
    dplyr::filter("QC_TYPE" == "RQC") %>%
    dplyr::select("ANALYSIS_ID") %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::right_join(d_annot$annot_responsecurves, by = c("ANALYSIS_ID"="ANALYSIS_ID"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_responsecurves_template)


  data@annot_batch_info <- data@annot_analyses %>%
    dplyr::group_by(.data$BATCH_NO) %>%
    dplyr::summarise(
      BATCH_ID = .data$BATCH_ID[1],
      id_batch_start = dplyr::first(.data$RUN_ID_ANNOT),
      id_batch_end = dplyr::last(.data$RUN_ID_ANNOT)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(.data$id_batch_start)%>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_batch_info_template)

  data@dataset_orig <- data@dataset_orig %>%
    dplyr::bind_rows(pkg.env$dataset_templates$dataset_orig_template)

  data@dataset <- data@dataset_orig  %>% dplyr::select(-dplyr::any_of(c("FEATURE_ID"))) %>%
    dplyr::left_join(data@annot_analyses  %>% dplyr::select("ANALYSIS_ID", "QC_TYPE", "BATCH_ID"), by = c("ANALYSIS_ID")) %>%
    dplyr::right_join(d_annot$annot_features %>% dplyr::select(dplyr::any_of(c("FEATURE_NAME", "NORM_ISTD_FEATURE_NAME", "isISTD", "FEATURE_ID", "isQUANTIFIER"))), by = c("FEATURE_NAME"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$dataset_orig_template)
  stopifnot(methods::validObject(data))
  data
}



