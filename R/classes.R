pkg.env <- new.env()

# Data structure templates
pkg.env$dataset_templates <- list(
  dataset_orig_template = tibble::tibble(
    "DATAFILE_NAME" = character(),
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
    "REPLICATE" = integer(),
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
    "isQUANTIFIER" = logical(),
    "VALID_INTEGRATION" = logical(),
    "FEATURE_RESPONSE_FACTOR" = numeric(),
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
#' @slot analysis_type Analysis type, one of "lipidomics", "metabolomics", "externalcalib", "others"
#' @slot dataset_orig Original imported analysis data. Required fields:
#' @slot dataset Processed analysis data. Required fields:
#' @slot dataset_QC_filtered Processed analysis data. Required fields:
#' @slot annot_analyses Annotation of analyses/runs
#' @slot annot_features Annotation of measured features.
#' @slot annot_istd Annotation of Internal Standard concentrations.
#' @slot annot_responsecurves Annotation of  Response curves (RQC). Required fields
#' @slot annot_studysamples Annotation of study samples. Required fields:
#' @slot annot_batch_info Annotation of batches. Required fields:
#' @slot metrics_qc QC information for each measured feature
#' @slot status_processing Status within the data processing workflow
#' @slot is_istd_normalized Flag if data has been ISTD normalized
#' @slot is_quantitated Flag if data has been quantitated using ISTD and sample amount
#' @slot is_drift_corrected Flag if data has been drift corrected
#' @slot is_batch_corrected Flag if data has been batch corrected
#' @slot is_isotope_corr Flag if one or more features have been isotope corrected
#' @export
#'
#' @importFrom utils tail
#' @importFrom tibble tibble

setClass("MidarExperiment",
         slots = c(
           analysis_type = "character",
           dataset_orig = "tbl_df",
           dataset = "tbl_df",
           dataset_QC_filtered = "tbl_df",
           annot_analyses = "tbl_df",
           annot_features = "tbl_df",
           annot_istd = "tbl_df",
           annot_responsecurves= "tbl_df",
           annot_studysamples = "tbl_df",
           annot_batch_info = "tbl_df",
           metrics_qc = "tbl_df",
           status_processing = "character",
           is_istd_normalized = "logical",
           is_quantitated = "logical",
           is_drift_corrected = "logical",
           is_batch_corrected = "logical",
           is_isotope_corr = "logical"
         ),
         prototype = list(
           analysis_type = "",
           dataset_orig = pkg.env$dataset_templates$annot_analyses_template,
           dataset = pkg.env$dataset_templates$annot_analyses_template,
           dataset_QC_filtered = pkg.env$dataset_templates$annot_analyses_template,
           annot_analyses = pkg.env$dataset_templates$annot_analyses_template,
           annot_features = pkg.env$dataset_templates$annot_features_template,
           annot_istd = pkg.env$dataset_templates$annot_istd_template,
           annot_responsecurves= pkg.env$dataset_templates$annot_responsecurves_template,
           annot_studysamples = tibble::tibble(),
           annot_batch_info = tibble::tibble(),
           metrics_qc = tibble::tibble(),
           status_processing = "No Data",
           is_istd_normalized = FALSE,
           is_quantitated = FALSE,
           is_drift_corrected = FALSE,
           is_batch_corrected = FALSE,
           is_isotope_corr = FALSE
         )
)

#' Constructor for the MidarExperiment object.
#' @importFrom methods new
#' @param analysis_type Analysis type, one of "lipidomics", "metabolomics", "externalcalib", "others"
#' @return `MidarExperiment` object
#' @export
MidarExperiment <- function(analysis_type = "") {
  methods::new("MidarExperiment", analysis_type = analysis_type)
}

#' Set analysis type
#' @description
#' Set the analysis type, i.e. "lipidomics", "metabolomics, "quantitative"
#' @param x MidarExperiment object
#' @return A character string
#' @export
setGeneric("analysis_type", function(x) standardGeneric("analysis_type"))


#' Get analysis type
#' @description
#' Get the analysis type defined for the MidarExperiment object
#' @param x MidarExperiment object
#' @param value Analysis type, one of "lipidomics", "metabolomics, "quantitative"
#' @export
setGeneric("analysis_type<-", function(x, value) standardGeneric("analysis_type<-"))

#' Get `analysis_type`
#' @param x MidarExperiment object
#' @return A character string
setMethod("analysis_type", "MidarExperiment", function(x) x@analysis_type)

#' Set `analysis_type`
#' @param x MidarExperiment object
#' @param value Analysis type, one of "lipidomics", "metabolomics, "quantitative"
setMethod("analysis_type<-", "MidarExperiment", function(x, value) {
  x@analysis_type <- value
  x
})


check_integrity <-  function(object, excl_unannotated_analyses) {
  #browser()
  if (nrow(object@dataset_orig) > 0 & nrow(object@annot_analyses) > 0) {
    d_xy <- length(setdiff(object@dataset_orig$DATAFILE_NAME %>% unique(), object@annot_analyses$DATAFILE_NAME))
    d_yx <- length(setdiff(object@annot_analyses$DATAFILE_NAME,object@dataset_orig$DATAFILE_NAME %>% unique()))
    if(d_xy > 0){
      if(d_xy == length(object@dataset_orig$DATAFILE_NAME %>% unique())) stop("Error: None of the measurements/samples have matching metadata . Please check data and metadata files.")
        if(!excl_unannotated_analyses){
          stop(glue::glue("Error: {d_xy} of {object@dataset_orig$DATAFILE_NAME %>% unique() %>% length()} measurements have no matching metadata."))
        if (d_xy < 50)
          writeLines(glue::glue("No metadata present for: {paste0(setdiff(object@dataset_orig$DATAFILE_NAME %>% unique(), object@annot_analyses$DATAFILE_NAME), collapse = ", ")} measurements."))
        else
          print("No metadata present for: Too many (> 50) to display")
        } else {
            writeLines(crayon::yellow(glue::glue("! Note: {d_xy} of {object@dataset_orig$DATAFILE_NAME %>% unique() %>% length()} measurements without matching metadata were excluded.")))
        }
      } else if(d_yx > 0) {
        stop(glue::glue("{d_yx} of {object@annot_analyses$DATAFILE_NAME %>% length()} sample metadata are not found in the measurement data."))
      } else {
      object@status_processing <- "DataMetadataLoaded"
      TRUE
      }
  }
}

##' @importFrom methods setValidity
##'
#methods::setValidity("MidarExperiment")


get_status_flag <- function(x) if_else(x, crayon::green$bold('\u2713'), crayon::red$bold('\u2717'))




#' Getter for specific slots of an MidarExperiments object
#'
#' $ syntax can be used to as a shortcut for getting specific variables and results from a MidarExperiment object
#' @return Value with a variable or a tibble
#' @param x MidarExperiment object
#' @param name MidarExperiment slot
#' @examples
#' mexp = MidarExperiment()
#' mexp$analysis_type
#' mexp$annot_analyses
#' @importFrom methods slot
#' @export
setMethod(f = "$",
  signature = c("MidarExperiment"),
  definition = function(x,name) {

  # check for other struct slots
  valid=c('analysis_type','dataset','annot_analyses', 'annot_features', 'annot_istd', 'metrics_qc', 'annot_batches', 'dataset_QC_filtered', 'is_istd_normalized')
  if (!name %in% valid) stop('"', name, '" is not valid for this object: ', class(x)[1])
  methods::slot(x,name)
  }
)



setMethod("show", "MidarExperiment", function(object) {
  cat("\n", is(object)[[1]], "\n",
      "\n",
      "  Data: ", "\n",
      "  \u2022 Samples: ", length(unique(object@dataset$ANALYSIS_ID)), "\n",
      "  \u2022 Features:  ", length(unique(object@dataset$FEATURE_NAME)), "\n",
      "\n",
      "  Metadata: ", "\n",
      "  \u2022 Sample annotation: ", get_status_flag(nrow(object@annot_analyses) > 0), "\n",
      "  \u2022 Feature annotation: ", get_status_flag(nrow(object@annot_features) > 0), "\n",
      "  \u2022 Internal standard annotation: ", get_status_flag(nrow(object@annot_istd) > 0), "\n",
      "  \u2022 Response curves annotation: ", get_status_flag(nrow(object@annot_responsecurves) > 0), "\n",
      "  \u2022 Study samples annotation: ", get_status_flag(nrow(object@annot_studysamples) > 0), "\n",
      "\n",
      "  Processing status: ",object@status_processing, "\n",
      "\n",
      "  Processing: ", "\n",
      "  \u2022 ISTD normalized: ", get_status_flag(object@is_istd_normalized), "\n",
      "  \u2022 ISTD quantitated: ", get_status_flag(object@is_quantitated), "\n",
      "  \u2022 Drift corrected: ", get_status_flag(object@is_drift_corrected), "\n",
      "  \u2022 Batch corrected: ", get_status_flag(object@is_batch_corrected), "\n",
      "  \u2022 Interference (isotope) corrected: ", get_status_flag(object@is_isotope_corr), "\n"
  )
})


#' @title Reads an Agilent MassHunter Quant CSV file
#' @description
#' Imports a .csv file with `Agilent MassHunter Quantitative Analysis` results.
#' Samples should be in rows, features/compounds in columns and must contain either peak areas, peak heights or intensities.
#' Additional columns, such as RT (rentention time), FWHM, PrecursorMZ, and CE will be imported and available from the `MidarExperiment` object for downstream analyses.
#'
#' When more than one file is provided, all files are imported and merged into one raw dataset. This can be useful, e.g. when importing datasets that are pre-processing in blocks resulting in different files.
#' Each Datafile/Feature pair needs to be unique within and across data files, presense of replicates return an error.
#'
#' @param data MidarExperiment object
#' @param file_dir_names One or more file names with path or folder path. When a folder name is given, all *.csv files in this folder will be read.
#'
#' @importFrom methods validObject
#' @importFrom fs is_dir path_tidy file_exists dir_ls
#'
#' @return MidarExperiment object
#' @examples
#' csvfile <- system.file("extdata", "Example_MHQuant_1.csv", package = "midar")
#' mexp <- MidarExperiment()
#' mexp <- read_masshunter_csv(mexp, csvfile)
#' mexp

#' @export

read_masshunter_csv <- function(data, file_dir_names) {
  #data@dataset_orig <- import_masshunter_csv(filenames, silent = FALSE)

  if(!fs::is_dir(file_dir_names))
    file_paths <- fs::path_tidy(file_dir_names)
  else
    file_paths <- fs::dir_ls(file_dir_names, glob = "*.csv")

  if(!all(fs::file_exists(file_paths))) stop("One or more given files do not exist. Please check file paths.")
  if(any(duplicated(file_paths))) stop("One or more given files are replicated. Please check file paths.")

  d_temp <- file_paths  |>
    map_dfr(import_masshunter_csv, .id = "DATA_SOURCE")

  # Test if ANALYSIS_IDs (=DATAFILE_NAME), FEATURE_NAMEs, and values are replicated
  if(nrow(d_temp) > nrow(d_temp |> distinct(.data$DATAFILE_NAME, .data$SOURCE_FEATURE_NAME, .keep_all = FALSE))){
    has_duplicated_id <- TRUE
    if(nrow(d_temp) > nrow(d_temp |> distinct(.data$DATAFILE_NAME, .data$SOURCE_FEATURE_NAME, .keep_all = TRUE)))
      has_duplicated_id_values <- TRUE
    else
      has_duplicated_id_values <- FALSE
  } else {
    has_duplicated_id <- FALSE
  }

  if(has_duplicated_id)
    if(has_duplicated_id_values)
      stop(glue::glue("Dataset(s) contains replicated reportings (analysis and feature pairs) with identical values. Please check dataset(s)."))
    else
      stop(glue::glue("Dataset(s) contains replicated reportings (analysis and feature pairs) with different values.Please check dataset(s)."))



  data@dataset_orig <- d_temp

  data@dataset_orig <- data@dataset_orig %>% dplyr::rename(Intensity = "Area")
  check_integrity(data, excl_unannotated_analyses = FALSE)
  #stopifnot(methods::validObject(data))
  #stopifnot(methods::validObject(data))
  data@status_processing <- "Raw Data"
  data
}


# ##########
# TODO: This below has to be changed, mapping of metadata to data should be a distinct function
# ##########
#' @title Import metadata from the MSOrganizer template (.XLM)
#' @param data MidarExperiment object
#' @param filename file name and path
#' @param excl_unannotated_analyses Exclude analyses (samples) that have no matching metadata
#' @return MidarExperiment object
#' @export
#'

read_msorganizer_xlm <- function(data, filename, excl_unannotated_analyses = FALSE) {
  d_annot <- import_msorganizer_xlm(filename)

  data@annot_analyses <- data@dataset_orig %>%
    dplyr::select("DATAFILE_NAME") %>%
    dplyr::distinct() %>%
    dplyr::right_join(d_annot$annot_analyses, by = c("DATAFILE_NAME" = "DATAFILE_NAME"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_analyses_template)

  data@annot_features <- data@dataset_orig %>%
    dplyr::select("SOURCE_FEATURE_NAME") %>%
    dplyr::distinct() %>%
    dplyr::right_join(d_annot$annot_features, by = c("SOURCE_FEATURE_NAME"="SOURCE_FEATURE_NAME"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_features_template)

  data@annot_istd <- data@annot_features %>%
    dplyr::select("QUANT_ISTD_FEATURE_NAME") %>%
    dplyr::distinct() %>%
    dplyr::left_join(d_annot$annot_istd, by = c("QUANT_ISTD_FEATURE_NAME"="QUANT_ISTD_FEATURE_NAME"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$annot_istd_template)


  data@annot_responsecurves <- data@annot_analyses %>%
    dplyr::filter(.data$QC_TYPE == "RQC") %>%
    dplyr::select("ANALYSIS_ID", "DATAFILE_NAME") %>%
    dplyr::distinct() %>%
    dplyr::ungroup() %>%
    dplyr::right_join(d_annot$annot_responsecurves, by = c("DATAFILE_NAME"="DATAFILE_NAME"), keep = FALSE) %>%
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
    dplyr::inner_join(data@annot_analyses  %>% dplyr::select("ANALYSIS_ID", "DATAFILE_NAME", "QC_TYPE", "SAMPLE_ID", "REPLICATE", "VALID_ANALYSIS", "BATCH_ID"), by = c("DATAFILE_NAME")) %>%
    dplyr::inner_join(d_annot$annot_features %>% filter(.data$VALID_INTEGRATION) |>  dplyr::select(dplyr::any_of(c("FEATURE_NAME", "FEATURE_NAME", "NORM_ISTD_FEATURE_NAME", "isISTD", "SOURCE_FEATURE_NAME", "FEATURE_ID", "isQUANTIFIER", "VALID_INTEGRATION", "FEATURE_RESPONSE_FACTOR", "INTERFERING_FEATURE", "INTERFERANCE_PROPORTION"))),
                      by = c("SOURCE_FEATURE_NAME"), keep = FALSE) %>%
    dplyr::bind_rows(pkg.env$dataset_templates$dataset_orig_template) |>
    mutate(Corrected_Interference = FALSE)
  #stopifnot(methods::validObject(data, excl_nonannotated_analyses))
  check_integrity(data, excl_unannotated_analyses = excl_unannotated_analyses)
  data@status_processing <- "Annotated Raw Data"


  writeLines(crayon::green(glue::glue("\u2713 Metadata successfully associated with {length(data@dataset$ANALYSIS_ID %>% unique())} samples and {length(data@dataset$FEATURE_NAME %>% unique())} features.")))
  data
}



