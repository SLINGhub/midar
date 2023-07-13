#' Read and convert an Agilent MassHunter Quant CSV result file
#'
#' @param filename File path of MassHunter Quant CSV file
#' @param silent Suppress messages
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom glue glue
#' @importFrom tidyselect vars_select_helpers
#' @return A tibble in the long format
#' @export

import_masshunter_csv <- function(filename, silent = FALSE) {

  if(!silent) print(glue::glue("Reading [{basename(filename)}] ..."))
  # if (shiny::isRunning())
  #   incProgress(1 / length(n_datafiles), detail = paste0(", basename(file)))
  #
  # Read Agilent MassHunter Quant Export file (CSV)
  suppressWarnings(suppressMessages(
    datWide <-
      readr::read_csv(
        file = filename,
        col_names = FALSE,
        na = c("#N/A", "NULL"),
        trim_ws = TRUE,
        col_types = readr::cols(.default = "c"),
        locale = readr::locale(encoding = 'ISO-8859-1'), num_threads = 4,progress = TRUE
      )))
  warnings_datWide <- readr::problems(datWide)



  # Remove text that is not required and remove dot chars that interfere later with the conversion wide to long
  # ToDo: Convert to tidyverse functions
  datWide[1, ] <-
    lapply(datWide[1, ], function(y)
      gsub(" Results", "", y))
  datWide[1, ] <-
    lapply(datWide[1, ], function(y)
      gsub(" Method", "", y))
  datWide[2, ] <- lapply(datWide[2, ], function(y)
    gsub("\\. ", "", y))
  datWide[2, ] <- lapply(datWide[2, ], function(y)
    gsub("\\.", "", y))
  datWide[2, ] <- lapply(datWide[2, ], function(y)
    gsub("/", "", y))


  datWide <- datWide |> dplyr::add_row(.after = 1)
  datWide[1, ] <- tibble::tibble(A = datWide[1,] |> unlist() |> dplyr::na_if("")) |>  tidyr::fill("A") |> unlist() |> as.list()

  datWide[1, ] <- replace(datWide[1,], stringr::str_detect(string = datWide[1,] , pattern = "AA"),"")
  datWide[2, ] <- replace(datWide[1,], !stringr::str_detect(string = datWide[1,] , pattern = "AA"),"")
  datWide[1, ] <- tibble::tibble(A = datWide[1,] |> unlist() |> dplyr::na_if("")) |>  tidyr::fill("A") |> unlist() |> as.list()

  # Concatenate rows
  datWide[1, ] <- paste(datWide[1, ], datWide[2, ], sep = " ") |> stringr::str_squish() |> as.list()
  datWide <- datWide[-2, ]


  # Fill in transition name in empty columns (different parameters of the same transition)
  # Dealing with exported "Qualifier" results: Adds the corresponding quantifier name (the first column of the group) with the Tag "QUAL" and the transition in front e.g. "Sph d18:0 [QUAL 302.3>266.2]"


  # Remove columns with no column names, as they are undefined (seems to be only Outlier Summary and Quantitation Message Summary)
  datWide <- datWide |> dplyr::select(-tidyselect::where(~ .x[2] == ""))

  # Concatenate rows containing parameters + transitions to the form parameter.sample and parameter.transition headers. This will later be converted with reshape()
  colnames(datWide) <- paste(datWide[2, ], datWide[1, ], sep = "\t")
  datWide <- datWide[-1:-2, ]


  # Rename some known column header names and remove columns that are not needed or not informative.

  new_colnames = c(
    "DATAFILE_NAME" = "Data File\tSample",
    "SampleName" = "Name\tSample",
    "AcqTimeStamp" = "AcqDate-Time\tSample",
    "SampleType" = "Type\tSample",
    "VialPosition" = "Pos\tSample",
    "InjVolume" = "Vol\tSample",
    "Dilution" = "Dil\tSample",
    "SampleGroup" = "Sample Group\tSample",
    "InstrumentName" = "Instrument Name\tSample",
    "InstrumentType" = "Instrument Type\tSample",
    "AcqMethodFile" = "Acq. Method File\tSample",
    "AcqMethodPath" = "Acq. Method Path\tSample",
    "DataPath" = "Data Path\tSample",
    "QuantitationMessage" = "Quantitation Message\tSample",
    "NameCompound" = "Name\tCompound"
  )

  datWide <- datWide %>% dplyr::rename(dplyr::any_of(new_colnames))

  if("QuantitationMessage" %in% names(datWide))  stop("Field 'Quantitation Message' currently not supported: Please re-export your data in MH without this field.")
  if("NameCompound" %in% names(datWide)) stop("Compound table format is currently not supported. Please re-export your data in MH with compounds as columns.")
  if(! "DATAFILE_NAME" %in% names(datWide)) stop("Unknown format or corrupt data file. Please try re-export from MH.")
  if (nrow(warnings_datWide)> 0) stop("Unknown format, or data file is corrupt. Please try re-export from MH.")
  # Remove ".Sample" from remaining sample description headers and remove known unused columns
  datWide <-
    datWide[, !(names(datWide) %in% c("NA\tSample", "Level\tSample", "Sample"))]
  names(datWide) <- sub("\\\tSample", "", names(datWide))

  keep.cols <- names(datWide) %in% c("")
  datWide <- datWide [!keep.cols]



  # Transform wide to long (tidy) table format
  # ------------------------------------------



  datWide <- datWide %>%
    dplyr::mutate(
      dplyr::across(.cols = dplyr::any_of(c("AcqTimeStamp")),
                    .fns = \(x) lubridate::parse_date_time(x, c("mdy_HM", "dmy_HM", "ymd_HM", "ydm_HM", "mdy_HM %p", "dmy_HM %p", "ymd_HM %p", "ydm_HM %p"))),
      dplyr::across(.cols = dplyr::any_of(c("SampleName")),
                    .fns = stringr::str_squish)
    )

  # obtain list with column names of the columns that define transition values (e.g. "RT Cer d16:0/18:0"). Delimuter is currently tab (\t)
  param_transition_names <-
    colnames(datWide[, -1:-tail(grep("\\\t", colnames(datWide), invert =  TRUE), 1)])


  # Obtain long table of all param-transition combinations, split param and compund name and then spread values of different param as columns
  datLong <- datWide %>%
    dplyr::mutate(RUN_ID = dplyr::row_number(), .before = 1) |>
    tidyr::pivot_longer(cols = tidyselect::all_of(param_transition_names), names_pattern = "(.*)\t(.*)$", names_to = c("Param", "SOURCE_FEATURE_NAME")) %>%
    tidyr::pivot_wider(names_from = "Param" ,values_from = "value")

  # Convert types of knows parameters and fields in the data set
  # ------------------------------------------------------------

  datLong <- datLong %>%
    dplyr::mutate(
      dplyr::across(.cols = dplyr::any_of(c("RT", "Area", "Height","FWHM","Width","SN","IntStart","IntEnd",
                                            "Symmetry","InjVolume", "Precursor Ion", "Product Ion", "Collision Energy")),
                    .fns = \(x) as.numeric(stringr::str_replace(x, ",", "."))),
      dplyr::across(.cols = dplyr::any_of(c("MI")),
                    .fns = as.logical),
      dplyr::across(.cols = dplyr::any_of(c("Ion Polarity")),
                    .fns = as.factor)
    )

  new_colnames <- c(DataName = "SampleName", PRECURSOR_MZ = "Precursor Ion", PRODUCT_MZ = "Product Ion",
                    CollisionEnergy = "Collision Energy", IonPolarity = "Ion Polarity")

  datLong <- datLong %>%
    dplyr::rename(dplyr::any_of(new_colnames)) %>%
    dplyr::mutate(DATAFILE_NAME = stringr::str_replace(.data$DATAFILE_NAME, "\\.d", ""), .before = "DataName") %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  if(!silent) {
    writeLines(crayon::green(glue::glue("\u2713 Imported {length(unique(datLong$DATAFILE_NAME))} samples with {length(unique(datLong$SOURCE_FEATURE_NAME))} transitions. \n")))
  }
  datLong
}


#' Reads a wide CSV file with Feature Intensities
#'
#' @param file File name and path of the MassHunter Quant CSV file
#' @param field Peak parameter (e.g. Area, RT)
#' @param silent Suppress messages
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom tidyselect vars_select_helpers
#' @return A tibble in the long format
#' @export

import_masshunter_csv_wide <- function(file, field, silent = FALSE) {
  sample_def_cols = c(
    "DATAFILE_NAME",
    "SampleName",
    "AcqTimeStamp",
    "SampleType",
    "VialPosition",
    "InjVolume",
    "SampleType"
  )
  d <- import_masshunter_csv(file, silent) %>%
    dplyr::select(tidyselect::any_of(sample_def_cols), "SOURCE_FEATURE_NAME", {{field}})

  d %>% tidyr::pivot_wider(names_from = "SOURCE_FEATURE_NAME", values_from = {{field}})
}



#' Reads a long CSV file with Feature Intensities
#'
#' @param file File name and path of a plain long-format CSV file
#' @param field Peak parameter (e.g. Area, RT)
#' @param silent Suppress messages
#'
#' @return A tibble in the long format
#' @export
read_long_table_csv <- function(file, field, silent = FALSE) {

  sample_def_cols = c(
    "DATAFILE_NAME",
    "SampleName",
    "AcqTimeStamp",
    "SampleType",
    "VialPosition",
    "InjVolume",
    "SampleType",
    "RUN_ID",
    "ANALYSIS_ID",
    "SOURCE_FEATURE_NAME",
    "Area",
    "RT",
    "FWHM",
    "SN"
  )

  d <- readr::read_csv(file, col_names = TRUE, trim_ws = TRUE, progress = TRUE) %>%
    dplyr::select(tidyselect::any_of(sample_def_cols), "SOURCE_FEATURE_NAME", Intensity = {{field}})

  d
}


#' Reads a wide CSV or XLSX sheet with Feature Values
#'
#' @param data MidarExperiment()
#' @param file File name and path of a plain wide-format CSV or XLS file
#' @param field Peak parameter (e.g. Area, RT)
#' @param sheet Sheet name
#' @param silent Suppress messages
#' @importFrom rlang :=
#'
#' @return A tibble in the long format
#' @export
#'
read_table_wide <- function(data, file, field, sheet = "", silent = FALSE) {

  if (!field %in% c("Intensity", "normIntensity", "Concentration", "RT", "FWHM")) stop("Field can only be one of: Intensity, normIntensity, Concentration or RT")

  var_field <- rlang::ensym(field)

  ext <- fs::path_ext(file)
  browser
  if(ext == "csv")
    d <- readr::read_csv(file, col_names = TRUE, trim_ws = TRUE, progress = FALSE, na = c("n/a", "N/A", "NA", "na", "ND", "N.D.", "n.d."), col_types = "cn")
  else if(ext == "xls" | ext == "xlsx"){
    if(sheet == "") stop("Please define sheet name via the `sheet` parameter")
    d <- readxl::read_excel(path = file, sheet = sheet, trim_ws = TRUE, progress = TRUE, na = c("n/a", "N/A"))
  }
  else
    stop("Invalid file format. Only CSV, XLS and XLSX supported.")

  d <- d |> dplyr::mutate(RUN_ID = row_number(),.before = 1)

  d <- d |> tidyr::pivot_longer(cols = -1:-2, names_to = "SOURCE_FEATURE_NAME" , values_to = field) |>
    dplyr::rename(ATAFILE_NAME = 2) |>
    mutate({{field}} := as.numeric(!!var_field))

  data@dataset_orig <- d
  data

}

#' Read and convert an Agilent MassHunter Quant CSV result file
#'
#' @param filename File path of MassHunter Quant CSV file
#' @param use_mrmkit_normdata use raw or MRMkit-normalized data
#' @param silent Suppress messages
#' @importFrom stats na.omit setNames
#' @importFrom utils tail
#' @importFrom glue glue
#' @importFrom readr read_csv
#' @importFrom tidyselect vars_select_helpers
#' @return A tibble in the long format
#' @export
read_MRMkit_raw_area_CSV<- function(filename, use_mrmkit_normdata = FALSE, silent = FALSE) {

  sep <- if_else(fs::path_ext(filename) == "csv", ",", "\t")

  d_mrmkit_raw <- readr::read_delim(filename, delim = sep, col_types = readr::cols(.default = "c"),
                                       col_names = TRUE, trim_ws = TRUE)

  # Extract MRMkit's "QC" info
  d_mrmkit_featureinfo <- d_mrmkit_raw %>%
    dplyr::filter(.data$name %in% c("Q1", "Q3", "RT", "D-ratio")) %>%
    tidyr::pivot_longer(-.data$name, names_to = "SOURCE_FEATURE_NAME", values_to = "value") %>%
    tidyr::pivot_wider(names_from = "name" ,values_from = "value") %>%
    dplyr::mutate(SOURCE_FEATURE_NAME = dplyr::if_else(stringr::str_detect(.data$SOURCE_FEATURE_NAME, "RT"), stringr::str_squish(stringr::str_extract(.data$SOURCE_FEATURE_NAME, ".*(?= RT)")),.data$SOURCE_FEATURE_NAME)) %>%
    dplyr::rename(PRECURSOR_MZ = .data$Q1, PRODUCT_MZ = .data$Q3) %>%
    dplyr::mutate(dplyr::across(dplyr::any_of(c("PRECURSOR_MZ", "PRODUCT_MZ", "RT")), as.numeric))


  d_mrmkit_data <- d_mrmkit_raw %>%
    dplyr::filter(!.data$name %in% c("Q1", "Q3", "RT", "D-ratio")) %>%
    dplyr::mutate(RUN_ID = dplyr::row_number(), .before = .data$name) %>%
    tidyr::pivot_longer(-.data$RUN_ID:-.data$name, names_to = "SOURCE_FEATURE_NAME", values_to = "Intensity") %>%
    dplyr::rename(DATAFILE_NAME  = .data$name) %>%
    dplyr::mutate(DATAFILE_NAME = stringr::str_remove(.data$DATAFILE_NAME, stringr::regex("\\.mzML|\\.d|\\.raw|\\.wiff|\\.lcd", ignore_case = TRUE))) |>
    dplyr::mutate(DATAFILE_NAME = stringr::str_squish(.data$DATAFILE_NAME)) |>
    dplyr::mutate(SOURCE_FEATURE_NAME = dplyr::if_else(stringr::str_detect(.data$SOURCE_FEATURE_NAME, "RT"), stringr::str_squish(stringr::str_extract(.data$SOURCE_FEATURE_NAME, ".*(?= RT)")),.data$SOURCE_FEATURE_NAME)) %>%
    dplyr::mutate(dplyr::across(.data$Intensity, as.numeric)) %>%
    dplyr::left_join(d_mrmkit_featureinfo, by = "SOURCE_FEATURE_NAME") %>%
    dplyr::relocate(.data$Intensity, .after = dplyr::last_col()) %>%
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  if(use_mrmkit_normdata)
    d_mrmkit_data |> filter(str_detect(.data$DATAFILE_NAME, "^norm_"))
  else
    d_mrmkit_data |> filter(!str_detect(.data$DATAFILE_NAME, "^norm_"))

}



#' Reads a long CSV file with Feature Intensities
#'
#' @param data MidarExperiment object
#' @param raw_area_file file path of MRMkit csv output with raw peak areas
#' @param final_results_file file path of MRMkit csv output with final processed normalized peak areas
#' @param silent No comments printed
#' @return A tibble in the long format
#' @export
read_MRMkit_results <- function(data, raw_area_file, final_results_file, silent = FALSE) {

  sample_def_cols = c(
    "DATAFILE_NAME",
    "SampleName",
    "AcqTimeStamp",
    "SampleType",
    "VialPosition",
    "InjVolume",
    "SampleType",
    "RUN_ID",
    "ANALYSIS_ID",
    "SOURCE_FEATURE_NAME",
    "Area",
    "RT",
    "FWHM",
    "SN"
  )

  d_raw <- read_MRMkit_raw_area_CSV(raw_area_file)
  data@dataset_orig <- d_raw

  sep <- if_else(fs::path_ext(final_results_file) == "csv", ",", "\t")

  if (final_results_file !=""){
    d_results <- readr::read_delim(final_results_file, delim = sep, col_names = TRUE, trim_ws = TRUE, progress = TRUE, show_col_types = FALSE) %>%
      dplyr::filter(!is.na(.data$type))%>%
      dplyr::mutate(filename = stringr::str_squish(str_remove(.data$filename, stringr::regex("\\.mzML", ignore_case = TRUE)))) %>%
      dplyr::rename(DATAFILE_NAME = .data$filename) %>%
      #dplyr::mutate(RUN_ID = row_number(), .before = ANALYSIS_ID) %>%
      dplyr::select(-tidyselect::any_of(c("batch", "type"))) %>%
      tidyr::pivot_longer(-.data$DATAFILE_NAME, names_to = "SOURCE_FEATURE_NAME", values_to = "normIntensity")

    data@dataset <- data@dataset_orig %>% left_join(d_results, by=c("DATAFILE_NAME", "SOURCE_FEATURE_NAME"))

  }
  data@status_processing <- "Raw Data"
  data
}



