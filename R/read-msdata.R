#' @title Imports an Agilent MassHunter Quant CSV file
#' @description
#' Imports .csv file(s) exported from `Agilent MassHunter (MH) Quantitative Anal  ysis` containing peak integration results.
#' Samples must be in rows, features/compounds in columns and must contain either peak areas, peak heights or intensities.
#' Additional columns, such as RT (retention time), FWHM, PrecursorMZ, and CE will be imported and will available from the `MidarExperiment` object for downstream analyses.
#'
#' When a patch to folder containing MH .csv files is provided, all files are imported and merged into one raw dataset. This is useful, e.g. when importing datasets that are pre-processing in blocks resulting in different files.
#' Each feature and raw data file pair must only occur once within and across all .csv source data files, duplicated return an error.
#'
#' @param data MidarExperiment object
#' @param path One or more file names with path, or a folder path, which case all *.csv files in this folder will be read.
#' @param file_format File format of the MassHunter export. One of "csv", "xls"
#'
#' @importFrom methods validObject
#' @importFrom fs is_dir path_tidy file_exists dir_ls
#'
#' @return MidarExperiment object
#' @examples
#' path_csvfile <- system.file("extdata", "Example_MHQuant_1.csv", package = "midar")
#' mexp <- MidarExperiment()
#' mexp <- rawdata_import_agilent(mexp, path_csvfile)
#' mexp

#' @export

rawdata_import_agilent <- function(data, path, file_format = "csv", expand_qualifier_names = TRUE, silent = FALSE) {
  if (file_format == "csv") {
    import_analysis(data, path, "read_masshunter_csv", "*.csv", expand_qualifier_names = expand_qualifier_names, silent = silent)
  } else {
    stop(glue::glue("This function currently only supports MH exports in the '*.csv' format, '{file_format}' is not supported)"))
  }
}


import_analysis <- function(data, path, import_function, file_ext, ...) {

  if (!fs::is_dir(path)) {
    file_paths <- fs::path_tidy(path)
  } else {
    file_paths <- fs::dir_ls(path, glob = file_ext)
  }


  if (!all(fs::file_exists(file_paths))) stop("One or more given files do not exist. Please check file paths.")
  if (any(duplicated(file_paths))) stop("One or more given files are replicated. Please check file paths.")

  names(file_paths) <- file_paths
  args <- list(...)

  d_temp <- file_paths |>
    purrr::map_dfr(.f = \(x) do.call(what = import_function, append(x, args)), .id = "data_source")


  # Test if analysis_ids, feature_names, and values are replicated
  if (nrow(d_temp) > nrow(d_temp |> distinct(.data$analysis_id, .data$feature_name, .keep_all = FALSE))) {
    has_duplicated_id <- TRUE
    if (nrow(d_temp) > nrow(d_temp |> distinct(.data$analysis_id, .data$feature_name, .keep_all = TRUE))) {
      has_duplicated_id_values <- TRUE
    } else {
      has_duplicated_id_values <- FALSE
    }
  } else {
    has_duplicated_id <- FALSE
  }

  if (has_duplicated_id) {
    if (has_duplicated_id_values) {
      stop(glue::glue("Dataset(s) contains replicated reportings (analysis and feature pairs) with identical intensity values. Please check dataset(s)."))
    } else {
      stop(glue::glue("Dataset(s) contains replicated reportings (analysis and feature pairs) with different intensity values.Please check dataset(s)."))
    }
  }

  data@dataset_orig <- d_temp

  data@dataset_orig <- data@dataset_orig |> dplyr::rename(feature_intensity = "feature_area")
  # TODO: excl_unannotated_analyses below

  check_integrity(data, excl_unannotated_analyses = FALSE)
  # stopifnot(methods::validObject(data))
  data@status_processing <- "Raw Data"
  data
}



#' Reads and parses one Agilent MassHunter Quant CSV result file
#'
#' @param path File path of MassHunter Quant CSV file
#' @param silent Suppress messages
#' @param expand_qualifier_names If TRUE, original qualifier names will be renamed by adding the quantifier name in front and placing qualifier name into square brackets(e.g. `Qualifier (422.3 -> 113.0)` ransition names of quantifier will be added to qualifier names
#' @return A tibble with the parse results in the long format


read_masshunter_csv <- function(path, expand_qualifier_names = TRUE, silent = FALSE) {
  # if(!silent) print(glue::glue("Reading [{basename(path)}] ..."))
  # if (shiny::isRunning())
  #   incProgress(1 / length(n_datafiles), detail = paste0(", basename(file)))
  #
  # Read Agilent MassHunter Quant Export file (CSV)

  suppressWarnings(suppressMessages(
    datWide <-
      readr::read_csv(
        file = path,
        col_names = FALSE,
        na = c("#N/A", "NULL", "NA"),
        trim_ws = TRUE,
        col_types = readr::cols(.default = "c"),
        locale = readr::locale(encoding = "UTF-8"),
        progress = TRUE
      )
  ))
  warnings_datWide <- readr::problems(datWide)

  # Remove text that is not required and remove dot chars that interfere later with the conversion wide to long
  # TODO: Convert to tidyverse functions
  datWide[2, ] <- lapply(datWide[2, ], \(y) gsub("\\. ", "", y))
  datWide[2, ] <- lapply(datWide[2, ], \(y) gsub("\\.", "", y))
  datWide[2, ] <- lapply(datWide[2, ], \(y) gsub("/", "", y))



  datWide <- datWide |> dplyr::add_row(.after = 1)

  if (datWide[1, 1] != "Sample") stop("Error parsing this file. It may in unsupported format, e.g. with features/analytes in rows, or corrupt. Please try re-export from Masshunter with samples in rows and features/analytes in columns.")

  feature_name_tbl <- tibble::tibble(`_prefixXXX_` = datWide[1, ] |> unlist() |> dplyr::na_if(""))

  # if parameter set, then use prefix the feature name and modify the qualifier name
  if (!expand_qualifier_names) {
    feature_name_tbl <- feature_name_tbl |>
      tidyr::fill("_prefixXXX_")
  } else {
    feature_name_tbl <- feature_name_tbl |>
      mutate(temp = .data$`_prefixXXX_`) |>
      tidyr::fill("temp") |>
      mutate(temp = if_else(!str_detect(.data$temp, "Qualifier \\(") & expand_qualifier_names, "", .data$temp)) |>
      mutate(`_prefixXXX_` = if_else(str_detect(.data$`_prefixXXX_`, "Qualifier \\(") & expand_qualifier_names, NA_character_, .data$`_prefixXXX_`)) |>
      tidyr::fill("_prefixXXX_") |>
      mutate(temp = str_replace(.data$temp, "Qualifier \\(", "[QUAL ")) |>
      mutate(temp = str_replace(.data$temp, "\\)", "]")) |>
      mutate(`_prefixXXX_` = paste0(.data$`_prefixXXX_`, " ", .data$temp)) |>
      select(!"temp")
  }


  datWide[1, ] <- feature_name_tbl |>
    unlist() |>
    as.list()

  datWide[1, ] <- replace(datWide[1, ], stringr::str_detect(string = datWide[1, ], pattern = "_prefixXXX_"), "")
  datWide[2, ] <- replace(datWide[1, ], !stringr::str_detect(string = datWide[1, ], pattern = "_prefixXXX_"), "")

  datWide[1, ] <- tibble::tibble(`_prefixXXX_` = datWide[1, ] |> unlist() |> dplyr::na_if("")) |>
    tidyr::fill("_prefixXXX_") |>
    unlist() |>
    as.list()

  # Concatenate rows
  datWide[1, ] <- paste(datWide[1, ], datWide[2, ], sep = " ") |>
    stringr::str_squish() |>
    as.list()
  datWide <- datWide[-2, ]


  # Fill in transition name in empty columns (different parameters of the same transition)
  # Dealing with exported "Qualifier" results: Adds the corresponding quantifier name (the first column of the group) with the Tag "QUAL" and the transition in front e.g. "Sph d18:0 [QUAL 302.3>266.2]"

  # Remove columns with no column names, as they are undefined (seems to be only Outlier Summary and Quantitation Message Summary)
  datWide <- datWide |> dplyr::select(-tidyselect::where(~ .x[2] == ""))

  # Concatenate rows containing parameters + transitions to the form parameter.sample and parameter.transition headers. This will later be converted with reshape()
  colnames(datWide) <- paste(datWide[2, ], datWide[1, ], sep = "\t")
  datWide <- datWide[-1:-2, ]

  # prefix columns from Method and Results Group
  datWide <- datWide |>
    dplyr::rename_with(.fn = ~ stringr::str_c("Method_", .x), .cols = dplyr::ends_with(" Method")) |>
    dplyr::rename_with(.fn = ~ stringr::str_c("Results_", .x), .cols = dplyr::ends_with(" Results")) |>
    dplyr::rename_with(.fn = ~ stringr::str_remove(.x, " Method$| Results$"))

  if (expand_qualifier_names) datWide <- datWide |> dplyr::rename_with(.fn = ~ stringr::str_replace(.x, " Method \\[QUAL| Results \\[QUAL", " [QUAL"))


  # Rename some known column header names and remove columns that are not needed or not informative.

  new_colnames <- c(
    "raw_data_filename" = "Data File\tSample",
    "sample_name" = "Name\tSample",
    "acquisition_time_stamp" = "AcqDate-Time\tSample",
    "sample_type" = "Type\tSample",
    "vial_position" = "Pos\tSample",
    "inj_volume" = "Vol\tSample",
    "dilution_factor" = "Dil\tSample",
    "sample_group" = "Sample Group\tSample",
    "instrument_name" = "Instrument Name\tSample",
    "instrument_type" = "Instrument Type\tSample",
    "acq_method_file" = "AcqMethod File\tSample",
    "acq_method_path" = "AcqMethod Path\tSample",
    "data_file_path" = "Data Path\tSample",
    "message_quantitation" = "Quantitation Message\tSample",
    "comment" = "Comment\tSample",
    "completed" = "Completed\tSample",
    "compound_name" = "Name\tCompound"
  )

  datWide <- datWide |> dplyr::rename(dplyr::any_of(new_colnames))

  # TODO: check if all these are caught
  if ("message_quantitation" %in% names(datWide)) stop("Field 'Quantitation Message' currently not supported: Please re-export your data in MH without this field.")
  if ("compound_name" %in% names(datWide)) stop("Compound table format is currently not supported. Please re-export your data in MH with compounds as columns.")
  if (!"raw_data_filename" %in% names(datWide)) stop("Error parsing this Masshunter .csv file. The file may be in an unsupported format or corrupt. Please try re-export from Masshunter.")
  if (nrow(warnings_datWide) > 0) stop("Unknown format, or data file is corrupt. Please try re-export from MH.")
  # Remove ".Sample" from remaining sample description headers and remove known unused columns
  datWide <-
    datWide[, !(names(datWide) %in% c("NA\tSample", "Level\tSample", "Sample"))]
  names(datWide) <- sub("\\\tSample", "", names(datWide))

  keep.cols <- names(datWide) %in% c("")
  datWide <- datWide[!keep.cols]



  # Transform wide to long (tidy) table format

  datWide <- datWide |>
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::any_of(c("acquisition_time_stamp")),
        .fns = \(x) lubridate::parse_date_time(x, c("mdy_HM", "dmy_HM", "ymd_HM", "ydm_HM", "mdy_HM %p", "dmy_HM %p", "ymd_HM %p", "ydm_HM %p"))
      ),
      dplyr::across(
        .cols = dplyr::any_of(c("sample_name")),
        .fns = stringr::str_squish
      )
    )

  # obtain list with column names of the columns that define transition values (e.g. "RT Cer d16:0/18:0"). Delimuter is currently tab (\t)
  param_transition_names <-
    colnames(datWide[, -1:-tail(grep("\\\t", colnames(datWide), invert = TRUE), 1)])


  # Obtain long table of all param-transition combinations, split param and compund name and then spread values of different param as columns
  datLong <- datWide |>
    dplyr::mutate(file_run_id = dplyr::row_number(), .before = 1) |>
    tidyr::pivot_longer(cols = tidyselect::all_of(param_transition_names), names_pattern = "(.*)\t(.*)$", names_to = c("Param", "feature_name")) |>
    tidyr::pivot_wider(names_from = "Param", values_from = "value")

  # Convert types of knows parameters and fields in the data set

  new_numeric_colnames <- c(
    "feature_rt" = "Results_RT",
    "method_target_rt" = "Method_RT",
    "feature_response" = "Results_Resp",
    "feature_height" = "Results_Height",
    "feature_area" = "Results_Area",
    "feature_response" = "Results_Resp",
    "feature_accuracy" = "Results_Accuracy",
    "feature_conc_calc" = "Results_Calc Conc",
    "feature_conc_final" = "Results_Final Conc",
    "feature_fwhm" = "Results_FWHM",
    "feature_width" = "Results_Width",
    "feature_sn_ratio" = "Results_SN",
    "feature_int_start" = "Results_IntStart",
    "feature_int_end" = "Results_IntEnd",
    "feature_symetry" = "Results_Symmetry",
    "method_precursor_mz" = "Method_Precursor Ion",
    "method_product_mz" = "Method_Product Ion",
    "method_collision_energy" = "Method_Collision Energy",
    "method_fragmentor" = "Method_Fragmentor",
    "method_multiplier" = "Method_Multiplier",
    "method_noise_raw_signal" = "Method_Noise of Raw Signal",
    "inj_volume" = "inj_volume"
  )

  new_int_colnames <- c(
    "method_time_segment" = "Method_TS"
  )

  new_logical_colnames <- c(
    "feature_manual_integration" = "Results_MI"
  )

  new_factor_colnames <- c(
    "method_polarity" = "Method_Ion Polarity"
  )

  new_char_colnames <- c(
    "method_compound_group" = "Method_CmpdGroup",
    "method_compound_id" = "Method_ID",
    "method_integration_method" = "Method_Int",
    "method_integration_parameters" = "Method_IntParms",
    "method_peak_smoothing" = "Method_Smoothing",
    "method_peak_smoothing_gauss_width" = "Method_Smoothing Gaussian Width",
    "method_peak_smoothing_function_width" = "Method_Smoothing Function Width",
    "method_transition" = "Method_Transition",
    "method_ion_source" = "Method_Ion Source",
    "method_type" = "Method_Type",
    "method_noise_algorithm" = "Method_Noise Alg"
  )

  datLong <- datLong |>
    dplyr::rename(dplyr::any_of(c(new_numeric_colnames, new_int_colnames, new_logical_colnames, new_factor_colnames, new_char_colnames)))

  # replace , with . for german systems
  datLong <- datLong |>
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::any_of(names(new_numeric_colnames)),
        .fns = \(x) as.numeric(stringr::str_replace(x, ",", "."))
      ),
      dplyr::across(
        .cols = dplyr::any_of(names(new_int_colnames)),
        .fns = as.integer
      ),
      dplyr::across(
        .cols = dplyr::any_of(names(new_logical_colnames)),
        .fns = as.logical
      ),
      dplyr::across(
        .cols = dplyr::any_of(names(new_char_colnames)),
        .fns = as.character
      ),
      dplyr::across(
        .cols = dplyr::any_of(names(new_factor_colnames)),
        .fns = as.factor
      )
    ) |>
    dplyr::mutate(raw_data_filename = stringr::str_replace(.data$raw_data_filename, "\\.d", "")) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish)) |>
    dplyr::relocate("sample_name", .after = "raw_data_filename")


  datLong <- datLong |>
    mutate(
      integration_qualifier = ((expand_qualifier_names & str_detect(.data$feature_name, " \\[QUAL ")) |
        (!expand_qualifier_names & str_detect(.data$feature_name, "^Qualifier \\("))),
      .after = "feature_name"
    )

    # Defines the field from the raw data to used as analysis id
  datLong <- datLong |>
    mutate(analysis_id = raw_data_filename, .before = 1)


  if (!silent) {
    if (!any(datLong$integration_qualifier)) {
      cli_alert_success(cli::col_green("Imported {length(unique(datLong$raw_data_filename))} samples with {length(unique(datLong$feature_name))} features \n"))
    } else if (any(datLong$integration_qualifier)) {
      cli_alert_success(cli::col_green("Imported {length(unique(datLong$raw_data_filename))} samples with {length(unique(datLong$feature_name))} features ({length(unique(datLong$feature_name[!datLong$integration_qualifier]))} quantifiers, {length(unique(datLong$feature_name[!datLong$integration_qualifier]))} qualifiers) \n"))
    }
  }

  datLong
}

#' Reads MRMkit results
#'
#' @param data MidarExperiment object
#' @param path file path of MRMkit csv output with raw peak areas
#' @return A tibble in the long format
#' @export
import_mrmkit <- function(data, path) {
  import_analysis(data, path, "read_mrmkit_result", "*.tsv|*.csv")
}


#' Reads a long CSV file with Feature Intensities
#'
#' @param path File name of the MRMkit result file (*.tsv or *.csv)
#' @param use_normalized_data Import raw peak areas or normalized peak areas from the file
#' @param silent No comments printed
#' @return A tibble in the long format
#' @export
read_mrmkit_result <- function(path, use_normalized_data = FALSE, silent = FALSE) {

  ext_file <- tolower(fs::path_ext(path))

  sep <- case_when(
    ext_file == "csv" ~ "csv",
    ext_file == "tsv" ~ "\t",
    TRUE ~ NA_character_
  )

  if (is.na(sep)) stop("Data file type/extension not supported. Please re-export results from MRMkit.")

  d_mrmkit_raw <- readr::read_delim(path,
    delim = sep, col_types = readr::cols(.default = "c"),
    col_names = TRUE, trim_ws = TRUE
  )


  # Extract MRMkit's "QC" info
  d_mrmkit_featureinfo <- d_mrmkit_raw |>
    dplyr::filter(.data$name %in% c("Q1", "Q3", "RT", "D-ratio")) |>
    tidyr::pivot_longer(-.data$name, names_to = "feature_name", values_to = "value") |>
    tidyr::pivot_wider(names_from = "name", values_from = "value") |>
    dplyr::mutate(feature_name = dplyr::if_else(stringr::str_detect(.data$feature_name, "RT"), stringr::str_squish(stringr::str_extract(.data$feature_name, ".*(?= RT)")), .data$feature_name)) |>
    dplyr::rename(precursor_mz = .data$Q1, product_mz = .data$Q1) |>
    dplyr::mutate(dplyr::across(dplyr::any_of(c("precursor_mz", "product_mz", "RT")), as.numeric))


  d_mrmkit_data <- d_mrmkit_raw |>
    dplyr::filter(!.data$name %in% c("Q1", "Q3", "RT", "D-ratio")) |>
    dplyr::mutate(file_run_id = dplyr::row_number(), .before = "name") |>
    tidyr::pivot_longer(-.data$file_run_id:-.data$name, names_to = "feature_name", values_to = "value") |>
    dplyr::rename(raw_data_filename = .data$name) |>
    dplyr::mutate(raw_data_filename = stringr::str_remove(.data$raw_data_filename, stringr::regex("\\.mzML$|\\.d$|\\.raw$|\\.wiff$|\\.wiff2$|\\.lcd$", ignore_case = TRUE))) |>
    dplyr::mutate(raw_data_filename = stringr::str_squish(.data$raw_data_filename)) |>
    dplyr::mutate(feature_name = dplyr::if_else(stringr::str_detect(.data$feature_name, "RT"), stringr::str_squish(stringr::str_extract(.data$feature_name, ".*(?= RT)")), .data$feature_name)) |>
    dplyr::mutate(dplyr::across(.data$value, as.numeric)) |>
    dplyr::left_join(d_mrmkit_featureinfo, by = "feature_name") |>
    dplyr::relocate(.data$value, .after = dplyr::last_col()) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  if (use_normalized_data) {
    d_mrmkit_data <- d_mrmkit_data |>
      filter(str_detect(.data$raw_data_filename, "^norm_")) |>
      dplyr::rename(norm_intensity = "value")
  } else {
    d_mrmkit_data <- d_mrmkit_data |>
      filter(!str_detect(.data$raw_data_filename, "^norm_")) |>
      dplyr::rename(feature_area = "value")
  }

  names(d_mrmkit_data) <- tolower(names(d_mrmkit_data))

  # sep <- if_else(fs::path_ext(final_results_file) == "csv", ",", "\t")

  # if (final_results_file !=""){
  #   d_results <- readr::read_delim(final_results_file, delim = sep, col_names = TRUE, trim_ws = TRUE, progress = TRUE, show_col_types = FALSE) |>
  #     dplyr::filter(!is.na(.data$type))|>
  #     dplyr::mutate(path = stringr::str_squish(str_remove(.data$path, stringr::regex("\\.mzML$|\\.d$|\\.raw$|\\.wiff$|\\.wiff2$|\\.lcd$", ignore_case = TRUE)))) |>
  #     dplyr::rename(raw_data_filename = .data$path) |>
  #     #dplyr::mutate(run_id = row_number(), .before = ANALYSIS_ID) |>
  #     dplyr::select(-tidyselect::any_of(c("batch", "type"))) |>
  #     tidyr::pivot_longer(-.data$raw_data_filename, names_to = "feature_name", values_to = "feature_norm_intensity")
  #
  #   data@dataset <- data@dataset_orig |> left_join(d_results, by=c("raw_data_filename", "feature_name"))
  #
  # }

  if (!silent) {
    txt <- glue::glue("\u2713 Imported {length(unique(d_mrmkit_data$raw_data_filename))} samples with {length(unique(d_mrmkit_data$feature_name))} features \n")
    writeLines(col_green(txt))
  }

  d_mrmkit_data
}


#' Reads a generic wide .csv or .xls sheet with analysis results data
#'
#' @description
#' Imports a tables with the column being the analysis identifier (e.g. sample name) and the following colums being measured features. The values in the table
#'
#'
#' @param file File name and path of a plain wide-format CSV or XLS file
#' @param value_type Parameter type of the feature, i.e. "area", "height", "intensity", "norm_area", "norm_height", "norm_intensity", "feature_conc", "rt", "fwhm" or "width"
#' @param sheet Sheet name in case an Excel file (.xls, .xlsx, .xlsm) is imported
#' @param silent Suppress messages

#' @return A tibble in the long format
#' @export
#'
read_analysisresult_table <- function(file, value_type = c("area", "height", "intensity", "norm_area", "norm_height", "norm_intensity", "feature_conc", "rt", "fwhm", "width"), sheet = "", silent = FALSE) {
  value_type <- match.arg(value_type)

  value_type <- dplyr::case_match(
    value_type,
    "area" ~ "feature_area",
    "height" ~ "feature_height",
    "intensity" ~ "feature_intensity",
    "norm_area" ~ "feature_norm_area",
    "norm_height" ~ "feature_norm_height",
    "norm_intensity" ~ "feature_norm_intensity",
    "feature_conc" ~ "feature_conc",
    "rt" ~ "feature_rt",
    "fwhm" ~ "feature_fwhm",
    "width" ~ "feature_width",
    .default = NULL
  )

  var_value_type <- rlang::ensym(value_type)
  ext <- fs::path_ext(file)

  if (ext == "csv") {
    d <- readr::read_csv(file, col_names = TRUE, trim_ws = TRUE, progress = FALSE, na = c("n/a", "N/A", "NA", "na", "ND", "N.D.", "n.d."), col_types = "cn", locale = readr::locale(encoding = "UTF-8"))
  } else if (ext == "xls" || ext == "xlsx" || ext == "xlm" || ext == "xlmx") {
    if (sheet == "") stop("Please define sheet name via the `sheet` parameter")
    # nms <- names(readxl::read_excel(path = file, sheet = sheet, trim_ws = TRUE, progress = TRUE, na = c("n/a", "N/A", "NA"), n_max = 0))
    # ct <- c("text", rep("numeric",length(nms)-1))
    # d <- readxl::read_excel(path = file, sheet = sheet, trim_ws = TRUE, progress = TRUE, na = c("n/a", "N/A", "NA"), col_types = ct)
    d <- .read_generic_excel(file, sheet)
  } else {
    stop("Invalid file format. Only CSV, XLS and XLSX supported.")
  }

  d |>
    dplyr::mutate(run_id = dplyr::row_number(), .before = 1) |>
    tidyr::pivot_longer(cols = -1:-2, names_to = "feature_name", values_to = value_type) |>
    dplyr::rename(analysis_id = 2) |>
    mutate({{ value_type }} := as.numeric(!!var_value_type))
}




#' Reads a long CSV file with Feature Intensities
#'
#' @param file File name and path of a plain long-format CSV file
#' @param analysis_id_col Column to be used as analysis_id
#' @param feature_name_col Column to be used feature_name
#' @param silent Suppress messages
#'
#' @return A tibble in the long format
#' @export
rawdata_import_mrmkit <- function(file, analysis_id_col = NULL, feature_name_col = NULL, silent = FALSE) {
  analysis_inf_cols <- c(
    "analysis_id",
    "raw_data_filename",
    "feature_name",
    "sample_name",
    "acquisition_time_stamp",
    "sample_type",
    "vial_position",
    "inj_volume",
    "sample_type",
    "run_id"
  )

  quant_cols <- c(
    "feature_area",
    "feature_height",
    "feature_response",
    "feature_intensity",
    "feature_rt",
    "feature_fwhm",
    "feature_sn_ratio"
  )

  d <- readr::read_csv(file, col_names = TRUE, trim_ws = TRUE, locale = readr::locale(encoding = "UTF-8"), progress = TRUE)

  if (!is.null(analysis_id_col) && analysis_id_col %in% names(d)) stop("Analysis Id column not defined. ")

  dplyr::select(tidyselect::any_of(.data$sample_def_cols), "feature_name", feature_intensity = {{ .data$field }})

  d
}




#' @title internal method to read excel sheets
#' @param path excel file name
#' @return tibble table
#' @noRd
.read_generic_excel <- function(path, sheetname) {
  nms <- names(readxl::read_excel(path = path, sheet = sheetname, trim_ws = TRUE, progress = TRUE, na = c("n/a", "N/A", "NA"), n_max = 0))
  n_col <- length(nms) - 1
  c_numeric <- rep("numeric", n_col)
  ct <- c("text", c_numeric)
  cat("c_numeric contains:", ct, "\n")
  readxl::read_excel(path = path, sheet = sheetname, trim_ws = TRUE, progress = TRUE, na = c("n/a", "N/A", "NA"), col_types = ct)
}


#' @title internal method to read csv files
#' @param path csv file name
#' @return tibble table
#' @noRd
# TODO remove this test function
.test_mult <- function(a, b) {
  nms <- c(1, 2, 3)
  n_col <- length(nms)
  rep("numeric", n_col)
}
