
#' @title Imports Agilent MassHunter Quantitative Analysis CSV files
#' @description
#' Imports .csv files exported from Agilent MassHunter Quantitative Analysis software, containing peak integration results.
#' The input files must have anlyses (samples) in rows, features/compounds in columns, and either peak areas, peak heights, or response as the values.
#' Additional columns, such as retention time (RT), full-width at half-maximum (FWHM), precursor m/z (PrecursorMZ), and collision energy (CE), will also be imported and made available in the `MidarExperiment` object for downstream analyses.
#'
#' When a directory path is provided, all matching .csv files in that directory will be imported and merged into a single dataset. This is useful when importing datasets that were pre-processed in blocks, resulting in multiple files.
#' Each unique combination of feature and raw data file must only occur once across all source data files. Duplicate combinations will result in an error.
#'
#' @param data MidarExperiment object
#' @param path One or more file paths, or a directory path (in which case all matching files will be imported)
#' @param import_metadata Logical, whether to extract and add metadata from the analysis result file
#' @param expand_qualifier_names Logical, whether to add the quantifier name in front of the qualifier name (the latter only has the m/z transition values)
#' @param silent Logical, whether to suppress most notifications
#' @return MidarExperiment object with the imported data
#' @examples
#' mexp <- MidarExperiment()

#' file_path = system.file("extdata", "Example_MHQuant_1.csv", package = "midar")
#'
#' mexp <- import_data_masshunter(
#'   data = mexp,
#'   path = file_path,
#'   import_metadata = TRUE,
#'   expand_qualifier_names = TRUE)
#'
#' print(mexp)
#'
#' @export

import_data_masshunter <- function(data = NULL, path, import_metadata, expand_qualifier_names = TRUE, silent = FALSE) {
  check_data(data)
  #if (fs::path_ext(path) == "csv") {
    data <- import_data_main(data, path, "parse_masshunter_csv", "*.csv", expand_qualifier_names = expand_qualifier_names, silent = silent)
    data <- set_intensity_var(data, variable_name = NULL, auto_select = TRUE, warnings = TRUE, "feature_area", "feature_response", "feature_height")
    if (import_metadata) data <- import_metadata_from_data(data, qc_type_column_name = "sample_type")
    data
  #} else {
  #  cli::cli_abort(glue::glue("This function currently only supports MH exports in the '*.csv' format, '{file_format}' is not supported)"))
  #}
}

#' @title Imports MRMkit peak integration results
#' @description
#' Imports tabular data files (*.tsv) generated from `MRMkit` containing peak integration results.
#' The input files must be in a long format with columns for the raw data file name, feature ID, peak intensity, and other parameters.
#' Additional information, such as retention time, FWHM, precursor/product MZ, and CE will also be imported and made available in the `MidarExperiment` object for downstream analyses.
#'
#' When a directory path is provided, all matching files in that directory will be imported and merged into a single dataset.
#' This is useful when importing datasets that were pre-processed in blocks, resulting in multiple files.
#' Each unique combination of feature and raw data file must only occur once across all source data files.
#' Duplicate combinations will result in an error.
#'
#' @param data MidarExperiment object
#' @param path One or more file paths, or a directory path (in which case all matching files will be imported)
#' @param import_metadata Logical, whether to import additional metadata columns (e.g., `batch_id`, `qc_type`)
#' @param silent Logical, whether to suppress most notifications
#' @return MidarExperiment object with the imported data
#' @examples
#' mexp <- MidarExperiment()
#'
#' file_path = system.file("extdata", "sPerfect_MRMkit.tsv", package = "midar")
#'
#' mexp <- import_data_mrmkit(
#'   data = mexp,
#'   path = file_path,
#'   import_metadata = TRUE)
#
#' print(mexp)
#' @export
import_data_mrmkit <- function(data = NULL, path, import_metadata, silent = FALSE) {
  check_data(data)
  data <- import_data_main(data = data, path = path, import_function = "parse_mrmkit_result", file_ext = "*.tsv|*.csv", silent = FALSE)
  data <- set_intensity_var(data, variable_name = NULL, auto_select = TRUE, warnings = TRUE, "feature_area", "feature_height")

  if (import_metadata) data <- import_metadata_from_data(data, qc_type_column_name = "sample_type")
  data
}

#' @title Import plain analysis results
#' @description
#' Imports .csv file(s) reprensenting analysis results in a wide format with samples in rows and features in columns.
#' The first column must contain the analysis identifier, subsequent columns can be metadata columns, followed by value of different features.
#' The table format must be in the long format with columns for the raw data file name, feature ID, and the peak intensity and other parameters.
#' Additional information, such as retention time, FWHM, precursor/product MZ, and CE will also be imported and will available from the `MidarExperiment` object for downstream analyses.
#'
#' When a path to a folder is provided, all .scv files in that folder will be imported and merged into one raw dataset. This is useful, e.g. when importing datasets that were pre-processed in blocks resulting in different files.
#' Each feature and raw data file pair must only occur once within and across all .csv source data files, duplicated return an error.
#'
#' @param data MidarExperiment object
#' @param path One or more file names with path, or a folder path, which case all *.csv files in this folder will be read.
#' @param variable_name Variable type representing the values in the table. Must be one of "intensity", "norm_intensity", "conc", "area", "height", "response")
#' @param analysis_id_col Column to be used as analysis_id. `NA` (default) used 'analysis_id' if present, or the first column if it contains unique values.
#' @param import_metadata Import additional metadata columns (e.g. batch ID, sample type) and add to the `MidarExperiment` object.
#' Only following metadata column names are supported: "qc_type", "batch_id", "is_quantifier", "is_istd", "run_seq_num", "precursor_mz", "product_mz", "collision_energy"
# #' @param silent Suppress notifications
#' @return MidarExperiment object
#' @examples
#' file_path <- system.file("extdata", "plain_wide_dataset.csv", package = "midar")
#'
#' mexp <- MidarExperiment()
#'
#' mexp <- import_data_csv(
#'   data = mexp,
#'   path = file_path,
#'  variable_name = "conc",
#'  import_metadata = TRUE)
#'
#' print(mexp)
#'
#' @export

import_data_csv <- function(data = NULL, path, variable_name, analysis_id_col = NA, import_metadata) {
  check_data(data)
  data <- import_data_main(data = data, path = path, import_function = "parse_plain_csv", file_ext = "*.csv", silent = FALSE, variable_name, analysis_id_col, import_metadata)
  data <- set_intensity_var(data, variable_name = paste0("feature_", str_remove(variable_name, "feature_")), auto_select = TRUE, warnings = TRUE, "feature_area", "feature_height", "feature_conc")
  if (import_metadata) data <- import_metadata_from_data(data, qc_type_column_name = "qc_type")

  data
}



import_data_main <- function(data = NULL, path, import_function, file_ext, silent, ...) {
  check_data(data)
  if (!fs::is_dir(path)) {
    file_paths <- fs::path_tidy(path)
  } else {
    file_paths <- fs::dir_ls(path, glob = file_ext)
  }

  if (!all(fs::file_exists(file_paths))) cli::cli_abort("One or more given files do not exist. Please check file paths.")
  if (any(duplicated(file_paths))) cli::cli_abort("One or more given files are replicated. Please check file paths.")

  names(file_paths) <- file_paths
  args <- list(...)

  d_raw <- file_paths |>
    purrr::map_dfr(.f = \(x) do.call(what = import_function, append(x, args)), .id = "data_source") |>
    relocate("analysis_id", "data_source")

  d_runorder <- d_raw |>
    select("analysis_id") |>
    distinct() |>
    mutate(run_seq_num = row_number())

  d_raw <- d_raw |>
    left_join(d_runorder, by = "analysis_id") |>
    relocate("run_seq_num", .before = 1)

  # VERIFY DATA, i.e. analysis_ids, feature_ids, and values are replicated ===
  ## which can be result of multiple imports of the same/overlapping data or due to parsing error

  if (nrow(d_raw) > nrow(d_raw |> distinct(.data$analysis_id, .data$feature_id, .keep_all = FALSE))) {
    has_duplicated_id <- TRUE
    if (nrow(d_raw) > nrow(d_raw |> distinct(.data$analysis_id, .data$feature_id, .keep_all = TRUE))) {
      has_duplicated_id_values <- TRUE
    } else {
      has_duplicated_id_values <- FALSE
    }
  } else {
    has_duplicated_id <- FALSE
  }

  if (has_duplicated_id) {
    if (has_duplicated_id_values) {
      cli::cli_abort(glue::glue("Imported data contains replicated reportings (analysis and feature pairs) with **identical** intensity values. Please check imported dataset(s)."))
    } else {
      cli::cli_abort(glue::glue("Imported data contains replicated reportings (analysis and feature pairs) with **different intensity values. Please check imported dataset(s)."))
    }
  }

  data@dataset_orig <- dplyr::bind_rows(pkg.env$table_templates$dataset_orig_template, d_raw)



  # TODO: excl_unmatched_analyses below

  check_integrity(data, excl_unmatched_analyses = FALSE)
  # stopifnot(methods::validObject(data))

  if (!silent) {
    if (!any(data@dataset_orig$integration_qualifier)) {
      cli_alert_success(cli::col_green("Imported {length(unique(data@dataset_orig$analysis_id))} analyses with {length(unique(data@dataset_orig$feature_id))} features"))
    } else if (any(data@dataset_orig$integration_qualifier)) {
      cli_alert_success(cli::col_green("Imported {length(unique(data@dataset_orig$analysis_id))} analyses with {length(unique(data@dataset_orig$feature_id))} features ({length(unique(data@dataset_orig$feature_id[!data@dataset_orig$integration_qualifier]))} quantifiers, {length(unique(data@dataset_orig$feature_id[!data@dataset_orig$integration_qualifier]))} qualifiers)"))
    }
  }

  data@status_processing <- "Raw Data"

  data
}



#' Reads and parses one Agilent MassHunter Quant CSV result file
#'
#' @param path File path of MassHunter Quant CSV file
#' @param silent Suppress messages
#' @param expand_qualifier_names If TRUE, original qualifier names will be renamed by adding the quantifier name in front and placing qualifier name into square brackets(e.g. `Qualifier (422.3 -> 113.0)` transition names of quantifier will be added to qualifier names
#' @return A tibble with the parse results in the long format
#' @examples
#' file_path = system.file("extdata", "Example_MHQuant_1.csv", package = "midar")
#'
#' tbl <- parse_masshunter_csv(
#'   path = file_path,
#'   expand_qualifier_names = TRUE)
#'
#' head(tbl)
#' @export

parse_masshunter_csv <- function(path, expand_qualifier_names = TRUE, silent = FALSE) {
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

  if (datWide[1, 1] != "Sample") cli::cli_abort("Error parsing this file. It may in unsupported format, e.g. with features/analytes in rows, or corrupt. Please try re-export from Masshunter with samples in rows and features/analytes in columns.")

  feature_id_tbl <- tibble::tibble(`_prefixXXX_` = datWide[1, ] |> unlist() |> dplyr::na_if(""))

  # if parameter set, then use prefix the feature name and modify the qualifier name
  if (!expand_qualifier_names) {
    feature_id_tbl <- feature_id_tbl |>
      tidyr::fill("_prefixXXX_")
  } else {
    feature_id_tbl <- feature_id_tbl |>
      mutate(temp = .data$`_prefixXXX_`) |>
      tidyr::fill("temp") |>
      mutate(temp = if_else(!str_detect(.data$temp, "Qualifier \\(") & expand_qualifier_names, "", .data$temp)) |>
      mutate(`_prefixXXX_` = if_else(str_detect(.data$`_prefixXXX_`, "Qualifier \\(") & expand_qualifier_names, NA_character_, .data$`_prefixXXX_`)) |>
      tidyr::fill("_prefixXXX_") |>
      mutate(temp = str_replace(.data$temp, "Qualifier \\(", "[QUAL ")) |>
      mutate(temp = str_replace(.data$temp, "\\)", "]")) |>
      mutate(`_prefixXXX_` = paste0(.data$`_prefixXXX_`, " ", .data$temp)) |>
      dplyr::select(!"temp")
  }


  datWide[1, ] <- feature_id_tbl |>
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
    "sample_level" = "Level\tSample",
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
  if ("message_quantitation" %in% names(datWide)) cli::cli_abort("Field 'Quantitation Message' currently not supported: Please re-export your data in MH without this field.")
  if ("compound_name" %in% names(datWide)) cli::cli_abort("Compound table format is currently not supported. Please re-export your data in MH with compounds as columns.")
  if (!"raw_data_filename" %in% names(datWide)) cli::cli_abort("The `Data File` column is missing or the file format is unsupported format. Please try re-export with `Data File` from Masshunter.")
  if (nrow(warnings_datWide) > 0) cli::cli_abort("Unknown format, or data file is corrupt. Please try re-export from MH.")
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
        .fns = \(x) lubridate::parse_date_time(x, c("mdy_HM", "dmy_HM", "ymd_HM", "ydm_HM", "mdy_HM %p", "dmy_HM %p", "ymd_HM %p", "ydm_HM %p", "%y-%m-%d %H:%M:%S"), quiet=TRUE)
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
    dplyr::mutate(file_run_seq_num = dplyr::row_number(), .before = 1) |>
    tidyr::pivot_longer(cols = tidyselect::all_of(param_transition_names), names_pattern = "(.*)\t(.*)$", names_to = c("Param", "feature_id")) |>
    tidyr::pivot_wider(names_from = "Param", values_from = "value")


  # Convert types of knows parameters and fields in the data set

  new_numeric_colnames <- c(
    "feature_rt" = "Results_RT",
    "method_target_rt" = "Method_RT",
    "feature_response" = "Results_Resp",
    "feature_height" = "Results_Height",
    "feature_area" = "Results_Area",
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
    "method_conc_expected" = "Method_Exp Conc",
    "method_fragmentor" = "Method_Fragmentor",
    "method_multiplier" = "Method_Multiplier",
    "method_noise_raw_signal" = "Method_Noise of Raw Signal",
    "inj_volume" = "inj_volume"
  )



  new_int_colnames <- c("method_time_segment" = "Method_TS")
  new_logical_colnames <- c("feature_manual_integration" = "Results_MI")
  new_factor_colnames <- c("method_polarity" = "Method_Ion Polarity")

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
      integration_qualifier = ((expand_qualifier_names & str_detect(.data$feature_id, " \\[QUAL ")) |
        (!expand_qualifier_names & str_detect(.data$feature_id, "^Qualifier \\("))),
      .after = "feature_id"
    )

    # Defines the field from the raw data to used as analysis id
  datLong <- datLong |>
    mutate(analysis_id = .data$raw_data_filename, .before = 1)

  datLong
}


#' Parses MRMkit peak integration results into a tibble
#'
#' @param path File name of the MRMkit result file (*.tsv or *.csv)
# #' @param use_normalized_data Import raw peak areas or normalized peak areas from the file
#' @param silent No comments printed
#' @return A tibble in the long format
#' @examples
#'
#' file_path = system.file("extdata", "sPerfect_MRMkit.tsv", package = "midar"),
#'
#' tbl <- parse_mrmkit_result(
#'   path = file_path)
#'
#' head(tbl)
#' @export

# TODO: remove support for norm intensity
parse_mrmkit_result <- function(path, silent = FALSE) {

  ext_file <- tolower(fs::path_ext(path))

  sep <- case_when(
    ext_file == "csv" ~ "csv",
    ext_file == "tsv" ~ "\t",
    TRUE ~ NA_character_
  )
  if (is.na(sep)) cli::cli_abort("Data file type/extension not supported.")

  d_mrmkit_raw <- readr::read_delim(path,
                                    delim = sep, col_types = readr::cols(.default = "c"),
                                    col_names = TRUE, trim_ws = TRUE)


  d_mrmkit_data <- d_mrmkit_raw |>
  dplyr::select(
    "raw_data_filename",
    "acquisition_time_stamp" = "time_stamp",
    "sample_type" = "sample_type",
    "batch_id" = "batch",
    "featurename" = "feature_name",
    "istd_feature_id" = "internal_standard",
    "feature_rt" = "rt_apex",
    "feature_area" = "area",
    #"feature_norm_intensity" = "area_normalized",
    "feature_height" = "height",
    "feature_fwhm" = "FWHM",
    "feature_int_start" = "rt_int_start",
    "feature_int_end" = "rt_int_end",
    "method_polarity" = "polarity",
    "method_precursor_mz" = "precursor_mz",
    "method_product_mz" = "product_mz",
    "method_collision_energy" = "collision_energy",
  )

  # extract and pivot to improve speed
  d_sample <- d_mrmkit_data |>
    dplyr::select("raw_data_filename", "acquisition_time_stamp", "batch_id") |>
    dplyr::distinct() |>
    dplyr::mutate(analysis_id = stringr::str_remove(.data$raw_data_filename, stringr::regex("\\.mzML$|\\.d$|\\.raw$|\\.wiff$|\\.wiff2$|\\.lcd$", ignore_case = TRUE))) |>
    dplyr::mutate(acquisition_time_stamp = lubridate::parse_date_time(.data$acquisition_time_stamp, c("mdy_HM", "dmy_HM", "ymd_HM", "ydm_HM", "mdy_HM %p", "dmy_HM %p", "ymd_HM %p", "ydm_HM %p", "%y-%m-%d %H:%M:%S"), quiet=TRUE)) |>
    dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))

  d_feature <- d_mrmkit_data |>
    dplyr::select(c("featurename", "istd_feature_id")) |>
    distinct() |>
    mutate(feature_id = stringr::str_squish(.data$featurename),
           istd_feature_id = stringr::str_squish(.data$istd_feature_id))

  # finalize data
  d_mrmkit_data <- d_mrmkit_data |>
    dplyr::mutate(dplyr::across(tidyselect::starts_with("feature_"), as.numeric)) |>
    dplyr::mutate(dplyr::across(tidyselect::ends_with("_mz"), as.numeric)) |>
    dplyr::mutate(dplyr::across(tidyselect::ends_with("_energy"), as.numeric)) |>
    dplyr::mutate(integration_qualifier = FALSE) |>
    dplyr::select(-c("acquisition_time_stamp",  "istd_feature_id", "batch_id")) |>
    left_join(d_sample, by = "raw_data_filename") |>
    left_join(d_feature, by = "featurename", keep = FALSE) |>
    dplyr::select(-"featurename") |>
    relocate("analysis_id", "raw_data_filename", "acquisition_time_stamp", "sample_type", "batch_id", "feature_id", "istd_feature_id", "integration_qualifier")

  # if (!use_normalized_data) {
  #   #d_mrmkit_data <- d_mrmkit_data |> mutate(feature_norm_intensity = NA_real_)
  #   d_mrmkit_data <- d_mrmkit_data |> select(-feature_norm_intensity)
  # }

  d_mrmkit_data

}



#' Reads a long CSV file with Feature Intensities
#'
#' @param path path name and path of a plain long-format CSV file
#' @param variable_name Name of the variable representing the values in the table. Must be one of "intensity", "norm_intensity", "conc", "area", "height", "response")
#' @param analysis_id_col Column to be used as analysis_id
#' @param import_metadata Import additional metadata columns (e.g. batch ID, sample type) and add to the `MidarExperiment` object
#' @return A tibble in the long format
#' @examples
#' file_path <- system.file("extdata", "plain_wide_dataset.csv", package = "midar")
#'
#' tbl <- parse_plain_csv(
#'  path = file_path,
#'  variable_name = "conc",
#'  analysis_id_col = "analysis_id"m
#'  import_metadata = TRUE)
#'
#' head(tbl)
#' @export
 parse_plain_csv <- function(path, variable_name, analysis_id_col = NA, import_metadata) {

  if(fs::path_ext(path) != "csv") cli::cli_abort("Only csv files are currently supported.")

  variable_name <- str_remove(variable_name, "feature_")
  rlang::arg_match(variable_name, c("area", "height", "intensity", "response", "conc", "conc_raw", "rt", "fwhm"))
  variable_name <- stringr::str_c("feature_", variable_name)
  variable_name_sym = rlang::sym(variable_name)

  if(import_metadata){
    analysis_metadata_cols <- c(
      "run_seq_num",
      "qc_type",
      "batch_id",
      "is_quantifier",
      "is_istd",
      "precursor_mz",
      "product_mz",
      "collision_energy"
    )
  } else {
    analysis_metadata_cols = NULL
  }

  d <- readr::read_csv(path, col_names = TRUE, trim_ws = TRUE, locale = readr::locale(encoding = "UTF-8"), col_types = "n", progress = FALSE)

  if (!is.null(analysis_id_col) & !is.na(analysis_id_col)){
    if (!analysis_id_col %in% names(d)) {
      cli::cli_abort("No column with the name `{analysis_id_col}` found in the data file.")
    }
  } else {
    if ("analysis_id" %in% names(d) & (length(unique(d$analysis_id)) == nrow(d))){
      analysis_id_col <- "analysis_id"
    } else {
      if (all(is.character(d[,1]))) {
        analysis_id_col <- "analysis_id"
        cli::cli_alert_info("The first column was used as `analysis_id`.")
      } else {
        cli::cli_abort("No column with the name `analysis_id` or other potential column with unique text identifiers found in the data file.")
      }
    }
  }

  if(import_metadata){
    analysis_metadata_cols <- c(
      "run_seq_num",
      "qc_type",
      "batch_id",
      "is_quantifier",
      "is_istd",
      "precursor_mz",
      "product_mz",
      "collision_energy"
    )
    if(!any(analysis_metadata_cols %in% names(d))) {
      cli::cli_abort("No metadata column names found. Please refer to `??import_data_plain` for more information on support metadata column names.")
    }
  } else {
    analysis_metadata_cols = NULL
  }
  analysis_id_col_sym = rlang::sym(analysis_id_col)

  # Verify all data columns are numeric

  val_cols <- setdiff(names(d), c(analysis_id_col, analysis_metadata_cols))
  analysis_cols <- setdiff(names(d), c(val_cols))

  if(!all(purrr::map_lgl(d[val_cols], is.numeric))) {
    if(import_metadata)
      cli::cli_abort("All columns with feature values must be numeric. Check your data, or if metadata is present, set `import_metadata = TRUE`.")
    else
      cli::cli_abort("All columns with feature values must be numeric. Check your data, and metadata column names, see `??import_data_plain`")
  }

  d_long <- d |>
    tidyr::pivot_longer(cols = -all_of(c(analysis_cols)), names_to = "feature_id", values_to = variable_name) |>
    dplyr::rename(analysis_id = !!analysis_id_col_sym) |>
    dplyr::mutate(across(any_of(c("analysis_id", "batch_id")), as.character)) |>
    dplyr::mutate(integration_qualifier = FALSE)

  d_long
  }





#' Reads a generic wide .csv or .xls sheet with analysis results data
#'
#' @description
#' Imports a tables with the column being the analysis identifier (e.g. sample name) and the following colums being measured features. The values in the table
#'
#'
#' @param path File name and path of a plain wide-format CSV or XLS file
#' @param value_type Parameter type of the feature, i.e. "area", "height", "intensity", "norm_area", "norm_height", "norm_intensity", "feature_conc", "rt", "fwhm" or "width"
#' @param sheet Sheet name in case an Excel file (.xls, .xlsx, .xlsm) is imported
#' @param silent Suppress messages

#' @return A tibble in the long format
#' @noRd
read_data_table <- function(path, value_type = c("area", "height", "intensity", "norm_area", "norm_height", "norm_intensity", "feature_conc", "rt", "fwhm", "width"), sheet = "", silent = FALSE) {
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
  ext <- fs::path_ext(path)

  if (ext == "csv") {
    d <- readr::read_csv(path, col_names = TRUE, trim_ws = TRUE, progress = FALSE,
                         na = c("n/a", "N/A", "NA", "na", "ND", "N.D.", "n.d."),
                         col_types = "cn", locale = readr::locale(encoding = "UTF-8"))
  } else if (ext == "xls" || ext == "xlsx" || ext == "xlm" || ext == "xlmx") {
    if (sheet == "") cli::cli_abort("Please define sheet name via the `sheet` parameter")
    # nms <- names(readxl::read_excel(path = path, sheet = sheet, trim_ws = TRUE, progress = TRUE, na = c("n/a", "N/A", "NA"), n_max = 0))
    # ct <- c("text", rep("numeric",length(nms)-1))
    # d <- readxl::read_excel(path = path, sheet = sheet, trim_ws = TRUE, progress = TRUE, na = c("n/a", "N/A", "NA"), col_types = ct)
    d <- .read_generic_excel(path, sheet)
  } else {
    cli::cli_abort("Invalid file format. Only CSV, XLS and XLSX supported.")
  }

  d |>
    dplyr::mutate(run_seq_num = dplyr::row_number(), .before = 1) |>
    tidyr::pivot_longer(cols = -1:-2, names_to = "feature_id", values_to = value_type) |>
    dplyr::rename(analysis_id = 2) |>
    mutate({{ value_type }} := as.numeric(!!var_value_type))
}







#' @title internal method to read excel sheets
#' @param path excel file name
#' @return tibble table
#' @noRd
.read_generic_excel <- function(path, sheetname) {
  data <- openxlsx2::read_xlsx(file = path,
                              sheet = sheetname,
                              col_names = TRUE,
                              skip_empty_rows = TRUE,
                              skip_empty_cols = TRUE,
                              detect_dates = TRUE,
                              na.strings = c("n/a", "N/A", "NA"),
                              rows = 1)|>
    mutate(across(where(is.character), str_trim))

  nms <- names(data)

  n_col <- length(nms) - 1
  c_numeric <- rep("numeric", n_col)
  ct <- c("character", c_numeric)
  data <- openxlsx2::read_xlsx(file = path,
                              sheet = sheetname,
                              col_names = TRUE,
                              detect_dates = TRUE,
                              skip_empty_rows = TRUE,
                              skip_empty_cols = TRUE,
                              na.strings = c("n/a", "N/A", "NA")
                              #types = ct
                              ) |>
    mutate(across(where(is.character), str_trim))
  data
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



# # for MRMkit wide (previous) result format
# parse_mrmkit_result_wide <- function(path, use_normalized_data = FALSE, silent = FALSE) {
#
#   ext_file <- tolower(fs::path_ext(path))
#
#   sep <- case_when(
#     ext_file == "csv" ~ "csv",
#     ext_file == "tsv" ~ "\t",
#     TRUE ~ NA_character_
#   )
#
#   if (is.na(sep)) cli::cli_abort("Data file type/extension not supported. Please re-export results from MRMkit.")
#
#   d_mrmkit_raw <- readr::read_delim(path,
#                                     delim = sep, col_types = readr::cols(.default = "c"),
#                                     col_names = TRUE, trim_ws = TRUE
#   )
#
#
#   # Extract MRMkit's "QC" info
#   d_mrmkit_featureinfo <- d_mrmkit_raw |>
#     dplyr::filter(.data$name %in% c("Q1", "Q3", "RT", "D-ratio")) |>
#     tidyr::pivot_longer(-.data$name, names_to = "feature_id", values_to = "value") |>
#     tidyr::pivot_wider(names_from = "name", values_from = "value") |>
#     dplyr::mutate(feature_id = dplyr::if_else(stringr::str_detect(.data$feature_id, "RT"), stringr::str_squish(stringr::str_extract(.data$feature_id, ".*(?= RT)")), .data$feature_id)) |>
#     dplyr::rename(precursor_mz = .data$Q1, product_mz = .data$Q1) |>
#     dplyr::mutate(dplyr::across(dplyr::any_of(c("precursor_mz", "product_mz", "RT")), as.numeric))
#
#
#   d_mrmkit_data <- d_mrmkit_raw |>
#     dplyr::filter(!.data$name %in% c("Q1", "Q3", "RT", "D-ratio")) |>
#     dplyr::mutate(file_run_seq_num = dplyr::row_number(), .before = "name") |>
#     tidyr::pivot_longer(-.data$file_run_seq_num:-.data$name, names_to = "feature_id", values_to = "value") |>
#     dplyr::rename(raw_data_filename = .data$name) |>
#     dplyr::mutate(raw_data_filename = stringr::str_remove(.data$raw_data_filename, stringr::regex("\\.mzML$|\\.d$|\\.raw$|\\.wiff$|\\.wiff2$|\\.lcd$", ignore_case = TRUE))) |>
#     dplyr::mutate(raw_data_filename = stringr::str_squish(.data$raw_data_filename)) |>
#     dplyr::mutate(feature_id = dplyr::if_else(stringr::str_detect(.data$feature_id, "RT"), stringr::str_squish(stringr::str_extract(.data$feature_id, ".*(?= RT)")), .data$feature_id)) |>
#     dplyr::mutate(dplyr::across(.data$value, as.numeric)) |>
#     dplyr::left_join(d_mrmkit_featureinfo, by = "feature_id") |>
#     dplyr::relocate(.data$value, .after = dplyr::last_col()) |>
#     dplyr::mutate(dplyr::across(tidyselect::where(is.character), stringr::str_squish))
#
#   if (use_normalized_data) {
#     d_mrmkit_data <- d_mrmkit_data |>
#       filter(str_detect(.data$raw_data_filename, "^norm_")) |>
#       dplyr::rename(norm_intensity = "value")
#   } else {
#     d_mrmkit_data <- d_mrmkit_data |>
#       filter(!str_detect(.data$raw_data_filename, "^norm_")) |>
#       dplyr::rename(feature_area = "value")
#   }
#
#   names(d_mrmkit_data) <- tolower(names(d_mrmkit_data))
#
#   # sep <- if_else(fs::path_ext(final_results_file) == "csv", ",", "\t")
#
#   # if (final_results_file !=""){
#   #   d_results <- readr::read_delim(final_results_file, delim = sep, col_names = TRUE, trim_ws = TRUE, progress = TRUE, show_col_types = FALSE) |>
#   #     dplyr::filter(!is.na(.data$type))|>
#   #     dplyr::mutate(path = stringr::str_squish(str_remove(.data$path, stringr::regex("\\.mzML$|\\.d$|\\.raw$|\\.wiff$|\\.wiff2$|\\.lcd$", ignore_case = TRUE)))) |>
#   #     dplyr::rename(raw_data_filename = .data$path) |>
#   #     #dplyr::mutate(run_seq_num = row_number(), .before = ANALYSIS_ID) |>
#   #     dplyr::select(-tidyselect::any_of(c("batch", "type"))) |>
#   #     tidyr::pivot_longer(-.data$raw_data_filename, names_to = "feature_id", values_to = "feature_norm_intensity")
#   #
#   #   data@dataset <- data@dataset_orig |> left_join(d_results, by=c("raw_data_filename", "feature_id"))
#   #
#   # }
#
#   if (!silent) {
#     cli_alert_success(cli::col_green("Imported {length(unique(d_mrmkit_data$raw_data_filename))} analyses with {length(unique(d_mrmkit_data$feature_id))} features \n"))
#   }
#
#   d_mrmkit_data
# }
