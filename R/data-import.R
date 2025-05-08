
#' @title Import Agilent MassHunter Quantitative Analysis CSV files
#' @description Imports .csv files exported from Agilent MassHunter Quantitative
#' Analysis software, containing peak integration results. The input files must
#' have anlyses (samples) in rows, features/compounds in columns, and either
#' peak areas, peak heights, or response as the values. Additional columns, such
#' as retention time (RT), full-width at half-maximum (FWHM), precursor m/z
#' (PrecursorMZ), and collision energy (CE), will also be imported and made
#' available in the `MidarExperiment` object for downstream analyses.
#'
#' When a directory path is provided, all matching .csv files in that directory
#' will be imported and merged into a single dataset. This is useful when
#' importing datasets that were pre-processed in blocks, resulting in multiple
#' files. Each unique combination of feature and raw data file must only occur
#' once across all source data files. Duplicate combinations will result in an
#' error.
#'
#' @param data MidarExperiment object
#' @param path One or more file paths, or a directory path (in which case all
#'   matching files will be imported)
#' @param import_metadata Logical, whether to extract and add metadata from the
#'   analysis result file
#' @param expand_qualifier_names Logical, whether to add the quantifier name in
#'   front of the qualifier name (the latter only has the m/z transition values)
#' @param conc_column Which concentration field of the masshunter data to use, in
#' case "Calc. Conc." and "Final. Conc." are present.  Default is "conc_final".
#' Must be one of "conc_calc" or "conc_final" (default).
#' @param silent Logical, whether to suppress most notifications
#' @return MidarExperiment object with the imported data
#' @examples
#' mexp <- MidarExperiment()

#' file_path = system.file("extdata", "MHQuant_demo.csv", package = "midar")
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

import_data_masshunter <- function(data = NULL,
                                   path,
                                   import_metadata = TRUE,
                                   expand_qualifier_names = TRUE,
                                   conc_column = "conc_final",
                                   silent = FALSE) {
  check_data(data)
  rlang::arg_match(conc_column, c("conc_calc", "conc_final"))
  #if (fs::path_ext(path) == "csv") {
  data <- import_data_main(data, path, "parse_masshunter_csv", "*.csv", expand_qualifier_names = expand_qualifier_names, silent = silent, conc_column = conc_column)
  data <- set_intensity_var(data, variable_name = NULL, auto_select = TRUE, warnings = TRUE, "feature_area", "feature_response", "feature_height")
  if (import_metadata) data <- import_metadata_from_data(data, qc_type_column_name = "sample_type")
  data

}

#' @title Import Skyline Peak Integration Results
#'
#' @description
#' This function imports tabular data files (*.csv) exported from `Skyline`,
#' containing peak integration results.
#'
#' @details
#' In Skyline, transitions are defined by
#' the `Molecule Name` and corresponding precursor and product m/z values,
#' rather than any identifier. When importing this data, the `feature_id` is
#' generated using the `Molecule Name` with either precursor/product names or
#' m/z values, unless `Molecule Name` uniquely identifies the features (refer to
#' the `transition_id_columns` argument below).
#'
#' The following supported columns from Skyline can be imported:
#'
#' | Skyline Column Name    | MiDAR Column Name      | Required  |
#' |------------------------|------------------------|-----------|
#' | `Replicate Name`       | `analysis_id`          | Yes       |
#' | `Molecule List Name`   | `feature_class`        | No        |
#' | `Molecule Name`        | `feature_id`           | Yes       |
#' | `Precursor Name`       | `feature_id`           | Yes*      |
#' | `Product Name`         | `feature_id`           | Yes*      |
#' | `Precursor Mz`         | `method_precursor_mz`  | Yes*      |
#' | `Product Mz`           | `method_product_mz`    | Yes*      |
#' | `Area`                 | `feature_area`         | Yes       |
#' | `Retention Time`       | `feature_rt`           | No        |
#'
#' *Requirements for these columns are described in `transition_id_columns`.
#'
#' To export results from Skyline, use the 'Molecule Transition Results' format
#' and include `Replicate Name`, `Molecule Name`, and either
#' `Precursor Mz`/`Product Mz` or `Precursor Name`/`Product Name` columns.
#' At least one feature variable, such as `Area` or `RT`, must also be exported.
#'
#' @param data A `MidarExperiment` object.
#' @param path One or more file paths, or a directory path from which all
#' matching files will be imported.
#' @param transition_id_columns A character vector specifying the columns that
#' define the transition (precursor and product) to use for unique `feature_id`
#' generation. Options are "name", "mz", or "none". If "none", `feature_id` is
#' derived from `Molecule Name` and `Precursor Name` or `Product Name`. If "mz",
#' `feature_id` is based on `Precursor Mz` and `Product Mz`. Using "none" will
#' result in `feature_id` being a copy of `Molecule Name`, with an error raised
#' if it is not unique for each transition.
#' @param import_metadata Logical; whether to import additional metadata columns
#' (e.g., precursor/product m/z values).
#' @param silent Logical; whether to suppress most notifications.
#' @return A `MidarExperiment` object containing the imported data.
#' @examples
#' mexp <- MidarExperiment()
#' file_path <- system.file("extdata", "Skyline_MoleculeTransitionResults.csv", package = "midar")
#' mexp <- import_data_skyline(
#'   data = mexp,
#'   path = file_path,
#'   transition_id_columns = "mz",
#'   import_metadata = TRUE
#' )
#' print(mexp)
#' @export
import_data_skyline <- function(data = NULL, path,
                                transition_id_columns = c("name", "mz", "none"),
                                import_metadata = TRUE, silent = FALSE) {
  check_data(data)
  rlang::arg_match(transition_id_columns, c("name", "mz", "none"))
  data <- import_data_main(data = data, path = path, import_function = "parse_skyline_result", file_ext = "*.tsv|*.csv", silent = FALSE, transition_id_columns = transition_id_columns)
  data <- set_intensity_var(data, variable_name = NULL, auto_select = TRUE, warnings = TRUE, "feature_area", "feature_height")

  if (import_metadata) data <- import_metadata_from_data(data, qc_type_column_name = "sample_type")
  data
}


#' @title Import MRMkit peak integration results
#' @description
#' Imports tabular data files (*.tsv) generated from `MRMkit` containing peak
#' integration results. The input files must be in a long format with columns
#' for the raw data file name, feature ID, peak intensity, and other arguments
#' Additional information, such as retention time, FWHM, precursor/product m/z,
#' and CE will also be imported and made available in the `MidarExperiment`
#' object for downstream analyses. Concentrations will also be imported if present,
#' whereby when both `Calc. Conc.` and `Final Conc.` are present in the files,
#' the argument `conc_column` can be used to specify which concentration field as
#' 'concentration' to use in downstream analyses.
#'
#' When a directory path is provided, all matching files in that directory will
#' be imported and merged into a single dataset. This is useful when importing
#' datasets that were pre-processed in blocks, resulting in multiple files. Each
#' unique combination of feature and raw data file must only occur once across
#' all source data files. Duplicate combinations will result in an error.
#'
#' @param data MidarExperiment object
#' @param path One or more file paths, or a directory path (in which case all matching files will be imported)
#' @param import_metadata Logical, whether to import additional metadata columns (e.g., `batch_id`, `qc_type`)
#' @param silent Logical, whether to suppress most notifications
#' @return MidarExperiment object with the imported data
#' @examples
#' mexp <- MidarExperiment()
#'
#' file_path = system.file("extdata", "MRMkit_demo.tsv", package = "midar")
#'
#' mexp <- import_data_mrmkit(
#'   data = mexp,
#'   path = file_path,
#'   import_metadata = TRUE)
#
#' print(mexp)
#' @export
import_data_mrmkit <- function(data = NULL, path, import_metadata = TRUE, silent = FALSE) {
  check_data(data)
  data <- import_data_main(data = data, path = path, import_function = "parse_mrmkit_result", file_ext = "*.tsv|*.csv", silent = FALSE)
  data <- set_intensity_var(data, variable_name = NULL, auto_select = TRUE, warnings = TRUE, "feature_area", "feature_height")

  if (import_metadata) data <- import_metadata_from_data(data, qc_type_column_name = "qc_type")
  data
}

#' (Depreciated) Import Wide CSV Files
#' @details
#' This function is deprecated. Please use [import_data_csv_wide()] instead.

#'
#' @param data MidarExperiment object
#' @param path One or more file names with path, or a folder path, which case all *.csv files in this folder will be read.
#' @param variable_name Variable type representing the values in the table. Must be one of "intensity", "norm_intensity", "conc", "area", "height", "response")
#' @param analysis_id_col Column to be used as analysis_id. `NA` (default) used 'analysis_id' if present, or the first column if it contains unique values.
#' @param import_metadata Import additional metadata columns (e.g. batch ID, sample type) and add to the `MidarExperiment` object.
#' Only following metadata column names are supported: "qc_type", "batch_id", "is_quantifier", "is_istd", "analysis_order"
#' @param first_feature_column Column number of the first column representing the feature values
#' @param na_strings A character vector of strings which are to be interpreted as NA values. Blank fields are also considered to be missing values.
# #' @param silent Su ppress notifications
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

import_data_csv <- function(data = NULL, path, variable_name, analysis_id_col = NA, import_metadata  = TRUE, first_feature_column = NA, na_strings = "NA") {
  cli::cli_alert_warning(col_yellow("The function import_data_csv is deprecated. Please use import_data_csv_wide instead."))
  import_data_csv_wide(data = data, path = path, variable_name = variable_name, analysis_id_col = analysis_id_col, import_metadata = import_metadata, first_feature_column = first_feature_column, na_strings = na_strings)
}

#' Import Analysis Results from Plain Wide-Format CSV Files
#'
#' @description
#' Imports analysis result data from wide-format `.csv` files, where each row corresponds
#' to a unique analysis-feature pair and columns contain analysis- or feature-specific variables.
#'
#' @details
#' The dataset must include two identifier columns: `"analysis_id"` and `"feature_id"`,
#' with each pair of values unique across the table. Additionally, the table must contain
#' at least one feature variable column, such as `"area"`, `"height"`, `"intensity"`,
#' `"norm_intensity"`, `"response"`, `"conc"`, `"rt"`, or `"fwhm"`. Some downstream
#' functions may require specific columns among these to be present.
#'
#' The `variable_name` argument specifies the data type represented in the table, which
#' must be one of: `"area"`, `"height"`, `"intensity"`, `"norm_intensity"`, `"response"`,
#' `"conc"`, `"conc_raw"`, `"rt"`, or `"fwhm"`.
#'
#' If there is no column named `analysis_id`, it will be inferred from the first column,
#' provided it contains unique values.
#'
#' When `import_metadata` is set to `TRUE`, the following metadata columns will be imported
#' if present:
#' \itemize{
#'   \item \code{analysis_order}
#'   \item \code{qc_type}
#'   \item \code{batch_id}
#'   \item \code{is_quantifier}
#' }
#'
#' To prevent additional non-metadata columns from being misinterpreted as features,
#' use the `first_feature_column` parameter to specify the column where feature data starts.
#'
#' If a directory path is provided to `path`, all `.csv` files in that directory will be
#' processed and merged into a single dataset. This facilitates handling datasets split
#' into multiple files during preprocessing. Ensure each feature and raw data file pair
#' appears only once to avoid duplication errors.
#'
#' The `na_strings` parameter allows specifying character strings to be interpreted as
#' missing values (NA). Blank fields are also treated as missing.
#'
#' @param data A \code{MidarExperiment} object.
#' @param path A file path or vector of file paths, or a directory path. If a directory is
#'   provided, all `.csv` files within it will be read.
#' @param variable_name A character string specifying the variable type contained in the
#'   data. Must be one of `"intensity"`, `"norm_intensity"`, `"conc"`, `"area"`, `"height"`, or `"response"`.
#' @param analysis_id_col The column name or index to be used as `analysis_id`. Defaults to `NA`,
#'   in which case `"analysis_id"` is used if present; otherwise, the first column is used if it contains unique values.
#' @param import_metadata Logical indicating whether to import additional metadata columns
#'   (e.g., batch ID, sample type) into the `MidarExperiment` object. Supported metadata columns are:
#'   `"qc_type"`, `"batch_id"`, `"is_quantifier"`, `"is_istd"`, and `"analysis_order"`.
#' @param first_feature_column Integer indicating the column number where feature value columns start.
#' @param na_strings Character vector of strings to interpret as NA values. Blank fields are also treated as NA.
#'
#' @return A \code{MidarExperiment} object containing the imported dataset.
#'
#' @examples
#' file_path <- system.file("extdata", "plain_wide_dataset.csv", package = "midar")
#' mexp <- MidarExperiment()
#' mexp <- import_data_csv_wide(
#'   data = mexp,
#'   path = file_path,
#'   variable_name = "conc",
#'   import_metadata = TRUE
#' )
#' print(mexp)
#'
#' @export

import_data_csv_wide <- function(data = NULL, path, variable_name, analysis_id_col = NA, import_metadata  = TRUE, first_feature_column = NA, na_strings = "NA") {
  check_data(data)
  data <- import_data_main(data = data, path = path, import_function = "parse_plain_wide_csv", file_ext = "*.csv", silent = FALSE, variable_name = variable_name, analysis_id_col = analysis_id_col, import_metadata = import_metadata, first_feature_column =first_feature_column, na_strings = na_strings)
  data <- set_intensity_var(data, variable_name = paste0("feature_", str_remove(variable_name, "feature_")), auto_select = FALSE, warnings = FALSE, "feature_area", "feature_height", "feature_conc")

  if (import_metadata)
    data <- import_metadata_from_data(data, qc_type_column_name = "qc_type")
  else{
    #set_analysis_order_analysismetadata(data, order_by = "default")
    data <- link_data_metadata(data)
  }

  data
}

#' Import Analysis Results from Long Format CSV Files
#'
#' @description
#' This function imports analysis results from CSV files in long table format,
#' where each row represents a unique observation of a feature-value pair for
#' an analysis (sample), along with associated feature variables and other
#' metadata. See "Details" below for more information on using this function.
#'
#' @details
#' When no column mapping is provided via the `column_mapping` argument, the
#' function will automatically detect and import columns with the following
#' names:
#'
#' | CSV Column Name      | MiDAR Internal Name   | Required  |
#' |----------------------|-----------------------|-----------|
#' | `analysis_id`        | `analysis_id`         | Yes       |
#' | `feature_id`         | `feature_id`          | Yes       |
#' | `qc_type`            | `qc_type`             | No        |
#' | `sample_id`          | `sample_id`           | No        |
#' | `batch_id`           | `batch_id`            | No        |
#' | `istd_feature_id`    | `istd_feature_id`     | No        |
#' | `feature_class`      | `feature_class`       | No        |
#' | `analyte_id`         | `analyte_id`          | No        |
#' | `precursor_mz`       | `method_precursor_mz` | No        |
#' | `product_mz`         | `method_product_mz`   | No        |
#' | `area`               | `feature_area`        | No        |
#' | `height`             | `feature_height`      | No        |
#' | `intensity`          | `feature_intensity`   | No        |
#' | `rt`                 | `feature_rt`          | No        |
#' | `fwhm`               | `feature_fwhm`        | No        |
#' | `width`              | `feature_width`       | No        |
#'
#' Detection of these columns is case-insensitive. Additionally, if feature
#' variable columns use the internal naming convention with prefixes
#' "feature_" or "method_" (e.g. `feature_area` instead of `area`), the
#' function will detect and import them automatically.
#'
#' To import data with different column names, provide a named vector mapping
#' CSV column names to the internal column names used by `MidarExperiment`.
#' The mapping should be in the format:
#' `c("analysis_id" = "[CSV column name for analysis]", "feature_id" = "[CSV
#' column name for feature]", ...)`, where the right-hand side refers to the
#' exact column name in the CSV file header. Columns matching internal names
#' do not require mapping and will be imported automatically. The mapping is
#' case-insensitive.
#'
#' Note that the dataset must contain an analysis identifier, either as an
#' `analysis_id` column or via a mapped column.
#'
#' The function processes all CSV files in the specified directory or the
#' given file(s), combining them into a single dataset. This supports datasets
#' split across multiple files during preprocessing. Each feature and raw data
#' file pair should appear only once to avoid duplication.
#'
#' The `na_strings` parameter allows specifying character strings that should
#' be interpreted as `NA`, ensuring proper handling of missing values.
#'
#' @param data A `MidarExperiment` object to which the imported data will be
#'   added.
#' @param path A single file path, multiple file paths, or a directory path. If
#'   a directory is provided, all `*.csv` files within will be imported.
#' @param column_mapping A named character vector mapping internal column names
#'   to CSV column names. Should include keys such as `"analysis_id"`,
#'   `"feature_id"`, and feature variable names. If `NULL` (default), the
#'   function attempts automatic detection.
#' @param import_metadata Logical indicating whether to import additional
#'   metadata columns (e.g., batch ID, sample type) into the
#'   `MidarExperiment` object. Supported metadata column names include
#'   `"qc_type"`, `"batch_id"`, `"is_quantifier"`, `"is_istd"`, and
#'   `"analysis_order"`.
#' @param na_strings Character vector of strings to interpret as missing values
#'   (`NA`). Blank fields are always treated as missing.
#' @param warn_unrecognized_columns Logical indicating whether to issue a
#'   warning when unknown columns are encountered in the dataset.
#' @param silent Logical indicating whether to suppress most notifications and
#'   messages.
#' @return A `MidarExperiment` object containing the imported data.
#'
#' @examples
#' file_path <- system.file("extdata", "plain_long_dataset.csv", package = "midar")
#' mexp <- MidarExperiment()
#'
#' # Define the column mapping; right side is the CSV column name
#' col_map <- c(
#'   "analysis_id" = "raw_data_filename",
#'   "qc_type" = "qc_type",
#'   "feature_id" = "feature_id",
#'   "feature_class" = "feature_class",
#'   "istd_feature_id" = "istd_feature_id",
#'   "feature_rt" = "rt",
#'   "feature_area" = "area"
#' )
#'
#' mexp <- import_data_csv_long(
#'   data = mexp,
#'   path = file_path,
#'   column_mapping = col_map,
#'   import_metadata = TRUE
#' )
#'
#' print(mexp)
#'
#' @export

import_data_csv_long <- function(data = NULL, path, import_metadata = TRUE, column_mapping = NULL, na_strings = "NA", warn_unrecognized_columns = TRUE, silent = FALSE) {
  check_data(data)
  data <- import_data_main(data = data, path = path, import_function = "parse_plain_long_csv", file_ext = "*.csv", silent = silent, column_mapping = column_mapping, na_strings = na_strings, warn_unrecognized_columns = warn_unrecognized_columns)
  data <- set_intensity_var(data, variable_name = NULL, auto_select = TRUE, warnings = FALSE, "feature_area", "feature_height", "feature_conc")

  if (import_metadata)
    data <- import_metadata_from_data(data, qc_type_column_name = "qc_type")
  else {
    data <- link_data_metadata(data)
  }

  data
}


import_data_main <- function(data = NULL, path, import_function, file_ext, na_strings, silent, ...) {
  check_data(data)
  if (all(!fs::is_dir(path))) {
    file_paths <- fs::path_tidy(path)
  } else {
    file_paths <- fs::dir_ls(path, glob = file_ext)
  }

  if (!all(fs::file_exists(file_paths))) cli::cli_abort(col_red("One or more given files do not exist. Please verify file paths."))
  if (any(duplicated(file_paths))) cli::cli_abort(col_red("One or more given files are duplicated. Please verify file paths."))

  names(file_paths) <- file_paths
  args <- list(...)

  d_raw <- file_paths |>
    purrr::map_dfr(.f = \(x) do.call(what = import_function, args = append(x, args)), .id = "data_source") |>
    relocate("analysis_id", "data_source")


  if(!"analysis_order" %in% names(d_raw)) {
    d_runorder <- d_raw |>
      select("analysis_id") |>
      distinct() |>
      mutate(analysis_order = row_number())

    d_raw <- d_raw |>
      left_join(d_runorder, by = "analysis_id") |>
      relocate("analysis_order", .before = 1)
  }



  # VERIFY DATA, i.e. analysis_ids, feature_ids, and values are replicated ===
  ## which can be result of multiple imports of the same/overlapping data or due to parsing error
  n_idpairs_distinct <- d_raw |>
    select("analysis_id", "feature_id") |>
    distinct(.keep_all = FALSE) |>
    nrow()
  if (nrow(d_raw) > n_idpairs_distinct) {
    has_duplicated_id <- TRUE

    n_idvalpairs_distinct <- d_raw |>
      select("analysis_order", "analysis_id", "feature_id",
        any_of(c("method_precursor_mz", "method_product_mz", "feature_class", "method_collision_energy", "method_polarity", "method_conc_expected",
          "feature_area", "feature_rt", "feature_fwhm", "feature_intensity", "feature_height", "feature_conc", "feauture_norm_intensity"))) |>
      distinct(.keep_all = FALSE) |>
      nrow()

    if (n_idvalpairs_distinct == n_idpairs_distinct ) {
      has_duplicated_values <- TRUE
    } else {
      has_duplicated_values <- FALSE
    }
  } else {
    has_duplicated_id <- FALSE
  }

  if (has_duplicated_id) {
    if (has_duplicated_values) {
      cli::cli_abort(glue::glue("Imported data contains duplicated reportings (analysis and feature pairs) with {cli::style_italic('identical')} feature variable values. Please verify imported dataset(s)."))
    } else {
      cli::cli_abort(glue::glue("Imported data contains duplicated reportings (analysis and feature pairs) with {cli::style_italic('different')} feature variable values. Please verify imported dataset(s)."))
    }
  }

  data@dataset_orig <- dplyr::bind_rows(pkg.env$table_templates$dataset_orig_template, d_raw)

  # TODO: excl_unmatched_analyses below

  #check_integrity_analyses(data, excl_unmatched_analyses = TRUE, silent = TRUE)
  # stopifnot(methods::validObject(data))

  if (!silent) {
    if (!any(data@dataset_orig$integration_qualifier)) {
      cli_alert_success(cli::col_green("Imported {length(unique(data@dataset_orig$analysis_id))} analyses with {length(unique(data@dataset_orig$feature_id))} features"))
    } else if (any(data@dataset_orig$integration_qualifier)) {
      cli_alert_success(cli::col_green("Imported {length(unique(data@dataset_orig$analysis_id))} analyses with {length(unique(data@dataset_orig$feature_id))} features ({length(unique(data@dataset_orig$feature_id[!data@dataset_orig$integration_qualifier]))} quantifiers, {length(unique(data@dataset_orig$feature_id[!data@dataset_orig$integration_qualifier]))} qualifiers)"))
    }
  }

  data@feature_intensity_var = ""
  data@status_processing <- "Raw data imported"

  data
}



#' Reads and parses one Agilent MassHunter Quant CSV result file
#'
#' @param path File path of MassHunter Quant CSV file
#' @param silent Suppress messages
#' @param expand_qualifier_names If TRUE, original qualifier names will be renamed by adding the quantifier name in front and placing qualifier name into square brackets(e.g. `Qualifier (422.3 -> 113.0)` transition names of quantifier will be added to qualifier names
#' @param conc_column Which concentration field of the masshunter data to use, in
#' case "Calc. Conc." and "Final. Conc." are present.  Default is "conc_final".
#' @return A tibble with the parse results in the long format
#' @examples
#' file_path = system.file("extdata", "MHQuant_demo.csv", package = "midar")
#'
#' tbl <- parse_masshunter_csv(
#'   path = file_path,
#'   expand_qualifier_names = TRUE)
#'
#' head(tbl)
#' @export

parse_masshunter_csv <- function(path, expand_qualifier_names = TRUE, silent = FALSE, conc_column = "conc_final") {

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

  if (!any(c("Data File", "Name", "Area", "RT", "Calc. Conc.", "Final Conc", "Height", "Resp") %in% unlist(datWide[3,]))) cli::cli_abort("Data file is in an unsupported or corrupted format. Please try re-export your data in MH with compounds as columns.")

  if ("Compound Method" == datWide[[1,1]]) cli::cli_abort(col_red("Compound table format is currently not supported. Please re-export your data in MH with compounds as columns."))

  if (!any(c("Data File") %in% unlist(datWide[3,]))) cli::cli_abort(col_red("'Data File' column is required and used as a unique identifier, but is missing or the file is in an unsupported/corrupt format. Please re-export from MassHunter with 'Data File' included, samples in rows, and features in columns."))

  if(datWide[1, 1] != "Sample") datWide[1, 1] <- "Sample"

  feature_id_tbl <- dplyr::tibble(`_prefixXXX_` = datWide[1, ] |> unlist() |> dplyr::na_if(""))

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

  datWide[1, ] <- dplyr::tibble(`_prefixXXX_` = datWide[1, ] |> unlist() |> dplyr::na_if("")) |>
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
  datWide <- datWide |> dplyr::select(-where(~ .x[2] == ""))

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
    dplyr::mutate(file_analysis_order = dplyr::row_number(), .before = 1) |>
    tidyr::pivot_longer(cols = all_of(param_transition_names), names_pattern = "(.*)\t(.*)$", names_to = c("Param", "feature_id")) |>
    tidyr::pivot_wider(names_from = "Param", values_from = "value")


  # Convert types of known parameters and fields in the data set

  new_numeric_colnames <- c(
    "feature_rt" = "Results_RT",
    "feature_response" = "Results_Resp",
    "feature_height" = "Results_Height",
    "feature_area" = "Results_Area",
    "feature_accuracy" = "Results_Accuracy",
    "feature_conc_calc" = "Results_CalcConc",
    "feature_conc_final" = "Results_Final Conc",
    "feature_ionratio" = "Results_Ratio",
    "feature_fwhm" = "Results_FWHM",
    "feature_width" = "Results_Width",
    "feature_sn_ratio" = "Results_SN",
    "feature_int_start" = "Results_IntStart",
    "feature_int_end" = "Results_IntEnd",
    "feature_symetry" = "Results_Symmetry",
    "feature_istd_responseratio" = "Results_ISTD Resp. Ratio",
    "method_target_rt" = "Method_RT",
    "method_conc_expected" = "Method_ExpConc",
    "method_precursor_mz" = "Method_Precursor Ion",
    "method_product_mz" = "Method_Product Ion",
    "method_collision_energy" = "Method_Collision Energy",
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
    dplyr::mutate(dplyr::across(where(is.character), stringr::str_squish)) |>
    dplyr::relocate(any_of("sample_name"), .after = "raw_data_filename")


  if("feature_conc_calc" %in% names(datLong) || "feature_conc_final" %in% names(datLong)){
    if("feature_conc_calc" %in% names(datLong) && "feature_conc_final" %in% names(datLong)){
      if(conc_column == "conc_calc" && "feature_conc_calc" %in% names(datLong)){
        datLong <- datLong |>
          dplyr::mutate(feature_conc = .data$feature_conc_calc)
      } else if(conc_column == "conc_final" & "feature_conc_final" %in% names(datLong)){
        datLong <- datLong |>
          dplyr::mutate(feature_conc = .data$feature_conc_final)
      }
    } else if("feature_conc_calc" %in% names(datLong)){
      datLong <- datLong |>
        dplyr::mutate(feature_conc = .data$feature_conc_calc)
    } else if("feature_conc_final" %in% names(datLong)){
      datLong <- datLong |>
        dplyr::mutate(feature_conc = .data$feature_conc_final)
    }
  }

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


#' Parses skyline peak integration results into a tibble
#'
#' @param path File name of the MRMkit result file (*.tsv or *.csv)
# #' @param use_normalized_data Import raw peak areas or normalized peak areas from the file
#' @param na_strings A character vector of strings to be interpreted as NA values.
#' @param silent No comments printed
#' @param ... Additional arguments passed to the function. Currently only `transition_id_columns`, which is used for
#' Skyline-like data files. The column mapped to `feature_id` will be appended with following info: if  `transition_id_columns = "name"`, the function will use the `method_precursor_name` and `method_product_name` columns to create unique feature IDs.
#' If `transition_id_columns = "mz"`, the function will use the `method_precursor_mz` and `method_product_mz` columns to create unique feature IDs.

#' @return A tibble in the long format
#' @examples
#'
#' file_path = system.file("extdata", "MRMkit_demo.tsv", package = "midar")
#'
#' tbl <- parse_mrmkit_result(path = file_path)
#'
#' head(tbl)
#' @export

# TODO: remove support for norm intensity
parse_skyline_result <- function(path, na_strings, silent = FALSE, ...) {


  col_map <- c(
    "analysis_id"= "Replicate Name",
    "feature_id" = "Molecule Name",
    "feature_class" = "Molecule List Name",
    "feature_rt" = "Retention Time",
    "feature_area" = "Area",
    "method_precursor_name" = "Precursor Name",
    "method_product_name" = "Product Name",
    "method_precursor_mz" = "Precursor Mz",
    "method_product_mz" = "Product Mz"
  )

  parse_plain_long_csv(path = path,
                       na_strings  = na_strings,
                       silent = silent,
                       column_mapping = col_map,
                       warn_unrecognized_columns = FALSE,
                       ...)
}

#' Parses MRMkit peak integration results into a tibble
#'
#' @param path File name of the MRMkit result file (*.tsv or *.csv)
# #' @param use_normalized_data Import raw peak areas or normalized peak areas from the file
#' @param silent No comments printed
#' @return A tibble in the long format
#' @examples
#'
#' file_path = system.file("extdata", "MRMkit_demo.tsv", package = "midar")
#'
#' tbl <- parse_mrmkit_result(path = file_path)
#'
#' head(tbl)
#' @export

# TODO: remove support for norm intensity
parse_mrmkit_result <- function(path, silent = FALSE) {

  col_map <- c(
      "raw_data_filename"= "raw_data_filename",
      "acquisition_time_stamp" = "time_stamp",
      "qc_type" = "sample_type",
      "batch_id" = "batch",
      "feature_id" = "feature_name",
      "istd_feature_id" = "internal_standard",
      "feature_rt" = "rt_apex",
      "feature_area" = "area",
      #"feature_norm_intensity" = "area_normalized",
      "feature_height" = "height",
      "feature_fwhm" = "FWHM",
      #"feature_width" = "width",
      "feature_int_start" = "rt_int_start",
      "feature_int_end" = "rt_int_end",
      "method_polarity" = "polarity",
      "method_precursor_mz" = "precursor_mz",
      "method_product_mz" = "product_mz",
      "method_collision_energy" = "collision_energy"
    )

  d_raw <- parse_plain_long_csv(path = path,
                       column_mapping = col_map,
                       warn_unrecognized_columns = FALSE,
                       silent = silent)

  d_raw <- d_raw |> mutate(feature_width = .data$feature_int_end - .data$feature_int_start)

  d_raw
}


#' Parses a plain long CSV file
#'
#' Parses a CSV table with analysis/samples and feature pairs in rows,
#' columns representing the feature variables.
#'
#' @param path File name (*.tsv or *.csv)
#' @param na_strings A character vector of strings to be interpreted as NA values.
#' Blank fields are also considered to be missing values.
#' @param column_mapping A named character vector with the mapping of the columns in the input file to the columns in the output table.
#' If NULL (default), the function will use the default mapping. See [import_data_csv_long()] for details.
#' @param warn_unrecognized_columns Logical indicating whether to issue a warning for unknown columns in the dataset.
#' @param ... Additional arguments passed to the function. Currently only `transition_id_columns`, which is used for
#' Skyline-like data files. The column mapped to `feature_id` will be appended with following info: if  `transition_id_columns = "name"`, the function will use the `method_precursor_name` and `method_product_name` columns to create unique feature IDs.
#' If `transition_id_columns = "mz"`, the function will use the `method_precursor_mz` and `method_product_mz` columns to create unique feature IDs.
#' @param silent No comments printed
#' @return A tibble in the long format
#' @examples
#'
#' file_path = system.file("extdata", "plain_long_dataset.csv", package = "midar")
#'
#' tbl <- parse_plain_long_csv(path = file_path)
#'
#' head(tbl)
#' @export
parse_plain_long_csv <- function(path, na_strings = "NA", silent = FALSE, column_mapping = NULL, warn_unrecognized_columns = TRUE, ...) {

  args <- list(...)

  ext_file <- tolower(fs::path_ext(path))
  sep <- case_when(
    ext_file == "csv" ~ ",",
    ext_file == "tsv" ~ "\t",
    TRUE ~ NA_character_
  )
  if (is.na(sep)) cli::cli_abort(col_red("Data file type/extension not supported."))

  d_raw <- readr::read_delim(path,
                                    delim = sep, col_types = readr::cols(.default = "c"),
                                    col_names = TRUE, trim_ws = TRUE)

  d_raw <- d_raw |> dplyr::rename_with(tolower)


  if(is.null(column_mapping)){
    column_mapping <- c(
      "analysis_id" = "analysis_id",
      "raw_data_filename" = "raw_data_filename",
      "acquisition_time_stamp" = "time_stamp",
      "qc_type" = "qc_type",
      "batch_id" = "batch_id",
      "feature_id" = "feature_id",
      "feature_class" = "feature_class",
      "istd_feature_id" = "istd_feature_id",
      "qc_type" = "qc_type",
      "feature_rt" = "rt",
      "feature_area" = "area",
      "feature_height" = "height",
      "feature_norm_intensity" = "norm_area",
      "feature_height" = "height",
      "feature_fwhm" = "fwhm",
      "feature_width" = "width",
      "feature_int_start" = "int_start",
      "feature_int_end" = "int_end",
      "method_polarity" = "polarity",
      "method_precursor_mz" = "precursor_mz",
      "method_product_mz" = "product_mz",
      "method_collision_energy" = "collision_energy"
    )

    d_raw <- d_raw |>
      dplyr::rename_with(
        .fn = ~ ifelse(
          str_detect(., "_id|_class|_name"),
          .,
          str_replace_all(., "^(method_|feature_)", "")
        )
      )
  }
  column_mapping <- tolower(column_mapping)
  names(column_mapping) <- tolower(names(column_mapping))

  column_mapping_filt <- column_mapping[column_mapping %in% names(d_raw)]
  d_raw <- d_raw |>
    rename(!!!column_mapping_filt)


  # -------- SKYLINE-Like format -----
  #Check if input is of Skyline format (analysis id defined by molecule name and ion names/mz)
  # TODO: A bit of a workaround... needs to be shifted, wrapped in a separate function called here
  if ("transition_id_columns" %in% names(args)){

    if(!"analysis_id" %in% names(d_raw) ){
      cli_abort(col_red("The `Replicate Name` column is missing in the data file. Please ensure it is included in the Skyline export."))
    }

    if(!"feature_id" %in% names(d_raw) ){
      cli_abort(col_red("The `Molecule Name` column is missing in the data file. Please ensure it is included in the Skyline export."))
    }

    if (args$transition_id_columns == "name") {
      if("method_precursor_name" %in% names(d_raw) &&
         "method_product_name" %in% names(d_raw) &&
         !any(is.na(d_raw$method_precursor_name)) &&
         !any(is.na(d_raw$method_product_name))){
        d_raw <- d_raw |>
          unite("feature_id", c("feature_id", "method_precursor_name", "method_product_name"), sep = "_", remove = FALSE)
      } else {
        cli_abort(col_red(
          "`Precursor Name` and/or `Product Name` columns are missing or contain no/missing values. Ensure these are included in the Skyline export, or modify `transition_id_columns` argument"
        ))
      }
    } else if (args$transition_id_columns == "mz") {
      if("method_precursor_mz" %in% names(d_raw) &&
         "method_product_mz" %in% names(d_raw) &&
         !any(is.na(d_raw$method_precursor_mz)) &&
         !any(is.na(d_raw$method_product_mz)))
         {
        d_raw <- d_raw |>
          unite("feature_id", c("feature_id", "method_precursor_mz", "method_product_mz"), sep = "_", remove = FALSE)
    } else {
      cli_abort(col_red(
        "`Precursor Mz` and/or `Product Mz` columns are missing or contain no/missing values. Ensure these are included in the Skyline export, or modify `transition_id_columns` argument"
      ))
      }
    } else if (args$transition_id_columns == "none") {
      # Do nothing, as Molecule Name was mapped to feature_id
    }

    # Check for duplicate feature IDs
    duplicate_rows <- d_raw |>
      group_by(.data$feature_id, .data$analysis_id) |>
      filter(n() > 1)

    if (nrow(duplicate_rows) > 0) {
      if(args$transition_id_columns == "none")
        cli::cli_abort(col_red("`Molecule Name` is not unique identifier for each transition. To generate unique feature IDs, please set the `transition_id_columns` argument to either `name` or `mz`."))
      else
        cli::cli_abort(col_red("Feature IDs are not unique even with precursor/product ion details added. Please check Precursor/Product Name/Mz columns in the data file."))
    }
  }
  # ---



  # Conditionally add 'analysis_id' if not present but raw_data_filename exists
  if ("raw_data_filename" %in% names(d_raw) && !("analysis_id" %in% names(d_raw))) {
    d_raw <- d_raw |> dplyr::mutate(analysis_id = .data$raw_data_filename)
  }


  # Check required IDs
  missing_ids <- setdiff("analysis_id", names(d_raw))
  if (length(missing_ids) > 0) {
    cli::cli_abort(col_red("Required `analysis_id` column is missing. Please map column names via the `column_mapping` argument or verify data file."))
  }

  d_raw <- d_raw |>
    dplyr::mutate(
      analysis_id = stringr::str_remove(
        .data$analysis_id,
        stringr::regex("\\.mzML$|\\.d$|\\.raw$|\\.wiff$|\\.wiff2$|\\.lcd$", ignore_case = TRUE)
      )
    )

  # 2. Check for unmapped columns

  unmapped <- setdiff(names(d_raw), c("analysis_id", names(column_mapping)))
  if (length(unmapped) > 0) {
    if(warn_unrecognized_columns){
      cli::cli_alert_warning(col_yellow("Following unrecognized columns present in the data and were ignored: {.val {unmapped}}."))
      cli::cli_alert_warning(col_yellow("Use argument `column_mapping` to define column mapping."))
    }
    d_raw <- d_raw |> select(-!!unmapped)
  }



  # 3. Check for presence of at least one quant/intensity column
  intensity_cols <- c(
    "feature_rt", "feature_area", "feature_norm_intensity", "feature_conc", "feature_intensity", "feature_response",
    "feature_height", "feature_fwhm", "feature_width", "feature_symmetry", "feature_int_start", "feature_int_end"
  )
  if (!any(intensity_cols %in% names(d_raw))) {
    cli::cli_abort("No feature variable column found. PLease column mapping via `column_mapping` and data source")
  }

  # extract and pivot to improve speed
  d_sample <- d_raw |>
    dplyr::select(any_of(c("analysis_id", "raw_data_filename", "qc_type", "acquisition_time_stamp", "batch_id"))) |>
    dplyr::distinct()


  # Conditionally parse acquisition_time_stamp if it exists
  if ("acquisition_time_stamp" %in% names(d_sample)) {
    d_sample <- d_sample |>
      dplyr::mutate(
        acquisition_time_stamp = lubridate::parse_date_time(
          .data$acquisition_time_stamp,
          c("mdy_HM", "dmy_HM", "ymd_HM", "ydm_HM",
            "mdy_HM %p", "dmy_HM %p", "ymd_HM %p", "ydm_HM %p",
            "%y-%m-%d %H:%M:%S"),
          quiet = TRUE
        )
      )
  }



  missing_ids <- setdiff("feature_id", names(d_raw))
  if (length(missing_ids) > 0) {
    cli::cli_abort(col_red("Required `feature_id` column is missing. Please map column names via the `column_mapping` argument or verify data file."))
  }

  # Always squish all character columns that exist
  d_sample <- d_sample |>
    dplyr::mutate(dplyr::across(where(is.character), stringr::str_squish))

  d_feature <- d_raw |>
    dplyr::select(any_of(c("feature_id", "feature_class", "istd_feature_id"))) |>
    dplyr::distinct() |>
    dplyr::mutate(
      dplyr::across(
        any_of(c("feature_id", "istd_feature_id")),
        ~ stringr::str_squish(as.character(.))
      )
    )


  # finalize data



  d_raw_final <- d_raw |>
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("feature_") & !dplyr::matches(c("_name", "_id", "_class")),
        as.numeric
      )
    ) |>
    dplyr::mutate(dplyr::across(ends_with("_mz"), as.numeric)) |>
    dplyr::mutate(dplyr::across(ends_with("_energy"), as.numeric)) |>
    dplyr::mutate(integration_qualifier = FALSE) |>
    dplyr::select(-any_of(c("acquisition_time_stamp",  "raw_data_filename", "qc_type", "feature_class", "istd_feature_id", "batch_id"))) |>
    left_join(d_sample, by = "analysis_id") |>
    left_join(d_feature, by = "feature_id", keep = FALSE) |>
    relocate(any_of(c("analysis_id", "raw_data_filename", "acquisition_time_stamp", "qc_type", "batch_id", "feature_id", "istd_feature_id","feature_class", "integration_qualifier")))

  # if (!use_normalized_data) {
  #   #d_mrmkit_data <- d_mrmkit_data |> mutate(feature_norm_intensity = NA_real_)
  #   d_mrmkit_data <- d_mrmkit_data |> select(-feature_norm_intensity)
  # }

  d_raw_final

}



#' Parses a plain wide CSV file
#'
#' Parses a CSV table with analysis/samples in rows, features values in columns.
#'
#' @param path path name and path of a plain long-format CSV file
#' @param variable_name Name of the variable representing the values in the table. Must be one of "intensity", "norm_intensity", "conc", "area", "height", "response")
#' @param analysis_id_col Column to be used as analysis_id
#' @param import_metadata Import additional metadata columns (e.g. batch ID, sample type) and add to the `MidarExperiment` object
#' @param first_feature_column Column number of the first column representing the feature values
#' @param na_strings A character vector of strings which are to be interpreted as NA values. Blank fields are also considered to be missing values.
#' @return A tibble in the long format
#' @examples
#' file_path <- system.file("extdata", "plain_wide_dataset.csv", package = "midar")
#'
#' tbl <- parse_plain_wide_csv(
#'  path = file_path,
#'  variable_name = "conc",
#'  analysis_id_col = "analysis_id",
#'  import_metadata = TRUE)
#'
#' head(tbl)
#' @export
 parse_plain_wide_csv <- function(path, variable_name, analysis_id_col = NA, import_metadata = TRUE, first_feature_column = NA, na_strings = "NA") {

  if(fs::path_ext(path) != "csv") cli::cli_abort(col_red("Only csv files are currently supported."))

  variable_name <- str_remove(variable_name, "feature_")
  rlang::arg_match(variable_name, c("area", "height", "intensity", "norm_intensity", "response", "conc", "conc_raw", "rt", "fwhm"))
  variable_name <- stringr::str_c("feature_", variable_name)
  variable_name_sym = rlang::sym(variable_name)


  d <- readr::read_csv(path, col_names = TRUE, trim_ws = TRUE, locale = readr::locale(encoding = "UTF-8"), col_types = "?", na = na_strings, progress = FALSE, name_repair = "minimal")



  if (!is.null(analysis_id_col) & !is.na(analysis_id_col)){
    # a column name was provided as analysis_id
    if (!analysis_id_col %in% names(d)) {
      cli::cli_abort(col_red("No column with the name `{analysis_id_col}` found in the data file."))
    }
  } else {
    # no column name provided, try to find analysis_id first
    if ("analysis_id" %in% names(d)){
      analysis_id_col <- "analysis_id"
    } else {
      # no analysis_id column found
      cli::cli_abort(col_red("Column `analysis_id` not found in imported data. Please set the column to be used as `analysis_id` using the `analysis_id_col` argument"))
    }
  }


  n_dup <- sum(duplicated(names(d)))
  if(n_dup > 0) {
    cli::cli_abort(col_red("{n_dup} duplicated column name(s) present in the data file. Please verify the data."))
  }

  analysis_metadata_cols <- c(
    "analysis_order",
    "qc_type",
    "batch_id",
    "is_quantifier",
    "is_istd"
  )

  analysis_metadata_cols <- intersect(analysis_metadata_cols,names(d))
  if(length(analysis_metadata_cols) > 0) {
    txt <- glue::glue_collapse(analysis_metadata_cols, sep = ", ")
    if(import_metadata){
      if("analysis_order" %in% analysis_metadata_cols) {
        if(!is.numeric(d$analysis_order))
          cli::cli_abort(col_red("Column `analysis_order` must contain unique numbers. Please verify your data or use `import_metadata = FALSE`"))
        if(any(duplicated(d$analysis_order)))
          cli::cli_abort(col_red("`analysis_order` contains duplicated values. Please verify your data."))
      }

      d <- d |>
        mutate(
          across(
            .cols = dplyr::any_of(names(c("qc_type", "batch_id"))),
            .fns = \(x) str_squish(as.character(x))
          )
        ) |>
        mutate(
          across(
            .cols = dplyr::any_of(names(c("is_quantifier", "is_istd"))),
            .fns = \(x) as.logical(str_squish(as.character(x)))
          )
        ) |>
        mutate(
          across(
            .cols = dplyr::any_of(names(c("analysis_order"))),
            .fns = \(x) as.numeric(str_squish(as.character(x)))
          )
        )

       cli::cli_alert_info(col_green("Metadata column(s) '{txt}' imported. To ignore, set `import_metadata = FALSE`"))


    } else {

      d <- d |> select(-all_of(analysis_metadata_cols))
      cli::cli_alert_info(cli::col_grey("Metadata column(s) '{txt}' found and ignored. To import the metadata, set `import_metadata = TRUE`"))

    }
  } else {
    analysis_metadata_cols = NULL
  }





  analysis_id_col_sym = rlang::sym(analysis_id_col)

  d[[analysis_id_col]] <- factor(d[[analysis_id_col]])


  n_dup <- sum(duplicated(d[[analysis_id_col]]))
  if(n_dup > 0) {
    cli::cli_abort(col_red("{n_dup} duplicated `analysis_id` present in the data file. Please verify the data."))
  }


  # Verify all data columns are numeric

  if(!is.na(first_feature_column)) {
    if(!is.numeric(first_feature_column)) {
      if(!first_feature_column %in% names(d)) {
        cli::cli_abort(col_red("Column `{first_feature_column}` not found in the data file. Please very set argument `first_feature_column`"))
      } else {
        first_feature_column_ix <- which(names(d) == first_feature_column)
      }
    } else {
      if(first_feature_column > ncol(d) | first_feature_column < 2) {
        cli::cli_abort(col_red("Column index set via `first_feature_column` out of range. Please verify the set value or use the column name."))
      } else {
        first_feature_column_ix <- first_feature_column
      }
    }
    sel_feature_names <- d[, first_feature_column_ix:ncol(d)]|>  names()
  } else {
    sel_feature_names <- d  |>  names()
  }



  val_cols <- setdiff(sel_feature_names, c(analysis_id_col, analysis_metadata_cols))
  analysis_cols <- setdiff(names(d), c(val_cols))






  if(!all(purrr::map_lgl(d[val_cols], is.numeric))) {
      cli::cli_abort(col_red("All columns with feature values must be numeric. Check your data, or if metadata is present, use `first_feature_column` to set start of feature data. If missing values are present as text, use `na_strings` to set the values to be treated as missing."))
  }

  d_long <- d |>
    tidyr::pivot_longer(cols = -all_of(c(analysis_cols)), names_to = "feature_id", values_to = variable_name) |>
    dplyr::rename(analysis_id = !!analysis_id_col_sym) |>
    dplyr::mutate(across(any_of(c("analysis_id", "batch_id", "batch_id")), as.character)) |>
    dplyr::mutate(integration_qualifier = FALSE)

  d_long
  }





#' #' Reads a generic wide .csv or .xls sheet with analysis results data
#' #'
#' #' @description
#' #' Imports a tables with the column being the analysis identifier (e.g. sample name) and the following colums being measured features. The values in the table
#' #'
#' #'
#' #' @param path File name and path of a plain wide-format CSV or XLS file
#' #' @param value_type Parameter type of the feature, i.e. "area", "height", "intensity", "norm_area", "norm_height", "norm_intensity", "feature_conc", "rt", "fwhm" or "width"
#' #' @param sheet Sheet name in case an Excel file (.xls, .xlsx, .xlsm) is imported
#' #' @param silent Suppress messages
#'
#' #' @return A tibble in the long format
#' #' @noRd
#' read_data_table <- function(path, value_type = c("area", "height", "intensity", "norm_area", "norm_height", "norm_intensity", "feature_conc", "rt", "fwhm", "width"), sheet = "", silent = FALSE) {
#'   value_type <- match.arg(value_type)
#'
#'   value_type <- dplyr::case_match(
#'     value_type,
#'     "area" ~ "feature_area",
#'     "height" ~ "feature_height",
#'     "intensity" ~ "feature_intensity",
#'     "norm_area" ~ "feature_norm_area",
#'     "norm_height" ~ "feature_norm_height",
#'     "norm_intensity" ~ "feature_norm_intensity",
#'     "feature_conc" ~ "feature_conc",
#'     "rt" ~ "feature_rt",
#'     "fwhm" ~ "feature_fwhm",
#'     "width" ~ "feature_width",
#'     .default = NULL
#'   )
#'
#'   var_value_type <- rlang::ensym(value_type)
#'   ext <- fs::path_ext(path)
#'
#'   if (ext == "csv") {
#'     d <- readr::read_csv(path, col_names = TRUE, trim_ws = TRUE, progress = FALSE,
#'                          na = c("n/a", "N/A", "NA", "na", "ND", "N.D.", "n.d."),
#'                          col_types = "cn", locale = readr::locale(encoding = "UTF-8"))
#'   } else if (ext == "xls" || ext == "xlsx" || ext == "xlm" || ext == "xlmx") {
#'     if (sheet == "") cli::cli_abort("Please define sheet name via the `sheet` parameter")
#'     d <- .read_generic_excel(path, sheet)
#'   } else {
#'     cli::cli_abort("Invalid file format. Only CSV, XLS and XLSX supported.")
#'   }
#'
#'   d |>
#'     dplyr::mutate(analysis_order = dplyr::row_number(), .before = 1) |>
#'     tidyr::pivot_longer(cols = -1:-2, names_to = "feature_id", values_to = value_type) |>
#'     dplyr::rename(analysis_id = 2) |>
#'     mutate({{ value_type }} := as.numeric(!!var_value_type))
#' }
#'




#'
#'
#' #' @title internal method to read excel sheets
#' #' @param path excel file name
#' #' @return tibble table
#' #' @noRd
#' .read_generic_excel <- function(path, sheetname) {
#'   data <- openxlsx2::read_xlsx(file = path,
#'                               sheet = sheetname,
#'                               col_names = TRUE,
#'                               skip_empty_rows = TRUE,
#'                               skip_empty_cols = TRUE,
#'                               detect_dates = TRUE,
#'                               na_strings = c("n/a", "N/A", "NA"),
#'                               rows = 1)|>
#'     mutate(across(where(is.character), str_trim))
#'
#'   nms <- names(data)
#'
#'   n_col <- length(nms) - 1
#'   c_numeric <- rep("numeric", n_col)
#'   ct <- c("character", c_numeric)
#'   data <- openxlsx2::read_xlsx(file = path,
#'                               sheet = sheetname,
#'                               col_names = TRUE,
#'                               detect_dates = TRUE,
#'                               skip_empty_rows = TRUE,
#'                               skip_empty_cols = TRUE,
#'                               na_strings = c("n/a", "N/A", "NA")
#'                               #types = ct
#'                               ) |>
#'     mutate(across(where(is.character), str_trim))
#'   data
#' }





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
#     dplyr::mutate(file_analysis_order = dplyr::row_number(), .before = "name") |>
#     tidyr::pivot_longer(-.data$file_analysis_order:-.data$name, names_to = "feature_id", values_to = "value") |>
#     dplyr::rename(raw_data_filename = .data$name) |>
#     dplyr::mutate(raw_data_filename = stringr::str_remove(.data$raw_data_filename, stringr::regex("\\.mzML$|\\.d$|\\.raw$|\\.wiff$|\\.wiff2$|\\.lcd$", ignore_case = TRUE))) |>
#     dplyr::mutate(raw_data_filename = stringr::str_squish(.data$raw_data_filename)) |>
#     dplyr::mutate(feature_id = dplyr::if_else(stringr::str_detect(.data$feature_id, "RT"), stringr::str_squish(stringr::str_extract(.data$feature_id, ".*(?= RT)")), .data$feature_id)) |>
#     dplyr::mutate(dplyr::across(.data$value, as.numeric)) |>
#     dplyr::left_join(d_mrmkit_featureinfo, by = "feature_id") |>
#     dplyr::relocate(.data$value, .after = dplyr::last_col()) |>
#     dplyr::mutate(dplyr::across(where(is.character), stringr::str_squish))
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
#   #     #dplyr::mutate(analysis_order = row_number(), .before = ANALYSIS_ID) |>
#   #     dplyr::select(-any_of(c("batch", "type"))) |>
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
