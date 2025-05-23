% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data-import.R
\name{import_data_csv_wide}
\alias{import_data_csv_wide}
\title{Import Analysis Results from Plain Wide-Format CSV Files}
\usage{
import_data_csv_wide(
  data = NULL,
  path,
  variable_name,
  analysis_id_col = NA,
  import_metadata = TRUE,
  first_feature_column = NA,
  na_strings = "NA"
)
}
\arguments{
\item{data}{A \code{MidarExperiment} object.}

\item{path}{A file path or vector of file paths, or a directory path. If a directory is
provided, all \code{.csv} files within it will be read.}

\item{variable_name}{A character string specifying the variable type contained in the
data. Must be one of \code{"intensity"}, \code{"norm_intensity"}, \code{"conc"}, \code{"area"}, \code{"height"}, or \code{"response"}.}

\item{analysis_id_col}{The column name or index to be used as \code{analysis_id}. Defaults to \code{NA},
in which case \code{"analysis_id"} is used if present; otherwise, the first column is used if it contains unique values.}

\item{import_metadata}{Logical indicating whether to import additional metadata columns
(e.g., batch ID, sample type) into the \code{MidarExperiment} object. Supported metadata columns are:
\code{"qc_type"}, \code{"batch_id"}, \code{"is_quantifier"}, \code{"is_istd"}, and \code{"analysis_order"}.}

\item{first_feature_column}{Integer indicating the column number where feature value columns start.}

\item{na_strings}{Character vector of strings to interpret as NA values. Blank fields are also treated as NA.}
}
\value{
A \code{MidarExperiment} object containing the imported dataset.
}
\description{
Imports analysis result data from wide-format \code{.csv} files, where each row corresponds
to a unique analysis-feature pair and columns contain analysis- or feature-specific variables.
}
\details{
The dataset must include two identifier columns: \code{"analysis_id"} and \code{"feature_id"},
with each pair of values unique across the table. Additionally, the table must contain
at least one feature variable column, such as \code{"area"}, \code{"height"}, \code{"intensity"},
\code{"norm_intensity"}, \code{"response"}, \code{"conc"}, \code{"rt"}, or \code{"fwhm"}. Some downstream
functions may require specific columns among these to be present.

The \code{variable_name} argument specifies the data type represented in the table, which
must be one of: \code{"area"}, \code{"height"}, \code{"intensity"}, \code{"norm_intensity"}, \code{"response"},
\code{"conc"}, \code{"conc_raw"}, \code{"rt"}, or \code{"fwhm"}.

If there is no column named \code{analysis_id}, it will be inferred from the first column,
provided it contains unique values.

When \code{import_metadata} is set to \code{TRUE}, the following metadata columns will be imported
if present:
\itemize{
\item \code{analysis_order}
\item \code{qc_type}
\item \code{batch_id}
\item \code{is_quantifier}
}

To prevent additional non-metadata columns from being misinterpreted as features,
use the \code{first_feature_column} parameter to specify the column where feature data starts.

If a directory path is provided to \code{path}, all \code{.csv} files in that directory will be
processed and merged into a single dataset. This facilitates handling datasets split
into multiple files during preprocessing. Ensure each feature and raw data file pair
appears only once to avoid duplication errors.

The \code{na_strings} parameter allows specifying character strings to be interpreted as
missing values (NA). Blank fields are also treated as missing.
}
\examples{
file_path <- system.file("extdata", "plain_wide_dataset.csv", package = "midar")
mexp <- MidarExperiment()
mexp <- import_data_csv_wide(
  data = mexp,
  path = file_path,
  variable_name = "conc",
  import_metadata = TRUE
)
print(mexp)

}
