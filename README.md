# MiDAR <img src="man/figures/logo.svg" align="right" height="139"/>

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) [![Codecov test coverage](https://codecov.io/gh/SLINGhub/midar/branch/master/graph/badge.svg)](https://app.codecov.io/gh/SLINGhub/midar?branch=master)

<!-- badges: end -->

`MiDAR` is an R package to reproducibly manage, post-process, visualize, apply quality control, and analyze small-molecule mass spectrometry (MS) datasets, e.g. from targeted lipidomics and metabolomics experiments.

`MiDAR` is tailored to handle different analytical designs, data types and data processing strategies. As such, this package provides functions to import data files from from different commercial and open-source tools. Data processing functions include, internal standard and sample amount-based normalization, quantification, as well as drift and batch corrections. Quality control (QC) functions provide QC metrics and plots of raw and processed data, and QC-based feature filtering.  

Datasets and processing steps are tracked, and can be saved as `MiDAR` S4 class object (RDS) files, Excel, PowerPoint and interactive HTML-based reports, enabling sharing and reporting the of all data, metadata and applied data processing steps.


## Installation

`{midar}` is currently only available via GitHub:

``` r
if (!require("pak")) install.packages("pak")
pak::pkg_install("SLINGhub/midar")
```

## Example

``` r
library(midar)

#' 
#' mexp <- MidarExperiment()
#' mexp <- 
#' mexp

# Get paths of example files included with this package
file_path <- system.file("extdata", "sPerfect_MRMkit.tsv", package = "midar", mustWork = TRUE)

# Create a MidarExperiment object (S4)
mexp <- MidarExperiment()

# Load data and available metadata from MRMkit output file
mexp <- data_import_mrmkit(data = mexp, path = file_path, use_metadata = TRUE)


# Normalize and quantitate each feature by internal standards
mexp <- calc_normalize_by_istd(mexp)
mexp <- calc_quant_by_istd(mexp)

# Get QC metrics for each feature
mexp <- qc_calc_metrics(mexp)

# Filter features according to QC criteria
mexp <- qc_filter_features(data = mexp,
                        min_cv_conc_bqc = 30,
                        min_intensity_bqc = 100,
                        min_signal_blank_ratio = 5,
                        min_response_rsquare = 0.8,
                        rqc_curve_used_for_filt = 1)
```

## Contributor Code of Conduct

Please note that the midar project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
