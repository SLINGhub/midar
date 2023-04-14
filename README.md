# MiDAR <img src="man/figures/logo.svg" align="right" height="139"/>

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) [![Codecov test coverage](https://codecov.io/gh/SLINGhub/midar/branch/master/graph/badge.svg)](https://app.codecov.io/gh/SLINGhub/midar?branch=master)

<!-- badges: end -->

`MiDAR` is an R package to manage, post-process, inspect and analyze small-molecule mass spectrometry (MS) raw datasets. Such datasets include pre-processed datasets from targeted lipidomics and metabolomics experiments.

`MiDAR` is tailored to handle a range of different real-life analytical designs, data types and processing strategies. One focus of MiDAR is on using data formats from diferent MS instrument types and data pre-processing softwares, and to collected all required metadata in structures tidy format.

`MiDAR` provides functions are provided to inspect and visualize the data sets, before and after processing steps, to inspect quality of the data, and effects of data processing steps. A basic set of exploratory data analysis and vizualization options is also provided. Finally, this package includes function to report the processed datasets, quality control results, and processing details to Excel, PowerPoint and interactive HTML-based reports.

## Installation

MiDAR is currently only available via installation from GitHub:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("SLINGhub/midar")
```

## Example

``` r
library(midar)

# Get paths of example files included with this package
masshunter_file <- system.file("extdata", "Example_MHQuant_1.csv", package = "midar", mustWork = TRUE)
metadata_file <- system.file("extdata", "Example_Metadata_1.xlsm", package = "midar", mustWork = TRUE)

# Create a MidarExperiment object (S4)
mexp <- MidarExperiment()

# Load data and metadata
mexp <- loadMasshunterCSV(data = mexp, filename = masshunter_file)
mexp <- loadMSOrganizerXLM(data = mexp, filename = metadata_file)

# Normalize and quantitate each feature by internal standards
mexp <- normalize_by_istd(mexp)
mexp <- quantitate_by_istd(mexp)

# Get QC metrics for each feature
mexp <- calculate_qc_metrics(mexp)

# Filter features according to QC criteria
mexp <- apply_qc_filter(data = mexp,
                        CV_BQC_max = 30,
                        Intensity_BQC_min = 100,
                        SB_RATIO_min = 5,
                        R2_min = 0.8,
                        RQC_CURVE = 1)
```

## Contributor Code of Conduct

Please note that the midar project is released with a [Contributor Code of Conduct](https://contributor-covenant.org/version/2/0/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.
