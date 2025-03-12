---
title: "MiDAR: Quantitative Mass Spectrometry Data Processing and Quality Control "
---

# MiDAR <img src="man/figures/logo.png" align="right" height="139"/>

<!-- badges: start -->

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental) [![Codecov test coverage](https://codecov.io/gh/SLINGhub/midar/branch/master/graph/badge.svg)](https://app.codecov.io/gh/SLINGhub/midar?branch=master)

<!-- badges: end -->

MiDAR is an R package for reproducible post-processing, quality control, and reporting of quantitative small-molecule mass spectrometry (MS) data.It offers a complete workflow, allowing users to import data, apply normalization and quantification methods, perform isotope correction, and address drift and batch effects.Additionally, MiDAR supports feature filtering, sharing curated datasets in various formats, and generating quality control plots and metrics to evaluate analytical data quality and the effects of post-processing steps.

The intended users of MiDAR are analytical and bioinformatics scientists. Its core tools are accessible to users with a basic understanding of R or coding, while also allowing for more advanced customisations and access to all data.The package supports analysts in annotating, inspecting, and processing their data independently, and facilitates the sharing of detailed annotated datasets and processing information for collaboration with colleagues, including bioinformaticians, for further analyses.

MiDAR emphasises fully documented, reproducible data processing workflows and also serves as a software library for creating automated data processing pipelines.
## Installation

In the console run:

``` r
if (!require("pak")) install.packages("pak")
pak::pkg_install("SLINGhub/midar")
```

## Example workflow

``` r
library(midar)

# Create a MidarExperiment object
myexp <- MidarExperiment()

# Load data and available metadata
myexp <- import_data_mrmkit(myexp, path = "data.tsv")
myexp <- import_metadata_msorganizer(myexp, path = "metadata.csv")


# Normalize and quantitate each feature by internal standards
myexp <- normalize_by_istd(myexp)
myexp <- quantify_by_istd(myexp)

# Filter features according to QC criteria
mexp <- filter_features_qc(
  data = mexp, 
  max.cv.conc.bqc = 30,
  min.min.min.min.signalblank.median.spl.pblk = 3,
)

# Export concentration data
myexp <- save_dataset_csv(
  myexp, 
  path = "mydata.csv", 
  variable = "norm_intensity", 
  filter_data = FALSE)
```
