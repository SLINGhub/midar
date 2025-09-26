## ----pkgs, include=FALSE------------------------------------------------------
library(midar)
knitr::opts_chunk$set(
  collapse = TRUE,
  warning = FALSE,
  comment = "#>"
)

## ----eval=FALSE---------------------------------------------------------------
# if (!require("pak")) install.packages("pak")
# pak::pkg_install("SLINGhub/midar")

## -----------------------------------------------------------------------------
# Path of example files included with this package
dir_path <- system.file("extdata",  package = "midar")

# Create a MidarExperiment object
mexp <- MidarExperiment()

# Load data and metadata
mexp <- import_data_mrmkit(mexp,
                           path = file.path(dir_path, "MRMkit_demo.tsv"),
                           import_metadata = TRUE)

mexp <- import_metadata_analyses(mexp, 
                                 path = file.path(dir_path, "MRMkit_AnalysesAnnot.csv"), 
                                 ignore_warnings = T)
mexp <- import_metadata_istds(mexp, 
                              path = file.path(dir_path, "MRMkit_ISTDconc.csv"))

# Normalize and quantitate features by internal standards
mexp <- normalize_by_istd(mexp)
mexp <- quantify_by_istd(mexp)

# Drift and batch-effect correction
mexp <- correct_drift_cubicspline(mexp, variable = "conc", ref_qc_types = "BQC")
mexp <- correct_batch_centering(mexp, variable = "conc", ref_qc_types = "BQC")

# Plot analysis time-trends of final concentrations of each feature 
plot_runscatter(mexp,
                variable = "conc",
                qc_types = c("BQC", "TQC", "SPL", "PBLK", "SBLK"),
                cap_outliers = TRUE,
                output_pdf = FALSE,
                path = "./output/runscatter_istd.pdf")

# Set features QC-filter criteria   
 mexp <- filter_features_qc(mexp,
                            max.cv.conc.bqc = 25,
                            min.signalblank.median.spl.pblk = 3,
                            include_qualifier = FALSE,
                            include_istd = FALSE)
 
# Save concentration data
 save_dataset_csv( mexp, 
                   path = "mydata.csv", 
                   variable = "conc", 
                   filter_data = TRUE)

