library(tidyverse)
library(midar)
library(usethis)
library(here)


# Create an unprocessed lipidomics dataset with data and metadata as a MidarExperiment object

lipidomics_dataset <- MidarExperiment(analysis_type = "lipidomics")

lipidomics_dataset <- import_data_mrmkit(data = lipidomics_dataset,
                                         path = here("data-raw/MRMkit_example.tsv"),
                                         import_metadata = T)

lipidomics_dataset <- import_metadata_midarxlm(lipidomics_dataset,
                                 here("data-raw/MRMkit_example_metadata.xlsm"),
                                 excl_unmatched_analyses = T)


usethis::use_data(lipidomics_dataset, overwrite = TRUE,compress = "xz")

save_dataset_csv(data = lipidomics_dataset, path = here("data-raw/MRMkit_example.csv"),
                 variable = "intensity", filter_data = FALSE, add_qctype = TRUE, qc_types = c("SPL", "BQC"))


#mexp <- normalize_by_istd(lipidomics_dataset)
#mexp <- quantify_by_istd(mexp <- normalize_by_istd(lipidomics_dataset))

# Create an unprocessed quantitative LC-MS (with calibration curve and QC levels)
# with data and metadata as a MidarExperiment object
mexp <- MidarExperiment()

mexp <- import_data_masshunter(data = mexp,
                                         path = here("data-raw/QuantLCMS_Example_MassHunter.csv"),
                                         import_metadata = T)

quant_lcms_dataset <- import_metadata_midarxlm(mexp,
                                               here("data-raw/QuantLCMS_Example_Metadata.xlsm"),
                                               excl_unmatched_analyses = T, ignore_warnings = T)

usethis::use_data(quant_lcms_dataset, overwrite = TRUE,compress = "xz")
