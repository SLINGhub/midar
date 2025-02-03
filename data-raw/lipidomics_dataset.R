library(tidyverse)
library(midar)
library(usethis)
library(here)


# Create an unprocessed dataset with data and metadata as a MidarExperiment object

lipidomics_dataset <- MidarExperiment(analysis_type = "lipidomics")

lipidomics_dataset <- import_data_mrmkit(data = lipidomics_dataset,
                                         path = here("data-raw/MRMkit_example.tsv"),
                                         import_metadata = T)

lipidomics_dataset <- import_metadata_midarxlm(lipidomics_dataset,
                                 here("data-raw/MRMkit_example_metadata.xlsm"),
                                 excl_unmatched_analyses = T)

usethis::use_data(lipidomics_dataset, overwrite = TRUE,compress = "xz")
