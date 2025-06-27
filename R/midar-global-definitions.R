
# Initializing
pkg.env <- new.env()
setOldClass(c("tbl_df", "tbl", "data.frame")) # allow S4 to see the S3 tbl_df


# Data structure templates
pkg.env$table_templates <- list(
  dataset_orig_template = dplyr::tibble(
    "analysis_id" = character(),
    "raw_data_filename" = character(),
    "acquisition_time_stamp" = as.Date(character()),
    "feature_id" = character(),
  ),
  dataset_template = dplyr::tibble(
    "analysis_order" = integer(),
    "analysis_id" = character(),
    "acquisition_time_stamp" = as.Date(character()),
    "qc_type" = factor(),
    "batch_id" = character(),
    "specimen" = character(),
    "sample_id" = character(),
    "replicate_no" = integer(),
    "feature_id" = character(),
    "feature_class" = character(),
    "feature_label" = character(),
    "is_istd" = logical(),
    "is_quantifier" = logical(),
    "feature_intensity" = numeric()

  ),
  annot_analyses_template = dplyr::tibble(
    "analysis_order" = integer(),
    "analysis_id" = character(),
    "sample_id" = character(),
    "qc_type" = factor(),
    "batch_id" = character(),
    "replicate_no" = integer(),
    "specimen" = character(),
    "sample_amount" = numeric(),
    "sample_amount_unit" = character(),
    "istd_volume" = numeric(),
    "valid_analysis" = logical(),
    "annot_order_num" = integer(),
    "remarks" = character()
  ),
  annot_features_template = dplyr::tibble(
    "feature_id" = character(),
    "feature_class" = character(),
    "analyte_id" = character(),
    "chem_formula" = character(),
    "molecular_weight" = numeric(),
    "is_istd" = logical(),
    "istd_feature_id" = character(),
    "quant_istd_feature_id" = character(),
    "is_quantifier" = logical(),
    "valid_feature" = logical(),
    "response_factor" = numeric(),
    "interference_feature_id" = character(),
    "interference_contribution" = numeric(),
    "curve_fit_model" = character(),
    "curve_fit_weighting" = character(),
    "remarks" = character()
  ),
  annot_istds_template = dplyr::tibble(
    "istd_feature_id" = character(),
    "quant_istd_feature_id" = character(),
    "istd_conc_nmolar" = numeric()
  ),
  annot_responsecurves_template = dplyr::tibble(
    "analysis_id" = character(),
    "curve_id" = character(),
    "analyzed_amount" = numeric()
  ),
  annot_qcconcentrations_template = dplyr::tibble(
    "sample_id" = character(),
    "analyte_id" = character(),
    "concentration" = numeric(),
    "concentration_unit" = character(),
    "include_in_analysis" = logical()
  ),
  annot_batch_info_template = dplyr::tibble(
    "batch_id" = character(),
    "batch_no" = numeric(),
    "id_batch_start" = numeric(),
    "id_batch_end" = numeric()
  ),
  parameters_processing_template = dplyr::tibble(
    "parameter_name" = character()
  )
)


pkg.env$qc_type_annotation <- list(
  qc_type_levels = c(
    "SBLK", "TBLK", "UBLK", "HQC", "MQC", "LQC",  "QC", "PBLK", "CAL","EQA", "PQC", "TQC", "BQC", "RQC", "EQC", "NIST",
    "LTR",  "SPL", "SST", "MBLK"
  ),
  qc_type_col = c(
    "SBLK" = "#1854f9",
    "TBLK" = "#db0202",
    "UBLK" = "#de21de",
    "BQC" = "#db0202",
    "TQC" = "#1854f9",
    "PQC" = "#f99f18",
    "LQC" = "#f27507",
    "MQC" = "#f27507",
    "HQC" = "#f27507",
    "QC" = "#f27507",
    "CAL" = "#4575b4",
    "RQC" = "#96a4ff",
    "EQC" = "#513c3c",
    "NIST" = "#002e6b",
    "LTR" = "#880391",
    "EQA" = "#880391",
    "PBLK" = "#216651",
    "SPL" = "#8e9b9e",
    "SST" = "#bafc03",
    "MBLK" = "black"
  ),
  qc_type_fillcol = c(
    "SBLK" = "#f891ff",
    "TBLK" = "#fffb03",
    "UBLK" = "#c1bd04",
    "BQC" = "#db0202",
    "TQC" = "#1854f9",
    "PQC" = "#f99f18",
    "LQC" = "#f5c969",
    "MQC" = "#f5c969",
    "HQC" = "#f5c969",
    "QC" = "#f5c969",
    "EQA" = "#de21de",
    "CAL" = "white",
    "RQC" = "#688ff9",
    "EQC" = "NA",
    "NIST" = "#cce2ff",
    "LTR" = "#880391",
    "PBLK" = "#e4f2c4",
    "SPL" = "NA",
    "SST" = "#aaaeaf",
    "MBLK" = "black"
  ),
  qc_type_shape = c(
    "SBLK" = 23,
    "TBLK" = 23,
    "UBLK" = 23,
    "BQC" = 16,
    "TQC" = 25,
    "PQC" = 23,
    "LQC" = 25,
    "MQC" = 23,
    "HQC" = 24,
    "QC" = 23,
    "EQA" = 22,
    "CAL" = 21,
    "RQC" = 6,
    "EQC" = 24,
    "NIST" = 23,
    "LTR" = 23,
    "PBLK" = 23,
    "SPL" = 21,
    "SST" = 10,
    "MBLK" = 10
  )
)
