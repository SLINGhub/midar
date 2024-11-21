
# Initializing
pkg.env <- new.env()
setOldClass(c("tbl_df", "tbl", "data.frame")) # allow S4 to see the S3 tbl_df


# Data structure templates
pkg.env$table_templates <- list(
  dataset_orig_template = tibble::tibble(
    "analysis_id" = character(),
    "raw_data_filename" = character(),
    "acquisition_time_stamp" = as.Date(character()),
    "feature_id" = character(),
  ),
  dataset_template = tibble::tibble(
    "run_seq_num" = integer(),
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
  annot_analyses_template = tibble::tibble(
    "run_seq_num" = integer(),
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
  annot_features_template = tibble::tibble(
    "feature_id" = character(),
    "feature_class" = character(),
    "is_istd" = logical(),
    "istd_feature_id" = character(),
    "quant_istd_feature_id" = character(),
    "is_quantifier" = logical(),
    "valid_feature" = logical(),
    "response_factor" = numeric(),
    "interference_feature_id" = character(),
    "interference_proportion" = numeric(),
    "remarks" = character()
  ),
  annot_istds_template = tibble::tibble(
    "istd_feature_id" = character(),
    "quant_istd_feature_id" = character(),
    "istd_conc_nmolar" = numeric()
  ),
  annot_responsecurves_template = tibble::tibble(
    "analysis_id" = character(),
    "curve_id" = character(),
    "analyzed_amount" = numeric()
  ),
  annot_batch_info_template = tibble::tibble(
    "batch_id" = character(),
    "batch_no" = numeric(),
    "id_batch_start" = numeric(),
    "id_batch_end" = numeric()
  ),
  parameters_processing_template = tibble::tibble(
    "parameter_name" = character()
  )
)


pkg.env$qc_type_annotation <- list(
  qc_type_levels = c(
    "SBLK", "TBLK", "UBLK", "PQC", "TQC", "BQC", "RQC", "EQC", "NIST",
    "LTR", "PBLK", "SPL", "SST", "MBLK"
  ),
  qc_type_col = c(
    "SBLK" = "#1854f9",
    "TBLK" = "#db0202",
    "UBLK" = "#de21de",
    "BQC" = "#db0202",
    "TQC" = "#1854f9",
    "PQC" = "#f99f18",
    "RQC" = "#96a4ff",
    "EQC" = "#513c3c",
    "NIST" = "#002e6b",
    "LTR" = "#880391",
    "PBLK" = "#08c105",
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
    "RQC" = "#688ff9",
    "EQC" = "NA",
    "NIST" = "#cce2ff",
    "LTR" = "#880391",
    "PBLK" = "#08c105",
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
    "PQC" = 25,
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
