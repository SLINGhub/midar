#' Calculate Quality Control (QC) Metrics for Features
#'
#' @description
#' Computes various quality control (QC) metrics for each feature in a
#' `MidarExperiment` object. Metrics are derived from different sample
#' types and can be computed either across the full dataset or as medians
#' of batch-wise calculations.
#'
#' @details
#'
#' **Batch-wise calculations**:
#' The function computes the following QC metrics for each feature and for
#' different QC sample types (e.g., SPL, TQC, BQC, PBLK, NIST, LTR)
#'
#' The format for the metrics is standardized as `metric_name_qc_type`, where
#' `qc_type` refers to the specific QC sample type for which the metric is
#' calculated. For example: `intensity_min_spl` refers to the minimum intensity

#' Statistics of normalized intensities , external calibration, and response
#' curves can be included by setting the relevant arguments
#'  (`include_norm_intensity_stats`, `include_conc_stats`,
#'  `include_response_stats`, `include_calibration_results`) to `TRUE`.
#'
#'  **Note** when corresponding underlying processed data is not available,
#'  the function will not raise an error but will return `NA` values for the
#'  respective metrics. This, however, does not apply for the optinal metrics
#'  mentioned above. For these cases an error will be raised if the underlying
#'  data is missing.
#'
#' If `use_batch_medians = TRUE`, batch-specific QC statistics are computed
#' first, and then the median of these values is returned for each feature.
#' However, response curve and calibration statistics are calculated per
#' curve, irrespective of batches and `use_batch_medians` settings.
#'
#' The calculated metrics are stored in the `metrics_qc` table of the
#' `MidarExperiment` objects and comprises following details
#'
#' - **Feature details**:
#'  Specific feature information extracted from the feature metadata tanle,
#'  such as feature class, associated ISTD, quantifier status.
#'
#' - **Feature MS Method Information** (if method variables are available in the analysis data).
#'   Extracts and summarizes method-related variables for each feature. If multiple
#'   values exist for the same feature, these will be concatenated into a string.
#'   The latter would indicate inconsistent analysis conditions.
#'   - `precursor_mz`: The m/z value of the precursor ion(s),
#'   - `product_mz`: The m/z value of the product ion(s), concatenated if multiple values exist for the same feature.
#'   - `collision_energy`: The collision energy used for fragmentation, concatenated if multiple values exist exist for the same feature.

#' - **Missing Value Metrics**:
#'   - `missing_intensity_prop_spl`: Proportion of missing intensities for the SPL sample type.
#'   - `missing_norm_intensity_prop_spl`: Proportion of missing normalized intensities for SPL samples.
#'   - `missing_conc_prop_spl`: Proportion of missing concentration values for SPL samples.
#'   - `na_in_all`: Indicator if a feature has all missing intensities across all samples
#'
#' - **Retention Time (RT) Metrics**: Requires that retention tim data are available.
#'   - `rt_min_*`: Minimum retention time across different QC sample types (e.g., SPL, BQC, TQC).
#'   - `rt_max_*`: Maximum retention time across different QC sample types.
#'   - `rt_median_*`: Median retention time for specific QC sample types like PBLK, SPL, BQC, TQC, etc.
#'
#' - **Intensity Metrics**:
#'   - `intensity_min_*`: Minimum intensity value for features across different QC sample types such as SPL, TQC, BQC, etc.
#'   - `intensity_max_*`: Maximum intensity values across sample types.
#'   - `intensity_median_*`: Median intensity for various QC sample types.
#'   - `intensity_cv_*`: Coefficient of variation (CV) of intensity values for specific QC types.
#'   - `sb_ratio_*`: Signal-to-blank ratios such as the ratio of intensity values for SPL vs PBLK, UBLK, or SBLK.
#'   - `intensity_q10_*`: The 10th percentile of intensity values for the SPL sample type.
#'
#' - **Normalized Intensity Metrics** (only if `include_norm_intensity_stats = TRUE`):
#'  Requires that raw intensities  were normalized, see [normalize_by_istd()]
#'  for details.
#'   - `norm_intensity_cv_*`: Coefficient of variation (CV) of normalized intensities for QC sample types like TQC, BQC, SPL, etc.
#'
#' - **Concentration Metrics** (only if `include_conc_stats = TRUE`):
#' Requires that concentration were calculated, see [quantify_by_istd()] or
#' [quantify_by_calibration()] for details.
#'   - `conc_median_*`: Median concentration values for different QC sample types like TQC, BQC, SPL, NIST, and LTR.
#'   - `conc_cv_*`: Coefficient of variation (CV) for concentration values.
#'   - `conc_dratio_sd_*`: The ratio of standard deviations of concentration between BQC or TQC and SPL samples.
#'   - `conc_dratio_mad_*`: The ratio of median absolute deviations (MAD) between BQC or TQC and SPL concentrations.
#'
#' - **Response Curve Metrics** (if `include_response_stats = TRUE`):
#'   Calculates response curve statistics for each feature and each curve
#'   (where `#` refers to the curve identifier). Requires that response curves
#'   are defined in the data. See [get_response_curve_stats()] for additional details.
#'   - `r2_rqc_#`: R-squared value of the linear regression for the response
#'     curve, representing the goodness of fit.
#'   - `slopenorm_rqc_#`: Normalized slope of the linear regression for the
#'     response curve, indicating the relationship between the response and
#'     concentration.
#'   - `y0norm_rqc_#`: Normalized intercept of the linear regression for the
#'     response curve, representing the baseline or starting value.
#'
#' - **External Calibration Results** Incorporates external calibration results,
#'  if `include_calibration_results = TRUE` and calibration curves are defined
#'  in the data:
#'   - `fit_model`: The regression model used for curve fitting.
#'   - `fit_weighting`: The weighting method applied during curve fitting.
#'   - `lowest_cal`: The lowest nonzero calibration concentration.
#'   - `highest_cal`: The highest calibration concentration.
#'   - `r.squared`: R-squared value indicating the goodness of fit.
#'   - `coef_a`:
#'     - For **linear fits**, this represents the slope of the regression line.
#'     - For **quadratic fits**, this represents the coefficient of the quadratic term (`x²`).
#'   - `coef_b`:
#'     - For **linear fits**, this represents the intercept of the regression line.
#'     - For **quadratic fits**, this represents the coefficient of the linear term (`x`).
#'   - `coef_c`:
#'     - Only present for **quadratic fits**, representing the intercept of the regression equation.
#'     - Set to `NA` for linear fits.
#'   - `sigma`: The residual standard error of the regression model.
#'   - `reg_failed`: Boolean flag indicating if regression fitting failed.
#'
#' @param data A `MidarExperiment` object containing data and metadata, whereby
#' data needs to be normalized and quantitated for specific QC metrics, such as
#' statistics based on normalized intensities and concentrations.
#' @param use_batch_medians Logical, whether to compute QC metrics using the
#'   median of batch-wise derived values instead of the full dataset. Default is
#'   FALSE.
#' @param include_norm_intensity_stats Logical. If `NA` (default), statistics on
#' normalized intensity values are included if the data is available. If `TRUE`,
#' they are always calculated, raising an error if data is missing.
#' @param include_conc_stats Logical. If `NA` (default), concentration-related
#' statistics are included if concentration data is available. If `TRUE`,
#' they are always calculated, raising an error if data is missing.
#' @param include_response_stats Logical. If `NA` (default), response curve statistics
#' are included if the required data is available. If `TRUE`, they are always
#' calculated, raising an error if data is missing.

#' @param include_calibration_results Logical, whether to incorporate external
#'   calibration results into the QC metrics table if available. Default is TRUE.
#'
#' @return A `MidarExperiment` object with an updated `metrics_qc` table
#'   containing computed QC metrics for each feature.
#' @export
calc_qc_metrics <- function(
    data = NULL,
    use_batch_medians = FALSE,
    include_norm_intensity_stats = NA,
    include_conc_stats = NA,
    include_response_stats = NA,
    include_calibration_results = NA
) {
  # Check if the input data is valid
  check_data(data)

  # If the analysis type is lipidomics, extract lipid class names
   # TODO: remove later when fixed
  if (!is.na(data@analysis_type) && tolower(data@analysis_type) == "lipidomics")
    data@dataset <- parse_lipid_feature_names(data@dataset,
                                 add_chain_composition = FALSE,
                                 use_as_feature_class = "lipid_class_lcb",
                                 add_transition_names = FALSE
                                 )

  # Select relevant feature information from the dataset
  d_feature_info <- data@annot_features |>
    select("valid_feature", "feature_id", "feature_class", "is_istd",
           "is_quantifier", "istd_feature_id", "quant_istd_feature_id",
           "response_factor")

  # Define method variables and template
  method_var <- c("method_precursor_mz", "method_product_mz", "method_collision_energy")
  d_method_template <- tibble("precursor_mz" = NA_character_,
                              "product_mz" = NA_character_,
                              "collision_energy" = NA_character_)

  # Check if the method variables exist in the dataset
  if (any(c(method_var) %in% names(data@dataset_orig))) {
    # Summarize method information for each feature
    d_method_info <- data@dataset_orig |>
      select(c("feature_id", any_of(method_var))) |>
      summarise(
        .by = "feature_id",
        precursor_mz = stringr::str_c(unique(.data$method_precursor_mz), collapse = "; "),
        product_mz = stringr::str_c(unique(.data$method_product_mz), collapse = "; "),
        collision_energy = stringr::str_c(unique(.data$method_collision_energy), collapse = "; ")
      ) |>
      ungroup() |>
      bind_rows(d_method_template) |>
      dplyr::mutate(across(where(is.character) & !c("feature_id"),
                           ~ as.numeric(suppressWarnings(as.numeric(.)))))
  } else {
    d_method_info <- d_feature_info |> select("feature_id") |> distinct() |>
      bind_rows(d_method_template)
  }


  # Summarize missing value statistics for different QC types
  d_stats_missingval <- data@dataset |>
    dplyr::filter(.data$qc_type %in% c("SPL", "NIST", "LTR", "BQC", "TQC", "PBLK",
                                       "SBLK", "UBLK", "IBLK", "CAL", "STD", "LQC",
                                       "MQC", "UQC")) |>
    dplyr::summarise(
      .by = "feature_id",
      missing_intensity_prop_spl = if ("feature_intensity" %in% names(data@dataset))
        sum(is.na(.data$feature_intensity[.data$qc_type == "SPL"])) /
        length(.data$feature_intensity[.data$qc_type == "SPL"]) else NA_real_,
      missing_norm_intensity_prop_spl = if ("feature_norm_intensity" %in% names(data@dataset))
        sum(is.na(.data$feature_norm_intensity[.data$qc_type == "SPL"])) /
        length(.data$feature_norm_intensity[.data$qc_type == "SPL"]) else NA_real_,
      missing_conc_prop_spl = if ("feature_conc" %in% names(data@dataset))
        sum(is.na(.data$feature_conc[.data$qc_type == "SPL"])) /
        length(.data$feature_conc[.data$qc_type == "SPL"]) else NA_real_,
      na_in_all = if ("feature_intensity" %in% names(data@dataset))
        all(is.na(.data$feature_intensity)) else NA_real_
    )

  # Set grouping variable depending on whether batch medians are used
  if (use_batch_medians)
    grp <- c("feature_id", "batch_id")
  else
    grp <- c("feature_id")

  # Select relevant variables needed for statistics
  d_stats_var <- data@dataset |>
    select(any_of(c("batch_id", "feature_id", "qc_type", "feature_rt",
                    "feature_intensity", "feature_conc", "feature_norm_intensity"))) |>
    filter(.data$qc_type != "RQC") |>
    mutate(qc_type = factor(.data$qc_type),
           batch_id = factor(.data$batch_id))

     # Using dtplyr to improve speed (2x), group by batch still 5x times slower than no batch grouping
     d_stats_var <- suppressMessages(dtplyr::lazy_dt(d_stats_var))

     # Summarize intensity statistics by group, feature (see before0)
     d_stats_var_final <- d_stats_var |> select(all_of(grp)) |> distinct()

     if ("feature_rt" %in% names(data@dataset)) {
       d_stats_var_rt <- d_stats_var |>
         summarise(
           .by = all_of(grp),

           rt_min_SPL = safe_min(.data$feature_rt[.data$qc_type == "SPL"], na.rm = TRUE),
           rt_max_SPL = safe_max(.data$feature_rt[.data$qc_type == "SPL"], na.rm = TRUE),
           rt_min_BQC = safe_min(.data$feature_rt[.data$qc_type == "BQC"], na.rm = TRUE),
           rt_min_TQC = safe_min(.data$feature_rt[.data$qc_type == "TQC"], na.rm = TRUE),
           rt_median_PBLK = median(.data$feature_rt[.data$qc_type == "PBLK"], na.rm = TRUE),
           rt_median_SPL = median(.data$feature_rt[.data$qc_type == "SPL"], na.rm = TRUE),
           rt_median_BQC = median(.data$feature_rt[.data$qc_type == "BQC"], na.rm = TRUE),
           rt_median_TQC = median(.data$feature_rt[.data$qc_type == "TQC"], na.rm = TRUE),
           rt_median_NIST = median(.data$feature_rt[.data$qc_type == "NIST"], na.rm = TRUE),
           rt_median_LTR = median(.data$feature_rt[.data$qc_type == "LTR"], na.rm = TRUE)
          )
       d_stats_var_final <- d_stats_var_final |> left_join(d_stats_var_rt, by = grp)
     }

     if ("feature_intensity" %in% names(data@dataset)) {
       d_stats_var_int <- d_stats_var |>
         summarise(
           .by = all_of(grp),

         intensity_min_SPL = safe_min(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
         intensity_max_SPL = safe_max(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),

         intensity_min_BQC = safe_min(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
         intensity_min_TQC = safe_min(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
         intensity_max_BQC = safe_max(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
         intensity_max_TQC = safe_max(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),

         intensity_median_PBLK = median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE),
         intensity_median_UBLK = median(.data$feature_intensity[.data$qc_type == "UBLK"], na.rm = TRUE),
         intensity_median_SBLK = median(.data$feature_intensity[.data$qc_type == "SBLK"], na.rm = TRUE),
         intensity_median_SPL = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
         intensity_median_BQC = median(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
         intensity_median_TQC = median(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
         intensity_median_NIST = median(.data$feature_intensity[.data$qc_type == "NIST"], na.rm = TRUE),
         intensity_median_LTR = median(.data$feature_intensity[.data$qc_type == "LTR"], na.rm = TRUE),

         intensity_cv_TQC = cv(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
         intensity_cv_BQC = cv(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
         intensity_cv_SPL = cv(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
         intensity_cv_LTR = cv(.data$feature_intensity[.data$qc_type == "LTR"], na.rm = TRUE),
         intensity_cv_NIST = cv(.data$feature_intensity[.data$qc_type == "NIST"], na.rm = TRUE),

         # Calculate quantiles within summarise
         intensity_q10_SPL = quantile(.data$feature_intensity[.data$qc_type == "SPL"], probs = 0.1, na.rm = TRUE, names = FALSE)
       ) |>
       mutate(
         sb_ratio_q10_pbk = .data$intensity_q10_SPL / .data$intensity_median_PBLK,
         sb_ratio_pblk = .data$intensity_median_SPL / .data$intensity_median_PBLK,
         sb_ratio_ublk = .data$intensity_median_SPL / .data$intensity_median_UBLK,
         sb_ratio_sblk = .data$intensity_median_SPL / .data$intensity_median_SBLK
       )
       d_stats_var_final <- d_stats_var_final |> left_join(d_stats_var_int, by = grp)
  }

    # table to collect all data
    if (!"feature_norm_intensity" %in% names(data@dataset)) {
      if (isTRUE(include_norm_intensity_stats)) {
        cli::cli_abort(col_red("Normalized intensity data is missing. Please normalize the data first using `normalize_by_*()` functions."))
      }
    } else if(is.na(include_norm_intensity_stats) || include_norm_intensity_stats){

      d_stats_var_norm_int <-  d_stats_var |>
        summarise(
         .by = all_of(grp),
          norm_intensity_cv_TQC = cv(.data$feature_norm_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
          norm_intensity_cv_BQC = cv(.data$feature_norm_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
          norm_intensity_cv_SPL = cv(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          norm_intensity_cv_LTR = cv(.data$feature_norm_intensity[.data$qc_type == "LTR"], na.rm = TRUE),
          norm_intensity_cv_NIST = cv(.data$feature_norm_intensity[.data$qc_type == "NIST"], na.rm = TRUE),
          normint_dratio_sd_bqc = sd(.data$feature_norm_intensity[.data$qc_type == "BQC"], na.rm = TRUE) / sd(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          normint_dratio_sd_tqc = sd(.data$feature_norm_intensity[.data$qc_type == "TQC"], na.rm = TRUE) / sd(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          normint_dratio_mad_bqc = mad(.data$feature_norm_intensity[.data$qc_type == "BQC"], na.rm = TRUE) / mad(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          normint_dratio_mad_tqc = mad(.data$feature_norm_intensity[.data$qc_type == "TQC"], na.rm = TRUE) / mad(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE)
         )
      d_stats_var_final <- d_stats_var_final |> left_join(d_stats_var_norm_int, by = grp)

    }


    if (!"feature_conc" %in% names(data@dataset)) {

      if (isTRUE(include_conc_stats)) {
        cli::cli_abort(col_red("Concentration data is missing. Please quantify the data first using `quantify_by_*()` functions."))
      }
    } else if(is.na(include_conc_stats) || include_conc_stats) {
      d_stats_var_conc <-  d_stats_var |>
        summarise(
          .by = all_of(grp),
           conc_median_TQC = median(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE),
           conc_median_BQC = median(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE),
           conc_median_SPL = median(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE),
           conc_median_NIST = median(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE),
           conc_median_LTR = median(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE),

           conc_cv_TQC = cv(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE),
           conc_cv_BQC = cv(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE),
           conc_cv_SPL = cv(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE),
           conc_cv_NIST = cv(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE),
           conc_cv_LTR = cv(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE),
           conc_dratio_sd_bqc = sd(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE) / sd(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE),
           conc_dratio_sd_tqc = sd(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE) / sd(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE),
           conc_dratio_mad_bqc = mad(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE) / mad(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE),
           conc_dratio_mad_tqc = mad(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE) / mad(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE)
         )
      d_stats_var_final <- d_stats_var_final |> left_join(d_stats_var_conc, by = grp)
    }



     # If batch medians are requested, calculate the median of all columns (except ID columns) for each feature
     if (use_batch_medians) {
       d_stats_var_final <- d_stats_var_final |>
         summarise(across(-ends_with("_id"), ~ median(.x, na.rm = TRUE)),
                   .by = "feature_id")
     }

     # Convert the processed data back to a regular data.frame from lazy_dt
     d_stats_var_final <- as.data.frame(d_stats_var_final)

     # Identify unique feature IDs present in the dataset
     features_in_dataset <- unique(data@dataset$feature_id)

     # Combine feature details and calculated metrics into a final tibble
     ## `in_data` is a logical flag indicating if the feature is present in both
     ## data and metadata
     data@metrics_qc <- tibble(
       "feature_id" = sort(union(unique(data@dataset_orig$feature_id),
                                 unique(data@annot_features$feature_id)))) |>
       mutate(in_data = .data$feature_id %in% features_in_dataset) |>
       left_join(d_feature_info, by = "feature_id") |>
       left_join(d_method_info, by = "feature_id") |>
       left_join(d_stats_missingval, by = "feature_id") |>
       left_join(d_stats_var_final, by = "feature_id") |>
       relocate("feature_id", "feature_class", "in_data", "valid_feature",
                "is_quantifier", "precursor_mz", "product_mz", "collision_energy")


     # If response curve statistics are to be included and RQC data is available,
     if (!"feature_intensity" %in% names(data@dataset) ) {
       if(isTRUE(include_response_stats)) {
         cli::cli_abort(col_red("Response curve data is missing. Please calculate response curves first using `calculate_response_curve()` function."))
       }
      } else if(is.na(include_response_stats) || include_response_stats){
         d_rqc_stats <- get_response_curve_stats(
           data,
           with_staturation_stats = FALSE,
           limit_to_rqc = TRUE,
           silent_invalid_data = if(isTRUE(include_response_stats)) FALSE else TRUE)
         if(!is.null(d_rqc_stats)){
           data@metrics_qc <- data@metrics_qc |>
             dplyr::left_join(d_rqc_stats, by = "feature_id")

         }
       }


     # If calibration results should be included and metrics_calibration is not
     # empty, join calibration data

     #TODOTODO: IMPORTANT temp solution- align with response curve function/error checking
     if (is.na(include_calibration_results) || include_calibration_results) {

      if(nrow(data@metrics_calibration) > 0){
         data@metrics_qc <- data@metrics_qc |>
           dplyr::left_join(data@metrics_calibration |>
                              dplyr::rename_with(~ str_replace(., "_cal_1", "_cal")) |>
                              dplyr::rename_with(~ str_replace(., "cal_cal", "cal")) |>
                              select(-"is_quantifier", -"fit_cal")
                            , by = "feature_id")
      } else {
        if(isTRUE(include_calibration_results)){
          cli::cli_abort(col_red("Calibration metrics are missing. Please calculate calibration results first using `calculate_calibration()` function."))
         }
      }
    }

     # Return the updated data object with the calculated QC metrics
     data
}


#' Feature Filtering Based on QC Criteria
#'
#' @description
#' Filters a dataset based on quality control (QC) criteria, including intensity,
#' coefficient of variation (CV), signal-to-blank ratios, D-ratio, response curve properties,
#' and proportion of missing values. Criteria apply to different QC types (BQC, TQC) and
#' measurement variables (concentration, intensity, and normalized intensity).
#'
#' To clear all existing filters, run `filter_features_qc()` with `clear_existing = TRUE`
#' and without any filtering criteria.
#' @details
#' This function implements filtering criteria based on quality control (QC) samples
#' and additional analytical parameters, following recommendations outlined by
#' Broadhurst et al. (2018). The implemented criteria evaluate data quality through
#' analysis of QC samples, blanks, and study samples.
#'
#' @references
#' Broadhurst, D., Goodacre, R., Reinke, S. N., Kuligowski, J., Wilson, I. D.,
#' Lewis, M. R., & Dunn, W. B. (2018). Guidelines and considerations for the use
#' of system suitability and quality control samples in mass spectrometry assays
#' applied in clinical studies. *Metabolomics*, 14(6), 72.
#' \doi{10.1007/s11306-018-1367-3}
#'
#' @param data MidarExperiment object.
#' @param clear_existing Logical. If `TRUE`, replaces any existing filters; if `FALSE`, adds new filters on top of existing ones. Default is `TRUE`.
#' @param use_batch_medians Logical. If `TRUE`, uses batch-wise median QC values for filtering. Default is `FALSE`.
#' @param include_qualifier Logical. If `TRUE`, includes qualifier features in the filtering process.
#' @param include_istd Logical. If `TRUE`, includes internal standards (ISTDs) in the filtering process.
#' @param features.to.keep A vector of feature identifiers to retain, even if they do not meet the filtering criteria.
#' @param max.prop.missing.intensity.spl Maximum proportion of missing intensity values among study samples (SPL). Default is `NA`.
#' @param max.prop.missing.normintensity.spl Maximum proportion of missing normalized intensity values among study samples (SPL). Default is `NA`.
#' @param max.prop.missing.conc.spl Maximum proportion of missing concentration values among study samples (SPL). Default is `NA`.
#' @param min.intensity.lowest.bqc Minimum intensity of the lowest BQC sample. Default is `NA`.
#' @param min.intensity.lowest.tqc Minimum intensity of the lowest TQC sample. Default is `NA`.
#' @param min.intensity.lowest.spl Minimum intensity of the lowest study sample (SPL). Default is `NA`.
#' @param min.intensity.median.bqc Minimum median intensity of BQC samples. Default is `NA`.
#' @param min.intensity.median.tqc Minimum median intensity of TQC samples. Default is `NA`.
#' @param min.intensity.median.spl Minimum median intensity of study samples (SPL). Default is `NA`.
#' @param min.intensity.highest.spl Minimum intensity of the highest intensity study sample (SPL). Default is `NA`.
#' @param min.signalblank.median.spl.pblk Minimum signal-to-blank ratio for SPL samples and PBLK. Default is `NA`.
#' @param min.signalblank.median.spl.ublk Minimum signal-to-blank ratio for SPL samples and UBLK. Default is `NA`.
#' @param min.signalblank.median.spl.sblk Minimum signal-to-blank ratio for SPL samples and SBLK. Default is `NA`.
#' @param max.cv.intensity.bqc Maximum CV for intensity in BQC samples. Default is `NA`.
#' @param max.cv.intensity.tqc Maximum CV for intensity in TQC samples. Default is `NA`.
#' @param max.cv.normintensity.bqc Maximum CV for normalized intensity in BQC samples. Default is `NA`.
#' @param max.cv.normintensity.tqc Maximum CV for normalized intensity in TQC samples. Default is `NA`.
#' @param max.cv.conc.bqc Maximum CV for concentration in BQC samples. Default is `NA`.
#' @param max.cv.conc.tqc Maximum CV for concentration in TQC samples. Default is `NA`.
#' @param response.curves.selection Select specific response curves by ID. Default is `NA`.
#' @param response.curves.summary Define the method to summarize multiple response curves. Default is `NA`.
#' @param min.rsquare.response Minimum R-squared value for the response curves. Default is `NA`.
#' @param min.slope.response Minimum slope for the response curve. Default is `NA`.
#' @param max.slope.response Maximum slope for the response curve. Default is `NA`.
#' @param max.yintercept.response Maximum y-intercept of the response curve. Default is `NA`.
#' @param max.dratio.sd.conc.bqc Maximum allowed D-ratio (SD of BQC / SD of SPL) using standard deviation for BQC samples. Default is `NA`.
#' @param max.dratio.sd.conc.tqc Maximum allowed D-ratio (SD of TQC / SD of SPL) using standard deviation for TQC samples. Default is `NA`.
#' @param max.dratio.mad.conc.bqc Maximum allowed D-ratio (MAD of BQC / MAD of SPL) using mean absolute deviation for BQC samples. Default is `NA`.
#' @param max.dratio.mad.conc.tqc Maximum allowed D-ratio (MAD of TQC / MAD of SPL) using mean absolute deviation for TQC samples. Default is `NA`.
#' @param max.dratio.sd.normint.bqc Maximum allowed D-ratio (SD of normalized intensity in BQC / SD of SPL) using standard deviation. Default is `NA`.
#' @param max.dratio.sd.normint.tqc Maximum allowed D-ratio (SD of normalized intensity in TQC / SD of SPL) using standard deviation. Default is `NA`.
#' @param max.dratio.mad.normint.bqc Maximum allowed D-ratio (MAD of normalized intensity in BQC / MAD of SPL) using mean absolute deviation. Default is `NA`.
#' @param max.dratio.mad.normint.tqc Maximum allowed D-ratio (MAD of normalized intensity in TQC / MAD of SPL) using mean absolute deviation. Default is `NA`.
#'
#' @return The input MidarExperiment object with the feature filtering criteria applied.

#' @export
filter_features_qc <- function(data = NULL,
                                    clear_existing = TRUE,
                                    use_batch_medians = FALSE,
                                    include_qualifier,
                                    include_istd,
                                    features.to.keep = NA,
                                    max.prop.missing.intensity.spl = NA,
                                    max.prop.missing.normintensity.spl = NA,
                                    max.prop.missing.conc.spl = NA,
                                    min.intensity.lowest.bqc = NA,
                                    min.intensity.lowest.tqc = NA,
                                    min.intensity.lowest.spl = NA,
                                    min.intensity.median.bqc = NA,
                                    min.intensity.median.tqc = NA,
                                    min.intensity.median.spl = NA,
                                    min.intensity.highest.spl = NA,
                                    min.signalblank.median.spl.pblk = NA,
                                    min.signalblank.median.spl.ublk = NA,
                                    min.signalblank.median.spl.sblk = NA,
                                    max.cv.intensity.bqc = NA,
                                    max.cv.intensity.tqc = NA,
                                    max.cv.normintensity.bqc = NA,
                                    max.cv.normintensity.tqc = NA,
                                    max.cv.conc.bqc = NA,
                                    max.cv.conc.tqc = NA,
                                    response.curves.selection = NA,
                                    response.curves.summary = NA,
                                    min.rsquare.response = NA,
                                    min.slope.response = NA,
                                    max.slope.response = NA,
                                    max.yintercept.response = NA,
                                    max.dratio.sd.conc.bqc = NA,
                                    max.dratio.sd.conc.tqc = NA,
                                    max.dratio.mad.conc.bqc = NA,
                                    max.dratio.mad.conc.tqc = NA,
                                    max.dratio.sd.normint.bqc = NA,
                                    max.dratio.sd.normint.tqc = NA,
                                    max.dratio.mad.normint.bqc = NA,
                                    max.dratio.mad.normint.tqc = NA) {


  check_data(data)

  if (missing(include_qualifier))
    cli::cli_abort(col_red("Argument {.arg include_qualifier} is missing. Please specify whether qualifier features should be included in the filtered data."))

  if (missing(include_istd))
    cli::cli_abort(col_red("Argument {.arg include_istd} is missing. Please specify whether internal standards should be included in the filtered data."))

   # Check if response curve ID is defined when r2 is set. TODO: verify need and extend
  if (all(is.na(response.curves.selection))) {
    if(!is.na(min.rsquare.response) | !is.na(min.slope.response) | !is.na(max.slope.response) |!is.na(max.yintercept.response) | !is.na(response.curves.summary))
      cli::cli_abort(cli::col_red("No response curves selected. Please set the curves using `response.curves.selection`, or remove `response.___` arguments to proceed without response filters."))
  } else {
    if(is.na(min.rsquare.response) & is.na(min.slope.response) & is.na(max.slope.response) & is.na(max.yintercept.response)) {
      cli::cli_abort(cli::col_red("No response filters were defined. Please set the appropriate `response._x_` arguments, remove `response.curves.selection`, or set it to `NA`."))
    } else {
      if (length(unique(response.curves.selection)) > 1){
        if(is.na(response.curves.summary)){
          cli::cli_abort(cli::col_red("Please set `response.curves.summary` to define curve how the results from different curves should be summarized for filtered, Must be either either 'mean', 'median', 'best' or 'worst'."))
        } else {
          rlang::arg_match(response.curves.summary, c("mean", "median", "best", "worst"))
          if (nrow(data@annot_responsecurves) == 0)
            cli::cli_abort(cli::col_red("No response curves are defined in the metadata. Please either remove `response.curves.selection` and any response filters, or reprocess with updated metadata"))
        }
      }
    }
  }

  # Check which criteria categories were defined
  arg_names <- names(as.list(match.call()))
    intensity_criteria_defined <- any(str_detect(arg_names, "[^\\.]intensity|signalblank"))
  norm_intensity_criteria_defined <- any(str_detect(arg_names, "norm_intensity"))
  conc_criteria_defined <- any(str_detect(arg_names, "conc"))
  resp_criteria_defined <- any(str_detect(arg_names, "response"))


  message("Calculating feature QC metrics - please wait...")
  data_local = calc_qc_metrics(data,
                         use_batch_medians = use_batch_medians,
                         include_norm_intensity_stats = data@is_istd_normalized,
                         include_conc_stats = data@is_quantitated,
                         include_response_stats = resp_criteria_defined)


  # Check if feature_ids defind with features.to.keep are present in the dataset
  if (!all(is.na(features.to.keep))) {
    keepers_not_defined <- setdiff(features.to.keep, unique(data@dataset$feature_id))
    txt <- glue::glue_collapse(keepers_not_defined, sep = ", ", last = ", and ")
    if (length(keepers_not_defined) > 0) cli::cli_abort(col_red("Following features defined via `features.to.keep` are not present in this dataset: {txt}"))
  }


  # Save QC filter criteria to MidarExperiment object
  # TODO: fix some of the param below ie. features.to.keep

  # Store the QC filter criteria in MidarExperiment object
  data@parameters_processing <- data@parameters_processing |>
    mutate(
      max.prop.missing.intensity.spl  = max.prop.missing.intensity.spl,
      max.prop.missing.normintensity.spl  = max.prop.missing.normintensity.spl,
      max.prop.missing.conc.spl  = max.prop.missing.conc.spl,
      min.intensity.lowest.bqc = min.intensity.lowest.bqc,
      min.intensity.lowest.tqc = min.intensity.lowest.tqc,
      min.intensity.lowest.spl = min.intensity.lowest.spl,
      min.intensity.median.bqc = min.intensity.median.bqc,
      min.intensity.median.tqc = min.intensity.median.tqc,
      min.intensity.median.spl = min.intensity.median.spl,
      min.intensity.highest.spl = min.intensity.highest.spl,
      max.cv.conc.bqc = max.cv.conc.bqc,
      max.cv.conc.tqc = max.cv.conc.tqc,
      max.cv.intensity.bqc = max.cv.intensity.bqc,
      max.cv.intensity.tqc = max.cv.intensity.tqc,
      max.cv.normintensity.bqc = max.cv.normintensity.bqc,
      max.cv.normintensity.tqc = max.cv.normintensity.tqc,
      max.dratio.sd.conc.bqc = max.dratio.sd.conc.bqc,
      max.dratio.sd.conc.tqc = max.dratio.sd.conc.tqc,
      max.dratio.mad.conc.bqc = max.dratio.mad.conc.bqc,
      max.dratio.mad.conc.tqc = max.dratio.mad.conc.tqc,
      max.dratio.sd.normint.bqc = max.dratio.sd.normint.bqc,
      max.dratio.sd.normint.tqc = max.dratio.sd.normint.tqc,
      min.signalblank.median.spl.pblk = min.signalblank.median.spl.pblk,
      min.signalblank.median.spl.ublk = min.signalblank.median.spl.ublk,
      min.signalblank.median.spl.sblk = min.signalblank.median.spl.sblk,
      response.curves.selection = list(response.curves.selection),
      min.rsquare.response = min.rsquare.response,
      min.slope.response = min.slope.response,
      max.slope.response = max.slope.response,
      max.yintercept.response = max.yintercept.response,
      include_qualifier = include_qualifier,
      include_istd = include_istd,
      features.to.keep = paste(features.to.keep, collapse = ", ")
    )

  metrics_qc_local <- data_local@metrics_qc
  ##tictoc::tic()

  metrics_qc_local <- metrics_qc_local |>
    mutate(
      pass_lod = comp_lgl_vec(
        list(
          compare_values(metrics_qc_local, "intensity_min_BQC", min.intensity.lowest.bqc, ">"),
          compare_values(metrics_qc_local, "intensity_min_TQC", min.intensity.lowest.tqc, ">"),
          compare_values(metrics_qc_local, "intensity_min_SPL", min.intensity.lowest.spl, ">"),
          compare_values(metrics_qc_local, "intensity_median_BQC", min.intensity.median.bqc, ">"),
          compare_values(metrics_qc_local, "intensity_median_TQC", min.intensity.median.tqc, ">"),
          compare_values(metrics_qc_local, "intensity_median_SPL", min.intensity.median.spl, ">"),
          compare_values(metrics_qc_local, "intensity_max_SPL", min.intensity.highest.spl, ">")
        ),
        .operator = "AND"
      ),
      filter_lod = !(is.na(min.intensity.lowest.bqc) & is.na(min.intensity.lowest.tqc) & is.na(min.intensity.lowest.spl) & is.na(min.intensity.median.bqc) & is.na(min.intensity.median.tqc) & is.na(min.intensity.median.spl) & is.na(min.intensity.highest.spl)),

      pass_sb = comp_lgl_vec(
        list(
          compare_values(metrics_qc_local, "sb_ratio_pblk", min.signalblank.median.spl.pblk, ">") | .data$is_istd & !include_istd,
          compare_values(metrics_qc_local, "sb_ratio_ublk", min.signalblank.median.spl.ublk, ">") | .data$is_istd & !include_istd,
          compare_values(metrics_qc_local, "sb_ratio_sblk", min.signalblank.median.spl.sblk, ">") | .data$is_istd & !include_istd
        ),
        .operator = "AND"
      ),
      filter_sb = !(is.na(min.signalblank.median.spl.pblk) & is.na(min.signalblank.median.spl.ublk) &is.na(min.signalblank.median.spl.sblk)),

      pass_cva = comp_lgl_vec(
        list(
          compare_values(metrics_qc_local, "conc_cv_BQC", max.cv.conc.bqc, "<"),
          compare_values(metrics_qc_local, "conc_cv_TQC", max.cv.conc.tqc, "<"),
          compare_values(metrics_qc_local, "norm_intensity_cv_BQC", max.cv.normintensity.bqc, "<"),
          compare_values(metrics_qc_local, "norm_intensity_cv_TQC", max.cv.normintensity.tqc, "<"),
          compare_values(metrics_qc_local, "intensity_cv_BQC", max.cv.intensity.bqc, "<"),
          compare_values(metrics_qc_local, "intensity_cv_TQC", max.cv.intensity.tqc, "<")
        ),
        .operator = "AND"
      ),
      filter_cva = !(is.na(max.cv.conc.bqc) &
                       is.na(max.cv.conc.tqc) &
                       is.na(max.cv.normintensity.bqc) &
                       is.na(max.cv.normintensity.tqc) &
                       is.na(max.cv.intensity.bqc) &
                       is.na(max.cv.intensity.tqc)),

      pass_dratio = comp_lgl_vec(
        list(
          compare_values(metrics_qc_local, "conc_dratio_sd_bqc", max.dratio.sd.conc.bqc, "<"),
          compare_values(metrics_qc_local, "conc_dratio_sd_tqc", max.dratio.sd.conc.tqc, "<"),
          compare_values(metrics_qc_local, "conc_dratio_mad_bqc", max.dratio.mad.conc.bqc, "<"),
          compare_values(metrics_qc_local, "conc_dratio_mad_tqc", max.dratio.mad.conc.tqc, "<"),
          compare_values(metrics_qc_local, "normint_dratio_sd_bqc", max.dratio.sd.normint.bqc, "<"),
          compare_values(metrics_qc_local, "normint_dratio_sd_tqc", max.dratio.sd.normint.tqc, "<"),
          compare_values(metrics_qc_local, "normint_dratio_mad_bqc", max.dratio.mad.normint.bqc, "<"),
          compare_values(metrics_qc_local, "normint_dratio_mad_tqc", max.dratio.mad.normint.tqc, "<")
        ),
        .operator = "AND"
      ),
      filter_dratio = !(is.na(max.dratio.sd.conc.bqc) &
                       is.na(max.dratio.sd.conc.tqc) &
                       is.na(max.dratio.mad.conc.bqc) &
                       is.na(max.dratio.mad.conc.tqc) &
                       is.na(max.dratio.sd.normint.bqc) &
                       is.na(max.dratio.sd.normint.tqc) &
                       is.na(max.dratio.mad.normint.bqc) &
                       is.na(max.dratio.mad.normint.tqc)),


      pass_missingval = comp_lgl_vec(
        list(
          compare_values(metrics_qc_local, "missing_intensity_prop_spl", max.prop.missing.intensity.spl, "<="),
          compare_values(metrics_qc_local, "missing_norm_intensity_prop_spl", max.prop.missing.normintensity.spl, "<="),
          compare_values(metrics_qc_local, "missing_conc_prop_spl", max.prop.missing.conc.spl, "<=")
        ),
        .operator = "AND"
      ),
      filter_missingval = !(is.na(max.prop.missing.intensity.spl) &
                          is.na(max.prop.missing.normintensity.spl) &
                          is.na(max.prop.missing.conc.spl))
    )

  ##tictoc::toc()
  # Check if linearity criteria are defined
  metrics_qc_local<- metrics_qc_local |> mutate(pass_linearity = NA)

  if (resp_criteria_defined){
    if (is.numeric(response.curves.selection)) {

      rqc_r2_col_names <- names(metrics_qc_local)[which(stringr::str_detect(names(metrics_qc_local), "r2_rqc"))]
      rqc_r2_col <-  rqc_r2_col_names[response.curves.selection]
      rqc_slope_col_names <- names(metrics_qc_local)[which(stringr::str_detect(names(metrics_qc_local), "slopenorm_rqc"))]
      rqc_slope_col <-  rqc_slope_col_names[response.curves.selection]
      rqc_y0_col_names <- names(metrics_qc_local)[which(stringr::str_detect(names(metrics_qc_local), "y0norm_rqc"))]
      rqc_y0_col <-  rqc_y0_col_names[response.curves.selection]

      if (any(is.na(rqc_r2_col))) cli::cli_abort(cli::col_red("The specified response curve index exceeds the available range. There are only {length(rqc_r2_col_names)} response curves in the dataset. Please adjust the indices set via `response.curves.selection`"))

    } else if(is.character(response.curves.selection)){
        rqc_r2_col <- paste0("r2_rqc_", response.curves.selection)
        rqc_slope_col <- paste0("slopenorm_rqc_", response.curves.selection)
        rqc_y0_col <- paste0("y0norm_rqc_", response.curves.selection)
        missing_curves <- response.curves.selection[!(rqc_r2_col %in% names(metrics_qc_local))]
        if (length(missing_curves) > 0) {
          cli::cli_abort(cli::col_red("The following response curves are not defined in the metadata: {paste(missing_curves, collapse=', ')}. Please adjust the curve ids or ensure correct identifiers."))
        }
    } else {
      cli::cli_abort(cli::col_red("The `response.curves.selection` must be specified as either numeric indices or identifiers provided as strings, corresponding to the response curve(s)."))
    }

    # Summarize metrics across curves (columns) based on specified criteria

    # Determine the function to apply based on response.curves.summary
    fun_r2 <- case_match(
      response.curves.summary,
      "worst" ~ "safe_min",
      "best" ~ "safe_max",
      .default = response.curves.summary
    )

    fun_slope_min <- case_match(
      response.curves.summary,
      c("worst") ~ "safe_min",
      c("best") ~ "safe_max",
      .default = response.curves.summary
    )

    fun_slope_max <- case_match(
      response.curves.summary,
      c("worst") ~ "safe_max",
      c("best") ~ "safe_min",
      .default = response.curves.summary
    )

    fun_y0 <- case_match(
      response.curves.summary,
      c("worst") ~ "safe_max",
      c("best") ~ "safe_min",
      .default = response.curves.summary
    )

    # Calculate summary metrics using purrr::pmap_dbl
    metrics_qc_local <- metrics_qc_local |>
      mutate(
        rqc_r2__sum__ = purrr::pmap_dbl(
          across(all_of(rqc_r2_col)), ~ do.call(fun_r2, list(c(...), na.rm = TRUE))
        ),
        rqc_slope__sum__min__ = purrr::pmap_dbl(
          across(all_of(rqc_slope_col)), ~ do.call(fun_slope_min, list(c(...), na.rm = TRUE))
        ),
        rqc_slope__sum__max__ = purrr::pmap_dbl(
          across(all_of(rqc_slope_col)), ~ do.call(fun_slope_max, list(c(...), na.rm = TRUE))
        ),
        rqc_y0__sum__ = purrr::pmap_dbl(
          across(all_of(rqc_y0_col)), ~ do.call(fun_y0, list(c(...), na.rm = TRUE))
        )
      )

    # Check if columns exist before mutating pass_linearity
    if (all(rqc_r2_col %in% names(metrics_qc_local))) {
      metrics_qc_local <- metrics_qc_local |>
        mutate(
          pass_linearity = if_else(
            !is.na(.data$rqc_r2__sum__) | (!is.na(min.slope.response) & !is.na(max.slope.response) & !is.na(max.yintercept.response)),
            (.data$rqc_r2__sum__ > min.rsquare.response | is.na(min.rsquare.response)) &
              (.data$rqc_slope__sum__min__ > min.slope.response | is.na(min.slope.response)) &
              (.data$rqc_slope__sum__max__ <= max.slope.response | is.na(max.slope.response)) &
              (.data$rqc_y0__sum__ < max.yintercept.response | is.na(max.yintercept.response)),
            NA
          ),
          filter_linearity = !(is.na(min.rsquare.response) &
                                is.na(min.slope.response) &
                                is.na(max.slope.response) &
                                is.na(max.yintercept.response))
        )
    }

  }


  #TODO: change name of pass_qualifer to a better name
  metrics_qc_local <- metrics_qc_local |>
    mutate(
      pass_istd = !.data$is_istd | include_istd,
      pass_qualifier = .data$is_quantifier | include_qualifier,
      pass_featureskeep = .data$feature_id %in% features.to.keep,
      filter_istd = TRUE,
      filter_qualifier = TRUE,
      filter_featureskeep = !all(is.na(features.to.keep)) && length(features.to.keep) > 0
    )



  # Check if filter has been previously set and if it should be overwritten

  if(!clear_existing && all("all_filter_pass" %in% names(data@metrics_qc))){

    metrics_old <- data@metrics_qc |>
      select(
        any_of(c(
          "feature_id",
          "batch_id",
          qc_pass_before = "all_filter_pass",
          pass_lod_before = "pass_lod",
          pass_sb_before = "pass_sb",
          pass_cva_before = "pass_cva",
          pass_linearity_before = "pass_linearity",
          pass_dratio_before =  "pass_dratio",
          pass_missingval_before = "pass_missingval",
          pass_istd_before = "pass_istd",
          pass_qualifier_before = "pass_qualifier",
          pass_featureskeep_before = "pass_featureskeep",
          filter_lod_before  = "filter_lod",
          filter_sb_before  = "filter_sb",
          filter_cva_before  = "filter_cva",
          filter_dratio_before  = "filter_dratio",
          filter_linearity_before  = "filter_linearity",
          filter_missingval_before  = "filter_missingval",
          filter_istd_before = "filter_istd",
          filter_qualifier_before = "filter_qualifier",
          filter_featureskeep_before = "filter_featureskeep"
        ))
      )

    metrics_qc_local <- metrics_qc_local |> full_join(metrics_old, by = "feature_id")

    prev_filters <- list()

    if(all(metrics_old$filter_missingval_before )) {
      if (all(metrics_qc_local$filter_missingval)) {
        prev_filters <- append(prev_filters, "Missing Values")
      } else {
       metrics_qc_local$pass_missingval <- metrics_qc_local$pass_missingval_before
       metrics_qc_local$filter_missingval <- metrics_qc_local$filter_missingval_before
      }
  }


    if(all(metrics_qc_local$filter_lod) && all(metrics_old$filter_lod_before )){
      prev_filters <- append(prev_filters, "Min-Intensity")
    } else {
      metrics_qc_local$pass_lod <- metrics_qc_local$pass_lod_before
      metrics_qc_local$filter_lod <- metrics_qc_local$filter_lod_before
    }

    if(all(metrics_old$filter_sb_before )){
      if(all(metrics_qc_local$filter_sb)){
        prev_filters <- append(prev_filters, "Signal-to-Blank")
      } else {
        metrics_qc_local$pass_sb <- metrics_qc_local$pass_sb_before
        metrics_qc_local$filter_sb <- metrics_qc_local$filter_sb_before
      }
    }

    if(all(metrics_old$filter_cva_before )) {
      if(all(metrics_qc_local$filter_cva)){
        prev_filters <- append(prev_filters, "%CV")
      } else {
        metrics_qc_local$pass_cva <- metrics_qc_local$pass_cva_before
        metrics_qc_local$filter_cva <- metrics_qc_local$filter_cva_before
      }
    }

    if(all(metrics_old$filter_dratio_before )) {
      if(all(metrics_qc_local$filter_dratio)) {
        prev_filters <- append(prev_filters, "D-ratio")
      } else {
        metrics_qc_local$pass_dratio <- metrics_qc_local$pass_dratio_before
        metrics_qc_local$filter_dratio <- metrics_qc_local$filter_dratio_before
      }
    }

    if ("filter_linearity_before" %in% names(metrics_qc_local)) {
      if(all(metrics_old$filter_linearity_before )) {
        if(all(metrics_qc_local$filter_linearity)){
          prev_filters <- append(prev_filters, "Linearity")
        } else {
          metrics_qc_local$pass_linearity <- metrics_qc_local$pass_linearity_before
          metrics_qc_local$filter_linearity <- metrics_qc_local$filter_linearity_before
        }
      }
    }

    # below filter are always defined (never NA). When they are changed in new filter application, they are always replaced a
    # and a message appears


    if ("filter_istd" %in% names(metrics_qc_local)) {
      if(!all(metrics_qc_local$pass_istd == metrics_old$pass_istd_before)){
        prev_filters <- append(prev_filters, "ISTD")
      }
    }

    if ("filter_qualifier" %in% names(metrics_qc_local)) {
        if(!all(metrics_qc_local$pass_qualifier == metrics_old$pass_qualifier_before)){
          prev_filters <- append(prev_filters, "Qualifier")
        }
    }

    if ("filter_featureskeep" %in% names(metrics_qc_local)) {
        if(!all(metrics_qc_local$pass_featureskeep == metrics_old$pass_featureskeep_before)){
          prev_filters <- append(prev_filters, "Keepers")
      }
    }


    if(length(prev_filters) > 0){
      cli::cli_alert_warning(cli::col_yellow(glue::glue("Replaced following previously defined QC filters: {glue::glue_collapse(prev_filters, sep = ', ', last = ', and ')}")))
    }
  }

  metrics_qc_local <- metrics_qc_local |>
      mutate(
        all_qc_filter_pass =
          (
            (is.na(.data$pass_lod) | .data$pass_lod) &
              (is.na(.data$pass_sb) | .data$pass_sb) &
              (is.na(.data$pass_cva) | .data$pass_cva) &
              (is.na(.data$pass_linearity) | .data$pass_linearity) &
              (is.na(.data$pass_dratio) | .data$pass_dratio) &
              (is.na(.data$pass_missingval) | .data$pass_missingval) &
              (is.na(.data$valid_feature) | .data$valid_feature) &
              (is.na(.data$pass_istd) | .data$pass_istd) &
              (is.na(.data$pass_qualifier) | .data$pass_qualifier)
          )
      )

  n_featured_forcedkeep <-  intersect(metrics_qc_local[!metrics_qc_local$all_qc_filter_pass,]$feature_id, features.to.keep)
  if(length(n_featured_forcedkeep) > 0){
    cli::cli_alert_warning(cli::col_yellow(glue::glue("The following features were forced to be retained despite not meeting filtering criteria: {glue::glue_collapse(n_featured_forcedkeep, sep = ', ', last = ', and ')}")))
  }

  metrics_qc_local <- metrics_qc_local |>
    mutate(
      all_filter_pass = .data$all_qc_filter_pass | (is.na(.data$pass_featureskeep) | .data$pass_featureskeep)
    )

    #TODO: deal with invalid integrations (as defined by user in metadata)
    d_filt <- metrics_qc_local |>
      filter(.data$all_filter_pass)

    d_metrics_temp <- metrics_qc_local

    if(!include_qualifier) d_metrics_temp <- d_metrics_temp |> filter(.data$is_quantifier)
    if(!include_istd) d_metrics_temp <- d_metrics_temp |> filter(!.data$is_istd)

    n_all_quant <- get_feature_count(data, is_istd = include_istd && NA, is_quantifier = TRUE)
    n_all_qual <- get_feature_count(data, is_istd = include_istd && NA, is_quantifier = FALSE)

    n_filt_quant <- nrow(d_metrics_temp |>  filter(.data$in_data, .data$is_quantifier, .data$all_filter_pass))
    n_filt_qual <- nrow(d_metrics_temp |>  filter(.data$in_data, !.data$is_quantifier, .data$all_filter_pass))

    if(!clear_existing && ("all_filter_pass" %in% names(data@metrics_qc))){
      n_filt_quant_before <- nrow(d_metrics_temp |>  filter(.data$in_data, .data$is_quantifier, .data$qc_pass_before))
      n_filt_qual_before <- nrow(d_metrics_temp |>  filter(.data$in_data,!.data$is_quantifier, .data$qc_pass_before))
      qc_pass_prev <- sum(d_metrics_temp$qc_pass_before, na.rm = TRUE)
    }

    n_istd_quant <- get_feature_count(data, is_istd = TRUE, is_quantifier = TRUE)
    n_istd_qual <- get_feature_count(data, is_istd = TRUE, is_quantifier = FALSE)

    if (!include_istd) {
      n_filt_quant <- nrow(d_filt |>  filter(.data$in_data,!.data$is_istd, .data$is_quantifier))
      n_filt_qual <- nrow(d_filt |>  filter(.data$in_data, !.data$is_istd, !.data$is_quantifier))

      if(!clear_existing && all("all_filter_pass" %in% names(data@metrics_qc))){
        n_filt_quant_before <- nrow(d_metrics_temp|>  filter(.data$in_data,!.data$is_istd, .data$is_quantifier, .data$qc_pass_before))
        n_filt_qual_before <- nrow(d_metrics_temp |>  filter(.data$in_data,!.data$is_istd, !.data$is_quantifier, .data$qc_pass_before))
      }
    }

    qc_pass_now <- sum(metrics_qc_local$all_filter_pass, na.rm = TRUE)


  filter_cleared <- !any(str_detect(arg_names[!arg_names %in% c("include_istd", "include_qualifier")], "\\."))


  if(include_qualifier){
    if(!clear_existing && all("all_filter_pass" %in% names(data@metrics_qc))){
        cli::cli_alert_success(cli::col_green("\rFeature QC filters were updated: {n_filt_quant} (before {n_filt_quant_before}) of {n_all_quant} quantifier and {n_filt_qual} of {n_all_qual} qualifier features meet QC criteria ({if_else(!include_istd, 'not including the', 'including the')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)"))
    } else {
      if (filter_cleared){
        cli::cli_alert_success(cli::col_green("\r{.strong {.emph Cleared}} all feature QC filters! All {n_all_quant} quantifier and all {n_all_qual} qualifier features are now selected ({if_else(!include_istd, 'not including the', 'including the')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)"))
      } else {
        cli::cli_alert_success(cli::col_green("\rNew feature QC filters were defined: {n_filt_quant} of {n_all_quant} quantifier and {n_filt_qual} of {n_all_qual} qualifier features meet QC criteria ({if_else(!include_istd, 'not including the', 'including the')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)"))
      }
    }
  } else {
    if(!clear_existing && all("all_filter_pass" %in% names(data@metrics_qc))){
      cli::cli_alert_success(cli::col_green("\rFeature QC filters were updated: {n_filt_quant} (before {n_filt_quant_before}) of {n_all_quant} quantifier features meet QC criteria ({if_else(!include_istd, 'not including the', 'including the')} {n_istd_quant} quantifier ISTD features)"))
    } else {
      if(!filter_cleared)
        cli::cli_alert_success(cli::col_green("\rNew feature QC filters were defined: {n_filt_quant} of {n_all_quant} quantifier features meet QC criteria ({if_else(!include_istd, 'not including the', 'including the')} {n_istd_quant} quantifier ISTD features)."))
      else
        cli::cli_alert_success(cli::col_green("\r{.strong {.emph Cleared}} all feature QC filters! All {n_all_quant} quantifier features are now selected ({if_else(!include_istd, 'not including the', 'including the')} {n_istd_quant} quantifier ISTD features)."))
      }
  }

  # TODO: cleanup
  #if (!include_qualifier) d_filt <- d_filt |> filter(.data$is_quantifier)
  #if (!include_istd) d_filt <- d_filt |> filter(!.data$is_istd)

  data@is_filtered <- TRUE
  data@metrics_qc <- metrics_qc_local |> select(-dplyr::ends_with("before"))


  data@dataset_filtered <- data@dataset |>
    dplyr::right_join(d_filt |> filter(.data$all_filter_pass) |> dplyr::select("feature_id"), by = "feature_id")

  data
}
