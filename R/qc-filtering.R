#' Calculate feature quality control (QC) metrics
#' @description
#' For each feature, different QC metrics are calculated, either across the
#' full analysis or as medians of batch-wise calculations
#'
#' @param data MidarExperiment object
#' @param batch_medians Use median of batch-wise derived QC values. Default is FALSE.
#' @param with_norm_intensity Include normalized intensity statistics of features. Default is TRUE.
#' @param with_conc Include concentration statistics of features. Default is TRUE.
#' @param with_linearity Include linearity statistics of response curves. This will increase the calculation time. Default is TRUE.
#' @return MidarExperiment object
#' @export
qc_calc_metrics <- function(data = NULL, batch_medians, with_norm_intensity = TRUE, with_conc = TRUE, with_linearity = TRUE ) {
  check_data(data)
   # TODO: remove later when fixed
  if (tolower(data@analysis_type) == "lipidomics") data <- lipidomics_get_lipid_class_names(data)
        # All features defined in the metadata
      d_feature_info <- data@annot_features |>
        select("valid_feature", "feature_id", "feature_class", "is_istd", "is_quantifier", "istd_feature_id", "quant_istd_feature_id", "response_factor")

      # All features MS method parameters
      method_var <- c("method_precursor_mz", "method_product_mz", "method_collision_energy")
      d_method_template <- tibble("precursor_mz" = NA_character_, "product_mz" = NA_character_, "collision_energy" = NA_character_)
      if(any(c(method_var) %in% names(data@dataset_orig))){

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
          dplyr::mutate(across(where(\(x) {is.numeric(suppressWarnings(as.numeric(x)))}) & !.data$feature_id, as.numeric)) #TODO: more elegant

      } else {
        d_method_info <- d_method_template
      }

      d_stats_missingval <- data@dataset |>
        dplyr::filter(.data$qc_type %in% c("SPL", "NIST", "LTR", "BQC", "TQC", "PBLK", "SBLK", "UBLK", "IBLK", "CAL", "STD", "LQC", "MQC", "UQC")) |>
        dplyr::summarise(
          .by = "feature_id",
          missing_intensity_prop_spl = if ("feature_intensity" %in% names(data@dataset )) sum(is.na(.data$feature_intensity[.data$qc_type == "SPL"]))/length(.data$feature_intensity[.data$qc_type == "SPL"]) else NA_real_,
          missing_norm_intensity_prop_spl = if ("feature_norm_intensity" %in% names(data@dataset )) sum(is.na(.data$feature_norm_intensity[.data$qc_type == "SPL"]))/length(.data$feature_norm_intensity[.data$qc_type == "SPL"]) else NA_real_,
          missing_conc_prop_spl = if ("feature_conc" %in% names(data@dataset )) sum(is.na(.data$feature_conc[.data$qc_type == "SPL"]))/length(.data$feature_conc[.data$qc_type == "SPL"]) else NA_real_,
          na_in_all = if ("feature_intensity" %in% names(data@dataset )) all(is.na(.data$feature_intensity)) else NA_real_
        )

      if (batch_medians) grp <- c("feature_id", "batch_id") else grp <- c("feature_id")
      d_stats_var <- data@dataset |>
        select(any_of(c("batch_id", "feature_id", "qc_type","feature_intensity", "feature_conc", "feature_norm_intensity"))) |>
        filter(.data$qc_type != "RQC") |>
        mutate(qc_type = factor(.data$qc_type),
               batch_id = factor(.data$batch_id))

     # Using dtplyr to improve speed (2x), group by batch still 5x times slower than no batch grouping
     d_stats_var <- suppressMessages(dtplyr::lazy_dt(d_stats_var))
     d_stats_var_int <-  d_stats_var |>
        summarise(
          .by = grp,
          intensity_min_SPL = safe_min(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          intensity_max_SPL = safe_max(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          intensity_min_BQC = safe_min(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
          intensity_min_TQC = safe_min(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
          intensity_median_PBLK = median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE),
          intensity_median_SPL = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          intensity_median_BQC = median(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
          intensity_median_TQC = median(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
          intensity_median_NIST = median(.data$feature_intensity[.data$qc_type == "NIST"], na.rm = TRUE),
          intensity_median_LTR = median(.data$feature_intensity[.data$qc_type == "LTR"], na.rm = TRUE),

          intensity_cv_TQC = cv(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE)  * 100,
          intensity_cv_BQC = cv(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE)  * 100,
          intensity_cv_SPL = cv(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE)  * 100,
          intensity_cv_LTR = cv(.data$feature_intensity[.data$qc_type == "LTR"], na.rm = TRUE)  * 100,
          intensity_cv_NIST = cv(.data$feature_intensity[.data$qc_type == "NIST"], na.rm = TRUE)  * 100,

          sb_ratio_q10_pbk = quantile(.data$feature_intensity[.data$qc_type == "SPL"], probs = 0.1, na.rm = TRUE, names = FALSE) / median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE, names = FALSE),
          sb_ratio_pblk = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE, names = FALSE) / median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE, names = FALSE),
          sb_ratio_ublk = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE, names = FALSE) / median(.data$feature_intensity[.data$qc_type == "UBLK"], na.rm = TRUE, names = FALSE),
          sb_ratio_sblk = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE, names = FALSE) / median(.data$feature_intensity[.data$qc_type == "SBLK"], na.rm = TRUE, names = FALSE),
        )

     d_stats_var_final <- d_stats_var_int

     if(data@is_istd_normalized){
       d_stats_var_norm_int <-  d_stats_var |>
         summarise(
           .by = grp,
            norm_intensity_cv_TQC = cv(.data$feature_norm_intensity[.data$qc_type == "TQC"], na.rm = TRUE) * 100,
            norm_intensity_cv_BQC = cv(.data$feature_norm_intensity[.data$qc_type == "BQC"], na.rm = TRUE) * 100,
            norm_intensity_cv_SPL = cv(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
            norm_intensity_cv_TQC = cv(.data$feature_norm_intensity[.data$qc_type == "TQC"], na.rm = TRUE)  * 100,
            norm_intensity_cv_BQC = cv(.data$feature_norm_intensity[.data$qc_type == "BQC"], na.rm = TRUE)  * 100,
          )
       d_stats_var_final <- d_stats_var_final |> left_join(d_stats_var_norm_int, by = grp)

     }

     if(data@is_quantitated){
        d_stats_var_conc <-  d_stats_var |>
          summarise(
            .by = grp,
             conc_median_TQC = median(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE),
             conc_median_BQC = median(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE),
             conc_median_SPL = median(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE),
             conc_median_NIST = median(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE),
             conc_median_LTR = median(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE),

             conc_cv_TQC = cv(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE) * 100,
             conc_cv_BQC = cv(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE)  * 100,
             conc_cv_SPL = cv(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE)  * 100,
             conc_cv_NIST = cv(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE)  * 100,
             conc_cv_LTR = cv(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE)  * 100,
             conc_dratio_sd_bqc_conc = sd(.data$feature_conc[.data$qc_type == "BQC"]) / sd(.data$feature_conc[.data$qc_type == "SPL"]),
             conc_dratio_sd_tqc_conc = sd(.data$feature_conc[.data$qc_type == "TQC"]) / sd(.data$feature_conc[.data$qc_type == "SPL"]),
             conc_dratio_sd_bqc_normint = sd(.data$feature_norm_intensity[.data$qc_type == "BQC"]) / sd(.data$feature_norm_intensity[.data$qc_type == "SPL"]),
             conc_dratio_sd_tqc_normint = sd(.data$feature_norm_intensity[.data$qc_type == "TQC"]) / sd(.data$feature_norm_intensity[.data$qc_type == "SPL"]),
             conc_dratio_mad_bqc_conc = mad(.data$feature_conc[.data$qc_type == "BQC"]) / mad(.data$feature_conc[.data$qc_type == "SPL"]),
             conc_dratio_mad_tqc_conc = mad(.data$feature_conc[.data$qc_type == "TQC"]) / mad(.data$feature_conc[.data$qc_type == "SPL"]),
             conc_dratio_mad_bqc_normint = mad(.data$feature_norm_intensity[.data$qc_type == "BQC"]) / mad(.data$feature_norm_intensity[.data$qc_type == "SPL"]),
             conc_dratio_mad_tqc_normint = mad(.data$feature_norm_intensity[.data$qc_type == "TQC"]) / mad(.data$feature_norm_intensity[.data$qc_type == "SPL"])
      )
        d_stats_var_final <- d_stats_var_final |> left_join(d_stats_var_conc, by = grp)
     }


    if (batch_medians){
      d_stats_var_final <- d_stats_var_final |>
        summarise(across(-ends_with("_id"), ~ median(.x, na.rm = TRUE)), .by= "feature_id")
    }
    # Convert to lazy_dt to data.frame
     d_stats_var_final <-  as.data.frame(d_stats_var_final)

    features_in_dataset <- unique(data@dataset$feature_id)

    # Combine details of features and metrics into one tibble
    ## `in_data` is a logical flag indicating if the feature is present in the dataset (= data and metadata)

    data@metrics_qc <-
      tibble("feature_id" = sort(union(unique(data@dataset_orig$feature_id),
                                     unique(data@annot_features$feature_id)))) |>
      mutate(in_data = .data$feature_id %in% features_in_dataset) |>
      left_join(d_feature_info, by = "feature_id") |>
      left_join(d_method_info, by = "feature_id") |>
      left_join(d_stats_missingval, by = "feature_id") |>
      left_join(d_stats_var_final, by = "feature_id") |>
      relocate("feature_id", "feature_class", "in_data", "valid_feature", "is_quantifier", "precursor_mz", "product_mz", "collision_energy")

  if (with_linearity & "RQC" %in% data@dataset$qc_type) {
    d_rqc_stats <- get_response_curve_stats(data,  with_staturation_stats = FALSE, limit_to_rqc = TRUE)
    data@metrics_qc <- data@metrics_qc |>
      dplyr::left_join(d_rqc_stats, by = "feature_id")
  }

  data
}


#' Linear regression statistics of response curves

#' @details
#' When using `with_staturation_stats` then the `lancer` package needs to be installed
#'

#' @param data MidarExperiment object
#' @param with_staturation_stats Include statistics and classification from the `lancer` package.
#' @param limit_to_rqc If TRUE (default) then only includes RQC qc type
#' @return Tibble with linear regression stats of all curves in a wide format
#' @export

get_response_curve_stats <- function(data = NULL, with_staturation_stats = FALSE, limit_to_rqc = FALSE) {
  check_data(data)
  get_lm_results <- function(data){

    dt <- data

    # Scale to max x and max y (1 = max)
    dt$x_scaled <- (dt$analyzed_amount) / (max(dt$analyzed_amount))
    dt$y_scaled <- (dt$feature_intensity) / (max(dt$feature_intensity))
    tryCatch(
      {
        res <- lm(y_scaled ~ x_scaled, data = dt, na.action = na.exclude)
        r.squared <- summary(res)$r.squared
        slope <- res$coefficients[[2]]
        intercept <- res$coefficients[1]
        return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], r.squared = r.squared , slope = slope, intercept = intercept))

      },
      error = function(e) {
        return(list(feature_id = dt$feature_id[1], curve_id = dt$curve_id[1], r.squared = NA_real_ , slope = NA_real_, intercept = NA_real_))

        }
    )

  }


  d_stats  <- data@dataset |>
    select("analysis_id", "feature_id", "feature_intensity") |>
    dplyr::inner_join(data@annot_responsecurves, by = "analysis_id") |>
    dplyr::filter(!all(is.na(.data$feature_intensity))) |>
    dplyr::group_split(.data$feature_id, .data$curve_id)

    d_stats <- map(d_stats, function(x) get_lm_results(x))

    d_stats <- d_stats |> bind_rows() |>
      dplyr::mutate(slopenorm = .data$slope,
                    y0norm = .data$intercept) |>
      dplyr::select("feature_id", "curve_id", r2 = "r.squared", "slopenorm", "y0norm") |>
      tidyr::pivot_wider(names_from = "curve_id", values_from = c("r2", "slopenorm", "y0norm"), names_prefix = "rqc_")


  if (with_staturation_stats){
    print("with_staturation_stats")
    if (!requireNamespace("lancer", quietly = TRUE)) {
      stop(
        "Package `lancer` must be installed when `with_staturation_stats = TRUE`. It is available from `https://github.com/SLINGhub/lancer`",
        call. = FALSE
      )
    }
    d_stats_lancer <- data@dataset |>
      select("analysis_id", "feature_id", "feature_intensity")
      dplyr::inner_join(data@annot_responsecurves, by = "analysis_id") |>
      dplyr::group_by(.data$feature_id, .data$curve_id) |>
      dplyr::filter(!all(is.na(.data$feature_intensity))) |>
      tidyr::nest() |>
      mutate(
        lancer_raw = map(data, \(x) lancer::summarise_curve_data(x, "analyzed_amount", "feature_intensity")),
        lancer = map(.data$lancer_raw, \(x) lancer::evaluate_linearity(x))
      ) |>
      select(-"lancer_raw") |>
      tidyr::unnest(c("lancer")) |>
      dplyr::select("feature_id", "curve_id", "r_corr", class_wf2 = "wf2_group",  "pra_linear", "mandel_p_val", "concavity") |>
      tidyr::pivot_wider(names_from = "curve_id", values_from = c("r_corr", "class_wf2", "pra_linear", "mandel_p_val", "concavity"), names_prefix = "rqc_") |>
      ungroup()

    d_stats <- d_stats |> left_join(d_stats_lancer, by = c("feature_id"))
  }
  d_stats
}
#' Feature Filtering Based on Quality Control Criteria
#'
#' @description
#' Filters a dataset based on quality control (QC) criteria, including intensity,
#' coefficient of variation (CV), signal-to-blank ratios, D-ratio, response curve properties,
#' and proportion of missing values. Criteria apply to different QC types (BQC, TQC) and
#' measurement variables (concentration, intensity, and normalized intensity).
#'
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
#' @param replace_existing Logical. If `TRUE`, replaces any existing filters; if `FALSE`, adds new filters on top of existing ones. Default is `TRUE`.
#' @param batch_medians Logical. If `TRUE`, uses batch-wise median QC values for filtering. Default is `FALSE`.
#' @param qualifier.include Logical. If `TRUE`, includes qualifier features in the filtering process. Default is `FALSE`.
#' @param istd.include Logical. If `TRUE`, includes internal standards (ISTDs) in the filtering process. Default is `FALSE`.
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
#' @param response.curves.select Select specific response curves by ID. Default is `NA`.
#' @param response.curves.summary Define the method to summarize multiple response curves. Default is `NA`.
#' @param min.rsquare.response Minimum R-squared value for the response curves. Default is `NA`.
#' @param min.slope.response Minimum slope for the response curve. Default is `NA`.
#' @param max.slope.response Maximum slope for the response curve. Default is `NA`.
#' @param max.yintercept.response Maximum y-intercept of the response curve. Default is `NA`.
#' @param max.dratio.sd.bqc Maximum allowed D-ratio (SD of BQC / SD of SPL) using standard deviation for BQC samples. Default is `NA`.
#' @param max.dratio.sd.tqc Maximum allowed D-ratio (SD of TQC / SD of SPL) using standard deviation for TQC samples. Default is `NA`.
#' @param max.dratio.mad.bqc Maximum allowed D-ratio (MAD of BQC / MAD of SPL) using mean absolute deviation for BQC samples. Default is `NA`.
#' @param max.dratio.mad.tqc Maximum allowed D-ratio (MAD of TQC / MAD of SPL) using mean absolute deviation for TQC samples. Default is `NA`.
#' @param max.dratio.sd.normint.bqc Maximum allowed D-ratio (SD of normalized intensity in BQC / SD of SPL) using standard deviation. Default is `NA`.
#' @param max.dratio.sd.normint.tqc Maximum allowed D-ratio (SD of normalized intensity in TQC / SD of SPL) using standard deviation. Default is `NA`.
#' @param max.dratio.mad.normint.bqc Maximum allowed D-ratio (MAD of normalized intensity in BQC / MAD of SPL) using mean absolute deviation. Default is `NA`.
#' @param max.dratio.mad.normint.tqc Maximum allowed D-ratio (MAD of normalized intensity in TQC / MAD of SPL) using mean absolute deviation. Default is `NA`.
#'
#' @return The input MidarExperiment object with the feature filtering criteria applied.

#' @export
filter_features_qc <- function(data = NULL,
                                    replace_existing = TRUE,
                                    batch_medians = FALSE,
                                    qualifier.include = FALSE,
                                    istd.include = FALSE,
                                    features.to.keep = NULL,
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
                                    response.curves.select = NA,
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
   # Check if RQC curve ID is defined when r2 is set
  if (all(is.na(response.curves.select))) {
    if(!is.na(min.rsquare.response) | !is.na(min.slope.response) | !is.na(max.slope.response) |!is.na(max.yintercept.response) | !is.na(response.curves.summary))
      cli::cli_abort(cli::col_red("No response curves selected. Please set the curves using `response.curves.select`, or remove `response.___` arguments to proceed without response filters."))
  } else {
    if(is.na(min.rsquare.response) & is.na(min.slope.response) & is.na(max.slope.response) & is.na(max.yintercept.response)) {
      cli::cli_abort(cli::col_red("No response filters were defined. Please set the appropriate `response.___` arguments, remove `response.curves.select`, or set it to NA."))
    } else {
      if (length(unique(response.curves.select)) > 1){
        if(is.na(response.curves.summary)){
          cli::cli_abort(cli::col_red("Please define `response.curves.summary` to define how the results from different curves should be summarized for filtered, Must be either o either 'mean', 'median', 'best' or 'worst'."))
        } else {
          rlang::arg_match(response.curves.summary, c("mean", "median", "best", "worst"))
          if (nrow(data@annot_responsecurves) == 0)
            cli::cli_abort(cli::col_red("No response curves are defined in the metadata. Please either remove `response.curves.select` and any response filters, or reprocess with updated metadata"))
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


  cat("Calculating feature QC metrics - please wait...")
  data_local = qc_calc_metrics(data,
                         batch_medians = batch_medians,
                         with_norm_intensity = data@is_istd_normalized,
                         with_conc = data@is_quantitated,
                         with_linearity = resp_criteria_defined)


  # Check if feature_ids defind with features.to.keep are present in the dataset
  if (!is.null(features.to.keep)) {
    keepers_not_defined <- setdiff(features.to.keep, unique(data@dataset$feature_id))
    txt <- glue::glue_collapse(keepers_not_defined, sep = ", ", last = ", and ")
    if (length(keepers_not_defined) > 0) cli::cli_abort(glue::glue("Following features defined via `features.to.keep` are not present in this dataset: {txt}"))
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
      response.curves.select_used_for_filt = list(response.curves.select),
      min.rsquare.response = min.rsquare.response,
      min.slope.response = min.slope.response,
      max.slope.response = max.slope.response,
      max.yintercept.response = max.yintercept.response,
      qualifier.include = qualifier.include,
      istd.include = istd.include,
      features.to.keep = NA
    )

  metrics_qc_local <- data_local@metrics_qc
  ##tictoc::tic()
  metrics_qc_local <- metrics_qc_local |>
    #rowwise() |>  #TODO make it vevtorized
    mutate(
      pass_lod = comp_lgl_vec(
        list(comp_val(dplyr::cur_data(), "intensity_min_BQC", min.intensity.lowest.bqc, ">"),
        comp_val(dplyr::cur_data(), "intensity_min_TQC", min.intensity.lowest.tqc, ">"),
        comp_val(dplyr::cur_data(), "intensity_min_SPL", min.intensity.lowest.spl, ">"),
        comp_val(dplyr::cur_data(), "intensity_median_BQC", min.intensity.median.bqc, ">"),
        comp_val(dplyr::cur_data(), "intensity_median_TQC", min.intensity.median.tqc, ">"),
        comp_val(dplyr::cur_data(), "intensity_median_SPL", min.intensity.median.spl, ">"),
        comp_val(dplyr::cur_data(), "intensity_max_SPL", min.intensity.highest.spl, ">")),
        .operator = "AND"),

      pass_sb = comp_lgl_vec(
        list(comp_val(dplyr::cur_data(), "sb_ratio_pblk", min.signalblank.median.spl.pblk, ">") | .data$is_istd & !istd.include,
        comp_val(dplyr::cur_data(), "sb_ratio_ublk", min.signalblank.median.spl.ublk, ">") | .data$is_istd & !istd.include,
        comp_val(dplyr::cur_data(), "sb_ratio_sblk", min.signalblank.median.spl.sblk, ">")),
        .operator = "AND"),

      pass_cva = comp_lgl_vec(
        list(
        comp_val(dplyr::cur_data(), "conc_cv_BQC", max.cv.conc.bqc, "<"),
        comp_val(dplyr::cur_data(), "conc_cv_TQC", max.cv.conc.tqc, "<"),
        comp_val(dplyr::cur_data(), "norm_intensity_cv_BQC",  max.cv.normintensity.bqc, "<"),
        comp_val(dplyr::cur_data(), "norm_intensity_cv_TQC", max.cv.normintensity.tqc, "<"),
        comp_val(dplyr::cur_data(), "intensity_cv_BQC",  max.cv.intensity.bqc, "<"),
        comp_val(dplyr::cur_data(), "intensity_cv_TQC", max.cv.intensity.tqc, "<")
        ),
        .operator = "AND"),

      pass_dratio = comp_lgl_vec(
        list(
          comp_val(dplyr::cur_data(), "conc_dratio_sd_bqc_conc", max.dratio.sd.conc.bqc, "<"),
          comp_val( dplyr::cur_data(), "conc_dratio_sd_tqc_conc", max.dratio.sd.conc.tqc, "<"),
          comp_val( dplyr::cur_data(), "conc_dratio_mad_bqc", max.dratio.mad.conc.bqc, "<"),
          comp_val( dplyr::cur_data(), "conc_dratio_mad_tqc", max.dratio.mad.conc.tqc, "<"),
          comp_val(dplyr::cur_data(), "conc_dratio_sd_bqc_normint", max.dratio.sd.normint.bqc, "<"),
          comp_val( dplyr::cur_data(), "conc_dratio_sd_tqc_normint", max.dratio.sd.normint.tqc, "<"),
          comp_val( dplyr::cur_data(), "conc_dratio_mad_bqc_normint", max.dratio.mad.normint.bqc, "<"),
          comp_val( dplyr::cur_data(), "conc_dratio_mad_tqc_normint", max.dratio.mad.normint.tqc, "<")
          ),
        .operator = "AND"),

      pass_missingval = comp_lgl_vec(
        list(comp_val(dplyr::cur_data(), "missing_intensity_prop_spl", max.prop.missing.intensity.spl, "<="),
        comp_val(dplyr::cur_data(), "missing_norm_intensity_prop_spl", max.prop.missing.normintensity.spl, "<="),
        comp_val(dplyr::cur_data(), "missing_conc_prop_spl", max.prop.missing.conc.spl, "<=")),
        .operator = "AND")
    )
  ##tictoc::toc()
  # Check if linearity criteria are defined
  metrics_qc_local<- metrics_qc_local |> mutate(pass_linearity = NA)

  if (resp_criteria_defined){
    #browser()
    if (is.numeric(response.curves.select)) {

      rqc_r2_col_names <- names(metrics_qc_local)[which(stringr::str_detect(names(metrics_qc_local), "r2_rqc"))]
      rqc_r2_col <-  rqc_r2_col_names[response.curves.select]
      rqc_slope_col_names <- names(metrics_qc_local)[which(stringr::str_detect(names(metrics_qc_local), "slopenorm_rqc"))]
      rqc_slope_col <-  rqc_slope_col_names[response.curves.select]
      rqc_y0_col_names <- names(metrics_qc_local)[which(stringr::str_detect(names(metrics_qc_local), "y0norm_rqc"))]
      rqc_y0_col <-  rqc_y0_col_names[response.curves.select]

      if (any(is.na(rqc_r2_col))) cli::cli_abort(cli::col_red("The specified RQC curve index exceeds the available range. There are only {length(rqc_r2_col_names)} RQC curves in the dataset. Please adjust the indices set via `response.curves.select`"))

    } else if(is.character(response.curves.select)){
        rqc_r2_col <- paste0("r2_rqc_", response.curves.select)
        rqc_slope_col <- paste0("slopenorm_rqc_", response.curves.select)
        rqc_y0_col <- paste0("y0norm_rqc_", response.curves.select)
        missing_curves <- response.curves.select[!(rqc_r2_col %in% names(metrics_qc_local))]
        if (length(missing_curves) > 0) {
          cli::cli_abort(cli::col_red("The following response curves are not defined in the metadata: {paste(missing_curves, collapse=', ')}. Please adjust the curve ids or ensure correct identifiers."))
        }
    } else {
      cli::cli_abort(cli::col_red("The `response.curves.select` must be specified as either numeric indices or identifiers provided as strings, corresponding to the RQC curve(s)."))
    }

    # Summarize metrics across curves (columns) based on specified criteria

    # Determine the function to apply based on response.curves.summary
    fun_r2 <- case_match(
      response.curves.summary,
      "worst" ~ "safe_min",
      "best" ~ "safe_max",
      .default = response.curves.summary
    )

    fun_slope <- case_match(
      response.curves.summary,
      c("worst") ~ "safe_min",
      c("best") ~ "safe_max",
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
        rqc_slope__sum__ = purrr::pmap_dbl(
          across(all_of(rqc_slope_col)), ~ do.call(fun_slope, list(c(...), na.rm = TRUE))
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
            !is.na(rqc_r2__sum__) | (!is.na(min.slope.response) & !is.na(max.slope.response) & !is.na(max.yintercept.response)),
            (rqc_r2__sum__ > min.rsquare.response | is.na(min.rsquare.response)) &
              (rqc_slope__sum__ > min.slope.response | is.na(min.slope.response)) &
              (rqc_slope__sum__ <= max.slope.response | is.na(max.slope.response)) &
              (rqc_y0__sum__ < max.yintercept.response | is.na(max.yintercept.response)),
            NA
          )
        )
    }

  }

  # Check if filter has been previously set and if it should be overwritten
  #browser()
  if(!replace_existing & nrow(data@metrics_qc) > 0){
    metrics_old <- data@metrics_qc |>
      select(
        any_of(c(
          "feature_id",
          "batch_id",
          qc_pass_before = "qc_pass",
          pass_lod_before = "pass_lod",
          pass_sb_before = "pass_sb",
          pass_cva_before = "pass_cva",
          pass_linearity_before = "pass_linearity",
          pass_dratio_before =  "pass_dratio",
          pass_missingval_before = "pass_missingval"
        ))
      )

    metrics_qc_local <- metrics_qc_local |> full_join(metrics_old, by = "feature_id")

    prev_filters <- list()
    if(!all(is.na(metrics_qc_local$pass_missingval)) & !all(is.na(metrics_qc_local$pass_missingval_before))) prev_filters <- append(prev_filters, "D-ratio")
    if(!all(is.na(metrics_qc_local$pass_lod)) & !all(is.na(metrics_qc_local$pass_lod_before))) prev_filters <- append(prev_filters, "Min-Intensity")
    if(!all(is.na(metrics_qc_local$pass_sb)) & !all(is.na(metrics_qc_local$pass_sb_before))) prev_filters <- append(prev_filters, "Signal-to-Blank")
    if(!all(is.na(metrics_qc_local$pass_linearity)) & !all(is.na(metrics_qc_local$pass_linearity_before))) prev_filters <- append(prev_filters, "Linearity")
    if(!all(is.na(metrics_qc_local$pass_cva)) & !all(is.na(metrics_qc_local$pass_cva_before))) prev_filters <- append(prev_filters, "%CV")
    if(!all(is.na(metrics_qc_local$pass_dratio)) & !all(is.na(metrics_qc_local$pass_dratio_before))) prev_filters <- append(prev_filters, "D-ratio")
    if(length(prev_filters) > 0){
      cli::cli_alert_warning(cli::col_yellow(glue::glue("Following previously set QC filters were replaced : {glue::glue_collapse(prev_filters, sep = ', ', last = ', and ')}")))
    }
     # Determine QC flags for each feature based on given criteria, excluding RQC

    metrics_qc_local <- metrics_qc_local |>
      mutate(
        pass_lod = if_else(is.na(.data$pass_lod), .data$pass_lod_before, .data$pass_lod),
        pass_sb = if_else(is.na(.data$pass_sb), .data$pass_sb_before, .data$pass_sb),
        pass_cva = if_else(is.na(.data$pass_cva), .data$pass_cva_before, .data$pass_cva),
        pass_linearity = if_else(is.na(.data$pass_linearity), .data$pass_linearity_before, .data$pass_linearity),
        pass_dratio = if_else(is.na(.data$pass_dratio), .data$pass_dratio_before, .data$pass_dratio),
        pass_missingval = if_else(is.na(.data$pass_missingval), .data$pass_missingval_before, .data$pass_missingval)
    )
  }

  metrics_qc_local <- metrics_qc_local |>
      mutate(
        qc_pass =
          (
            (is.na(.data$pass_lod) | .data$pass_lod) &
              (is.na(.data$pass_sb) | .data$pass_sb) &
              (is.na(.data$pass_cva) | .data$pass_cva) &
              (is.na(.data$pass_linearity) | .data$pass_linearity) &
              (is.na(.data$pass_dratio) | .data$pass_dratio) &
              (is.na(.data$pass_missingval) | .data$pass_missingval) &
              (is.na(.data$valid_feature) | .data$valid_feature)
          ) |
          (
            .data$feature_id %in% features.to.keep
          )
      )

    #TODO: deal with invalid integrations (as defined by user in metadata)
    d_filt <- metrics_qc_local |>
      filter(.data$qc_pass)

    if(!qualifier.include) d_metrics_temp <- metrics_qc_local |> filter(.data$is_quantifier)
    if(!istd.include) d_metrics_temp <- metrics_qc_local |> filter(!.data$is_istd)

    n_all_quant <- get_feature_count(data, isistd = FALSE, isquantifier = TRUE)
    n_all_qual <- get_feature_count(data, isistd = FALSE, isquantifier = FALSE)

    n_filt_quant <- nrow(d_metrics_temp |>  filter(.data$in_data, .data$is_quantifier, .data$qc_pass))
    n_filt_qual <- nrow(d_metrics_temp |>  filter(.data$in_data, !.data$is_quantifier, .data$qc_pass))

    if(!replace_existing & nrow(data@metrics_qc) > 0){
      n_filt_quant_before <- nrow(d_metrics_temp |>  filter(.data$in_data, .data$is_quantifier, .data$qc_pass_before))
      n_filt_qual_before <- nrow(d_metrics_temp |>  filter(.data$in_data,!.data$is_quantifier, .data$qc_pass_before))
      qc_pass_prev <- sum(d_metrics_temp$qc_pass_before, na.rm = TRUE)
    }

    n_istd_quant <- get_feature_count(data, isistd = TRUE, isquantifier = TRUE)
    n_istd_qual <- get_feature_count(data, isistd = TRUE, isquantifier = FALSE)

    if (!istd.include) {
      n_filt_quant <- nrow(d_filt |>  filter(.data$in_data,!.data$is_istd, .data$is_quantifier))
      n_filt_qual <- nrow(d_filt |>  filter(.data$in_data, !.data$is_istd, !.data$is_quantifier))

      if(!replace_existing & nrow(data@metrics_qc) > 0){
        n_filt_quant_before <- nrow(d_metrics_temp|>  filter(.data$in_data,!.data$is_istd, .data$is_quantifier, .data$qc_pass_before))
        n_filt_qual_before <- nrow(d_metrics_temp |>  filter(.data$in_data,!.data$is_istd, !.data$is_quantifier, .data$qc_pass_before))
      }
    }

    qc_pass_now <- sum(metrics_qc_local$qc_pass, na.rm = TRUE)

  if(qualifier.include){
    if(!replace_existing & nrow(data@metrics_qc) > 0){

      cli::cli_alert_success(cli::col_green(glue::glue("\rQC filter criteria were added to existing: {n_filt_quant} (before {n_filt_quant_before}) of {n_all_quant} quantifier and {n_filt_qual} of {n_all_qual} qualifier features meet QC criteria ({if_else(!istd.include, 'excluding the', 'including the')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)")))
    } else
      cli::cli_alert_success(cli::col_green(glue::glue("\rQC filter criteria were defined: {n_filt_quant} of {n_all_quant} quantifier and {n_filt_qual} of {n_all_qual} qualifier features meet QC criteria ({if_else(!istd.include, 'excluding the', 'including the')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)")))
  } else {
    if(!replace_existing & nrow(data@metrics_qc) > 0){
      cli::cli_alert_success(cli::col_green(glue::glue("\rQC filter criteria were added to existing: {n_filt_quant} (before {n_filt_quant_before}) of {n_all_quant} quantifier features meet QC criteria ({if_else(!istd.include, 'excluding the', 'including the')} {n_istd_quant} quantifier ISTD features)")))
    } else {
      cli::cli_alert_success(cli::col_green((glue::glue("\rNew QC filter criteria were defined: {n_filt_quant} of {n_all_quant} quantifier features meet QC criteria ({if_else(!istd.include, 'excluding the', 'including the')} {n_istd_quant} quantifier ISTD features)."))))
    }
  }
  if (!qualifier.include) d_filt <- d_filt |> filter(.data$is_quantifier)
  if (!istd.include) d_filt <- d_filt |> filter(!.data$is_istd)

  data@is_filtered <- TRUE
  data@metrics_qc <- metrics_qc_local |> select(-dplyr::ends_with("before"))


  data@dataset_filtered <- data@dataset |>
    dplyr::right_join(d_filt |> filter(.data$qc_pass) |> dplyr::select("feature_id"), by = "feature_id")

  data
}
