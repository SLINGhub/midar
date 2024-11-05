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
            norm_intensity_cv_SPL = cv(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE) * 100
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
             conc_dratio_sd_bqc = sd(.data$feature_conc[.data$qc_type == "BQC"]) / sd(.data$feature_conc[.data$qc_type == "SPL"]),
             conc_dratio_sd_tqc = sd(.data$feature_conc[.data$qc_type == "TQC"]) / sd(.data$feature_conc[.data$qc_type == "SPL"]),
             conc_dratio_mad_bqc = mad(.data$feature_conc[.data$qc_type == "BQC"]) / mad(.data$feature_conc[.data$qc_type == "SPL"]),
             conc_dratio_mad_tqc = mad(.data$feature_conc[.data$qc_type == "TQC"]) / mad(.data$feature_conc[.data$qc_type == "SPL"])
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


#' Filter dataset according to QC and other criteria
#' @description
#' Filter dataset according to QC parameter criteria, remove features that are internal standards (ISTDs) or not annotated as quantifier (optional).
#' Exclude features and analyses that were annotated as not valid in the metadata (valid_feature, valid_analysis).
#'
#' Missing qc parameter values, e.g. no mean BQC intensity, because no BQC values present in the data, are handled as following:
#' - If a filter is applied (e.g. `intensity.min.bqc.min`) but qc value is NA then this feature will be excluded.
#' - If a filter is no applied then no matter if qc value is defined or NA, no filtering will be applied
#'
#' Note: When `istd.include` is TRUE, then `signalblank.median.sblk.min` and `signalblank.median.sblk.min` are ignored for ISTDs, because these blanks are defined as containing ISTDs.
#'
#' @param data MidarExperiment object
#' @param overwrite Add/change previously defined qc filters or define a new set of filters. Default is `FALSE`.
#' @param batch_medians Use median of batch-wise derived QC values. Default is FALSE.
#' @param qualifier.include Include qualifier features
#' @param istd.include Include Internal Standards (ISTD). If set to FALSE, meaning ISTDs will be included, then `min_signal_blank_ratio` is ignored, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#
#' @param missing.intensity.spl.prop.max NA Proportion of missing raw intensities
#' @param missing.norm_intensity.spl.prop.max NA Proportion of missing normalized intensities
#' @param missing.conc.spl.prop.max NA Proportion of missing final concentrations
#' @param outlier.technical.exlude Remove samples classified as outliers
#' @param intensity.min.bqc.min Minimum median feature intensity of BQC
#' @param intensity.min.tqc.min Minimum median feature intensity of TQC
#' @param intensity.min.spl.min Minimum median feature intensity of study samples (SPL)
#' @param intensity.median.bqc.min Minimum median feature intensity of BQC
#' @param intensity.median.tqc.min Minimum median feature intensity of TQC
#' @param intensity.median.spl.min Minimum median feature intensity of study samples (SPL)
#' @param intensity.max.spl.min Minimum maximun feature intensity oof study samples (SPL)
#' @param cv.conc.bqc.max = Maximum %CV of BQC
#' @param cv.conc.tqc.max Maximum %CV of TQC
#' @param cv.intensity.bqc.min Maximum %CV of BQC
#' @param cv.intensity.tqc.min Maximum %CV of TQC

#' @param dratio.conc.bqc.sd.max D-ratio defined as CV_BQC/CV_SPL, based on sd
#' @param dratio.conc.tqc.sd.max D-ratio defined as CV_TQC/CV_SPL, based on sd
#' @param dratio.conc.bqc.mad.max D-ratio defined as CV_BQC/CV_SPL, based on mad
#' @param dratio.conc.tqc.mad.max D-ratio defined as CV_TQC/CV_SPL, based on mad
#' @param signalblank.median.pblk.min = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Process Blank (PBLK)
#' @param signalblank.median.ublk.min = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Unprocessed Blank (UBLK)
#' @param signalblank.median.sblk.min = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Solvent Blank (SBLK)

#' @param response.curve.id A single value or vector specifying the identifiers or indices of the response curves to be used for filtering. This can be either numeric indices or text identifiers.
#' @param response.curve.summary A string indicating how to summarize the metrics of multiple selected curves, if applicable. Must be either 'mean', 'median', 'best', or 'worst'.
#' @param response.rsquare.min The minimum acceptable R squared value for the response curve(s). Refer to the details section for more information.
#' @param response.slope.min The minimum normalized slope for the response curve(s). Refer to the details section for further explanation.
#' @param response.slope.max The maximum normalized slope for the response curve(s). Refer to the details section for additional details.
#' @param response.yintersect.max The maximum allowable normalized y0 intercept for the response curve(s).

#' @param features_to_keep A vector specifying the features that should be retained, regardless of whether they meet the specified quality control filtering criteria.

#' @details
#'
#' Response curves
#'
#' A normalized slope of 1 indicates a proportional relationship between the response and the sample amount. Negative slopes reflect an inverse relationship, where the response decreases as the sample amount increases.
#' The normalized y-intercept represents the value of the response when the sample amount is zero. A value of 1 corresponds to the maximum measured response on the curve.
#'
#' @return MidarExperiment object
#' @export

# TODO: Reporting of qc filters applied on NA data (currently returns FALSE= Exclude when qc value is NA)
## TODO:  Handle feature with all being NA in SPL or QC or all.
qc_apply_feature_filter <- function(data = NULL,
                            overwrite = TRUE,
                            batch_medians = FALSE,
                            qualifier.include = FALSE,
                            istd.include = FALSE,
                            missing.intensity.spl.prop.max  = NA,
                            missing.norm_intensity.spl.prop.max  = NA,
                            missing.conc.spl.prop.max  = NA,
                            intensity.min.bqc.min = NA,
                            intensity.min.tqc.min = NA,
                            intensity.min.spl.min = NA,
                            intensity.median.bqc.min = NA,
                            intensity.median.tqc.min = NA,
                            intensity.median.spl.min = NA,
                            intensity.max.spl.min = NA,
                            signalblank.median.pblk.min = NA,
                            signalblank.median.ublk.min = NA,
                            signalblank.median.sblk.min = NA,
                            cv.conc.bqc.max = NA,
                            cv.conc.tqc.max = NA,
                            cv.intensity.bqc.min = NA,
                            cv.intensity.tqc.min = NA,
                            dratio.conc.bqc.sd.max = NA,
                            dratio.conc.tqc.sd.max = NA,
                            dratio.conc.bqc.mad.max = NA,
                            dratio.conc.tqc.mad.max = NA,
                            response.curve.id = NA,
                            response.curve.summary = NA,
                            response.rsquare.min = NA,
                            response.slope.min = NA,
                            response.slope.max = NA,
                            response.yintersect.max = NA,
                            outlier.technical.exlude = FALSE,
                            features_to_keep = NULL) {

  check_data(data)
   # Check if RQC curve ID is defined when r2 is set
  if (all(is.na(response.curve.id))) {
    if(!is.na(response.rsquare.min) | !is.na(response.slope.min) | !is.na(response.slope.max) |!is.na(response.yintersect.max) | !is.na(response.curve.summary))
      cli::cli_abort(cli::col_red("No response curves selected. Please set the curves using `response.curve.id`, or remove `response.___` arguments to proceed without response filters."))
  } else {
    if(is.na(response.rsquare.min) & is.na(response.slope.min) & is.na(response.slope.max) & is.na(response.yintersect.max)) {
      cli::cli_abort(cli::col_red("No response filters were defined. Please set the appropriate `response.___` arguments, remove `response.curve.id`, or set it to NA."))
    } else {
      if (length(unique(response.curve.id)) > 1){
        if(is.na(response.curve.summary)){
          cli::cli_abort(cli::col_red("Please define `response.curve.summary` to define how the results from different curves should be summarized for filtered, Must be either o either 'mean', 'median', 'best' or 'worst'."))
        } else {
          rlang::arg_match(response.curve.summary, c("mean", "median", "best", "worst"))
          if (nrow(data@annot_responsecurves) == 0)
            cli::cli_abort(cli::col_red("No response curves are defined in the metadata. Please either remove `response.curve.id` and any response filters, or reprocess with updated metadata"))
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


  # Check if feature_ids defind with features_to_keep are present in the dataset
  if (!is.null(features_to_keep)) {
    keepers_not_defined <- setdiff(features_to_keep, unique(data@dataset$feature_id))
    txt <- glue::glue_collapse(keepers_not_defined, sep = ", ", last = ", and ")
    if (length(keepers_not_defined) > 0) cli::cli_abort(glue::glue("Following features defined via `features_to_keep` are not present in this dataset: {txt}"))
  }

  # Save QC filter criteria to MidarExperiment object
  # TODO: fix some of the param below ie. features_to_keep

  # Store the QC filter criteria in MidarExperiment object
  data@parameters_processing <- data@parameters_processing |>
    mutate(
      outlier.technical.exlude = outlier.technical.exlude,
      missing.intensity.spl.prop.max  = missing.intensity.spl.prop.max,
      missing.norm_intensity.spl.prop.max  = missing.norm_intensity.spl.prop.max,
      missing.conc.spl.prop.max  = missing.conc.spl.prop.max,
      intensity.min.bqc.min = intensity.min.bqc.min,
      intensity.min.tqc.min = intensity.min.tqc.min,
      intensity.min.spl.min = intensity.min.spl.min,
      intensity.median.bqc.min = intensity.median.bqc.min,
      intensity.median.tqc.min = intensity.median.tqc.min,
      intensity.median.spl.min = intensity.median.spl.min,
      intensity.max.spl.min = intensity.max.spl.min,
      cv.conc.bqc.max = cv.conc.bqc.max,
      cv.conc.tqc.max = cv.conc.tqc.max,
      cv.intensity.bqc.min = cv.intensity.bqc.min,
      cv.intensity.tqc.min = cv.intensity.tqc.min,
      dratio.conc.bqc.sd.max = dratio.conc.bqc.sd.max,
      dratio.conc.tqc.sd.max = dratio.conc.tqc.sd.max,
      dratio.conc.bqc.mad.max = dratio.conc.bqc.mad.max,
      dratio.conc.tqc.mad.max = dratio.conc.tqc.mad.max,
      signalblank.median.pblk.min = signalblank.median.pblk.min,
      signalblank.median.ublk.min = signalblank.median.ublk.min,
      signalblank.median.sblk.min = signalblank.median.sblk.min,
      response.curve.id_used_for_filt = list(response.curve.id),
      response.rsquare.min = response.rsquare.min,
      response.slope.min = response.slope.min,
      response.slope.max = response.slope.max,
      response.yintersect.max = response.yintersect.max,
      qualifier.include = qualifier.include,
      istd.include = istd.include,
      features_to_keep = NA
    )

  metrics_qc_local <- data_local@metrics_qc
  ##tictoc::tic()
  metrics_qc_local <- metrics_qc_local |>
    #rowwise() |>  #TODO make it vevtorized
    mutate(
      pass_lod = comp_lgl_vec(
        list(comp_val(dplyr::cur_data(), "intensity_min_BQC", intensity.min.bqc.min, ">"),
        comp_val(dplyr::cur_data(), "intensity_min_TQC", intensity.min.tqc.min, ">"),
        comp_val(dplyr::cur_data(), "intensity_min_SPL", intensity.min.spl.min, ">"),
        comp_val(dplyr::cur_data(), "intensity_median_BQC", intensity.median.bqc.min, ">"),
        comp_val(dplyr::cur_data(), "intensity_median_TQC", intensity.median.tqc.min, ">"),
        comp_val(dplyr::cur_data(), "intensity_median_SPL", intensity.median.spl.min, ">"),
        comp_val(dplyr::cur_data(), "intensity_max_SPL", intensity.max.spl.min, ">")),
        .operator = "AND"),

      pass_sb = comp_lgl_vec(
        list(comp_val(dplyr::cur_data(), "sb_ratio_pblk", signalblank.median.pblk.min, ">") | .data$is_istd & !istd.include,
        comp_val(dplyr::cur_data(), "sb_ratio_ublk", signalblank.median.ublk.min, ">") | .data$is_istd & !istd.include,
        comp_val(dplyr::cur_data(), "sb_ratio_sblk", signalblank.median.sblk.min, ">")),
        .operator = "AND"),

      pass_cva = comp_lgl_vec(
        list(comp_val(dplyr::cur_data(), "conc_cv_BQC", cv.conc.bqc.max, "<"),
        comp_val(dplyr::cur_data(), "conc_cv_TQC", cv.conc.tqc.max, "<"),
        comp_val(dplyr::cur_data(), "intensity_cv_BQC",  cv.intensity.bqc.min, "<"),
        comp_val(dplyr::cur_data(), "intensity_cv_TQC", cv.intensity.tqc.min, "<")),
        .operator = "AND"),

      pass_dratio = comp_lgl_vec(
        list(comp_val(dplyr::cur_data(), "conc_dratio_sd_bqc", dratio.conc.bqc.sd.max, "<"),
          comp_val( dplyr::cur_data(), "conc_dratio_sd_tqc", dratio.conc.tqc.sd.max, "<"),
          comp_val( dplyr::cur_data(), "conc_dratio_mad_bqc", dratio.conc.bqc.mad.max, "<"),
          comp_val( dplyr::cur_data(), "conc_dratio_mad_tqc", dratio.conc.tqc.mad.max, "<")),
        .operator = "AND"),

      pass_missingval = comp_lgl_vec(
        list(comp_val(dplyr::cur_data(), "missing_intensity_prop_spl", missing.intensity.spl.prop.max, "<="),
        comp_val(dplyr::cur_data(), "missing_norm_intensity_prop_spl", missing.norm_intensity.spl.prop.max, "<="),
        comp_val(dplyr::cur_data(), "missing_conc_prop_spl", missing.conc.spl.prop.max, "<=")),
        .operator = "AND")
    )
  ##tictoc::toc()
  # Check if linearity criteria are defined
  metrics_qc_local<- metrics_qc_local |> mutate(pass_linearity = NA)

  if (resp_criteria_defined){
    #browser()
    if (is.numeric(response.curve.id)) {

      rqc_r2_col_names <- names(metrics_qc_local)[which(stringr::str_detect(names(metrics_qc_local), "r2_rqc"))]
      rqc_r2_col <-  rqc_r2_col_names[response.curve.id]
      rqc_slope_col_names <- names(metrics_qc_local)[which(stringr::str_detect(names(metrics_qc_local), "slopenorm_rqc"))]
      rqc_slope_col <-  rqc_slope_col_names[response.curve.id]
      rqc_y0_col_names <- names(metrics_qc_local)[which(stringr::str_detect(names(metrics_qc_local), "y0norm_rqc"))]
      rqc_y0_col <-  rqc_y0_col_names[response.curve.id]

      if (any(is.na(rqc_r2_col))) cli::cli_abort(cli::col_red("The specified RQC curve index exceeds the available range. There are only {length(rqc_r2_col_names)} RQC curves in the dataset. Please adjust the indices set via `response.curve.id`"))

    } else if(is.character(response.curve.id)){
        rqc_r2_col <- paste0("r2_rqc_", response.curve.id)
        rqc_slope_col <- paste0("slopenorm_rqc_", response.curve.id)
        rqc_y0_col <- paste0("y0norm_rqc_", response.curve.id)
        missing_curves <- response.curve.id[!(rqc_r2_col %in% names(metrics_qc_local))]
        if (length(missing_curves) > 0) {
          cli::cli_abort(cli::col_red("The following response curves are not defined in the metadata: {paste(missing_curves, collapse=', ')}. Please adjust the curve ids or ensure correct identifiers."))
        }
    } else {
      cli::cli_abort(cli::col_red("The `response.curve.id` must be specified as either numeric indices or identifiers provided as strings, corresponding to the RQC curve(s)."))
    }

    # Summarize metrics across curves (columns) based on specified criteria

    # Determine the function to apply based on response.curve.summary
    fun_r2 <- case_match(
      response.curve.summary,
      "worst" ~ "safe_min",
      "best" ~ "safe_max",
      .default = response.curve.summary
    )

    fun_slope <- case_match(
      response.curve.summary,
      c("worst") ~ "safe_min",
      c("best") ~ "safe_max",
      .default = response.curve.summary
    )

    fun_y0 <- case_match(
      response.curve.summary,
      c("worst") ~ "safe_max",
      c("best") ~ "safe_min",
      .default = response.curve.summary
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
            !is.na(rqc_r2__sum__) | (!is.na(response.slope.min) & !is.na(response.slope.max) & !is.na(response.yintersect.max)),
            (rqc_r2__sum__ > response.rsquare.min | is.na(response.rsquare.min)) &
              (rqc_slope__sum__ > response.slope.min | is.na(response.slope.min)) &
              (rqc_slope__sum__ <= response.slope.max | is.na(response.slope.max)) &
              (rqc_y0__sum__ < response.yintersect.max | is.na(response.yintersect.max)),
            NA
          )
        )
    }

  }

  # Check if filter has been previously set and if it should be overwritten
  #browser()
  if(!overwrite & nrow(data@metrics_qc) > 0){
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
            .data$feature_id %in% features_to_keep
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

    if(!overwrite & nrow(data@metrics_qc) > 0){
      n_filt_quant_before <- nrow(d_metrics_temp |>  filter(.data$in_data, .data$is_quantifier, .data$qc_pass_before))
      n_filt_qual_before <- nrow(d_metrics_temp |>  filter(.data$in_data,!.data$is_quantifier, .data$qc_pass_before))
      qc_pass_prev <- sum(d_metrics_temp$qc_pass_before, na.rm = TRUE)
    }

    n_istd_quant <- get_feature_count(data, isistd = TRUE, isquantifier = TRUE)
    n_istd_qual <- get_feature_count(data, isistd = TRUE, isquantifier = FALSE)

    if (!istd.include) {
      n_filt_quant <- nrow(d_filt |>  filter(.data$in_data,!.data$is_istd, .data$is_quantifier))
      n_filt_qual <- nrow(d_filt |>  filter(.data$in_data, !.data$is_istd, !.data$is_quantifier))

      if(!overwrite & nrow(data@metrics_qc) > 0){
        n_filt_quant_before <- nrow(d_metrics_temp|>  filter(.data$in_data,!.data$is_istd, .data$is_quantifier, .data$qc_pass_before))
        n_filt_qual_before <- nrow(d_metrics_temp |>  filter(.data$in_data,!.data$is_istd, !.data$is_quantifier, .data$qc_pass_before))
      }
    }

    qc_pass_now <- sum(metrics_qc_local$qc_pass, na.rm = TRUE)

  if(qualifier.include){
    if(!overwrite & nrow(data@metrics_qc) > 0){

      cli::cli_alert_success(cli::col_green(glue::glue("\rQC filter criteria were added to existing: {n_filt_quant} (before {n_filt_quant_before}) of {n_all_quant} quantifier and {n_filt_qual} of {n_all_qual} qualifier features meet QC criteria ({if_else(!istd.include, 'excluding the', 'including the')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)")))
    } else
      cli::cli_alert_success(cli::col_green(glue::glue("\rQC filter criteria were defined: {n_filt_quant} of {n_all_quant} quantifier and {n_filt_qual} of {n_all_qual} qualifier features meet QC criteria ({if_else(!istd.include, 'excluding the', 'including the')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)")))
  } else {
    if(!overwrite & nrow(data@metrics_qc) > 0){
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
    dplyr::right_join(metrics_qc_local |> filter(.data$qc_pass) |> dplyr::select("feature_id"), by = "feature_id") #|>
    #filter( !(.data$outlier_technical & outlier.technical.exlude)) #TODO

  data
}
