#' Calculate feature quality control (QC) metrics
#' @description
#' For each feature, different QC metrics are calculated, either across the
#' full analysis or as medians of batch-wise calculations
#'
#' @param data MidarExperiment object
#' @param batchwise_median Use median of batch-wise derived QC values.
#' @return MidarExperiment object
#' @export
qc_calculate_metrics <- function(data, batchwise_median ) {
  # if(!(c("feature_norm_intensity") %in% names(data@dataset))) warning("No normali is not normalized")

  # TODO: remove later when fixed
  if (tolower(data@analysis_type) == "lipidomics") data <- lipidomics_get_lipid_class_names(data)
        # All features defined in the metadata
      d_feature_info <- data@annot_features |>
        #filter(.data$valid_feature) |>
        select("valid_feature", "feature_id", "feature_class", "is_quantifier", "is_istd", "norm_istd_feature_id", "quant_istd_feature_id", "response_factor") |>
        mutate(is_istd = feature_id == norm_istd_feature_id)

      # All features MS method parameters
      method_var <- c("method_precursor_mz", "method_product_mz", "method_collision_energy")
      d_method_template <- tibble("precursor_mz" = NA_character_, "product_mz" = NA_character_, "collision_energy" = NA_character_)
      if(any(c(method_var) %in% names(data@dataset_orig))){
        d_method_info <- data@dataset_orig |>
          select(c("feature_id", any_of(method_var))) |>
          group_by(feature_id) |>
          summarise(
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
        dplyr::filter(.data$qc_type %in% c("SPL", "NIST", "LTR", "BQC", "TQC", "PBLK", "SBLK", "UBLK", "IBLK", "CAL", "STD", "LQC", "MQC", "UQC")) %>%
        dplyr::group_by(.data$feature_id) %>%
        dplyr::summarise(
          missing_intensity_prop_spl = sum(is.na(.data$feature_intensity[.data$qc_type == "SPL"]))/length(.data$feature_intensity[.data$qc_type == "SPL"]),
          missing_normintensity_prop_spl = sum(is.na(.data$feature_norm_intensity[.data$qc_type == "SPL"]))/length(.data$feature_norm_intensity[.data$qc_type == "SPL"]),
          missing_conc_prop_spl = sum(is.na(.data$feature_conc[.data$qc_type == "SPL"]))/length(.data$feature_conc[.data$qc_type == "SPL"]),
          na_in_all_spl = all(is.na(.data$feature_conc[.data$qc_type == "SPL"]))
        )

      if (batchwise_median) grp <- c("feature_id", "batch_id") else grp <- c("feature_id")

      d_stats_var <- data@dataset |>
        dplyr::group_by(dplyr::pick(grp)) |>
        summarise(
          Int_min_SPL = min_val(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          Int_min_BQC = min_val(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
          Int_min_TQC = min_val(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
          Int_max_SPL = max_val(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          Int_median_PBLK = median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE),
          Int_median_SPL = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
          Int_median_BQC = median(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
          Int_median_TQC = median(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
          Int_median_NIST = median(.data$feature_intensity[.data$qc_type == "NIST"], na.rm = TRUE),
          Int_median_LTR = median(.data$feature_intensity[.data$qc_type == "LTR"], na.rm = TRUE),

          conc_median_TQC = median(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE),
          conc_median_BQC = median(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE),
          conc_median_SPL = median(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE),
          conc_median_NIST = median(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE),
          conc_median_LTR = median(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE),

          SB_Ratio_q10_pbk = quantile(.data$feature_intensity[.data$qc_type == "SPL"], probs = 0.1, na.rm = TRUE, names = FALSE) / median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE, names = FALSE),
          SB_Ratio_pblk = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE, names = FALSE) / median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE, names = FALSE),
          SB_Ratio_ublk = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE, names = FALSE) / median(.data$feature_intensity[.data$qc_type == "UBLK"], na.rm = TRUE, names = FALSE),
          SB_Ratio_sblk = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE, names = FALSE) / median(.data$feature_intensity[.data$qc_type == "SBLK"], na.rm = TRUE, names = FALSE),

          Int_CV_TQC = sd(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE) / mean(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE) * 100,
          Int_CV_BQC = sd(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE) / mean(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE) * 100,
          Int_CV_SPL = sd(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE) / mean(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
          normInt_CV_TQC = sd(.data$feature_norm_intensity[.data$qc_type == "TQC"], na.rm = TRUE) / mean(.data$feature_norm_intensity[.data$qc_type == "TQC"], na.rm = TRUE) * 100,
          normInt_CV_BQC = sd(.data$feature_norm_intensity[.data$qc_type == "BQC"], na.rm = TRUE) / mean(.data$feature_norm_intensity[.data$qc_type == "BQC"], na.rm = TRUE) * 100,
          normInt_CV_SPL = sd(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE) / mean(.data$feature_norm_intensity[.data$qc_type == "SPL"], na.rm = TRUE) * 100,

          conc_CV_TQC = sd(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE) / mean(.data$feature_conc[.data$qc_type == "TQC"], na.rm = TRUE) * 100,
          conc_CV_BQC = sd(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE) / mean(.data$feature_conc[.data$qc_type == "BQC"], na.rm = TRUE) * 100,
          conc_CV_SPL = sd(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE) / mean(.data$feature_conc[.data$qc_type == "SPL"], na.rm = TRUE) * 100,
          conc_CV_NIST = sd(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE) / mean(.data$feature_conc[.data$qc_type == "NIST"], na.rm = TRUE) * 100,
          conc_CV_LTR = sd(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE) / mean(.data$feature_conc[.data$qc_type == "LTR"], na.rm = TRUE) * 100,
          conc_dratio_cv_bqc = .data$conc_CV_BQC / .data$conc_CV_SPL,
          conc_dratio_cv_tqc = .data$conc_CV_TQC / .data$conc_CV_SPL
    )

   if (batchwise_median){
     d_stats_var <- d_stats_var |>
        group_by(pick("feature_id")) |>
        summarise(across(is.numeric, median, na.rm = TRUE))
   }


  data@metrics_qc <-
    tibble("feature_id" = sort(union(unique(data@dataset_orig$feature_id),
                                     unique(data@annot_features$feature_id)))) |>
    left_join(d_feature_info, by = "feature_id") |>
    left_join(d_method_info, by = "feature_id") |>
    left_join(d_stats_missingval, by = "feature_id") |>
    left_join(d_stats_var, by = "feature_id") |>
    relocate(feature_id, feature_class, valid_feature, is_quantifier, precursor_mz, product_mz, collision_energy)

  if ("RQC" %in% data@dataset$qc_type) {
    d_rqc_stats <- get_response_curve_stats(data,  limit_to_rqc)
    data@metrics_qc <- data@metrics_qc |>
      dplyr::left_join(d_rqc_stats, by = "feature_id")
  }
  data
}


#' Linear regression statistics of response curves

#' @param data MidarExperiment object
#' @param with_staturation_stats Include PPA and Mandel stats (from {lancer})
#' @param limit_to_rqc If TRUE (default) then only includes RQC qc type
#' @return Tibble with linear regression stats of all curves in a wide format
#' @export

get_response_curve_stats <- function(data, with_staturation_stats = FALSE, limit_to_rqc = FALSE) {
  model <- as.formula("feature_intensity ~ relative_sample_amount")

  browser()
  d_stats  <- data@dataset %>%
    dplyr::inner_join(data@annot_responsecurves, by = "analysis_id") %>%
    dplyr::group_by(.data$feature_id, .data$rqc_series_id) %>%
    dplyr::filter(!all(is.na(.data$feature_intensity))) %>%
    tidyr::nest() %>%
    mutate(
      models = purrr::map(data, function(x) lm(model, data = x, na.action = na.exclude)),
      mandel = map(data, \(x) lancer::calculate_mandel(x, "relative_sample_amount", "feature_intensity")),
      pra = map(data, \(x) lancer::calculate_pra_linear(x, "relative_sample_amount", "feature_intensity")),
      stats = purrr::map(.data$models, function(x) broom::glance(x)),
      model = purrr::map(.data$models, function(x) {
        broom::tidy(x) |>
          select(.data$term, .data$estimate) |>
          pivot_wider(names_from = "term", values_from = "estimate")
      })
    ) %>%
    tidyr::unnest(c("stats", "model", "pra", "mandel")) %>%

    dplyr::mutate(y0rel = .data$`(Intercept)` / .data$relative_sample_amount) |>
    dplyr::select("feature_id", "rqc_series_id", r2 = "r.squared", y0rel = "y0rel", "mandel_stats", "mandel_p_val", "pra") %>%
    tidyr::pivot_wider(names_from = "rqc_series_id", values_from = c("r2", "y0rel", "mandel_stats", "mandel_p_val", "pra"), names_prefix = "rqc_") |>
    ungroup()


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
#' Note: When `istds.include` is TRUE, then `signalblank.median.sblk.min` and `signalblank.median.sblk.min` are ignored for ISTDs, because these blanks are defined as containing ISTDs.
#'
#' @param data MidarExperiment object
#' @param qualifier.include Include qualifier features
#' @param istds.include Include Internal Standards (ISTD). If set to FALSE, meaning ISTDs will be included, then `min_signal_blank_ratio` is ignored, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#
#' @param missing.intensity.spl.prop.max NA Proportion of missing raw intensities
#' @param missing.normintensity.spl.prop.max NA Proportion of missing normalized intensities
#' @param missing.conc.spl.prop.max NA Proportion of missing final concentrations
#' @param outlier.technical.exlude Remove samples classified as outliers
#' @param intensity.min.bqc.min Minimum median signal intensity of BQC
#' @param intensity.min.tqc.min Minimum median signal intensity of TQC
#' @param intensity.min.spl.min Minimum median signal intensity of study samples (SPL)
#' @param intensity.median.bqc.min Minimum median signal intensity of BQC
#' @param intensity.median.tqc.min Minimum median signal intensity of TQC
#' @param intensity.median.spl.min Minimum median signal intensity of study samples (SPL)
#' @param intensity.max.spl.min Minimum maximun signal intensity oof study samples (SPL)
#' @param cv.conc.bqc.max = Maximum %CV of BQC
#' @param cv.conc.tqc.max Maximum %CV of TQC
#' @param cv.intensity.bqc.min Maximum %CV of BQC
#' @param cv.intensity.tqc.min Maximum %CV of TQC

#' @param dratio.conc.bqc.max D-ratio defined as CV_BQC/CV_SPL
#' @param dratio.conc.tqc.max D-ratio defined as CV_TQC/CV_SPL
#' @param signalblank.median.pblk.min = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Process Blank (PBLK)
#' @param signalblank.median.ublk.min = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Unprocessed Blank (UBLK)
#' @param signalblank.median.sblk.min = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Solvent Blank (SBLK)
#' @param response.rsquare.min = Minimum r squared of RQC curve defined under `response.curve.id`
#' @param response.yintersect.rel.max = Minimum relative y0 intersect, whereby 1 refers to a `relative_sample_amount` of 100%. Used to filter for curves that have a good r2 but are flat or even have a negative slope.
#' @param response.curve.id Name of RQC curve as string, or index number of curve to use for filtering (first curve is 1)
#' @param features_to_keep Features that must be kept, even if they did not meet the given QC criteria
#' @return MidarExperiment object
#' @export

# TODO: Reporting of qc filters applied on NA data (currently returns FALSE= Exclude when qc value is NA)
apply_qc_filter <- function(data,
                            qualifier.include = FALSE,
                            istds.include = FALSE,
                            missing.intensity.spl.prop.max  = NA,
                            missing.normintensity.spl.prop.max  = NA,
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
                            dratio.conc.bqc.max = NA,
                            dratio.conc.tqc.max = NA,
                            response.rsquare.min = NA,
                            response.yintersect.rel.max = NA,
                            response.curve.id = NA,
                            outlier.technical.exlude = FALSE,
                            features_to_keep = NULL) {
  if ((!is.na(response.rsquare.min)) & is.na(response.curve.id) & nrow(data@annot_responsecurves) > 0) cli::cli_abort("RQC Curve ID not defined! Please set response.curve.id parameter or set response.rsquare.min to NA if you which not to filter based on RQC r2 values.")
  if (((!is.na(response.rsquare.min)) | !is.na(response.curve.id)) & nrow(data@annot_responsecurves) == 0) cli::cli_abort("No RQC curves were defined in the metadata. Please reprocess with updated metadata, or to ignore linearity filtering, remove or set response.curve.id and response.rsquare.min to NA")


  if (nrow(data@metrics_qc) == 0) {
    cli::cli_abort("QC info has not yet been calculated. Please apply 'calculate_qc_metrics' first.")
  }
  if (!is.null(features_to_keep)) {
    keepers_not_defined <- setdiff(features_to_keep, unique(data@dataset$feature_id))
    txt <- glue::glue_collapse(keepers_not_defined, sep = ", ", last = ", and ")
    if (length(keepers_not_defined) > 0) cli::cli_abort(glue::glue("Following defined in features_to_keep are not present in this dataset: {txt}"))
  }
  if(nrow(data@metrics_qc) == 0)
    data <- calculate_qc_metrics(data)

  # if (is.na(intensity.min.bqc.min)) intensity.min.bqc.min <- NA
  # if (is.na(intensity.min.tqc.min)) intensity.min.tqc.min <- NA
  # if (is.na(intensity.min.spl.min)) intensity.min.spl.min <- NA
  #
  # if (is.na(intensity.median.bqc.min)) intensity.median.bqc.min <- NA
  # if (is.na(intensity.median.tqc.min)) intensity.median.tqc.min <- NA
  # if (is.na(intensity.median.spl.min)) intensity.median.spl.min <- NA
  #
  # if (is.na(intensity.max.spl.min)) intensity.max.spl.min <- NA
  #
  # if (is.na(cv.conc.bqc.max)) cv.conc.bqc.max <- NA
  # if (is.na(cv.conc.tqc.max)) cv.conc.tqc.max <- NA
  # if (is.na(cv.intensity.bqc.min)) cv.intensity.bqc.min <- NA
  # if (is.na(cv.intensity.tqc.min)) cv.intensity.tqc.min <- NA
  # if (is.na(dratio.conc.bqc.max)) dratio.conc.bqc.max <- NA
  # if (is.na(dratio.conc.tqc.max)) dratio.conc.tqc.max <- NA
  # if (is.na(signalblank.median.pblk.min)) signalblank.median.pblk.min <- NA
  # if (is.na(signalblank.median.ublk.min)) signalblank.median.ublk.min <- NA
  # if (is.na(signalblank.median.sblk.min)) signalblank.median.sblk.min <- NA
  # if (is.na(response.rsquare.min)) response.rsquare.min <- NA
  # if (is.na(response.yintersect.rel.max)) response.yintersect.rel.max <- NA
  # if (is.na(missingval.spl.prop.max)) missingval.spl.prop.max <- NA


  # Save QC filters
  # TODO: fix some of the param below ie. features_to_keep
  data@parameters_processing <- data@parameters_processing |>
    mutate(
      outlier.technical.exlude = outlier.technical.exlude,
      missing.intensity.spl.prop.max  = missing.intensity.spl.prop.max,
      missing.normintensity.spl.prop.max  = missing.normintensity.spl.prop.max,
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
      dratio.conc.bqc.max = dratio.conc.bqc.max,
      dratio.conc.tqc.max = dratio.conc.tqc.max,
      signalblank.median.pblk.min = signalblank.median.pblk.min,
      signalblank.median.ublk.min = signalblank.median.ublk.min,
      signalblank.median.sblk.min = signalblank.median.sblk.min,
      response.curve.id_used_for_filt = response.curve.id,
      response.rsquare.min = response.rsquare.min,
      response.yintersect.rel.max = response.yintersect.rel.max,
      qualifier.include = qualifier.include,
      istds.include = istds.include,
      features_to_keep = NA
    )


  # Function to used to compare qc values with criteria and deal with NA
  # Behaviour:
  # value is NA , threshold NA -> NA
  # value is num , threshold NA -> NA
  # value is NA , threshold is Num -> FALSE
  # value is num , threshold is num -> TRUE/FALSE

  # TODO: Add to function description,
  # TODO: make this function public for user to build own?

  comp_val <- function(val, threshold, operator) {
    if (is.na(val))
      if (is.na(threshold)) NA else FALSE
    else
      if (is.na(threshold)) NA else get(operator)(val, threshold)
  }

  # Get the result of AND/OR of boolean elements of vectors, return NA when all NA
  comp_lgl_vec <- function(lgl_vec, .operator){

    v <- lgl_vec[!is.na(lgl_vec)]
    if (length(v) == 0) return(NA)
    if (.operator == "AND")
      all(v)
    else eif (.operator == "OR")
      any(v)

  }

  data@metrics_qc <- data@metrics_qc |>
    rowwise() |>
    mutate(
      pass_lod = comp_lgl_vec(
        c(comp_val(.data$Int_min_BQC, intensity.min.bqc.min, ">"),
        comp_val(.data$Int_min_TQC, intensity.min.tqc.min, ">"),
        comp_val(.data$Int_min_SPL, intensity.min.spl.min, ">"),
        comp_val(.data$Int_median_BQC, intensity.median.bqc.min, ">"),
        comp_val(.data$Int_median_TQC, intensity.median.tqc.min, ">"),
        comp_val(.data$Int_median_SPL, intensity.median.spl.min, ">"),
        comp_val(.data$Int_max_SPL, intensity.max.spl.min, ">")),
        .operator = "AND"),

      pass_sb = comp_lgl_vec(
        c(comp_val(.data$SB_Ratio_pblk, signalblank.median.pblk.min, ">") | .data$is_istd & !istds.include,
        comp_val(.data$SB_Ratio_ublk, signalblank.median.ublk.min, ">") | .data$is_istd & !istds.include,
        comp_val(.data$SB_Ratio_sblk, signalblank.median.sblk.min, ">")),
        .operator = "AND"),


      pass_cva = comp_lgl_vec(
        c(comp_val(.data$conc_CV_BQC, cv.conc.bqc.max, "<"),
        comp_val(.data$conc_CV_TQC, cv.conc.tqc.max, "<"),
        comp_val(.data$Int_CV_BQC,  cv.intensity.bqc.min, "<"),
        comp_val(.data$Int_CV_TQC, cv.intensity.tqc.min, "<")),
        .operator = "AND"),


      pass_dratio = comp_lgl_vec(
        c(comp_val(.data$conc_dratio_cv_bqc, dratio.conc.bqc.max, "<"),
        comp_val( .data$conc_dratio_cv_tqc, dratio.conc.tqc.max, "<")),
        .operator = "AND"),

      pass_missingval = comp_lgl_vec(
        c(comp_val(.data$missing_intensity_prop_spl, missing.intensity.spl.prop.max, "<="),
        comp_val(.data$missing_normintensity_prop_spl, missing.normintensity.spl.prop.max, "<="),
        comp_val(.data$missing_conc_prop_spl, missing.conc.spl.prop.max, "<=")),
        .operator = "AND"),

      #pass_no_na = !(.data$na_in_all_spl)
    )

  if (is.numeric(response.curve.id)) {
    rqc_r2_col_names <- names(data@metrics_qc)[which(stringr::str_detect(names(data@metrics_qc), "r2_rqc"))]
    rqc_r2_col <- rqc_r2_col_names[response.curve.id]
    rqc_y0_col_names <- names(data@metrics_qc)[which(stringr::str_detect(names(data@metrics_qc), "y0rel_rqc"))]
    rqc_y0_col <- rqc_y0_col_names[response.curve.id]

    if (is.na(rqc_r2_col)) cli::cli_abort(glue::glue("RQC curve index exceeds the {length(rqc_r2_col_names)} present RQC curves in the dataset. Please check `response.curve.id` value."))
  } else {
    rqc_r2_col <- paste0("r2_rqc_", response.curve.id)
    rqc_y0_col <- paste0("y0_rqc_", response.curve.id)
  }


  if (rqc_r2_col %in% names(data@metrics_qc)) {
    data@metrics_qc <- data@metrics_qc |>
      mutate(
        pass_linearity = if_else(is.na(!!ensym(rqc_r2_col)), NA, !!ensym(rqc_r2_col) > response.rsquare.min &
          !!ensym(rqc_y0_col) < response.yintersect.rel.max)
      )
  } else {
    data@metrics_qc <- data@metrics_qc |>
      mutate(
        pass_linearity = NA)
  }


    data@metrics_qc <- data@metrics_qc |>
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


    d_filt <- data@metrics_qc |> filter(qc_pass)

    n_all_quant <- nrow(data@metrics_qc |> filter(.data$is_quantifier))
    n_all_qual <- nrow(data@metrics_qc |> filter(!.data$is_quantifier))
    n_filt_quant <- nrow(d_filt |>  filter(.data$is_quantifier))
    n_filt_qual <- nrow(d_filt |>  filter(!.data$is_quantifier))

    n_istd_quant <- nrow(data@metrics_qc |> filter(is_istd, .data$is_quantifier))
    n_istd_qual <- nrow(data@metrics_qc |> filter(is_istd, !.data$is_quantifier))

    if (!istds.include) {
      n_all_quant <- n_all_quant - n_istd_quant
      n_all_qual <- n_all_qual - n_istd_qual
      n_filt_quant <- nrow(d_filt |>  filter(!.data$is_istd, .data$is_quantifier))
      n_filt_qual <- nrow(d_filt |>  filter(!.data$is_istd, !.data$is_quantifier))
    }



  if(qualifier.include)
    cli::cli_alert_success(cli::col_green(glue::glue("QC filtering applied: {n_filt_quant} of {n_all_quant} quantifier and {n_filt_qual} of {n_all_qual} qualifier features passed QC criteria ({if_else(!istds.include, 'excluding the', 'including the')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)")))
  else
    cli::cli_alert_success(cli::col_green((glue::glue("QC filtering applied: {n_filt_quant} of {n_all_quant}  quantifier features passed QC criteria ({if_else(!istds.include, 'excluding the', 'including the')} {n_istd_quant} quantifier ISTD features)."))))

  if (!qualifier.include) d_filt <- d_filt |> filter(.data$is_quantifier)
  if (!istds.include) d_filt <- d_filt |> filter(!.data$is_istd)

  data@dataset_filtered <- data@dataset %>%
    dplyr::right_join(d_filt |> dplyr::select("feature_id"), by = "feature_id") #|>
    #filter( !(.data$outlier_technical & outlier.technical.exlude)) #TODO
  data
}
