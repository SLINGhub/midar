#' Calculate QC metrics for each feature
#'
#' @param data MidarExperiment object
#' @importFrom glue glue
#' @importFrom readr write_csv
#' @importFrom tidyr pivot_wider nest unnest

#' @importFrom purrr map
#' @importFrom broom glance
#' @importFrom dplyr summarise
#' @return MidarExperiment object
#' @export
calculate_qc_metrics <- function(data) {
  # if(!(c("feature_norm_intensity") %in% names(data@dataset))) warning("No normali is not normalized")


  # TODO: remove later when fixed
  if (tolower(data@analysis_type) == "lipidomics") data <- lipidomics_get_lipid_class_names(data)

  ds1 <- data@dataset %>%
    dplyr::filter(.data$qc_type %in% c("SPL", "NIST", "LTR", "BQC", "TQC", "PBLK", "SBLK", "UBLK", "IBLK", "CAL", "STD", "LQC", "MQC", "UQC")) %>%
    dplyr::group_by(.data$feature_name, .data$feature_class) %>%
    dplyr::summarise(
      # PrecursorMz = paste0(unique(.data$precursor_mz), collapse = ","),
      # ProductMz = paste0(unique(.data$product_mz), collapse = ","),
      valid_integration = unique(.data$valid_integration),
      missing_prop_spl = sum(is.na(.data$feature_intensity[.data$qc_type == "SPL"]))/length(.data$feature_intensity[.data$qc_type == "SPL"]),
      na_in_all_spl = all(is.na(.data$feature_conc[.data$qc_type == "SPL"])),
      is_quantifier = unique(.data$is_quantifier),
      is_istd = unique(.data$is_istd),
      norm_istd = unique(.data$norm_istd_feature_name),
      quant_istd = unique(.data$quant_istd_feature_name),
      feature_response_factor = unique(.data$feature_response_factor),
      Int_min_SPL = min_val(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
      Int_min_BQC = min_val(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
      Int_min_TQC = min_val(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
      Int_median_PBLK = median(.data$feature_intensity[.data$qc_type == "PBLK"], na.rm = TRUE),
      Int_median_SPL = median(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
      Int_median_BQC = median(.data$feature_intensity[.data$qc_type == "BQC"], na.rm = TRUE),
      Int_median_TQC = median(.data$feature_intensity[.data$qc_type == "TQC"], na.rm = TRUE),
      Int_median_NIST = median(.data$feature_intensity[.data$qc_type == "NIST"], na.rm = TRUE),
      Int_median_LTR = median(.data$feature_intensity[.data$qc_type == "LTR"], na.rm = TRUE),
      Int_max_SPL = max_val(.data$feature_intensity[.data$qc_type == "SPL"], na.rm = TRUE),
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


  data@metrics_qc <- ds1


  if ("RQC" %in% data@dataset$qc_type) {
    model <- as.formula("feature_intensity ~ relative_sample_amount")
    ds2 <- data@dataset %>%
      dplyr::filter(.data$qc_type %in% c("RQC")) %>%
      dplyr::full_join(data@annot_responsecurves, by = "analysis_id") %>%
      dplyr::group_by(.data$feature_name, .data$rqc_series_id) %>%
      dplyr::filter(!all(is.na(.data$feature_intensity))) %>%
      tidyr::nest() %>%
      mutate(
        models = purrr::map(data, function(x) lm(model, data = x, na.action = na.exclude)),
        # mandel = map(data, \(x) DCVtestkit::calculate_mandel(x, "relative_sample_amount", "feature_intensity")),
        # ppa = map(data, \(x) DCVtestkit::calculate_pra_linear(x, "relative_sample_amount", "feature_intensity")),
        stats = purrr::map(.data$models, function(x) broom::glance(x)),
        model = purrr::map(.data$models, function(x) {
          broom::tidy(x) |>
            select(.data$term, .data$estimate) |>
            pivot_wider(names_from = "term", values_from = "estimate")
        })
      ) %>%
      tidyr::unnest(c("stats", "model")) %>%
      dplyr::mutate(y0rel = .data$`(Intercept)` / .data$relative_sample_amount) |>
      dplyr::select("feature_name", "rqc_series_id", r2 = "r.squared", y0rel = "y0rel") %>%
      tidyr::pivot_wider(names_from = "rqc_series_id", values_from = c("r2", "y0rel"), names_prefix = "rqc_") |>
      ungroup()

    data@metrics_qc <- data@metrics_qc %>%
      dplyr::left_join(ds2, by = "feature_name") |>
      ungroup()
  }

  data
}





#' Filter dataset according to QC and other criteria
#' @description
#' Filter dataset according to QC parameter criteria, remove features that are internal standards (ISTDs) or not annotated as quantifier (optional).
#' Exclude features and analyses that were annotated as not valid in the metadata (valid_integration, valid_analysis).
#'
#' Note: When `istds.include` is TRUE, then `min_signal_blank_ratio` is ignored for ISTDs, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#'
#' @param data MidarExperiment object
#' @param qualifier.include Include qualifier features
#' @param istds.include Include Internal Standards (ISTD). If set to FALSE, meaning ISTDs will be included, then `min_signal_blank_ratio` is ignored, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#
#' @param missingval.spl.prop.max NA Proportion of missing values, default is 0.00
#' @param outlier.technical.exlude Remove samples classified as outliers
#' @param intensity.min.bqc.min Minimum median signal intensity of BQC
#' @param intensity.min.tqc.min Minimum median signal intensity of TQC
#' @param intensity.min.spl.min Minimum median signal intensity of study samples (SPL)
#' @param intensity.median.bqc.min Minimum median signal intensity of BQC
#' @param intensity.median.tqc.min Minimum median signal intensity of TQC
#' @param intensity.median.spl.min Minimum median signal intensity of study samples (SPL)
#' @param intensity_spl_max_min Minimum maximun signal intensity oof study samples (SPL)
#' @param cv.conc.bqc.min = Maximum %CV of BQC
#' @param cv.conc.tqc.min Maximum %CV of TQC
#' @param cv.conc.tqc.min = Maximum %CV of TQC
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

# TODO default
apply_qc_filter <- function(data,
                            qualifier.include = FALSE,
                            istds.include = FALSE,
                            missingval.spl.prop.max = 1.00,
                            intensity.min.bqc.min = NA,
                            intensity.min.tqc.min = NA,
                            intensity.min.spl.min = NA,
                            intensity.median.bqc.min = NA,
                            intensity.median.tqc.min = NA,
                            intensity.median.spl.min = NA,
                            intensity_spl_max_min = NA,
                            signalblank.median.pblk.min = NA,
                            signalblank.median.ublk.min = NA,
                            signalblank.median.sblk.min = NA,
                            cv.conc.bqc.min = NA,
                            cv.conc.tqc.min = NA,
                            cv.intensity.bqc.min = NA,
                            cv.intensity.tqc.min = NA,
                            dratio.conc.bqc.max = NA,
                            dratio.conc.tqc.max = NA,
                            response.rsquare.min = NA,
                            response.yintersect.rel.max = NA,
                            response.curve.id = NA,
                            outlier.technical.exlude = FALSE,
                            features_to_keep = NULL) {
  if ((!is.na(response.rsquare.min)) & is.na(response.curve.id) & nrow(data@annot_responsecurves) > 0) stop("RQC Curve ID not defined! Please set response.curve.id parameter or set response.rsquare.min to NA if you which not to filter based on RQC r2 values.")
  if (((!is.na(response.rsquare.min)) | !is.na(response.curve.id)) & nrow(data@annot_responsecurves) == 0) stop("No RQC curves were defined in the metadata. Please reprocess with updated metadata, or to ignore linearity filtering, remove or set response.curve.id and response.rsquare.min to NA")


  if (nrow(data@metrics_qc) == 0) {
    stop("QC info has not yet been calculated. Please apply 'calculate_qc_metrics' first.")
  }
  if (!is.null(features_to_keep)) {
    keepers_not_defined <- setdiff(features_to_keep, unique(data@dataset$feature_name))
    txt <- glue::glue_collapse(keepers_not_defined, sep = ", ", last = ", and ")
    if (length(keepers_not_defined) > 0) stop(glue::glue("Following defined in features_to_keep are not present in this dataset: {txt}"))
  }
  if(nrow(data@metrics_qc) == 0)
    data <- calculate_qc_metrics(data)

  if (is.na(intensity.min.bqc.min)) intensity.min.bqc.min <- -Inf
  if (is.na(intensity.min.tqc.min)) intensity.min.tqc.min <- -Inf
  if (is.na(intensity.min.spl.min)) intensity.min.spl.min <- -Inf

  if (is.na(intensity.median.bqc.min)) intensity.median.bqc.min <- -Inf
  if (is.na(intensity.median.tqc.min)) intensity.median.tqc.min <- -Inf
  if (is.na(intensity.median.spl.min)) intensity.median.spl.min <- -Inf

  if (is.na(intensity_spl_max_min)) intensity_spl_max_min <- -Inf

  if (is.na(cv.conc.bqc.min)) cv.conc.bqc.min <- Inf
  if (is.na(cv.conc.tqc.min)) cv.conc.tqc.min <- Inf
  if (is.na(cv.intensity.bqc.min)) cv.intensity.bqc.min <- Inf
  if (is.na(cv.intensity.tqc.min)) cv.intensity.tqc.min <- Inf
  if (is.na(dratio.conc.bqc.max)) dratio.conc.bqc.max <- Inf
  if (is.na(dratio.conc.tqc.max)) dratio.conc.tqc.max <- Inf
  if (is.na(signalblank.median.pblk.min)) signalblank.median.pblk.min <- -Inf
  if (is.na(signalblank.median.ublk.min)) signalblank.median.ublk.min <- -Inf
  if (is.na(signalblank.median.sblk.min)) signalblank.median.sblk.min <- -Inf
  if (is.na(response.rsquare.min)) response.rsquare.min <- -Inf
  if (is.na(response.yintersect.rel.max)) response.yintersect.rel.max <- Inf
  if (is.na(missingval.spl.prop.max)) missingval.spl.prop.max <- Inf


  # TODO: fix some of the param below ie. features_to_keep
  data@parameters_processing <- data@parameters_processing |>
    mutate(
      outlier.technical.exlude = outlier.technical.exlude,
      intensity.min.bqc.min = intensity.min.bqc.min,
      intensity.min.tqc.min = intensity.min.tqc.min,
      intensity.min.spl.min = intensity.min.spl.min,
      intensity.median.bqc.min = intensity.median.bqc.min,
      intensity.median.tqc.min = intensity.median.tqc.min,
      intensity.median.spl.min = intensity.median.spl.min,
      intensity_spl_max_min = intensity_spl_max_min,
      cv.conc.bqc.min = cv.conc.bqc.min,
      cv.conc.tqc.min = cv.conc.tqc.min,
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
      missingval.spl.prop.max = missingval.spl.prop.max,
      qualifier.include = qualifier.include,
      istds.include = istds.include,
      features_to_keep = NA
    )


  data@metrics_qc <- data@metrics_qc |>
    mutate(
      pass_lod =
        (is.na(.data$Int_min_BQC) | .data$Int_min_BQC > intensity.min.bqc.min) &
        (is.na(.data$Int_min_TQC) | .data$Int_min_TQC > intensity.min.tqc.min) &
        (is.na(.data$Int_min_SPL) | .data$Int_min_SPL > intensity.min.spl.min) &
        (is.na(.data$Int_median_SPL) | .data$Int_median_BQC > intensity.median.bqc.min) &
        (is.na(.data$Int_median_TQC) | .data$Int_median_TQC > intensity.median.tqc.min) &
        (is.na(.data$Int_median_SPL) | .data$Int_median_SPL > intensity.median.spl.min) &
        (is.na(.data$Int_max_SPL) | .data$Int_max_SPL > intensity_spl_max_min),

      pass_sb = (
        (is.na(.data$SB_Ratio_pblk) | .data$SB_Ratio_pblk > signalblank.median.pblk.min | .data$is_istd & !istds.include) &
        (is.na(.data$SB_Ratio_ublk) | .data$SB_Ratio_ublk > signalblank.median.ublk.min) &
        (is.na(.data$SB_Ratio_sblk) | .data$SB_Ratio_sblk > signalblank.median.sblk.min)),

      pass_cva =
        (is.na(.data$conc_CV_BQC) | .data$conc_CV_BQC < cv.conc.bqc.min) &
        (is.na(.data$conc_CV_TQC) | .data$conc_CV_TQC < cv.conc.tqc.min) &
        (is.na(.data$Int_CV_BQC) | .data$Int_CV_BQC < cv.intensity.bqc.min) &
        (is.na(.data$Int_CV_TQC) | .data$Int_CV_TQC < cv.intensity.tqc.min),

      pass_dratio =
        (is.na(.data$conc_dratio_cv_bqc) | .data$conc_dratio_cv_bqc < dratio.conc.bqc.max) &
        (is.na(.data$conc_dratio_cv_tqc) | .data$conc_dratio_cv_tqc < dratio.conc.tqc.max),

      pass_missingval =
        (is.na(.data$missing_prop_spl) |.data$missing_prop_spl <= missingval.spl.prop.max)

      #pass_no_na = !(.data$na_in_all_spl)
    )

  #pass_linearity = TRUE,

  if (is.numeric(response.curve.id)) {
    rqc_r2_col_names <- names(data@metrics_qc)[which(stringr::str_detect(names(data@metrics_qc), "r2_rqc"))]
    rqc_r2_col <- rqc_r2_col_names[response.curve.id]
    rqc_y0_col_names <- names(data@metrics_qc)[which(stringr::str_detect(names(data@metrics_qc), "y0rel_rqc"))]
    rqc_y0_col <- rqc_y0_col_names[response.curve.id]

    if (is.na(rqc_r2_col)) stop(glue::glue("RQC curve index exceeds the {length(rqc_r2_col_names)} present RQC curves in the dataset. Please check `response.curve.id` value."))
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
            (is.na(.data$valid_integration) | .data$valid_integration)
          ) |
          (
            .data$feature_name %in% features_to_keep
          )
    )

    #TODO: deal with invalid integrations (as defined by user in metadata)


    d_filt <- data@metrics_qc |> filter(qc_pass)

    n_istd_quant <- nrow(data@metrics_qc |> filter(is_istd, .data$is_quantifier))
    n_istd_qual <- nrow(data@metrics_qc |> filter(is_istd, !.data$is_quantifier))

    n_all_quant <- nrow(data@metrics_qc |> filter(.data$is_quantifier))
    n_all_qual <- nrow(data@metrics_qc |> filter(!.data$is_quantifier))

    n_filt_quant <- nrow(d_filt |>  filter(.data$is_quantifier))
    n_filt_qual <- nrow(d_filt |>  filter(!.data$is_quantifier))

    if (!istds.include) {
      d_filt <- d_filt |> filter(!.data$is_istd)

      n_all_quant <- n_all_quant - n_istd_quant
      n_all_qual <- n_all_qual - n_istd_qual

      n_filt_quant <- n_filt_quant - n_istd_quant
      n_filt_qual <- n_filt_qual - n_istd_qual

    }

  if (!qualifier.include) d_filt <- d_filt |> filter(.data$is_quantifier)

  if(qualifier.include)
   writeLines(crayon::green(glue::glue("\u2713 QC filtering applied: {n_filt_quant} of {n_all_quant} quantifier and {n_filt_qual} of {n_all_qual} qualifier features passed QC criteria ({if_else(!istds.include, 'excluding', 'including')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)")))
  else
   writeLines(crayon::green(glue::glue("\u2713 QC filtering applied: {n_filt_quant} of {n_all_quant}  quantifier features passed QC criteria ({if_else(!istds.include, 'excluding', 'including')} {n_istd_quant} ISTDs).")))



  data@dataset_filtered <- data@dataset %>%
    dplyr::right_join(d_filt |> dplyr::select("feature_name"), by = "feature_name") |>
    filter(.data$valid_analysis, !(.data$outlier_technical & outlier.technical.exlude))
  data
}
