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
#' Note: When `istds_exclude` is FALSE, then `min_signal_blank_ratio` is ignored, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#'
#' @param data MidarExperiment object
#' @param missing_prop_spl_max NA Proportion of missing values, default is 0.00
#' @param outlier_technical_exlude Remove samples classified as outliers
#' @param intensity_bqc_min_min Minimum median signal intensity of BQC
#' @param intensity_tqc_min_min Minimum median signal intensity of TQC
#' @param intensity_spl_min_min Minimum median signal intensity of SPL
#' @param intensity_bqc_median_min Minimum median signal intensity of BQC
#' @param intensity_tqc_median_min Minimum median signal intensity of TQC
#' @param intensity_spl_median_min Minimum median signal intensity of SPL
#' @param intensity_spl_max_min Minimum maximun signal intensity of SPL
#' @param cv_conc_min_bqc = Maximum %CV of BQC
#' @param cv_conc_min_tqc Maximum %CV of TQC
#' @param cv_conc_min_tqc = Maximum %CV of TQC
#' @param cv_intensity_min_bqc Maximum %CV of BQC
#' @param cv_intensity_min_tqc Maximum %CV of TQC
#' @param intensity_min_bqc Minimum median signal intensity of TQC
#' @param intensity_min_tqc Minimum median signal intensity of TQC
#' @param intensity_min_spl Minimum median signal intensity of study samples (SPL)

#' @param dratio_conc_max_bqc D-ratio defined as CV_BQC/CV_SPL
#' @param dratio_conc_max_tqc D-ratio defined as CV_TQC/CV_SPL
#' @param signal_blank_min_pblk = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Process Blank (PBLK)
#' @param signal_blank_min_ublk = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Unprocessed Blank (UBLK)
#' @param signal_blank_min_sblk = Signal-to-Blank ratio. Calculated from the median of study samples and the median of the Solvent Blank (SBLK)
#' @param response_rsquare_min = Minimum r squared of RQC curve defined under `response_curve_name`
#' @param response_y0_rel_max = Minimum relative y0 intersect, whereby 1 refers to a `relative_sample_amount` of 100%. Used to filter for curves that have a good r2 but are flat or even have a negative slope.
#' @param response_curve_name Name of RQC curve as string, or index number of curve to use for filtering (first curve is 1)
#' @param qualifier_exclude Remove features where Quantifier is set to FALSE.
#' @param istds_exclude Remove Internal Standards (ISTD). If set to FALSE, meaning ISTDs will be included, then `min_signal_blank_ratio` is ignored, because the the S/B is based on Processed Blanks (PBLK) that contain ISTDs.
#' @param features_to_keep Features that must be kept, even if they did not meet the given QC criteria
#' @return MidarExperiment object
#' @export

# TODO default
apply_qc_filter <- function(data,
                            missing_prop_spl_max = 1.00,
                            intensity_bqc_min_min = NA,
                            intensity_tqc_min_min = NA,
                            intensity_spl_min_min = NA,
                            intensity_bqc_median_min = NA,
                            intensity_tqc_median_min = NA,
                            intensity_spl_median_min = NA,
                            intensity_spl_max_min = NA,
                            signal_blank_min_pblk = NA,
                            signal_blank_min_ublk = NA,
                            signal_blank_min_sblk = NA,
                            cv_conc_min_bqc = NA,
                            cv_conc_min_tqc = NA,
                            cv_intensity_min_bqc = NA,
                            cv_intensity_min_tqc = NA,
                            dratio_conc_max_bqc = NA,
                            dratio_conc_max_tqc = NA,
                            response_rsquare_min = NA,
                            response_y0_rel_max = NA,
                            response_curve_name = NA,
                            outlier_technical_exlude = FALSE,
                            qualifier_exclude = TRUE,
                            istds_exclude = TRUE,
                            features_to_keep = NULL) {
  if ((!is.na(response_rsquare_min)) & is.na(response_curve_name) & nrow(data@annot_responsecurves) > 0) stop("RQC Curve ID not defined! Please set response_curve_name parameter or set response_rsquare_min to NA if you which not to filter based on RQC r2 values.")
  if (((!is.na(response_rsquare_min)) | !is.na(response_curve_name)) & nrow(data@annot_responsecurves) == 0) stop("No RQC curves were defined in the metadata. Please reprocess with updated metadata, or to ignore linearity filtering, remove or set response_curve_name and response_rsquare_min to NA")


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




  if (is.na(intensity_bqc_min_min)) intensity_bqc_min_min <- -Inf
  if (is.na(intensity_tqc_min_min)) intensity_tqc_min_min <- -Inf
  if (is.na(intensity_spl_min_min)) intensity_spl_min_min <- -Inf

  if (is.na(intensity_bqc_median_min)) intensity_bqc_median_min <- -Inf
  if (is.na(intensity_tqc_median_min)) intensity_tqc_median_min <- -Inf
  if (is.na(intensity_spl_median_min)) intensity_spl_median_min <- -Inf

  if (is.na(intensity_spl_max_min)) intensity_spl_max_min <- -Inf

  if (is.na(cv_conc_min_bqc)) cv_conc_min_bqc <- Inf
  if (is.na(cv_conc_min_tqc)) cv_conc_min_tqc <- Inf
  if (is.na(cv_intensity_min_bqc)) cv_intensity_min_bqc <- Inf
  if (is.na(cv_intensity_min_tqc)) cv_intensity_min_tqc <- Inf
  if (is.na(dratio_conc_max_bqc)) dratio_conc_max_bqc <- Inf
  if (is.na(dratio_conc_max_tqc)) dratio_conc_max_tqc <- Inf
  if (is.na(signal_blank_min_pblk)) signal_blank_min_pblk <- -Inf
  if (is.na(signal_blank_min_ublk)) signal_blank_min_ublk <- -Inf
  if (is.na(signal_blank_min_sblk)) signal_blank_min_sblk <- -Inf
  if (is.na(response_rsquare_min)) response_rsquare_min <- -Inf
  if (is.na(response_y0_rel_max)) response_y0_rel_max <- Inf
  if (is.na(missing_prop_spl_max)) missing_prop_spl_max <- Inf


  # TODO: fix some of the param below ie. features_to_keep
  data@parameters_processing <- data@parameters_processing |>
    mutate(
      outlier_technical_exlude = outlier_technical_exlude,
      intensity_bqc_min_min = intensity_bqc_min_min,
      intensity_tqc_min_min = intensity_tqc_min_min,
      intensity_spl_min_min = intensity_spl_min_min,
      intensity_bqc_median_min = intensity_bqc_median_min,
      intensity_tqc_median_min = intensity_tqc_median_min,
      intensity_spl_median_min = intensity_spl_median_min,
      intensity_spl_max_min = intensity_spl_max_min,
      cv_conc_min_bqc = cv_conc_min_bqc,
      cv_conc_min_tqc = cv_conc_min_tqc,
      cv_intensity_min_bqc = cv_intensity_min_bqc,
      cv_intensity_min_tqc = cv_intensity_min_tqc,
      dratio_conc_max_bqc = dratio_conc_max_bqc,
      dratio_conc_max_tqc = dratio_conc_max_tqc,
      signal_blank_min_pblk = signal_blank_min_pblk,
      signal_blank_min_ublk = signal_blank_min_ublk,
      signal_blank_min_sblk = signal_blank_min_sblk,
      response_curve_name_used_for_filt = response_curve_name,
      response_rsquare_min = response_rsquare_min,
      response_y0_rel_max = response_y0_rel_max,
      missing_prop_spl_max = missing_prop_spl_max,
      qualifier_exclude = qualifier_exclude,
      istds_exclude = istds_exclude,
      features_to_keep = NA
    )


  data@metrics_qc <- data@metrics_qc |>
    mutate(
      pass_lod =
        (is.na(.data$Int_min_BQC) | .data$Int_min_BQC > intensity_bqc_min_min) &
        (is.na(.data$Int_min_TQC) | .data$Int_min_TQC > intensity_tqc_min_min) &
        (is.na(.data$Int_min_SPL) | .data$Int_min_SPL > intensity_spl_min_min) &
        (is.na(.data$Int_median_SPL) | .data$Int_median_BQC > intensity_bqc_median_min) &
        (is.na(.data$Int_median_TQC) | .data$Int_median_TQC > intensity_tqc_median_min) &
        (is.na(.data$Int_median_SPL) | .data$Int_median_SPL > intensity_spl_median_min) &
        (is.na(.data$Int_max_SPL) | .data$Int_max_SPL > intensity_spl_max_min),

      pass_sb = (
        (is.na(.data$SB_Ratio_pblk) | .data$SB_Ratio_pblk > signal_blank_min_pblk | .data$is_istd & !istds_exclude) &
        (is.na(.data$SB_Ratio_ublk) | .data$SB_Ratio_ublk > signal_blank_min_ublk) &
        (is.na(.data$SB_Ratio_sblk) | .data$SB_Ratio_sblk > signal_blank_min_sblk)),

      pass_cva =
        (is.na(.data$conc_CV_BQC) | .data$conc_CV_BQC < cv_conc_min_bqc) &
        (is.na(.data$conc_CV_TQC) | .data$conc_CV_TQC < cv_conc_min_tqc) &
        (is.na(.data$Int_CV_BQC) | .data$Int_CV_BQC < cv_intensity_min_bqc) &
        (is.na(.data$Int_CV_TQC) | .data$Int_CV_TQC < cv_intensity_min_tqc),

      pass_dratio =
        (is.na(.data$conc_dratio_cv_bqc) | .data$conc_dratio_cv_bqc < dratio_conc_max_bqc) &
        (is.na(.data$conc_dratio_cv_tqc) | .data$conc_dratio_cv_tqc < dratio_conc_max_tqc),

      pass_missingval =
        (is.na(.data$missing_prop_spl) |.data$missing_prop_spl <= missing_prop_spl_max)

      #pass_no_na = !(.data$na_in_all_spl)
    )

  #pass_linearity = TRUE,

  if (is.numeric(response_curve_name)) {
    rqc_r2_col_names <- names(data@metrics_qc)[which(stringr::str_detect(names(data@metrics_qc), "r2_rqc"))]
    rqc_r2_col <- rqc_r2_col_names[response_curve_name]
    rqc_y0_col_names <- names(data@metrics_qc)[which(stringr::str_detect(names(data@metrics_qc), "y0rel_rqc"))]
    rqc_y0_col <- rqc_y0_col_names[response_curve_name]

    if (is.na(rqc_r2_col)) stop(glue::glue("RQC curve index exceeds the {length(rqc_r2_col_names)} present RQC curves in the dataset. Please check `response_curve_name` value."))
  } else {
    rqc_r2_col <- paste0("r2_rqc_", response_curve_name)
    rqc_y0_col <- paste0("y0_rqc_", response_curve_name)
  }


  if (rqc_r2_col %in% names(data@metrics_qc)) {
    data@metrics_qc <- data@metrics_qc |>
      mutate(
        pass_linearity = if_else(is.na(!!ensym(rqc_r2_col)), NA, !!ensym(rqc_r2_col) > response_rsquare_min &
          !!ensym(rqc_y0_col) < response_y0_rel_max)
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

    if (istds_exclude) {
      d_filt <- d_filt |> filter(!.data$is_istd)

      n_all_quant <- n_all_quant - n_istd_quant
      n_all_qual <- n_all_qual - n_istd_qual

      n_filt_quant <- n_filt_quant - n_istd_quant
      n_filt_qual <- n_filt_qual - n_istd_qual

    }







  if (qualifier_exclude) d_filt <- d_filt |> filter(.data$is_quantifier)








  if(!qualifier_exclude)
   writeLines(crayon::green(glue::glue("\u2713 QC filtering applied: {n_filt_quant} of {n_all_quant} quantifier and {n_filt_qual} of {n_all_qual} qualifier features passed QC criteria ({if_else(istds_exclude, 'excluding', 'including')} {n_istd_quant} quantifier and {n_istd_qual} qualifier ISTD features)")))
  else
   writeLines(crayon::green(glue::glue("\u2713 QC filtering applied: {n_filt_quant} of {n_all_quant}  quantifier features passed QC criteria ({if_else(istds_exclude, 'excluding', 'including')} {n_istd_quant} ISTDs).")))



  data@dataset_filtered <- data@dataset %>%
    dplyr::right_join(d_filt |> dplyr::select("feature_name"), by = "feature_name") |>
    filter(.data$valid_analysis, !(.data$outlier_technical & outlier_technical_exlude))
  data
}
