
#' Linear Regression Statistics of Response Curves
#'
#' This function calculates linear regression statistics (R-squared, slope, and intercept)
#' for each response curve in the provided `MidarExperiment` object. Optionally, it can include
#' additional statistics from the `lancer` package (if installed) when `with_staturation_stats` is set to TRUE.
#'
#' @param data A `MidarExperiment` object containing the dataset and response curve annotations.
#' @param with_staturation_stats Logical, if TRUE, include additional statistics from the `lancer` package.
#'   Note: The `lancer` package must be installed when this argument is set to TRUE.
#' @param limit_to_rqc Logical, if TRUE (default), only include rows with `qc_type == "RQC"`.
#' @param silent_invalid_data Logical, if TRUE suppresses raising an error when
#' required data or metadata are missing, or there is a mismatch between them.

#'
#' @return A tibble with linear regression statistics (`r2`, `slopenorm`, `y0norm`) for each curve,
#'   or `NULL` if no data matches the criteria.
#'
#' @export
get_response_curve_stats <- function(data = NULL,
                                     with_staturation_stats = FALSE,
                                     limit_to_rqc = FALSE,
                                     silent_invalid_data = FALSE) {
  check_data(data)
  get_lm_results <- function(tbl){

    dt <- tbl

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

  d_stats  <- data@dataset

  if (nrow(data@annot_responsecurves) == 0) {
    if(!silent_invalid_data) cli::cli_abort(col_red("No response curve metadata found. Please import the corresponding metadata first."))
    return(NULL)
  }

  if (limit_to_rqc){
    d_stats <- d_stats |>
      dplyr::filter(.data$qc_type == "RQC")
    if (nrow(d_stats) == 0) {
      if(!silent_invalid_data) cli::cli_abort(col_red("No analyses/samples of QC type `RQC` found. Please verify the analysis metadata."))
      return(NULL)
    }
  }

  d_stats  <- d_stats |>
    select("analysis_id", "feature_id", "feature_intensity") |>
    dplyr::inner_join(data@annot_responsecurves, by = "analysis_id")

  n_diff <- setdiff(data@annot_responsecurves$analysis_id, d_stats$analysis_id)

  if(length(n_diff) > 0){
    if(!silent_invalid_data) cli::cli_abort(col_red("One or more analysis IDs in the response curve metadata do not match the dataset. Please verify your metadata."))
    return(NULL)
  }

  d_stats  <- d_stats |>
    dplyr::filter(!all(is.na(.data$feature_intensity))) |>
    dplyr::group_split(.data$feature_id, .data$curve_id)

  d_stats <- map(d_stats, function(x) get_lm_results(x))

  d_stats <- d_stats |> bind_rows() |>
    dplyr::mutate(slopenorm = .data$slope,
                  y0norm = .data$intercept) |>
    dplyr::select("feature_id", "curve_id", r2 = "r.squared", "slopenorm", "y0norm") |>
    tidyr::pivot_wider(names_from = "curve_id", values_from = c("r2", "slopenorm", "y0norm"), names_prefix = "rqc_")

  if (with_staturation_stats){
    if (!rlang::is_installed("lancer")) {
      cli::cli_abort( # nocov start
        c(
          "{.strong Package `lancer`} must be installed when `with_saturation_stats = TRUE`.",
          "It is available from {.url https://github.com/SLINGhub/lancer}.",
          "Install it using {.code pak::pkg_install(\"SLINGhub/lancer\")}."
        ),
        class = "missing_package_error"
      ) # nocov end
    }

    d_stats_lancer <- data@dataset |>
      select("analysis_id", "feature_id", "feature_intensity") |>
      inner_join(data@annot_responsecurves, by = "analysis_id") |>
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
