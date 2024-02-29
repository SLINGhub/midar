
#' Get list of analyses classified as technical outliers
#'
#' @description
#' Retrieves IDs os analysis technical outliers, i.e, analyses with systematic differences in multiple features.
#' This is based on SD or MAD fences from the mean or median, respectively, from the principal components of a PCA, or of the relative mean abundance (RMA)
#'
#' @param data MidarExperiment object
#' @param variable Feature variable used for outlier detection
#' @param summarize_fun Function used to summarize the features, either "pca" based on PCA, or "rma" based on mean relative abundance (RMA) of all features
#' @param outlier_detection Outlier detection method, either based on SD ("sd") or MAD ("mad")
#' @param fence_multiplicator Multiplicator for SD or MAD, respectively.
#' @param pca_component PCA component to be used
#' @param log_transform Log-transform data for outlier detection
#' @return MidarExperiment object
#'
#' @export
#'
#' @return ggplot2 object
#' @export

analysis_outlier_detection <- function(data,
                                  variable = c("feature_intensity", "feature_norm_intensity", "feature_conc"),
                                  summarize_fun = c("pca", "rma"),
                                  outlier_detection = c("sd", "mad"),
                                  fence_multiplicator,
                                  pca_component,
                                  log_transform = TRUE,
                                  print_outliers = TRUE)
{


  variable <- rlang::arg_match(variable)
  summarize_fun <- rlang::arg_match(summarize_fun)
  outlier_detection <- rlang::arg_match(outlier_detection)

  variable_s <- rlang::sym(variable)
  pc_x = rlang::sym(paste0(".fittedPC",pca_component))

  if(summarize_fun == "rma") stop("Relative Mean Abundance has not yet been implemented. Please use 'pca'")

  d_wide <- data@dataset_QC_filtered  |>
    filter(.data$qc_type %in% c("BQC", "TQC", "SPL")) |>
    filter(!.data$is_istd) |>
    dplyr::select("analysis_id", "qc_type", "batch_id", "feature_name", {{variable}})

  d_filt <- d_wide |>
    tidyr::pivot_wider(id_cols = "analysis_id", names_from = "feature_name", values_from = {{variable}})

  m_raw <- d_filt  |>
    tibble::column_to_rownames("analysis_id") |>
    dplyr::select(where(~!any(is.na(.)))) |>
    as.matrix()

  if(log_transform) m_raw <- log2(m_raw)
  pca_res <- prcomp(m_raw, scale = TRUE, center = TRUE)

  d_metadata <- d_wide  %>% dplyr::select("analysis_id", "qc_type", "batch_id") |> dplyr::distinct()
  pca_annot <- pca_res |> broom::augment(d_metadata)

  if(summarize_fun == "sd"){
    d_outlier <- pca_annot |> filter(.data$.fittedPC1 > (mean(.data$.fittedPC1 ) + fence_multiplicator * sd(.data$.fittedPC1 )) |
                                     .data$.fittedPC1 < (mean(.data$.fittedPC1 ) - fence_multiplicator * sd(.data$.fittedPC1 )))

  } else {
    d_outlier <- pca_annot |> filter(.data$.fittedPC1 > (median(.data$.fittedPC1 ) + fence_multiplicator * mad(.data$.fittedPC1 )) |
                                       .data$.fittedPC1 < (median(.data$.fittedPC1 ) - fence_multiplicator * mad(.data$.fittedPC1 )))
  }


  if(nrow(d_outlier) > 0) {
    data@dataset <- data@dataset |>
      mutate(outlier_technical = analysis_id %in% d_outlier$analysis_id,
             outlier_technical_note = if_else(outlier_technical, glue::glue("PCA, {fence_multiplicator} x {summarize_fun}"), NA_character_))
    data@has_outliers_tech <- TRUE
    data@excl_outliers_tech <- TRUE
    data@dataset_QC_filtered |> filter(run_id < 0)
  }

  if(print_outliers){
    writeLines(glue::glue_collapse(d_outlier$analysis_id, sep = ", ", width = 80, last = ", and "))
   }

  writeLines(crayon::green(glue::glue("\u2713 {nrow(d_outlier)} analyses/samples were classified as technical outlier(s). Please (re)apply 'apply_qc_filter()' and use 'clear_outlier()' to clear all outlier classifications. \n")))

  data

}

#' Clear all analysis/sample outlier classifications
#'
#' @description MidarExperiment object
#' @param data  description
#' @return MidarExperiment object
#'
#' @export
clear_outlier <- function(data){
  data@dataset$outlier_technical <- FALSE
  data@dataset$outlier_technical_note <- NA_character_
  data@has_outliers_tech <- FALSE
  data@excl_outliers_tech <- FALSE
  writeLines(crayon::green(glue::glue("All analysis/sample outlier classifications were cleared. Please reapply 'apply_qc_filter()'. \n")))
  data
}
