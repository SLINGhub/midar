#' Plot Retention Time versus Chain Length and Saturation
#'
#' Generates  scatter plots of retention time (RT) versus either chain
#' length, degree of saturation (double bonds), or equivalent carbon number (ECN)
#' of lipid features of diffent feature classes. This
#' visualization can be useful in identifying annotation (peak picking) errors in reversed-phase (RP)-LC lipidomics dataset
#' arising from isotopic,
#' isobaric, isomeric, or unknown interferences.
#'
#' The retention time can be either plotted against the total number of carbon atoms with
#' the total number of double bonds as curves, or opposite, with the total double bond
#' number as x axis and the total number of carbon atoms as curves. Alternatively,
#' the retention time can be plotted against the ECN, which is calculated as
#' \eqn{ECN = C_{total} - ecn_k \times DB_{total}}, where \eqn{ecn_k} is a constant that may need
#' to be adjusted to the specific chromatographic properties. The default value is
#' \eqn{ecn_k = 1.5}.
#'
#' @param data A `MidarExperiment` object
#' @param x_var Variable to use for the x-axis. One ofEither "total_c", "total_db" or "ecn".
#' @param qc_types A character vector of QC types to include in the plot. If `NA`, all
#' @param outliers_highlight Whether to highlight potential outliers in the plot. Default is `TRUE`.
#' @param outlier_residual_min Minimum value for the residuals to be considered an outlier (default is `0.15`). The value corresponds to the RT difference betweem the
#' fitted line and the median RT of the feature. The value is used to flag outliers.
#' @param outlier_print Whether to print the features that are flagged as potential outliers to the console. Default is `TRUE`.
#' @param ecn_k Constant for ECN calculation (ECN = C - ecn_k* DB), see Details. Default is `1.5`.
#' @param robust_regression Whether to use robust regression, which is less sensitive to outlier (default is `TRUE`).
#' @param include_qualifier Whether to include qualifier features.
#' @param cols_page Number of facet columns, representing different feature classes, shown per page (default is `5`).
#' @param point_size Size of the data points. Default is 2
#' @param point_transparency Alpha transparency of the data point. Default is 0.9
#' @param line_transparency Alpha transparency of the regression lines. Default is 0.9
#' @param base_font_size Base font size for the plot.
#'
#' @return A `ggplot2` object representing faceted scatter plots
#'
#' @export

plot_rt_vs_chain <- function(
  data = NULL,
  x_var = c("total_c", "total_db", "ecn"),
  qc_types = NA,
  outliers_highlight = TRUE,
  outlier_residual_min = 0.15,
  outlier_print = TRUE,
  ecn_k = 1.5,
  include_qualifier = FALSE,
  robust_regression = TRUE,
  cols_page = 5,
  point_size = 2,
  point_transparency = 0.9,
  line_transparency = 0.5,
  base_font_size = 8
) {
  # -- Data validation --

  if (inherits(data, "MidarExperiment")) {
    check_data(data)
  } else if (inherits(data, "data.frame")) {
    if (length(unique(data[[1]])) < nrow(data)) {
      cli::cli_abort(col_red(
        "The first column has one or more replicated feature identifier. Please check your data."
      ))
    }

    if (!is.numeric(data[[2]])) {
      cli::cli_abort(col_red(
        "All values in the second column must be numeric, as they are interpreted as retention times. Please check your data."
      ))
    }
    data <- MidarExperiment()
    data@dataset <- tibble(
      analysis_id = "NA",
      qc_type = "SPL",
      feature_id = data[, 1],
      feature_rt = data[, 2]
    )
  } else {
    cli::cli_abort(
      "The input data must be a `MidarExperiment` object or a `data.frame`."
    )
  }

  x_var <- rlang::arg_match(x_var, c("total_c", "total_db", "ecn"))

  d_lipid_info <- data@dataset |>
    dplyr::select(
      "analysis_id",
      "feature_id",
      "qc_type",
      "is_istd",
      "is_quantifier",
      "feature_class",
      "feature_rt"
    ) |>
    filter(.data$qc_type %in% qc_types) |>
    group_by(.data$feature_id) |>
    summarise(
      rt_median = median(.data$feature_rt, na.rm = TRUE),
      lipid_class = dplyr::first(.data$feature_class),
      is_istd = dplyr::first(.data$is_istd),
      is_quantifier = dplyr::first(.data$is_quantifier),
      .groups = "drop"
    ) |>
    parse_lipid_feature_names() |>
    ungroup()

  # -- Determine x and grouping variables --
  if (x_var == "ecn") {
    d_lipid_info <- d_lipid_info |>
      dplyr::mutate(
        ecn = .data$total_c - ecn_k * .data$total_db,
        total_grp = 0,
      )
    x_var <- "ecn"
    group_var <- "total_grp"
    col_var <- "total_db"
    x_title <- glue::glue(
      "Equivalent Carbon Number (ECN) of the chains (ecn_k = {ecn_k})"
    )
  } else if (x_var == "total_c") {
    x_var <- "total_c"
    group_var <- "total_db"
    col_var <- "total_db"
    x_title <- "Total Carbon Number (C) of the chains"
  } else if (x_var == "total_db") {
    x_var <- "total_db"
    group_var <- "total_c"
    col_var <- "total_c"
    x_title <- "Total Double Bond Number (DB) of the chains"
  }

  # -- Remove qualifiers if requested --
  if (!include_qualifier) {
    d_lipid_info <- d_lipid_info |> filter(.data$is_quantifier)
  }

  d_lipid_info$lipid_class_lcb <- factor(d_lipid_info$lipid_class_lcb)

  # -- Fit models and predict fitted values per group --
  fit_models <- function(data, x_var, group_var, remove_outlier_order) {
    if (remove_outlier_order) {
      data <- data |>
        mutate(
          rt_median = if_else(.data$is_outlier_order, NA_real_, .data$rt_median)
        )
    }

    fitted_lines <- data |>
      group_by(across(all_of(group_var)), .data$lipid_class_lcb) |>
      nest() |>
      mutate(
        model = purrr::map(
          data,
          ~ {
            if (n_distinct(.x[[x_var]]) >= 3) {
              n <- nrow(.x)
              weights <- rep(1, n)
              if (n >= 1) {
                weights[1] <- 0.8
              } # Set first value
              if (n >= 2) {
                weights[n] <- 0.6
              } # Set last value (if n >= 2)
              # Robust regression with weights
              tryCatch(
                withCallingHandlers(
                  if (robust_regression) {
                    MASS::rlm(
                      rt_median ~ poly(.x[[x_var]], 2, raw = FALSE),
                      data = .x,
                      weights = weights
                    )
                  } else {
                    lm(rt_median ~ poly(.x[[x_var]], 2, raw = FALSE), data = .x)
                  },
                  warning = function(w) invokeRestart("muffleWarning")
                ),
                error = function(e) NULL,
                warning = function(w) NULL
              )
            } else if (n_distinct(.x[[x_var]]) == 2) {
              tryCatch(
                withCallingHandlers(
                  lm(rt_median ~ .x[[x_var]], data = .x),
                  warning = function(w) invokeRestart("muffleWarning")
                ),
                error = function(e) NULL,
                warning = function(w) NULL
              )
            } else {
              NULL
            }
          }
        )
      ) |>
      filter(!purrr::map_lgl(.data$model, is.null)) |>
      mutate(
        pred = purrr::map2(
          .data$model,
          data,
          ~ suppressWarnings(broom::augment(.x, newdata = .y))
        )
      ) |>
      unnest("pred")
    fitted_lines |> ungroup()
  }

  # -- NEW: Yes detection using MAD on regression residuals --
  # Merge fitted lines with original data to calculate residuals
  d_plot <- d_lipid_info

  avg_rt_median <- d_plot |>
    group_by(.data$lipid_class_lcb, !!sym(x_var), !!sym(group_var)) |>
    summarise(avg_rt = mean(.data$rt_median, na.rm = TRUE)) |> # Calculate average rt_median for the current group
    ungroup() |> # Ungroup to access the entire dataset
    arrange(.data$lipid_class_lcb) |> # Sort by lipid_class to ensure proper order
    group_by(.data$lipid_class_lcb, !!sym(group_var)) |> # Regroup to access the current group
    mutate(
      #prev_avg_rt = dplyr::lag(avg_rt, default = 0),
      next_avg_rt = dplyr::lead(.data$avg_rt, default = Inf)
    ) |> # Access average from the next group
    ungroup() # Remove grouping

  d_plot <- d_plot |>
    left_join(avg_rt_median, by = c("lipid_class_lcb", x_var, group_var)) |>
    group_by(.data$lipid_class_lcb, !!sym(group_var)) |>
    arrange(desc(!!sym(x_var)), .by_group = TRUE, ) |>
    mutate(
      # Flag if current rt_median is less than or equal to the average rt_median of the previous x_var
      is_outlier_order = .data$rt_median >= .data$next_avg_rt,
    ) |>
    ungroup()

  fitted_lines <- fit_models(
    d_plot,
    x_var,
    group_var,
    remove_outlier_order = TRUE
  )

  d_plot <- d_plot |>
    left_join(
      fitted_lines |> select("feature_id", ".fitted"),
      by = "feature_id",
      "lipid_class_lcb"
    ) |>
    group_by(.data$lipid_class_lcb, !!sym(group_var)) |>
    arrange(desc(!!sym(x_var)), .by_group = TRUE, ) |>
    mutate(
      # Flag if current rt_median is less than or equal to the average rt_median of the previous x_var
      residual = .data$rt_median - .data$.fitted,
      is_outlier_resid = abs(.data$residual) > outlier_residual_min,
      is_outlier = .data$is_outlier_resid | .data$is_outlier_order,
      is_outlier_plot = outliers_highlight & .data$is_outlier
    ) |>
    ungroup()

  # For points: open symbol for non-outliers, filled for outliers
  d_plot <- d_plot |>
    mutate(
      # Shape 21 is circle (can be open or filled by fill aesthetic)
      point_shape = 21,
      #point_shape = ifelse(!is_outlier, 24, 21),
      point_fill = ifelse(
        .data$is_outlier,
        as.character(!!sym(col_var)),
        NA_character_
      ),
      #stroke =  ifelse(is_outlier, 0.7, 0.5)
    )

  # If you want a legend for outliers, add a variable:
  d_plot <- d_plot |>
    mutate(
      is_outlier_plot = ifelse(
        !is.na(.data$is_outlier_plot) & .data$is_outlier_plot,
        "Yes",
        "No"
      )
    )

  # Ensure grouping columns are factors for plotting
  d_plot[[group_var]] <- factor(d_plot[[group_var]])
  fitted_lines[[group_var]] <- factor(fitted_lines[[group_var]])
  d_plot[[col_var]] <- factor(d_plot[[col_var]])
  fitted_lines[[col_var]] <- factor(fitted_lines[[col_var]])

  # Colour palette
  alt_bright_dark_colors <- c(
    "#023557",
    "#40bce6",
    "#018c45",
    "#a1f073",
    "#E31A1C",
    "#FDBF6F",
    "#6a2185",
    "#bfaed4",
    "#FF7F00",
    "#AAAAAA",
    "#8B4513"
  )

  n_colors_needed <- d_plot |>
    tidyr::drop_na(!!sym(x_var)) |>
    dplyr::pull(!!sym(col_var)) |>
    unique() |>
    length()

  p <- ggplot(
    d_plot |> drop_na(!!sym(x_var), .data$rt_median),
    aes(x = !!sym(x_var), y = .data$rt_median, group = !!sym(group_var))
  )

  if (x_var == "ecn") {
    p <- p +
      geom_line(
        data = fitted_lines,
        aes(x = !!sym(x_var), y = .data$.fitted),
        color = "#afc7c7",
        inherit.aes = FALSE,
        linewidth = 0.7,
        alpha = line_transparency
      )
  } else {
    p <- p +
      geom_line(
        data = fitted_lines,
        aes(x = !!sym(x_var), y = .data$.fitted, color = !!sym(group_var)),
        inherit.aes = FALSE,
        linewidth = 0.7,
        alpha = line_transparency
      )
  }

  # Open points for non-outliers, filled for outliers using fill aesthetic
  p <- p +
    geom_point(
      aes(
        shape = .data$is_outlier_plot,
        fill = !!sym(col_var),
        color = !!sym(col_var)
      ),
      size = point_size,
      alpha = point_transparency,
      stroke = 0.5
    ) +
    facet_wrap(vars(.data$lipid_class_lcb), ncol = cols_page, scales = "free") +
    scale_x_continuous(
      breaks = function(x) {
        seq(
          floor(min(x, na.rm = TRUE) / 2) * 2,
          ceiling(max(x, na.rm = TRUE) / 2) * 2,
          by = 2
        )
      },
      expand = c(0.1, 0)
    ) +
    scale_y_continuous(
      expand = c(0.1, 0)
    )

  # Set colour and fill scales
  p <- p +
    scale_color_manual(
      values = rep(alt_bright_dark_colors, length.out = n_colors_needed)
    ) +
    scale_fill_manual(
      values = rep(alt_bright_dark_colors, length.out = n_colors_needed),
      na.value = "white"
    ) +
    # Show open symbol for non-outlier, filled for outlier in legend
    scale_shape_manual(values = c("No" = 1, "Yes" = 8)) +
    labs(x = x_title, y = "Median Retention Time", shape = "Outlier") +
    theme_bw(base_size = base_font_size) +
    theme(
      panel.grid.minor = element_blank()
    )

  if (!outliers_highlight) {
    p <- p + ggplot2::guides(shape = "none") # Do not show shape legend
  }

  if (outlier_print) {
    species <- d_plot |> filter(.data$is_outlier) |> pull(.data$feature_id)
    if (length(species) > 0) {
      cli::cli_alert_info(
        "The following features were flagged as potential annotation outliers: {glue::glue_collapse(species, sep = ', ')}"
      )
    } else {
      cli::cli_alert_success("No potential annotation outliers were detected.")
    }
  }
  p
}
