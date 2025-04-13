#' Calculate Average Molecular Weight from Chemical Formulas
#'
#' Calculates the average molecular weight of one or more chemical formulas, based on the natural isotopic distribution of elements.
#' The calculation uses the \pkg{enviPat} package to retrieve isotopic masses and abundances, and computes the weighted mean of the isotopic distribution.
#'
#' This function calculates the *average molecular weight* (not monoisotopic mass), which reflects the average of the isotopic distribution found in nature.
#'
#' Isotopes can be specified explicitly in the formula. Atomic mass numbers for isotopes must be enclosed in square brackets (e.g., `[13]C` for carbon-13).
#' Deuterium must be written as `D` instead of `[3]H`.
#'
#' @param formula A character vector of one or more chemical formulas to process.
#'
#' @return A numeric vector of average molecular weights, one for each formula.
#'
#' @details
#' The function uses the \pkg{enviPat} package to validate and parse chemical formulas, calculate isotopic patterns,
#' and determine the average molecular weight based on weighted means of isotopic abundances.
#'
#' @references
#' Loos, M., Gerber, C., Corona, F., Hollender, J., & Singer, H. (2015).
#' Accelerated Isotope Fine Structure Calculation Using Pruned Transition Trees.
#' *Analytical Chemistry*, 87(11), 5738–5744. \doi{10.1021/acs.analchem.5b00941}
#'
#' @examples
#' calc_average_molweight(c("C6H12O6", "[13]C6H12O6", "C8H10N4O2", "D2O"))
#'
#' @export



calc_average_molweight <- function(formula) {

  # Use enviPat to calculate the molecular weight

  if(length(formula) == 0)
    cli_abort(col_red("No chemical formula provided. Please provide on or more valid chemical formula."))

  rlang::check_installed("enviPat")
  # isotopes was obtained via data(isotopes, package = "enviPat") and saved as internal dataset
  formula_checked <- enviPat::check_chemform(isotopes = isotopes, chemforms = formula)

  if(any(formula_checked$warning))
    cli_abort(col_red("Following invalid chemical formula defined: {glue::glue_collapse(formula_checked[formula_checked$warning, ]$new_formula, sep = ', ')}. Please verify feature metadata."))

  # Calculate the molecular weight
  pattern <- enviPat::isopattern(
    isotopes = isotopes,
    chemforms  = formula,
    threshold = 0.0001,
    plotit = FALSE,
    verbose = FALSE,
    charge = FALSE,
    emass = 0.00054858,
    algo = 1
  )

  # Calculate average mass by mean of the isotopic distribution
  weighted_means <- purrr::map_dbl(pattern, function(mat) {
    df <- as.data.frame(mat)
    stats::weighted.mean(df[[1]], df[[2]])
  })

  unname(weighted_means)
}
