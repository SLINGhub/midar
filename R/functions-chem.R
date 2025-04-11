# Function to calculate molecular weight from a chemical formula

calc_averagemolmass <- function(formula) {

  # Use enviPat to calculate the molecular weight
  rlang::check_installed("enviPat")
  data(isotopes, package = "enviPat")
  formula_checked <- enviPat::check_chemform(isotopes = isotopes, chemforms = formula)

  if(any(formula_checked$warning))
    cli_abort(col_red("Following invalid chemical formula defined: {glue::glue_collapse(formula_checked[formula_checked$warning, ]$new_formula, sep = ', ')}. Please verify feature metadata."))

  # Calculate the molecular weight


  pattern <- enviPat::isopattern(
    isotopes = isotopes,
    chemforms  = formula,
    threshold=0.0001,
    plotit=FALSE,
    verbose = FALSE,
    charge=FALSE,
    emass=0.00054858,
    algo=1
  )

  # Calculate average mass by mean of the isotopic distribution
  weighted_means <- purrr::map_dbl(pattern, function(mat) {
    df <- as.data.frame(mat)
    weighted.mean(df[[1]], df[[2]])
  })

  unname(weighted_means)
}
