# Classes and associated methods for addtional additional structures


# Define a new class for assertr_tibble

as_assertr_tibble <- function(x, ...) {
  if (!inherits(x, "data.frame")) {
    stop("x must be a data.frame")
  }
  x = as_tibble(x)
  class(x) <- c("assertr_tibble", class(x))
  x
}

# Define a new print methods for assertr_tibble
#' @export
tbl_sum.assertr_tibble <- function(x, ...) {
#c("Metadata" = "Errors and Warnings")
  "-------------------------------------------------------------------------------------------"
}
#' @export
tbl_format_header.assertr_tibble <- function(x, setup, ...) {
  cli::col_black(names(setup$tbl_sum), "-", setup$tbl_sum)
}

#' @export
tbl_format_footer.assertr_tibble <- function(x, setup, ...) {
  #cli::style_italic(paste0  "-----------------------------------------------------------------------------------------------------"nThe table has ", setup$rows_total, " rows in total."))
  cli::style_italic(paste0("--------------------------------------------------------------------------------------------\nE = Error, W = Warning, W* = Supressed Warning, N = Note", "\n--------------------------------------------------------------------------------------------"))
  #cli::style_italic(paste0("E = Error, W = Warning, W* = Warning Ignored, N = Note", "\n--------------------------------------------------------------------------------------------"))

  }

#' @export
ctl_new_pillar.assertr_tibble <- function(controller, x, width, ..., title = NULL) {
  out <- NextMethod()
  pillar::new_pillar(list(
    title = out$title,
    #type = out$type,
    data = out$data
  ))
}
