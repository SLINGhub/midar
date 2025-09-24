# Currently, ggvenn throws a warning "Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.". 
# This function suppresses that specific warning
ggvenn_nowarning <- function(...) {
  withCallingHandlers(
    ggvenn::ggvenn(...),
    warning = function(w) {
      if (grepl("Using `size` aesthetic for lines was deprecated", w$message)) {
        invokeRestart("muffleWarning")
      }
    }
  )
}
