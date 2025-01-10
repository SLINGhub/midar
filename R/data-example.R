#' @title Load an example MidarExperiment dataset
#' @description
#' Load an example MidarExperiment dataset. Dataset 1 is a small dataset (Burla et al, 2024, see below) and Dataset 2 a larger dataset (Tan et al, 2022).See Details below.
#'
#' @param data MidarExperiment object, optional. Data will be overwritten if provided.
#' @param dataset Dataset type. Either 1 or 2. Default is 1.
#'
#' @return MidarExperiment object
#' @examples
#' myexp <- MidarExperiment()
#' myexp <- data_load_example(myexp)
#' myexp

#' @export
data_load_example <- function(data = NULL, dataset = 1){

  if(dataset != 1) stop("Only dataset 1 is currently available")
  data <- lipidomics_dataset
  check_data(data)
  data
}

