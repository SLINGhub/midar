% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lipidomics.R
\name{get_lipid_class_names}
\alias{get_lipid_class_names}
\title{Get lipid class, species and transition names}
\usage{
get_lipid_class_names(
  data = NULL,
  use_as_feature_class = "lipid_class",
  add_transition_names = FALSE
)
}
\arguments{
\item{data}{MidarExperiment object}

\item{use_as_feature_class}{Set feature_class to lipid_class}

\item{add_transition_names}{add transition name and transition group, based on information in square brackets in feature_id}
}
\value{
MidarExperiment object
}
\description{
This function retrieves lipid class, species and transition names from the \code{feature_id} column and adds them as columns to the dataset.
}
