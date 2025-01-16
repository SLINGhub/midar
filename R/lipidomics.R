# TODO: Replace with RGOSLIN !



###################################################################################################################
# Determine Lipid and Transition Class Name from analyte name
###################################################################################################################
# Retrieves "lipid class" from the analyte name and adds it a factor column. Lipid class is defined as group of lipids sharing same head group and modications but with different chain lengths and saturations.
# lipid_class_base are for example all ceramides (including deoxy etc) and PC (including PC-O etc)
# Different transitions of the same lipid class are indicatec in transition class

get_analyte_name <- function(transition_name) {
  analyte_name <- str_trim(str_replace(transition_name, "\\[.*?\\]", ""))
  analyte_name <- str_trim(str_replace(analyte_name, "  ", " "))
  analyte_name <- str_trim(str_replace(analyte_name, "\\+", ""))
  analyte_name <- str_trim(str_replace(analyte_name, "\\[.*?\\]", ""))
  analyte_name <- str_trim(str_replace(analyte_name, "  ", " "))
  return(analyte_name)
}




#' Get lipid class, species and transition names
#'
#' This function retrieves lipid class, species and transition names from the `feature_id` column and adds them as columns to the dataset.
#'
#' @param data MidarExperiment object
#' @param use_as_feature_class Set feature_class to lipid_class
#' @param add_transition_names add transition name and transition group, based on information in square brackets in feature_id
#' @return MidarExperiment object
#' @export


get_lipid_class_names <- function(data = NULL, use_as_feature_class = "lipid_class", add_transition_names = FALSE) {
  check_data(data)
  use_as_feature_class_s <- rlang::sym(use_as_feature_class)

  dat <- data@dataset |>
    select("feature_id") |>
    unique()
  data@dataset <- data@dataset |> dplyr::select(!any_of(c("analyte_name", "lipid_class", "lipid_class_by_lcb", "lipid_class_base", "transition_name", "transition_group")))



  # cat("Retrieving lipid class/transition names...", fill = FALSE)
  dat_temp <- dat |> mutate(lipid_class_base = (str_squish(str_extract(.data$feature_id, "[A-z0-9]+[[:blank:]]*"))))
  dat_temp <- dat_temp |> mutate(lipid_class = (str_squish(str_extract(.data$feature_id, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m])"))))

  # add a "-", except for between name and sphingoid base
  dat_temp <- dat_temp |> mutate(lipid_class = str_replace(.data$lipid_class, "([[:blank:]]+)([^d|t|m]{1})", "-\\2"))
  dat_temp <- dat_temp |> mutate(lipid_class_by_lcb = str_extract(.data$feature_id, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m][0-9\\:]{4})"))
  # Add  transition name  (defined by flanking "[]")
  dat_temp <- dat_temp |> mutate(transition_name = str_extract(.data$feature_id, "(?<=\\[).*?(?=\\])"))
  dat_temp <- dat_temp |> mutate(transition_name = ifelse(is.na(.data$transition_name), "", .data$transition_name))

  # Add  transition class (adds measured product class in square brackets after the analyte class name). In case of NL indicated with [-...] (e.g.  TG 48:2 [-16:1]), [NL] is indicated (TG [NL]). Maybe useful for normalization
  dat_temp <- dat_temp |> mutate(transition_class = ifelse(grepl("\\[\\-|\\[NL", .data$feature_id), "M>M-NL", str_replace(paste("", str_trim(str_extract(.data$feature_id, "\\[[a-zA-Z0-9\\-\\: ]*\\]"))), " NA", .data$transition_name)))

  dat_temp <- dat_temp |> mutate(analyte_name = get_analyte_name(.data$feature_id))

  # Convert the transition names for each lipid_class to a number (e.g. for Cer d18:1 "M>SphB","M>SphB-H2O"  becomes 1,2:  for for Cer m18:1 "M-H2O>SphB", "M-H2O>SphB" also becomes 1, 2)
  dat_temp <- dat_temp |>
    group_by(.data$analyte_name) |>
    # ToDo: replace superseeded do function
    dplyr::do(
      mutate(.data, transition_group = match(.$transition_name, levels(as.factor(dat_temp[dat_temp$lipid_class_by_lcb == unique(.$lipid_class_by_lcb), ]$transition_name))))
    ) |>
    mutate(transition_group = if_else(is.na(.data$transition_group), 1, .data$transition_group)) |>
    ungroup()

  data@dataset <- data@dataset |>
    dplyr::left_join(dat_temp, by = c("feature_id")) |>
    dplyr::relocate("analyte_name", "lipid_class", "lipid_class_by_lcb", "lipid_class_base", "transition_name", "transition_group", .after = "feature_id") |>
    mutate(feature_class = !!use_as_feature_class_s) |>
    ungroup()


  if (!add_transition_names) data@dataset <- data@dataset |> select(!any_of(c("transition_group", "transition_name")))




  data
}
