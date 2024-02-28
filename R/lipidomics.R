###################################################################################################################
# Determine Lipid and Transition Class Name from Compound name
###################################################################################################################
# Retrieves "lipid class" from the compound name and adds it a factor column. Lipid class is defined as group of lipids sharing same head group and modications but with different chain lengths and saturations.
# lipidClassBase are for example all ceramides (including deoxy etc) and PC (including PC-O etc)
# Different transitions of the same lipid class are indicatec in transition class

get_CompoundName <- function(transition_name){
  compound_name <- str_trim(str_replace(transition_name, "\\[.*?\\]",""))
  compound_name <- str_trim(str_replace(compound_name, "  "," "))
  compound_name <- str_trim(str_replace(compound_name, "\\+",""))
  compound_name <- str_trim(str_replace(compound_name, "\\[.*?\\]",""))
  compound_name <- str_trim(str_replace(compound_name, "  "," "))
  return(compound_name)
}


#' Retrieve lipid name, lipid class and transition from feauture names
#' @param data MidarExperiment object
#' @importFrom dplyr do
#' @return MidarExperiment object
#' @export


add_lipid_class_transition <- function(data){
  dat <- data |> select(.data$feature_name) |> unique()

  #cat("Retrieving lipid class/transition names...", fill = FALSE)
  dat_temp <- dat %>%  mutate(LipidClassBase = (str_squish(str_extract(.data$feature_name, "[A-z0-9]+[[:blank:]]*"))))
  dat_temp <- dat_temp %>%  mutate(LipidClass = (str_squish(str_extract(.data$feature_name, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m])"))))

  # add a "-", except for between name and sphingoid base
  dat_temp <- dat_temp %>% mutate(LipidClass = str_replace(.data$LipidClass, "([[:blank:]]+)([^d|t|m]{1})", "-\\2"))
  dat_temp <- dat_temp %>% mutate(LipidClassSL = str_extract(.data$feature_name, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m][0-9\\:]{4})"))
  # Add  transition name  (defined by flanking "[]")
  dat_temp <- dat_temp %>%  mutate(TransitionName = str_extract(.data$feature_name, "(?<=\\[).*?(?=\\])"))
  dat_temp <- dat_temp %>%  mutate(TransitionName = ifelse(is.na(.data$TransitionName),"",.data$TransitionName))

  # Add  transition class (adds measured product class in square brackets after the compound class name). In case of NL indicated with [-...] (e.g.  TG 48:2 [-16:1]), [NL] is indicated (TG [NL]). Maybe useful for normalization
  dat_temp <- dat_temp %>%  mutate(transitionClass =  ifelse(grepl("\\[\\-|\\[NL", .data$feature_name),"M>M-NL",str_replace(paste("", str_trim(str_extract(.data$feature_name, "\\[[a-zA-Z0-9\\-\\: ]*\\]")))," NA", .data$TransitionName)))

  dat_temp <- dat_temp %>%  mutate(CompoundName = get_CompoundName(.data$feature_name))

  # Convert the transition names for each lipidClass to a number (e.g. for Cer d18:1 "M>SphB","M>SphB-H2O"  becomes 1,2:  for for Cer m18:1 "M-H2O>SphB", "M-H2O>SphB" also becomes 1, 2)
  dat_temp <- dat_temp %>%
    group_by(.data$CompoundName) %>%
    #ToDo: replace superseeded do function
    dplyr::do(
      mutate(.data, TransitionGroup = match(.$TransitionName , levels(as.factor(dat_temp[dat_temp$LipidClassSL == unique(.$LipidClassSL),]$TransitionName))))
    ) %>%
    mutate(TransitioGroup = ifelse(is.na(.data$TransitionGroup), 1, .data$TransitionGroup)) |>
    ungroup()

  data <- data |>
    right_join(dat_temp, by = c("feature_name")) |>
    relocate("CompoundName", "LipidClass", "LipidClassSL", "LipidClassBase", "TransitionName","TransitionGroup", .after = .data$feature_name)

  data
}
