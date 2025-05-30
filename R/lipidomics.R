# TODO: Replace with RGOSLIN !



# Determine Lipid and Transition Class Name from analyte name
# Retrieves "lipid class" from the analyte name and adds it a factor column. Lipid class is defined as group of lipids sharing same head group and modications but with different chain lengths and saturations.
# lipid_class_base are for example all ceramides (including deoxy etc) and PC (including PC-O etc)
# Different transitions of the same lipid class are indicatec in transition class

get_analyte_id <- function(transition_name, remove_nl_transitions) {
  #analyte_id <- str_trim(str_replace(transition_name, "\\[.*?\\]", ""))
  analyte_id <- remove_bracket_substrings(transition_name, remove_nl_transitions= remove_nl_transitions)

  analyte_id <- str_trim(str_replace(analyte_id, "  ", " "))
  analyte_id <- str_trim(str_replace(analyte_id, "\\+", ""))
  #analyte_id <- str_trim(str_replace(analyte_id, "\\[.*?\\]", ""))
  analyte_id <- str_trim(str_replace(analyte_id, "  ", " "))
  analyte_id <- str_trim(str_replace(analyte_id, "\\[ISTD\\]", ""))
  analyte_id <- str_trim(str_replace(analyte_id, "\\[IS\\]", ""))
  analyte_id <- str_trim(str_replace(analyte_id, "\\(ISTD\\)", ""))
  analyte_id <- str_trim(str_replace(analyte_id, "\\(IS\\)", ""))
  analyte_id <- str_trim(str_replace(analyte_id, "ISTD", ""))
  analyte_id <- str_trim(str_replace(analyte_id, "IS", ""))
  analyte_id <- str_replace(analyte_id, "SM (?=\\d{2}:\\d{1})", "SM d")
  analyte_id <- str_replace(analyte_id, "(?<=m\\d{2}:\\d{1})\\/", ";O/")
  analyte_id <- str_replace(analyte_id, "(?<=d\\d{2}:\\d{1})\\/", ";O2/")
  analyte_id <- str_replace(analyte_id, "(?<=t\\d{2}:\\d{1})\\/", ";O3/")
  analyte_id <- str_replace(analyte_id, "(?<=m\\d{2}:\\d{1})$", ";O")
  analyte_id <- str_replace(analyte_id, "(?<=d\\d{2}:\\d{1})$", ";O2")
  analyte_id <- str_replace(analyte_id, "(?<=t\\d{2}:\\d{1})$", ";O3")
  analyte_id <- str_replace(analyte_id, " [dmt]{1}(?=\\d{2}:\\d{1})", " ")
  analyte_id <- str_replace(analyte_id, "^AC", "CAR")
  analyte_id <- str_replace(analyte_id, "^Sph", "SPB")
  analyte_id <- str_replace(analyte_id, "^S1P", "SPBP")
  analyte_id <- str_replace(analyte_id, "^dhCer", "Cer")
  analyte_id <- str_replace(analyte_id, "^Cer1P", "CerP")
  analyte_id <- str_replace(analyte_id, "^S1P", "SPBP")
  analyte_id <- str_replace(analyte_id, "^Hex1Cer", "HexCer")
  analyte_id <- str_replace(analyte_id, "15\\-MHDA", "16:0;15Me")
  analyte_id <- str_replace(analyte_id, "\\-OH", ";OH")
  analyte_id <- str_replace(analyte_id, "\\-Me", ";Me")
  analyte_id <- str_replace(analyte_id, "^COH", "Chol")
  analyte_id <- str_replace(analyte_id, " a ", " (a) ")
  analyte_id <- str_replace(analyte_id, " b", " (b) ")
  analyte_id <- str_replace(analyte_id, " b", " (c) ")
  analyte_id <- str_replace(analyte_id, " ab", " (ab) ")
  analyte_id <- str_replace(analyte_id, " bc", " (bc) ")
  analyte_id <- str_remove(analyte_id, "\\|.*$")

  analyte_id <- remove_leading_round_brackets(analyte_id)

  return(str_squish(analyte_id))
}




#' Get lipid class, species and transition names
#'
#' This function retrieves lipid class, species and transition names from the `feature_id` column and adds them as columns to the dataset.
#'
#' @param tbl A data frame containing a `feature_id` column
#' @param use_as_feature_class Set feature_class to lipid_class
#' @param add_transition_names add transition name and transition group, based on information in square brackets in feature_id
#' @param add_chain_composition add total_c and total_db to the dataset
#' @return MidarExperiment object


parse_lipid_feature_names <- function(tbl, use_as_feature_class = "lipid_class_lcb", add_transition_names = FALSE, add_chain_composition = TRUE) {
  #check_installed("rgoslin")
  use_as_feature_class_s <- rlang::sym(use_as_feature_class)

  dat <- tbl |>
    select("feature_id") |>
    unique()

  dat_temp <- dat

  # cat("Retrieving lipid class/transition names...", fill = FALSE)
  #dat_temp <- dat |> mutate(lipid_class_base = (str_squish(str_extract(.data$feature_id, "[A-z0-9]+[[:blank:]]*"))))
  #dat_temp <- dat_temp |> mutate(lipid_class = (str_squish(str_extract(.data$feature_id, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m])"))))

  # add a "-", except for between name and sphingoid base
  #dat_temp <- dat_temp |> mutate(lipid_class = str_replace(.data$lipid_class, "([[:blank:]]+)([^d|t|m]{1})", "-\\2"))
  #dat_temp <- dat_temp |> mutate(lipid_class_by_lcb = str_extract(.data$feature_id, "[A-z0-9]+[[:blank:]]*([A-Z]{1}|[d|t|m][0-9\\:]{4})"))

  # Clean transition names

  dat_temp <- dat_temp |> mutate(feature_id = str_squish(str_replace(.data$feature_id, "\\s*\\(\\-H2O\\)?\\s*", " [M-H2O]")))

  # Add  transition name  (defined by flanking "[]")
  dat_temp <- dat_temp |> mutate(transition_name = str_extract(.data$feature_id, "(?<=\\[).*?(?=\\])"))
  dat_temp <- dat_temp |> mutate(transition_name = ifelse(is.na(.data$transition_name), "", .data$transition_name))

  # Add  transition class (adds measured product class in square brackets after the analyte class name). In case of NL indicated with [-...] (e.g.  TG 48:2 [-16:1]), [NL] is indicated (TG [NL]). Maybe useful for normalization
  dat_temp <- dat_temp |> mutate(transition_class = ifelse(grepl("\\[\\-|\\[NL", .data$feature_id), "M>M-NL", str_replace(paste("", str_trim(str_extract(.data$feature_id, "\\[[a-zA-Z0-9\\-\\: ]*\\]"))), " NA", .data$transition_name)))

  dat_temp <- dat_temp |> mutate(analyte_id = get_analyte_id(.data$feature_id, remove_nl_transitions = FALSE))

  dat_temp <- covert_lyso_pl(dat_temp)
  dat_temp <- convert_triglycerides(dat_temp)
  dat_temp <- dat_temp |> mutate(analyte_id = normalize_isotope_labels(.data$analyte_id))

  # # Convert the transition names for each lipid_class to a number (e.g. for Cer d18:1 "M>SphB","M>SphB-H2O"  becomes 1,2:  for for Cer m18:1 "M-H2O>SphB", "M-H2O>SphB" also becomes 1, 2)
  # dat_temp <- dat_temp |>
  #   group_by(.data$analyte_id) |>
  #   # ToDo: replace superseeded do function
  #   dplyr::do(
  #     mutate(.data, transition_group_id = match(.$transition_name, levels(as.factor(dat_temp[dat_temp$lipid_class_by_lcb == unique(.$lipid_class_by_lcb), ]$transition_name))))
  #   ) |>
  #   mutate(transition_group_id = if_else(is.na(.data$transition_group_id), 1, .data$transition_group_id)) |>
  #   ungroup()

  dat_temp_rgoslin <- dat_temp |>
    mutate(lipid_name = str_trim(str_replace(.data$analyte_id, "\\[.*?\\]", ""))) |>
    mutate(lipid_name = str_replace(.data$lipid_name, "\\s*\\(.*?\\)", "")) |>
   mutate(lipid_name = str_replace(.data$lipid_name, "^(\\S+\\s+\\S+).*", "\\1"))   #Remove all text after 2nd space
  d_goslin <- suppressMessages(rgoslin::parseLipidNames(dat_temp_rgoslin$lipid_name))

  goslin_err <- d_goslin |>  filter(.data$Grammar == "NOT_PARSEABLE") |> pull("Original.Name")

  goslin_err_features <- dat_temp_rgoslin |> filter(.data$lipid_name %in% goslin_err) |> pull("feature_id")
  if(length(goslin_err_features) > 0){
    cli::cli_alert_warning(col_yellow("The names of {length(goslin_err)} of {nrow(dat_temp)} feature_id(s) could not be parsed: {glue::glue_collapse(goslin_err_features, sep = ', ')}."))
  }

  d_goslin <- d_goslin |>
    #filter(Grammar != "NOT_PARSEABLE") |>
    mutate(lipid_class_lcb = .data$Extended.Species.Name) |>
    mutate(Normalized.Name = if_else(str_detect(.data$Original.Name, " \\-P| P\\-"), str_replace(.data$Normalized.Name, " O\\-", " P-"), .data$Normalized.Name)) |>
    mutate(lipid_class_lcb = if_else(.data$Lipid.Maps.Category == "SP" & str_detect(.data$Normalized.Name, "\\;O[1-3]"),
                                     stringr::str_c(.data$Extended.Species.Name, str_extract(.data$Normalized.Name,"\\;O[1-3]")),.data$lipid_class_lcb)) |>
    mutate(lipid_class_lcb = if_else(.data$Extended.Species.Name == "TG" & str_detect(.data$Normalized.Name, "TG O\\-"),
                                     "TG-O" ,.data$lipid_class_lcb)) |>
    mutate(lipid_class_lcb = if_else(.data$Extended.Species.Name == "DG" & str_detect(.data$Normalized.Name, "DG O\\-"),
                                     "DG-O" ,.data$lipid_class_lcb)) |>
    mutate(lipid_class_lcb = if_else(.data$Extended.Species.Name == "SPBP",
                                     "SPBP" ,.data$lipid_class_lcb)) |>
    mutate(lipid_class_lcb = if_else(.data$Extended.Species.Name == "ST 27:1",
                                     "CE" ,.data$lipid_class_lcb)) |>
    mutate(lipid_class_lcb = if_else(.data$Extended.Species.Name == "ST 27:1;O",
                                     "Chol" ,.data$lipid_class_lcb)) |>
    ungroup() |>
    select(analyte_name = "Normalized.Name", lipid_class = "Extended.Species.Name", "lipid_class_lcb", lipid_class_base = "Lipid.Maps.Category", "chem_formula" = "Sum.Formula", total_c = "Total.C", total_db = "Total.DB")


  dat_temp <- dat_temp |>
    dplyr::bind_cols(d_goslin) |>
    group_by(.data$analyte_id) |>
    mutate(
      transition_group_id = match(
        .data$transition_name,
        unique(.data$transition_name[.data$lipid_class_lcb == dplyr::first(.data$lipid_class_lcb)])
      ),
      transition_group_id = if_else(is.na(.data$transition_group_id), 1L, .data$transition_group_id)
    ) |>
    ungroup()

  if (!add_transition_names) dat_temp<- dat_temp|> select(!any_of(c("transition_group_id", "transition_name")))
  if (!add_chain_composition) dat_temp<- dat_temp|> select(!any_of(c("total_c", "total_db")))

  dat <- tbl |>
    dplyr::select(-any_of(c("Original.Name", "analyte_id", "lipid_class", "lipid_class_lcb", "lipid_class_base", "transition_name", "transition_class", "analyte_name", "transition_group_id",  "chem_formula", "total_c", "total_db"))) |>
    dplyr::left_join(dat_temp, by = c("feature_id")) |>
    dplyr::relocate(any_of(c("lipid_class", "lipid_class_lcb", "lipid_class_base", "transition_name", "transition_group_id", "chem_formula", "total_c", "total_db")), .after = "feature_id") |>
    mutate(feature_class = !!use_as_feature_class_s) |>
    ungroup()

  dat |> ungroup()
}


convert_triglycerides <- function(dat){
  d_tg <- dat |>
    select("analyte_id") |>
    distinct() |>
    mutate(new_nl_name = .data$analyte_id) |>
    mutate(new_nl_name = str_replace(.data$new_nl_name, "DG O-", "DG-O ")) |>
    mutate(new_nl_name = str_replace(.data$new_nl_name, "TG O-", "TG-O ")) |>
    filter(str_detect(analyte_id,
                      "^(DG|DG-O|TG|TG-O)\\s\\d{2}:\\d{1,2}\\s.*?(\\[-\\d{2}:\\d\\]|\\[NL-\\d{2}:\\d\\]|\\[NL\\s\\d{2}:\\d\\])")) |>


    tidyr::separate(.data$new_nl_name, into = c("class", "sum_chain", "nl_chain", "isomer"), sep = " ", remove = FALSE, fill = "right") |>
    mutate(isomer_temp = .data$nl_chain) |>
    mutate(nl_chain = if_else(str_detect(.data$nl_chain, "\\[NL-|\\[\\-"), .data$nl_chain, .data$isomer)) |>
    mutate(isomer = if_else(!str_detect(.data$isomer_temp, "\\[NL-|\\[\\-"), .data$isomer_temp, NA)) |>
    select(-"isomer_temp") |>
    mutate(nl_chain = str_replace(.data$nl_chain, "\\[NL-|\\[\\-|\\]", "")) |>
    mutate(nl_chain = str_replace(.data$nl_chain, "\\]$", "")) |>
    tidyr::separate(.data$sum_chain, into = c("sum_chain_c", "sum_chain_db"), sep = ":", remove = FALSE, convert = TRUE) |>
    tidyr::separate(.data$nl_chain, into = c("nl_chain_c", "nl_chain_db"), sep = ":", remove = FALSE, convert = TRUE ) |>
    mutate(second_c = .data$sum_chain_c - .data$nl_chain_c,
           second_db = .data$sum_chain_db - .data$nl_chain_db) |>
    mutate(new_nl_name = str_squish(glue::glue("{class} {nl_chain_c}:{nl_chain_db}_{second_c}:{second_db} {isomer}"))) |>
    mutate(new_nl_name = str_replace(.data$new_nl_name, " NA", "")) |>
    mutate(new_nl_name = str_replace(.data$new_nl_name, "TG-O ", "TG O-")) |>
    mutate(new_nl_name = str_replace(.data$new_nl_name, "DG-O ", "DG O-"))



  dat <- dat |>
    left_join(d_tg, by = "analyte_id") |>
    mutate(analyte_id = if_else(
      str_detect(.data$analyte_id, "^TG"),
        if_else(str_detect(.data$analyte_id, "^[DT]G\\s*\\d{2}(:\\d{1,2})?(\\s+\\S+)?\\s*\\[(NL)?-"),
                .data$new_nl_name,
                str_remove(.data$analyte_id, " \\[SIM\\]")),
      .data$analyte_id
    )) |>
    select("feature_id", "analyte_id", "transition_name", "transition_class")


  dat
}

covert_lyso_pl <- function(dat){

  d_lyso <- dat |>
    filter(str_detect(.data$analyte_id, "^LPC|^LPE|^LPI|^LPS")) |>
    mutate(analyte_id_temp = str_replace(.data$analyte_id, "(?i)(\\[sn1\\] |-sn1)", " (sn1)")) |>
    mutate(analyte_id_temp = str_replace(.data$analyte_id_temp, "(?i)(\\[sn2\\] |-sn2)", " (sn1)")) |>
    select("feature_id", "analyte_id", "analyte_id_temp") |>
    distinct() |>
    mutate(analyte_id_temp = str_remove(.data$analyte_id_temp, "\\(")) |>
    mutate(analyte_id_temp = str_remove(.data$analyte_id_temp, "\\)")) |>
    tidyr::separate(.data$analyte_id_temp, into = c("class", "chain", "isomer1", "isomer2"), fill = "right", sep = " ", remove = FALSE) |>
    mutate(analyte_new = case_when(
      (.data$isomer1 == "a" | .data$isomer1 == "sn2") & !is.na(.data$isomer1) ~ glue::glue("{class} 0:0/{chain} ({isomer2})"),
      (.data$isomer1 == "b" | .data$isomer1 == "sn1") & !is.na(.data$isomer1) ~ glue::glue("{class} {chain}/0:0 ({isomer2})"),
      .default = glue::glue("{class} {chain} ({isomer2})"))) |>
    mutate(analyte_new = str_remove(.data$analyte_new, fixed(" (NA)")))  |>
    select("feature_id", "analyte_new", "analyte_id")

  dat <- dat |>
    left_join(d_lyso, by = c("feature_id", "analyte_id")) |>
    mutate(analyte_id =
             if_else(
              str_detect(.data$analyte_id, "^LPC|^LPE|^LPI|^LPS"),
              str_squish(.data$analyte_new),
              .data$analyte_id
    )) |>
    select(-"analyte_new")

  dat

}



remove_bracket_substrings <- function(x, remove_nl_transitions) {
  # Helper: remove all brackets + surrounding space
  remove_all_brackets <- function(str) {
    gsub("\\s*\\[[^]]*\\]\\s*", " ", str)
  }

  # Helper: remove brackets not starting with NL- or -
  remove_selective_brackets <- function(str) {
    gsub("\\s*\\[(?!NL-|-)\\S[^]]*\\]", " ", str, perl = TRUE)
  }

  # Apply to each element
  sapply(x, function(str) {
    cleaned <- if (grepl("^(TG|DG)", str) & !remove_nl_transitions) {
      remove_selective_brackets(str)
    } else {
      remove_all_brackets(str)
    }
    str_squish(cleaned)  # final cleanup: trim and collapse multiple spaces
  }, USE.NAMES = FALSE)
}


normalize_isotope_labels <- function(x) {
  x <- stringr::str_squish(x)

  # Convert dash-based comma-separated isotope labels into space-based form
  x_dash <- stringr::str_replace_all(
    x,
    "(?<=\\d)-((?:[dD]\\d{1,2}|\\d{1,2}[A-Z]{1}\\d{0,2})(?:,(?:[dD]\\d{1,2}|\\d{1,2}[A-Z]{1}\\d{0,2}))*)",
    " \\1"
  )
  # Now split by comma into space-separated labels
  x_dash <- stringr::str_replace_all(x_dash, ",", " ")

  # Match isotope labels
  matches <- stringr::str_match_all(
    x_dash,
    "(?<=\\s|^|$)([dD]\\d{1,2}|\\d{1,2}[A-Z]{1}\\d{0,2})(?=\\s|$|$)"
  )

  sapply(seq_along(matches), function(i) {
    found <- unique(na.omit(matches[[i]][, 2]))
    if (length(found) == 0) return(x[i])

    # Normalize labels
    norm <- sapply(found, function(s) {
      if (stringr::str_detect(s, "^[dD]\\d{1,2}$")) {
        toupper(s)  # D7
      } else if (stringr::str_detect(s, "^\\d{1,2}[A-Z]{1}\\d{1,2}$")) {
        s
      } else if (stringr::str_detect(s, "^\\d{1,2}[A-Z]{1}$")) {
        paste0(s, "1")
      } else {
        s
      }
    })
    norm_block <- paste0("[", paste(norm, collapse = ","), "]")

    # Remove original isotope labels from text
    cleaned <- x_dash[i]
    for (pat in found) {
      pat_escaped <- stringr::str_replace_all(pat, "([$$$$])", "\\\\\\1")
      cleaned <- stringr::str_remove_all(
        cleaned,
        paste0("(?<=\\s|^)[$$]?", pat_escaped, "[$$]?(?=\\s|$)")
      )
    }
    cleaned <- stringr::str_squish(cleaned)

    # Extract all square bracket blocks (with preceding spaces)
    after_core <- paste(stringr::str_extract_all(cleaned, "\\[.*?\\]")[[1]], collapse = "")
    # Remove all square bracket blocks (with preceding spaces) to get core
    core <- stringr::str_trim(stringr::str_remove_all(cleaned, "\\[.*?\\]"))

    # Return with normalized isotope block after core, then after_core
    paste0(core, norm_block, " ", after_core)
  })
}




remove_leading_round_brackets <- function(x) {
  str_replace(
    x,
    "^\\s*([^\\s()]+)\\s*\\(\\s*([^()]+?)\\s*\\)",
    "\\1 \\2"
  )
}

