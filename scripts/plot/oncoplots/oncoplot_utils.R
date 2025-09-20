library(dplyr)

shorten_consequence <- function(data) {
  data <- data |>
    mutate(across(everything(), ~ gsub("missense_variant", "Missense", .))) |>
    mutate(across(everything(), ~ gsub("stop_gained", "Nonsense", .))) |>
    mutate(across(everything(), ~ gsub("splice_acceptor_variant", "Splice site", .))) |>
    mutate(across(everything(), ~ gsub("splice_donor_variant", "Splice site", .))) |>
    mutate(across(everything(), ~ gsub("frameshift_variant", "Frameshift", .))) |>
    mutate(across(everything(), ~ gsub("inframe_deletion", "Inframe indel", .))) |>
    mutate(across(everything(), ~ gsub("inframe_insertion", "Inframe indel", .))) |>
    mutate(across(everything(), ~ gsub("stop_lost", "Stop codon lost", .))) |>
    mutate(across(everything(), ~ gsub("protein_altering_variant", "Inframe indel", .))) |>
    mutate(across(everything(), ~ gsub("start_lost", "Start codon lost", .)))

  return(data)
}

get_order <- function(dataframe, variable) {
  order <- dataframe |>
    group_by({{ variable }}) |>
    count({{ variable }}) |>
    arrange(n, desc({{ variable }})) |>
    pull({{ variable }})
  return(order)
}

add_missing_samples <- function(sample_list, dataframe, label) {
  colnames_df <- colnames(dataframe)
  missing_samples <- setdiff(sample_list, unique(dataframe[[colnames_df[2]]]))
  assign("missing_samples", missing_samples, envir = .GlobalEnv)

  filler <- dataframe[2, 1][[colnames_df[1]]]

  missing_samples_df <- tibble(
    !!colnames_df[1] := filler,
    !!colnames_df[2] := missing_samples,
    !!colnames_df[3] := label
  )

  dataframe <- bind_rows(dataframe, missing_samples_df)

  return(dataframe)
}

remove_duplicates <- function(dataframe) {
  colnames_df <- colnames(dataframe)

  duplicates <- dataframe |>
    group_by(!!sym(colnames_df[2]), !!sym(colnames_df[1])) |>
    filter(duplicated(!!sym(colnames_df[2])) | duplicated(!!sym(colnames_df[2]), fromLast = TRUE))
  assign("duplicates", duplicates, envir = .GlobalEnv)

  dataframe_unique <- dataframe |>
    group_by(!!sym(colnames_df[2]), !!sym(colnames_df[1])) |>
    mutate(is_duplicate = duplicated(!!sym(colnames_df[2])) | duplicated(!!sym(colnames_df[2]), fromLast = TRUE)) |>
    mutate(!!sym(colnames_df[3]) := if_else(is_duplicate, "Multi hit", !!sym(colnames_df[3]))) |>
    filter(!is_duplicate | row_number() == 1) |>
    select(-is_duplicate)

  return(dataframe_unique)
}

expand_dataframe <- function(dataframe, label) {
  colnames_df <- colnames(dataframe)

  all_combinations <- expand_grid(
    !!sym(colnames_df[2]) := unique(dataframe[[colnames_df[2]]]),
    !!sym(colnames_df[1]) := unique(dataframe[[colnames_df[1]]])
  )

  variants_expanded <- all_combinations |>
    left_join(dataframe, by = c(colnames_df[1], colnames_df[2])) |>
    mutate(!!sym(colnames_df[3]) := ifelse(is.na(!!sym(colnames_df[3])), label, !!sym(colnames_df[3])))

  return(variants_expanded)
}
